// ABT #1540: thermal network of coils whose sections are stacked along the column (CONTIGUOUS
// sections in the winding window) with a tape lying across the window between them.
//
// Fixtures (tests/testData), all stored designs loaded by the web app:
//   thermal_island_220kw_litz_contiguous_12121.json   E 96/42/26, litz P/S/P/S/P, zero-thickness tapes.
//   thermal_island_litz_rect_contiguous_12121.json    ER 35/20/11, litz primary and one rectangular
//                                                     strip per secondary section, zero-thickness tapes.
//   thermal_island_pq3540_litz_rect_contiguous_1212121.json  PQ 35/40, litz + strips, 25 um tapes.
// The first two are VALID layouts (no copper inside copper, everything inside the window) that the
// network refused with "disconnected node(s)" before the fix.

#include "physical_models/Temperature.h"
#include "constructive_models/Mas.h"
#include "processors/MagneticSimulator.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;
using json = nlohmann::json;

namespace {

const std::string kLitzE96 = "thermal_island_220kw_litz_contiguous_12121.json";
const std::string kLitzAndStripEr35 = "thermal_island_litz_rect_contiguous_12121.json";
const std::string kLitzAndStripPq3540 = "thermal_island_pq3540_litz_rect_contiguous_1212121.json";

// "Touching" for the geometric checks of these tests: the network's own tolerance (0.1 um).
constexpr double kTouching = 1e-7;

json load_abt1540_fixture(const std::string& name) {
    auto path = std::filesystem::path(__FILE__).parent_path() / "testData" / name;
    std::ifstream file(path);
    if (!file.good()) {
        throw std::runtime_error("cannot open test fixture " + path.string());
    }
    json data;
    file >> data;
    return data;
}

// What the web does with a loaded MAS file: mas_autocomplete(mas, false, {}) (WebLibMKF
// `mas_autocomplete` binding, called by loadMasIntoApp); the stored layout is kept.
OpenMagnetics::Mas autocomplete_like_the_web(const json& masJson) {
    OpenMagnetics::Mas mas(masJson);
    return OpenMagnetics::mas_autocomplete(mas, false, json::object());
}

// Gives every tape lying across the window (the insulation sections between the stacked
// conduction sections) a real thickness, and moves each section below a tape down by that much,
// so the layout stays a valid one. `interference` is how far the tape and the turns on either
// side of it are drawn into each other: stored coordinates carry nanometres of rounding, and
// 10 nm is far inside the 0.1 um the network takes as touching.
void give_tapes_across_window_a_thickness(json& masJson, double thickness, double interference) {
    auto& coil = masJson["magnetic"]["coil"];
    const double step = thickness - 2 * interference;
    std::map<std::string, double> shiftOfSection;
    size_t tapesAbove = 0;
    double previousY = std::numeric_limits<double>::max();
    for (auto& section : coil["sectionsDescription"]) {
        double y = section["coordinates"][1];
        if (!(y < previousY)) {
            throw std::runtime_error("fixture sections are not stacked from the top down");
        }
        previousY = y;
        double shift = -static_cast<double>(tapesAbove) * step;
        if (section["type"] == "insulation") {
            shift -= step / 2;
            section["dimensions"][1] = thickness;
            ++tapesAbove;
        }
        shiftOfSection[section["name"]] = shift;
        section["coordinates"][1] = y + shift;
    }
    if (tapesAbove == 0) {
        throw std::runtime_error("fixture has no tape across the window");
    }
    for (auto& layer : coil["layersDescription"]) {
        const std::string sectionName = layer["section"];
        layer["coordinates"][1] = layer["coordinates"][1].get<double>() + shiftOfSection.at(sectionName);
        if (layer["type"] == "insulation" && layer["orientation"] == "contiguous") {
            layer["dimensions"][1] = thickness;
        }
    }
    for (auto& turn : coil["turnsDescription"]) {
        const std::string sectionName = turn["section"];
        turn["coordinates"][1] = turn["coordinates"][1].get<double>() + shiftOfSection.at(sectionName);
    }
}

struct BuiltNetwork {
    std::vector<ThermalNetworkNode> nodes;
    std::vector<ThermalResistanceElement> resistances;
    std::string error;  // what() of the solve, empty when it solved
    ThermalResult result;

    bool has_conduction_edge(size_t a, size_t b) const {
        for (const auto& r : resistances) {
            if (r.type != HeatTransferType::CONDUCTION) continue;
            if ((r.nodeFromId == a && r.nodeToId == b) || (r.nodeFromId == b && r.nodeToId == a)) return true;
        }
        return false;
    }
    // Face of node `a` on its conduction edge to node `b`; throws when there is none.
    ThermalNodeFace face_towards(size_t a, size_t b) const {
        for (const auto& r : resistances) {
            if (r.type != HeatTransferType::CONDUCTION) continue;
            if (r.nodeFromId == a && r.nodeToId == b) return r.quadrantFrom;
            if (r.nodeFromId == b && r.nodeToId == a) return r.quadrantTo;
        }
        throw std::runtime_error("no conduction edge between " + nodes[a].name + " and " + nodes[b].name);
    }
    std::vector<size_t> of_part(ThermalNodePartType part) const {
        std::vector<size_t> indices;
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].part == part) indices.push_back(i);
        }
        return indices;
    }
    double left(size_t i) const { return nodes[i].physicalCoordinates[0] - nodes[i].dimensions.width / 2; }
    double right(size_t i) const { return nodes[i].physicalCoordinates[0] + nodes[i].dimensions.width / 2; }
    double bottom(size_t i) const { return nodes[i].physicalCoordinates[1] - nodes[i].dimensions.height / 2; }
    double top(size_t i) const { return nodes[i].physicalCoordinates[1] + nodes[i].dimensions.height / 2; }
    double mean_turn_temperature(size_t windingIndex) const {
        double sum = 0;
        size_t count = 0;
        for (size_t i : of_part(ThermalNodePartType::TURN)) {
            if (nodes[i].windingIndex.value() != windingIndex) continue;
            sum += result.nodeTemperatures.at(nodes[i].name);
            ++count;
        }
        if (count == 0) throw std::runtime_error("no turn nodes of winding " + std::to_string(windingIndex));
        return sum / static_cast<double>(count);
    }
    std::pair<double, double> turn_temperature_range() const {
        double lowest = std::numeric_limits<double>::max();
        double highest = std::numeric_limits<double>::lowest();
        for (size_t i : of_part(ThermalNodePartType::TURN)) {
            double t = result.nodeTemperatures.at(nodes[i].name);
            lowest = std::min(lowest, t);
            highest = std::max(highest, t);
        }
        return {lowest, highest};
    }
};

// Builds and solves the network with the stated losses: `lossPerWinding[w]` watts in winding w,
// shared between its turns by their length; no dependence on the loss models.
BuiltNetwork build_network(OpenMagnetics::Magnetic magnetic, double ambientTemperature, double coreLosses,
                           const std::vector<double>& lossPerWinding) {
    auto coil = magnetic.get_coil();
    auto turns = coil.get_turns_description().value();
    std::vector<double> lengthPerWinding(lossPerWinding.size(), 0.0);
    for (auto& turn : turns) {
        lengthPerWinding.at(coil.get_winding_index_by_name(turn.get_winding())) += turn.get_length();
    }
    std::vector<WindingLossesPerElement> perTurn;
    double windingLosses = 0;
    for (auto& turn : turns) {
        size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        OhmicLosses ohmic;
        ohmic.set_losses(lossPerWinding.at(windingIndex) * turn.get_length() / lengthPerWinding.at(windingIndex));
        ohmic.set_origin(ResultOrigin::SIMULATION);
        ohmic.set_method_used("fixed loss per winding, shared by turn length");
        WindingLossesPerElement element;
        element.set_name(turn.get_name());
        element.set_ohmic_losses(ohmic);
        perTurn.push_back(element);
        windingLosses += ohmic.get_losses();
    }
    WindingLossesOutput windingLossesOutput;
    windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
    windingLossesOutput.set_method_used("fixed loss per winding, shared by turn length");
    windingLossesOutput.set_winding_losses(windingLosses);
    windingLossesOutput.set_winding_losses_per_turn(perTurn);

    TemperatureConfig config;
    config.ambientTemperature = ambientTemperature;
    config.coreLosses = coreLosses;
    config.windingLosses = windingLosses;
    config.windingLossesOutput = windingLossesOutput;
    config.plotSchematic = false;

    BuiltNetwork network;
    Temperature temperature(magnetic, config);
    try {
        network.result = temperature.calculateTemperatures();
    }
    catch (const std::exception& e) {
        network.error = e.what();
        if (network.error.empty()) network.error = "(exception without a message)";
    }
    network.nodes = temperature.getNodes();
    network.resistances = temperature.getResistances();
    return network;
}

BuiltNetwork build_fixture_network(const json& masJson, double coreLosses, const std::vector<double>& lossPerWinding) {
    auto autocompleted = autocomplete_like_the_web(masJson);
    return build_network(autocompleted.get_magnetic(), 40.0, coreLosses, lossPerWinding);
}

// Nodes of `network` that the tapes lying across the window became: wider than thick.
std::vector<size_t> tapes_across_window(const BuiltNetwork& network) {
    std::vector<size_t> tapes;
    for (size_t i : network.of_part(ThermalNodePartType::INSULATION_LAYER)) {
        if (network.nodes[i].dimensions.width > network.nodes[i].dimensions.height) tapes.push_back(i);
    }
    return tapes;
}

} // namespace

TEST_CASE("ThermalNetwork ABT 1540: a zero-thickness tape across the window is not a node", "[temperature][abt-1540][smoke-test]") {
    // D2. Both fixtures draw every insulation layer with zero thickness: the wraps between the
    // layers of a section as [0, span] and the tapes between the stacked sections as [span, 0].
    // Zero thickness is the winder's boundary marker, not material, so none of them is a node.
    // Only dimensions[0] used to be tested, so the tapes across the window became nodes without
    // height (L_2, L_5, ... in the ticket's error text).
    settings.reset();
    auto fixture = GENERATE(kLitzE96, kLitzAndStripEr35);
    INFO(fixture);
    auto masJson = load_abt1540_fixture(fixture);
    size_t zeroThicknessTapesAcrossWindow = 0;
    for (auto& layer : masJson["magnetic"]["coil"]["layersDescription"]) {
        if (layer["type"] != "insulation") continue;
        double thickness = layer["dimensions"][layer["orientation"] == "contiguous" ? 1 : 0];
        REQUIRE(thickness == 0.0);
        if (layer["orientation"] == "contiguous") ++zeroThicknessTapesAcrossWindow;
    }
    REQUIRE(zeroThicknessTapesAcrossWindow == 4);

    auto network = build_fixture_network(masJson, 1.0, {1.0, 1.0});
    CHECK(network.of_part(ThermalNodePartType::INSULATION_LAYER).empty());
    CHECK(network.of_part(ThermalNodePartType::TURN).size() == masJson["magnetic"]["coil"]["turnsDescription"].size());
    CHECK(network.error == "");
}

TEST_CASE("ThermalNetwork ABT 1540: turns resting on a tape across the window keep their edge to each other", "[temperature][abt-1540][smoke-test]") {
    // D1. The same two layouts with 50 um tapes between the sections. The tape under a row of
    // turns runs along the whole row; it is not between two neighbours of that row, and they
    // conduct to each other like any two touching turns of adjacent layers. The old test took the
    // tape as standing between them (its radial span overlaps theirs) and removed that edge.
    settings.reset();
    auto fixture = GENERATE(kLitzE96, kLitzAndStripEr35);
    INFO(fixture);
    auto masJson = load_abt1540_fixture(fixture);
    give_tapes_across_window_a_thickness(masJson, 50e-6, 10e-9);
    auto network = build_fixture_network(masJson, 1.0, {1.0, 1.0});

    auto tapes = tapes_across_window(network);
    REQUIRE(tapes.size() == 4);
    auto turns = network.of_part(ThermalNodePartType::TURN);
    size_t neighboursOnATape = 0;
    for (size_t a = 0; a < turns.size(); ++a) {
        for (size_t b = a + 1; b < turns.size(); ++b) {
            size_t first = turns[a], second = turns[b];
            // Same row, side by side and touching.
            if (std::abs(network.nodes[first].physicalCoordinates[1] - network.nodes[second].physicalCoordinates[1]) > kTouching) continue;
            double radialGap = std::max(network.left(first), network.left(second)) - std::min(network.right(first), network.right(second));
            if (std::abs(radialGap) > kTouching) continue;
            // Resting on (or under) a tape.
            bool onATape = false;
            for (size_t tape : tapes) {
                if (std::abs(network.bottom(first) - network.top(tape)) < kTouching ||
                    std::abs(network.top(first) - network.bottom(tape)) < kTouching) {
                    onATape = true;
                }
            }
            if (!onATape) continue;
            ++neighboursOnATape;
            INFO(network.nodes[first].name << " | " << network.nodes[second].name);
            CHECK(network.has_conduction_edge(first, second));
        }
    }
    // Counted on the stored drawings. E 96: 8 rows touch a tape; the primary rows have 3 turns
    // (2 neighbour pairs), the secondary rows 2 (their outer layer is shorter): 4 x 2 + 4 x 1.
    // ER 35: 4 primary rows of 4 litz layers touch a tape (3 pairs each); the strips are alone.
    CHECK(neighboursOnATape == 12);
    CHECK(network.error == "");
}

TEST_CASE("ThermalNetwork ABT 1540: a tape across the window conducts to the turns on both of its faces", "[temperature][abt-1540][smoke-test]") {
    // D3. With a real thickness the tape is a node, and the heat between two stacked sections has
    // to cross it: every turn resting on the tape gets an edge from its bottom face to the tape's
    // top face, every turn under it from its top face to the tape's bottom face, and the two
    // sections have no direct turn-to-turn edge across it. Only layers to the left and right of a
    // turn used to be candidates, so the tape had no turn edge at all.
    settings.reset();
    auto fixture = GENERATE(kLitzE96, kLitzAndStripEr35);
    INFO(fixture);
    auto masJson = load_abt1540_fixture(fixture);
    give_tapes_across_window_a_thickness(masJson, 50e-6, 10e-9);
    auto network = build_fixture_network(masJson, 1.0, {1.0, 1.0});

    auto tapes = tapes_across_window(network);
    REQUIRE(tapes.size() == 4);
    auto turns = network.of_part(ThermalNodePartType::TURN);
    for (size_t tape : tapes) {
        INFO("tape " << network.nodes[tape].name);
        std::vector<size_t> above, below;
        for (size_t turn : turns) {
            bool overTheTape = network.right(turn) > network.left(tape) && network.left(turn) < network.right(tape);
            if (!overTheTape) continue;
            if (std::abs(network.bottom(turn) - network.top(tape)) < kTouching) above.push_back(turn);
            if (std::abs(network.top(turn) - network.bottom(tape)) < kTouching) below.push_back(turn);
        }
        REQUIRE(!above.empty());
        REQUIRE(!below.empty());
        for (size_t turn : above) {
            INFO("turn on the tape: " << network.nodes[turn].name);
            REQUIRE(network.has_conduction_edge(turn, tape));
            CHECK(network.face_towards(turn, tape) == ThermalNodeFace::TANGENTIAL_RIGHT);  // turn's bottom
            CHECK(network.face_towards(tape, turn) == ThermalNodeFace::TANGENTIAL_LEFT);   // tape's top
        }
        for (size_t turn : below) {
            INFO("turn under the tape: " << network.nodes[turn].name);
            REQUIRE(network.has_conduction_edge(turn, tape));
            CHECK(network.face_towards(turn, tape) == ThermalNodeFace::TANGENTIAL_LEFT);   // turn's top
            CHECK(network.face_towards(tape, turn) == ThermalNodeFace::TANGENTIAL_RIGHT);  // tape's bottom
        }
        for (size_t upper : above) {
            for (size_t lower : below) {
                INFO(network.nodes[upper].name << " | " << network.nodes[lower].name);
                CHECK_FALSE(network.has_conduction_edge(upper, lower));
            }
        }
        // The tape's top and bottom lie against turns (or the still air of the window between
        // spread turns): they only conduct. Whatever sees the room is one of its thin edges.
        for (const auto& r : network.resistances) {
            if (r.nodeFromId != tape && r.nodeToId != tape) continue;
            if (r.type == HeatTransferType::CONDUCTION) continue;
            auto tapeFace = r.nodeFromId == tape ? r.quadrantFrom : r.quadrantTo;
            INFO("edge to " << network.nodes[r.nodeFromId == tape ? r.nodeToId : r.nodeFromId].name);
            CHECK((tapeFace == ThermalNodeFace::RADIAL_INNER || tapeFace == ThermalNodeFace::RADIAL_OUTER));
        }
        // The tape's thin edges are not faces a turn bears on.
        for (size_t turn : turns) {
            if (!network.has_conduction_edge(turn, tape)) continue;
            auto tapeFace = network.face_towards(tape, turn);
            INFO("turn " << network.nodes[turn].name);
            CHECK((tapeFace == ThermalNodeFace::TANGENTIAL_LEFT || tapeFace == ThermalNodeFace::TANGENTIAL_RIGHT));
        }
    }
    CHECK(network.error == "");
}

TEST_CASE("ThermalNetwork ABT 1540: a strip conducts to every turn lying on it", "[temperature][abt-1540][smoke-test]") {
    // D4. ER 35/20/11: each secondary section is one 3.645 x 0.709 mm strip, with four 0.746 mm
    // litz bundles of the primary lying on it and four under it (zero-thickness tapes: they
    // touch). Each bundle bears on the strip with its bottom (top) face. The network measured
    // every turn as a circle of half its SMALLER side, so the strip ended 0.35 mm from its own
    // centre and only the bundles next to that centre were its neighbours; and it read the
    // facing faces from the centre offsets, which made a bundle standing off the strip's centre
    // a side neighbour.
    settings.reset();
    auto masJson = load_abt1540_fixture(kLitzAndStripEr35);
    auto network = build_fixture_network(masJson, 1.0, {1.0, 1.0});

    auto turns = network.of_part(ThermalNodePartType::TURN);
    size_t strips = 0;
    size_t bundlesOnStrips = 0;
    for (size_t strip : turns) {
        if (network.nodes[strip].crossSectionalShape != TurnCrossSectionalShape::RECTANGULAR) continue;
        ++strips;
        for (size_t bundle : turns) {
            if (network.nodes[bundle].crossSectionalShape != TurnCrossSectionalShape::ROUND) continue;
            double x = network.nodes[bundle].physicalCoordinates[0];
            if (x < network.left(strip) || x > network.right(strip)) continue;
            bool onTop = std::abs(network.bottom(bundle) - network.top(strip)) < kTouching;
            bool under = std::abs(network.top(bundle) - network.bottom(strip)) < kTouching;
            if (!onTop && !under) continue;
            ++bundlesOnStrips;
            INFO(network.nodes[bundle].name << " on " << network.nodes[strip].name);
            REQUIRE(network.has_conduction_edge(bundle, strip));
            CHECK(network.face_towards(strip, bundle) == (onTop ? ThermalNodeFace::TANGENTIAL_LEFT : ThermalNodeFace::TANGENTIAL_RIGHT));
            CHECK(network.face_towards(bundle, strip) == (onTop ? ThermalNodeFace::TANGENTIAL_RIGHT : ThermalNodeFace::TANGENTIAL_LEFT));
        }
    }
    CHECK(strips == 2);
    CHECK(bundlesOnStrips == 16);  // 2 strips x (4 bundles on top + 4 under)
    CHECK(network.error == "");
}

TEST_CASE("ThermalNetwork ABT 1540: stacked sections reach temperatures that follow the losses", "[temperature][abt-1540][smoke-test]") {
    // The two valid layouts, as stored (zero-thickness tapes: the sections touch) and with 50 um
    // tapes, at 40 C ambient with stated losses instead of the loss models' (the stored operating
    // points give kilowatts of copper loss; that is not this network's business).
    //
    // What must hold whatever the model's detail:
    //  - the network solves: every node has a path to ambient;
    //  - everything is above ambient;
    //  - the core is no hotter than the whole loss leaving through its bounding box alone at
    //    5 W/m2K, the bottom of still-air convection plus radiation (the real surface is larger:
    //    winding overhang, window faces), and no cooler than the whole loss leaving through three
    //    times that box at 25 W/m2K, more than convection and radiation reach at these
    //    temperatures;
    //  - the hot spot is no further above ambient than twice that upper bound: the way from a
    //    turn to the core's surface is taken as no worse than the whole exterior film. Loose on
    //    purpose; it catches a winding left hanging on one thin path, not a 20 % error;
    //  - the winding that dissipates is the hotter one, both ways round;
    //  - doubling every loss heats every turn further;
    //  - a 50 um film between the sections moves the hot spot by less than 5 % of its rise: per
    //    turn it is a fraction of a K/W in a contact of 1.4 to 5.6 K/W.
    // Seen at 1 W core + 1 W primary + 1 W secondary. E 96/42/26 (five stacks, box 616 cm2):
    // turns 45.7 to 49.9 C, core 45.1 C; core bounds 40.6 to 49.7 C, hot spot bound 59.5 C.
    // ER 35/20/11: turns 81.2 to 116.7 C, core 80.6 C, its two strips the hot spot: each carries
    // 0.5 W and sheds it only through the eight bundles lying on it.
    settings.reset();
    auto fixture = GENERATE(kLitzE96, kLitzAndStripEr35);
    INFO(fixture);
    const double ambient = 40.0;
    auto stored = load_abt1540_fixture(fixture);
    auto taped = stored;
    give_tapes_across_window_a_thickness(taped, 50e-6, 10e-9);

    auto magnetic = autocomplete_like_the_web(stored).get_magnetic();
    auto core = magnetic.get_core().get_processed_description().value();
    const double boxArea = 2 * (core.get_width() * core.get_height() + core.get_width() * core.get_depth() +
                                core.get_height() * core.get_depth());
    INFO("core bounding box area " << boxArea * 1e4 << " cm2");
    REQUIRE(boxArea > 0);

    auto even = build_fixture_network(stored, 1.0, {1.0, 1.0});
    auto primaryOnly = build_fixture_network(stored, 1.0, {2.0, 0.0});
    auto secondaryOnly = build_fixture_network(stored, 1.0, {0.0, 2.0});
    auto doubled = build_fixture_network(stored, 2.0, {2.0, 2.0});
    auto evenTaped = build_fixture_network(taped, 1.0, {1.0, 1.0});
    for (const BuiltNetwork* network : {&even, &primaryOnly, &secondaryOnly, &doubled, &evenTaped}) {
        REQUIRE(network->error == "");
        REQUIRE(network->result.converged);
        for (const auto& [name, temperature] : network->result.nodeTemperatures) {
            INFO(name);
            REQUIRE(std::isfinite(temperature));
            CHECK(temperature >= ambient);
        }
    }

    const double totalLosses = 3.0;
    const double hottestCore = ambient + totalLosses / (5.0 * boxArea);
    const double coolestCore = ambient + totalLosses / (25.0 * 3.0 * boxArea);
    const double hottestSpot = ambient + 2 * (hottestCore - ambient);
    for (const BuiltNetwork* network : {&even, &primaryOnly, &secondaryOnly, &evenTaped}) {
        auto [lowest, highest] = network->turn_temperature_range();
        INFO("turns " << lowest << " to " << highest << " C, core " << network->result.averageCoreTemperature
             << " C; core bounds " << coolestCore << " to " << hottestCore << " C, hot spot bound " << hottestSpot << " C");
        CHECK(lowest > ambient);
        CHECK(highest < hottestSpot);
        CHECK(network->result.maximumTemperature < hottestSpot);
        CHECK(network->result.averageCoreTemperature > coolestCore);
        CHECK(network->result.averageCoreTemperature < hottestCore);
    }

    CHECK(primaryOnly.mean_turn_temperature(0) > primaryOnly.mean_turn_temperature(1));
    CHECK(secondaryOnly.mean_turn_temperature(1) > secondaryOnly.mean_turn_temperature(0));
    CHECK(primaryOnly.mean_turn_temperature(0) > secondaryOnly.mean_turn_temperature(0));
    CHECK(secondaryOnly.mean_turn_temperature(1) > primaryOnly.mean_turn_temperature(1));

    for (size_t turn : even.of_part(ThermalNodePartType::TURN)) {
        const std::string& name = even.nodes[turn].name;
        INFO(name);
        CHECK(doubled.result.nodeTemperatures.at(name) > even.result.nodeTemperatures.at(name));
    }

    const double rise = even.result.maximumTemperature - ambient;
    INFO("hot spot " << even.result.maximumTemperature << " C as stored, " << evenTaped.result.maximumTemperature << " C with 50 um tapes");
    CHECK(std::abs(evenTaped.result.maximumTemperature - even.result.maximumTemperature) < 0.05 * rise);
}

TEST_CASE("ThermalNetwork ABT 1540: stored stacked-section designs simulate as the web loads them", "[temperature][abt-1540]") {
    // mas_autocomplete + MagneticSimulator::simulate, the web's path for a loaded file, on the
    // two valid layouts that threw "disconnected node(s)" and on the PQ 35/40 file, which
    // simulated before the fix and must keep doing so.
    settings.reset();
    auto fixture = GENERATE(kLitzE96, kLitzAndStripEr35, kLitzAndStripPq3540);
    INFO(fixture);
    auto autocompleted = autocomplete_like_the_web(load_abt1540_fixture(fixture));
    OpenMagnetics::MagneticSimulator magneticSimulator;
    OpenMagnetics::Mas simulated;
    REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(autocompleted.get_inputs(), autocompleted.get_magnetic()));
    REQUIRE(simulated.get_outputs().size() == autocompleted.get_inputs().get_operating_points().size());
    for (size_t index = 0; index < simulated.get_outputs().size(); ++index) {
        const auto& output = simulated.get_outputs()[index];
        REQUIRE(output.get_temperature());
        double ambient = autocompleted.get_inputs().get_operating_points()[index].get_conditions().get_ambient_temperature();
        double maximumTemperature = output.get_temperature()->get_maximum_temperature();
        INFO("operating point " << index << ": ambient " << ambient << " C, maximum " << maximumTemperature << " C");
        CHECK(std::isfinite(maximumTemperature));
        CHECK(maximumTemperature > ambient);
    }
    settings.reset();
}

namespace {

// The web builder's re-wind of a loaded design (WebLibMKF wind binding): same windings, bobbin
// and section pattern, wound again even when the result does not fit the window.
OpenMagnetics::Magnetic rewind_like_the_builder(const OpenMagnetics::Mas& autocompleted, const std::vector<double>& proportions,
                                                const std::vector<size_t>& pattern, double tapeThickness) {
    json coilJson;
    to_json(coilJson, autocompleted.get_magnetic().get_coil());
    settings.set_coil_wind_even_if_not_fit(true);
    settings.set_coil_delimit_and_compact(true);
    settings.set_coil_include_additional_coordinates(true);
    auto windings = std::vector<OpenMagnetics::Winding>(coilJson["functionalDescription"]);
    OpenMagnetics::Coil coil;
    coil.set_bobbin_from_json(coilJson["bobbin"]);
    coil.set_functional_description(windings);
    coil.preload_margins(std::vector<std::vector<double>>(pattern.size(), std::vector<double>{0, 0}));
    coil.set_core_columns(autocompleted.get_magnetic().get_core().get_processed_description()->get_columns());
    coil.wind_by_sections(proportions, pattern, 1);
    for (auto& section : coil.get_sections_description_conduction()) {
        coil.set_layers_orientation(WindingOrientation::OVERLAPPING, section.get_name());
        coil.set_turns_alignment(CoilAlignment::SPREAD, section.get_name());
    }
    coil.set_interlayer_insulation(0.0, std::nullopt, std::nullopt, false);
    coil.set_intersection_insulation(tapeThickness, 1, std::nullopt, std::nullopt, false);
    coil.wind(proportions, pattern, 1);
    if (!coil.get_turns_description()) {
        throw std::runtime_error("the re-wind produced no turns");
    }
    coil.delimit_and_compact();
    auto magnetic = autocompleted.get_magnetic();
    magnetic.set_coil(coil);
    return magnetic;
}

} // namespace

TEST_CASE("ThermalNetwork ABT 1540: layouts drawn with copper inside copper still get a connected network", "[temperature][abt-1540][smoke-test]") {
    // Three layouts whose turns are drawn into each other, as the web shows them:
    //  - PQ 35/40 as stored: litz bundles 0.17 mm into each other, strips 0.07 mm. It solved
    //    before ABT #1540 and must keep solving;
    //  - the same design re-wound by the builder: sections taller than the window, turns up to
    //    0.69 mm into each other; "8 disconnected nodes" before;
    //  - PQ 65/60 of ABT #1487, wound on load: 0.87 mm; "2 disconnected nodes" before.
    // The network takes two outlines drawn into each other as touching (it always has) and says
    // nothing about the overlap; these designs get a temperature like any other. Bounds as loose
    // as in the valid layouts' test: above ambient, hot spot under twice the rise of the whole
    // loss leaving through the core's bounding box at 5 W/m2K, the lossy winding the hotter one.
    settings.reset();
    auto layout = GENERATE(std::string("PQ 35/40 stored"), std::string("PQ 35/40 re-wound"), std::string("PQ 65/60"));
    INFO(layout);
    OpenMagnetics::Magnetic magnetic;
    if (layout == "PQ 65/60") {
        magnetic = autocomplete_like_the_web(load_abt1540_fixture("abt1487_pq65_sp_contiguous.json")).get_magnetic();
    }
    else {
        auto autocompleted = autocomplete_like_the_web(load_abt1540_fixture(kLitzAndStripPq3540));
        magnetic = autocompleted.get_magnetic();
        if (layout == "PQ 35/40 re-wound") {
            magnetic = rewind_like_the_builder(autocompleted, {0.6636, 0.3364}, {0, 1, 0, 1, 0, 1, 0}, 25e-6);
            settings.reset();
        }
    }
    const double ambient = 40.0;
    auto core = magnetic.get_core().get_processed_description().value();
    const double boxArea = 2 * (core.get_width() * core.get_height() + core.get_width() * core.get_depth() +
                                core.get_height() * core.get_depth());
    const double hottestSpot = ambient + 2 * 3.0 / (5.0 * boxArea);
    // Winding 0 is the primary in the PQ 35/40 file and in the PQ 65/60 file.
    auto even = build_network(magnetic, ambient, 1.0, {1.0, 1.0});
    auto primaryOnly = build_network(magnetic, ambient, 1.0, {2.0, 0.0});
    auto secondaryOnly = build_network(magnetic, ambient, 1.0, {0.0, 2.0});
    for (const BuiltNetwork* network : {&even, &primaryOnly, &secondaryOnly}) {
        REQUIRE(network->error == "");
        REQUIRE(network->result.converged);
        for (const auto& [name, temperature] : network->result.nodeTemperatures) {
            INFO(name);
            REQUIRE(std::isfinite(temperature));
            CHECK(temperature >= ambient);
        }
        auto [lowest, highest] = network->turn_temperature_range();
        INFO("turns " << lowest << " to " << highest << " C, core " << network->result.averageCoreTemperature
             << " C, hot spot bound " << hottestSpot << " C");
        CHECK(lowest > ambient);
        CHECK(highest < hottestSpot);
    }
    CHECK(primaryOnly.mean_turn_temperature(0) > primaryOnly.mean_turn_temperature(1));
    CHECK(secondaryOnly.mean_turn_temperature(1) > secondaryOnly.mean_turn_temperature(0));
}
