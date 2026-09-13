// ABT #1175 (WP6): multi-chamber (split) bobbins - MAS-RFC 0014 part A.
//
// Vendor facts used as references (datasheets, not MKF output):
//   TDK RM8 B65812N1008D002: 2 sections, centre flange 0.56 mm, total winding width 5.5 mm.
//   TDK E20/10/6 B66206K1106T002: 2 sections of 6.1 mm, partition 0.65 mm.
//   Miles Platts ETD34 FD9645: 5 chambers of 2.2 mm, 4 partitions of 0.6 mm, 14.8 mm span
//     (5 x 2.2 + 4 x 0.6 = 13.4 mm between the end flanges; the 14.8 mm span adds the two 0.7 mm
//     end flanges).
//   Miles Platts ETD49 5-chamber: 5 equal chambers of 5.6 mm in 32.7 mm (walls (32.7 - 28) / 4).
//   Norwe EE 20/2k (09991-106), from the MAS catalogue: 12.5 mm between flanges, 0.6 mm wall.
#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/LeakageInductance.h"
#include "processors/Inputs.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <numbers>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

const std::string twoChamberE20 = "Bobbin EE 20 horizontal longer-creepage 8-pin 3.81mm 2-chamber (Norwe 09991-106)";
const std::string twoChamberEtd49 = "Bobbin ETD 49 horizontal lr 2-chamber 20-pin (Norwe 90565-186)";
const std::string untranscribedEtd19 = "Bobbin ETD 19 horizontal lr 2-chamber 10-pin (Norwe N0002-186)";

json functional_description_of(const std::string& bobbinName) {
    auto bobbin = find_bobbin_by_name(bobbinName);
    json functionalDescription;
    to_json(functionalDescription, bobbin.get_functional_description().value());
    return functionalDescription;
}

// A catalogue row's functional description with its chamber labels replaced.
json with_chambers(json functionalDescription, int64_t numberChambers, const std::vector<double>& chamberWidths,
                   const std::vector<double>& wallThicknesses) {
    auto& dimensions = functionalDescription["dimensions"];
    for (auto it = dimensions.begin(); it != dimensions.end();) {
        const std::string key = it.key();
        bool chamberLabel = key.size() > 1 && (key[0] == 'c' || key[0] == 'w') &&
                            std::all_of(key.begin() + 1, key.end(), ::isdigit);
        it = chamberLabel ? dimensions.erase(it) : std::next(it);
    }
    for (size_t index = 0; index < chamberWidths.size(); ++index) {
        dimensions["c" + std::to_string(index + 1)] = json{{"nominal", chamberWidths[index]}};
    }
    for (size_t index = 0; index < wallThicknesses.size(); ++index) {
        dimensions["w" + std::to_string(index + 1)] = json{{"nominal", wallThicknesses[index]}};
    }
    functionalDescription["numberChambers"] = numberChambers;
    return functionalDescription;
}

OpenMagnetics::Bobbin bobbin_from(json functionalDescription, const std::string& name) {
    json bobbinJson;
    bobbinJson["name"] = name;
    bobbinJson["functionalDescription"] = functionalDescription;
    return OpenMagnetics::Bobbin(bobbinJson);
}

void check_stack(OpenMagnetics::Bobbin bobbin, const std::vector<double>& chamberWidths, const std::vector<double>& wallThicknesses) {
    auto processed = bobbin.get_processed_description().value();
    auto windows = processed.get_winding_windows();
    REQUIRE(windows.size() == chamberWidths.size());
    REQUIRE(processed.get_dividers());
    auto dividers = processed.get_dividers().value();
    REQUIRE(dividers.size() == wallThicknesses.size());
    double stack = 0;
    for (auto width : chamberWidths) stack += width;
    for (auto thickness : wallThicknesses) stack += thickness;
    double top = stack / 2;
    for (size_t index = 0; index < windows.size(); ++index) {
        INFO("chamber " << index);
        CHECK(windows[index].get_column() == windows[0].get_column());
        CHECK_THAT(windows[index].get_coordinates().value()[0], WithinAbs(windows[0].get_coordinates().value()[0], 1e-12));
        CHECK_THAT(windows[index].get_width().value(), WithinAbs(windows[0].get_width().value(), 1e-12));
        CHECK_THAT(windows[index].get_height().value(), WithinAbs(chamberWidths[index], 1e-12));
        CHECK_THAT(windows[index].get_coordinates().value()[1], WithinAbs(top - chamberWidths[index] / 2, 1e-12));
        top -= chamberWidths[index];
        if (index < dividers.size()) {
            CHECK_THAT(dividers[index].get_thickness(), WithinAbs(wallThicknesses[index], 1e-12));
            CHECK_THAT(dividers[index].get_coordinates()[1], WithinAbs(top - wallThicknesses[index] / 2, 1e-12));
            CHECK_FALSE(dividers[index].get_height());
            top -= wallThicknesses[index];
        }
    }
}

json flyback_coil_json(const json& bobbin, int64_t primaryTurns, const std::string& primaryWire,
                       int64_t secondaryTurns, const std::string& secondaryWire, bool secondaryFirst = false) {
    json coilJson;
    coilJson["bobbin"] = bobbin;
    json primary = {{"name", "Primary"}, {"numberTurns", primaryTurns}, {"numberParallels", 1},
                    {"isolationSide", "primary"}, {"wire", primaryWire}};
    json secondary = {{"name", "Secondary"}, {"numberTurns", secondaryTurns}, {"numberParallels", 1},
                      {"isolationSide", "secondary"}, {"wire", secondaryWire}};
    coilJson["functionalDescription"] = secondaryFirst ? json::array({secondary, primary}) : json::array({primary, secondary});
    return coilJson;
}

}  // namespace

TEST_CASE("A two-chamber catalogue bobbin splits into two chambers sharing one column (ABT #1175)",
          "[constructive-model][bobbin][chambers]") {
    SECTION("E 20/10/6, Norwe 09991-106: 12.5 mm between flanges, one 0.6 mm wall") {
        auto bobbin = find_bobbin_by_name(twoChamberE20);
        REQUIRE(bobbin.get_processed_description());
        check_stack(bobbin, {0.00595, 0.00595}, {0.0006});
        auto windows = bobbin.get_processed_description()->get_winding_windows();
        // The two chambers hold the winding length between the flanges minus the wall.
        CHECK_THAT(windows[0].get_height().value() + windows[1].get_height().value(), WithinAbs(0.0125 - 0.0006, 1e-12));
        CHECK(bobbin.get_number_chambers() == 2);
        CHECK(bobbin.are_windows_chambers_of_same_column(0, 1));
        CHECK(bobbin.get_dividers_between_windows(0, 1) == std::vector<size_t>{0});
        CHECK(bobbin.get_dividers_between_windows(1, 0) == std::vector<size_t>{0});
    }

    SECTION("ETD 49, Norwe 90565-186: the wall bounding a chamber is the thinner of flange and divider") {
        auto bobbin = find_bobbin_by_name(twoChamberEtd49);
        check_stack(bobbin, {0.0159, 0.0159}, {0.001});
        // Flange (h1 - h2) / 2 = (35.2 - 32.2) / 2 = 1.5 mm; the centre wall is 1 mm.
        CHECK_THAT(bobbin.get_processed_description()->get_wall_thickness(), WithinAbs(0.0015, 1e-9));
        CHECK_THAT(bobbin.get_column_and_wall_thickness(0).second, WithinAbs(0.001, 1e-12));
        CHECK_THAT(bobbin.get_column_and_wall_thickness(1).second, WithinAbs(0.001, 1e-12));
        CHECK_THROWS_AS(bobbin.get_column_and_wall_thickness(2), InvalidInputException);
    }

    SECTION("a plain former keeps one window, no dividers and its flange as its wall") {
        auto bobbin = find_bobbin_by_name("Bobbin ETD 49");
        CHECK(bobbin.get_processed_description()->get_winding_windows().size() == 1);
        CHECK_FALSE(bobbin.get_processed_description()->get_dividers());
        CHECK(bobbin.get_number_chambers() == 1);
        CHECK_THAT(bobbin.get_column_and_wall_thickness(0).second, WithinAbs(0.0015, 1e-9));
        CHECK_THROWS_AS(bobbin.get_column_and_wall_thickness(1), InvalidInputException);
    }
}

TEST_CASE("Chamber stacks built from vendor drawings (ABT #1175)", "[constructive-model][bobbin][chambers]") {
    SECTION("TDK E20/10/6 B66206K1106T002: 2 x 6.1 mm, 0.65 mm partition") {
        auto bobbin = bobbin_from(with_chambers(functional_description_of(twoChamberE20), 2, {0.0061, 0.0061}, {0.00065}), "TDK E20 2k");
        check_stack(bobbin, {0.0061, 0.0061}, {0.00065});
        CHECK_THAT(bobbin.get_processed_description()->get_dividers().value()[0].get_coordinates()[1], WithinAbs(0, 1e-15));
    }
    SECTION("TDK RM8 B65812N1008D002: 5.5 mm winding width, 0.56 mm centre flange") {
        double chamber = (0.0055 - 0.00056) / 2;
        auto bobbin = bobbin_from(with_chambers(functional_description_of("Bobbin RM 8"), 2, {chamber, chamber}, {0.00056}), "TDK RM8 2k");
        check_stack(bobbin, {chamber, chamber}, {0.00056});
    }
    SECTION("Miles Platts ETD34 FD9645: 5 x 2.2 mm chambers, 4 x 0.6 mm partitions") {
        std::vector<double> chambers(5, 0.0022);
        std::vector<double> walls(4, 0.0006);
        auto bobbin = bobbin_from(with_chambers(functional_description_of("Bobbin ETD 34"), 5, chambers, walls), "Miles Platts ETD34 5k");
        check_stack(bobbin, chambers, walls);
        // Plus the two 0.7 mm end flanges, the stack is the catalogue's 14.8 mm span.
        CHECK_THAT(5 * 0.0022 + 4 * 0.0006 + 2 * 0.0007, WithinAbs(0.0148, 1e-12));
        CHECK(bobbin.get_dividers_between_windows(0, 4) == std::vector<size_t>{0, 1, 2, 3});
        CHECK(bobbin.get_dividers_between_windows(3, 1) == std::vector<size_t>{1, 2});
    }
    SECTION("Miles Platts ETD49 5-chamber: 5 x 5.6 mm in 32.7 mm") {
        std::vector<double> chambers(5, 0.0056);
        std::vector<double> walls(4, (0.0327 - 5 * 0.0056) / 4);
        auto bobbin = bobbin_from(with_chambers(functional_description_of("Bobbin ETD 49"), 5, chambers, walls), "Miles Platts ETD49 5k");
        check_stack(bobbin, chambers, walls);
        CHECK_THAT(bobbin.get_processed_description()->get_winding_windows()[2].get_coordinates().value()[1], WithinAbs(0, 1e-15));
    }
}

TEST_CASE("A chamber count without the walls is never split by invention (ABT #1175)", "[constructive-model][bobbin][chambers]") {
    SECTION("the catalogue keeps an untranscribed chamber former as functional data") {
        auto names = get_bobbin_names();
        REQUIRE(std::find(names.begin(), names.end(), untranscribedEtd19) != names.end());
        auto bobbin = find_bobbin_by_name(untranscribedEtd19);
        REQUIRE(bobbin.get_functional_description());
        CHECK(bobbin.get_functional_description()->get_number_chambers() == 2);
        CHECK_FALSE(bobbin.get_processed_description());
        CHECK_THROWS_WITH(bobbin.process_data(), ContainsSubstring("carry no 'c1'"));
    }
    SECTION("a design that uses it gets the refusal") {
        auto coilJson = flyback_coil_json(untranscribedEtd19, 20, "Round 0.25 - Grade 1", 5, "Round 0.5 - Grade 1");
        OpenMagnetics::Coil coil(coilJson, false);
        CHECK_THROWS_WITH(coil.resolve_bobbin(), ContainsSubstring("carry no 'c1'"));
    }
    auto base = functional_description_of(twoChamberE20);
    SECTION("a missing wall label throws") {
        auto description = with_chambers(base, 2, {0.00595, 0.00595}, {});
        CHECK_THROWS_WITH(bobbin_from(description, "no wall"), ContainsSubstring("no 'w1'"));
    }
    SECTION("a label beyond the declared count throws") {
        auto description = with_chambers(base, 2, {0.004, 0.004, 0.004}, {0.0003});
        CHECK_THROWS_WITH(bobbin_from(description, "extra chamber"), ContainsSubstring("'c3'"));
    }
    SECTION("a stack longer than the flanges throws") {
        auto description = with_chambers(base, 2, {0.0080, 0.0080}, {0.0006});
        CHECK_THROWS_WITH(bobbin_from(description, "too long"), ContainsSubstring("stack to"));
    }
    SECTION("a non-positive label throws") {
        auto description = with_chambers(base, 2, {0.00595, 0.00595}, {0});
        CHECK_THROWS_WITH(bobbin_from(description, "zero wall"), ContainsSubstring("must be positive"));
    }
}

TEST_CASE("A flyback wound into a two-chamber bobbin puts the primary in chamber 0 and the secondary in chamber 1 (ABT #1175)",
          "[constructive-model][coil][chambers]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    for (bool secondaryFirst : {false, true}) {
        INFO("secondary declared first: " << secondaryFirst);
        auto coilJson = flyback_coil_json(twoChamberE20, 40, "Round 0.25 - Grade 1", 10, "Round 0.5 - Grade 1", secondaryFirst);
        OpenMagnetics::Coil coil(coilJson, false);
        REQUIRE(coil.wind());
        auto bobbin = coil.resolve_bobbin();
        auto windows = bobbin.get_processed_description()->get_winding_windows();
        const auto range1 = coil.get_sections_description().value();
        for (const auto& section : range1) {
            if (section.get_type() != ElectricalType::CONDUCTION) {
                FAIL_CHECK("section " << section.get_name() << " is not a conduction section: no tape belongs between chambers");
                continue;
            }
            size_t windowIndex = coil.resolve_section_winding_window_index(section);
            const auto& window = windows[windowIndex];
            double top = window.get_coordinates().value()[1] + window.get_height().value() / 2;
            double bottom = window.get_coordinates().value()[1] - window.get_height().value() / 2;
            CHECK(section.get_coordinates()[1] + section.get_dimensions()[1] / 2 <= top + 1e-9);
            CHECK(section.get_coordinates()[1] - section.get_dimensions()[1] / 2 >= bottom - 1e-9);
            if (section.get_partial_windings()[0].get_winding() == "Primary") {
                CHECK(windowIndex == 0);
            }
            else {
                CHECK(windowIndex == 1);
            }
        }
        const auto range2 = coil.get_turns_description().value();
        for (const auto& turn : range2) {
            if (turn.get_winding() == "Primary") {
                CHECK(turn.get_coordinates()[1] > 0.0003);
            }
            else {
                CHECK(turn.get_coordinates()[1] < -0.0003);
            }
        }
    }
    settings.reset();
}

TEST_CASE("A windingWindow on the winding overrides the default chamber (ABT #1175)", "[constructive-model][coil][chambers]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto coilJson = flyback_coil_json(twoChamberE20, 40, "Round 0.25 - Grade 1", 10, "Round 0.5 - Grade 1");
    coilJson["functionalDescription"][0]["windingWindow"] = 1;
    coilJson["functionalDescription"][1]["windingWindow"] = 0;
    OpenMagnetics::Coil coil(coilJson, false);
    REQUIRE(coil.wind());
    const auto range3 = coil.get_turns_description().value();
    for (const auto& turn : range3) {
        CHECK((turn.get_winding() == "Primary" ? turn.get_coordinates()[1] < 0 : turn.get_coordinates()[1] > 0));
    }

    SECTION("more isolation sides than chambers is refused, not piled into one chamber") {
        auto threeSides = flyback_coil_json(twoChamberE20, 20, "Round 0.25 - Grade 1", 10, "Round 0.25 - Grade 1");
        threeSides["functionalDescription"].push_back(json{{"name", "Tertiary"}, {"numberTurns", 5}, {"numberParallels", 1},
                                                           {"isolationSide", "tertiary"}, {"wire", "Round 0.25 - Grade 1"}});
        OpenMagnetics::Coil crowded(threeSides, false);
        CHECK_THROWS_WITH(crowded.wind(), ContainsSubstring("only 2 chambers"));
    }
    settings.reset();
}

TEST_CASE("Chambers on one column get no second crossings; lateral windows keep theirs (ABT #1175)",
          "[constructive-model][coil][chambers][multi-column]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_include_additional_coordinates(true);
    auto coilJson = flyback_coil_json(twoChamberE20, 40, "Round 0.25 - Grade 1", 10, "Round 0.5 - Grade 1");
    OpenMagnetics::Coil coil(coilJson, false);
    REQUIRE(coil.wind());
    auto turns = coil.get_turns_description().value();
    REQUIRE(turns.size() == 50);
    for (const auto& turn : turns) {
        CHECK(turn.get_coordinates()[0] > 0);
        CHECK_FALSE(turn.get_additional_coordinates());
    }
    settings.reset();
}

TEST_CASE("A divider thicker than the required DTI is the barrier between two chambers (ABT #1175)",
          "[constructive-model][insulation][chambers]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    // IEC 62368-1, reinforced insulation, 230 V mains: a real solid-insulation requirement.
    DimensionWithTolerance altitude;
    altitude.set_maximum(2000);
    DimensionWithTolerance mainSupplyVoltage;
    mainSupplyVoltage.set_nominal(230);
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(100000, 0.001, 25, WaveformLabel::SINUSOIDAL, 400, 0.5, 0, {4.0});
    auto designRequirements = inputs.get_design_requirements();
    designRequirements.set_insulation(OpenMagneticsTesting::get_quick_insulation_requirements(
        altitude, Cti::GROUP_I, IsolationClass::REINFORCED, mainSupplyVoltage, OvervoltageCategory::II, PollutionDegree::PD2,
        {InsulationStandards::IEC_623681}));
    inputs.set_design_requirements(designRequirements);
    double requiredDti = InsulationCoordinator().calculate_distance_through_insulation(inputs);
    INFO("required DTI " << requiredDti);
    REQUIRE(requiredDti > 0);
    REQUIRE(requiredDti <= 0.0006);

    auto coilJson = flyback_coil_json(twoChamberE20, 40, "Round 0.25 - Grade 1", 10, "Round 0.5 - Grade 1");
    OpenMagnetics::Coil coil(coilJson, false);
    coil.set_inputs(inputs);
    REQUIRE(coil.wind());

    for (auto pair : {std::pair<size_t, size_t>{0, 1}, std::pair<size_t, size_t>{1, 0}}) {
        auto interface = coil.get_coil_section_interface(pair.first, pair.second);
        REQUIRE(interface);
        REQUIRE(interface->get_barrier_divider_index());
        CHECK(interface->get_barrier_divider_index().value() == 0);
        CHECK(interface->get_total_margin_tape_distance() == 0);
        CHECK(interface->get_number_layers_insulation() == 0);
        CHECK(interface->get_solid_insulation_thickness() == 0);
    }
    const auto range4 = coil.get_sections_description().value();
    for (const auto& section : range4) {
        INFO(section.get_name());
        CHECK(section.get_type() == ElectricalType::CONDUCTION);
        auto margin = OpenMagnetics::Coil::resolve_margin(section);
        CHECK(margin[0] == 0);
        CHECK(margin[1] == 0);
    }

    SECTION("the same windings in ONE window of the same former get tape and margin") {
        auto bobbin = find_bobbin_by_name(twoChamberE20);
        auto processed = bobbin.get_processed_description().value();
        auto window = processed.get_winding_windows()[0];
        window.set_height(0.0125);
        window.set_coordinates(std::vector<double>{window.get_coordinates().value()[0], 0});
        window.set_area(window.get_width().value() * 0.0125);
        processed.set_winding_windows({window});
        processed.set_dividers(std::nullopt);
        json singleWindowBobbin;
        singleWindowBobbin["name"] = "E20 single window";
        singleWindowBobbin["processedDescription"] = processed;
        OpenMagnetics::Coil single(flyback_coil_json(singleWindowBobbin, 40, "Round 0.25 - Grade 1", 10, "Round 0.5 - Grade 1"), false);
        single.set_inputs(inputs);
        single.wind();
        bool anyInsulationSection = false;
        const auto range5 = single.get_sections_description().value();
        for (const auto& section : range5) {
            anyInsulationSection = anyInsulationSection || section.get_type() == ElectricalType::INSULATION;
        }
        CHECK(anyInsulationSection);
        auto interface = single.get_coil_section_interface(0, 1);
        REQUIRE(interface);
        CHECK_FALSE(interface->get_barrier_divider_index());
    }

    SECTION("a wall thinner than the DTI between two isolation sides is refused") {
        auto bobbin = find_bobbin_by_name(twoChamberE20);
        auto processed = bobbin.get_processed_description().value();
        auto dividers = processed.get_dividers().value();
        dividers[0].set_thickness(requiredDti / 4);
        processed.set_dividers(dividers);
        json thinWallBobbin;
        thinWallBobbin["name"] = "E20 thin wall";
        thinWallBobbin["processedDescription"] = processed;
        OpenMagnetics::Coil thin(flyback_coil_json(thinWallBobbin, 40, "Round 0.25 - Grade 1", 10, "Round 0.5 - Grade 1"), false);
        thin.set_inputs(inputs);
        CHECK_THROWS_WITH(thin.wind(), ContainsSubstring("no wall between them is a barrier"));
    }
    settings.reset();
}

TEST_CASE("A series junction across a wall reserves the divider's crossing slot (ABT #1175)", "[constructive-model][coil][chambers]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto bobbin = find_bobbin_by_name(twoChamberE20);
    auto processed = bobbin.get_processed_description().value();
    auto window = processed.get_winding_windows()[0];

    auto seriesCoil = [&](const CoreBobbinProcessedDescription& description) {
        json bobbinJson;
        bobbinJson["name"] = "E20 2k";
        bobbinJson["processedDescription"] = description;
        json coilJson = flyback_coil_json(bobbinJson, 30, "Round 0.25 - Grade 1", 30, "Round 0.25 - Grade 1");
        // Two halves of one primary: same side, one per chamber, joined by a wire junction "J" that
        // is not a pin of the bobbin.
        coilJson["functionalDescription"][1]["name"] = "Primary B";
        coilJson["functionalDescription"][1]["isolationSide"] = "primary";
        coilJson["functionalDescription"][0]["windingWindow"] = 0;
        coilJson["functionalDescription"][1]["windingWindow"] = 1;
        coilJson["functionalDescription"][0]["connections"] = json::array({json{{"pinName", "1"}, {"end", "start"}},
                                                                          json{{"pinName", "J"}, {"end", "finish"}}});
        coilJson["functionalDescription"][1]["connections"] = json::array({json{{"pinName", "J"}, {"end", "start"}},
                                                                          json{{"pinName", "2"}, {"end", "finish"}}});
        return OpenMagnetics::Coil(coilJson, false);
    };

    SECTION("through the slot") {
        auto slotted = processed;
        auto dividers = slotted.get_dividers().value();
        DividerCrossingSlot slot;
        slot.set_width(0.001);
        slot.set_depth(0.0015);
        slot.set_angle(90);
        dividers[0].set_crossing_slot(slot);
        slotted.set_dividers(dividers);
        auto coil = seriesCoil(slotted);
        REQUIRE(coil.wind());
        auto crossovers = coil.get_chamber_crossovers();
        REQUIRE(crossovers.size() == 1);
        const auto& crossover = crossovers[0];
        CHECK(crossover.fromWinding == "Primary");
        CHECK(crossover.toWinding == "Primary B");
        CHECK(crossover.dividerIndex == 0);
        CHECK(crossover.throughSlot);
        REQUIRE(crossover.slotAngle);
        CHECK(crossover.slotAngle.value() == 90);
        double wireOuterWidth = coil.get_wires()[0].get_maximum_outer_width();
        double flangeReach = window.get_coordinates().value()[0] + window.get_width().value() / 2;
        CHECK_THAT(crossover.coordinates[0], WithinAbs(flangeReach - 0.0015 + wireOuterWidth / 2, 1e-12));
        CHECK_THAT(crossover.coordinates[1], WithinAbs(0, 1e-12));
        CHECK_THAT(crossover.dimensions[1], WithinAbs(0.0006, 1e-12));
        REQUIRE(crossover.waypoints.size() == 2);
        CHECK_THAT(crossover.waypoints[0][1], WithinAbs(0.0003, 1e-12));
        CHECK_THAT(crossover.waypoints[1][1], WithinAbs(-0.0003, 1e-12));
        CHECK_THAT(crossover.routedLength, WithinAbs(0.0006, 1e-12));
    }

    SECTION("a full-height wall with no slot has no route") {
        auto coil = seriesCoil(processed);
        CHECK_THROWS_WITH(coil.wind(), ContainsSubstring("has no crossingSlot"));
    }

    SECTION("a slot narrower than the wire is refused") {
        auto slotted = processed;
        auto dividers = slotted.get_dividers().value();
        DividerCrossingSlot slot;
        slot.set_width(0.0001);
        slot.set_depth(0.0015);
        dividers[0].set_crossing_slot(slot);
        slotted.set_dividers(dividers);
        auto coil = seriesCoil(slotted);
        CHECK_THROWS_WITH(coil.wind(), ContainsSubstring("crossingSlot, which is"));
    }

    SECTION("a wall that stops short of the flanges is crossed over its rim") {
        auto shortWall = processed;
        auto dividers = shortWall.get_dividers().value();
        dividers[0].set_height(window.get_width().value() / 2);
        shortWall.set_dividers(dividers);
        auto coil = seriesCoil(shortWall);
        REQUIRE(coil.wind());
        auto crossovers = coil.get_chamber_crossovers();
        REQUIRE(crossovers.size() == 1);
        CHECK_FALSE(crossovers[0].throughSlot);
        double columnSurface = window.get_coordinates().value()[0] - window.get_width().value() / 2;
        double wireOuterWidth = coil.get_wires()[0].get_maximum_outer_width();
        CHECK_THAT(crossovers[0].coordinates[0], WithinAbs(columnSurface + window.get_width().value() / 2 + wireOuterWidth / 2, 1e-12));
    }

    SECTION("a junction on a pin of the bobbin is two leads, not a crossover") {
        auto pinned = processed;
        MAS::Pin pin;
        pin.set_name("J");
        pin.set_shape(PinShape::ROUND);
        pin.set_type(PinDescriptionType::THT);
        pin.set_dimensions({0.0006, 0.0006, 0.003});
        pin.set_coordinates(std::vector<double>{0, 0.008, -0.004});
        pinned.set_pins(std::vector<MAS::Pin>{pin});
        auto coil = seriesCoil(pinned);
        REQUIRE(coil.wind());
        CHECK(coil.get_chamber_crossovers().empty());
    }
    settings.reset();
}

// ETD 49 (Norwe 90565-186: 2 x 15.9 mm chambers, 1 mm wall, 8 mm deep), 144 + 144 turns of 0.8 mm so each
// winding fills ~86 % of the chamber depth - McLyman's 1D model assumes the sections span the build.
// Not the E 20 record: the E-family processor leaves columnDepth at 0, which shortens every turn of a
// catalogue E bobbin by 4 x depth (reported separately); a round ETD column has no such defect.
TEST_CASE("Leakage of a flyback rises when its windings go side by side into two chambers (ABT #1175)",
          "[physical-model][leakage-inductance][chambers]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const int64_t primaryTurns = 144;
    const int64_t secondaryTurns = 144;
    const std::string primaryWire = "Round 0.80 - Grade 1";
    const std::string secondaryWire = "Round 0.80 - Grade 1";
    const double frequency = 100000;
    auto core = OpenMagnetics::Core::create_quick_core("ETD 49/25/16", "3C97", OpenMagnetics::Core::create_ground_gapping(0.0002, 3));

    auto chamberedCoil = OpenMagnetics::Coil(flyback_coil_json(twoChamberEtd49, primaryTurns, primaryWire, secondaryTurns, secondaryWire), false);
    REQUIRE(chamberedCoil.wind());
    OpenMagnetics::Magnetic chambered;
    chambered.set_core(core);
    chambered.set_coil(chamberedCoil);
    double chamberedLeakage = LeakageInductance().calculate_leakage_inductance(chambered, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();

    auto processed = find_bobbin_by_name(twoChamberEtd49).get_processed_description().value();
    auto window = processed.get_winding_windows()[0];
    window.set_height(0.0328);
    window.set_coordinates(std::vector<double>{window.get_coordinates().value()[0], 0});
    window.set_area(window.get_width().value() * 0.0328);
    processed.set_winding_windows({window});
    processed.set_dividers(std::nullopt);
    json singleWindowBobbin;
    singleWindowBobbin["name"] = "ETD 49 single window";
    singleWindowBobbin["processedDescription"] = processed;
    auto singleCoil = OpenMagnetics::Coil(flyback_coil_json(singleWindowBobbin, primaryTurns, primaryWire, secondaryTurns, secondaryWire), false);
    REQUIRE(singleCoil.wind());
    OpenMagnetics::Magnetic single;
    single.set_core(core);
    single.set_coil(singleCoil);
    double singleLeakage = LeakageInductance().calculate_leakage_inductance(single, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();

    // McLyman, Transformer and Inductor Design Handbook, eq. 17-5, side-by-side sections:
    // L = 4 pi MLT N^2 / b (sum c + sum a / 3) 1e-9 H, lengths in cm; a = axial width of each
    // winding's copper, c = the axial gap between them, b = the radial build the leakage flux
    // crosses (the chamber depth), MLT the mean turn length, N the primary turns.
    auto turns = chamberedCoil.get_turns_description().value();
    double primaryTop = -1, primaryBottom = 1, secondaryTop = -1, secondaryBottom = 1, lengthSum = 0;
    for (const auto& turn : turns) {
        double y = turn.get_coordinates()[1];
        double half = turn.get_dimensions().value()[1] / 2;
        if (turn.get_winding() == "Primary") {
            primaryTop = std::max(primaryTop, y + half);
            primaryBottom = std::min(primaryBottom, y - half);
        }
        else {
            secondaryTop = std::max(secondaryTop, y + half);
            secondaryBottom = std::min(secondaryBottom, y - half);
        }
        lengthSum += turn.get_length();
    }
    double meanTurnLengthCm = lengthSum / turns.size() * 100;
    double primaryWidthCm = (primaryTop - primaryBottom) * 100;
    double secondaryWidthCm = (secondaryTop - secondaryBottom) * 100;
    double gapCm = (primaryBottom - secondaryTop) * 100;
    double buildCm = window.get_width().value() * 100;
    double mcLyman = 4 * std::numbers::pi * meanTurnLengthCm * primaryTurns * primaryTurns / buildCm *
                     (gapCm + (primaryWidthCm + secondaryWidthCm) / 3) * 1e-9;

    INFO("chambered " << chamberedLeakage << " H, single window " << singleLeakage << " H, McLyman 17-5 " << mcLyman
         << " H (MLT " << meanTurnLengthCm << " cm, a " << primaryWidthCm << " + " << secondaryWidthCm << " cm, c " << gapCm
         << " cm, b " << buildCm << " cm)");
    std::cout << "[ABT #1175 leakage] chambered " << chamberedLeakage << " H, single window " << singleLeakage
              << " H, McLyman 17-5 " << mcLyman << " H (MLT " << meanTurnLengthCm << " cm, a " << primaryWidthCm << " + "
              << secondaryWidthCm << " cm, c " << gapCm << " cm, b " << buildCm << " cm)" << std::endl;
    CHECK(chamberedLeakage > singleLeakage);
    CHECK_THAT(chamberedLeakage, WithinRel(mcLyman, 0.15));
    settings.reset();
}
