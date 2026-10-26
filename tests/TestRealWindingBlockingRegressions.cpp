// Real-winding connection blocking: designs MVB++ builds as FEM-grade CAD must keep getting
// blocking applied by magnetic_autocomplete. MVB++ refuses a real-winding build whose blocking was
// not applied (magnetic_autocomplete_safe), so a regression here turns its fixtures red.

#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "processors/ResultsAnalyzer.h"
#include "support/Logger.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "support/Settings.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <fstream>
#include <set>
#include <optional>
#include <source_location>

#include "TestingUtils.h"

using json = nlohmann::json;

namespace {

json load_test_data(const std::string& name) {
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("testData")
                           .append(name));
    REQUIRE(file.good());
    return json::parse(file);
}

}  // namespace

// ABT #1194: MVB++'s cm37 fixture (E16, Primary + Secondary, 22 turns x 2 parallels each,
// interleaved). Built exactly as MVB++'s magnetic_autocomplete_safe builds it: Core and Coil from
// the stored MAS (the coil NOT re-wound by its constructor), then magnetic_autocomplete with real
// winding geometry on and no fit relaxations.
TEST_CASE("Real winding: the interleaved N-filar E16 (MVB++ cm37) gets connection blocking",
          "[constructive-model][coil][real-winding][abt1194]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    const auto masJson = load_test_data("abt1194_interleaved_nfilar_e16.json");
    const auto& magneticJson = masJson.at("magnetic");

    settings.set_coil_use_real_winding_geometry(true);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
    REQUIRE(enriched.get_coil().is_real_winding_blocking_applied());
    settings.reset();
}

// ABT #1860: MVB++'s current_transformer_complete, a one-turn pass-through primary (one Round 6.0
// turn, 4.19 mm, subtending 46.5 deg of T 25.3/14.8/10) beside 100 turns of Round 0.5. Built the
// way magnetic_autocomplete_safe builds it, for both toroid mountings. The contiguous sector
// apportioning (ABT #1734) counted every ring as holding a turn at any angle, so the primary
// claimed 0 deg, got a 0.29 deg layer (fill factor 158), and magnetic_autocomplete refused the
// real wind before any connection blocking ran.
TEST_CASE("Real winding: a one-turn pass-through toroid primary (MVB++ current_transformer_complete) gets connection blocking",
          "[constructive-model][coil][real-winding][abt1860]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    const auto magneticJson = json::parse(R"({"core": {"name": "Toroid T 25.3/14.8/10 N30", "functionalDescription": {"type": "toroidal", "material": "N30", "shape": "T 25.3/14.8/10", "gapping": [], "numberStacks": 1}}, "coil": {"bobbin": "Basic", "functionalDescription": [{"name": "Primary", "numberTurns": 1, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 6.0 - Single Build"}, {"name": "Secondary", "numberTurns": 100, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 0.5 - Grade 1"}]}})");
    for (auto mounting : {MAS::OrientationEnum::VERTICAL, MAS::OrientationEnum::HORIZONTAL}) {
        settings.reset();
        settings.set_coil_use_real_winding_geometry(true);
        settings.set_toroid_mounting(mounting);
        OpenMagnetics::Core core(magneticJson.at("core"));
        OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::Magnetic enriched;
        REQUIRE_NOTHROW(enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{}));
        auto& enrichedCoil = enriched.get_mutable_coil();
        INFO(enrichedCoil.get_last_fit_failure());
        CHECK(enrichedCoil.is_real_winding_blocking_applied());
        CHECK(enrichedCoil.are_sections_and_layers_fitting());
    }
    settings.reset();
}

namespace {

// The parallels each conduction section of `winding` carries, in wound order.
std::vector<std::set<size_t>> parallels_per_section(const OpenMagnetics::Coil& coil, const std::string& winding) {
    std::vector<std::set<size_t>> parallelsPerSection;
    const auto sections = coil.get_sections_description().value();
    for (const auto& section : sections) {
        if (section.get_type() != MAS::ElectricalType::CONDUCTION ||
            section.get_partial_windings()[0].get_winding() != winding) {
            continue;
        }
        std::set<size_t> carried;
        const auto proportion = section.get_partial_windings()[0].get_parallels_proportion();
        for (size_t parallelIndex = 0; parallelIndex < proportion.size(); ++parallelIndex) {
            if (proportion[parallelIndex] > 0) {
                carried.insert(parallelIndex);
            }
        }
        parallelsPerSection.push_back(carried);
    }
    return parallelsPerSection;
}

// Primary 20 turns x 2 parallels + Secondary 20 x 1, interleaved twice on a PQ 40/40: the pinned
// bifilar layout (Test_Real_Geometry_Bifilar_Interleaved). Unstated, the heuristic picks
// consecutive turns (parallels == sections), one parallel per section, for the ideal and the real
// wind alike.
OpenMagnetics::Coil bifilar_interleaved_pq40() {
    return OpenMagneticsTesting::get_quick_coil({20, 20}, {2, 1}, "PQ 40/40", 2);
}

}  // namespace

// ABT #1487: the crossing stations are cross-sections, not turns, and must not change how a
// winding's parallels are shared out between its sections. 2 turns x 3 parallels in 3 interleaved
// sections is one parallel per section (fewer turns than sections: the heuristic has no other
// choice); with one station per layer the counter read 3 turns == 3 sections, chose "one turn per
// section" and redistribute_section_turns_for_blocking then spread every parallel over every
// section. Each parallel spanned every layer, was charged a station in each, and the raise grew
// without bound (the web's PQ 65/60 with an edge-wound 3 x 4 secondary threw "turn blocking and
// the per-layer crossing stations do not converge ... 76 layers").
TEST_CASE("Real winding keeps one parallel per section when there are fewer turns than sections (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    auto build = []() {
        return OpenMagneticsTesting::get_quick_coil({20, 2}, {1, 3}, "PQ 40/40", 3,
                                                    MAS::WindingOrientation::OVERLAPPING,
                                                    MAS::WindingOrientation::OVERLAPPING,
                                                    MAS::CoilAlignment::SPREAD, MAS::CoilAlignment::CENTERED,
                                                    {OpenMagnetics::find_wire_by_name("Round 0.5 - Grade 1"),
                                                     OpenMagnetics::find_wire_by_name("Round 1.00 - Grade 1")});
    };

    settings.reset();
    auto ideal = build();
    REQUIRE(ideal.wind());
    const auto secondary = ideal.get_functional_description()[1].get_name();
    const auto idealSplit = parallels_per_section(ideal, secondary);
    REQUIRE(idealSplit == std::vector<std::set<size_t>>{{0}, {1}, {2}});

    settings.set_coil_use_real_winding_geometry(true);
    auto real = build();
    REQUIRE_NOTHROW(real.wind());
    REQUIRE(real.is_real_winding_blocking_applied());
    CHECK(parallels_per_section(real, secondary) == idealSplit);
    settings.reset();
}

// ABT #1487 (Alf, 2026-09-29): "These are two different ways of winding that must BOTH be supported.
// The user can choose one or the other." A stated winding style is honoured by real winding in
// both directions, on the same design whose unstated layout is side by side.
TEST_CASE("Real winding honours a stated winding style both ways (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    const std::string primary = "winding 0";
    const std::vector<std::set<size_t>> sideBySide{{0, 1}, {0, 1}};
    const std::vector<std::set<size_t>> onePerSection{{0}, {1}};

    settings.reset();
    {
        auto ideal = bifilar_interleaved_pq40();
        REQUIRE(parallels_per_section(ideal, primary) == onePerSection);
    }
    settings.set_coil_use_real_winding_geometry(true);
    {
        // Unstated: the heuristic's choice, the same one the ideal wind made.
        auto unstated = bifilar_interleaved_pq40();
        REQUIRE(unstated.is_real_winding_blocking_applied());
        CHECK(parallels_per_section(unstated, primary) == onePerSection);
    }
    for (auto [style, expected] : {std::pair{MAS::WindingStyle::WIND_BY_CONSECUTIVE_TURNS, onePerSection},
                                   std::pair{MAS::WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS, sideBySide}}) {
        INFO("stated style " << (style == MAS::WindingStyle::WIND_BY_CONSECUTIVE_TURNS ? "turns" : "parallels"));
        auto stated = bifilar_interleaved_pq40();
        stated.preload_winding_style_overrides({{primary, style}});
        REQUIRE_NOTHROW(stated.wind(2));
        REQUIRE(stated.is_real_winding_blocking_applied());
        CHECK(parallels_per_section(stated, primary) == expected);
    }
    settings.reset();
}

// A stored design keeps the way its parallels were wound through magnetic_autocomplete's
// real-winding re-wind (the path MVB++ and the web take). A winding-studio choice is transient and
// never saved, so the stored sections' parallelsProportion is the only record of it: a design
// stored side by side must not come back one parallel per section because the heuristic, asked
// again, would pick that. Stored both ways from the SAME design, each through the re-wind.
TEST_CASE("magnetic_autocomplete keeps a stored design's parallel layout through the real-winding re-wind (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    const std::string primary = "winding 0";
    for (auto [style, expected] : {std::pair{MAS::WindingStyle::WIND_BY_CONSECUTIVE_TURNS,
                                             std::vector<std::set<size_t>>{{0}, {1}}},
                                   std::pair{MAS::WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS,
                                             std::vector<std::set<size_t>>{{0, 1}, {0, 1}}}}) {
        INFO("stored " << (style == MAS::WindingStyle::WIND_BY_CONSECUTIVE_TURNS ? "one parallel per section" : "side by side"));
        settings.reset();
        auto stored = bifilar_interleaved_pq40();   // ideal wind, then the user's choice, as stored
        stored.preload_winding_style_overrides({{primary, style}});
        REQUIRE(stored.wind(2));
        REQUIRE(parallels_per_section(stored, primary) == expected);
        json coilJson;
        to_json(coilJson, stored);

        settings.set_coil_use_real_winding_geometry(true);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(OpenMagneticsTesting::get_quick_core("PQ 40/40", OpenMagneticsTesting::get_ground_gap(0.001)));
        magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));   // the override did not survive
        REQUIRE(magnetic.get_coil().get_winding_style_overrides().empty());
        auto enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
        REQUIRE(enriched.get_coil().is_real_winding_blocking_applied());
        CHECK(parallels_per_section(enriched.get_coil(), primary) == expected);
    }
    settings.reset();
}

namespace {

// The PQ 65/60 field design of ABT #1487: seven CONTIGUOUS (axially stacked) sections S,P,S,P,S,P,S.
// Secondary 3 turns x 4 parallels of edge-wound 10.085 x 0.885 mm rectangular wire, one whole
// parallel per section; Primary 45 turns x 1 parallel of 1.865 mm litz, split in series 1/3 per
// section. Stored section heights: secondary 2.918 mm, primary 7.591 mm. `edit` changes the
// stored MAS before it goes through magnetic_autocomplete with real winding on.
OpenMagnetics::Magnetic pq65_through_real_winding(const std::function<void(json&)>& edit) {
    auto masJson = load_test_data("abt1487_pq65_sp_contiguous.json");
    edit(masJson);
    OpenMagnetics::Mas mas(masJson);
    OpenMagnetics::Settings::GetInstance().set_coil_use_real_winding_geometry(true);
    return OpenMagnetics::magnetic_autocomplete(OpenMagnetics::Magnetic(masJson.at("magnetic")), json{}, mas.get_inputs());
}

json& pq65_winding(json& masJson, const std::string& name) {
    for (auto& winding : masJson["magnetic"]["coil"]["functionalDescription"]) {
        if (winding["name"] == name) {
            return winding;
        }
    }
    throw std::runtime_error("no winding " + name);
}

void set_pq65_stored_heights(json& masJson, double secondaryHeight, double primaryHeight) {
    for (auto& section : masJson["magnetic"]["coil"]["sectionsDescription"]) {
        if (section["type"] != "conduction") {
            continue;
        }
        const std::string name = section["name"];
        section["dimensions"][1] = name.rfind("Secondary", 0) == 0 ? secondaryHeight : primaryHeight;
    }
}

struct SectionFrame {
    std::map<std::string, std::pair<double, double>> axialExtent;   // conduction section -> {low, high}
    std::map<std::string, std::string> sectionOfLayer;
    std::map<std::string, std::string> sectionOfTurn;
    double copperOuterFace = std::numeric_limits<double>::lowest();
};

SectionFrame section_frame(const OpenMagnetics::Coil& coil) {
    SectionFrame frame;
    const auto sectionsCopy = coil.get_sections_description().value();
    for (const auto& section : sectionsCopy) {
        if (section.get_type() != MAS::ElectricalType::CONDUCTION) {
            continue;
        }
        frame.axialExtent[section.get_name()] = {section.get_coordinates()[1] - section.get_dimensions()[1] / 2,
                                                 section.get_coordinates()[1] + section.get_dimensions()[1] / 2};
    }
    const auto layersCopy = coil.get_layers_description().value();
    for (const auto& layer : layersCopy) {
        if (layer.get_type() != MAS::ElectricalType::CONDUCTION) {
            continue;
        }
        frame.sectionOfLayer[layer.get_name()] = layer.get_section().value();
        frame.copperOuterFace = std::max(frame.copperOuterFace, layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2);
    }
    const auto turnsCopy = coil.get_turns_description().value();
    for (const auto& turn : turnsCopy) {
        frame.sectionOfTurn[turn.get_name()] = turn.get_section().value();
    }
    return frame;
}

// Every terminal reservation of `winding` (its lead, stub and the squeezes it charges) stays at the
// height of the section its terminal turn lies in, and squeezes only that section's layers.
void check_terminals_local(const std::vector<OpenMagnetics::ConnectionReservedSpace>& spaces, const SectionFrame& frame,
                           const std::string& winding) {
    size_t checked = 0;
    for (const auto& space : spaces) {
        if (!space.isTerminal || space.winding != winding) {
            continue;
        }
        // An entrance's reservations name the turn they reach (toTurn), an exit's the turn they leave.
        const std::string& terminalTurn = space.toTurn.empty() ? space.fromTurn : space.toTurn;
        INFO(winding << " p" << space.parallel << " terminal space of turn '" << terminalTurn << "' in '" << space.section
                     << "' layer '" << space.layer << "' at y " << space.coordinates[1] * 1e3 << " mm");
        REQUIRE(frame.sectionOfTurn.contains(terminalTurn));
        const std::string& section = frame.sectionOfTurn.at(terminalTurn);
        CHECK(space.section == section);
        const auto [low, high] = frame.axialExtent.at(section);
        CHECK(space.coordinates[1] - space.dimensions[1] / 2 >= low - 1e-9);
        CHECK(space.coordinates[1] + space.dimensions[1] / 2 <= high + 1e-9);
        if (space.kind == OpenMagnetics::ConnectionKind::LAYER_SQUEEZE) {
            CHECK(frame.sectionOfLayer.at(space.layer) == section);
        }
        ++checked;
    }
    REQUIRE(checked > 0);
}

}  // namespace

// ABT #1487 (Alf, 2026-09-29), rules R1-R4 on the PQ 65/60 S-P-S-P-S-P-S design, its primary cut to
// 18 turns (6 per section): 24 turns genuinely overflow radially, 7 layers of the 1.865 mm litz
// being 13.06 mm of build in the 12.05 mm window (an overflow wind() reports, tests below). R1: a section
// holding a whole parallel (every secondary section) keeps that parallel's terminals local. R4: the
// series-split primary's own start and finish are local too, in its first and last section. R2: the
// primary's links between its sections run OUTSIDE the winding build, at a lane (exitSlot) of their
// own, squeezing only the two sections they join. R3: the build (26.6 mm, the secondary's edge-wound
// width) leaves 0.86 mm to the window's outer boundary and the 1.865 mm litz link does not fit there;
// it runs where the PQ's outer legs leave the window open, within the core outline.
TEST_CASE("Real winding on axially stacked sections: local terminals, outside series links (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto magnetic = pq65_through_real_winding([](json& masJson) { pq65_winding(masJson, "Primary")["numberTurns"] = 18; });
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto frame = section_frame(coil);
    const auto spaces = coil.get_connection_reserved_spaces();
    const auto routes = coil.get_connection_layout().routes;

    check_terminals_local(spaces, frame, "Secondary");   // R1
    check_terminals_local(spaces, frame, "Primary");     // R4

    const auto window = coil.resolve_bobbin().get_processed_description()->get_winding_windows()[0];
    const double windowOuterFace = window.get_coordinates().value()[0] + window.get_width().value() / 2;
    std::vector<double> terminalSlots;
    for (const auto& route : routes) {
        if ((route.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE || route.kind == OpenMagnetics::ConnectionKind::TERMINAL_EXIT) &&
            route.exitSlot) {
            terminalSlots.push_back(route.exitSlot.value());
        }
    }
    REQUIRE(!terminalSlots.empty());
    size_t outsideLinks = 0;
    for (const auto& route : routes) {
        if (route.kind != OpenMagnetics::ConnectionKind::EDGE_CONTINUATION || !route.exitSlot) {
            continue;
        }
        ++outsideLinks;
        INFO("outside link " << route.winding << " p" << route.parallel << " " << route.fromTurn << " -> " << route.toTurn);
        CHECK(route.winding == "Primary");
        const std::set<std::string> joined{frame.sectionOfTurn.at(route.fromTurn), frame.sectionOfTurn.at(route.toTurn)};
        CHECK(joined.size() == 2);
        // Outboard of every section's copper, and (R3) beyond the window here.
        // Radial out, axial outside the build, radial in -- plus a stub at either end when the
        // turn does not already sit on its section-edge row (a turn on the row has none).
        REQUIRE(route.waypoints.size() >= 4);
        double linkX = std::numeric_limits<double>::lowest();
        for (const auto& waypoint : route.waypoints) {
            linkX = std::max(linkX, waypoint[0]);
        }
        const double linkWidth = OpenMagnetics::resolve_dimensional_values(
            coil.resolve_wire(coil.get_winding_index_by_name("Primary")).get_maximum_outer_width());
        CHECK(linkX - linkWidth / 2 >= frame.copperOuterFace - 1e-9);
        CHECK(linkX + linkWidth / 2 > windowOuterFace);
        // A lane of its own, at a different angle from every terminal.
        for (double slot : terminalSlots) {
            CHECK(std::abs(route.exitSlot.value() - slot) > 1e-6);
        }
        // It squeezes, and reserves, nothing outside the two sections it joins.
        for (const auto& space : spaces) {
            if (space.winding != route.winding || space.fromTurn != route.fromTurn || space.toTurn != route.toTurn) {
                continue;
            }
            CHECK(joined.contains(space.section));
            if (space.kind == OpenMagnetics::ConnectionKind::LAYER_SQUEEZE) {
                CHECK(joined.contains(frame.sectionOfLayer.at(space.layer)));
            }
        }
    }
    CHECK(outsideLinks == 2);   // Primary section 0 -> 1 -> 2
    settings.reset();
}

// ABT #1487 R5 (Alf, 2026-09-29): a STORED design's section heights -- a customer may have corrected
// them -- survive magnetic_autocomplete's real-winding re-wind. Only a section whose crossing stations
// need more height than it was given grows, by exactly that need, and the height comes from its
// neighbouring sections; everything else stays as stored. Primary here is 0.5 mm round wire, so the
// primary sections have height to give.
TEST_CASE("magnetic_autocomplete keeps stored heights of axially stacked sections through the real-winding re-wind (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    auto heights = [](const OpenMagnetics::Coil& coil) {
        std::vector<std::pair<std::string, double>> result;
        const auto sectionsCopy = coil.get_sections_description().value();
        for (const auto& section : sectionsCopy) {
            if (section.get_type() == MAS::ElectricalType::CONDUCTION) {
                result.push_back({section.get_name(), section.get_dimensions()[1]});
            }
        }
        return result;
    };
    auto roundPrimary = [](json& masJson) { pq65_winding(masJson, "Primary")["wire"] = "Round 0.5 - Grade 1"; };

    SECTION("heights that hold their crossings are kept as stored") {
        settings.reset();
        auto magnetic = pq65_through_real_winding([&](json& masJson) {
            roundPrimary(masJson);
            set_pq65_stored_heights(masJson, 0.0036, 0.0068);
        });
        REQUIRE(magnetic.get_coil().is_real_winding_blocking_applied());
        const auto kept = heights(magnetic.get_coil());
        REQUIRE(kept.size() == 7);
        for (const auto& [name, height] : kept) {
            INFO(name);
            CHECK(std::abs(height - (name.rfind("Secondary", 0) == 0 ? 0.0036 : 0.0068)) < 1e-9);
        }
    }
    SECTION("a section too short for its crossings grows by exactly their need, from its neighbours") {
        settings.reset();
        auto magnetic = pq65_through_real_winding([&](json& masJson) {
            roundPrimary(masJson);
            set_pq65_stored_heights(masJson, 0.002918, 0.007591);
        });
        auto& coil = magnetic.get_mutable_coil();
        REQUIRE(coil.is_real_winding_blocking_applied());
        // 3 turns of one parallel in one layer (the 10.085 mm edge-wound wire fills the width): 3 turns
        // plus that layer's crossing station, one 0.885 mm wire height each.
        const double wireHeight = OpenMagnetics::resolve_dimensional_values(
            coil.resolve_wire(coil.get_winding_index_by_name("Secondary")).get_maximum_outer_height());
        const double need = 4 * wireHeight;
        const double grow = need - 0.002918;
        REQUIRE(grow > 0);
        // S0 and S3 take all from their one neighbour, S1 and S2 half from each of theirs.
        const std::vector<std::pair<std::string, double>> expected{
            {"Secondary section 0", need}, {"Primary section 0", 0.007591 - grow - grow / 2},
            {"Secondary section 1", need}, {"Primary section 1", 0.007591 - grow / 2 - grow / 2},
            {"Secondary section 2", need}, {"Primary section 2", 0.007591 - grow / 2 - grow},
            {"Secondary section 3", need}};
        const auto kept = heights(coil);
        REQUIRE(kept.size() == expected.size());
        for (size_t index = 0; index < expected.size(); ++index) {
            INFO(expected[index].first);
            CHECK(kept[index].first == expected[index].first);
            CHECK(std::abs(kept[index].second - expected[index].second) < 1e-9);
        }
    }
    settings.reset();
}

namespace {

// The PQ 65/60 field design of ABT #1487 wound directly, as the web's first wind does: no stored
// sections, the core set for the outside links' room (R3), the secondary wound by consecutive turns.
OpenMagnetics::Coil pq65_direct_coil(std::optional<int64_t> primaryTurns) {
    auto masJson = load_test_data("abt1487_pq65_sp_contiguous.json");
    if (primaryTurns) {
        pq65_winding(masJson, "Primary")["numberTurns"] = primaryTurns.value();
    }
    auto coilJson = masJson["magnetic"]["coil"];
    coilJson.erase("sectionsDescription");
    OpenMagnetics::Coil coil(coilJson, false);
    OpenMagnetics::Core core(masJson["magnetic"]["core"]);
    coil.set_core_geometry(core);
    coil.preload_winding_style_overrides({{"Secondary", MAS::WindingStyle::WIND_BY_CONSECUTIVE_TURNS}});
    return coil;
}

void check_every_declared_turn_present(const OpenMagnetics::Coil& coil) {
    std::set<std::string> names;
    const auto turnsCopy = coil.get_turns_description().value();
    for (const auto& turn : turnsCopy) {
        names.insert(turn.get_name());
    }
    for (const auto& winding : coil.get_functional_description()) {
        for (int64_t parallel = 0; parallel < winding.get_number_parallels(); ++parallel) {
            for (int64_t turn = 0; turn < winding.get_number_turns(); ++turn) {
                const std::string name = winding.get_name() + " parallel " + std::to_string(parallel) + " turn " + std::to_string(turn);
                INFO(name);
                CHECK(names.contains(name));
            }
        }
    }
}

}  // namespace

// ABT #1487 (owner decision, Alf 2026-09-29): a real-winding layout that does not fit is never a
// throw. The field design (primary 45 turns of 1.865 mm litz, 15 per section) wound directly needs
// more layers than the 12.05 mm window takes: the overlapping layers grow out radially, every turn
// is laid out, and wind() says it does not fit, naming how far and which way.
TEST_CASE("Real winding lays out an overflowing axially stacked design in full and reports it (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto coil = pq65_direct_coil(std::nullopt);
    settings.set_coil_use_real_winding_geometry(true);
    CHECK_FALSE(coil.wind({1, 0, 1, 0, 1, 0, 1}, 1));
    REQUIRE(coil.get_turns_description());
    check_every_declared_turn_present(coil);
    const std::string failure = coil.get_last_fit_failure();
    INFO(failure);
    CHECK(failure.find("real winding does not fit") != std::string::npos);
    CHECK(failure.find("radially outward") != std::string::npos);
    // It grows out radially, not along the column: no turn passes the window's top or bottom.
    const auto window = coil.resolve_bobbin().get_processed_description()->get_winding_windows()[0];
    const double windowTop = window.get_coordinates().value()[1] + window.get_height().value() / 2;
    const double windowBottom = window.get_coordinates().value()[1] - window.get_height().value() / 2;
    const auto turnsCopy = coil.get_turns_description().value();
    for (const auto& turn : turnsCopy) {
        INFO(turn.get_name());
        const double halfHeight = turn.get_dimensions().value()[1] / 2;
        CHECK(turn.get_coordinates()[1] + halfHeight <= windowTop + 1e-6);
        CHECK(turn.get_coordinates()[1] - halfHeight >= windowBottom - 1e-6);
    }
    settings.reset();
}

// ABT #1487 (R5 + owner decision): the field design as stored (secondary sections 2.918 mm, primary
// 7.591 mm) through magnetic_autocomplete. The secondary sections need 3.54 mm for their crossings
// and the 1.865 mm litz primary sections have almost nothing to give: the stack grows along the
// column past the window, is laid out all the same, and the fit failure says so.
TEST_CASE("magnetic_autocomplete reports a stored stack that outgrows the window along the column (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto magnetic = pq65_through_real_winding([](json& masJson) { set_pq65_stored_heights(masJson, 0.002918, 0.007591); });
    const auto& coil = magnetic.get_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    REQUIRE(coil.get_turns_description());
    check_every_declared_turn_present(coil);
    const std::string failure = coil.get_last_fit_failure();
    INFO(failure);
    CHECK(failure.find("real winding does not fit") != std::string::npos);
    CHECK(failure.find("grows") != std::string::npos);
    CHECK(failure.find("along the column past its stored height") != std::string::npos);
    settings.reset();
}

namespace {

// Every conduction layer, and every turn, lies within the height of its own section.
void check_layers_within_own_sections(const OpenMagnetics::Coil& coil) {
    const auto frame = section_frame(coil);
    // The sections themselves stay stacked: none reaches into another's height.
    std::vector<std::pair<double, double>> extents;
    for (const auto& [name, extent] : frame.axialExtent) {
        extents.push_back(extent);
    }
    std::sort(extents.begin(), extents.end());
    for (size_t index = 1; index < extents.size(); ++index) {
        INFO("section spanning " << extents[index - 1].first * 1e3 << " .. " << extents[index - 1].second * 1e3 << " mm");
        CHECK(extents[index - 1].second <= extents[index].first + 1e-9);
    }
    const auto layersCopy = coil.get_layers_description().value();
    size_t checkedLayers = 0;
    for (const auto& layer : layersCopy) {
        if (layer.get_type() != MAS::ElectricalType::CONDUCTION) {
            continue;
        }
        INFO(layer.get_name());
        const auto [low, high] = frame.axialExtent.at(layer.get_section().value());
        CHECK(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 >= low - 1e-9);
        CHECK(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 <= high + 1e-9);
        ++checkedLayers;
    }
    CHECK(checkedLayers >= 7);
    const auto turnsCopy = coil.get_turns_description().value();
    for (const auto& turn : turnsCopy) {
        INFO(turn.get_name());
        const auto [low, high] = frame.axialExtent.at(turn.get_section().value());
        const double halfHeight = turn.get_dimensions().value()[1] / 2;
        CHECK(turn.get_coordinates()[1] - halfHeight >= low - 1e-9);
        CHECK(turn.get_coordinates()[1] + halfHeight <= high + 1e-9);
    }
}

}  // namespace

// ABT #1487, two layout bugs axially stacked sections exposed, on the field design cut to a 15-turn
// primary, which FITS wound directly (at 18 turns the outside links already find no room).
// (i) get_connection_reserved_spaces counted a layer of ANOTHER section as crossed by a link between
// two layers of one section whenever its centre fell radially between them (the wide edge-wound secondary, centred in the window, sat between the primary's two
// layers): every primary layer step squeezed all four secondary sections and the wind threw "No
// layers in section". (ii) align_blocked_layer_turns spread a blocked layer's turns over the whole
// window height rather than its own section's, laying 34.566 mm layers across the other sections.
TEST_CASE("Real winding on axially stacked sections keeps links and blocked layers in their own section (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto coil = pq65_direct_coil(15);
    settings.set_coil_use_real_winding_geometry(true);
    const bool fits = coil.wind({1, 0, 1, 0, 1, 0, 1}, 1);
    INFO(coil.get_last_fit_failure());
    REQUIRE(fits);
    REQUIRE(coil.is_real_winding_blocking_applied());
    check_every_declared_turn_present(coil);
    const auto frame = section_frame(coil);

    // (i) a link squeezes only layers of the sections its two turns lie in -- never another
    // section's layer at another height (the secondary's layer, centred at x 21.48 mm, lies
    // radially between the primary's section-1 layers at 20.54 and 22.41 mm).
    auto turnOf = [](const std::string& name) { return name.substr(0, name.find('_')); };
    size_t squeezes = 0;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (space.isTerminal || space.kind != OpenMagnetics::ConnectionKind::LAYER_SQUEEZE) {
            continue;
        }
        INFO(space.winding << " " << space.fromTurn << " -> " << space.toTurn << " squeezes '" << space.layer << "'");
        REQUIRE(frame.sectionOfTurn.contains(turnOf(space.fromTurn)));
        REQUIRE(frame.sectionOfTurn.contains(turnOf(space.toTurn)));
        const std::set<std::string> joined{frame.sectionOfTurn.at(turnOf(space.fromTurn)), frame.sectionOfTurn.at(turnOf(space.toTurn))};
        CHECK(joined.contains(frame.sectionOfLayer.at(space.layer)));
        ++squeezes;
    }
    CHECK(squeezes > 0);

    // (ii) every layer, and every turn in it, stays within its own section's height.
    check_layers_within_own_sections(coil);
    settings.reset();

    // (ii) again where lead depths block layer slots: the stored field heights with a 0.5 mm round
    // primary (it fits: the secondary sections grow from their primary neighbours, R5), re-wound
    // with real winding. The bug spread primary section 2's blocked layers over 34.566 mm.
    auto magnetic = pq65_through_real_winding([](json& masJson) {
        pq65_winding(masJson, "Primary")["wire"] = "Round 0.5 - Grade 1";
        set_pq65_stored_heights(masJson, 0.002918, 0.007591);
    });
    REQUIRE(magnetic.get_coil().is_real_winding_blocking_applied());
    CHECK(magnetic.get_coil().get_last_fit_failure().empty());
    check_layers_within_own_sections(magnetic.get_coil());
    settings.reset();
}

// ABT #1487 (Alf, 2026-09-29): "the input terminal should be on top, on the same side the section
// is, as it doesn't make sense it is in the bottom and wind the layers up". On axially stacked
// sections every conductor's START is local (R1, R4), so its section starts at the edge facing the
// window end the section lies on and winds away from it: the start turn sits on that edge's row and
// its lead leaves by that edge. Here the primary's first section is the top-most primary section,
// so its start turn is on that section's top row and its lead leaves at the top; every secondary
// section holds a whole parallel and starts at the edge towards its own window end.
TEST_CASE("Real winding on axially stacked sections starts each local start at the section edge facing its window end (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto magnetic = pq65_through_real_winding([](json& masJson) { pq65_winding(masJson, "Primary")["numberTurns"] = 18; });
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto frame = section_frame(coil);
    const auto turnsCopy = coil.get_turns_description().value();
    const double windowCenterY =
        coil.resolve_bobbin().get_processed_description()->get_winding_windows()[0].get_coordinates().value()[1];

    // Precondition: the primary's first section is its top-most one.
    auto turnNamed = [&](const std::string& name) -> const auto& {
        for (const auto& turn : turnsCopy) {
            if (turn.get_name() == name) {
                return turn;
            }
        }
        FAIL("no turn " << name);
        throw std::runtime_error("unreachable");
    };
    const auto& primaryStart = turnNamed("Primary parallel 0 turn 0");
    const std::string primaryStartSection = frame.sectionOfTurn.at(primaryStart.get_name());
    for (const auto& [section, extent] : frame.axialExtent) {
        if (section.rfind("Primary", 0) == 0) {
            REQUIRE(extent.second <= frame.axialExtent.at(primaryStartSection).second + 1e-9);
        }
    }

    const size_t secondaryParallels = coil.get_number_parallels(coil.get_winding_index_by_name("Secondary"));
    std::vector<std::pair<std::string, int64_t>> starts{{"Primary", 0}};
    for (size_t parallel = 0; parallel < secondaryParallels; ++parallel) {
        starts.push_back({"Secondary", int64_t(parallel)});
    }
    for (const auto& [winding, parallel] : starts) {
        const auto& start = turnNamed(winding + " parallel " + std::to_string(parallel) + " turn 0");
        const std::string section = frame.sectionOfTurn.at(start.get_name());
        const auto [low, high] = frame.axialExtent.at(section);
        const bool sectionFacesTop = (low + high) / 2 >= windowCenterY;
        double edgeRowY = sectionFacesTop ? std::numeric_limits<double>::lowest() : std::numeric_limits<double>::max();
        for (const auto& turn : turnsCopy) {
            if (turn.get_section().value() == section && turn.get_winding() == winding) {
                edgeRowY = sectionFacesTop ? std::max(edgeRowY, turn.get_coordinates()[1]) : std::min(edgeRowY, turn.get_coordinates()[1]);
            }
        }
        INFO(winding << " p" << parallel << " starts at '" << start.get_name() << "' in '" << section << "' [" << low * 1e3
                     << ", " << high * 1e3 << "] mm, turn y " << start.get_coordinates()[1] * 1e3 << " mm, edge row y "
                     << edgeRowY * 1e3 << " mm, facing " << (sectionFacesTop ? "top" : "bottom"));
        CHECK(std::abs(start.get_coordinates()[1] - edgeRowY) < 1e-6);

        size_t entranceRoutes = 0;
        for (const auto& route : coil.get_connection_layout().routes) {
            if (route.kind != OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE || route.winding != winding ||
                route.parallel != parallel) {
                continue;
            }
            ++entranceRoutes;
            for (const auto& waypoint : route.waypoints) {
                INFO("entrance waypoint y " << waypoint[1] * 1e3 << " mm");
                // It leaves by the section's edge facing its window end: never towards the other half.
                CHECK((sectionFacesTop ? waypoint[1] >= start.get_coordinates()[1] - 1e-9
                                       : waypoint[1] <= start.get_coordinates()[1] + 1e-9));
                CHECK((sectionFacesTop ? waypoint[1] >= (low + high) / 2 : waypoint[1] <= (low + high) / 2));
            }
        }
        CHECK(entranceRoutes == 1);
    }
    settings.reset();
}

namespace {

// The turns painted with the fit-problem ring in an SVG (the CSS rule for the class does not count).
size_t count_fit_problem_rings(const std::string& svg) {
    size_t count = 0;
    const std::string needle = "class=\"fit_problem_turn\"";
    for (size_t position = svg.find(needle); position != std::string::npos; position = svg.find(needle, position + 1)) {
        ++count;
    }
    return count;
}

std::string paint_turns_with_fit_problems(OpenMagnetics::Magnetic magnetic, const std::string& fileName) {
    const auto outFile = std::filesystem::path{std::source_location::current().file_name()}
                             .parent_path().append("..").append("output").append(fileName);
    std::filesystem::remove(outFile);
    OpenMagnetics::Painter painter(outFile);
    painter.paint_core(magnetic);
    painter.paint_coil_turns(magnetic);
    painter.paint_winding_fit_problems(magnetic);
    return painter.export_svg();
}

}  // namespace

// The web's 2D view draws a real-winding layout that does not fit as it was wound and marks WHERE it
// overflows: Painter::paint_winding_fit_problems outlines the winding window and rings every turn
// whose copper leaves it. The field design as stored grows along the column past the window; a
// design that fits (0.5 mm round primary) has no ring at all.
TEST_CASE("paint_winding_fit_problems rings the turns an unfit real winding leaves outside the window (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487][painter]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    {
        auto unfit = pq65_through_real_winding([](json& masJson) { set_pq65_stored_heights(masJson, 0.002918, 0.007591); });
        REQUIRE_FALSE(unfit.get_coil().get_last_fit_failure().empty());
        const auto svg = paint_turns_with_fit_problems(unfit, "abt1487_fit_problems_unfit.svg");
        CHECK(svg.find("class=\"fit_problem_window\"") != std::string::npos);
        CHECK(count_fit_problem_rings(svg) > 0);
    }
    settings.reset();
    {
        auto fitting = pq65_through_real_winding([](json& masJson) {
            pq65_winding(masJson, "Primary")["wire"] = "Round 0.5 - Grade 1";
            set_pq65_stored_heights(masJson, 0.002918, 0.007591);
        });
        REQUIRE(fitting.get_coil().get_last_fit_failure().empty());
        const auto svg = paint_turns_with_fit_problems(fitting, "abt1487_fit_problems_fitting.svg");
        CHECK(svg.find("class=\"fit_problem_window\"") != std::string::npos);
        CHECK(count_fit_problem_rings(svg) == 0);
    }
    settings.reset();
}

// wind_magnetic_coil_as_described is the re-wind magnetic_autocomplete shares with the web's
// display paths, and those decide what to draw from its result: it must hand back wind()'s false
// for a layout that does not fit (with the reason in get_last_fit_failure), and true when it fits.
TEST_CASE("wind_magnetic_coil_as_described returns false for a real winding that does not fit (ABT #1487)",
          "[constructive-model][coil][real-winding][abt1487]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto rewind = [&](const std::function<void(json&)>& edit) {
        auto masJson = load_test_data("abt1487_pq65_sp_contiguous.json");
        edit(masJson);
        OpenMagnetics::Mas mas(masJson);
        settings.set_coil_use_real_winding_geometry(true);
        auto magnetic = OpenMagnetics::magnetic_autocomplete(OpenMagnetics::Magnetic(masJson.at("magnetic")), json{}, mas.get_inputs());
        REQUIRE(OpenMagnetics::magnetic_coil_needs_winding(magnetic));
        const bool fits = OpenMagnetics::wind_magnetic_coil_as_described(magnetic, json{}, mas.get_inputs());
        return std::make_pair(fits, magnetic.get_coil().get_last_fit_failure());
    };
    {
        const auto [fits, failure] = rewind([](json& masJson) { set_pq65_stored_heights(masJson, 0.002918, 0.007591); });
        INFO(failure);
        CHECK_FALSE(fits);
        // The reason is named (here a turn of the grown stack lies outside the window).
        CHECK_FALSE(failure.empty());
    }
    settings.reset();
    {
        const auto [fits, failure] = rewind([](json& masJson) {
            pq65_winding(masJson, "Primary")["wire"] = "Round 0.5 - Grade 1";
            set_pq65_stored_heights(masJson, 0.002918, 0.007591);
        });
        INFO(failure);
        CHECK(fits);
        CHECK(failure.empty());
    }
    settings.reset();
}

namespace {

// Autocompletes a stored design the way the web does. With freshWind the stored turns and layers are
// dropped, so the coil is wound from its sections even with real winding off (a stored layout is
// otherwise taken as it is and nothing is wound).
// With dropSections too, the coil states no layout at all and MKF plans it (ABT #1496: such a plan
// that does not fit throws, where a stated one is laid out and warned about).
OpenMagnetics::Magnetic autocomplete_stored_design(const std::string& name, bool freshWind = false, bool dropSections = false) {
    const auto masJson = load_test_data(name);
    json magneticJson = masJson.contains("magnetic") ? masJson.at("magnetic") : masJson;
    if (freshWind) {
        magneticJson["coil"].erase("turnsDescription");
        magneticJson["coil"].erase("layersDescription");
    }
    if (dropSections) {
        magneticJson["coil"].erase("sectionsDescription");
    }
    std::optional<OpenMagnetics::Inputs> inputs;
    if (masJson.contains("inputs")) {
        inputs = OpenMagnetics::Inputs(masJson.at("inputs"), false);
    }
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagnetics::Core(magneticJson.at("core")));
    magnetic.set_coil(OpenMagnetics::Coil(magneticJson.at("coil"), false));
    return OpenMagnetics::magnetic_autocomplete(magnetic, json{}, inputs);
}

size_t count_conduction_layers_of(OpenMagnetics::Coil& coil, const std::string& windingName) {
    size_t count = 0;
    const auto layers = coil.get_layers_by_type(MAS::ElectricalType::CONDUCTION);
    for (const auto& layer : layers) {
        if (layer.get_partial_windings()[0].get_winding() == windingName) {
            count++;
        }
    }
    return count;
}

}  // namespace

// ABT #1496: LLC_trafo_first (ETD 59, Primary 24 + Secondary 8 + Tertiary 8 turns of a 2.64 mm
// wire, 11.475 mm wide window). The equal split gives each winding 3.81 mm and the primary needs
// two 16-turn layers (5.28 mm), which try_rewind cannot hand out because the equal split leaves
// nothing unused; the winder then crammed all 24 turns into ONE layer, turn 0 10.6 mm below the
// window floor, and autocomplete handed that on. 5.28 + 2.64 + 2.64 mm fit side by side.
TEST_CASE("Real winding: autocomplete re-splits a window the equal split filled, so every turn is inside (ABT #1496)",
          "[constructive-model][coil][real-winding][abt1496]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto magnetic = autocomplete_stored_design("LLC_trafo_first.json");
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    INFO(coil.get_last_fit_failure());
    CHECK(coil.are_turns_inside_winding_window());
    CHECK(coil.are_sections_and_layers_fitting());
    CHECK(count_conduction_layers_of(coil, "Primary") == 2);
    settings.reset();
}

// ABT #1496: with real winding on, every re-wind after the ideal one (the crossing raise, the turn
// blocking) starts from the split that fit, not from the caller's equal one.
// catalogue_missing_turns (E 20/10/6): the SEC's 1.048 mm litz fits in the 0.46 of the window the
// re-split gives it, not in the equal third (1.02 mm). turns_misaligned (E 30/15/7): the 3.754 mm
// primary wire needs 0.74 of the window; half of it is 2.52 mm.
TEST_CASE("Real winding: re-winds start from the split that fit (ABT #1496)",
          "[constructive-model][coil][real-winding][abt1496]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    for (const std::string name : {"catalogue_missing_turns.json", "turns_misaligned.json"}) {
        settings.reset();
        settings.set_coil_use_real_winding_geometry(true);
        INFO(name);
        auto magnetic = autocomplete_stored_design(name);
        auto& coil = magnetic.get_mutable_coil();
        REQUIRE(coil.get_turns_description());
        INFO(coil.get_last_fit_failure());
        CHECK(coil.are_turns_inside_winding_window());
        CHECK(coil.are_sections_and_layers_fitting());
    }
    settings.reset();
}

// ABT #1496: an interior section that ends on a layer boundary below its fair share must not leave
// the outer section more than its whole layers hold. abt1194_interleaved_nfilar_e16: the Secondary's
// inner section stopped at one layer, the outer one was handed 40 conductors for its 2 x 17 rows and
// its turns came out 229 um past the window.
TEST_CASE("Real winding: an interior section takes its own whole layers when the outer ones cannot hold the rest (ABT #1496)",
          "[constructive-model][coil][real-winding][abt1496]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto magnetic = autocomplete_stored_design("abt1194_interleaved_nfilar_e16.json");
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    INFO(coil.get_last_fit_failure());
    CHECK(coil.are_turns_inside_winding_window());
    CHECK(coil.are_sections_and_layers_fitting());
    settings.reset();
}

// ABT #1496 (Alf, 2026-10-10: a layout MKF plans that does not fit THROWS, real winding on and off,
// with the fit-failure reason). bug_capacitance_error: 12 turns of a 0.935 mm wire, 2 per layer in
// a 2.8 mm tall window, need 6 layers (5.61 mm) in a 4.505 mm wide one -- no split can fit it. It
// used to come back wound past the window, silently. Its stated sections are dropped too, so MKF
// plans the whole layout.
TEST_CASE("magnetic_autocomplete throws COIL_WINDING_ERROR for a winding that does not fit its window (ABT #1496)",
          "[constructive-model][coil][real-winding][magnetic-autocomplete][abt1496]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    for (const bool realWinding : {false, true}) {
        settings.reset();
        settings.set_coil_use_real_winding_geometry(realWinding);
        INFO("real winding " << realWinding);
        bool threw = false;
        try {
            autocomplete_stored_design("bug_capacitance_error.json", /*freshWind=*/true, /*dropSections=*/true);
        }
        catch (const OpenMagnetics::OpenMagneticsException& e) {
            threw = true;
            INFO(e.what());
            CHECK(e.code() == OpenMagnetics::ErrorCode::COIL_WINDING_ERROR);
            CHECK(std::string(e.what()).find("does not fit its window") != std::string::npos);
        }
        CHECK(threw);
    }
    settings.reset();
}

namespace {

// ABT #2073: a winding refused because its copper does not fit must say so, with the numbers, before the
// position of the first turn that escapes. Returns the COIL_WINDING_ERROR text of autocompleting the design.
std::string autocomplete_refusal(const json& core, const json& functionalDescription) {
    json coil = {{"bobbin", "Basic"}, {"functionalDescription", functionalDescription}};
    try {
        OpenMagnetics::magnetic_autocomplete(OpenMagnetics::Magnetic(json{{"core", core}, {"coil", coil}}));
    }
    catch (const OpenMagnetics::OpenMagneticsException& e) {
        CHECK(e.code() == OpenMagnetics::ErrorCode::COIL_WINDING_ERROR);
        return e.what();
    }
    return "";
}

json ungapped_core(const std::string& shape, const std::string& material) {
    return {{"name", shape + " " + material},
            {"functionalDescription", {{"type", "two-piece set"}, {"shape", shape}, {"material", material},
                                       {"gapping", json::array()}, {"numberStacks", 1}}}};
}

}  // namespace

// Bartecka, Klos, Paska (Energies 2024, 17, 55): RM 14, 10 turns of Litz 180x0.1 and 4 of Litz 400x0.1, round
// bundles. The primary needs two 1.932 mm layers and the secondary one 2.880 mm layer, 6.744 mm side by side, and the
// RM 14 former is 6.0 mm across: no split of the window holds them, and the refusal used to name only a turn
// 635 um past the window, which read like a winder that never tried a second layer.
TEST_CASE("Autocomplete refusal names the whole layers that do not fit side by side (ABT #2073)",
          "[constructive-model][coil][magnetic-autocomplete][abt2073]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    json functionalDescription = json::array();
    functionalDescription.push_back({{"name", "Primary"}, {"numberTurns", 10}, {"numberParallels", 1}, {"isolationSide", "primary"},
                                     {"wire", "Litz 180x0.1 - Grade 1 - Unserved"}});
    functionalDescription.push_back({{"name", "Secondary"}, {"numberTurns", 4}, {"numberParallels", 1}, {"isolationSide", "secondary"},
                                     {"wire", {{"name", "Litz 400x0.1 bare"}, {"type", "litz"}, {"strand", "Round 0.1 - Grade 1"},
                                               {"numberConductors", 400}, {"coating", {{"type", "bare"}}}}}});
    const std::string what = autocomplete_refusal(ungapped_core("RM 14", "N97"), functionalDescription);
    INFO(what);
    CHECK(what.find("the whole layers of the windings need") != std::string::npos);
    CHECK(what.find("it gives 6.000 mm") != std::string::npos);
    CHECK(what.find("Primary section 0: 2 layers of up to 9 turns, 3.864 mm") != std::string::npos);
    CHECK(what.find("as laid out, Primary parallel 0 turn 0") != std::string::npos);
    settings.reset();
}

// Power Integrations DER-984: PQ 32/20 with 26 turns of served Litz 200x#42, one winding. Its turns need more layers
// than the section is across, so wind_by_rectangular_layers crams the rest along the column; the refusal says how many
// layers the turns need and how many fit.
TEST_CASE("Autocomplete refusal of a single winding names the layers its turns need (ABT #2073)",
          "[constructive-model][coil][magnetic-autocomplete][abt2073]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    json functionalDescription = json::array();
    functionalDescription.push_back({{"name", "Boost inductor"}, {"numberTurns", 26}, {"numberParallels", 1}, {"isolationSide", "primary"},
                                     {"wire", {{"name", "Litz 200x#42 served"}, {"type", "litz"}, {"strand", "Round 42.0 - Single Build"},
                                               {"numberConductors", 200}, {"coating", {{"type", "served"}, {"numberLayers", 1}}}}}});
    const std::string what = autocomplete_refusal(ungapped_core("PQ 32/20", "PC95"), functionalDescription);
    INFO(what);
    CHECK(what.find("section 'Boost inductor section 0' needs") != std::string::npos);
    CHECK(what.find("room for") != std::string::npos);
    CHECK(what.find("as laid out, Boost inductor parallel 0 turn 0") != std::string::npos);
    settings.reset();
}

namespace {

// Autocompletes `name` with the log collector on and returns the magnetic and every WARNING the
// Utils module logged on the way.
std::pair<OpenMagnetics::Magnetic, std::vector<std::string>> autocomplete_collecting_warnings(const std::string& name,
                                                                                            bool freshWind) {
    auto& logger = OpenMagnetics::Logger::getInstance();
    logger.disableCollector();
    logger.enableCollector(OpenMagnetics::LogLevel::WARNING);
    std::vector<std::string> warnings;
    OpenMagnetics::Magnetic magnetic;
    try {
        magnetic = autocomplete_stored_design(name, freshWind);
    }
    catch (...) {
        logger.disableCollector();
        throw;
    }
    for (const auto& record : logger.drainCollected()) {
        if (record.level == OpenMagnetics::LogLevel::WARNING && record.moduleOfOrigin == "Utils") {
            warnings.push_back(record.message);
        }
    }
    logger.disableCollector();
    return {magnetic, warnings};
}

bool any_contains(const std::vector<std::string>& messages, const std::string& needle) {
    return std::any_of(messages.begin(), messages.end(),
                       [&](const std::string& message) { return message.find(needle) != std::string::npos; });
}

}  // namespace

// ABT #1496 (Alf, 2026-10-10): a layout the design STATES (sections, layers or turns of its own) is
// the customer's, so an unfit one is laid out in full as stated, for the views to ring the turns
// outside the window, and the reason is logged as a WARNING for the caller's log collector -- real
// winding on and off. bug_capacitance_error states its sections; the 6 layers it needs (5.61 mm) do
// not fit the 4.505 mm window.
TEST_CASE("magnetic_autocomplete lays an unfit stated layout out in full and warns why (ABT #1496)",
          "[constructive-model][coil][real-winding][magnetic-autocomplete][abt1496]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    for (const bool realWinding : {false, true}) {
        settings.reset();
        settings.set_coil_use_real_winding_geometry(realWinding);
        INFO("real winding " << realWinding);
        auto [magnetic, warnings] = autocomplete_collecting_warnings("bug_capacitance_error.json", /*freshWind=*/true);
        for (const auto& warning : warnings) {
            UNSCOPED_INFO(warning);
        }
        CHECK(any_contains(warnings, "The stated winding layout is laid out as stated, but the winding does not fit its window"));
        auto& coil = magnetic.get_mutable_coil();
        REQUIRE(coil.get_turns_description());
        // Real winding lays each (layer, parallel) out with one extra, length-free crossing
        // station (ABT #685), so the turns stated are the turns that carry length.
        const auto laidOut = coil.get_turns_description().value();
        size_t turnsLaidOut = std::count_if(laidOut.begin(), laidOut.end(),
                                            [](const auto& turn) { return turn.get_length() > 0; });
        size_t turnsStated = 0;
        for (const auto& winding : coil.get_functional_description()) {
            turnsStated += winding.get_number_turns() * winding.get_number_parallels();
        }
        CHECK(turnsLaidOut == turnsStated);
        CHECK_FALSE((coil.are_turns_inside_winding_window() && coil.are_sections_and_layers_fitting()));
    }
    settings.reset();
}

// ABT #2071 (Alf, 2026-10-26): an unfit wind whose blocked fill cannot close OVERFLOWS ALONG THE TURN
// AXIS. bug_capacitance_error (EL 18, 11 turns of 0.935 mm litz, stated sections, real winding): a
// 2.8 mm layer holds 2 stations, the entrance lead's radial corridor takes 1 of them in every layer it
// crosses and the other is that layer's closing station, so no layer the winding spills into holds a
// turn. Those layers grow taller than the window instead, at the wire's pitch and clear of the
// corridor, to 3 stations each (corridor + closing station + the 2 turns a layer of this window
// holds): the fill closes in the 6 layers the ideal wind needs, blocking applied, every stated turn
// laid out, and the copper past the window's top reported as TURNS_OUTSIDE_WINDOW.
TEST_CASE("Real winding: an unfit wind whose layers cannot close beside the lead corridor grows them past the window",
          "[constructive-model][coil][real-winding][magnetic-autocomplete][abt-2071]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto [magnetic, warnings] = autocomplete_collecting_warnings("bug_capacitance_error.json", /*freshWind=*/true);
    for (const auto& warning : warnings) {
        UNSCOPED_INFO(warning);
    }
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    CHECK(coil.is_real_winding_blocking_applied());
    CHECK_FALSE(coil.get_connection_reserved_spaces().empty());

    const auto laidOut = coil.get_turns_description().value();
    size_t turnsLaidOut = std::count_if(laidOut.begin(), laidOut.end(),
                                        [](const auto& turn) { return turn.get_length() > 0; });
    size_t turnsStated = 0;
    for (const auto& winding : coil.get_functional_description()) {
        turnsStated += winding.get_number_turns() * winding.get_number_parallels();
    }
    CHECK(turnsStated == 11);
    CHECK(turnsLaidOut == turnsStated);

    // 6 layers, as the ideal wind needs: one station per layer plus the 11 turns = 17 stations, layer
    // 0 (no corridor, no closing station) with 2 and every other layer with 3, one more than the
    // window holds.
    size_t conductionLayers = 0;
    std::map<std::string, size_t> stationsPerLayer;
    const auto layersLaidOut = coil.get_layers_description().value();
    for (const auto& layer : layersLaidOut) {
        if (layer.get_type() == MAS::ElectricalType::CONDUCTION) {
            ++conductionLayers;
        }
    }
    for (const auto& turn : laidOut) {
        REQUIRE(turn.get_layer());
        stationsPerLayer[turn.get_layer().value()]++;
    }
    CHECK(conductionLayers == 6);
    CHECK(laidOut.size() == 17);
    size_t grownLayers = 0;
    for (const auto& [layerName, stations] : stationsPerLayer) {
        INFO(layerName << ": " << stations << " stations");
        CHECK(stations <= 3);
        grownLayers += stations == 3 ? 1 : 0;
    }
    CHECK(grownLayers == 5);

    // Exact geometry: no two stations overlap, and the copper past the window is reported, on its
    // top (the layers grew past the edge the corridor does not cross).
    CHECK(coil.judge_copper_overlaps().empty());
    const auto outside = coil.judge_turns_outside_window();
    CHECK(std::any_of(outside.begin(), outside.end(), [](const OpenMagnetics::TurnOutsideWindow& turn) {
        return turn.side == OpenMagnetics::TurnOutsideWindow::Side::TOP;
    }));
    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    size_t outsideReported = 0;
    for (const auto& finding : OpenMagnetics::ResultsAnalyzer::analyze_results(mas)) {
        outsideReported += finding.code == OpenMagnetics::ResultFindingCode::TURNS_OUTSIDE_WINDOW ? 1 : 0;
    }
    CHECK(outsideReported == 1);
    settings.reset();
}

// ABT #1496: a stated layout kept as it is, with no wind (planar never re-winds), is judged the same
// way. OM Oyang paper example: its stored turns lie outside the winding window, and it used to come
// back from autocomplete as it was, with nothing said.
TEST_CASE("magnetic_autocomplete warns when a stated layout it keeps as is does not fit (ABT #1496)",
          "[constructive-model][coil][planar][magnetic-autocomplete][abt1496]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto [magnetic, warnings] = autocomplete_collecting_warnings("OM Oyang paper example.json", /*freshWind=*/false);
    for (const auto& warning : warnings) {
        UNSCOPED_INFO(warning);
    }
    REQUIRE(magnetic.get_mutable_coil().is_planar());
    CHECK(any_contains(warnings, "The stated winding layout is kept as stated, but it does not fit its window"));
    settings.reset();
}

// ABT #1496: autocomplete never hands back an unfit layout without saying so. bug_nan_losses_per_turn
// (real winding) came back judged fitting by wind() and reporting itself unfit right after: the
// lead charge, applied after the verdict, counted Secondary 1 Half 2's 3.5 mm lead rows at each
// edge inside its section, which delimit_and_compact had shrunk to the 9.06 mm its layers keep
// clear of them (section filling factor 1.21). Measured inside the section, the leads take none of
// it and the layout fits, with nothing to warn about.
TEST_CASE("Real winding: autocomplete hands back a layout that fits, or says why not (ABT #1496)",
          "[constructive-model][coil][real-winding][magnetic-autocomplete][abt1496]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto [magnetic, warnings] = autocomplete_collecting_warnings("bug_nan_losses_per_turn.json", /*freshWind=*/false);
    for (const auto& warning : warnings) {
        UNSCOPED_INFO(warning);
    }
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    INFO(coil.get_last_fit_failure());
    CHECK_FALSE(any_contains(warnings, "does not fit its window"));
    CHECK(coil.are_turns_inside_winding_window());
    CHECK(coil.are_sections_and_layers_fitting());
    settings.reset();
}

// ABT #1496 (Alf, 2026-10-10: sections of one winding may differ in width, only when equal ones do
// not fit). bug_3d_foil_2kw_ih (E 55/28/21, 8.75 mm wide window), Secondary - Primary foil -
// Secondary: with real winding on, the Secondary is raised from 36 to 37 turns x 3 parallels of a
// 1.205 mm litz, which split 19 + 18 per parallel: 57 and 54 physical turns at 28 per layer, so one
// section needs 3 layers and the other 2. Equal sections claim 2 x 3 layers (2 x 3.61 mm + the
// 1.6 mm foil, past the window); unequal ones fit (3.61 + 2.41 + 1.6 mm). The design itself, with
// real winding on, then stops at the inner foil's terminal lead, which crosses the outer Secondary
// as a run as tall as the sheet (see the ABT #1496 report); so the split is pinned here with the 37
// turns stated and real winding off, where no lead is laid. MKF plans the layout (the stated
// sections are dropped) in the file's own Secondary, Primary, Secondary order.
TEST_CASE("The sections of one winding take the width of their own whole layers when equal ones do not fit (ABT #1496)",
          "[constructive-model][coil][abt1496]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(false);
    const auto masJson = load_test_data("bug_3d_foil_2kw_ih.json");
    json magneticJson = masJson.at("magnetic");
    magneticJson["coil"].erase("turnsDescription");
    magneticJson["coil"].erase("layersDescription");
    magneticJson["coil"].erase("sectionsDescription");
    REQUIRE(magneticJson["coil"]["functionalDescription"][1]["name"] == "Secondary");
    magneticJson["coil"]["functionalDescription"][1]["numberTurns"] = 37;
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagnetics::Core(magneticJson.at("core")));
    magnetic.set_coil(OpenMagnetics::Coil(magneticJson.at("coil"), false));
    auto done = OpenMagnetics::magnetic_autocomplete(magnetic, json{{"interleavingPattern", {1, 0, 1}}},
                                                     OpenMagnetics::Inputs(masJson.at("inputs"), false));
    auto& coil = done.get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    INFO(coil.get_last_fit_failure());
    CHECK(coil.are_turns_inside_winding_window());
    CHECK(coil.are_sections_and_layers_fitting());
    size_t layersInSection0 = 0;
    size_t layersInSection1 = 0;
    for (const auto& layer : coil.get_layers_by_type(MAS::ElectricalType::CONDUCTION)) {
        if (layer.get_section() == "Secondary section 0") layersInSection0++;
        if (layer.get_section() == "Secondary section 1") layersInSection1++;
    }
    CHECK(layersInSection0 + layersInSection1 == 5);
    CHECK(layersInSection0 != layersInSection1);
    const auto section0 = coil.get_section_by_name("Secondary section 0");
    const auto section1 = coil.get_section_by_name("Secondary section 1");
    CHECK(section0.get_dimensions()[0] != Catch::Approx(section1.get_dimensions()[0]).epsilon(1e-3));
    settings.reset();
}

// ABT #1496 (Alf, 2026-10-11): a foil sheet is not routed out across the outer windings as its own
// lead. A winding shop joins a separate conductor to the foil end (a wire soldered to it, a strip
// welded across it) and takes that out radially, or lets the end leave axially through the end
// face. MAS states it on the connection (leadConductor, exit). The real-winding lead model reserves
// the stated conductor's own cross-section, and a foil end that has to cross outer layers with
// neither stated throws: the sheet (a run as tall as the foil) is no size anybody can reserve.
// 17_cllc_xfmr_e5528_3c92a with its second winding made 4 turns of Foil 0.2: one sheet per layer,
// so its start, on the innermost of its layers, crosses the other three.
namespace {

OpenMagnetics::Magnetic foil_lead_magnetic(const std::vector<MAS::ConnectionElement>& foilConnections) {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    // Laid out ideally (it fits), then read under the real-winding lead model.
    auto mas = OpenMagneticsTesting::mas_loader(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/')) +
                                                "/../MAS/examples/17_cllc_xfmr_e5528_3c92a.json");
    auto magneticIn = mas.get_magnetic();
    auto& windings = magneticIn.get_mutable_coil().get_mutable_functional_description();
    windings[1].set_wire("Foil 0.2");
    windings[1].set_number_turns(4);
    windings[1].set_number_parallels(1);
    if (foilConnections.empty()) {
        windings[1].set_connections(std::nullopt);
    }
    else {
        windings[1].set_connections(foilConnections);
    }
    magneticIn.get_mutable_coil().set_turns_description(std::nullopt);
    magneticIn.get_mutable_coil().set_layers_description(std::nullopt);
    magneticIn.get_mutable_coil().set_sections_description(std::nullopt);
    auto magnetic = OpenMagnetics::magnetic_autocomplete(magneticIn);
    settings.set_coil_use_real_winding_geometry(true);
    return magnetic;
}

MAS::ConnectionElement foil_end(MAS::End end, std::optional<std::string> leadConductorName, std::optional<MAS::LeadExit> exit) {
    MAS::ConnectionElement connection;
    connection.set_end(end);
    if (leadConductorName) {
        // Given like a winding's own wire: here by its catalogue name.
        connection.set_lead_conductor(MAS::WireDataOrNameUnion(leadConductorName.value()));
    }
    connection.set_exit(exit);
    return connection;
}

// The same foil winding wound by the real winder itself.
OpenMagnetics::Magnetic foil_lead_magnetic_wound_real(const std::vector<MAS::ConnectionElement>& foilConnections) {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto mas = OpenMagneticsTesting::mas_loader(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/')) +
                                                "/../MAS/examples/17_cllc_xfmr_e5528_3c92a.json");
    auto magneticIn = mas.get_magnetic();
    auto& windings = magneticIn.get_mutable_coil().get_mutable_functional_description();
    windings[1].set_wire("Foil 0.2");
    windings[1].set_number_turns(4);
    windings[1].set_number_parallels(1);
    windings[1].set_connections(foilConnections);
    magneticIn.get_mutable_coil().set_turns_description(std::nullopt);
    magneticIn.get_mutable_coil().set_layers_description(std::nullopt);
    magneticIn.get_mutable_coil().set_sections_description(std::nullopt);
    settings.set_coil_use_real_winding_geometry(true);
    return OpenMagnetics::magnetic_autocomplete(magneticIn);
}

}  // namespace

TEST_CASE("Real winding: a foil end crossing outer layers with no stated lead conductor throws (ABT #1496)",
          "[constructive-model][coil][real-winding][foil][abt1496]") {
    auto magnetic = foil_lead_magnetic({});
    auto& coil = magnetic.get_mutable_coil();
    std::string message;
    try {
        coil.get_connection_reserved_spaces();
    }
    catch (const OpenMagnetics::CoilException& e) {
        CHECK(e.code() == OpenMagnetics::ErrorCode::COIL_WINDING_ERROR);
        message = e.what();
    }
    INFO(message);
    CHECK(message.find("sheet is not routed out as its own lead") != std::string::npos);
    CHECK(message.find("connections[].leadConductor") != std::string::npos);
    OpenMagnetics::Settings::GetInstance().reset();
}

TEST_CASE("Real winding: a foil end's lead reserves its stated lead conductor, not the sheet (ABT #1496)",
          "[constructive-model][coil][real-winding][foil][abt1496]") {
    const std::string leadName = "Round 1.00 - Grade 1";
    auto magnetic = foil_lead_magnetic({foil_end(MAS::End::START, leadName, std::nullopt),
                                        foil_end(MAS::End::FINISH, leadName, std::nullopt)});
    auto& coil = magnetic.get_mutable_coil();
    const double leadOuterDiameter = OpenMagnetics::Wire(OpenMagnetics::find_wire_by_name(leadName)).get_maximum_outer_width();
    const auto foilName = coil.get_functional_description()[1].get_name();
    const double foilHeight = coil.resolve_wire(1).get_maximum_outer_height();
    REQUIRE(foilHeight > 5 * leadOuterDiameter);
    size_t squeezes = 0;
    size_t terminalRects = 0;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (space.winding != foilName || !space.isTerminal) {
            continue;
        }
        ++terminalRects;
        INFO("kind " << int(space.kind) << " layer '" << space.layer << "' " << space.dimensions[0] * 1e3 << " x "
                     << space.dimensions[1] * 1e3 << " mm");
        if (space.kind == OpenMagnetics::ConnectionKind::LAYER_SQUEEZE) {
            ++squeezes;
            // The crossed layer loses the lead conductor's own slot, not a sheet-tall band.
            CHECK(std::abs(space.dimensions[1] - leadOuterDiameter) < 1e-9);
        }
        else {
            CHECK(std::min(space.dimensions[0], space.dimensions[1]) <= leadOuterDiameter + 1e-9);
        }
    }
    CHECK(terminalRects > 0);
    CHECK(squeezes > 0);
    OpenMagnetics::Settings::GetInstance().reset();
}

TEST_CASE("Real winding: a foil end's lead conductor named outside the wire database throws (ABT #1496)",
          "[constructive-model][coil][real-winding][foil][abt1496]") {
    // Wherever the lead is first read (the layout or the lead model), the unknown name throws.
    std::string message;
    try {
        auto magnetic = foil_lead_magnetic({foil_end(MAS::End::START, std::string("No such lead wire"), std::nullopt),
                                            foil_end(MAS::End::FINISH, std::string("No such lead wire"), std::nullopt)});
        magnetic.get_mutable_coil().get_connection_reserved_spaces();
    }
    catch (const OpenMagnetics::WireNotFoundException& e) {
        message = e.what();
    }
    INFO(message);
    CHECK(message.find("lead conductor of winding") != std::string::npos);
    CHECK(message.find("No such lead wire") != std::string::npos);
    OpenMagnetics::Settings::GetInstance().reset();
}

TEST_CASE("Real winding: a foil end that exits axially crosses no outer layer (ABT #1496)",
          "[constructive-model][coil][real-winding][foil][abt1496]") {
    auto magnetic = foil_lead_magnetic({foil_end(MAS::End::START, std::nullopt, MAS::LeadExit::AXIAL),
                                        foil_end(MAS::End::FINISH, std::nullopt, MAS::LeadExit::AXIAL)});
    auto& coil = magnetic.get_mutable_coil();
    const auto foilName = coil.get_functional_description()[1].get_name();
    std::vector<OpenMagnetics::ConnectionRoute> routes;
    auto spaces = coil.get_connection_reserved_spaces(&routes);
    size_t terminalRects = 0;
    for (const auto& space : spaces) {
        if (space.winding != foilName || !space.isTerminal) {
            continue;
        }
        ++terminalRects;
        CHECK(space.kind != OpenMagnetics::ConnectionKind::LAYER_SQUEEZE);
    }
    CHECK(terminalRects == 2);
    size_t terminalRoutes = 0;
    for (const auto& route : routes) {
        if (route.winding != foilName || (route.kind != OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE &&
                                          route.kind != OpenMagnetics::ConnectionKind::TERMINAL_EXIT)) {
            continue;
        }
        ++terminalRoutes;
        // Along its own column: one straight leg, constant radial position.
        REQUIRE(route.waypoints.size() == 2);
        CHECK(std::abs(route.waypoints[0][0] - route.waypoints[1][0]) < 1e-12);
    }
    CHECK(terminalRoutes == 2);
    OpenMagnetics::Settings::GetInstance().reset();
}

namespace {

// ABT #1423 judges each terminal route of a really wound coil for a tangential exit: a bend wider
// than the connection face's half straight cannot leave radially. Each foil terminal route is
// judged here with a 20 mm bend, wider than the E 55/28/21 column's 8.2 mm half straight (the
// coil's other leads keep the bend they were planned for). Returns, per route, whether it left
// tangentially, or the refusal's text when it could not.
std::vector<std::string> judge_foil_leads_with_a_wide_bend(OpenMagnetics::Coil& coil) {
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto foilName = coil.get_functional_description()[1].get_name();
    const auto turns = coil.get_turns_description().value();
    const auto routes = coil.get_connection_layout().routes;
    std::vector<std::string> verdicts;
    for (const auto& route : routes) {
        if (route.winding != foilName || (route.kind != OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE &&
                                          route.kind != OpenMagnetics::ConnectionKind::TERMINAL_EXIT)) {
            continue;
        }
        auto wide = route;
        wide.plannedBendRadius = 0.020;
        try {
            verdicts.push_back(coil.plan_tangent_departure(wide, routes, turns) ? "tangent" : "radial");
        }
        catch (const OpenMagnetics::CoilException& e) {
            verdicts.push_back(e.what());
        }
    }
    return verdicts;
}

}  // namespace

TEST_CASE("Real winding: a foil end's stated lead conductor is judged for a tangential exit like any lead (ABT #1496)",
          "[constructive-model][coil][real-winding][foil][abt1496]") {
    // A foil used to be exempt -- its sheet ends in a tab, not in a bent wire -- but a STATED lead
    // conductor is a bent wire like any other (Alf, 2026-10-11), so its bend is judged: with a bend
    // that cannot fit, it leaves tangentially or the wind is refused for the tangent reason. A
    // radial exit with an unfittable bend is the defect.
    auto magnetic = foil_lead_magnetic_wound_real({foil_end(MAS::End::START, std::string("Round 1.00 - Grade 1"), std::nullopt),
                                                   foil_end(MAS::End::FINISH, std::string("Round 1.00 - Grade 1"), std::nullopt)});
    const auto verdicts = judge_foil_leads_with_a_wide_bend(magnetic.get_mutable_coil());
    CHECK(verdicts.size() == 2);
    for (const auto& verdict : verdicts) {
        INFO(verdict);
        CHECK((verdict == "tangent" || verdict.find("must leave tangentially") != std::string::npos));
    }
    OpenMagnetics::Settings::GetInstance().reset();
}

TEST_CASE("Real winding: a foil end with no stated lead is not judged for a tangential exit (ABT #1496)",
          "[constructive-model][coil][real-winding][foil][abt1496]") {
    // The other half of the rule: a foil with no stated lead leaves by its sheet's tab (here
    // axially, crossing nothing), so even an unfittable bend sends nothing tangentially.
    auto magnetic = foil_lead_magnetic_wound_real({foil_end(MAS::End::START, std::nullopt, MAS::LeadExit::AXIAL),
                                                   foil_end(MAS::End::FINISH, std::nullopt, MAS::LeadExit::AXIAL)});
    const auto verdicts = judge_foil_leads_with_a_wide_bend(magnetic.get_mutable_coil());
    CHECK(verdicts.size() == 2);
    for (const auto& verdict : verdicts) {
        CHECK(verdict == "radial");
    }
    OpenMagnetics::Settings::GetInstance().reset();
}

TEST_CASE("Real winding: a foil end's axial lead lies on the sheet and leaves past the flange (ABT #1907)",
          "[constructive-model][coil][real-winding][foil][abt1907]") {
    // ABT #1907 (Alf, 2026-10-18): a foil end that exits axially on a stated lead conductor is that
    // conductor laid TANGENT to the sheet's outer face along the sheet's height, running on along the
    // turn axis through the end face to one lead diameter past the flange's outer face. The reserved
    // run starts at the sheet's edge, not at the sheet's centre inside its copper, and a winding
    // end's leads stand side by side along the face, one outer diameter each, in parallel order.
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    const std::string leadName = "Round 0.5 - Grade 1";
    auto mas = OpenMagneticsTesting::mas_loader(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/')) +
                                                "/../MAS/examples/17_cllc_xfmr_e5528_3c92a.json");
    auto magneticIn = mas.get_magnetic();
    auto& windings = magneticIn.get_mutable_coil().get_mutable_functional_description();
    windings[1].set_wire("Foil 0.2");
    windings[1].set_number_turns(1);
    windings[1].set_number_parallels(4);
    windings[1].set_connections(std::vector<MAS::ConnectionElement>{foil_end(MAS::End::START, leadName, MAS::LeadExit::AXIAL),
                                                                    foil_end(MAS::End::FINISH, leadName, MAS::LeadExit::AXIAL)});
    magneticIn.get_mutable_coil().set_turns_description(std::nullopt);
    magneticIn.get_mutable_coil().set_layers_description(std::nullopt);
    magneticIn.get_mutable_coil().set_sections_description(std::nullopt);
    settings.set_coil_use_real_winding_geometry(true);
    auto magnetic = OpenMagnetics::magnetic_autocomplete(magneticIn);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.is_real_winding_blocking_applied());
    const auto foilName = coil.get_functional_description()[1].get_name();
    auto sheet = coil.resolve_wire(1);
    REQUIRE(sheet.get_type() == MAS::WireType::FOIL);
    auto lead = OpenMagnetics::Wire(OpenMagnetics::find_wire_by_name(leadName));
    const double leadCopper = lead.get_maximum_conducting_width();
    const double leadOuter = lead.get_maximum_outer_width();
    const double sheetThickness = sheet.get_maximum_conducting_width();
    const double sheetHeight = sheet.get_maximum_conducting_height();
    const double wall = coil.resolve_bobbin().get_processed_description().value().get_wall_thickness();
    REQUIRE(wall > 0);
    const auto turns = coil.get_turns_description().value();
    std::vector<OpenMagnetics::ConnectionRoute> routes;
    const auto spaces = coil.get_connection_reserved_spaces(&routes);
    std::map<OpenMagnetics::ConnectionKind, std::map<int64_t, double>> alongFaceByEnd;
    size_t axialRoutes = 0;
    for (const auto& route : routes) {
        if (route.winding != foilName || (route.kind != OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE &&
                                          route.kind != OpenMagnetics::ConnectionKind::TERMINAL_EXIT)) {
            continue;
        }
        ++axialRoutes;
        REQUIRE(route.axialLead);
        const auto& axial = route.axialLead.value();
        const bool entrance = route.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE;
        const std::string turnName = entrance ? route.toTurn : route.fromTurn;
        const auto turn = std::find_if(turns.begin(), turns.end(), [&](const auto& t) { return t.get_name() == turnName; });
        REQUIRE(turn != turns.end());
        const double sheetX = turn->get_coordinates()[0];
        const double sheetY = turn->get_coordinates()[1];
        INFO(turnName << ": layerAxis " << axial.layerAxis * 1e3 << " mm, sheet at " << sheetX * 1e3 << " mm");
        // Tangent, copper on copper: no gap, no overlap.
        CHECK_THAT(axial.layerAxis - sheetX, Catch::Matchers::WithinAbs((sheetThickness + leadCopper) / 2, 1e-9));
        CHECK_THAT(axial.conductingDiameter, Catch::Matchers::WithinAbs(leadCopper, 1e-12));
        CHECK_THAT(axial.outerDiameter, Catch::Matchers::WithinAbs(leadOuter, 1e-12));
        // The attach is the sheet's whole height, ending at the edge the lead leaves by.
        const double outward = axial.windowEdge > sheetY ? 1.0 : -1.0;
        CHECK_THAT(axial.attachTo - sheetY, Catch::Matchers::WithinAbs(outward * sheetHeight / 2, 1e-9));
        CHECK_THAT(axial.attachFrom - sheetY, Catch::Matchers::WithinAbs(-outward * sheetHeight / 2, 1e-9));
        // Through the flange, one lead diameter past its outer face.
        CHECK_THAT(axial.flangeOuterFace - axial.windowEdge, Catch::Matchers::WithinAbs(outward * wall, 1e-12));
        CHECK_THAT(axial.end - axial.flangeOuterFace, Catch::Matchers::WithinAbs(outward * leadOuter, 1e-12));
        CHECK_THAT(axial.length, Catch::Matchers::WithinAbs(std::abs(axial.end - axial.attachFrom), 1e-12));
        CHECK_FALSE(route.exitSlot);
        CHECK_FALSE(route.tangentDeparture);
        // The projection on the window: the lead's own column, from the attach's far end to the edge.
        REQUIRE(route.waypoints.size() == 2);
        for (const auto& waypoint : route.waypoints) {
            CHECK_THAT(waypoint[0], Catch::Matchers::WithinAbs(axial.layerAxis, 1e-12));
        }
        alongFaceByEnd[route.kind][route.parallel] = axial.alongFace;
        // The weld strip (Alf, 2026-10-19): the copper that makes lead and sheet one conductor.
        // Sized by MKF's rule: sheet centre to lead centreline, sqrt(3)/2 of the lead's copper
        // diameter wide, along the whole attach.
        const auto& weld = axial.weld;
        INFO("weld: layer " << weld.layerFrom * 1e3 << " -> " << weld.layerTo * 1e3 << " mm, width " << weld.width * 1e3
                            << " mm, turn axis " << weld.from * 1e3 << " -> " << weld.to * 1e3 << " mm");
        CHECK_THAT(weld.layerFrom, Catch::Matchers::WithinAbs(sheetX, 1e-9));
        CHECK_THAT(weld.layerTo, Catch::Matchers::WithinAbs(axial.layerAxis, 1e-12));
        CHECK_THAT(weld.width, Catch::Matchers::WithinAbs(std::sqrt(3.0) / 2 * leadCopper, 1e-9));
        CHECK_THAT(weld.from, Catch::Matchers::WithinAbs(axial.attachFrom, 1e-12));
        CHECK_THAT(weld.to, Catch::Matchers::WithinAbs(axial.attachTo, 1e-12));
        // It lies ACROSS the contact line (the sheet's outer face, at the lead's centreline) ...
        const double contactLayer = sheetX + sheetThickness / 2;
        CHECK(weld.layerFrom < contactLayer - 1e-6);
        CHECK(weld.layerTo > contactLayer + 1e-6);
        // ... inside the sheet's copper and the lead's copper on the layer axis ...
        CHECK(weld.layerFrom >= sheetX - sheetThickness / 2 - 1e-12);
        CHECK(weld.layerTo <= axial.layerAxis + leadCopper / 2 + 1e-12);
        // ... on this parallel's own sheet along the face, inside the lead's own column ...
        CHECK(axial.alongFace - weld.width / 2 >= double(route.parallel) * leadOuter - 1e-9);
        CHECK(weld.width < leadCopper);
        // ... and within the sheet's height and the lead's extent on the turn axis.
        const double weldLow = std::min(weld.from, weld.to);
        const double weldHigh = std::max(weld.from, weld.to);
        CHECK(weldHigh - weldLow > 0);
        CHECK(weldLow >= sheetY - sheetHeight / 2 - 1e-9);
        CHECK(weldHigh <= sheetY + sheetHeight / 2 + 1e-9);
        CHECK(weldLow >= std::min(axial.attachFrom, axial.end) - 1e-12);
        CHECK(weldHigh <= std::max(axial.attachFrom, axial.end) + 1e-12);
        // The reserved run starts at the sheet's edge: it never lies inside the sheet's copper.
        size_t reserved = 0;
        for (const auto& space : spaces) {
            if (!space.isTerminal || space.winding != foilName || space.parallel != route.parallel || space.kind != route.kind) {
                continue;
            }
            ++reserved;
            const double low = space.coordinates[1] - space.dimensions[1] / 2;
            const double high = space.coordinates[1] + space.dimensions[1] / 2;
            CHECK_THAT(space.coordinates[0], Catch::Matchers::WithinAbs(axial.layerAxis, 1e-9));
            CHECK((high <= sheetY - sheetHeight / 2 + 1e-9 || low >= sheetY + sheetHeight / 2 - 1e-9));
        }
        CHECK(reserved == 1);
    }
    CHECK(axialRoutes == 8);
    for (const auto& [kind, byParallel] : alongFaceByEnd) {
        REQUIRE(byParallel.size() == 4);
        for (const auto& [parallel, alongFace] : byParallel) {
            CHECK_THAT(alongFace, Catch::Matchers::WithinAbs((double(parallel) + 0.5) * leadOuter, 1e-9));
        }
    }
    settings.reset();
}
