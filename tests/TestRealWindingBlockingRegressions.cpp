// Real-winding connection blocking: designs MVB++ builds as FEM-grade CAD must keep getting
// blocking applied by magnetic_autocomplete. MVB++ refuses a real-winding build whose blocking was
// not applied (magnetic_autocomplete_safe), so a regression here turns its fixtures red.

#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "support/Settings.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <limits>
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
