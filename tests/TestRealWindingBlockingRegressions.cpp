// Real-winding connection blocking: designs MVB++ builds as FEM-grade CAD must keep getting
// blocking applied by magnetic_autocomplete. MVB++ refuses a real-winding build whose blocking was
// not applied (magnetic_autocomplete_safe), so a regression here turns its fixtures red.

#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include "support/Settings.h"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <set>
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
