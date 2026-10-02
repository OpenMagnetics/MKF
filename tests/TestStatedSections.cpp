// Stated sections: a coil that arrives with a sectionsDescription (and no turns) has said where
// each section is, how big it is and which turns it holds. These guard that the winder keeps what
// was stated, and that copper which cannot fit a stated section is an error, never a silent
// overflow. Fixtures: a WE-FC choke on an ET 20 former, one 11.8 x 1.85 mm window holding four
// 1.8 mm sections at y = +5.0, +2.4, -2.4, -5.0 mm (0.8 / 3.0 / 0.8 mm apart along the column).
//   we_fc_et20_four_spaced_sections.json          30 + 30 turns of 0.4076 mm wire, 15 per section:
//                                                  4 layers of 4 fit 1.85 x 1.8 mm.
//   we_fc_et20_four_spaced_sections_overfull.json 20 + 20 turns of 0.5025 mm wire, 10 per section:
//                                                  3 layers of 3 fit, so 10 cannot.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <source_location>
#include <fstream>
#include "support/Settings.h"
#include "support/Utils.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "TestingUtils.h"
#include "json.hpp"

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::ContainsSubstring;

namespace {

json load_stated_sections_fixture(const std::string& name) {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "stated_sections/" + name);
    std::ifstream file(path);
    REQUIRE(file.good());
    return json::parse(file);
}

struct SettingsReset {
    ~SettingsReset() { settings.reset(); }
};

// Every turn of a stated section lies inside that section's stated rect.
void check_turns_inside_stated_sections(const OpenMagnetics::Coil& coil, const std::vector<Section>& stated) {
    REQUIRE(coil.get_turns_description());
    for (const auto& section : stated) {
        const double y0 = section.get_coordinates()[1] - section.get_dimensions()[1] / 2;
        const double y1 = section.get_coordinates()[1] + section.get_dimensions()[1] / 2;
        const double x0 = section.get_coordinates()[0] - section.get_dimensions()[0] / 2;
        const double x1 = section.get_coordinates()[0] + section.get_dimensions()[0] / 2;
        size_t turnsInSection = 0;
        const auto turns = coil.get_turns_description().value();
        for (const auto& turn : turns) {
            if (turn.get_section().value() != section.get_name()) {
                continue;
            }
            ++turnsInSection;
            const auto& c = turn.get_coordinates();
            const auto d = turn.get_dimensions().value();
            INFO(turn.get_name() << " at y=" << c[1] * 1e3 << " mm; section '" << section.get_name()
                 << "' spans y " << y0 * 1e3 << ".." << y1 * 1e3 << " mm");
            CHECK(c[1] - d[1] / 2 >= y0 - 1e-9);
            CHECK(c[1] + d[1] / 2 <= y1 + 1e-9);
            CHECK(c[0] - d[0] / 2 >= x0 - 1e-9);
            CHECK(c[0] + d[0] / 2 <= x1 + 1e-9);
        }
        CHECK(turnsInSection > 0);
    }
}

}  // namespace

TEST_CASE("Test_Autocomplete_Keeps_Stated_Section_Rects", "[constructive-model][coil][stated-sections][smoke-test]") {
    SettingsReset reset;
    auto magneticJson = load_stated_sections_fixture("we_fc_et20_four_spaced_sections.json");
    auto stated = magneticJson["coil"]["sectionsDescription"].get<std::vector<Section>>();
    OpenMagnetics::Magnetic magnetic(magneticJson);

    auto result = magnetic_autocomplete(magnetic);

    auto wound = result.get_coil().get_sections_description_conduction();
    REQUIRE(wound.size() == stated.size());
    for (const auto& section : stated) {
        auto found = std::find_if(wound.begin(), wound.end(), [&](const Section& s) { return s.get_name() == section.get_name(); });
        REQUIRE(found != wound.end());
        INFO(section.get_name());
        for (size_t axis = 0; axis < 2; ++axis) {
            CHECK(std::abs(found->get_coordinates()[axis] - section.get_coordinates()[axis]) < 1e-12);
            CHECK(std::abs(found->get_dimensions()[axis] - section.get_dimensions()[axis]) < 1e-12);
        }
    }
    check_turns_inside_stated_sections(result.get_coil(), stated);
    for (const auto& section : stated) {
        size_t count = 0;
        const auto turns = result.get_coil().get_turns_description().value();
        for (const auto& turn : turns) {
            count += turn.get_section().value() == section.get_name();
        }
        CHECK(count == 15);
    }
}

TEST_CASE("Test_Autocomplete_Keeps_Stated_Odd_Turn_Split", "[constructive-model][coil][stated-sections][smoke-test]") {
    SettingsReset reset;
    auto magneticJson = load_stated_sections_fixture("we_fc_et20_four_spaced_sections.json");
    // 29 turns per winding: the first section of each winding states 15, the second 14.
    for (auto& winding : magneticJson["coil"]["functionalDescription"]) {
        winding["numberTurns"] = 29;
    }
    std::map<std::string, size_t> statedTurns;
    for (auto& section : magneticJson["coil"]["sectionsDescription"]) {
        bool first = section["name"].get<std::string>().ends_with("section 0");
        section["partialWindings"][0]["parallelsProportion"] = std::vector<double>{(first ? 15.0 : 14.0) / 29.0};
        statedTurns[section["name"]] = first ? 15 : 14;
    }
    auto stated = magneticJson["coil"]["sectionsDescription"].get<std::vector<Section>>();
    OpenMagnetics::Magnetic magnetic(magneticJson);

    auto result = magnetic_autocomplete(magnetic);

    check_turns_inside_stated_sections(result.get_coil(), stated);
    for (const auto& [name, expected] : statedTurns) {
        size_t count = 0;
        const auto turns = result.get_coil().get_turns_description().value();
        for (const auto& turn : turns) {
            count += turn.get_section().value() == name;
        }
        INFO(name);
        CHECK(count == expected);
    }
}

// A stated layout that does not fit is reported as every other unfit wind is (ABT #930: the coil comes
// back, the reason in get_last_fit_failure), with its turns laid out in the stated sections, never
// re-stacked over the window to make them fit.
TEST_CASE("Test_Autocomplete_Reports_When_Stated_Sections_Overflow", "[constructive-model][coil][stated-sections][smoke-test]") {
    SettingsReset reset;
    auto magneticJson = load_stated_sections_fixture("we_fc_et20_four_spaced_sections_overfull.json");
    auto stated = magneticJson["coil"]["sectionsDescription"].get<std::vector<Section>>();
    OpenMagnetics::Magnetic magnetic(magneticJson);
    OpenMagnetics::Magnetic result;
    REQUIRE_NOTHROW(result = magnetic_autocomplete(magnetic));
    const auto failure = result.get_coil().get_last_fit_failure();
    INFO(failure);
    CHECK_THAT(failure, ContainsSubstring("overflows") && ContainsSubstring("needs 2.010 mm") && ContainsSubstring("has 1.800 mm"));
    auto wound = result.get_coil().get_sections_description_conduction();
    REQUIRE(wound.size() == stated.size());
    for (size_t index = 0; index < stated.size(); ++index) {
        INFO(stated[index].get_name());
        for (size_t axis = 0; axis < 2; ++axis) {
            CHECK(std::abs(wound[index].get_coordinates()[axis] - stated[index].get_coordinates()[axis]) < 1e-12);
            CHECK(std::abs(wound[index].get_dimensions()[axis] - stated[index].get_dimensions()[axis]) < 1e-12);
        }
    }
}

TEST_CASE("Test_Autocomplete_Reports_When_Margins_Do_Not_Fit", "[constructive-model][coil][stated-sections][smoke-test]") {
    SettingsReset reset;
    auto magneticJson = load_stated_sections_fixture("we_fc_et20_four_spaced_sections.json");
    // 0.4 + 1.5 mm of margin across a 1.85 mm wide contiguous section: no room for copper.
    std::vector<std::vector<double>> margins{{0, 0.4e-3}, {0.4e-3, 1.5e-3}, {1.5e-3, 0.4e-3}, {0.4e-3, 0}};
    for (size_t i = 0; i < margins.size(); ++i) {
        magneticJson["coil"]["sectionsDescription"][i]["margin"] = margins[i];
    }
    OpenMagnetics::Magnetic magnetic(magneticJson);
    OpenMagnetics::Magnetic result;
    REQUIRE_NOTHROW(result = magnetic_autocomplete(magnetic));
    CHECK_FALSE(result.get_coil().get_turns_description());
    // ABT #930: a wind with no turns names why -- the margin the winder holds for a section takes
    // the whole 1.85 mm width of the window.
    CHECK_THAT(result.get_coil().get_last_fit_failure(),
               ContainsSubstring("mm of margin across a winding window 1.850 mm wide") &&
               ContainsSubstring("leaves no room for one turn"));
}

TEST_CASE("Test_Wind_By_Layers_Throws_On_Overfull_Stated_Section", "[constructive-model][coil][stated-sections][smoke-test]") {
    SettingsReset reset;
    auto magneticJson = load_stated_sections_fixture("we_fc_et20_four_spaced_sections_overfull.json");
    OpenMagnetics::Coil coil(magneticJson["coil"], false);
    REQUIRE_THROWS_WITH(coil.wind_by_layers(),
                        ContainsSubstring("wind_by_layers") && ContainsSubstring("'L1 section 0'") &&
                        ContainsSubstring("needs 2.010 mm") && ContainsSubstring("has 1.800 mm"));
}

TEST_CASE("Test_Wind_By_Turns_Throws_On_Overfull_Stated_Layer", "[constructive-model][coil][stated-sections][smoke-test]") {
    SettingsReset reset;
    auto magneticJson = load_stated_sections_fixture("we_fc_et20_four_spaced_sections_overfull.json");
    OpenMagnetics::Coil coil(magneticJson["coil"], false);
    settings.set_coil_wind_even_if_not_fit(true);   // the caller's opt-in: lay the over-full layers out
    REQUIRE(coil.wind_by_layers());
    settings.set_coil_wind_even_if_not_fit(false);
    REQUIRE_THROWS_WITH(coil.wind_by_turns(),
                        ContainsSubstring("wind_by_turns") && ContainsSubstring("'L1 section 0'") &&
                        ContainsSubstring("needs 2.010 mm") && ContainsSubstring("has 1.800 mm"));
}

TEST_CASE("Test_Fitting_Check_On_Stated_Sections_Without_Filling_Factor", "[constructive-model][coil][stated-sections][smoke-test]") {
    SettingsReset reset;
    {
        auto magneticJson = load_stated_sections_fixture("we_fc_et20_four_spaced_sections.json");
        OpenMagnetics::Coil coil(magneticJson["coil"], false);
        REQUIRE(coil.wind_by_layers());
        REQUIRE(coil.wind_by_turns());
        REQUIRE_FALSE(coil.get_sections_description().value()[0].get_filling_factor());
        CHECK(coil.are_sections_and_layers_fitting());
        check_turns_inside_stated_sections(coil, magneticJson["coil"]["sectionsDescription"].get<std::vector<Section>>());
    }
    {
        auto magneticJson = load_stated_sections_fixture("we_fc_et20_four_spaced_sections_overfull.json");
        OpenMagnetics::Coil coil(magneticJson["coil"], false);
        settings.set_coil_wind_even_if_not_fit(true);
        REQUIRE(coil.wind_by_layers());
        REQUIRE(coil.wind_by_turns());
        CHECK_FALSE(coil.are_sections_and_layers_fitting());
    }
}

