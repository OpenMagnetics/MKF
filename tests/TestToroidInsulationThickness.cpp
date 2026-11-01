// The thickness of a toroid's insulation, read back from a coil the way the web reads it: a Coil
// built from the stored JSON without winding (WebLibMKF get_insulation_section_thickness /
// get_insulation_layer_thickness).
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <source_location>
#include <fstream>
#include "constructive_models/Coil.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "json.hpp"

using json = nlohmann::json;
using Catch::Matchers::WithinRel;
using Catch::Matchers::ContainsSubstring;

// The web CMC wizard's Path A coil as the coil adviser returned it (T 10.5/5.5/20 in R10KZ, 2 x 6
// turns of Litz 135x0.071, contiguous sectors). Its two walls between the sectors are stored as
// polar sections and layers of [2.75 mm radial height, 0.520873 deg]: one 25 um tape, the angle
// being wound_distance_to_angle(25 um, 2.75 mm). The web turned that angle back into a thickness
// itself, with radialHeight * sin(angle) and the angle in degrees taken as radians, which gave
// 1.3685 mm. It then wound walls 55 times thicker than the advised ones, and the winder rightly
// refused the coil: "6 turns of a 1.2695 mm wire need 3 rings, 3.8085 mm deep, in a window of
// 2.75 mm radial height". The thickness comes from MKF, the inverse of the angle MKF made.
TEST_CASE("Test_Toroid_Insulation_Wall_Thickness_Read_Back_From_A_Stored_Coil", "[constructive-model][coil][toroidal][insulation][bug]") {
    OpenMagnetics::settings.reset();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "web_adviser_flows/cmc_path_a_r10kz_advised_coil.json");
    std::ifstream file(path);
    REQUIRE(file.good());
    OpenMagnetics::Coil coil(json::parse(file), false);

    const double tapeThickness = 25e-6;
    for (const auto& sectionName : {"Insulation between Line and Neutral section 1", "Insulation between Neutral and Line section 3"}) {
        INFO(sectionName);
        CHECK_THAT(coil.get_insulation_section_thickness(sectionName), WithinRel(tapeThickness, 1e-9));
    }
    for (const auto& layerName : {"Insulation between Line and Neutral section 1 insulation layer 0", "Insulation between Neutral and Line section 3 insulation layer 0"}) {
        INFO(layerName);
        CHECK_THAT(coil.get_insulation_layer_thickness(std::string(layerName)), WithinRel(tapeThickness, 1e-9));
    }

    // A name that is not an insulation section has no thickness. It used to sum no layers and
    // answer 0 m, which the web would have wound as no wall at all.
    CHECK_THROWS_WITH(coil.get_insulation_section_thickness("Insulation between Line and Neutral section 9"), ContainsSubstring("Not found section"));
    CHECK_THROWS_WITH(coil.get_insulation_section_thickness("Line section 0"), ContainsSubstring("is not an insulation section"));
}
