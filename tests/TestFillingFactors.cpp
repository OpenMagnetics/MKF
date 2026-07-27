// ABT #245: areaFillingFactor must be the area fraction and nothing else.
//
// It used to be max(largest per-layer factor, wound area / available area), so a coil
// with one section that could not hold its turns reported 43264% "area fill" for a
// winding whose real areal fill was 2.35%, and the builder printed that as a percentage
// next to a "winding does not fit" banner. The overfill is now reported separately as
// maxLayerFillingFactor and the verdict as windingFits.
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <fstream>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

namespace {
constexpr auto kCorruptToroidFixture =
    "/home/alf/OpenMagnetics/WebFrontend/tests/fixtures/toroidal_stale_pin_corrupt_t402416.json";
}

TEST_CASE("Test_Filling_Factors_Area_Fraction_Is_Not_An_Overfill_Ratio", "[coil][filling-factor][smoke-test]") {
    settings.reset();
    std::ifstream file(kCorruptToroidFixture);
    if (!file.good()) {
        SKIP("WebFrontend fixture not available in this checkout");
    }
    json masJson = json::parse(file);
    json magneticJson = masJson.contains("magnetic") ? masJson["magnetic"] : masJson;

    // As loaded: T 40/24/16 x4 stacks, 21x2 + 21x2, with a corrupted wind whose
    // secondary section spans 0.0048 deg and cannot hold its turns.
    OpenMagnetics::Magnetic corrupted(magneticJson);
    auto corruptedFactors = corrupted.get_mutable_coil().calculate_filling_factor();

    // The verdict is preserved...
    CHECK_FALSE(corruptedFactors.windingFits);
    CHECK(corruptedFactors.maxLayerFillingFactor > 1);
    // ...but the area fraction is an area fraction: this winding is a few percent of
    // the 452 mm2 bore, not several thousand.
    CHECK(corruptedFactors.areaFillingFactor > 0);
    CHECK(corruptedFactors.areaFillingFactor < 1);

    // A clean re-wind of the same magnetic fits and reports single-digit percentages.
    json cleanJson = magneticJson;
    cleanJson["coil"].erase("turnsDescription");
    cleanJson["coil"].erase("layersDescription");
    cleanJson["coil"].erase("sectionsDescription");
    cleanJson["coil"].erase("groupsDescription");
    OpenMagnetics::Magnetic clean(cleanJson);
    auto cleanCoil = clean.get_mutable_coil();
    cleanCoil.wind();
    auto cleanFactors = cleanCoil.calculate_filling_factor();

    CHECK(cleanFactors.windingFits);
    CHECK(cleanFactors.maxLayerFillingFactor <= 1);
    CHECK_THAT(cleanFactors.areaFillingFactor, Catch::Matchers::WithinAbs(0.0235, 0.005));
    CHECK(cleanFactors.contiguousFillingFactor > 0);
    CHECK(cleanFactors.overlappingFillingFactor > 0);
}
