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
#include "Fixtures.h"

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

// A bobbin may legitimately carry no winding-window `area`: MAS marks it optional on
// windingWindowElement, which requires only width+height (rectangular) or
// angle+radialHeight (round). calculate_filling_factor used to read it as
// `windingWindows[0].get_area().value()`, so such a bobbin threw std::bad_optional_access
// — and because the WASM binding returned a bare what(), the browser JSON.parse'd the
// message and showed
//     SyntaxError: Unexpected token 'b', "bad_optional_access" is not valid JSON
// which names neither the bobbin nor the field.
//
// The fixture is the exact coil the MagneticBuilder sent for WE 744025006, a drumRing
// inductor whose catalogue MAS carries a processedDescription-only bobbin (no
// functionalDescription, so process_bobbin cannot fill the area in either — that is a
// legal bobbin, not a broken one). Its window is 0.3 mm x 2.08 mm with area and shape null.
TEST_CASE("Test_Filling_Factors_Bobbin_Without_Winding_Window_Area", "[coil][filling-factor][smoke-test]") {
    settings.reset();
    json coilJson = OpenMagneticsTesting::fixtures::get_json("coil-drumring-bobbin-without-winding-window-area");

    // Precondition: this is only a regression test while the fixture really omits the area.
    REQUIRE(coilJson["bobbin"]["processedDescription"]["windingWindows"][0]["area"].is_null());

    OpenMagnetics::Coil coil(coilJson, false);
    OpenMagnetics::Coil::FillingFactorsOutput factors;
    REQUIRE_NOTHROW(factors = coil.calculate_filling_factor());

    // The area is derived from the window's own dimensions with the same formula
    // create_quick_bobbin writes: 0.0003 m x 0.00208 m = 6.24e-7 m2. A wrong denominator
    // would still not throw, so pin the value rather than just the absence of a throw.
    const double expectedAvailableArea = 0.0002999999999999999 * 0.0020800000000000003;
    CHECK(factors.areaFillingFactor > 0);
    CHECK(factors.areaFillingFactor < 1);

    // Same coil, area stated explicitly: the derivation must agree with what MKF would
    // have written itself, so both paths give the same answer.
    json statedJson = coilJson;
    statedJson["bobbin"]["processedDescription"]["windingWindows"][0]["area"] = expectedAvailableArea;
    OpenMagnetics::Coil stated(statedJson, false);
    auto statedFactors = stated.calculate_filling_factor();
    CHECK_THAT(factors.areaFillingFactor,
               Catch::Matchers::WithinRel(statedFactors.areaFillingFactor, 1e-12));
}

// A drum's winding must end up ON the former, through the path a user actually travels.
//
// Alf's report was visual: on WE 74402500030, a drumRing whose 0.3 mm window carries one
// 0.191 mm layer, the turns rendered floating in the middle of the window with equal 0.053 mm
// gaps either side, when the drum post's surface is at 0.650 mm. A coil does not levitate off
// the bobbin it is wound on.
//
// The cause was in magnetic_autocomplete (ABT #998): the sections orientation was chosen from
// the CORE TYPE, so every non-two-piece core — drums, rods, molded, and the whole
// piece-and-plate family — was told to wind as a flat spiral. CONTIGUOUS means the turns
// advance RADIALLY, so MKF then took the section's radial position from the turns alignment
// (Coil.cpp:12817) and a CENTERED alignment put the copper mid-window. Every step after the
// orientation was behaving correctly, which is why compaction looked guilty for a while: it
// faithfully shrink-wrapped a section that was already in the wrong place.
//
// This asserts the END-TO-END property rather than any single step, because the intermediate
// contracts are exactly what was misunderstood: MKF is entitled to centre a genuine spiral,
// and a test saying otherwise would pin a contract MKF does not have. Feeding a catalogue
// magnetic that states NO orientation is the point — it exercises the choice the engine makes,
// which is where the defect lived.
TEST_CASE("Test_Drum_Winding_Sits_On_The_Former", "[coil][winding][smoke-test]") {
    settings.reset();
    json magneticJson = OpenMagneticsTesting::fixtures::get_json("magnetic-drumring-74402500030");

    // Precondition: the fixture must NOT state an orientation, or this tests the caller's
    // configuration instead of the engine's choice — and explicit configuration wins.
    REQUIRE(magneticJson["coil"]["bobbin"]["processedDescription"]["windingWindows"][0]
                .contains("sectionsOrientation") == false);

    // Through magnetic_autocomplete, NOT a bare coil.wind(): the orientation is chosen there,
    // so a test that winds the coil directly bypasses the defect entirely. Confirmed by
    // negative control — winding directly passes even with ABT #998 reverted, which is exactly
    // the "green either way" test that proves nothing.
    OpenMagnetics::Magnetic magnetic = OpenMagnetics::magnetic_autocomplete(
        OpenMagnetics::Magnetic(magneticJson));
    auto coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_turns_description());

    // By VALUE: get_turns_description() returns by value, so a reference bound through it
    // dangles at the end of the full expression (ABT #964, ABT #988).
    auto turns = coil.get_turns_description().value();
    double innermost = std::numeric_limits<double>::max();
    for (const auto& turn : turns) {
        innermost = std::min(innermost, turn.get_coordinates()[0] - turn.get_dimensions().value()[0] / 2);
    }

    // The drum post's surface: bobbin columnWidth, which for a zero-thickness former is the
    // core column half-width. The wire lies on it.
    const double formerSurface = 0.00065;
    CHECK_THAT(innermost, Catch::Matchers::WithinAbs(formerSurface, 1e-6));
}
