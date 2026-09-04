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

// A winding narrower than its window must sit ON the former, not float in the middle of it.
//
// WE 74402500030 is a drumRing whose 0.3 mm window carries one 0.191 mm layer. Wound, the
// turns land against the drum post at 0.650 mm — correct. Run delimit_and_compact over the
// result and the layer MOVES OUT to the window centre, 0.7045..0.8955, with equal 0.053 mm
// gaps on both sides. A coil does not levitate off the bobbin it is wound on.
//
// The section orientation is what exposes it: with OVERLAPPING the section starts at the
// inner edge and stays there; with CONTIGUOUS the section is created spanning the FULL window
// width and the winding inside it ends up centred.
//
// THIS TEST FAILS ON PURPOSE — it is a reproducer for an unfixed defect, not a regression
// guard, and it is committed failing so the signal is not lost. Traced this far:
//
//   wind_by_sections  creates the section full-window-width, centred      (0.800 / 0.300)
//   wind_by_layers    places the layer against the former                 (0.7455 / 0.191)
//   delimit_and_compact  shrink-wraps the section ONTO the layer          (0.7455 / 0.191)  <- correct
//   ...a later pass inside wind_inner moves BOTH to the window centre     (0.800  / 0.191)  <- the bug
//
// The last step is the one still unattributed. wind_inner runs wind_by_turns +
// delimit_and_compact up to four times (the real-winding re-wind paths at 4887/4922/4949/5035),
// and the layer is already at 0.800 by the second compaction, so a re-wind is re-placing it
// centred rather than against the former. Instrumenting the two compaction calls shows the
// first behaving correctly, which is what rules compaction itself out.
//
// Not patched speculatively: this is a physics path, and a guessed fix here would produce a
// plausible wrong number rather than an error — the failure mode this whole area keeps
// exhibiting.
TEST_CASE("Test_Compaction_Keeps_A_Narrow_Winding_On_The_Former", "[coil][compaction][smoke-test]") {
    settings.reset();
    json coilJson = OpenMagneticsTesting::fixtures::get_json("coil-drumring-contiguous-narrow-section");

    auto windAndReport = [&](bool compact) {
        settings.reset();
        settings.set_coil_delimit_and_compact(compact);
        OpenMagnetics::Coil coil(coilJson, false);
        coil.wind({1.0}, {0}, 1);
        if (compact) {
            coil.delimit_and_compact();
        }
        REQUIRE(coil.get_turns_description());
        double innermost = std::numeric_limits<double>::max();
        // By VALUE: get_turns_description() returns by value, so a reference bound through
        // it dangles the moment the full expression ends (ABT #964, ABT #988).
        auto turns = coil.get_turns_description().value();
        for (const auto& turn : turns) {
            innermost = std::min(innermost, turn.get_coordinates()[0] - turn.get_dimensions().value()[0] / 2);
        }
        return innermost;
    };

    const double formerSurface = 0.00065;   // bobbin columnWidth: the post the wire lies on

    // Without compaction the winding is already correct, which is what makes compaction the
    // suspect rather than the winder.
    CHECK_THAT(windAndReport(false), Catch::Matchers::WithinAbs(formerSurface, 1e-6));

    // With it, the winding must not have moved off the former.
    CHECK_THAT(windAndReport(true), Catch::Matchers::WithinAbs(formerSurface, 1e-6));
}
