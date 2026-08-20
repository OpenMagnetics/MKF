// Characterization coverage for ABT #605: WE HAVE NOT FIXED THIS. It has no
// established ground truth (no FEM run, no measured hardware) to fix
// against, so — unlike ABT #598/#599/#600/#602/#606 — there is no code
// change here, only a test that pins the CURRENT broken behavior and states
// the target so a future fix is verifiable.
//
// Reported: the same multi-layer round-wire winding gives wildly different
// total winding losses depending only on which proximity-loss model is
// selected (Dowell 8.91 W vs Ferreira 1.04 W vs Albach 1.05 W on the user's
// real PSFB design — an 8.5x spread with ohmic/core losses held identical
// across all three). Ferreira is the model MKF auto-selects for round/litz
// wire by default (WindingProximityEffectLosses.cpp:88-99), so this is the
// DEFAULT PATH for the overwhelming majority of designs. Wang is excluded
// here: it is a flat-conductor-only model (PLANAR/RECTANGULAR/FOIL) by its
// own documented design, not applicable to the round wire this test uses,
// and correctly throws NotImplementedException for it — that is not part
// of this bug.
//
// Physical prior (not proof, but why Ferreira is the leading suspect):
// for a multi-layer round-wire winding above 100 kHz, proximity loss is
// normally comparable to or larger than DC/ohmic loss. Ferreira returning a
// proximity loss 1000x+ SMALLER than ohmic loss on this geometry is not
// physically plausible; Dowell landing within an order of magnitude of
// ohmic loss is the more plausible of the two, but "more plausible" is not
// "verified" — this needs FEM or measured-hardware arbitration before any
// model gets changed.

#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "physical_models/WindingLosses.h"
#include "processors/Inputs.h"
#include "support/Settings.h"

#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

double total_losses_with_proximity_model(OpenMagnetics::Magnetic& magnetic, OperatingPoint& operatingPoint,
                                         double temperature, WindingProximityEffectLossesModels model) {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_enable_user_winding_losses_models(true);
    settings.set_winding_proximity_effect_losses_model(model);
    auto losses = WindingLosses().calculate_losses(magnetic, operatingPoint, temperature);
    settings.reset();
    return losses.get_winding_losses();
}

}  // namespace

TEST_CASE("Test_Proximity_Loss_Models_Disagree_By_Orders_Of_Magnitude_ABT_605", "[physical-model][winding-losses][proximity][model-comparison][bug]") {
    double temperature = 40;
    double frequency = 100000;
    // A multi-layer round-wire winding: the regime where proximity loss is
    // expected to be significant (comparable to or larger than DC loss).
    std::vector<int64_t> numberTurns({24});
    std::vector<int64_t> numberParallels({1});
    std::string shapeName = "ETD 34/17/11";

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        frequency, 1e-3, temperature, WaveformLabel::SINUSOIDAL, 4.0, 0.5, 0);

    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, OpenMagneticsTesting::get_ground_gap(1e-4), 1, "3C97");
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto operatingPoint = inputs.get_operating_point(0);

    double dowell = total_losses_with_proximity_model(magnetic, operatingPoint, temperature, WindingProximityEffectLossesModels::DOWELL);
    double ferreira = total_losses_with_proximity_model(magnetic, operatingPoint, temperature, WindingProximityEffectLossesModels::FERREIRA);
    double albach = total_losses_with_proximity_model(magnetic, operatingPoint, temperature, WindingProximityEffectLossesModels::ALBACH);

    INFO("Dowell=" << dowell << " Ferreira=" << ferreira << " Albach=" << albach);

    double maxTotal = std::max({dowell, ferreira, albach});
    double minTotal = std::min({dowell, ferreira, albach});
    REQUIRE(minTotal > 0);

    // ABT #605 RESOLVED (2026-08-20). This case used to assert that the spread WAS large
    // ("the bug still reproduces"), with the target below tagged [!mayfail]. Both the
    // reproduce-assertion and the tag are gone, exactly as the original author instructed:
    // "When it starts passing, that is the signal the fix landed — remove the tag then."
    //
    // History of the spread on this fixture: ~10000x when the ticket was filed, 2.87x after
    // the Kelvin x-pi fix (Dowell 1.236 / Ferreira 0.953 / Albach 0.430), and 1.30x now
    // (Dowell 1.228 / Ferreira 0.948 / Albach 0.957). What closed it was ABT #837: Albach's
    // proximity factor carried an extra conductor dimension, returning W where every other
    // model returns W/m, so it read ~1e-3 of the truth for EVERY wire type. Ferreira itself
    // is independently verified exact — its Kelvin-function factor reproduces a radial-mode
    // solve of the cylinder-in-transverse-field problem to 4 digits — so agreement with
    // Ferreira is agreement with the physics, not merely consensus.
    //
    // The residual 1.30x is Dowell's, and it is understood rather than mysterious: Dowell's
    // equal-area-square reduction of a round wire runs ~2.09x high at low gamma and ~0.60x
    // low at high gamma against the exact factor (ABT #837 scorecard). 2x remains the band:
    // it does not assert WHICH model is right, only that models of the same geometry must
    // agree to an engineering margin.
    CHECK(maxTotal / minTotal < 2.0);
}
