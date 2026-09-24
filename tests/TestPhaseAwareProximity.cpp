// Phase-aware proximity field (MAS excitation convention, 2026-09-24).
//
// Convention: every winding's voltage is in the common dot reference; primary-side windings
// (isolationSide "primary") are PASSIVE (+ into the dot), every other winding SOURCES (+ out of
// the dot). The physical dot-referenced current of winding k is c_k i_k, c_k = +1 primary side,
// -1 otherwise, and the inducing field of a harmonic is the COMPLEX sum of c_k I_k,h (phasor
// from the waveform) times the turn's unit-current field.
//
// These tests use only API that also exists on the amplitude-only code (origin/main before the
// change), so the same file demonstrates the old behaviour failing them.
#include <source_location>
#include "physical_models/MagneticField.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingProximityEffectLosses.h"
#include "processors/Inputs.h"
#include "support/Settings.h"
#include "json.hpp"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

// Golden values of the amplitude-only model, produced by origin/main 38af34a3 (printed by the
// tests below with %.17g). Single winding and exact antiphase must not move.
#define GOLDEN_SINGLE_ALBACH 0.071061297472233192
#define GOLDEN_SINGLE_BINNS 0.09268152574623513
#define GOLDEN_ANTIPHASE_ALBACH 0.18785344918292851
#define GOLDEN_ANTIPHASE_BINNS 0.20013109596372328
#define GOLDEN_ANTIPHASE_PIPELINE 0.37276841004685124

namespace {

const double frequency = 200000;
const double temperature = 25;

OperatingPointExcitation make_current_excitation(const std::vector<double>& time, const std::vector<double>& data) {
    Waveform waveform;
    waveform.set_data(data);
    waveform.set_time(time);
    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform, frequency);
    auto harmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledWaveform, frequency);
    SignalDescriptor current;
    current.set_waveform(waveform);
    current.set_harmonics(harmonics);
    current.set_processed(OpenMagnetics::Inputs::calculate_processed_data(harmonics, sampledWaveform, true));
    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    excitation.set_current(current);
    return excitation;
}

// peak * sin(w t + phase), sampled over one period (endpoint included).
OperatingPointExcitation sinusoidal_excitation(double peak, double phase) {
    const size_t numberPoints = 256;
    std::vector<double> time;
    std::vector<double> data;
    for (size_t i = 0; i <= numberPoints; ++i) {
        double t = static_cast<double>(i) / static_cast<double>(numberPoints) / frequency;
        time.push_back(t);
        data.push_back(peak * std::sin(2 * std::numbers::pi * frequency * t + phase));
    }
    return make_current_excitation(time, data);
}

// A current pulse train on a uniform grid of 256 intervals: `level` on samples
// [startSample, startSample + widthSamples) (modulo the period), zero elsewhere, linear between
// grid points. Defined on the grid so that a train started half a period later is an EXACT shift.
OperatingPointExcitation pulse_excitation(double level, size_t startSample, size_t widthSamples) {
    const size_t numberPoints = 256;
    std::vector<double> time;
    std::vector<double> data;
    for (size_t i = 0; i <= numberPoints; ++i) {
        time.push_back(static_cast<double>(i) / static_cast<double>(numberPoints) / frequency);
        size_t offset = (i % numberPoints + numberPoints - startSample) % numberPoints;
        data.push_back(offset < widthSamples ? level : 0.0);
    }
    return make_current_excitation(time, data);
}

OperatingPoint make_operating_point(std::vector<OperatingPointExcitation> excitations) {
    OperatingPoint operatingPoint;
    OperatingConditions conditions;
    conditions.set_ambient_temperature(temperature);
    operatingPoint.set_conditions(conditions);
    operatingPoint.set_excitations_per_winding(excitations);
    return operatingPoint;
}

OpenMagnetics::Magnetic make_magnetic(std::vector<int64_t> numberTurns, std::vector<std::string> isolationSides) {
    std::string shapeName = "ETD 34";
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, json::array(), 1, "3C97");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    json coilJson;
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    for (size_t i = 0; i < numberTurns.size(); ++i) {
        json windingJson;
        windingJson["name"] = "winding " + std::to_string(i);
        windingJson["numberTurns"] = numberTurns[i];
        windingJson["numberParallels"] = 1;
        windingJson["isolationSide"] = isolationSides[i];
        windingJson["wire"] = "Round 0.475 - Grade 1";
        coilJson["functionalDescription"].push_back(windingJson);
    }
    OpenMagnetics::Coil coil(coilJson, 1, WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                             CoilAlignment::CENTERED, CoilAlignment::CENTERED);
    coil.delimit_and_compact();
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    return magnetic;
}

double total_proximity_losses(const WindingLossesOutput& windingLossesOutput) {
    double total = 0;
    auto windingLossesPerTurn = windingLossesOutput.get_winding_losses_per_turn().value();
    for (auto& turnLosses : windingLossesPerTurn) {
        auto proximity = turnLosses.get_proximity_effect_losses();
        if (!proximity) {
            continue;
        }
        for (auto loss : proximity->get_losses_per_harmonic()) {
            total += loss;
        }
    }
    return total;
}

// Proximity losses of a given (real) field output: the loss the field would cause if it were
// the whole field of the harmonic.
double proximity_of_field(OpenMagnetics::Magnetic& magnetic, const OperatingPoint& operatingPoint, const std::vector<ComplexField>& fieldPerFrequency) {
    WindingWindowMagneticStrengthFieldOutput output;
    output.set_field_per_frequency(fieldPerFrequency);
    output.set_method_used("test");
    output.set_origin(ResultOrigin::SIMULATION);
    auto ohmic = WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, temperature);
    auto losses = WindingProximityEffectLosses::calculate_proximity_effect_losses(magnetic.get_coil(), temperature, ohmic, output);
    return total_proximity_losses(losses);
}

// Proximity losses as the production path computes them: the field model's default directions
// and phases, then the proximity models.
double proximity_default(OpenMagnetics::Magnetic& magnetic, const OperatingPoint& operatingPoint, MagneticFieldStrengthModels model) {
    MagneticField magneticField(model);
    auto fieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic);
    auto ohmic = WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, temperature);
    auto losses = WindingProximityEffectLosses::calculate_proximity_effect_losses(magnetic.get_coil(), temperature, ohmic, fieldOutput);
    return total_proximity_losses(losses);
}

// In-phase (real) field of the given windings only, with explicit signs.
std::vector<ComplexField> field_with_directions(OpenMagnetics::Magnetic& magnetic, const OperatingPoint& operatingPoint, MagneticFieldStrengthModels model, std::vector<int8_t> directions) {
    MagneticField magneticField(model);
    return magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, std::nullopt, directions).get_field_per_frequency();
}

std::vector<ComplexField> linear_combination(const std::vector<ComplexField>& a, double weightA, const std::vector<ComplexField>& b, double weightB) {
    REQUIRE(a.size() == b.size());
    auto result = a;
    for (size_t h = 0; h < a.size(); ++h) {
        REQUIRE(a[h].get_data().size() == b[h].get_data().size());
        for (size_t p = 0; p < a[h].get_data().size(); ++p) {
            auto& point = result[h].get_mutable_data()[p];
            point.set_real(weightA * a[h].get_data()[p].get_real() + weightB * b[h].get_data()[p].get_real());
            point.set_imaginary(weightA * a[h].get_data()[p].get_imaginary() + weightB * b[h].get_data()[p].get_imaginary());
        }
    }
    return result;
}

const std::vector<std::pair<MagneticFieldStrengthModels, std::string>> fieldModels = {
    {MagneticFieldStrengthModels::ALBACH, "ALBACH"},
    {MagneticFieldStrengthModels::BINNS_LAWRENSON, "BINNS_LAWRENSON"},
};

void prepare_settings() {
    settings.reset();
    settings.set_magnetic_field_include_fringing(false);
}

}  // namespace

TEST_CASE("Test_Phase_Aware_Single_Winding_Unchanged", "[physical-model][winding-losses][phase-aware]") {
    // Regression anchor: one winding. The golden values were produced by the amplitude-only
    // model (origin/main 38af34a3) and must be reproduced bit for bit (1e-12 relative).
    prepare_settings();
    auto magnetic = make_magnetic({12}, {"primary"});
    auto operatingPoint = make_operating_point({sinusoidal_excitation(2.0, 0.3)});
    const std::vector<double> golden = {GOLDEN_SINGLE_ALBACH, GOLDEN_SINGLE_BINNS};
    for (size_t m = 0; m < fieldModels.size(); ++m) {
        double losses = proximity_default(magnetic, operatingPoint, fieldModels[m].first);
        std::printf("[phase-aware] single %s %.17g\n", fieldModels[m].second.c_str(), losses);
        INFO(fieldModels[m].second);
        REQUIRE(losses > 0);
        CHECK_THAT(losses, Catch::Matchers::WithinRel(golden[m], 1e-12));
    }
    settings.reset();
}

TEST_CASE("Test_Phase_Aware_Antiphase_Transformer_Unchanged", "[physical-model][winding-losses][phase-aware]") {
    // Regression anchor: primary + secondary whose MAS currents are the SAME sinusoid. In the
    // source convention the secondary is + out of the dot, so physically the two MMFs oppose —
    // the exact antiphase the amplitude-only model assumed. Golden values from origin/main
    // 38af34a3, reproduced to 1e-12 relative. The whole WindingLosses pipeline is checked too.
    prepare_settings();
    auto magnetic = make_magnetic({8, 8}, {"primary", "secondary"});
    auto operatingPoint = make_operating_point({sinusoidal_excitation(2.0, 0.3), sinusoidal_excitation(2.0, 0.3)});
    const std::vector<double> golden = {GOLDEN_ANTIPHASE_ALBACH, GOLDEN_ANTIPHASE_BINNS};
    for (size_t m = 0; m < fieldModels.size(); ++m) {
        double losses = proximity_default(magnetic, operatingPoint, fieldModels[m].first);
        std::printf("[phase-aware] antiphase %s %.17g\n", fieldModels[m].second.c_str(), losses);
        INFO(fieldModels[m].second);
        REQUIRE(losses > 0);
        CHECK_THAT(losses, Catch::Matchers::WithinRel(golden[m], 1e-12));
    }

    WindingLosses windingLosses;
    auto pipeline = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
    std::printf("[phase-aware] antiphase pipeline %.17g\n", pipeline.get_winding_losses());
    CHECK_THAT(pipeline.get_winding_losses(), Catch::Matchers::WithinRel(GOLDEN_ANTIPHASE_PIPELINE, 1e-12));
    settings.reset();
}

TEST_CASE("Test_Phase_Aware_Shifted_Secondary_Superposition", "[physical-model][winding-losses][phase-aware]") {
    // The secondary current shifted by 45 and 90 degrees. The loss must change, and must equal
    // an independent evaluation: build H1 (primary alone) and H2 (secondary alone, unshifted,
    // source sign) from single-winding fields, form H = H1 + H2 exp(j phi) by hand, and take the
    // loss of |H|^2 = |Re H|^2 + |Im H|^2 as the loss of the real field Re H plus that of Im H
    // (the proximity models here are quadratic forms of the field).
    prepare_settings();
    auto magnetic = make_magnetic({8, 8}, {"primary", "secondary"});
    auto unshifted = make_operating_point({sinusoidal_excitation(2.0, 0.3), sinusoidal_excitation(2.0, 0.3)});

    for (auto& [model, modelName] : fieldModels) {
        auto primaryAlone = field_with_directions(magnetic, unshifted, model, {1, 0});
        auto secondaryAlone = field_with_directions(magnetic, unshifted, model, {0, -1});
        double antiphaseLosses = proximity_default(magnetic, unshifted, model);
        // Sanity of the construction itself: at phi = 0 the hand-built field is the model's.
        double superposedAtZero = proximity_of_field(magnetic, unshifted, linear_combination(primaryAlone, 1, secondaryAlone, 1));
        INFO(modelName);
        CHECK_THAT(superposedAtZero, Catch::Matchers::WithinRel(antiphaseLosses, 1e-9));

        for (double shiftDegrees : {45.0, 90.0}) {
            double shift = shiftDegrees * std::numbers::pi / 180;
            auto shifted = make_operating_point({sinusoidal_excitation(2.0, 0.3), sinusoidal_excitation(2.0, 0.3 + shift)});
            double losses = proximity_default(magnetic, shifted, model);

            auto inPhase = linear_combination(primaryAlone, 1, secondaryAlone, std::cos(shift));
            auto quadrature = linear_combination(primaryAlone, 0, secondaryAlone, std::sin(shift));
            double expected = proximity_of_field(magnetic, shifted, inPhase) + proximity_of_field(magnetic, shifted, quadrature);
            std::printf("[phase-aware] %s shift %.0f deg: model %.10g superposition %.10g antiphase %.10g\n",
                        modelName.c_str(), shiftDegrees, losses, expected, antiphaseLosses);
            INFO("shift " << shiftDegrees << " deg");
            CHECK_THAT(losses, Catch::Matchers::WithinRel(expected, 1e-9));
            // Less cancellation than exact antiphase: the loss must rise noticeably.
            CHECK(losses > 1.05 * antiphaseLosses);
        }
    }
    settings.reset();
}

TEST_CASE("Test_Phase_Aware_Push_Pull_Primary_Halves_Add", "[physical-model][winding-losses][phase-aware]") {
    // Push-pull, 4 windings: Primary Half 1/2 on the primary side, Secondary Half 1/2 on the
    // secondary side. In the dot-referenced passive convention Primary Half 1 conducts +I in the
    // first half-period and Primary Half 2 conducts -I (its energy flows in while the flux
    // ramps down) in the second: the half-period-shifted pulse trains below. At the fundamental
    // -pulse(t - T/2) has the SAME phasor as pulse(t), so both halves are +1 (primary side) with
    // equal phase and their fields ADD. The index rule (+1, -1, -1, -1) with amplitudes only made
    // them cancel. Secondaries carry no current here to isolate the primary halves. The second
    // train is an exact half-period shift of the first (see pulse_excitation).
    prepare_settings();
    auto magnetic = make_magnetic({6, 6, 3, 3}, {"primary", "primary", "secondary", "secondary"});
    auto operatingPoint = make_operating_point({pulse_excitation(3.0, 0, 96), pulse_excitation(-3.0, 128, 96),
                                                sinusoidal_excitation(0.0, 0.0), sinusoidal_excitation(0.0, 0.0)});

    for (auto& [model, modelName] : fieldModels) {
        // Independent reference: each half alone with its own physical sign (+1, passive,
        // primary side); the fields of the two halves summed as they physically are.
        auto half1 = field_with_directions(magnetic, operatingPoint, model, {1, 0, 0, 0});
        auto half2 = field_with_directions(magnetic, operatingPoint, model, {0, 1, 0, 0});
        // Every harmonic of the two pulse trains: odd harmonics add, even ones cancel; the
        // single-winding fields here are amplitude-signed, so compare at the fundamental only.
        std::vector<ComplexField> half1Fundamental = {half1[0]};
        std::vector<ComplexField> half2Fundamental = {half2[0]};
        double adding = proximity_of_field(magnetic, operatingPoint, linear_combination(half1Fundamental, 1, half2Fundamental, 1));
        double cancelling = proximity_of_field(magnetic, operatingPoint, linear_combination(half1Fundamental, 1, half2Fundamental, -1));

        MagneticField magneticField(model);
        auto fieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic);
        std::vector<ComplexField> modelFundamental = {fieldOutput.get_field_per_frequency()[0]};
        REQUIRE(modelFundamental[0].get_frequency() == frequency);
        double modelLosses = proximity_of_field(magnetic, operatingPoint, modelFundamental);
        std::printf("[phase-aware] push-pull %s fundamental: model %.10g adding %.10g cancelling %.10g\n",
                    modelName.c_str(), modelLosses, adding, cancelling);
        INFO(modelName);
        REQUIRE(adding > 2 * cancelling);
        CHECK_THAT(modelLosses, Catch::Matchers::WithinRel(adding, 1e-9));
    }
    settings.reset();
}
