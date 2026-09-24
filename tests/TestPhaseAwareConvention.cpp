// MAS excitation convention (2026-09-24): the direction and phase primitives behind the
// phase-aware proximity field (see TestPhaseAwareProximity.cpp for the loss-level checks).
#include "physical_models/MagneticField.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingProximityEffectLosses.h"
#include "processors/Inputs.h"
#include "support/CoilMesher.h"
#include "support/Settings.h"
#include "json.hpp"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

const double frequency = 100000;

// Sum of sinusoids peak_i * sin(w t + phase_i), sampled over one period (endpoint included).
SignalDescriptor sinusoid_signal(std::vector<std::pair<double, double>> terms) {
    const size_t numberPoints = 256;
    std::vector<double> time;
    std::vector<double> data;
    for (size_t i = 0; i <= numberPoints; ++i) {
        double t = static_cast<double>(i) / static_cast<double>(numberPoints) / frequency;
        time.push_back(t);
        double value = 0;
        for (auto [peak, phase] : terms) {
            value += peak * std::sin(2 * std::numbers::pi * frequency * t + phase);
        }
        data.push_back(value);
    }
    Waveform waveform;
    waveform.set_data(data);
    waveform.set_time(time);
    auto sampledWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform, frequency);
    auto harmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledWaveform, frequency);
    SignalDescriptor signal;
    signal.set_waveform(waveform);
    signal.set_harmonics(harmonics);
    signal.set_processed(OpenMagnetics::Inputs::calculate_processed_data(harmonics, sampledWaveform, true));
    return signal;
}

OperatingPointExcitation sinusoid(double peak, double phase) {
    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    excitation.set_current(sinusoid_signal({{peak, phase}}));
    return excitation;
}

OpenMagnetics::Magnetic make_magnetic(std::vector<int64_t> numberTurns, std::vector<std::string> isolationSides, json gapping = json::array()) {
    auto core = OpenMagneticsTesting::get_quick_core("ETD 34", gapping, 1, "3C97");
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

OperatingPoint make_operating_point(std::vector<OperatingPointExcitation> excitations) {
    OperatingPoint operatingPoint;
    OperatingConditions conditions;
    conditions.set_ambient_temperature(25);
    operatingPoint.set_conditions(conditions);
    operatingPoint.set_excitations_per_winding(excitations);
    return operatingPoint;
}

}  // namespace

TEST_CASE("Test_Phase_Aware_Direction_By_Isolation_Side", "[physical-model][phase-aware]") {
    auto magnetic = make_magnetic({4, 4, 2, 2}, {"primary", "primary", "secondary", "tertiary"});
    auto directions = CoilMesher::calculate_current_direction_per_winding(magnetic.get_coil());
    REQUIRE(directions == std::vector<int8_t>{1, 1, -1, -1});
    CHECK(CoilMesher::get_reference_winding_index(magnetic.get_coil()) == 0);

    // The reference winding is the FIRST primary-side one, wherever it sits.
    auto secondaryFirst = make_magnetic({2, 4}, {"secondary", "primary"});
    CHECK(CoilMesher::calculate_current_direction_per_winding(secondaryFirst.get_coil()) == std::vector<int8_t>{-1, 1});
    CHECK(CoilMesher::get_reference_winding_index(secondaryFirst.get_coil()) == 1);

    auto noPrimary = make_magnetic({2, 4}, {"secondary", "tertiary"});
    CHECK_THROWS(CoilMesher::get_reference_winding_index(noPrimary.get_coil()));
}

TEST_CASE("Test_Phase_Aware_Phase_From_Waveform", "[physical-model][phase-aware]") {
    auto magnetic = make_magnetic({4, 4}, {"primary", "secondary"});
    for (double shiftDegrees : {0.0, 25.0, 90.0, 135.0, -60.0}) {
        double shift = shiftDegrees * std::numbers::pi / 180;
        auto operatingPoint = make_operating_point({sinusoid(1.0, 0.7), sinusoid(3.0, 0.7 + shift)});
        auto phases = CoilMesher::calculate_current_phase_per_winding(magnetic.get_coil(), operatingPoint, {1});
        REQUIRE(phases.size() == 1);
        INFO("shift " << shiftDegrees);
        // The reference winding is the gauge: exactly zero.
        CHECK(phases[0][0] == 0.0);
        CHECK_THAT(std::remainder(phases[0][1] - shift, 2 * std::numbers::pi), Catch::Matchers::WithinAbs(0, 1e-9));
    }

    // The reference winding carries none of a harmonic: the gauge moves to the first winding
    // that does, and a zero-amplitude winding keeps phase 0.
    auto operatingPoint = make_operating_point({sinusoid(0.0, 0.0), sinusoid(3.0, 1.1)});
    auto phases = CoilMesher::calculate_current_phase_per_winding(magnetic.get_coil(), operatingPoint, {1});
    CHECK(phases[0][0] == 0.0);
    CHECK(phases[0][1] == 0.0);

    // An excited winding without a waveform has no phase: throw, never assume one.
    auto noWaveform = make_operating_point({sinusoid(1.0, 0.0), sinusoid(3.0, 0.2)});
    auto current = noWaveform.get_excitations_per_winding()[1].get_current().value();
    current.set_waveform(std::nullopt);
    noWaveform.get_mutable_excitations_per_winding()[1].set_current(current);
    CHECK_THROWS(CoilMesher::calculate_current_phase_per_winding(magnetic.get_coil(), noWaveform, {1}));
}

TEST_CASE("Test_Phase_Aware_Quadrature_Field", "[physical-model][magnetic-field][phase-aware]") {
    settings.reset();
    settings.set_magnetic_field_include_fringing(false);
    auto magnetic = make_magnetic({8, 8}, {"primary", "secondary"});
    for (auto model : {MagneticFieldStrengthModels::ALBACH, MagneticFieldStrengthModels::BINNS_LAWRENSON}) {
        MagneticField magneticField(model);
        // In phase (source convention -> physical antiphase): no quadrature field at all.
        auto antiphase = magneticField.calculate_magnetic_field_strength_field(make_operating_point({sinusoid(2.0, 0.3), sinusoid(2.0, 0.3)}), magnetic);
        REQUIRE(antiphase.get_quadrature_field_per_frequency().size() == antiphase.get_field_per_frequency().size());
        double inPhaseSquared = 0;
        double quadratureSquared = 0;
        for (size_t p = 0; p < antiphase.get_field_per_frequency()[0].get_data().size(); ++p) {
            auto& inPhase = antiphase.get_field_per_frequency()[0].get_data()[p];
            auto& quadrature = antiphase.get_quadrature_field_per_frequency()[0].get_data()[p];
            CHECK(inPhase.get_point() == quadrature.get_point());
            CHECK(inPhase.get_turn_index() == quadrature.get_turn_index());
            inPhaseSquared += inPhase.get_real() * inPhase.get_real() + inPhase.get_imaginary() * inPhase.get_imaginary();
            quadratureSquared += quadrature.get_real() * quadrature.get_real() + quadrature.get_imaginary() * quadrature.get_imaginary();
        }
        REQUIRE(inPhaseSquared > 0);
        CHECK(quadratureSquared <= 1e-24 * inPhaseSquared);

        // 90 degrees: the secondary's field moves entirely into the quadrature component, so the
        // quadrature field equals the secondary-alone field (up to sign).
        auto quadratureCase = magneticField.calculate_magnetic_field_strength_field(make_operating_point({sinusoid(2.0, 0.3), sinusoid(2.0, 0.3 + std::numbers::pi / 2)}), magnetic);
        auto secondaryAlone = magneticField.calculate_magnetic_field_strength_field(make_operating_point({sinusoid(2.0, 0.3), sinusoid(2.0, 0.3)}), magnetic, std::nullopt, std::vector<int8_t>{0, -1});
        auto primaryAlone = magneticField.calculate_magnetic_field_strength_field(make_operating_point({sinusoid(2.0, 0.3), sinusoid(2.0, 0.3)}), magnetic, std::nullopt, std::vector<int8_t>{1, 0});
        const auto& quadratureData = quadratureCase.get_quadrature_field_per_frequency()[0].get_data();
        const auto& inPhaseData = quadratureCase.get_field_per_frequency()[0].get_data();
        const auto& secondaryData = secondaryAlone.get_field_per_frequency()[0].get_data();
        const auto& primaryData = primaryAlone.get_field_per_frequency()[0].get_data();
        REQUIRE(quadratureData.size() == secondaryData.size());
        double scale = 0;
        for (auto& point : secondaryData) {
            scale = std::max(scale, std::hypot(point.get_real(), point.get_imaginary()));
        }
        REQUIRE(scale > 0);
        for (size_t p = 0; p < quadratureData.size(); ++p) {
            CHECK_THAT(std::abs(quadratureData[p].get_real()), Catch::Matchers::WithinAbs(std::abs(secondaryData[p].get_real()), 1e-9 * scale));
            CHECK_THAT(std::abs(quadratureData[p].get_imaginary()), Catch::Matchers::WithinAbs(std::abs(secondaryData[p].get_imaginary()), 1e-9 * scale));
            CHECK_THAT(inPhaseData[p].get_real(), Catch::Matchers::WithinAbs(primaryData[p].get_real(), 1e-9 * scale));
            CHECK_THAT(inPhaseData[p].get_imaginary(), Catch::Matchers::WithinAbs(primaryData[p].get_imaginary(), 1e-9 * scale));
        }
    }
    settings.reset();
}

namespace {

double proximity_of(OpenMagnetics::Magnetic& magnetic, const OperatingPoint& operatingPoint, const std::vector<ComplexField>& inPhase, const std::vector<ComplexField>& quadrature) {
    auto ohmic = WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, 25);
    WindingWindowMagneticStrengthFieldPhasorOutput output;
    output.set_field_per_frequency(inPhase);
    output.set_quadrature_field_per_frequency(quadrature);
    output.set_method_used("test");
    output.set_origin(ResultOrigin::SIMULATION);
    auto losses = WindingProximityEffectLosses::calculate_proximity_effect_losses(magnetic.get_coil(), 25, ohmic, output);
    double total = 0;
    auto perTurn = losses.get_winding_losses_per_turn().value();
    for (auto& turn : perTurn) {
        auto proximity = turn.get_proximity_effect_losses();
        if (proximity) {
            for (auto loss : proximity->get_losses_per_harmonic()) {
                total += loss;
            }
        }
    }
    return total;
}

// A gapped 8:8 transformer: primary = load + magnetizing current, secondary = load (source
// convention), and the magnetizing current supplied explicitly on the primary so the gap field
// magnitude does not depend on the winding currents.
OperatingPoint gapped_operating_point(double loadPhase, double magnetizingPeak, double magnetizingPhase) {
    auto primary = sinusoid(0, 0);
    primary.set_current(sinusoid_signal({{2.0, loadPhase}, {magnetizingPeak, magnetizingPhase}}));
    primary.set_magnetizing_current(sinusoid_signal({{0.8, 0.0}}));
    auto secondary = sinusoid(2.0, loadPhase);
    return make_operating_point({primary, secondary});
}

}  // namespace

TEST_CASE("Test_Phase_Aware_Gap_Fringing_Carries_Magnetizing_Phase", "[physical-model][magnetic-field][phase-aware]") {
    // The gap fringing field is the field of the magnetizing current i_m = sum c_k N_k i_k / N_r,
    // so it carries i_m's phase. Independent check: the model's phasor field must equal
    //     H_turns (fringing off)  +  F * exp(j theta_m)
    // where F is the gap field alone taken from an excitation whose magnetizing current is in
    // phase with the gauge (so it is purely real), and theta_m is computed here from the
    // analytic phasors, not by MKF.
    auto magnetic = make_magnetic({8, 8}, {"primary", "secondary"}, OpenMagneticsTesting::get_ground_gap(0.001));
    const double loadPhase = 0.3;
    const double magnetizingPeak = 0.8;
    const double magnetizingPhase = 0.3 - std::numbers::pi / 2;
    auto operatingPoint = gapped_operating_point(loadPhase, magnetizingPeak, magnetizingPhase);
    // Gauge: the primary current's phase. Phasors of A sin(wt + p) are A exp(j (p - pi/2)).
    std::complex<double> primaryPhasor = std::polar(2.0, loadPhase) + std::polar(magnetizingPeak, magnetizingPhase);
    std::complex<double> secondaryPhasor = std::polar(2.0, loadPhase);
    std::complex<double> magnetizing = (8.0 * primaryPhasor - 8.0 * secondaryPhasor) / 8.0;
    double theta = std::arg(magnetizing) - std::arg(primaryPhasor);

    // F alone: primary carrying only a current in phase with its own (gauge) phase; no turn
    // field (directions 0), so the whole field is the gap field, purely in phase.
    auto reference = make_operating_point({sinusoid(1.0, 0.0), sinusoid(1.0, 0.0)});
    {
        auto primary = reference.get_excitations_per_winding()[0];
        primary.set_magnetizing_current(sinusoid_signal({{0.8, 0.0}}));
        reference.get_mutable_excitations_per_winding()[0] = primary;
        auto secondary = reference.get_excitations_per_winding()[1];
        secondary.set_current(sinusoid_signal({{0.0, 0.0}}));
        reference.get_mutable_excitations_per_winding()[1] = secondary;
    }

    for (auto model : {MagneticFieldStrengthModels::ALBACH, MagneticFieldStrengthModels::BINNS_LAWRENSON}) {
        INFO("model " << static_cast<int>(model));
        settings.reset();
        settings.set_magnetic_field_strength_fringing_effect_model(MagneticFieldStrengthFringingEffectModels::ROSHEN);
        MagneticField magneticField(model, MagneticFieldStrengthFringingEffectModels::ROSHEN);

        settings.set_magnetic_field_include_fringing(true);
        auto full = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic);
        auto gapOnly = magneticField.calculate_magnetic_field_strength_field(reference, magnetic, std::nullopt, std::vector<int8_t>{0, 0});
        settings.set_magnetic_field_include_fringing(false);
        auto turnsOnly = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic);

        REQUIRE(full.get_field_per_frequency().size() == 1);
        const auto& inPhase = full.get_field_per_frequency()[0].get_data();
        const auto& quadrature = full.get_quadrature_field_per_frequency()[0].get_data();
        const auto& turnsInPhase = turnsOnly.get_field_per_frequency()[0].get_data();
        const auto& turnsQuadrature = turnsOnly.get_quadrature_field_per_frequency()[0].get_data();
        const auto& gap = gapOnly.get_field_per_frequency()[0].get_data();
        const auto& gapQuadrature = gapOnly.get_quadrature_field_per_frequency()[0].get_data();
        REQUIRE(inPhase.size() == gap.size());
        REQUIRE(inPhase.size() == turnsInPhase.size());
        double scale = 0;
        double gapScale = 0;
        for (size_t p = 0; p < inPhase.size(); ++p) {
            scale = std::max(scale, std::hypot(inPhase[p].get_real(), inPhase[p].get_imaginary()));
            gapScale = std::max(gapScale, std::hypot(gap[p].get_real(), gap[p].get_imaginary()));
            REQUIRE(gapQuadrature[p].get_real() == 0);
            REQUIRE(gapQuadrature[p].get_imaginary() == 0);
        }
        REQUIRE(gapScale > 1e-3 * scale);  // the gap field is not negligible here
        for (size_t p = 0; p < inPhase.size(); ++p) {
            CHECK_THAT(inPhase[p].get_real(), Catch::Matchers::WithinAbs(turnsInPhase[p].get_real() + gap[p].get_real() * std::cos(theta), 1e-9 * scale));
            CHECK_THAT(inPhase[p].get_imaginary(), Catch::Matchers::WithinAbs(turnsInPhase[p].get_imaginary() + gap[p].get_imaginary() * std::cos(theta), 1e-9 * scale));
            CHECK_THAT(quadrature[p].get_real(), Catch::Matchers::WithinAbs(turnsQuadrature[p].get_real() + gap[p].get_real() * std::sin(theta), 1e-9 * scale));
            CHECK_THAT(quadrature[p].get_imaginary(), Catch::Matchers::WithinAbs(turnsQuadrature[p].get_imaginary() + gap[p].get_imaginary() * std::sin(theta), 1e-9 * scale));
        }

        // The loss move is exactly the magnetizing-phase effect: the same field with the gap
        // term on the gauge phase (theta = 0, what the amplitude-only model did) gives the old
        // loss; with theta it gives the model's loss.
        std::vector<ComplexField> inPhaseAtZero = turnsOnly.get_field_per_frequency();
        for (size_t p = 0; p < inPhase.size(); ++p) {
            auto& point = inPhaseAtZero[0].get_mutable_data()[p];
            point.set_real(point.get_real() + gap[p].get_real());
            point.set_imaginary(point.get_imaginary() + gap[p].get_imaginary());
        }
        settings.set_magnetic_field_include_fringing(true);
        double modelLoss = proximity_of(magnetic, operatingPoint, full.get_field_per_frequency(), full.get_quadrature_field_per_frequency());
        double gaugePhaseLoss = proximity_of(magnetic, operatingPoint, inPhaseAtZero, turnsOnly.get_quadrature_field_per_frequency());
        std::printf("[phase-aware] fringing model %d theta_m %.4f rad: loss %.10g, with the gap field on the gauge phase %.10g\n",
                    static_cast<int>(model), theta, modelLoss, gaugePhaseLoss);
        CHECK(std::abs(modelLoss - gaugePhaseLoss) > 1e-3 * gaugePhaseLoss);
    }
    settings.reset();
}

TEST_CASE("Test_Phase_Aware_Gap_Fringing_Undefined_Magnetizing_Phase_Throws", "[physical-model][magnetic-field][phase-aware]") {
    // Winding currents that cancel exactly (N1 i1 = N2 i2 in the convention) leave the
    // magnetizing current without a phase while the gap still carries a field: throw.
    auto magnetic = make_magnetic({8, 8}, {"primary", "secondary"}, OpenMagneticsTesting::get_ground_gap(0.001));
    auto operatingPoint = gapped_operating_point(0.3, 0.0, 0.0);
    settings.reset();
    settings.set_magnetic_field_include_fringing(true);
    MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON, MagneticFieldStrengthFringingEffectModels::ROSHEN);
    CHECK_THROWS(magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic));
    settings.reset();
}
