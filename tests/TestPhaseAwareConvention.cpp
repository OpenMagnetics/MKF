// MAS excitation convention (2026-09-24): the direction and phase primitives behind the
// phase-aware proximity field (see TestPhaseAwareProximity.cpp for the loss-level checks).
#include "physical_models/MagneticField.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingProximityEffectLosses.h"
#include "processors/Inputs.h"
#include "processors/WaveformProcessor.h"
#include "support/CoilMesher.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include "json.hpp"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <source_location>
#include <limits>
#include <map>
#include <optional>
#include <algorithm>
#include <cctype>
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

// Peak primary voltage that carries a magnetizing current of peak magnetizingPeak on this
// magnetic: v = L_m di_m/dt, so |V| = w L_m |I_m| (L_m from MKF's magnetizing-inductance model).
double magnetizing_voltage_peak(OpenMagnetics::Magnetic& magnetic, double magnetizingPeak) {
    auto core = magnetic.get_core();
    if (!core.get_processed_description()) {
        core.process_data();
        core.process_gap();
    }
    double magnetizingInductance = MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(core, magnetic.get_coil())
                                       .get_magnetizing_inductance().get_nominal().value();
    REQUIRE(magnetizingInductance > 0);
    return 2 * std::numbers::pi * frequency * magnetizingInductance * magnetizingPeak;
}

// A gapped 8:8 transformer: primary = load + an extra term, secondary = load (source
// convention). With several windings the gap fringing field is the field of the magnetizing
// flux, which MKF takes from winding 0's voltage by Faraday (B = V / (j w N A_e)): the primary
// voltage is the one that carries the magnetizing current i_m = magnetizingPeak sin(w t +
// magnetizingPhase) through L_m, v = L_m di_m/dt = w L_m magnetizingPeak sin(w t +
// magnetizingPhase + pi/2), so the flux (and the fringing field) has i_m's phase. i_m itself is
// also recorded as the excitation's magnetizing current, consistent with that voltage.
OperatingPoint gapped_operating_point(double loadPhase, double extraPeak, double extraPhase, double magnetizingPeak, double magnetizingPhase, double voltagePeak) {
    auto primary = sinusoid(0, 0);
    primary.set_current(sinusoid_signal({{2.0, loadPhase}, {extraPeak, extraPhase}}));
    primary.set_magnetizing_current(sinusoid_signal({{magnetizingPeak, magnetizingPhase}}));
    primary.set_voltage(sinusoid_signal({{voltagePeak, magnetizingPhase + std::numbers::pi / 2}}));
    auto secondary = sinusoid(2.0, loadPhase);
    return make_operating_point({primary, secondary});
}

// The model's phasor field must equal H_turns (fringing off) + F exp(j theta_m), with F the gap
// field alone (taken where the fringing phase is 0, so it is purely real) and theta_m the phase
// of the magnetizing flux (Faraday on the primary voltage) relative to the primary current,
// computed here analytically.
// Returns {model loss, loss with the gap field on the gauge phase}.
std::pair<double, double> check_fringing_phase(OpenMagnetics::Magnetic& magnetic, const OperatingPoint& operatingPoint, double theta, double magnetizingPeak, MagneticFieldStrengthModels model) {
    // F alone: the magnetizing current in phase with the primary current (theta = 0) and no turn
    // field (directions 0). Same magnetizing amplitude (same primary voltage amplitude), so the
    // same gap field magnitude.
    auto reference = make_operating_point({sinusoid(1.0, 0.0), sinusoid(0.0, 0.0)});
    {
        auto primary = reference.get_excitations_per_winding()[0];
        primary.set_magnetizing_current(sinusoid_signal({{magnetizingPeak, 0.0}}));
        primary.set_voltage(sinusoid_signal({{magnetizing_voltage_peak(magnetic, magnetizingPeak), std::numbers::pi / 2}}));
        reference.get_mutable_excitations_per_winding()[0] = primary;
    }
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
        CHECK(std::abs(gapQuadrature[p].get_real()) <= 1e-12 * std::max(1.0, std::abs(gap[p].get_real())));
        CHECK(std::abs(gapQuadrature[p].get_imaginary()) <= 1e-12 * std::max(1.0, std::abs(gap[p].get_imaginary())));
    }
    REQUIRE(gapScale > 1e-3 * scale);  // the gap field is not negligible here
    for (size_t p = 0; p < inPhase.size(); ++p) {
        CHECK_THAT(inPhase[p].get_real(), Catch::Matchers::WithinAbs(turnsInPhase[p].get_real() + gap[p].get_real() * std::cos(theta), 1e-9 * scale));
        CHECK_THAT(inPhase[p].get_imaginary(), Catch::Matchers::WithinAbs(turnsInPhase[p].get_imaginary() + gap[p].get_imaginary() * std::cos(theta), 1e-9 * scale));
        CHECK_THAT(quadrature[p].get_real(), Catch::Matchers::WithinAbs(turnsQuadrature[p].get_real() + gap[p].get_real() * std::sin(theta), 1e-9 * scale));
        CHECK_THAT(quadrature[p].get_imaginary(), Catch::Matchers::WithinAbs(turnsQuadrature[p].get_imaginary() + gap[p].get_imaginary() * std::sin(theta), 1e-9 * scale));
    }

    // The loss move is exactly the magnetizing-phase effect: the same field with the gap term on
    // the gauge phase (what the amplitude-only model did) against the model's.
    std::vector<ComplexField> inPhaseAtZero = turnsOnly.get_field_per_frequency();
    for (size_t p = 0; p < inPhase.size(); ++p) {
        auto& point = inPhaseAtZero[0].get_mutable_data()[p];
        point.set_real(point.get_real() + gap[p].get_real());
        point.set_imaginary(point.get_imaginary() + gap[p].get_imaginary());
    }
    settings.set_magnetic_field_include_fringing(true);
    double modelLoss = proximity_of(magnetic, operatingPoint, full.get_field_per_frequency(), full.get_quadrature_field_per_frequency());
    double gaugePhaseLoss = proximity_of(magnetic, operatingPoint, inPhaseAtZero, turnsOnly.get_quadrature_field_per_frequency());
    settings.reset();
    return {modelLoss, gaugePhaseLoss};
}

}  // namespace

TEST_CASE("Test_Phase_Aware_Gap_Fringing_Carries_Magnetizing_Phase", "[physical-model][magnetic-field][phase-aware]") {
    // The gap fringing field is the field of the magnetizing flux, so it carries the magnetizing
    // current's phase, referred to the primary current (the gauge). With two windings MKF takes
    // that flux from the primary voltage (Faraday), which each case sets to the voltage carrying
    // its magnetizing current (gapped_operating_point). Phasors of A sin(wt + p) are
    // A exp(j (p - pi/2)); the -pi/2 cancels in every phase difference.
    auto magnetic = make_magnetic({8, 8}, {"primary", "secondary"}, OpenMagneticsTesting::get_ground_gap(0.001));
    const double loadPhase = 0.3;
    const double magnetizingPeak = 0.8;
    std::complex<double> secondaryPhasor = std::polar(2.0, loadPhase);
    const double voltagePeak = magnetizing_voltage_peak(magnetic, magnetizingPeak);

    for (auto model : {MagneticFieldStrengthModels::ALBACH, MagneticFieldStrengthModels::BINNS_LAWRENSON}) {
        INFO("model " << static_cast<int>(model));

        SECTION("magnetizing current unrelated to the winding currents") {
            // The primary's extra term (0.5 A at -1.0 rad) is NOT the magnetizing current
            // (0.8 A at 0.9 rad, carried by the primary voltage): the fringing phase must follow
            // the magnetizing flux, not the winding currents.
            auto operatingPoint = gapped_operating_point(loadPhase, 0.5, -1.0, magnetizingPeak, 0.9, voltagePeak);
            std::complex<double> primaryPhasor = std::polar(2.0, loadPhase) + std::polar(0.5, -1.0);
            double theta = 0.9 - std::arg(primaryPhasor);
            auto [modelLoss, gaugePhaseLoss] = check_fringing_phase(magnetic, operatingPoint, theta, magnetizingPeak, model);
            std::printf("[phase-aware] fringing model %d unrelated: theta_m %.4f rad, loss %.10g, gap field on the gauge phase %.10g\n",
                        static_cast<int>(model), theta, modelLoss, gaugePhaseLoss);
            CHECK(std::abs(modelLoss - gaugePhaseLoss) > 1e-3 * gaugePhaseLoss);
        }

        SECTION("winding currents satisfy Faraday: both phases agree") {
            // The primary carries load + i_m, with i_m equal to the magnetizing current: the phase
            // of sum c_k N_k i_k / N_r and that of the magnetizing current coincide.
            const double magnetizingPhase = loadPhase - std::numbers::pi / 2;
            auto operatingPoint = gapped_operating_point(loadPhase, magnetizingPeak, magnetizingPhase, magnetizingPeak, magnetizingPhase, voltagePeak);
            std::complex<double> primaryPhasor = std::polar(2.0, loadPhase) + std::polar(magnetizingPeak, magnetizingPhase);
            double thetaFromMagnetizingCurrent = magnetizingPhase - std::arg(primaryPhasor);
            double thetaFromWindingCurrents = std::arg((8.0 * primaryPhasor - 8.0 * secondaryPhasor) / 8.0) - std::arg(primaryPhasor);
            CHECK_THAT(std::remainder(thetaFromMagnetizingCurrent - thetaFromWindingCurrents, 2 * std::numbers::pi), Catch::Matchers::WithinAbs(0, 1e-12));
            auto [modelLoss, gaugePhaseLoss] = check_fringing_phase(magnetic, operatingPoint, thetaFromWindingCurrents, magnetizingPeak, model);
            std::printf("[phase-aware] fringing model %d Faraday: theta_m %.4f rad, loss %.10g, gap field on the gauge phase %.10g\n",
                        static_cast<int>(model), thetaFromWindingCurrents, modelLoss, gaugePhaseLoss);
        }

        SECTION("winding currents cancel (no magnetizing part): computes with the magnetizing current's phase") {
            // N1 i1 = N2 i2 exactly, as reflected-load operating points give. The fringing phase
            // does not come from the winding currents (their net MMF is zero) but from the primary
            // voltage, so this computes.
            auto operatingPoint = gapped_operating_point(loadPhase, 0.0, 0.0, magnetizingPeak, 0.9, voltagePeak);
            double theta = 0.9 - loadPhase;
            auto [modelLoss, gaugePhaseLoss] = check_fringing_phase(magnetic, operatingPoint, theta, magnetizingPeak, model);
            std::printf("[phase-aware] fringing model %d cancelling currents: theta_m %.4f rad, loss %.10g, gap field on the gauge phase %.10g\n",
                        static_cast<int>(model), theta, modelLoss, gaugePhaseLoss);
            CHECK(modelLoss > 0);
        }
    }
}

TEST_CASE("Test_Phase_Aware_Gap_Fringing_Needs_Magnetizing_Waveform", "[physical-model][magnetic-field][phase-aware]") {
    // With several windings the gap field comes only from the magnetizing flux, i.e. from the
    // primary voltage waveform (Faraday): without it, throw (no fallback to the currents or to the
    // derived magnetizing current).
    auto magnetic = make_magnetic({8, 8}, {"primary", "secondary"}, OpenMagneticsTesting::get_ground_gap(0.001));
    auto operatingPoint = gapped_operating_point(0.3, 0.0, 0.0, 0.8, 0.9, magnetizing_voltage_peak(magnetic, 0.8));
    auto primary = operatingPoint.get_excitations_per_winding()[0];
    auto voltage = primary.get_voltage().value();
    voltage.set_waveform(std::nullopt);
    primary.set_voltage(voltage);
    operatingPoint.get_mutable_excitations_per_winding()[0] = primary;
    settings.reset();
    settings.set_magnetic_field_include_fringing(true);
    MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON, MagneticFieldStrengthFringingEffectModels::ROSHEN);
    CHECK_THROWS(magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic));
    settings.reset();
}

// ABT #1880: a harmonic's phase comes from the same Fourier coefficient as its MAS amplitude
// (X_k = 2 c_k, WaveformProcessor::calculate_fourier_series), never from a second, resampled DFT.
namespace {

const double period1880 = 1.0 / frequency;

// Two-tone sine A1 sin(w t + p1) + A2 sin(2 w t + p2): closed axis (128 points over [0, T], the
// last equal to the first) or open axis (128 samples over [0, T)).
Waveform two_tone_1880(bool closed, double a1, double p1, double a2, double p2) {
    const size_t numberPoints = 128;
    const double step = period1880 / static_cast<double>(closed ? numberPoints - 1 : numberPoints);
    std::vector<double> time, data;
    for (size_t i = 0; i < numberPoints; ++i) {
        double t = static_cast<double>(i) * step;
        time.push_back(t);
        double wt = 2 * std::numbers::pi * frequency * t;
        data.push_back(a1 * std::sin(wt + p1) + a2 * std::sin(2 * wt + p2));
    }
    if (closed) {
        time.back() = period1880;
        data.back() = data.front();
    }
    Waveform waveform;
    waveform.set_time(time);
    waveform.set_data(data);
    return waveform;
}

// Triangle between 0 and 1, minimum at tau, maximum at tau + T/2:
//   x(t) = 1/2 - 4/pi^2 sum_{k odd} cos(k w (t - tau)) / k^2, so X_k = -4/(pi^2 k^2) exp(-j k w tau).
double triangle_value_1880(double t, double tau) {
    double u = std::fmod(t - tau + 2 * period1880, period1880);
    return u < period1880 / 2 ? 2 * u / period1880 : 2 * (period1880 - u) / period1880;
}

std::complex<double> triangle_phasor_1880(size_t k, double tau) {
    if (k == 0) {
        return 0.5;
    }
    if (k % 2 == 0) {
        return 0;
    }
    double kd = static_cast<double>(k);
    return -4 / (std::numbers::pi * std::numbers::pi * kd * kd) * std::polar(1.0, -kd * 2 * std::numbers::pi * frequency * tau);
}

// Closed: the four knots of the triangle. Open: 128 samples, tau on the sample grid, so the
// samples' linear interpolant IS the triangle.
Waveform triangle_1880(bool closed, double tau) {
    std::vector<double> time, data;
    if (closed) {
        time = {0, tau, tau + period1880 / 2, period1880};
    }
    else {
        for (size_t i = 0; i < 128; ++i) {
            time.push_back(static_cast<double>(i) * period1880 / 128);
        }
    }
    for (auto t : time) {
        data.push_back(triangle_value_1880(t, tau));
    }
    Waveform waveform;
    waveform.set_time(time);
    waveform.set_data(data);
    return waveform;
}

SignalDescriptor signal_1880(const Waveform& waveform) {
    SignalDescriptor signal;
    signal.set_waveform(waveform);
    auto sampled = OpenMagnetics::Inputs::calculate_sampled_waveform(waveform, frequency);
    signal.set_harmonics(OpenMagnetics::Inputs::calculate_harmonics_data(waveform, sampled, frequency));
    return signal;
}

OperatingPointExcitation excitation_1880(const Waveform& waveform) {
    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    excitation.set_current(signal_1880(waveform));
    return excitation;
}

// A painted SVG with its random CSS class names replaced by the style each one names, and the style block dropped:
// two paintings of the same field compare equal.
std::string normalized_svg_1880(const std::string& svg) {
    // Plain string scans: std::regex recurses per character and overflows the stack on a whole painting.
    const size_t styleBegin = svg.find("<style");
    const size_t styleEnd = svg.find("</style>");
    REQUIRE(styleBegin != std::string::npos);
    REQUIRE(styleEnd != std::string::npos);
    const std::string style = svg.substr(styleBegin, styleEnd - styleBegin);
    std::map<std::string, std::string> stylePerClass;
    for (size_t dot = style.find('.'); dot != std::string::npos; dot = style.find('.', dot + 1)) {
        size_t open = style.find('{', dot);
        size_t close = style.find('}', dot);
        if (open == std::string::npos || close == std::string::npos || open > close) {
            break;
        }
        std::string name = style.substr(dot + 1, open - dot - 1);
        name.erase(std::remove_if(name.begin(), name.end(), [](unsigned char c) { return std::isspace(c); }), name.end());
        std::string body = style.substr(open + 1, close - open - 1);
        body.erase(std::remove_if(body.begin(), body.end(), [](unsigned char c) { return std::isspace(c); }), body.end());
        stylePerClass[name] = body;
        dot = close;
    }
    const std::string body = svg.substr(0, styleBegin) + svg.substr(styleEnd + 8);
    std::string result;
    const std::string marker = "class=\"";
    size_t last = 0;
    for (size_t at = body.find(marker); at != std::string::npos; at = body.find(marker, last)) {
        size_t nameBegin = at + marker.size();
        size_t nameEnd = body.find('"', nameBegin);
        REQUIRE(nameEnd != std::string::npos);
        std::string name = body.substr(nameBegin, nameEnd - nameBegin);
        auto found = stylePerClass.find(name);
        result += body.substr(last, nameBegin - last) + (found == stylePerClass.end() ? name : found->second);
        last = nameEnd;
    }
    result += body.substr(last);
    return result;
}

}  // namespace

TEST_CASE("Test_Harmonic_Phase_Analytic_Waveforms", "[physical-model][phase-aware][abt-1880]") {
    settings.reset();
    // Phase-shifted two-tone sine: the interpolant of its points has, at k = 1 and 2, exactly the
    // sine's phase (p - pi/2 in the exp(+j w t) convention), closed knots and open samples alike.
    for (bool closed : {true, false}) {
        const double p1 = 0.7, p2 = -1.9;
        auto waveform = two_tone_1880(closed, 1.0, p1, 0.3, p2);
        for (auto [k, phase] : std::vector<std::pair<size_t, double>>{{1, p1 - std::numbers::pi / 2}, {2, p2 - std::numbers::pi / 2}}) {
            auto phasor = CoilMesher::calculate_harmonic_phasor(waveform, frequency, frequency * static_cast<double>(k));
            INFO("two-tone sine, " << (closed ? "closed" : "open") << " axis, harmonic " << k);
            CHECK_THAT(std::remainder(std::arg(phasor) - phase, 2 * std::numbers::pi), Catch::Matchers::WithinAbs(0, 1e-9));
        }
    }
    // Shifted triangle: the whole complex phasor of every harmonic, X_k = -4/(pi^2 k^2) exp(-j k w tau).
    for (bool closed : {true, false}) {
        const double tau = period1880 / 8;
        auto waveform = triangle_1880(closed, tau);
        for (size_t k = 0; k <= 9; ++k) {
            auto phasor = CoilMesher::calculate_harmonic_phasor(waveform, frequency, frequency * static_cast<double>(k));
            auto expected = triangle_phasor_1880(k, tau);
            INFO("triangle, " << (closed ? "closed" : "open") << " axis, harmonic " << k << ": " << phasor << " vs " << expected);
            CHECK(std::abs(phasor - expected) <= 1e-9 * std::abs(triangle_phasor_1880(1, tau)));
            if (k % 2 == 1) {
                CHECK_THAT(std::remainder(std::arg(phasor) - std::arg(expected), 2 * std::numbers::pi), Catch::Matchers::WithinAbs(0, 1e-9));
            }
        }
    }
    // A waveform whose axis does not start at 0: the phase is still referred to t = 0.
    {
        const double tau = period1880 / 8;
        // The same periodic triangle given over [-T/4, 3T/4].
        std::vector<double> shiftedTime = {-period1880 / 4, 0, tau, tau + period1880 / 2, 3 * period1880 / 4};
        std::vector<double> shiftedData;
        for (auto t : shiftedTime) {
            shiftedData.push_back(triangle_value_1880(t, tau));
        }
        Waveform shifted;
        shifted.set_time(shiftedTime);
        shifted.set_data(shiftedData);
        for (size_t k : {1, 3, 5}) {
            auto phasor = WaveformProcessor::calculate_harmonic_phasor(shifted, frequency, frequency * static_cast<double>(k));
            INFO("triangle over [-T/4, 3T/4], harmonic " << k);
            CHECK(std::abs(phasor - triangle_phasor_1880(k, tau)) <= 1e-9 * std::abs(triangle_phasor_1880(1, tau)));
        }
    }
}

TEST_CASE("Test_Harmonic_Phase_Is_Arg_Of_Harmonics_Coefficient", "[physical-model][phase-aware][abt-1880]") {
    settings.reset();
    // Two windings with content up to the Nyquist harmonic: a shifted triangle (odd harmonics
    // to k = 63) and a two-tone sine whose higher harmonics are rounding noise (c_k ~ 1e-17).
    // CoilMesher's phase of every listed harmonic is arg(X_k) of the coefficient the harmonics
    // function's amplitude comes from, bit for bit. The noise-level ones are zero both in the list
    // and in the waveform (below 1e-6 of the largest harmonic, WaveformProcessor::is_zero_harmonic):
    // they carry no current and get no phase (ABT #1880 reopen).
    auto magnetic = make_magnetic({4, 4}, {"primary", "secondary"});
    for (bool closed : {true, false}) {
        auto triangleWaveform = triangle_1880(closed, period1880 / 8);
        auto sineWaveform = two_tone_1880(closed, 2.0, 0.4, 0.5, 1.3);
        auto operatingPoint = make_operating_point({excitation_1880(triangleWaveform), excitation_1880(sineWaveform)});
        std::vector<Waveform> waveforms = {triangleWaveform, sineWaveform};
        std::vector<WaveformProcessor::FourierSeries> seriesPerWinding;
        for (auto& waveform : waveforms) {
            seriesPerWinding.push_back(WaveformProcessor::calculate_fourier_series(waveform, frequency, settings.get_inputs_number_points_sampled_waveforms()));
        }
        auto harmonics0 = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
        std::vector<size_t> indexes;
        for (size_t index = 0; index < harmonics0.get_frequencies().size(); ++index) {
            indexes.push_back(index);
        }
        REQUIRE(indexes.size() == 65);
        std::vector<std::optional<double>> gauge;
        std::vector<std::vector<double>> currentAmplitudes;
        auto phases = CoilMesher::calculate_current_phase_per_winding(magnetic.get_coil(), operatingPoint, indexes, &gauge, &currentAmplitudes);
        REQUIRE(phases.size() == indexes.size());
        REQUIRE(currentAmplitudes.size() == indexes.size());
        size_t zeroHarmonics = 0;
        size_t compared = 0;
        for (size_t position = 0; position < indexes.size(); ++position) {
            size_t index = indexes[position];
            std::vector<std::complex<double>> phasorPerWinding;
            std::vector<bool> excited;
            for (size_t windingIndex = 0; windingIndex < 2; ++windingIndex) {
                auto harmonics = operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_harmonics().value();
                double harmonicFrequency = harmonics.get_frequencies()[index];
                auto phasor = WaveformProcessor::calculate_harmonic_phasor(seriesPerWinding[windingIndex], frequency, harmonicFrequency);
                phasorPerWinding.push_back(phasor);
                const bool zero = harmonics.get_amplitudes()[index] == 0 ||
                                  WaveformProcessor::is_zero_harmonic(seriesPerWinding[windingIndex], harmonics.get_amplitudes()[index], phasor);
                excited.push_back(!zero);
                if (zero && harmonics.get_amplitudes()[index] != 0) {
                    ++zeroHarmonics;
                }
                CHECK(currentAmplitudes[position][windingIndex] == (zero ? 0.0 : harmonics.get_amplitudes()[index]));
                // The amplitude MAS lists is |X_k| of the same coefficient.
                if (excited.back()) {
                    INFO((closed ? "closed" : "open") << " axis, winding " << windingIndex << ", harmonic " << index);
                    CHECK(std::abs(std::abs(phasor) - harmonics.get_amplitudes()[index]) <= 4 * std::numeric_limits<double>::epsilon() * harmonics.get_amplitudes()[index]);
                }
            }
            if (!excited[0] && !excited[1]) {
                continue;
            }
            size_t gaugeWinding = excited[0] ? 0 : 1;
            double gaugePhase = std::arg(phasorPerWinding[gaugeWinding]);
            INFO((closed ? "closed" : "open") << " axis, harmonic " << index);
            REQUIRE(gauge[position].has_value());
            CHECK(gauge[position].value() == gaugePhase);
            for (size_t windingIndex = 0; windingIndex < 2; ++windingIndex) {
                double expected = excited[windingIndex] ? std::arg(phasorPerWinding[windingIndex]) - gaugePhase : 0.0;
                CHECK(phases[position][windingIndex] == expected);
                ++compared;
            }
        }
        // Triangle: DC and the 32 odd harmonics; sine: k = 1, 2 (k = 1 shared). Every other listed one is noise.
        CHECK(compared == 2 * 34);
        CHECK(zeroHarmonics > 0);
    }
}

TEST_CASE("Test_Harmonic_Phase_Inconsistent_Harmonics_Throw", "[physical-model][phase-aware][abt-1880]") {
    settings.reset();
    auto magnetic = make_magnetic({4, 4}, {"primary", "secondary"});
    auto sineWaveform = two_tone_1880(true, 2.0, 0.4, 0.0, 0.0);
    auto triangleWaveform = triangle_1880(true, period1880 / 8);

    // The harmonics list claims 0.5 A at 3 f; the pure sine it belongs to has none there (c_3 ~ 1e-17).
    {
        auto operatingPoint = make_operating_point({excitation_1880(triangleWaveform), excitation_1880(sineWaveform)});
        auto current = operatingPoint.get_excitations_per_winding()[1].get_current().value();
        auto harmonics = current.get_harmonics().value();
        REQUIRE(harmonics.get_frequencies()[3] == 3 * frequency);
        harmonics.get_mutable_amplitudes()[3] = 0.5;
        current.set_harmonics(harmonics);
        operatingPoint.get_mutable_excitations_per_winding()[1].set_current(current);
        CHECK_THROWS(CoilMesher::calculate_current_phase_per_winding(magnetic.get_coil(), operatingPoint, {3}));
        // The consistent list goes through.
        CHECK_NOTHROW(CoilMesher::calculate_current_phase_per_winding(magnetic.get_coil(), operatingPoint, {1}));
    }
    // A listed frequency that is not a harmonic of the excitation (1.5 f) has no coefficient.
    {
        auto operatingPoint = make_operating_point({excitation_1880(triangleWaveform), excitation_1880(sineWaveform)});
        auto current = operatingPoint.get_excitations_per_winding()[1].get_current().value();
        auto harmonics = current.get_harmonics().value();
        harmonics.get_mutable_frequencies()[1] = 1.5 * frequency;
        current.set_harmonics(harmonics);
        operatingPoint.get_mutable_excitations_per_winding()[1].set_current(current);
        CHECK_THROWS(CoilMesher::calculate_current_phase_per_winding(magnetic.get_coil(), operatingPoint, {1}));
    }
    // Beyond the series (k = 65 on a 128-point waveform) and off-grid frequencies throw.
    CHECK_THROWS(CoilMesher::calculate_harmonic_phasor(sineWaveform, frequency, 65 * frequency));
    CHECK_THROWS(CoilMesher::calculate_harmonic_phasor(sineWaveform, frequency, 2.5 * frequency));
}

TEST_CASE("Test_Harmonic_Phase_Zero_Harmonics_Contribute_Nothing", "[physical-model][phase-aware][abt-1880]") {
    // The web's test_wound_coil.json (ER 25/6/15, 15:4x3 turns): winding 1's current is a 3-sample symmetric triangle
    // at 100 kHz whose even harmonics are zero both in its MAS list (7.1e-15 A at 200 kHz) and in its waveform. A zero
    // harmonic carries no current whatever its phase: the phases, the mesh and the field plot must go through, and it
    // must contribute nothing (ABT #1880 reopen: plot_magnetic_field threw "its phase is undefined").
    settings.reset();
    clear_databases();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_wound_coil_zero_even_harmonics.json");
    // Parsed as the WASM's plot_magnetic_field does (Magnetic and OperatingPoint straight from JSON): through
    // OpenMagnetics::Inputs the operating point would be re-processed and its stored harmonics replaced.
    std::ifstream file(path);
    auto masJson = json::parse(file);
    OpenMagnetics::Magnetic magnetic(masJson["magnetic"]);
    OperatingPoint operatingPoint(masJson["inputs"]["operatingPoints"][0]);
    auto coil = magnetic.get_coil();

    auto current = operatingPoint.get_excitations_per_winding()[1].get_current().value();
    REQUIRE(current.get_waveform()->get_data() == std::vector<double>{-20, 20, -20});
    auto harmonics = current.get_harmonics().value();
    REQUIRE(harmonics.get_frequencies()[2] == 200000);
    // Not an exact zero: the case the reopen is about.
    REQUIRE(harmonics.get_amplitudes()[2] > 0);
    REQUIRE(harmonics.get_amplitudes()[2] < 1e-13);
    auto series = WaveformProcessor::calculate_waveform_fourier_series(current.get_waveform().value(), frequency, settings.get_inputs_number_points_sampled_waveforms());
    const double largest = WaveformProcessor::calculate_largest_harmonic_amplitude(series);
    // The scale is the fundamental. The stored list was computed by an older MKF (2e-4 off the exact series).
    REQUIRE_THAT(largest, Catch::Matchers::WithinRel(harmonics.get_amplitudes()[1], 1e-3));

    // Per-winding phases over every listed harmonic: the zero ones carry 0 A and a placeholder phase 0, never a gauge.
    std::vector<size_t> indexes;
    for (size_t index = 0; index < harmonics.get_frequencies().size(); ++index) {
        indexes.push_back(index);
    }
    std::vector<std::optional<double>> gauge;
    std::vector<std::vector<double>> currentAmplitudes;
    std::vector<std::vector<double>> phases;
    REQUIRE_NOTHROW(phases = CoilMesher::calculate_current_phase_per_winding(coil, operatingPoint, indexes, &gauge, &currentAmplitudes));
    std::vector<size_t> zeroIndexes;
    for (size_t position = 0; position < indexes.size(); ++position) {
        size_t index = indexes[position];
        auto phasor = WaveformProcessor::calculate_harmonic_phasor(series, frequency, harmonics.get_frequencies()[index]);
        INFO("winding 1, harmonic " << index << ": listed " << harmonics.get_amplitudes()[index] << " A, waveform " << std::abs(phasor) << " A");
        if (WaveformProcessor::is_zero_harmonic(series, harmonics.get_amplitudes()[index], phasor)) {
            zeroIndexes.push_back(index);
            CHECK(std::abs(phasor) < 1e-6 * largest);
            CHECK(currentAmplitudes[position][1] == 0);
            CHECK(phases[position][1] == 0);
        }
        else {
            CHECK(currentAmplitudes[position][1] == harmonics.get_amplitudes()[index]);
        }
    }
    // DC and the 31 even harmonics are zero; the 32 odd ones carry the triangle.
    std::vector<size_t> evenIndexes;
    for (size_t index = 0; index < indexes.size(); index += 2) {
        evenIndexes.push_back(index);
    }
    CHECK(zeroIndexes == evenIndexes);

    // Contributes nothing: the same operating point with winding 1's zero harmonics written as exact zeros gives the
    // identical mesh (every point, every harmonic), phases and painted field, bit for bit.
    auto exactZeroOperatingPoint = operatingPoint;
    {
        auto exactZeroCurrent = current;
        auto exactZeroHarmonics = harmonics;
        size_t written = 0;
        for (auto index : zeroIndexes) {
            if (exactZeroHarmonics.get_amplitudes()[index] != 0) {
                exactZeroHarmonics.get_mutable_amplitudes()[index] = 0;
                ++written;
            }
        }
        REQUIRE(written == zeroIndexes.size());
        exactZeroCurrent.set_harmonics(exactZeroHarmonics);
        exactZeroOperatingPoint.get_mutable_excitations_per_winding()[1].set_current(exactZeroCurrent);
    }
    auto mesh = CoilMesher().generate_mesh_inducing_coil_phasors(magnetic, operatingPoint, settings.get_harmonic_amplitude_threshold());
    auto exactZeroMesh = CoilMesher().generate_mesh_inducing_coil_phasors(magnetic, exactZeroOperatingPoint, settings.get_harmonic_amplitude_threshold());
    REQUIRE(mesh.fieldPerHarmonic.size() == exactZeroMesh.fieldPerHarmonic.size());
    CHECK(mesh.currentPhasePerHarmonicPerWinding == exactZeroMesh.currentPhasePerHarmonicPerWinding);
    CHECK(mesh.gaugePhasePerHarmonic == exactZeroMesh.gaugePhasePerHarmonic);
    CHECK(mesh.currentAmplitudePerHarmonicPerWinding == exactZeroMesh.currentAmplitudePerHarmonicPerWinding);
    // The mesher reads harmonic i of every winding, and winding 0 lists only [0, 100k, 300k]: at i = 2 winding 1's
    // 200 kHz zero is meshed beside winding 0's 300 kHz. That position must carry winding 1's zero.
    size_t positionsWithAZero = 0;
    for (size_t position = 0; position < mesh.fieldPerHarmonic.size(); ++position) {
        INFO("meshed harmonic " << position << " at " << mesh.fieldPerHarmonic[position].get_frequency() << " Hz");
        const auto& data = mesh.fieldPerHarmonic[position].get_data();
        const auto& exactZeroData = exactZeroMesh.fieldPerHarmonic[position].get_data();
        REQUIRE(data.size() == exactZeroData.size());
        size_t differing = 0;
        for (size_t pointIndex = 0; pointIndex < data.size(); ++pointIndex) {
            if (data[pointIndex].get_value() != exactZeroData[pointIndex].get_value()) {
                ++differing;
            }
        }
        CHECK(differing == 0);
        if (mesh.currentAmplitudePerHarmonicPerWinding[position][0] > 0 && mesh.currentAmplitudePerHarmonicPerWinding[position][1] == 0) {
            ++positionsWithAZero;
        }
    }
    CHECK(positionsWithAZero > 0);

    // plot_magnetic_field as the WASM calls it, with and without the gap fringing: it draws, and draws exactly what the
    // exact-zero list draws.
    for (bool fringing : {true, false}) {
        INFO("fringing " << fringing);
        settings.set_painter_include_fringing(fringing);
        std::filesystem::path emptyFilepath;
        std::string svg, exactZeroSvg;
        {
            Painter painter(emptyFilepath);
            REQUIRE_NOTHROW(painter.paint_magnetic_field(operatingPoint, magnetic));
            svg = painter.export_svg();
        }
        {
            Painter painter(emptyFilepath);
            REQUIRE_NOTHROW(painter.paint_magnetic_field(exactZeroOperatingPoint, magnetic));
            exactZeroSvg = painter.export_svg();
        }
        CHECK_FALSE(svg.empty());
        CHECK(normalized_svg_1880(svg) == normalized_svg_1880(exactZeroSvg));
    }
    settings.reset();
}
