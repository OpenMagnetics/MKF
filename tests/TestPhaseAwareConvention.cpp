// MAS excitation convention (2026-09-24): the direction and phase primitives behind the
// phase-aware proximity field (see TestPhaseAwareProximity.cpp for the loss-level checks).
#include "physical_models/MagneticField.h"
#include "physical_models/LeakageInductance.h"
#include "processors/Inputs.h"
#include "support/CoilMesher.h"
#include "support/Settings.h"
#include "json.hpp"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

const double frequency = 100000;

OperatingPointExcitation sinusoid(double peak, double phase) {
    const size_t numberPoints = 256;
    std::vector<double> time;
    std::vector<double> data;
    for (size_t i = 0; i <= numberPoints; ++i) {
        double t = static_cast<double>(i) / static_cast<double>(numberPoints) / frequency;
        time.push_back(t);
        data.push_back(peak * std::sin(2 * std::numbers::pi * frequency * t + phase));
    }
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

OpenMagnetics::Magnetic make_magnetic(std::vector<int64_t> numberTurns, std::vector<std::string> isolationSides) {
    auto core = OpenMagneticsTesting::get_quick_core("ETD 34", json::array(), 1, "3C97");
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
