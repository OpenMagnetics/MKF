#pragma once
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "processors/NgspiceRunner.h"
#include "processors/Inputs.h"
#include <MAS.hpp>

using namespace OpenMagnetics;
using namespace MAS;

namespace NgspiceTestHelpers {

/**
 * @brief Skip test if ngspice is not available
 */
inline void skip_if_ngspice_unavailable() {
    NgspiceRunner runner;
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
}

/**
 * @brief Check if ngspice is available
 */
inline bool is_ngspice_available() {
    NgspiceRunner runner;
    return runner.is_available();
}

/**
 * @brief Setup common converter input voltage
 */
inline DimensionWithTolerance setup_input_voltage(double nominal, double min, double max) {
    DimensionWithTolerance inputVoltage;
    inputVoltage.set_nominal(nominal);
    inputVoltage.set_minimum(min);
    inputVoltage.set_maximum(max);
    return inputVoltage;
}

/**
 * @brief Validate voltage waveform bounds
 */
inline void validate_voltage_bounds(const Waveform& waveform, double expectedMin, double expectedMax, 
                                    double tolerancePercent = 0.2) {
    auto data = waveform.get_data();
    double actualMax = *std::max_element(data.begin(), data.end());
    double actualMin = *std::min_element(data.begin(), data.end());
    
    CHECK(actualMax >= expectedMax * (1 - tolerancePercent));
    CHECK(actualMax <= expectedMax * (1 + tolerancePercent));
    CHECK(actualMin >= expectedMin * (1 - tolerancePercent));
    CHECK(actualMin <= expectedMin * (1 + tolerancePercent));
}

/**
 * @brief Validate current is within expected range
 */
inline void validate_current_range(const Waveform& waveform, double minExpected, double maxExpected) {
    auto data = waveform.get_data();
    double actualMax = *std::max_element(data.begin(), data.end());
    double actualMin = *std::min_element(data.begin(), data.end());
    
    CHECK(actualMax >= minExpected);
    CHECK(actualMax <= maxExpected);
    CHECK(actualMin >= 0);  // Current should not go negative for most converters
}

/**
 * @brief Validate operating point has required windings
 */
inline void validate_operating_point(const OperatingPoint& op, size_t expectedWindingCount) {
    REQUIRE(!op.get_excitations_per_winding().empty());
    REQUIRE(op.get_excitations_per_winding().size() >= expectedWindingCount);
}

/**
 * @brief Run simulation with fallback to analytical
 */
template<typename Converter>
inline std::vector<OperatingPoint> simulate_with_fallback(
    Converter& converter,
    const std::vector<double>& turnsRatios,
    double magnetizingInductance) {
    
    if (!is_ngspice_available()) {
        // Return analytical calculation
        return converter.process_operating_points(turnsRatios, magnetizingInductance);
    }
    
    auto result = converter.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
    
    if (result.empty()) {
        // Fallback to analytical
        return converter.process_operating_points(turnsRatios, magnetizingInductance);
    }
    
    return result;
}

/**
 * @brief Common test configuration for converters
 */
struct ConverterTestConfig {
    double inputVoltageNominal = 48.0;
    double inputVoltageMin = 36.0;
    double inputVoltageMax = 72.0;
    double outputVoltage = 12.0;
    double outputCurrent = 5.0;
    double switchingFrequency = 100000.0;
    double diodeVoltageDrop = 0.5;
    double efficiency = 0.9;
    double currentRippleRatio = 0.3;
};

/**
 * @brief Setup common converter test configuration
 */
template<typename Converter>
inline void setup_converter_common(Converter& converter, const ConverterTestConfig& config = {}) {
    converter.set_input_voltage(setup_input_voltage(
        config.inputVoltageNominal, config.inputVoltageMin, config.inputVoltageMax));
    converter.set_diode_voltage_drop(config.diodeVoltageDrop);
    converter.set_efficiency(config.efficiency);
    converter.set_current_ripple_ratio(config.currentRippleRatio);
}

} // namespace NgspiceTestHelpers
