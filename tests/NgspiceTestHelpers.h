#pragma once
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "processors/NgspiceRunner.h"
#include "processors/Inputs.h"
#include "converter_models/Topology.h"
#include <MAS.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

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

// ─────────────────────────────────────────────────────────────────────────────
// §5.1 converter-port DC-stream gate — shared helpers used by every
// per-topology Test_<Topo>_ConverterPortWaveforms TEST_CASE.
//
// The signals returned by simulate_and_extract_topology_waveforms() —
//   ConverterWaveforms.{input_voltage, input_current,
//                       output_voltages, output_currents}
// — represent the converter's *external* ports (DC source rail, DC
// filtered output rails), NOT the AC magnetic-component winding ports.
// They MUST therefore be DC (or as DC as the chosen output filter
// allows): mean ≈ nameplate, ripple bounded by Cout / load.
//
// A signal that violates this is almost certainly a winding-port
// (bipolar AC) signal smuggled into the converter-port stream — the
// inverse of the bug TestVoltSecondBalance.cpp catches on the
// winding-port stream.
//
// Bounds (loose enough to absorb startup transients, switching ripple,
// and ESR; tighten case-by-case as topologies stabilise):
//   input_voltage:   |mean − Vin_nom|/Vin_nom        < 1 %
//                    (max−min)/|mean|                < 1 %  (hard DC src)
//   output_voltage:  |mean − Vout_nom|/Vout_nom      < voutMeanTol
//                    (max−min)/|mean|                < 25 % (filtered DC)
//   output_current:  |mean − Iout_nom|/Iout_nom      < ioutMeanTol
// (No input_current bound — too sensitive to η model and startup inrush.)
// ─────────────────────────────────────────────────────────────────────────────
namespace ConverterPortChecks {

inline constexpr double kVinMeanTol    = 0.01;
inline constexpr double kVinRippleTol  = 0.01;
inline constexpr double kVoutMeanTol   = 0.10;
inline constexpr double kVoutRippleTol = 0.25;
inline constexpr double kIoutMeanTol   = 0.25;

inline double mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
}

inline double ripple_over_mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    auto [mn, mx] = std::minmax_element(v.begin(), v.end());
    double m = mean(v);
    if (std::fabs(m) < 1e-12) return 0.0;
    return (*mx - *mn) / std::fabs(m);
}

inline void check_dc_ports(const OpenMagnetics::ConverterWaveforms& w,
                           const std::string&             topoName,
                           size_t                         opIdx,
                           double                         vinNom,
                           const std::vector<double>&     voutNom,
                           const std::vector<double>&     ioutNom,
                           double                         voutMeanTol = kVoutMeanTol,
                           double                         ioutMeanTol = kIoutMeanTol) {
    const auto& vinData = w.get_input_voltage().get_data();
    REQUIRE(!vinData.empty());
    const double vinMean   = mean(vinData);
    const double vinRipple = ripple_over_mean(vinData);
    INFO(topoName << " OP " << opIdx
         << " input_voltage: mean=" << vinMean
         << " (nom " << vinNom << "), ripple/mean=" << vinRipple);
    CHECK(std::fabs(vinMean - vinNom) / vinNom < kVinMeanTol);
    CHECK(vinRipple < kVinRippleTol);

    const auto& voutWfs = w.get_output_voltages();
    REQUIRE(voutWfs.size() == voutNom.size());
    for (size_t i = 0; i < voutWfs.size(); ++i) {
        const auto& d = voutWfs[i].get_data();
        REQUIRE(!d.empty());
        const double m = mean(d);
        const double r = ripple_over_mean(d);
        INFO(topoName << " OP " << opIdx
             << " output_voltage[" << i << "]: mean=" << m
             << " (nom " << voutNom[i] << "), ripple/mean=" << r);
        CHECK(std::fabs(m - voutNom[i]) / voutNom[i] < voutMeanTol);
        CHECK(r < kVoutRippleTol);
    }

    const auto& ioutWfs = w.get_output_currents();
    REQUIRE(ioutWfs.size() == ioutNom.size());
    for (size_t i = 0; i < ioutWfs.size(); ++i) {
        const auto& d = ioutWfs[i].get_data();
        REQUIRE(!d.empty());
        const double m = mean(d);
        INFO(topoName << " OP " << opIdx
             << " output_current[" << i << "]: mean=" << m
             << " (nom " << ioutNom[i] << ")");
        CHECK(std::fabs(m - ioutNom[i]) / ioutNom[i] < ioutMeanTol);
    }
}

} // namespace ConverterPortChecks
