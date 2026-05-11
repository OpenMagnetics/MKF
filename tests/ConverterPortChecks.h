#pragma once
// ─────────────────────────────────────────────────────────────────────────────
// ConverterPortChecks.h — §5.1 converter-port DC-stream gate helper.
//
// Used by every per-topology Test_<Topo>_ConverterPortWaveforms TEST_CASE
// (see CONVERTER_MODELS_GOLDEN_GUIDE.md §5.1).
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
//
// This header deliberately does NOT pull in `using namespace OpenMagnetics;`
// or `using namespace MAS;` (unlike NgspiceTestHelpers.h, which is
// preserved for backwards compatibility) so it can be safely included
// from test files that need both namespaces visible without ambiguity
// (e.g. TestTopologyAsymmetricHalfBridge.cpp).
// ─────────────────────────────────────────────────────────────────────────────

#include <catch2/catch_test_macros.hpp>
#include "converter_models/Topology.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

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

inline double rms(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x * x;
    return std::sqrt(s / static_cast<double>(v.size()));
}

inline void check_dc_ports(const OpenMagnetics::ConverterWaveforms& w,
                           const std::string&                        topoName,
                           size_t                                    opIdx,
                           double                                    vinNom,
                           const std::vector<double>&                voutNom,
                           const std::vector<double>&                ioutNom,
                           double                                    voutMeanTol = kVoutMeanTol,
                           double                                    ioutMeanTol = kIoutMeanTol) {
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

// ─────────────────────────────────────────────────────────────────────────────
// PFC-specific port-checks (line-frequency input → DC bus output).
// The §5.1 DC-stream gate (`check_dc_ports`) does NOT apply to a PFC's input
// port because that port is fundamentally AC (rectified line).
//
// Bounds:
//   input_voltage   line-cycle MEAN ≈ 2·√2·Vrms / π ≈ 0.9·Vrms     (vinTol)
//                   line-cycle RMS  ≈ Vrms                          (vinTol)
//   output_voltage  MEAN            ≈ Vbus_nominal                  (voutMeanTol)
//   output_current  MEAN            ≈ Iload_nominal                 (ioutMeanTol)
//
// (No bound on input_current here — PF/THD are validated separately by
// PfcSimulationWaveforms.powerFactor / .currentThd in the smoke tests.)
// ─────────────────────────────────────────────────────────────────────────────
inline void check_pfc_ports(const OpenMagnetics::ConverterWaveforms& w,
                            const std::string&                        topoName,
                            size_t                                    opIdx,
                            double                                    vinRms,
                            double                                    voutNom,
                            double                                    ioutNom,
                            double                                    vinTol      = 0.05,
                            double                                    voutMeanTol = kVoutMeanTol,
                            double                                    ioutMeanTol = kIoutMeanTol) {
    const auto& vinData = w.get_input_voltage().get_data();
    REQUIRE(!vinData.empty());
    const double vinMean   = mean(vinData);
    const double vinRmsAct = rms(vinData);
    const double expectedVinMean = 2.0 * std::sqrt(2.0) * vinRms / M_PI; // ≈ 0.9·Vrms
    INFO(topoName << " OP " << opIdx
         << " input_voltage (rectified): mean=" << vinMean
         << " (expected " << expectedVinMean << "), rms=" << vinRmsAct
         << " (expected " << vinRms << ")");
    CHECK(std::fabs(vinMean   - expectedVinMean) / expectedVinMean < vinTol);
    CHECK(std::fabs(vinRmsAct - vinRms)          / vinRms          < vinTol);

    const auto& voutWfs = w.get_output_voltages();
    REQUIRE(voutWfs.size() == 1);
    const auto& voutData = voutWfs[0].get_data();
    REQUIRE(!voutData.empty());
    const double voutMean = mean(voutData);
    INFO(topoName << " OP " << opIdx
         << " output_voltage: mean=" << voutMean << " (nom " << voutNom << ")");
    CHECK(std::fabs(voutMean - voutNom) / voutNom < voutMeanTol);

    const auto& ioutWfs = w.get_output_currents();
    REQUIRE(ioutWfs.size() == 1);
    const auto& ioutData = ioutWfs[0].get_data();
    REQUIRE(!ioutData.empty());
    const double ioutMean = mean(ioutData);
    INFO(topoName << " OP " << opIdx
         << " output_current: mean=" << ioutMean << " (nom " << ioutNom << ")");
    CHECK(std::fabs(ioutMean - ioutNom) / ioutNom < ioutMeanTol);
}

} // namespace ConverterPortChecks
