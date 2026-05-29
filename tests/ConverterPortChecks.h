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
#include "MAS.hpp"
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

// Trapezoidal time-weighted integral of v(t) over the full time range.
// Required for ngspice waveforms whose time grid is non-uniform (denser
// around switching events) — a naive sum/N over-weights the dense regions
// and is wrong by an order of magnitude on power-stage waveforms.
inline double time_weighted_mean(const std::vector<double>& t,
                                 const std::vector<double>& v) {
    if (t.size() < 2 || t.size() != v.size()) return 0.0;
    double area = 0.0;
    for (size_t i = 1; i < t.size(); ++i) {
        area += 0.5 * (v[i] + v[i-1]) * (t[i] - t[i-1]);
    }
    return area / (t.back() - t.front());
}

inline double time_weighted_rms(const std::vector<double>& t,
                                const std::vector<double>& v) {
    if (t.size() < 2 || t.size() != v.size()) return 0.0;
    double area = 0.0;
    for (size_t i = 1; i < t.size(); ++i) {
        const double s2 = 0.5 * (v[i]*v[i] + v[i-1]*v[i-1]);
        area += s2 * (t[i] - t[i-1]);
    }
    return std::sqrt(area / (t.back() - t.front()));
}

inline double time_weighted_mean_product(const std::vector<double>& t,
                                         const std::vector<double>& a,
                                         const std::vector<double>& b) {
    if (t.size() < 2 || t.size() != a.size() || t.size() != b.size()) return 0.0;
    double area = 0.0;
    for (size_t i = 1; i < t.size(); ++i) {
        const double pcur  = a[i]   * b[i];
        const double pprev = a[i-1] * b[i-1];
        area += 0.5 * (pcur + pprev) * (t[i] - t[i-1]);
    }
    return area / (t.back() - t.front());
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

// ─────────────────────────────────────────────────────────────────────────────
// PFC switching-leg port-checks (Phase-6 PtP gate).
//
// Consumes the OperatingPoint produced by
// PowerFactorCorrection::simulate_with_ngspice_switching, which packs four
// diagnostic windings:
//
//   Winding 0  InputPort       V = vin_rect   I = i(L_boost)
//   Winding 1  PowerStage      V = vbus       I = i(L_boost)
//   Winding 2  VoltageLoop     V = vea        I = i_ref
//   Winding 3  CurrentLoop     V = vc_i       I = i_sense
//
// Applies the same DC/AC port semantics as `check_pfc_ports` but on the
// SPICE-simulated waveforms instead of the analytical ConverterWaveforms:
//
//   InputPort.voltage   line-cycle MEAN ≈ 0.9·Vrms                vinTol
//                       line-cycle RMS  ≈     Vrms                vinTol
//   PowerStage.voltage  MEAN            ≈ Vbus_nominal            voutMeanTol
//   Power balance       <vin·iL>        ≈ Pout_nominal            pinTol
//
// (We gate Pin = <vin·iL> rather than <iL> alone because PFC zero-crossing
// distortion can skew the inductor-current line-mean a few percent above
// the analytical 2·Iin_pk/π value while still preserving real power
// balance — <vin·iL> is the physically meaningful quantity.)
//
// The OperatingPoint MUST already be trimmed to a whole number of full
// line cycles (the default `trimToLastLineCycle = true` of
// simulate_with_ngspice_switching satisfies this).
// ─────────────────────────────────────────────────────────────────────────────
// `bipolarInput` selects the input-port semantics:
//   • false (bridged boost family): InputPort.voltage is the rectified line,
//     so its mean ≈ (2√2/π)·Vrms ≈ 0.9·Vrms and RMS ≈ Vrms.
//   • true  (bridgeless totem-pole): InputPort.voltage is the SIGNED AC line,
//     so its mean ≈ 0 and RMS ≈ Vrms. Power balance <vin·iL> is identical in
//     both cases — for the totem-pole both vin and iL are bipolar and in
//     phase, so their product is positive over the whole line cycle.
inline void check_pfc_switching_ports(const MAS::OperatingPoint&            op,
                                      const std::string&                    topoName,
                                      double                                vinRms,
                                      double                                voutNom,
                                      double                                poutNom,
                                      double                                vinTol      = 0.05,
                                      double                                voutMeanTol = kVoutMeanTol,
                                      double                                pinTol      = 0.05,
                                      bool                                  bipolarInput = false) {
    const auto& wnds = op.get_excitations_per_winding();
    REQUIRE(wnds.size() >= 2);

    auto findWinding = [&](const std::string& name) -> const MAS::OperatingPointExcitation* {
        for (const auto& w : wnds) {
            if (w.get_name().has_value() && w.get_name().value() == name) return &w;
        }
        return nullptr;
    };
    const auto* inputPort  = findWinding("InputPort");
    const auto* powerStage = findWinding("PowerStage");
    REQUIRE(inputPort  != nullptr);
    REQUIRE(powerStage != nullptr);

    // ---- Input port: rectified line ----------------------------------------
    // NOTE: MAS::OperatingPointExcitation::get_voltage() returns
    // std::optional<SignalDescriptor> *by value* — every chained
    // .get_xxx() call yields a fresh temporary, so we MUST copy each
    // optional into a named local before reaching the underlying data,
    // otherwise we read from destroyed temporaries.
    auto vinV = inputPort->get_voltage();
    REQUIRE(vinV.has_value());
    auto vinWfOpt = vinV->get_waveform();
    REQUIRE(vinWfOpt.has_value());
    const MAS::Waveform vinWf = vinWfOpt.value();
    auto vinTimeOpt = vinWf.get_time();
    REQUIRE(vinTimeOpt.has_value());
    const std::vector<double>& tvec    = vinTimeOpt.value();
    const std::vector<double>& vinData = vinWf.get_data();
    REQUIRE(!vinData.empty());
    REQUIRE(tvec.size() == vinData.size());
    const double vinMean       = time_weighted_mean(tvec, vinData);
    const double vinRmsAct     = time_weighted_rms (tvec, vinData);
    // Bridged input is rectified (mean ≈ 0.9·Vrms); bridgeless totem-pole
    // input is the signed line (mean ≈ 0). RMS ≈ Vrms in both cases.
    const double expectVinMean =
        bipolarInput ? 0.0 : 2.0 * std::sqrt(2.0) * vinRms / M_PI;
    INFO(topoName << " switching InputPort.voltage: mean=" << vinMean
         << " (expected " << expectVinMean << "), rms=" << vinRmsAct
         << " (expected " << vinRms << ")");
    if (bipolarInput)
        CHECK(std::fabs(vinMean) / vinRms < vinTol);   // mean ≈ 0 (signed line)
    else
        CHECK(std::fabs(vinMean - expectVinMean) / expectVinMean < vinTol);
    CHECK(std::fabs(vinRmsAct - vinRms)        / vinRms        < vinTol);

    // ---- Power-stage port: bus voltage --------------------------------------
    auto vbusV = powerStage->get_voltage();
    REQUIRE(vbusV.has_value());
    auto vbusWfOpt = vbusV->get_waveform();
    REQUIRE(vbusWfOpt.has_value());
    const MAS::Waveform vbusWf = vbusWfOpt.value();
    const std::vector<double>& vbusData = vbusWf.get_data();
    auto vbusTimeOpt = vbusWf.get_time();
    REQUIRE(vbusTimeOpt.has_value());
    const std::vector<double>& tvecVbus = vbusTimeOpt.value();
    REQUIRE(!vbusData.empty());
    const double vbusMean = time_weighted_mean(tvecVbus, vbusData);
    INFO(topoName << " switching PowerStage.voltage: mean=" << vbusMean
         << " (nom " << voutNom << ")");
    CHECK(std::fabs(vbusMean - voutNom) / voutNom < voutMeanTol);

    // ---- Power balance: <vin · iL> -----------------------------------------
    auto iLI = powerStage->get_current();
    REQUIRE(iLI.has_value());
    auto iLWfOpt = iLI->get_waveform();
    REQUIRE(iLWfOpt.has_value());
    const MAS::Waveform iLWf = iLWfOpt.value();
    const std::vector<double>& iLData = iLWf.get_data();
    REQUIRE(!iLData.empty());
    REQUIRE(iLData.size() == vinData.size());
    const double pinAvg = time_weighted_mean_product(tvec, vinData, iLData);
    INFO(topoName << " switching power balance <vin·iL>: " << pinAvg
         << " W (nom " << poutNom << ", tol ±" << 100.0*pinTol << " %)");
    CHECK(std::fabs(pinAvg - poutNom) / poutNom < pinTol);
}

} // namespace ConverterPortChecks
