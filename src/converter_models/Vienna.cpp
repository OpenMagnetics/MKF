// =====================================================================
// Vienna rectifier — Phase 1+2+3 (incremental).
//
// Phase 1: skeleton (constructor, run_checks, design-requirements stub).
// Phase 2: per-phase analytical solver at peak-of-line (single sample).
// Phase 3 (this file): adds samplingStrategy=fullLineCycle — each of the
//   three "Phase A/B/C" windings carries the full 50/60 Hz line cycle,
//   shifted by 120° per phase, so the wizard plot shows a recognisable
//   3-phase pattern instead of three identical peak-of-line snapshots.
//   The remaining Phase-3 items (viennaII, alternative switch types,
//   synchronousRectifier, phaseCount>1, peakOfLinePlusSectors) are
//   gated separately in run_checks and documented in VIENNA_PHASE3_PLAN.md.
//
// Deferred (throws on use, planned next):
//   - samplingStrategy=peakOfLinePlusSectors
//   - viennaII variant, backToBackMosfet / singleMosfetIn4DiodeBridge
//     switch types, synchronousRectifier, phaseCount > 1.
//   - Neutral-point voltage-ripple modelling (NP balancing controller).
//   - DC-bus capacitor sizing.
//   - SPICE netlist emission for the full 3-phase circuit.
//
// REFERENCES
//   [1] Kolar & Zach, PCIM 1994 (the original Vienna paper).
//   [2] TI TIDUCJ0B (TIDM-1000 user guide); TI TIDA-010257.
//   [3] ST UM2975 (STDES-VIENNARECT 15 kW).
//   [4] src/converter_models/VIENNA_PLAN.md - full design rationale.
// =====================================================================
#include "converter_models/Vienna.h"
#include "converter_models/Topology.h"
#include "support/Utils.h"
#include "processors/NgspiceRunner.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <string>

namespace OpenMagnetics {

// =====================================================================
// Constructors
// =====================================================================
// Validate the spec shape up-front. Without these guards, mis-shaped specs
// (e.g. a DimensionWithTolerance object passed where the schema requires a
// bare number) escape from deep inside the quicktype-generated from_json
// machinery as `[json.exception.type_error.302] type must be number, but
// is object` with no indication of which field is at fault. Field types
// per MAS/schemas/inputs/topologies/vienna.json.
static void validate_vienna_spec_shape(const json& j) {
    if (!j.contains("lineToLineVoltage")) {
        throw std::runtime_error("Vienna: 'lineToLineVoltage' is required");
    }
    if (!j.at("lineToLineVoltage").is_object()) {
        throw std::runtime_error(
            "Vienna: 'lineToLineVoltage' must be a DimensionWithTolerance "
            "object (e.g. {\"nominal\": 400.0}), not a bare number");
    }
    if (!j.contains("outputDcVoltage")) {
        throw std::runtime_error("Vienna: 'outputDcVoltage' is required");
    }
    if (!j.at("outputDcVoltage").is_number()) {
        throw std::runtime_error(
            "Vienna: 'outputDcVoltage' must be a bare number (e.g. 800.0), "
            "not a DimensionWithTolerance object");
    }
    if (!j.contains("switchingFrequency")) {
        throw std::runtime_error("Vienna: 'switchingFrequency' is required");
    }
    if (!j.at("switchingFrequency").is_number()) {
        throw std::runtime_error(
            "Vienna: 'switchingFrequency' must be a bare number");
    }
    if (!j.contains("operatingPoints")) {
        throw std::runtime_error("Vienna: 'operatingPoints' is required");
    }
    if (!j.at("operatingPoints").is_array() || j.at("operatingPoints").empty()) {
        throw std::runtime_error(
            "Vienna: 'operatingPoints' must be a non-empty array");
    }
}

Vienna::Vienna(const json& j) {
    validate_vienna_spec_shape(j);
    from_json(j, *static_cast<ViennaRectifier*>(this));
}

AdvancedVienna::AdvancedVienna(const json& j) {
    validate_vienna_spec_shape(j);
    from_json(j, *this);
}

// =====================================================================
// Static analytical helpers
// =====================================================================
double Vienna::compute_phase_peak_voltage(double V_LL_rms) {
    if (V_LL_rms <= 0)
        throw std::runtime_error("Vienna: lineToLineVoltage must be > 0");
    return std::sqrt(2.0) * V_LL_rms / std::sqrt(3.0);
}

double Vienna::compute_modulation_index(double V_phase_peak, double Vdc) {
    if (Vdc <= 0)
        throw std::runtime_error("Vienna: outputDcVoltage must be > 0");
    return V_phase_peak / (Vdc / 2.0);
}

double Vienna::compute_line_peak_current(double P, double V_phase_rms,
                                         double eff, double pf) {
    if (P <= 0)
        throw std::runtime_error("Vienna: input power must be > 0");
    if (V_phase_rms <= 0)
        throw std::runtime_error("Vienna: phase RMS voltage must be > 0");
    if (eff <= 0 || eff > 1.0)
        throw std::runtime_error("Vienna: efficiency must be in (0, 1]");
    if (pf <= 0 || pf > 1.0)
        throw std::runtime_error("Vienna: power factor must be in (0, 1]");
    return std::sqrt(2.0) * P / (3.0 * V_phase_rms * eff * pf);
}

double Vienna::compute_inductor_for_ripple(double V_phase_peak, double M,
                                           double Fsw, double DeltaI_pp_target) {
    if (Fsw <= 0)
        throw std::runtime_error("Vienna: switchingFrequency must be > 0");
    if (DeltaI_pp_target <= 0)
        throw std::runtime_error("Vienna: ripple target must be > 0");
    if (M <= 0 || M > 1.0)
        throw std::runtime_error("Vienna: modulation index out of (0, 1]");
    // Duty at peak-of-line = 1 - M (boost duty for V_phase_peak -> Vdc/2).
    double duty_at_peak = 1.0 - M;
    return V_phase_peak * duty_at_peak / (DeltaI_pp_target * Fsw);
}

double Vienna::compute_switch_rms(double I_pk, double M) {
    // Hartmann ETH dissertation 19755 (2011) §3.2, single-switch Vienna I at
    // unity power factor. Derived from
    //   I_sw_rms^2 = I_pk^2 * (1/(2*pi)) * integral_0^pi sin^2(t)*(1 - M*sin(t)) dt
    //              = I_pk^2 * (1/4 - 2*M/(3*pi))
    // Positive over the entire physical range (M <= 1, gated by over-modulation
    // check). NO silent clamp: throws if the closed form goes outside its domain.
    if (I_pk < 0)
        throw std::runtime_error("Vienna: I_pk must be >= 0");
    if (M < 0 || M > 1.0)
        throw std::runtime_error("Vienna: modulation index out of [0, 1]");
    double inner = 0.25 - 2.0 * M / (3.0 * M_PI);
    if (inner < 0)
        throw std::runtime_error(
            "Vienna: switch-RMS closed form yields negative argument at M=" +
            std::to_string(M) + " - outside Hartmann derivation domain");
    return I_pk * std::sqrt(inner);
}

double Vienna::compute_diode_avg(double I_pk) {
    if (I_pk < 0)
        throw std::runtime_error("Vienna: I_pk must be >= 0");
    return I_pk / M_PI;
}

double Vienna::compute_switch_rms_vienna_ii(double I_pk, double M) {
    if (I_pk < 0)
        throw std::runtime_error("Vienna II: I_pk must be >= 0");
    // Domain: M ≤ 3π/8 ≈ 1.178 — comfortably beyond the M ≤ 1 over-mod
    // gate. Inner term is positive for any valid M.
    double inner = 1.0 / 8.0 - M / (3.0 * M_PI);
    if (inner < 0)
        throw std::runtime_error(
            "Vienna II: switch-RMS closed form yields negative argument at M=" +
            std::to_string(M));
    return I_pk * std::sqrt(inner);
}

double Vienna::compute_switch_avg_vienna_ii(double I_pk, double M) {
    if (I_pk < 0)
        throw std::runtime_error("Vienna II: I_pk must be >= 0");
    // Always positive for the valid M ≤ 1 range (1/π ≈ 0.318 ≫ M/4 ≤ 0.25).
    return I_pk * (1.0 / M_PI - M / 4.0);
}

// =====================================================================
// run_checks
// =====================================================================
bool Vienna::run_checks(bool assertErrors) {
    _assertErrors = assertErrors;
    auto& ops = get_operating_points();
    if (ops.empty()) {
        if (assertErrors) throw std::runtime_error("Vienna: no operating points");
        return false;
    }
    if (get_switching_frequency() <= 0) {
        if (assertErrors) throw std::runtime_error("Vienna: switchingFrequency must be > 0");
        return false;
    }
    if (get_output_dc_voltage() <= 0) {
        if (assertErrors) throw std::runtime_error("Vienna: outputDcVoltage must be > 0");
        return false;
    }
    auto& vll = get_line_to_line_voltage();
    if (!vll.get_nominal().has_value() &&
        !(vll.get_minimum().has_value() && vll.get_maximum().has_value())) {
        if (assertErrors) throw std::runtime_error(
            "Vienna: lineToLineVoltage must specify nominal (or both min and max)");
        return false;
    }
    for (const auto& op : ops) {
        if (op.get_output_voltages().empty() || op.get_output_currents().empty()) {
            if (assertErrors) throw std::runtime_error(
                "Vienna: OP missing output voltage/current (Vdc and total Idc)");
            return false;
        }
    }
    // Variant gates — Phase 3 supports both viennaI and viennaII.
    // switchType gate fully lifted in Phase 3 Item 4. All three variants
    // are now analytically modelled:
    //   - tType                       → standard Vienna I single-switch leg
    //                                   (or Vienna II two-switch leg if
    //                                   viennaVariant=viennaII).
    //   - backToBackMosfet            → analytically identical to Vienna II
    //                                   (same physical two-switch arrangement,
    //                                   different schematic-symbol on BOM).
    //   - singleMosfetIn4DiodeBridge  → single MOSFET inside a 4-diode
    //                                   bridge per leg. Same Vienna I
    //                                   per-switch RMS / avg, but with 4
    //                                   forward diode drops in the
    //                                   conduction path (higher conduction
    //                                   loss — caller's loss model).
    if (get_synchronous_rectifier().value_or(false)) {
        if (assertErrors) throw std::runtime_error(
            "Vienna: synchronousRectifier is not modelled in Phase 1+2 "
            "(deferred to Phase 3).");
        return false;
    }
    if (get_phase_count().value_or(1) > 1) {
        if (assertErrors) throw std::runtime_error(
            "Vienna: phaseCount > 1 (interleaved channels) is not supported "
            "in Phase 1+2 (deferred to Phase 3).");
        return false;
    }
    // Phase 3 supports all three sampling strategies: peakOfLineOnly,
    // fullLineCycle, peakOfLinePlusSectors. No gate left here.
    return true;
}

// =====================================================================
// build_line_cycle_waveform — full-line-cycle envelope + switching ripple
// =====================================================================
//
// Direct answer to the wizard's "all 3 phases look identical" report: in
// Phase 1+2 the three windings carried the SAME switching-period snapshot
// at peak-of-line, because each phase was sampled at its own worst-case
// angle and that is the same waveform by symmetry. With fullLineCycle the
// three windings carry the full 50/60 Hz line cycle, shifted by 120° per
// phase, so A/B/C trace a recognisable 3-phase pattern in the wizard plot.
//
// Sampling note: 4096 samples per line cycle gives ~2 samples per
// switching cycle at the typical 100 kHz / 50 Hz ratio. The line envelope
// is therefore exact and the switching ripple is aliased — visible as a
// reduced-amplitude high-frequency texture on top of the sinusoid. Both
// content bands survive into the FFT-based harmonic decomposition the
// MagneticAdviser runs downstream, so the per-band core-loss attribution
// still works (line component from the envelope, switching component from
// the aliased ripple band). For high-fidelity switching-period view of
// the inductor at its stress point use samplingStrategy=peakOfLineOnly.
Waveform Vienna::build_line_cycle_waveform(
    LineCycleKind kind,
    double I_pk, double V_phase_peak, double Vdc,
    double L, double Fsw, double F_line,
    double phaseOffsetRad,
    size_t numSamples)
{
    if (numSamples < 2)
        throw std::runtime_error("Vienna::build_line_cycle_waveform: numSamples must be >= 2");
    if (F_line <= 0)
        throw std::runtime_error("Vienna::build_line_cycle_waveform: F_line must be > 0");
    if (Fsw <= F_line)
        throw std::runtime_error(
            "Vienna::build_line_cycle_waveform: Fsw must be > F_line");
    if (L <= 0)
        throw std::runtime_error("Vienna::build_line_cycle_waveform: L must be > 0");

    const double T_line  = 1.0 / F_line;
    const double T_sw    = 1.0 / Fsw;
    const double omega   = 2.0 * M_PI * F_line;
    const double Vhalf   = Vdc / 2.0;

    std::vector<double> time(numSamples);
    std::vector<double> data(numSamples);
    for (size_t i = 0; i < numSamples; ++i) {
        double t = static_cast<double>(i) / static_cast<double>(numSamples - 1) * T_line;
        time[i]  = t;

        double theta     = omega * t + phaseOffsetRad;
        double sinTheta  = std::sin(theta);
        double Vphase_t  = V_phase_peak * sinTheta;
        double Iavg_t    = I_pk * sinTheta;

        // Per-angle boost duty (1 − |Vphase|/Vhalf). Floors at 0 in the
        // (rare) over-modulated edge where |Vphase| > Vhalf; the M ≤ 1 gate
        // in process_design_requirements should prevent that, but be safe.
        double dutyAbs = 1.0 - std::abs(Vphase_t) / Vhalf;
        if (dutyAbs < 0.0) dutyAbs = 0.0;
        if (dutyAbs > 1.0) dutyAbs = 1.0;

        // Local switching-period ripple peak-to-peak.
        double dI_pp = std::abs(Vphase_t) * dutyAbs * T_sw / L;

        // Sub-sampled triangular ripple: tri(2π·Fsw·t) in [-1,+1], scaled
        // by dI_pp/2. Captures sign correctly so the aliased high-freq
        // band appears in the FFT.
        double tri = 0.0;
        {
            double swPhase = std::fmod(Fsw * t, 1.0);            // [0,1)
            tri = (swPhase < 0.5) ? (4.0 * swPhase - 1.0)
                                  : (3.0 - 4.0 * swPhase);
        }

        if (kind == LineCycleKind::CURRENT) {
            data[i] = Iavg_t + 0.5 * dI_pp * tri;
        }
        else {
            // Inductor voltage = +V_phase during the switch-ON segment,
            // V_phase − sign(Vphase)·Vhalf during OFF (so the cycle average
            // matches the boost equilibrium and the high-freq band is
            // visible). Use the same triangular phase to alternate.
            double V_on  = Vphase_t;
            double V_off = Vphase_t - ((Vphase_t >= 0) ? Vhalf : -Vhalf);
            data[i] = (tri >= 0) ? V_on : V_off;
        }
    }

    Waveform wf;
    wf.set_data(data);
    wf.set_time(time);
    return wf;
}

// =====================================================================
// process_design_requirements
//
// Computes:
//   M  = V_phase_peak / (Vdc/2)            modulation index (must be <= 1)
//   I_pk per phase from worst-case OP power
//   L  from currentRippleRatio * I_pk peak-to-peak target
// =====================================================================
DesignRequirements Vienna::process_design_requirements() {
    if (!run_checks(true)) {
        // run_checks(true) throws - this line is for the compiler.
        throw std::runtime_error("Vienna: run_checks failed");
    }

    auto& vll = get_line_to_line_voltage();
    double V_LL = vll.get_nominal().value_or(
        (vll.get_minimum().value_or(0) + vll.get_maximum().value_or(0)) / 2.0);
    if (V_LL <= 0)
        throw std::runtime_error("Vienna: lineToLineVoltage nominal must be > 0");

    double Vdc = get_output_dc_voltage();
    double Fsw = get_switching_frequency();
    double Lf  = get_line_frequency().value_or(50.0);
    if (Lf <= 0)
        throw std::runtime_error("Vienna: lineFrequency must be > 0");
    double eff = get_efficiency().value_or(0.97);
    double pf  = get_power_factor().value_or(0.99);
    double rippleRatio = get_current_ripple_ratio().value_or(0.25);
    if (rippleRatio <= 0 || rippleRatio >= 1.0)
        throw std::runtime_error("Vienna: currentRippleRatio must be in (0, 1)");

    double V_phase_peak = compute_phase_peak_voltage(V_LL);
    double V_phase_rms  = V_phase_peak / std::sqrt(2.0);
    if (Vdc <= std::sqrt(2.0) * V_LL)
        throw std::runtime_error(
            "Vienna: outputDcVoltage (" + std::to_string(Vdc) +
            " V) must exceed sqrt(2)*V_LL (" + std::to_string(std::sqrt(2.0) * V_LL) +
            " V) for boost-PFC operation.");

    double M = compute_modulation_index(V_phase_peak, Vdc);
    if (M > 1.0)
        throw std::runtime_error(
            "Vienna: modulation index M=" + std::to_string(M) +
            " > 1 (over-modulation). Increase outputDcVoltage above sqrt(2)*V_LL.");
    computedModulationIndex     = M;
    computedSwitchVoltageStress = Vdc / 2.0;

    // Worst-case I_pk over the operating-point list (largest output power).
    auto& ops = get_operating_points();
    double P_max = 0;
    for (const auto& op : ops) {
        double Vout = op.get_output_voltages()[0];
        double Iout = op.get_output_currents()[0];
        if (Vout <= 0 || Iout <= 0)
            throw std::runtime_error("Vienna: OP output voltage/current must be > 0");
        double P = Vout * Iout;
        if (P > P_max) P_max = P;
    }
    if (P_max <= 0)
        throw std::runtime_error("Vienna: worst-case OP power is non-positive");

    double I_pk = compute_line_peak_current(P_max, V_phase_rms, eff, pf);
    computedLinePeakCurrent = I_pk;

    // L from peak-to-peak ripple target (ripple expressed as fraction of I_pk).
    double DeltaI_pp_target = rippleRatio * I_pk;
    double L = compute_inductor_for_ripple(V_phase_peak, M, Fsw, DeltaI_pp_target);

    // User overrides take precedence (validation against published designs).
    if (userBoostInductance > 0) L = userBoostInductance;
    if (L <= 0)
        throw std::runtime_error("Vienna: derived boost inductance non-positive");
    computedBoostInductance = L;

    // ---- Build DesignRequirements ----
    DesignRequirements designRequirements;

    // Vienna has NO transformer turns ratio - the boost inductors are
    // single-winding (primary only). Emit a single 1:1 entry to keep the
    // downstream API consistent.
    DimensionWithTolerance nTol;
    nTol.set_nominal(1.0);
    designRequirements.get_mutable_turns_ratios().push_back(nTol);

    DimensionWithTolerance lTol;
    lTol.set_nominal(roundFloat(L, 10));
    designRequirements.set_magnetizing_inductance(lTol);

    // Three identical boost inductors, all on the primary side. The MKF
    // isolation framework expects per-output sides; for a non-isolated PFC
    // we emit a single primary side.
    std::vector<IsolationSide> isolationSides;
    isolationSides.push_back(get_isolation_side_from_index(0));
    designRequirements.set_isolation_sides(isolationSides);
    designRequirements.set_topology(Topologies::VIENNA_RECTIFIER_CONVERTER);

    return designRequirements;
}

// =====================================================================
// process_operating_points
//
// Dispatch on samplingStrategy:
//   peakOfLineOnly (default)    → 1 OP per input × peak-of-line snapshot
//   fullLineCycle (Phase 3)     → 1 OP per input × full line-cycle envelope
//   peakOfLinePlusSectors (P-3) → 6 OPs per input (DPWM sector mid-points)
//
// All three paths emit 3 windings per OP (Phase A/B/C). The first two
// match the legacy single-OP-per-input convention; the third produces
// six OPs per input so the adviser sizes worst-case across DPWM sectors.
// =====================================================================
std::vector<OperatingPoint> Vienna::process_operating_points(
    const std::vector<double>& /*turnsRatios*/, double /*magnetizingInductance*/)
{
    if (computedBoostInductance <= 0)
        process_design_requirements();

    const auto strategy =
        get_sampling_strategy().value_or(ViennaSamplingStrategy::PEAK_OF_LINE_ONLY);

    std::vector<OperatingPoint> result;
    auto& ops = get_operating_points();

    if (strategy == ViennaSamplingStrategy::PEAK_OF_LINE_PLUS_SECTORS) {
        // Six DPWM sector mid-points across one line cycle. Phase A's sector
        // angle θ_sec; phases B and C are at θ_sec ± 2π/3 inside each OP, so
        // every OP captures a single instant in the balanced-3-phase grid
        // and shows what each phase is doing at that instant.
        static const double sectorAngles[6] = {
            M_PI / 6.0,
            M_PI / 2.0,            // = peak-of-line for Phase A
            5.0 * M_PI / 6.0,
            7.0 * M_PI / 6.0,
            3.0 * M_PI / 2.0,
            11.0 * M_PI / 6.0,
        };
        static const char* sectorLabels[6] = {
            "30°", "90°", "150°", "210°", "270°", "330°",
        };
        for (size_t i = 0; i < ops.size(); ++i) {
            // Set the per-op diagnostic side-effects (lastSwitchRmsCurrent,
            // lastInductorPeakCurrent, etc.) by running the peak-of-line
            // solver once. emit_switching_period_op_at_line_angle does NOT
            // touch the last* fields, so without this call the sector
            // path would leave them stale from the previous input OP.
            (void) process_operating_point_for_input_voltage(ops[i]);

            for (size_t s = 0; s < 6; ++s) {
                auto op = emit_switching_period_op_at_line_angle(ops[i], sectorAngles[s]);
                op.set_name("OP " + std::to_string(i)
                             + " sec " + sectorLabels[s]);
                result.push_back(op);
            }
        }
        return result;
    }

    // peakOfLineOnly and fullLineCycle share the existing per-OP solver.
    for (size_t i = 0; i < ops.size(); ++i) {
        auto op = process_operating_point_for_input_voltage(ops[i]);
        op.set_name("OP " + std::to_string(i));
        result.push_back(op);
    }
    return result;
}

std::vector<OperatingPoint> Vienna::process_operating_points(Magnetic /*magnetic*/) {
    auto req = process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios())
        turnsRatios.push_back(resolve_dimensional_values(tr));
    double L = resolve_dimensional_values(req.get_magnetizing_inductance());
    return process_operating_points(turnsRatios, L);
}

// =====================================================================
// process_operating_point_for_input_voltage
//
// Per-OP solver. Samples ONLY at peak-of-line (Phase 1+2 scope).
//
// At each phase's peak-of-line, the boost inductor sees:
//   - Average current = I_pk
//   - Duty = 1 - M (boost from V_phase_peak to Vdc/2)
//   - V_L during ON  = +V_phase_peak
//   - V_L during OFF = V_phase_peak - Vdc/2 (negative)
//   - Peak-to-peak ripple = V_phase_peak * d * Tsw / L
//
// We synthesise three triangular waveforms (one per phase) at the same
// peak-of-line operating point. They are by design identical; the
// magnetic adviser uses any of them to size one inductor.
// =====================================================================
OperatingPoint Vienna::process_operating_point_for_input_voltage(
    const TopologyExcitation& viennaOpPoint)
{
    OperatingPoint operatingPoint;

    if (viennaOpPoint.get_output_voltages().empty() ||
        viennaOpPoint.get_output_currents().empty())
        throw std::runtime_error("Vienna: OP missing output voltage/current");

    double Fsw = get_switching_frequency();
    if (Fsw <= 0)
        throw std::runtime_error("Vienna: switchingFrequency must be > 0");

    auto& vll = get_line_to_line_voltage();
    double V_LL = vll.get_nominal().value_or(
        (vll.get_minimum().value_or(0) + vll.get_maximum().value_or(0)) / 2.0);
    double Vdc = get_output_dc_voltage();
    double V_phase_peak = compute_phase_peak_voltage(V_LL);
    double V_phase_rms  = V_phase_peak / std::sqrt(2.0);
    double M = compute_modulation_index(V_phase_peak, Vdc);
    if (M > 1.0)
        throw std::runtime_error(
            "Vienna: modulation index M=" + std::to_string(M) +
            " > 1 (over-modulation). Increase outputDcVoltage above sqrt(2)*V_LL.");

    double eff = get_efficiency().value_or(0.97);
    double pf  = get_power_factor().value_or(0.99);

    double Vout = viennaOpPoint.get_output_voltages()[0];
    double Iout = viennaOpPoint.get_output_currents()[0];
    double P    = Vout * Iout;
    if (P <= 0)
        throw std::runtime_error("Vienna: OP output power must be > 0");

    double I_pk = compute_line_peak_current(P, V_phase_rms, eff, pf);
    double L    = computedBoostInductance;
    if (L <= 0)
        throw std::runtime_error("Vienna: boost inductance not initialised");

    double duty_at_peak = 1.0 - M;
    double V_L_on       = V_phase_peak;          // switch closed: full phase voltage across L
    double V_L_off      = V_phase_peak - Vdc / 2.0;  // switch open: rectifier conducts to bus
    double V_L_pp       = V_L_on - V_L_off;      // = Vdc/2
    double Tsw          = 1.0 / Fsw;
    double DeltaI_pp    = V_L_on * duty_at_peak * Tsw / L;

    // Diagnostics
    lastInductorPeakCurrent      = I_pk + DeltaI_pp / 2.0;
    lastInductorRipplePeakToPeak = DeltaI_pp;
    lastDutyAtPeak               = duty_at_peak;
    lastSwitchVoltageStress      = Vdc / 2.0;
    // Variant-aware switch / diode diagnostics.
    //
    //   viennaVariant=viennaII       → 2 switches/leg, each conducts half
    //                                  the duty → use Vienna II forms.
    //   switchType=backToBackMosfet  → physically identical to Vienna II
    //                                  (two switches per leg). Treat as
    //                                  Vienna II for the diagnostics.
    //   switchType=singleMosfetIn4DiodeBridge → single switch per leg, both
    //                                  polarities via the bridge → Vienna I
    //                                  forms apply unchanged. The 4 diode
    //                                  drops show up in the caller's
    //                                  conduction-loss model, not in the
    //                                  per-switch RMS/avg numbers.
    const auto vvariant = get_vienna_variant().value_or(ViennaVariant::VIENNA_I);
    const auto stype    = get_switch_type().value_or(ViennaSwitchType::T_TYPE);
    const bool twoSwitchPerLeg = (vvariant == ViennaVariant::VIENNA_II) ||
                                  (stype    == ViennaSwitchType::BACK_TO_BACK_MOSFET);
    lastSwitchRmsCurrent = twoSwitchPerLeg
        ? compute_switch_rms_vienna_ii(I_pk, M)
        : compute_switch_rms(I_pk, M);
    lastDiodeAvgCurrent  = twoSwitchPerLeg
        ? compute_switch_avg_vienna_ii(I_pk, M)  // body-diode current
        : compute_diode_avg(I_pk);                // fast-rectifier current
    lastModulationIndex          = M;
    lastInputPower               = P / eff;

    // ---- Build three phase-inductor windings ----
    //
    // peakOfLineOnly (default, Phase 1+2): each winding is the switching-
    // period waveform at the worst-case stress point (its own peak-of-line);
    // all three are by definition identical. Good for adviser worst-case
    // sizing; poor for visualisation (wizard shows three identical traces).
    //
    // fullLineCycle (Phase 3): each winding is the full 50/60 Hz line cycle,
    // with the three phases offset by ±120°. Wizard sees a recognisable
    // 3-phase pattern; adviser sees both line-frequency and switching-
    // frequency content via the FFT.
    static const char* phaseNames[3] = { "Phase A", "Phase B", "Phase C" };
    const double F_line = get_line_frequency().value_or(50.0);
    const bool fullLineCycle =
        get_sampling_strategy().value_or(ViennaSamplingStrategy::PEAK_OF_LINE_ONLY)
            == ViennaSamplingStrategy::FULL_LINE_CYCLE;

    // Positive-sequence A-B-C: Phase A at 0, B at -120°, C at +120°. The
    // wizard plot then shows A leading B leading C, matching the grid
    // convention.
    static const double phaseOffsets[3] = {
        0.0,
        -2.0 * M_PI / 3.0,
        +2.0 * M_PI / 3.0,
    };

    for (int ph = 0; ph < 3; ++ph) {
        Waveform currentWaveform;
        Waveform voltageWaveform;
        double opFreq = Fsw;

        if (fullLineCycle) {
            currentWaveform = build_line_cycle_waveform(
                LineCycleKind::CURRENT, I_pk, V_phase_peak, Vdc,
                L, Fsw, F_line, phaseOffsets[ph]);
            voltageWaveform = build_line_cycle_waveform(
                LineCycleKind::VOLTAGE, I_pk, V_phase_peak, Vdc,
                L, Fsw, F_line, phaseOffsets[ph]);
            // The waveform's natural period is the line period; downstream
            // FFT analysis treats `frequency` as the fundamental of the
            // signal. Use F_line so the line-frequency component lands on
            // bin 1 of the harmonic decomposition.
            opFreq = F_line;
        }
        else {
            currentWaveform = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, DeltaI_pp, Fsw, duty_at_peak, I_pk, 0);
            voltageWaveform = Inputs::create_waveform(
                WaveformLabel::RECTANGULAR, V_L_pp, Fsw, duty_at_peak,
                V_L_off + V_L_pp / 2.0, 0);
        }

        auto excitation = complete_excitation(
            currentWaveform, voltageWaveform, opFreq, phaseNames[ph]);
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    OperatingConditions conditions;
    conditions.set_ambient_temperature(viennaOpPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    return operatingPoint;
}


// =====================================================================
// emit_switching_period_op_at_line_angle  (Phase 3, Item 2)
//
// Builds ONE switching-period operating point at a given line angle
// θ_sec for Phase A; Phases B and C are at θ_sec ± 2π/3. Each phase's
// local inductor sees its own instantaneous (V_phase, duty, ΔI_pp)
// because |sin| varies per phase at the same wall-clock instant.
//
// Used by samplingStrategy=peakOfLinePlusSectors (called 6× with the
// DPWM sector mid-points). The default peakOfLineOnly path still goes
// through process_operating_point_for_input_voltage at θ=π/2.
// =====================================================================
OperatingPoint Vienna::emit_switching_period_op_at_line_angle(
    const TopologyExcitation& viennaOpPoint,
    double sectorAngleRad)
{
    if (viennaOpPoint.get_output_voltages().empty() ||
        viennaOpPoint.get_output_currents().empty())
        throw std::runtime_error("Vienna: OP missing output voltage/current");

    double Fsw = get_switching_frequency();
    if (Fsw <= 0)
        throw std::runtime_error("Vienna: switchingFrequency must be > 0");

    auto& vll = get_line_to_line_voltage();
    double V_LL = vll.get_nominal().value_or(
        (vll.get_minimum().value_or(0) + vll.get_maximum().value_or(0)) / 2.0);
    double Vdc = get_output_dc_voltage();
    double V_phase_peak = compute_phase_peak_voltage(V_LL);
    double V_phase_rms  = V_phase_peak / std::sqrt(2.0);
    double M = compute_modulation_index(V_phase_peak, Vdc);
    if (M > 1.0)
        throw std::runtime_error(
            "Vienna: modulation index M=" + std::to_string(M) +
            " > 1 (over-modulation).");

    double eff = get_efficiency().value_or(0.97);
    double pf  = get_power_factor().value_or(0.99);

    double Vout = viennaOpPoint.get_output_voltages()[0];
    double Iout = viennaOpPoint.get_output_currents()[0];
    double P    = Vout * Iout;
    if (P <= 0)
        throw std::runtime_error("Vienna: OP output power must be > 0");

    double I_pk = compute_line_peak_current(P, V_phase_rms, eff, pf);
    double L    = computedBoostInductance;
    if (L <= 0)
        throw std::runtime_error("Vienna: boost inductance not initialised");

    const double Vhalf = Vdc / 2.0;
    const double Tsw   = 1.0 / Fsw;

    static const char*  phaseNames[3]   = { "Phase A", "Phase B", "Phase C" };
    static const double phaseOffsets[3] = {
        0.0, -2.0 * M_PI / 3.0, +2.0 * M_PI / 3.0,
    };

    OperatingPoint operatingPoint;
    for (int ph = 0; ph < 3; ++ph) {
        const double theta_ph   = sectorAngleRad + phaseOffsets[ph];
        const double sinTheta   = std::sin(theta_ph);
        const double Vphase_inst = V_phase_peak * sinTheta;
        const double Vphase_abs  = std::abs(Vphase_inst);
        const double I_avg_inst  = I_pk * sinTheta;

        // Local boost duty at this instant. Below the "off-time floor" of
        // |V_phase|/Vhalf the converter has nothing to do (sin near zero
        // crossing), so duty saturates at 0.
        double duty_inst = 1.0 - Vphase_abs / Vhalf;
        if (duty_inst < 0.0) duty_inst = 0.0;
        if (duty_inst > 1.0) duty_inst = 1.0;

        // Switching-period ripple. ΔI = V_phase·d·T_sw / L. Sign follows
        // the inductor-current sign so the triangular waveform is centred
        // on the local DC average I_avg_inst (which may itself be ±).
        const double dI_pp = Vphase_abs * duty_inst * Tsw / L;
        const double V_L_on  = Vphase_inst;
        const double V_L_off = Vphase_inst - ((Vphase_inst >= 0) ? Vhalf : -Vhalf);
        const double V_L_pp  = std::abs(V_L_on - V_L_off);
        const double V_L_off_mid = V_L_off + (V_L_on - V_L_off) / 2.0;

        Waveform currentWaveform = Inputs::create_waveform(
            WaveformLabel::TRIANGULAR, dI_pp, Fsw, duty_inst, I_avg_inst, 0);
        Waveform voltageWaveform = Inputs::create_waveform(
            WaveformLabel::RECTANGULAR, V_L_pp, Fsw, duty_inst, V_L_off_mid, 0);

        auto excitation = complete_excitation(
            currentWaveform, voltageWaveform, Fsw, phaseNames[ph]);
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    OperatingConditions conditions;
    conditions.set_ambient_temperature(viennaOpPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    return operatingPoint;
}


// =====================================================================
// AdvancedVienna - user supplies L directly.
// =====================================================================
DesignRequirements AdvancedVienna::process_design_requirements() {
    if (desiredBoostInductance.has_value() && desiredBoostInductance.value() > 0)
        set_user_boost_inductance(desiredBoostInductance.value());
    return Vienna::process_design_requirements();
}

// =====================================================================
// generate_ngspice_circuit  (Phase-1 SPICE: single-phase boost emulation)
//
// Emits ONE phase of the Vienna as a stand-alone boost converter at
// peak-of-line (line angle = pi/2 for the modelled phase). Inputs:
//   - Vphase  DC source        = V_phase_peak  (frozen at line peak)
//   - L_boost                  = computedBoostInductance
//   - Sense Vsrc Vph_sense     for inductor-current readout
//   - Voltage-controlled switch SBOOST + body-diode DSWBD (T-type leg)
//     gated at fsw, duty = 1 - M (boost duty at peak-of-line)
//   - Fast diode DBOOST from sw node to upper half-bus
//   - Cout (split-bus emulation: one cap = Vdc/2) with .ic
//   - Rload sized so I_avg(L) ~= I_pk (power-balance through one phase)
//
// The three Vienna phase inductors are identical by analytical
// symmetry; we duplicate this waveform across "Phase A/B/C" windings in
// the extractor.
// =====================================================================
std::string Vienna::generate_ngspice_circuit(
    const std::vector<double>& /*turnsRatios*/,
    double magnetizingInductance,
    size_t /*inputVoltageIndex*/,
    size_t operatingPointIndex)
{
    auto& ops = get_operating_points();
    if (ops.empty())
        throw std::runtime_error("Vienna generate_ngspice_circuit: no operating points");
    auto& vOp = ops[std::min(operatingPointIndex, ops.size() - 1)];

    double Fsw = get_switching_frequency();
    if (Fsw <= 0)
        throw std::runtime_error("Vienna generate_ngspice_circuit: switching frequency must be > 0");

    auto& vll = get_line_to_line_voltage();
    double V_LL = vll.get_nominal().value_or(
        (vll.get_minimum().value_or(0) + vll.get_maximum().value_or(0)) / 2.0);
    double Vdc  = get_output_dc_voltage();
    double V_phase_peak = compute_phase_peak_voltage(V_LL);
    double V_phase_rms  = V_phase_peak / std::sqrt(2.0);
    double M    = compute_modulation_index(V_phase_peak, Vdc);
    if (M > 1.0)
        throw std::runtime_error(
            "Vienna generate_ngspice_circuit: M=" + std::to_string(M) +
            " > 1 (over-modulation).");

    double eff = get_efficiency().value_or(0.97);
    double pf  = get_power_factor().value_or(0.99);

    double Vout_total = vOp.get_output_voltages()[0];
    double Iout_total = vOp.get_output_currents()[0];
    double P_total    = Vout_total * Iout_total;
    if (P_total <= 0)
        throw std::runtime_error("Vienna generate_ngspice_circuit: P must be > 0");

    double I_pk = compute_line_peak_current(P_total, V_phase_rms, eff, pf);
    double L    = (magnetizingInductance > 0) ? magnetizingInductance
                                              : computedBoostInductance;
    if (L <= 0)
        throw std::runtime_error("Vienna generate_ngspice_circuit: boost L not initialised");

    double period      = 1.0 / Fsw;
    double duty        = 1.0 - M;                 // boost duty at peak-of-line
    double tOn         = duty * period;
    double tEdge       = std::min(period * 0.005, 20e-9);
    if (tOn <= 2.0 * tEdge)
        throw std::runtime_error(
            "Vienna generate_ngspice_circuit: tOn=" + std::to_string(tOn) +
            " too small for ramp tEdge=" + std::to_string(tEdge));

    // Per-phase upper half-bus emulation. Average current into the cap (when
    // switch is OFF) equals I_pk*(1-d) = I_pk*M in steady state, so
    //   R_load_per_phase = (Vdc/2) / (I_pk * M).
    double V_half      = Vdc / 2.0;
    double Iload_phase = I_pk * M;
    if (Iload_phase <= 0)
        throw std::runtime_error("Vienna: Iload per phase non-positive");
    double Rload       = V_half / Iload_phase;

    int    numPeriodsTotal = numSteadyStatePeriods + numPeriodsToExtract;
    double simTime         = numPeriodsTotal * period;
    double startTime       = numSteadyStatePeriods * period;
    double maxStep         = period / 200.0;

    std::ostringstream c;
    c << "* Vienna PFC v0.1 (Phase-1 SPICE: single-phase boost emulation at peak-of-line)\n";
    c << "* Phase A at line peak (frozen DC). The three Vienna phase inductors\n";
    c << "* are identical by analytical symmetry; this waveform is replicated\n";
    c << "* across Phase A/B/C windings in the extractor.\n";
    c << "* V_LL=" << V_LL << "V  Vdc=" << Vdc << "V  Fsw=" << (Fsw/1e3) << "kHz\n";
    c << "* V_phase_peak=" << V_phase_peak << "V  M=" << M << "  d_peak=" << duty << "\n";
    c << "* L=" << (L*1e6) << "uH  I_pk=" << I_pk << "A  Rload_phase=" << Rload << " ohm\n\n";

    // DC bus probe (V_LL rail — single-phase emulation has no real 3-phase bus,
    // so we emit the V_phase_peak rail as the "input" for converter-port checks
    // and downstream Vin sanity gates use V_phase_peak, not V_LL).
    c << "Vphase vphase 0 " << V_phase_peak << "\n";
    c << "Rphase_dummy vphase 0 1Meg\n\n";

    // Boost inductor with primary current sense (IC= sets initial current)
    c << "Vph_sense vphase l_a 0\n";
    c << "Lboost l_a sw_node " << std::scientific << L << std::fixed
      << " IC=" << I_pk << "\n\n";

    // Switch model: T-type bidir leg behaves as a low-side switch from sw_node
    // to ground (neutral midpoint). Use SW1 with antiparallel ideal diode
    // (so the switch carries inductor current regardless of polarity).
    c << ".model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg\n";
    c << ".model DIDEAL D(IS=1e-12 RS=0.05)\n";
    c << ".model DBOOST D(Is=1e-8 N=0.01 RS=0.01)\n\n";

    c << "Vpwm pwm 0 PULSE(0 5 0 "
      << std::scientific << tEdge << " " << tEdge << " "
      << tOn << " " << period << std::fixed << ")\n";
    c << "Ssw sw_node 0 pwm 0 SW1\n";
    c << "Dsw_bd 0 sw_node DIDEAL\n";
    c << "Rsnub_sw sw_node 0 1k\n";
    c << "Csnub_sw sw_node 0 1n\n\n";

    // Boost diode to upper half-bus + cap + load
    c << "Dboost sw_node vdc_plus DBOOST\n";
    c << "Resr vdc_plus vdc_cap 0.05\n";
    c << "Cout vdc_cap 0 " << std::scientific << 47e-6 << std::fixed << "\n";
    c << "Rload vdc_cap 0 " << Rload << "\n\n";

    // Initial conditions: pre-charge bus to steady-state.
    c << ".ic v(vdc_cap)=" << V_half
      << " v(vdc_plus)="   << V_half << "\n\n";

    c << ".options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500\n";
    c << ".options METHOD=GEAR TRTOL=7\n\n";
    c << ".tran " << std::scientific << maxStep << " " << simTime
      << " " << startTime << " " << maxStep << " UIC\n\n";

    c << ".save v(vphase) i(Vphase) v(l_a) v(sw_node) v(vdc_plus) v(vdc_cap)"
      << " i(Vph_sense)\n\n";
    c << ".end\n";
    return c.str();
}

// =====================================================================
// simulate_and_extract_operating_points
//
// One netlist per (Vin x OP) combo. Vienna currently has a single
// lineToLineVoltage spec (no nominal/min/max sweep), so this is one
// netlist per OP. The single-phase-emulation waveform is mapped to all
// three "Phase A/B/C" windings (identical by analytical symmetry).
// =====================================================================
std::vector<OperatingPoint> Vienna::simulate_and_extract_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t numberOfPeriods)
{
    NgspiceRunner runner;
    if (!runner.is_available())
        return process_operating_points(turnsRatios, magnetizingInductance);

    int numPeriods = (numberOfPeriods > 0) ? static_cast<int>(numberOfPeriods)
                                           : get_num_periods_to_extract();
    int originalNumPeriodsToExtract = get_num_periods_to_extract();
    set_num_periods_to_extract(numPeriods);

    std::vector<OperatingPoint> result;
    auto& ops = get_operating_points();
    for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
        std::string netlist = generate_ngspice_circuit(
            turnsRatios, magnetizingInductance, 0, opIdx);
        double Fsw = get_switching_frequency();

        SimulationConfig config;
        config.frequency        = Fsw;
        config.extractOnePeriod = true;
        config.numberOfPeriods  = get_num_periods_to_extract();
        config.keepTempFiles    = false;

        auto simResult = runner.run_simulation(netlist, config);
        if (!simResult.success) {
            std::cerr << "Vienna sim failed: " << simResult.errorMessage
                      << ". Falling back to analytical." << std::endl;
            set_num_periods_to_extract(originalNumPeriodsToExtract);
            return process_operating_points(turnsRatios, magnetizingInductance);
        }

        // Map the same single-phase inductor waveform to all 3 windings.
        NgspiceRunner::WaveformNameMapping waveformMapping;
        std::vector<std::string> windingNames;
        std::vector<bool> flipCurrentSign;
        for (const char* nm : {"Phase A", "Phase B", "Phase C"}) {
            waveformMapping.push_back({{"voltage", "l_a"},
                                       {"current", "vph_sense#branch"}});
            windingNames.push_back(nm);
            flipCurrentSign.push_back(false);
        }

        OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
            simResult, waveformMapping, Fsw, windingNames,
            ops[opIdx].get_ambient_temperature(), flipCurrentSign);
        operatingPoint.set_name("OP " + std::to_string(opIdx));
        result.push_back(operatingPoint);
    }

    set_num_periods_to_extract(originalNumPeriodsToExtract);
    return result;
}

// =====================================================================
// simulate_and_extract_topology_waveforms
//
// §5.1 converter-port stream: for the single-phase emulation we report
// v(vphase) as input_voltage (the V_phase_peak rail). input_current is
// reconstructed from power balance (P_total / (eta * V_phase_peak * 3) per
// phase). output_voltage is v(vdc_cap) scaled by 2 to reconstitute the
// full Vdc bus; output_current = Vout/Rload_nom. Single-phase circuit so
// only one output stream.
// =====================================================================
std::vector<ConverterWaveforms> Vienna::simulate_and_extract_topology_waveforms(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t numberOfPeriods)
{
    NgspiceRunner runner;
    if (!runner.is_available())
        throw std::runtime_error("Vienna simulate_and_extract_topology_waveforms: ngspice not available");

    int originalNumPeriods = get_num_periods_to_extract();
    set_num_periods_to_extract(static_cast<int>(numberOfPeriods));

    std::vector<ConverterWaveforms> results;
    auto& ops = get_operating_points();
    for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
        std::string netlist = generate_ngspice_circuit(
            turnsRatios, magnetizingInductance, 0, opIdx);
        double Fsw = get_switching_frequency();

        SimulationConfig config;
        config.frequency        = Fsw;
        config.extractOnePeriod = true;
        config.numberOfPeriods  = numberOfPeriods;
        config.keepTempFiles    = false;

        auto simResult = runner.run_simulation(netlist, config);
        if (!simResult.success)
            throw std::runtime_error("Vienna simulation failed: " + simResult.errorMessage);

        std::map<std::string, size_t> nameToIndex;
        for (size_t i = 0; i < simResult.waveformNames.size(); ++i) {
            std::string lower = simResult.waveformNames[i];
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            nameToIndex[lower] = i;
        }
        auto getWf = [&](const std::string& name) -> Waveform {
            std::string lower = name;
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            auto it = nameToIndex.find(lower);
            return (it != nameToIndex.end()) ? simResult.waveforms[it->second] : Waveform();
        };

        ConverterWaveforms wf;
        wf.set_switching_frequency(Fsw);
        wf.set_operating_point_name("Vienna op. point " + std::to_string(opIdx));

        Waveform vphase = getWf("vphase");
        if (!vphase.get_data().empty()) wf.set_input_voltage(vphase);

        // Power-balance reconstruction of per-phase input current.
        double Vphase_local = 0.0;
        if (!vphase.get_data().empty()) {
            for (double v : vphase.get_data()) Vphase_local += v;
            Vphase_local /= vphase.get_data().size();
        }
        double Pout_total = ops[opIdx].get_output_voltages()[0] *
                            ops[opIdx].get_output_currents()[0];
        double effLocal   = get_efficiency().value_or(0.97);
        // Per-phase DC-equivalent input current at the frozen peak (one of three).
        const double Iin_dc = (Vphase_local > 0.0)
            ? (Pout_total / (effLocal * Vphase_local * 3.0)) : 0.0;
        Waveform iInWf = vphase;
        for (auto& v : iInWf.get_mutable_data()) v = Iin_dc;
        wf.set_input_current(iInWf);

        // Output rail: SPICE simulates upper half-bus (Vdc/2). Reconstitute
        // the full Vdc bus by scaling x2 (split-bus symmetric assumption).
        Waveform vCap = getWf("vdc_cap");
        if (!vCap.get_data().empty()) {
            Waveform vFull = vCap;
            for (auto& v : vFull.get_mutable_data()) v *= 2.0;
            wf.get_mutable_output_voltages().push_back(vFull);

            const double Vo_nom  = ops[opIdx].get_output_voltages()[0];
            const double Io_nom  = ops[opIdx].get_output_currents()[0];
            const double Rload_nom = (Io_nom > 0) ? (Vo_nom / Io_nom) : 100.0;
            Waveform ioutWf = vFull;
            for (auto& v : ioutWf.get_mutable_data()) v = v / Rload_nom;
            wf.get_mutable_output_currents().push_back(ioutWf);
        }
        results.push_back(wf);
    }
    set_num_periods_to_extract(originalNumPeriods);
    return results;
}

} // namespace OpenMagnetics
