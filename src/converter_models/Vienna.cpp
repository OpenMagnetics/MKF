// =====================================================================
// Vienna rectifier - Phase 1+2 implementation.
//
// Phase 1: skeleton (constructor, run_checks, design-requirements stub).
// Phase 2: per-phase analytical solver at peak-of-line (single sample).
//
// The three boost inductors of a balanced 3-phase Vienna are identical by
// design, and each sees its worst-case stress at its own peak-of-line
// instant. Phase 1+2 emits one OperatingPoint with three "windings" -
// Phase A / B / C inductor - each shaped as a triangular boost-inductor
// waveform at its peak-of-line operating point. The downstream magnetic
// adviser uses any of the three (they are equivalent) to size a single
// powder-core inductor; production would build three identical units.
//
// Deferred to Phase 3+ (throws on use):
//   - viennaII variant, backToBackMosfet / singleMosfetIn4DiodeBridge
//     switch types, synchronousRectifier, phaseCount > 1.
//   - peakOfLinePlusSectors / fullLineCycle sampling strategies.
//   - Neutral-point voltage-ripple modelling (NP balancing controller).
//   - DC-bus capacitor sizing.
//   - SPICE netlist emission.
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
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>

namespace OpenMagnetics {

// =====================================================================
// Constructors
// =====================================================================
Vienna::Vienna(const json& j) {
    from_json(j, *static_cast<ViennaRectifier*>(this));
}

AdvancedVienna::AdvancedVienna(const json& j) {
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
    // Variant gates - Phase 1+2 supports only the default.
    if (get_vienna_variant().has_value() &&
        get_vienna_variant().value() != ViennaVariant::VIENNA_I) {
        if (assertErrors) throw std::runtime_error(
            "Vienna: only viennaI is supported in Phase 1+2 (viennaII deferred to Phase 3).");
        return false;
    }
    if (get_switch_type().has_value() &&
        get_switch_type().value() != ViennaSwitchType::T_TYPE) {
        if (assertErrors) throw std::runtime_error(
            "Vienna: only switchType=tType is supported in Phase 1+2 "
            "(backToBackMosfet / singleMosfetIn4DiodeBridge deferred to Phase 3).");
        return false;
    }
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
    if (get_sampling_strategy().has_value() &&
        get_sampling_strategy().value() != ViennaSamplingStrategy::PEAK_OF_LINE_ONLY) {
        if (assertErrors) throw std::runtime_error(
            "Vienna: only samplingStrategy=peakOfLineOnly is implemented in "
            "Phase 1+2 (peakOfLinePlusSectors / fullLineCycle deferred to Phase 3).");
        return false;
    }
    return true;
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
// =====================================================================
std::vector<OperatingPoint> Vienna::process_operating_points(
    const std::vector<double>& /*turnsRatios*/, double /*magnetizingInductance*/)
{
    if (computedBoostInductance <= 0)
        process_design_requirements();

    std::vector<OperatingPoint> result;
    auto& ops = get_operating_points();
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
    lastSwitchRmsCurrent         = compute_switch_rms(I_pk, M);
    lastDiodeAvgCurrent          = compute_diode_avg(I_pk);
    lastModulationIndex          = M;
    lastInputPower               = P / eff;

    // ---- Build three identical phase-inductor windings ----
    // Each winding represents the phase inductor at its own peak-of-line.
    // Triangular ripple around the dc average I_pk; rectangular voltage
    // alternates +V_L_on (during d) and V_L_off (during 1-d).
    static const char* phaseNames[3] = { "Phase A", "Phase B", "Phase C" };
    for (int ph = 0; ph < 3; ++ph) {
        Waveform currentWaveform = Inputs::create_waveform(
            WaveformLabel::TRIANGULAR, DeltaI_pp, Fsw, duty_at_peak, I_pk, 0);
        Waveform voltageWaveform = Inputs::create_waveform(
            WaveformLabel::RECTANGULAR, V_L_pp, Fsw, duty_at_peak,
            V_L_off + V_L_pp / 2.0, 0);
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

} // namespace OpenMagnetics
