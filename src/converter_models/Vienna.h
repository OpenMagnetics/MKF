#pragma once
// Version: 0.1 - Phase 1+2 (skeleton + per-phase analytical solver at peak-of-line)
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"

namespace OpenMagnetics {
using namespace MAS;

/**
 * @brief Vienna rectifier (3-phase, 3-level, unidirectional boost-type PFC).
 *
 * Three identical line-frequency boost inductors (one per phase) carry pure-AC
 * line current at 120 deg phase shift. Switch blocking voltage = Vdc/2 (the
 * 3-level advantage). Originally proposed by Kolar & Zach, PCIM 1994; modern
 * SiC industrial implementations include TI TIDA-010257, ST STDES-VIENNARECT,
 * Microchip MSCSICPFC-REF5 and Infineon REF_11KW_PFC_SIC_QD.
 *
 * Phase 1+2 scope (this version):
 *   - Schema-driven configuration (lineToLineVoltage, outputDcVoltage,
 *     switchingFrequency, currentRippleRatio, viennaVariant, switchType,
 *     synchronousRectifier, phaseCount, samplingStrategy).
 *   - Single-channel Vienna I, T-type switches, diode rectifier only.
 *     Throws on every other variant (deferred to Phase 3+).
 *   - Sampling strategy: only `peakOfLineOnly` is implemented. Throws on
 *     `peakOfLinePlusSectors` and `fullLineCycle` (deferred to Phase 3+).
 *   - process_operating_points returns ONE OperatingPoint with THREE
 *     "windings" (Phase A inductor, Phase B inductor, Phase C inductor).
 *     Each winding represents the phase inductor at its OWN peak-of-line
 *     instant - the worst-case sizing point. The three inductors are
 *     identical by design (balanced 3-phase grid).
 *
 * KEY EQUATIONS (Kolar 1994; TI TIDUCJ0B; ST UM2975)
 *   V_phase_peak = sqrt(2) * V_LL / sqrt(3)
 *   M            = V_phase_peak / (Vdc / 2)         (modulation index, must be <= 1)
 *   I_pk         = sqrt(2) * P / (3 * V_phase_rms * eta * PF)
 *   d(peak)      = 1 - M                             (per-phase duty at peak-of-line)
 *   DeltaI_L     = V_phase_peak * d(peak) / (L * Fsw)   (peak-to-peak inductor ripple)
 *   L            = V_phase_peak * (1 - M) / (DeltaI_L_target * Fsw)
 *   Vsw          = Vdc / 2                          (switch blocking voltage, 3-level)
 *   I_sw_rms     = I_pk * sqrt(1/4 - 2*M/(3*pi))    [Hartmann ETH 19755 (2011)]
 *   I_D_avg      = I_pk / pi                         (per fast rectifier diode)
 *
 * REFERENCES
 *   [1] J. W. Kolar, F. C. Zach, "A Novel Three-Phase Three-Switch
 *       Three-Level Unity Power Factor PWM Rectifier", PCIM 1994.
 *   [2] J. W. Kolar, "Vienna Rectifier and Beyond", APEC 2018 plenary
 *       (ETH-PES).
 *   [3] T. Friedli, J. W. Kolar, "The Essence of Three-Phase PFC Rectifier
 *       Systems Part II", IEEE TPEL.
 *   [4] TI TIDUCJ0B (TIDM-1000 user guide); TI TIDA-010257 (10 kW Vienna
 *       PFC); ST UM2975 (STDES-VIENNARECT 15 kW); Microchip DS50002952B.
 *   [5] src/converter_models/VIENNA_PLAN.md - full design rationale.
 *
 * ACCURACY DISCLAIMER
 *   Phase 1+2 samples ONLY at peak-of-line. A full line-cycle envelope
 *   (THD validation, sector-weighted RMS, neutral-point ripple modelling)
 *   is deferred to Phase 3+. Three-phase grid is assumed balanced.
 *   Neutral-point voltage is assumed balanced by an outer control loop.
 *
 * DISAMBIGUATION
 *   Vienna vs 2-level 3-phase boost rectifier: 2-level uses 6 IGBTs +
 *     6 diodes, switch V-stress = Vdc (full bus). Different topology.
 *   Vienna vs 3L-ANPC: ANPC is bidirectional (V2G fast charging). Vienna
 *     is intrinsically unidirectional (diode bridge).
 *   Vienna I vs Vienna II: single switch per leg vs two switches per leg.
 *     Same topology, covered by `viennaVariant` flag. Phase 1+2 supports
 *     viennaI only.
 *   Vienna vs single-phase totem-pole PFC: 3-phase vs 1-phase. Single-
 *     phase is in PowerFactorCorrection.cpp.
 */
class Vienna : public MAS::ViennaRectifier, public Topology {
public:
    enum class Phase { A, B, C };

private:
    // Computed design values - filled by process_design_requirements().
    double computedBoostInductance     = 0;   // single value (all 3 phases identical) [H]
    double computedModulationIndex     = 0;   // M = V_phase_peak / (Vdc/2)
    double computedLinePeakCurrent     = 0;   // I_pk per phase [A]
    double computedSwitchVoltageStress = 0;   // = Vdc/2 [V]

    // User overrides for L (skipped MAS schema; settable via C++).
    double userBoostInductance = 0;           // 0 = unset

    // Diagnostics populated per operating point.
    mutable double lastInductorPeakCurrent  = 0.0;   // I_pk
    mutable double lastInductorRipplePeakToPeak = 0.0;
    mutable double lastDutyAtPeak           = 0.0;   // 1 - M
    mutable double lastSwitchVoltageStress  = 0.0;   // Vdc/2
    mutable double lastSwitchRmsCurrent     = 0.0;   // Kolar 1994
    mutable double lastDiodeAvgCurrent      = 0.0;   // I_pk/pi
    mutable double lastBoostDiodeAvgCurrent = 0.0;   // I_pk · M/4 (Phase 3)
    mutable double lastBoostDiodeRmsCurrent = 0.0;   // sync-rec MOSFET I_rms (Phase 3)
    mutable double lastModulationIndex      = 0.0;
    mutable double lastInputPower           = 0.0;

    // SPICE simulation tuning (mirrors Llc / Cllc / Src).
    int numPeriodsToExtract   = 5;
    int numSteadyStatePeriods = 20;

public:
    bool _assertErrors = false;

    Vienna(const json& j);
    Vienna() {}

    // Vienna switches connect each phase node to the neutral midpoint -
    // they route inductor current rather than setting a bridge-midpoint
    // voltage. So `is_bridge_topology()` stays false (default).

    MAS::Topology topology_kind() const override {
        return MAS::Topology::VIENNA_RECTIFIER_CONVERTER;
    }

    // A 3-phase Vienna uses three identical, separately-wound per-phase boost
    // inductors (La/Lb/Lc) — modeled here as N excitations in one operating
    // point. The from-converter adviser designs ONE of them and notes that N
    // identical inductors are required. See Topology::uses_identical_per_phase_inductors.
    bool uses_identical_per_phase_inductors() const override { return true; }

    // Computed-design accessors
    double get_computed_boost_inductance()      const { return computedBoostInductance; }
    double get_computed_modulation_index()      const { return computedModulationIndex; }
    double get_computed_line_peak_current()     const { return computedLinePeakCurrent; }
    double get_computed_switch_voltage_stress() const { return computedSwitchVoltageStress; }

    // User override (Lf)
    void   set_user_boost_inductance(double v) { userBoostInductance = v; }
    double get_user_boost_inductance() const   { return userBoostInductance; }

    // Diagnostics (last solved OP)
    double get_last_inductor_peak_current()      const { return lastInductorPeakCurrent; }
    double get_last_inductor_ripple_peak_to_peak() const { return lastInductorRipplePeakToPeak; }
    double get_last_duty_at_peak()               const { return lastDutyAtPeak; }
    double get_last_switch_voltage_stress()      const { return lastSwitchVoltageStress; }
    double get_last_switch_rms_current()         const { return lastSwitchRmsCurrent; }
    double get_last_diode_avg_current()          const { return lastDiodeAvgCurrent; }
    double get_last_boost_diode_avg_current()    const { return lastBoostDiodeAvgCurrent; }
    double get_last_boost_diode_rms_current()    const { return lastBoostDiodeRmsCurrent; }
    double get_last_modulation_index()           const { return lastModulationIndex; }
    double get_last_input_power()                const { return lastInputPower; }

    // Static analytical helpers (Phase 1+2)
    /** V_phase_peak = sqrt(2) * V_LL / sqrt(3). */
    static double compute_phase_peak_voltage(double V_LL_rms);
    /** Modulation index M = V_phase_peak / (Vdc/2). Must be <= 1 (no over-modulation). */
    static double compute_modulation_index(double V_phase_peak, double Vdc);
    /** I_pk = sqrt(2) * P / (3 * V_phase_rms * eta * PF). */
    static double compute_line_peak_current(double P, double V_phase_rms, double eff, double pf);
    /** L = V_phase_peak * (1 - M) / (DeltaI_pp_target * Fsw). */
    static double compute_inductor_for_ripple(double V_phase_peak, double M, double Fsw, double DeltaI_pp_target);
    /** Vienna I per-switch RMS via Hartmann ETH 19755 (2011) §3.2 closed form.
     *  Single switch per leg conducts both half-cycles through the diode
     *  bridge: I_sw_rms² = I_pk² · (1/4 − 2M/(3π)). */
    static double compute_switch_rms(double I_pk, double M);
    /** Per-diode average current = I_pk / pi (full-wave bridge). */
    static double compute_diode_avg(double I_pk);

    // ── Vienna II helpers (Phase 3, Item 3) ───────────────────────────────
    //
    // Vienna II splits each leg into TWO active switches (back-to-back
    // MOSFETs forming a bidirectional clamp). Each switch conducts only
    // during ONE polarity half-cycle of the phase current — exactly half
    // the duty of Vienna I's single-switch-per-leg arrangement. The fast-
    // rectifier diodes disappear; their job is taken by the body diode of
    // the OFF-polarity switch (see also synchronousRectifier, Item 5).
    //
    // Derivation (Friedli & Kolar, "Essence of 3-φ PFC Part II" §IV.B):
    //   I_sw_avg = I_pk · (1/(2π)) · ∫₀^π sin(t)·(1 − M·sin(t)) dt
    //            = I_pk · (1/π − M/4)
    //   I_sw_rms² = I_pk² · (1/(2π)) · ∫₀^π sin²(t)·(1 − M·sin(t)) dt
    //             = I_pk² · (1/8 − M/(3π))
    //
    // Each closed form is exactly HALF its Vienna I counterpart (same
    // integrand, half the integration range).

    /** Vienna II per-switch RMS — half the Vienna I duty.
     *  I_sw_rms² = I_pk² · (1/8 − M/(3π)). */
    static double compute_switch_rms_vienna_ii(double I_pk, double M);
    /** Vienna II per-switch average — half the Vienna I duty.
     *  I_sw_avg = I_pk · (1/π − M/4). */
    static double compute_switch_avg_vienna_ii(double I_pk, double M);

    // ── Boost-diode / sync-rectifier MOSFET currents (Phase 3, Item 5) ───
    //
    // The boost diode (Vienna I/II fast rectifier) conducts the inductor
    // current during the (1-d) portion of each switching cycle. With
    // d(θ) = 1 − M·|sin(θ)|, the diode duty per switching cycle is
    // M·|sin(θ)|, and at small switching-period ripple the instantaneous
    // diode current ≈ |i_L(θ)| = I_pk·|sin(θ)|.
    //
    // Integrating over the half-cycle where the diode is the active rect:
    //   I_diode_avg = I_pk · (1/(2π)) · ∫₀^π M·sin²(θ) dθ
    //               = I_pk · M / 4
    //   I_diode_rms² = I_pk² · (1/(2π)) · ∫₀^π M·sin³(θ) dθ
    //                = I_pk² · 2M / (3π)
    //
    // With synchronousRectifier=true these same closed forms apply to the
    // sync-rec MOSFET (the device that replaces the fast diode). The
    // conduction-loss model swaps from Vf·I_avg to Rds·I_rms² — Rds is
    // the caller's loss model, but the I_rms field is what they need.
    static double compute_boost_diode_avg(double I_pk, double M);
    static double compute_boost_diode_rms(double I_pk, double M);

    // Phase-3 fullLineCycle waveform builder.
    //
    // For samplingStrategy=fullLineCycle, each phase inductor sees its own
    // sinusoidal current envelope at the line frequency, offset by ±120°
    // from the other two phases. The instantaneous switching-frequency
    // ripple |ΔI(θ)| varies with the angle:
    //
    //   i_avg(t)   = I_pk · sin(ω_line·t + φ)
    //   V_phase(t) = V_peak · sin(ω_line·t + φ)
    //   d(t)       = 1 − |V_phase(t)| / (Vdc/2)
    //                (≥0 by construction since M ≤ 1)
    //   ΔI_pp(t)   = |V_phase(t)| · d(t) / (L · Fsw)
    //
    // Builds a Waveform with `numSamples` uniform time samples across one
    // line period, carrying the envelope plus a superposed triangular ripple
    // term. The envelope is exact; the ripple is sub-sampled (a typical
    // 100 kHz switching × 50 Hz line gives 2000 switching cycles per line
    // cycle, so the ripple visible at 1-4 samples/cycle is an aliased proxy
    // — adequate for plotting and for FFT-based harmonic analysis, since
    // both line-frequency and switching-frequency content remain present in
    // the spectrum).
    //
    // `phaseOffsetRad` shifts the line angle (0 = Phase A, -2π/3 = Phase B,
    // +2π/3 = Phase C in a positive-sequence A-B-C grid).
    //
    // `kind` picks the quantity: CURRENT returns inductor current; VOLTAGE
    // returns inductor terminal voltage.
    enum class LineCycleKind { CURRENT, VOLTAGE };
    static Waveform build_line_cycle_waveform(
        LineCycleKind kind,
        double I_pk, double V_phase_peak, double Vdc,
        double L, double Fsw, double F_line,
        double phaseOffsetRad,
        size_t numSamples = 4096);

    // Topology interface
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;

    /** Returns ONE OperatingPoint with THREE windings (Phase A/B/C boost inductor),
     *  each sampled at its own peak-of-line instant. */
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    /** Per-OP analytical solver. peakOfLineOnly path → switching-period
     *  triangular at sectorAngleRad = π/2; fullLineCycle path → full
     *  line-cycle envelope across all three windings.
     */
    OperatingPoint process_operating_point_for_input_voltage(
        const TopologyExcitation& viennaOpPoint);

    /** Switching-period OP at an arbitrary line angle. The angle is
     *  applied to Phase A; Phases B/C are at ±2π/3 from it. Used as the
     *  per-sector building block for samplingStrategy=peakOfLinePlusSectors
     *  (called for each of 6 sector mid-points) and as the implementation
     *  of peakOfLineOnly (sectorAngleRad = π/2).
     */
    OperatingPoint emit_switching_period_op_at_line_angle(
        const TopologyExcitation& viennaOpPoint,
        double sectorAngleRad);

    // ── SPICE simulation tuning ─────────────────────────────────────────────
    int  get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int v)  { numSteadyStatePeriods = v; }
    int  get_num_periods_to_extract()   const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int v)    { numPeriodsToExtract = v; }

    // ── SPICE buildout (Phase-1 SPICE: single-phase boost emulation) ────────
    //
    // FIXME-vienna-1 (single-phase emulation):
    //   The Phase-1 SPICE netlist models ONE phase as a stand-alone boost
    //   converter at its own peak-of-line operating point (frozen DC source
    //   = V_phase_peak, switch to neutral, diode to upper half-bus). The
    //   resulting inductor waveform is duplicated to all three windings
    //   (Phase A/B/C) by analytical symmetry — exactly matching the
    //   analytical model's assumption. Full 3-phase netlist with all three
    //   phases at frozen-DC line-angle voltages, T-type switches, split-bus
    //   caps + neutral-point modelling (per VIENNA_PLAN.md §A.6) is
    //   deferred to Phase 3+.
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex      = 0,
        size_t operatingPointIndex    = 0);

    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 1);

    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 1);
};


// ----------------------------------------------------------------------------
// AdvancedVienna: user supplies L directly (bypasses ripple-based derivation).
// ----------------------------------------------------------------------------
class AdvancedVienna : public Vienna {
private:
    std::optional<double> desiredBoostInductance;

public:
    AdvancedVienna() = default;
    ~AdvancedVienna() = default;

    AdvancedVienna(const json& j);

    DesignRequirements process_design_requirements() override;

    std::optional<double> get_desired_boost_inductance() const { return desiredBoostInductance; }
    void set_desired_boost_inductance(std::optional<double> v) { desiredBoostInductance = v; }
};


// ----------------------------------------------------------------------------
// JSON serialization for AdvancedVienna (Vienna itself uses MAS::ViennaRectifier
// generated from_json/to_json via the explicit-cast pattern in the ctor).
// ----------------------------------------------------------------------------
inline void from_json(const json& j, AdvancedVienna& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_line_frequency(get_stack_optional<double>(j, "lineFrequency"));
    x.set_line_to_line_voltage(j.at("lineToLineVoltage").get<DimensionWithTolerance>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<TopologyExcitation>>());
    x.set_output_dc_voltage(j.at("outputDcVoltage").get<double>());
    x.set_phase_count(get_stack_optional<int64_t>(j, "phaseCount"));
    x.set_power_factor(get_stack_optional<double>(j, "powerFactor"));
    x.set_sampling_strategy(get_stack_optional<ViennaSamplingStrategy>(j, "samplingStrategy"));
    x.set_switching_frequency(j.at("switchingFrequency").get<double>());
    x.set_switch_type(get_stack_optional<ViennaSwitchType>(j, "switchType"));
    x.set_synchronous_rectifier(get_stack_optional<bool>(j, "synchronousRectifier"));
    x.set_vienna_variant(get_stack_optional<ViennaVariant>(j, "viennaVariant"));

    x.set_desired_boost_inductance(get_stack_optional<double>(j, "desiredBoostInductance"));
}

inline void to_json(json& j, const AdvancedVienna& x) {
    j = json::object();
    j["currentRippleRatio"]    = x.get_current_ripple_ratio();
    j["efficiency"]            = x.get_efficiency();
    j["lineFrequency"]         = x.get_line_frequency();
    j["lineToLineVoltage"]     = x.get_line_to_line_voltage();
    j["operatingPoints"]       = x.get_operating_points();
    j["outputDcVoltage"]       = x.get_output_dc_voltage();
    j["phaseCount"]            = x.get_phase_count();
    j["powerFactor"]           = x.get_power_factor();
    j["samplingStrategy"]      = x.get_sampling_strategy();
    j["switchingFrequency"]    = x.get_switching_frequency();
    j["switchType"]            = x.get_switch_type();
    j["synchronousRectifier"]  = x.get_synchronous_rectifier();
    j["viennaVariant"]         = x.get_vienna_variant();
    j["desiredBoostInductance"]= x.get_desired_boost_inductance();
}

} // namespace OpenMagnetics
