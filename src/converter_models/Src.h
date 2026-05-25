#pragma once
// Version: 0.1 — Phase 1+2 (skeleton + FHA steady-state)
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"

namespace OpenMagnetics {
using namespace MAS;

/**
 * @brief Series Resonant Converter (SRC) — two-element Lr+Cr tank, no Lm branch.
 *
 * Inherits converter-level parameters from MAS::SeriesResonant and the Topology
 * interface, mirroring the Llc class structure.
 *
 * Phase 1+2 (this version): FHA-based steady-state only.
 *   - Steigerwald gain   M = π / [2·sqrt(1 + Q²·(Λ - 1/Λ)²)]   (half-bridge)
 *                        M = π /     sqrt(1 + Q²·(Λ - 1/Λ)²)   (full-bridge)
 *     where Λ = fsw/fr and the "M" maps Vbridge·n → Vout via the 4/π
 *     fundamental factor and 8·n²/π² FHA load reflection.
 *   - Sinusoidal iLr(t), vCr(t) waveform synthesis (above-resonance CCM).
 *   - Throws when below-resonance CCM or DCM is detected (deferred to future
 *     phases — see SRC_PLAN.md §11 Phase 3+).
 *
 * SRC has |M| ≤ 1 (step-down only). Above-resonance ZVS, below-resonance ZCS
 * (capacitive, hard-switching). Operating point is set by fsw alone; D is fixed
 * at 50 % for HB/FB (only FB phase-shift uses a D control).
 *
 * REFERENCES
 * [1] R. L. Steigerwald, "A Comparison of Half-Bridge Resonant Converter
 *     Topologies", IEEE Trans. Power Electronics, vol. 3, no. 2, pp. 174–182,
 *     Apr. 1988.
 * [2] M. K. Kazimierczuk, "Resonant Power Converters", Wiley, 2nd ed. 2011,
 *     Ch. 4 (Series Resonant).
 * [3] src/converter_models/SRC_PLAN.md — full design rationale and roadmap.
 */
class Src : public MAS::SeriesResonant, public Topology {
private:
    // Computed resonant-tank values (filled by process_design_requirements).
    double computedResonantInductance  = 0;   // Lr [H]
    double computedResonantCapacitance = 0;   // Cr [F]
    double computedResonantFrequency   = 0;   // fr = 1/(2π·sqrt(Lr·Cr)) [Hz]

    // User-supplied component overrides (skipped MAS schema; settable via C++).
    double userResonantInductance  = 0;       // 0 = unset
    double userResonantCapacitance = 0;       // 0 = unset

    // Diagnostics populated per operating point inside
    // process_operating_point_for_input_voltage().
    mutable double lastGainM         = 0.0;   // |Vout·n / (k_bridge·Vin)|, FHA
    mutable double lastNormalizedFsw = 0.0;   // Λ = fsw / fr
    mutable double lastIrPeak        = 0.0;   // primary tank current peak [A]
    mutable double lastVcrPeak       = 0.0;   // resonant cap peak voltage [V]
    mutable bool   lastIsAboveResonance = true;

    // ---- Per-OP diagnostic vectors (one entry per V_in solved in
    // process_operating_points). Cleared at top, pushed at end of each iter.
    mutable std::vector<std::string> perOpName;
    mutable std::vector<double>      perOpGainM;
    mutable std::vector<double>      perOpNormalizedFsw;
    mutable std::vector<double>      perOpIrPeak;
    mutable std::vector<double>      perOpVcrPeak;
    mutable std::vector<int>         perOpIsAboveResonance;

    // SPICE simulation tuning (mirrors Llc / Cllc).
    int numPeriodsToExtract   = 5;
    int numSteadyStatePeriods = 10;

    // Extra-component waveforms — one entry per operating point, populated by
    // process_operating_point_for_input_voltage and consumed by
    // get_extra_components_inputs. Cleared at the start of every
    // process_operating_points() call.
    mutable std::vector<Waveform> extraCapVoltageWaveforms;   // Vcr full-period
    mutable std::vector<Waveform> extraCapCurrentWaveforms;   // ILr full-period
    mutable std::vector<Waveform> extraIndVoltageWaveforms;   // VLr full-period
    mutable std::vector<Waveform> extraIndCurrentWaveforms;   // ILr full-period
    // CURRENT_DOUBLER: per output, two filter inductors Lo1/Lo2 each carrying
    // ≈ Iout/2 DC plus ±I_sec/4 ripple. Cleared on every process_operating_points.
    mutable std::vector<Waveform> extraLoCurrentWaveforms;    // ILo1 (CD)
    mutable std::vector<Waveform> extraLoVoltageWaveforms;    // VLo1 (CD)
    mutable std::vector<Waveform> extraLo2CurrentWaveforms;   // ILo2 (CD)
    mutable std::vector<Waveform> extraLo2VoltageWaveforms;   // VLo2 (CD)

public:
    bool _assertErrors = false;

    Src(const json& j);
    Src() {};

    bool is_bridge_topology() const override { return true; }

    MAS::Topologies topology_kind() const override {
        return MAS::Topologies::SERIES_RESONANT_CONVERTER;
    }

    // ── Computed-tank accessors ─────────────────────────────────────────────
    double get_computed_resonant_inductance()  const { return computedResonantInductance; }
    double get_computed_resonant_capacitance() const { return computedResonantCapacitance; }
    double get_computed_resonant_frequency()   const { return computedResonantFrequency; }

    // ── User overrides for Lr/Cr ────────────────────────────────────────────
    void set_user_resonant_inductance(double v)  { userResonantInductance  = v; }
    void set_user_resonant_capacitance(double v) { userResonantCapacitance = v; }
    double get_user_resonant_inductance()  const { return userResonantInductance; }
    double get_user_resonant_capacitance() const { return userResonantCapacitance; }

    // ── Diagnostics ─────────────────────────────────────────────────────────
    /** FHA gain |M| = (Vout·n)/(k_bridge·Vin) at last solved operating point. */
    double get_last_gain() const { return lastGainM; }
    /** Normalized switching frequency Λ = fsw/fr at last solved operating point. */
    double get_last_normalized_fsw() const { return lastNormalizedFsw; }
    /** Peak resonant-tank (primary winding) current at last solved operating point. */
    double get_last_ir_peak() const { return lastIrPeak; }
    /** Peak resonant-cap voltage at last solved operating point. */
    double get_last_vcr_peak() const { return lastVcrPeak; }
    /** True iff last operating point was above resonance (Λ ≥ 1). */
    bool   get_last_is_above_resonance() const { return lastIsAboveResonance; }

    const std::vector<std::string>& get_per_op_name()                 const { return perOpName; }
    const std::vector<double>&      get_per_op_gain_m()               const { return perOpGainM; }
    const std::vector<double>&      get_per_op_normalized_fsw()       const { return perOpNormalizedFsw; }
    const std::vector<double>&      get_per_op_ir_peak()              const { return perOpIrPeak; }
    const std::vector<double>&      get_per_op_vcr_peak()             const { return perOpVcrPeak; }
    const std::vector<int>&         get_per_op_is_above_resonance()   const { return perOpIsAboveResonance; }

    // ── Helpers ─────────────────────────────────────────────────────────────
    /** Returns 0.5 for half-bridge, 1.0 for full-bridge / phase-shift FB.
     *  k_bridge · Vin = peak fundamental of bridge midpoint square-wave / (4/π). */
    double get_bridge_voltage_factor() const;

    /** Returns user-provided fr if set; else 1/(2π·√(Lr·Cr)) once tank computed;
     *  else geometric mean of fsw_min, fsw_max as initial seed. Throws if none
     *  available. */
    double get_effective_resonant_frequency() const;

    /** Effective secondary rectifier type. Defaults to FULL_BRIDGE_DIODE when
     *  rectifierType is unset in the schema (matches Phase-1/2 behaviour). */
    SrcRectifierType get_effective_rectifier_type() const {
        return get_rectifier_type().value_or(SrcRectifierType::FULL_BRIDGE_DIODE);
    }

    /** True iff the effective secondary rectifier is center-tapped (CT). */
    bool is_center_tapped() const {
        return get_effective_rectifier_type() == SrcRectifierType::CENTER_TAPPED_DIODE;
    }

    /** Number of secondary windings produced per output. CT splits each output
     *  into two half-windings; FB and CD use a single secondary winding. */
    size_t windings_per_output() const { return is_center_tapped() ? 2 : 1; }

    // ── Topology interface ──────────────────────────────────────────────────
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;

    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    /** FHA-based steady-state per operating point. Above-resonance CCM only;
     *  throws on Λ < 1 (ZCS/below-resonance) or DCM. */
    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const TopologyExcitation& srcOpPoint,
        const std::vector<double>& turnsRatios);

    // ── SPICE simulation ────────────────────────────────────────────────────
    int  get_num_periods_to_extract()   const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int v)    { numPeriodsToExtract = v; }
    int  get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int v)  { numSteadyStatePeriods = v; }

    /** Emit a self-contained ngspice deck for one (Vin, OP) combination.
     *  Behavioural-PULSE bridge (no SW1/body-diode model); full-bridge diode
     *  rectifier; separate Lr (no "integrated" branch). HB and FB bridges
     *  supported. Mirrors Llc::generate_ngspice_circuit. */
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex   = 0,
        size_t operatingPointIndex = 0);

    /** Run ngspice across (Vin × OP) and return one OperatingPoint per cell. */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 0);

    /** Run ngspice and emit ConverterWaveforms for port-level diagnostics
     *  (input_voltage, input_current, output_voltages, output_currents). */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);

    /** Extra-component inputs: SRC has TWO discrete components beyond the
     *  transformer — the resonant capacitor Cr and the series inductor Lr.
     *  Returns a vector of length 2: { CAS::Inputs(resonantCapacitor),
     *  MAS::Inputs(seriesInductor) }. The Lr Inputs is sized to either the
     *  full computedResonantInductance (mode=IDEAL) or to the residual
     *  Lr_external = max(Lr - Llk, 0) once a physical magnetic is available
     *  (mode=REAL). Requires a prior process_operating_points() call so
     *  extraCap and extraInd waveforms are populated. */
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;
};


// ────────────────────────────────────────────────────────────────────────────
// AdvancedSrc: user supplies desired turns ratios & resonant components.
// ────────────────────────────────────────────────────────────────────────────
class AdvancedSrc : public Src {
private:
    std::vector<double>      desiredTurnsRatios;
    std::optional<double>    desiredResonantInductance;
    std::optional<double>    desiredResonantCapacitance;

public:
    AdvancedSrc() = default;
    ~AdvancedSrc() = default;

    AdvancedSrc(const json& j);

    DesignRequirements process_design_requirements() override;

    const std::vector<double>& get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double>& v) { desiredTurnsRatios = v; }

    std::optional<double> get_desired_resonant_inductance() const { return desiredResonantInductance; }
    void set_desired_resonant_inductance(std::optional<double> v) { desiredResonantInductance = v; }

    std::optional<double> get_desired_resonant_capacitance() const { return desiredResonantCapacitance; }
    void set_desired_resonant_capacitance(std::optional<double> v) { desiredResonantCapacitance = v; }
};


// ────────────────────────────────────────────────────────────────────────────
// JSON serialization for AdvancedSrc (Src itself uses MAS::SeriesResonant's
// generated from_json/to_json via the explicit-cast pattern in the ctor).
// ────────────────────────────────────────────────────────────────────────────
inline void from_json(const json& j, AdvancedSrc& x) {
    x.set_bridge_type(get_stack_optional<SrcBridgeType>(j, "bridgeType"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_isolated(get_stack_optional<bool>(j, "isolated"));
    x.set_max_switching_frequency(j.at("maxSwitchingFrequency").get<double>());
    x.set_min_switching_frequency(j.at("minSwitchingFrequency").get<double>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<TopologyExcitation>>());
    x.set_quality_factor(get_stack_optional<double>(j, "qualityFactor"));
    x.set_rectifier_type(get_stack_optional<SrcRectifierType>(j, "rectifierType"));
    x.set_resonant_capacitance(get_stack_optional<double>(j, "resonantCapacitance"));
    x.set_resonant_frequency(get_stack_optional<double>(j, "resonantFrequency"));
    x.set_series_inductance(get_stack_optional<double>(j, "seriesInductance"));
    x.set_use_synchronous_rectifier(get_stack_optional<bool>(j, "useSynchronousRectifier"));

    auto dtr = get_stack_optional<std::vector<double>>(j, "desiredTurnsRatios");
    if (dtr.has_value()) x.set_desired_turns_ratios(dtr.value());
    x.set_desired_resonant_inductance(get_stack_optional<double>(j, "desiredResonantInductance"));
    x.set_desired_resonant_capacitance(get_stack_optional<double>(j, "desiredResonantCapacitance"));
}

inline void to_json(json& j, const AdvancedSrc& x) {
    j = json::object();
    j["bridgeType"]                = x.get_bridge_type();
    j["efficiency"]                = x.get_efficiency();
    j["inputVoltage"]              = x.get_input_voltage();
    j["isolated"]                  = x.get_isolated();
    j["maxSwitchingFrequency"]     = x.get_max_switching_frequency();
    j["minSwitchingFrequency"]     = x.get_min_switching_frequency();
    j["operatingPoints"]           = x.get_operating_points();
    j["qualityFactor"]             = x.get_quality_factor();
    j["rectifierType"]             = x.get_rectifier_type();
    j["resonantCapacitance"]       = x.get_resonant_capacitance();
    j["resonantFrequency"]         = x.get_resonant_frequency();
    j["seriesInductance"]          = x.get_series_inductance();
    j["useSynchronousRectifier"]   = x.get_use_synchronous_rectifier();
    j["desiredTurnsRatios"]        = x.get_desired_turns_ratios();
    j["desiredResonantInductance"] = x.get_desired_resonant_inductance();
    j["desiredResonantCapacitance"]= x.get_desired_resonant_capacitance();
}

} // namespace OpenMagnetics
