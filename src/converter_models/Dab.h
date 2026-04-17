#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"

namespace OpenMagnetics {
using namespace MAS;

/**
 * @brief Dual Active Bridge (DAB) DC-DC Converter
 *
 * Inherits converter-level parameters from MAS::DualActiveBridge and
 * the Topology interface, exactly mirroring the LLC pattern.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW
 * =====================================================================
 *
 *   +--[Q1]--+--[Q3]--+             +--[Q5]--+--[Q7]--+
 *   |        |        |             |        |        |
 *  V1     bridge_a  V1           V2     bridge_c  V2
 *   |        |        |             |        |        |
 *   +--[Q2]--+--[Q4]--+             +--[Q6]--+--[Q8]--+
 *            |     L_series              |
 *            +---/\/\/---[Transformer]---+
 *                         Np : Ns
 *
 * Two full H-bridges connected through a high-frequency transformer
 * and a series inductance L (leakage or external coupling inductor).
 *
 * =====================================================================
 * MODULATION CONVENTION (Huang et al. 2018 / Rosano-Maniktala 3rd ed.)
 * =====================================================================
 *
 * A single-cell 2-bridge DAB has at most **three independent phase-shift
 * degrees of freedom**. We adopt the standard TPS naming:
 *
 *   D1 = primary-bridge intra-leg shift  (leg Q1/Q2 vs leg Q3/Q4, radians)
 *   D2 = secondary-bridge intra-leg shift (leg Q5/Q6 vs leg Q7/Q8, radians)
 *   D3 = inter-bridge (outer) shift       (primary Q1 vs secondary Q5, rad)
 *
 *        Internally the analytical equations below call D3 "phi" (φ).
 *        phi ≡ D3 — same variable, different name. The MAS schema
 *        exposes it as `innerPhaseShift3`.
 *
 * Modes (following Huang et al., Energies 11(9):2419, 2018):
 *   SPS — Single Phase Shift:   D1 = D2 = 0,  D3 free                (1 DOF)
 *   EPS — Extended Phase Shift: D1 > 0, D2 = 0, D3 free              (2 DOF)
 *   DPS — Dual Phase Shift:     D1 = D2 > 0 (symmetric), D3 free     (2 DOF)
 *   TPS — Triple Phase Shift:   D1, D2, D3 all independent           (3 DOF)
 *
 * References:
 *   [1] Huang, Wang, Wu, Liu, Wei — "General Analysis of Switching
 *       Modes in a DAB with Triple Phase Shift Modulation",
 *       Energies 11(9):2419, 2018.
 *   [2] Maniktala & Rosano — "Power Supplies A to Z", 3rd Edition,
 *       Chapter 2 (pp. 78–148). Introduces the D-based notation.
 *   [3] Calderon et al. — "Generalization of Bidirectional DAB
 *       Modulation Schemes under TPS Control", Energies 16(22):7577, 2023.
 *   [4] TI TIDA-010054 (tidues0e.pdf) — SPS reference design.
 *   [5] Demetriades PhD thesis — Chapter 6 (SPS waveform derivation).
 *
 * =====================================================================
 * KEY EQUATIONS (Single Phase Shift / SPS — D1 = D2 = 0, D3 ≡ phi)
 * =====================================================================
 *
 * Turns ratio:
 *   N = V1_nom / V2_nom  (primary-to-secondary)
 *
 * Voltage conversion ratio:
 *   d = N * V2 / V1   (d = 1 at nominal operating point)
 *
 * Power transfer (SPS):
 *   P = N * V1 * V2 * D3 * (pi - |D3|) / (2 * pi^2 * Fs * L)
 *
 * D3 for desired power:
 *   D3 = (pi/2) * (1 - sqrt(1 - 8*Fs*L*P / (N*V1*V2)))
 *
 * Inductor current (piecewise linear, referred to primary):
 *   Ibase = V1 / (2*pi*Fs*L)
 *   i1 = 0.5 * (2*D3 - (1-d)*pi) * Ibase
 *   i2 = 0.5 * (2*d*D3 + (1-d)*pi) * Ibase
 *
 * ZVS boundaries:
 *   Primary:   D3_zvs > (1 - 1/d) * pi/2
 *   Secondary: D3_zvs > (1 - d) * pi/2
 *
 * Waveform: 4 intervals per period, piecewise linear inductor current.
 *   Interval 1 (0 to t_D3):        slope = (V1 + N*V2)/L
 *   Interval 2 (t_D3 to Ts/2):     slope = (V1 - N*V2)/L
 *   Interval 3 (Ts/2 to Ts/2+t_D3):slope = -(V1 + N*V2)/L
 *   Interval 4 (Ts/2+t_D3 to Ts):  slope = -(V1 - N*V2)/L
 *
 * For EPS/DPS/TPS the waveform has up to 8 sub-intervals; see
 * compute_power_general() and dab_build_subintervals() in Dab.cpp.
 *
 * Transformer excitation (for MAS OperatingPoint):
 *   Primary:   V(t) = bipolar rectangular +/- V1
 *              I(t) = piecewise linear iL(t)
 *   Secondary: V(t) = bipolar rectangular +/- V2
 *              I(t) = N * iL(t)
 *
 * Magnetizing inductance:
 *   Lm is the transformer magnetizing inductance.
 *   Lm >> L (typically 10x-100x the series inductance).
 *   Lm determines the magnetizing current ripple.
 *   The magnetizing current is a triangular wave superimposed on iL.
 *   For MAS design: Lm_min = V1 / (4 * Fs * Delta_Im_max)
 */
class Dab : public MAS::DualActiveBridge, public Topology {
private:
    int numPeriodsToExtract = 2;
    int numSteadyStatePeriods = 3;

    // Computed design values (filled by process_design_requirements)
    double computedSeriesInductance = 0;       // L (leakage + external)
    double computedMagnetizingInductance = 0;  // Lm of transformer
    double computedDeadTime = 200e-9;          // Default 200 ns dead time
    double computedD3Rad = 0;                  // D3 (outer inter-bridge shift)
                                               // in radians, computed if the
                                               // user didn't specify one.

    // Diagnostic outputs from the per-OP analytical solver, populated inside
    // process_operating_point_for_input_voltage. These mirror the equivalent
    // accessors on the LLC model so callers can introspect ZVS margins, the
    // detected modulation regime, and the sub-interval timing structure.
    mutable int lastModulationType = 0;        // 0=SPS, 1=EPS, 2=DPS, 3=TPS
    mutable double lastZvsMarginPrimary = 0;   // D3 - D3_min,primary   (rad)
    mutable double lastZvsMarginSecondary = 0; // D3 - D3_min,secondary (rad)
    mutable double lastD3Rad = 0;              // Outer shift actually used (rad)
    mutable double lastVoltageConversionRatio = 1.0; // d = N·V₂/V₁
    mutable std::vector<double> lastSubIntervalTimes; // boundary angles in radians, [0, 2π]

public:
    bool _assertErrors = false;

    Dab(const json& j);
    Dab() {};

    // ---- Simulation tuning ----
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { numPeriodsToExtract = value; }
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { numSteadyStatePeriods = value; }

    // ---- Computed value accessors ----
    double get_computed_series_inductance() const { return computedSeriesInductance; }
    void set_computed_series_inductance(double value) { computedSeriesInductance = value; }

    double get_computed_magnetizing_inductance() const { return computedMagnetizingInductance; }
    void set_computed_magnetizing_inductance(double value) { computedMagnetizingInductance = value; }

    double get_computed_dead_time() const { return computedDeadTime; }
    void set_computed_dead_time(double value) { computedDeadTime = value; }

    /** Computed D3 (outer inter-bridge shift) in radians. Set by
     *  process_design_requirements when the user does not provide an
     *  explicit innerPhaseShift3 value. */
    double get_computed_d3_rad() const { return computedD3Rad; }
    /** Alias for get_computed_d3_rad(): SPS phase shift φ (rad). */
    double get_computed_phase_shift() const { return computedD3Rad; }

    // ---- Per-OP diagnostic accessors (populated by process_operating_point_for_input_voltage) ----
    /** Modulation type used for the last solved op point: 0=SPS, 1=EPS, 2=DPS, 3=TPS. */
    int get_last_modulation_type() const { return lastModulationType; }
    /** ZVS margin on the primary bridge: D3 - (1 - 1/d)·π/2 (rad). > 0 means ZVS achieved. */
    double get_last_zvs_margin_primary() const { return lastZvsMarginPrimary; }
    /** ZVS margin on the secondary bridge: D3 - (1 - d)·π/2 (rad). > 0 means ZVS achieved. */
    double get_last_zvs_margin_secondary() const { return lastZvsMarginSecondary; }
    /** Outer inter-bridge shift D3 actually used in the last solved op point (rad, signed). */
    double get_last_d3_rad() const { return lastD3Rad; }
    /** Voltage conversion ratio d = N·V₂/V₁ for the last solved op point. */
    double get_last_voltage_conversion_ratio() const { return lastVoltageConversionRatio; }
    /** Sub-interval boundary angles for the last solved op point (rad, [0, 2π]). */
    const std::vector<double>& get_last_sub_interval_times() const { return lastSubIntervalTimes; }

    // ---- Topology interface ----
    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;

    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    // Per-operating-point analytical waveforms
    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const DabOperatingPoint& dabOpPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // ---- DAB-specific calculations ----
    //
    // All functions below use the Huang et al. 2018 TPS convention:
    //   D1 = primary-bridge intra-leg shift (rad, from innerPhaseShift1)
    //   D2 = secondary-bridge intra-leg shift (rad, from innerPhaseShift2)
    //   D3 = inter-bridge (outer) shift (rad, from innerPhaseShift3)
    //
    // The SPS-only helpers (compute_power, compute_series_inductance,
    // compute_phase_shift, compute_switching_currents, compute_primary_rms_current,
    // check_zvs_*) take D3 directly because D1 = D2 = 0 in SPS. The "general"
    // helpers (compute_power_general, compute_phase_shift_general) take all
    // three shifts explicitly to cover EPS, DPS and TPS.

    /** SPS power transfer: P = N·V1·V2·D3·(π-|D3|)/(2π²·Fs·L). */
    static double compute_power(double V1, double V2, double N,
                                double D3, double Fs, double L);

    /** SPS: required series inductance L for desired power at given D3. */
    static double compute_series_inductance(double V1, double V2, double N,
                                            double D3, double Fs, double P);

    /** SPS: D3 (outer shift) for desired power with given inductance. */
    static double compute_phase_shift(double V1, double V2, double N,
                                      double Fs, double L, double P);

    /** Compute voltage conversion ratio d = N*V2/V1 */
    static double compute_voltage_ratio(double V1, double V2, double N);

    /** SPS: inductor current at the two switching instants i1, i2. */
    static void compute_switching_currents(double V1, double V2, double N,
                                           double D3, double Fs, double L,
                                           double& i1, double& i2);

    /** SPS: primary RMS current from (i1, i2, D3). */
    static double compute_primary_rms_current(double i1, double i2, double D3);

    /** SPS ZVS checks: D3 > (1-1/d)·π/2 on primary, D3 > (1-d)·π/2 on secondary. */
    static bool check_zvs_primary(double D3, double d);
    static bool check_zvs_secondary(double D3, double d);

    // ---- Modulation-shift accessors ----
    //
    // Read-out helpers that convert the MAS schema fields (degrees) to the
    // radians used by the analytical/SPICE engines. They also implement
    // the "DPS implicit symmetry" (D2 = D1 if the user only supplied D1).

    /** D1 — primary-bridge intra-leg shift in radians. */
    static double get_D1_rad(const DabOperatingPoint& op) {
        double deg = op.get_inner_phase_shift1().value_or(0.0);
        return deg * M_PI / 180.0;
    }

    /** D2 — secondary-bridge intra-leg shift in radians.
     *  DPS special case: if innerPhaseShift2 is absent, enforce D2 = D1
     *  (Huang 2018 / Rosano-Maniktala symmetric DPS).
     *  Other modes: returns innerPhaseShift2 if present, else 0. */
    static double get_D2_rad(const DabOperatingPoint& op) {
        if (op.get_inner_phase_shift2().has_value()) {
            return op.get_inner_phase_shift2().value() * M_PI / 180.0;
        }
        auto modTypeOpt = op.get_modulation_type();
        if (modTypeOpt.has_value() && modTypeOpt.value() == ModulationType::DPS) {
            return get_D1_rad(op);
        }
        return 0.0;
    }

    /** D3 — inter-bridge (outer) shift in radians. Returns the signed value
     *  directly from innerPhaseShift3; callers decide how to interpret "no
     *  user value" (typically fall back to computedD3Rad from the design
     *  step). Positive = forward power flow. */
    static double get_D3_rad(const DabOperatingPoint& op) {
        double deg = op.get_inner_phase_shift3().value_or(0.0);
        return deg * M_PI / 180.0;
    }

    /** Power transfer for general modulation (exact closed-form integration
     *  over the sub-intervals induced by D1, D2, D3). Used by EPS/DPS/TPS. */
    static double compute_power_general(double V1, double V2, double N,
                                         double D3, double D1, double D2,
                                         double Fs, double L);

    /** Binary-search for D3 given desired power and fixed (D1, D2). */
    static double compute_phase_shift_general(double V1, double V2, double N,
                                               double D1, double D2,
                                               double Fs, double L, double P);

    // ---- SPICE simulation ----
    //
    // The generated DAB netlist deliberately uses idealised voltage-controlled
    // switches (`SW1` model with hysteresis) rather than full MOSFET models.
    // The omitted effects are:
    //   • Finite Rds(on) variation with temperature
    //   • Coss / parasitic drain-source capacitance
    //   • Body-diode reverse recovery (trr, Qrr)
    //   • Gate-driver propagation delay
    //   • Switching-loss energy (Eon, Eoff)
    //
    // None of these effects flow through the magnetic component (the user's
    // transformer winding sees only iL, vL, and the bridge voltages, all of
    // which are dominated by the resonant tank dynamics, not the silicon).
    // MKF's purpose is magnetic-component design (winding loss, core loss,
    // flux density, leakage inductance budgeting), so the SPICE model is
    // tuned for accurate transformer waveforms — not for converter-level
    // efficiency, ZVS validation, or EMI prediction. The analytical model
    // (process_operating_point_for_input_voltage) and the SPICE model agree
    // to within ~5 % NRMSE on the primary winding current across all four
    // modulation modes (SPS / EPS / DPS / TPS), which is the relevant
    // accuracy bar for magnetic design.
    //
    // If you need realistic switching transitions (e.g. for ZVS margin
    // verification at light load), wrap the generated netlist with your own
    // MOSFET sub-circuit replacing each `S* … SW1` line.
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);

    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    /**
     * @brief Simulate and extract topology-level waveforms for converter validation
     * 
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @param numberOfPeriods Number of switching periods to simulate (default 2)
     * @return Vector of ConverterWaveforms extracted from simulation
     */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);
};


/**
 * @brief AdvancedDab: user supplies desired turns ratios & inductances directly
 */
class AdvancedDab : public Dab {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredMagnetizingInductance;
    std::optional<double> desiredSeriesInductance;

public:
    AdvancedDab() = default;
    ~AdvancedDab() = default;

    AdvancedDab(const json& j);

    Inputs process();

    const double& get_desired_magnetizing_inductance() const { return desiredMagnetizingInductance; }
    void set_desired_magnetizing_inductance(const double& value) { desiredMagnetizingInductance = value; }

    const std::vector<double>& get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double>& value) { desiredTurnsRatios = value; }

    std::optional<double> get_desired_series_inductance() const { return desiredSeriesInductance; }
    void set_desired_series_inductance(std::optional<double> value) { desiredSeriesInductance = value; }
};


// ---- JSON serialization ----
void from_json(const json& j, AdvancedDab& x);
void to_json(json& j, const AdvancedDab& x);

inline void from_json(const json& j, AdvancedDab& x) {
    // DualActiveBridge base fields
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<DabOperatingPoint>>());
    x.set_series_inductance(get_stack_optional<double>(j, "seriesInductance"));
    x.set_use_leakage_inductance(get_stack_optional<bool>(j, "useLeakageInductance"));

    // AdvancedDab extra fields
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_magnetizing_inductance(j.at("desiredMagnetizingInductance").get<double>());
    x.set_desired_series_inductance(get_stack_optional<double>(j, "desiredSeriesInductance"));
}

inline void to_json(json& j, const AdvancedDab& x) {
    j = json::object();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["operatingPoints"] = x.get_operating_points();
    j["seriesInductance"] = x.get_series_inductance();
    j["useLeakageInductance"] = x.get_use_leakage_inductance();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredMagnetizingInductance"] = x.get_desired_magnetizing_inductance();
    j["desiredSeriesInductance"] = x.get_desired_series_inductance();
}

} // namespace OpenMagnetics
