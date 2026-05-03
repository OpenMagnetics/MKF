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
 * @brief Phase-Shifted Half Bridge (PSHB) — Three-Level Pinheiro-Barbi NPC
 *        DC-DC Converter
 *
 * Inherits from MAS::PhaseShiftedHalfBridge and the Topology interface.
 *
 * =====================================================================
 * TOPOLOGY DISAMBIGUATION (READ THIS FIRST)
 * =====================================================================
 *
 * This model implements the **3-level NPC (Neutral Point Clamped)
 * Pinheiro-Barbi converter** (IEEE TPE 8(4) 1993):
 *   - Two stacked half-bridge legs sharing a split-capacitor input bus
 *   - 4 switches per leg + 2 clamp diodes per leg = 8 SiC + 4 diodes total
 *   - Each leg can output Vin, Vin/2 (clamped via the diodes to mid_cap),
 *     or 0 — i.e. THREE distinct levels per leg
 *   - Phase-shift modulation between the two legs controls the duty of
 *     the differential 3-level pulse: ±Vin/2 / 0
 *
 * It is **NOT** the asymmetric (complementary) half bridge of
 * Imbertson-Mohan 1993.  An AHB has only ONE leg, two complementary
 * switches with asymmetric duty (D / 1−D), a DC-blocking capacitor
 * in series with the primary, and a 2-level *asymmetric* primary
 * voltage [+(1−D)·Vin, −D·Vin].  AHB conversion ratio is
 * Vo = 2·D·(1−D)·Vin/n, which differs from this model's
 * Vo = (Vin/2)·D_eff/n.  See the planned `AsymmetricHalfBridge` class.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW (this model — Pinheiro-Barbi 1993, 3-level NPC)
 * =====================================================================
 *
 *           +Vin
 *            │
 *           C1 ─── mid_cap (= Vin/2)
 *            │       │   │
 *      ┌──Q1A──┐ ┌──Q1B──┐
 *      │   nA1 │ │   nB1 │
 *      │   ▲   │ │   ▲   │
 *      │   │ DA_h│   │ DB_h│   (clamp diode anodes at mid_cap, cathodes at nA1/nB1)
 *      │   ◄──┘ │ │   ◄──┘ │
 *      │       │ │       │
 *      ├──Q2A──┤ ├──Q2B──┤
 *      │  mid_A│ │  mid_B│   ← primary winding sits between mid_A and mid_B
 *      │       │ │       │     (with Lr in series)
 *      ├──Q3A──┤ ├──Q3B──┤
 *      │   nA2 │ │   nB2 │
 *      │   │ DA_l│   │ DB_l│   (clamp diode anodes at nA2/nB2, cathodes at mid_cap)
 *      │   ►──┐ │ │   ►──┐ │
 *      │       │ │       │
 *      └──Q4A──┘ └──Q4B──┘
 *           │
 *           0
 *
 * Per-leg states (e.g. leg-A):
 *   Vin output: Q1A + Q2A on, Q3A + Q4A off  → mid_A = Vin
 *   Vin/2 out:  Q2A + Q3A on, Q1A + Q4A off
 *               (D_high clamps nA1 to mid_cap when current flows from load
 *                to mid_A; D_low clamps nA2 to mid_cap when current flows
 *                the other way; either way mid_A ≈ Vin/2)
 *   0 output:   Q3A + Q4A on, Q1A + Q2A off  → mid_A = 0
 *
 * PWM scheme:
 *   - Inner switches Q2,Q3 run at 50 % complementary duty (Q2 first half-
 *     cycle, Q3 second half) — common to both legs.
 *   - Outer switches Q1,Q4 modulate the duty within each half:
 *       Leg-A: Q1 on during [0, t_act], Q4 on during [Thalf, Thalf+t_act]
 *       Leg-B: phase-shifted by φ_shift = (1−Deff)·Thalf
 *   - Differential mid_A − mid_B has 3-level shape ±Vin/2 / 0 with active
 *     interval = Deff·Thalf and freewheel interval = (1−Deff)·Thalf.
 *
 * =====================================================================
 * KEY EQUATIONS
 * =====================================================================
 *
 * Effective primary voltage amplitude:
 *   Vpri_pk = Vin / 2   (half of full-bridge)
 *
 * Output voltage (center-tapped rectifier):
 *   Vo = (Vin/2) * D_eff / n - Vd
 *
 * Output voltage (current doubler):
 *   Vo = (Vin/2) * D_eff / (2*n) - Vd
 *
 * Output voltage (full-bridge rectifier):
 *   Vo = (Vin/2) * D_eff / n - 2*Vd
 *
 * Turns ratio (center-tapped, for target D_eff at nominal Vin):
 *   n = (Vin_nom/2) * D_eff_nom / (Vo + Vd)
 *
 * Effective duty cycle:
 *   D_eff = phase_shift / 180  (phase_shift in degrees)
 *
 * Primary voltage waveform (3-level, same shape as PSFB but half amplitude):
 *   +(Vin/2)  during power transfer
 *   0         during freewheeling
 *   -(Vin/2)  during opposite power transfer
 *   0         during opposite freewheeling
 *
 * Primary current: same shape as PSFB, but currents are higher for
 *   the same output power because reflected load current = Io/n and
 *   n is smaller (due to half the primary voltage).
 *
 * Magnetizing inductance:
 *   Im_peak = (Vin/2) * D_eff / (4 * Fs * Lm)
 *   Lm = (Vin/2) * D_eff / (4 * Fs * Im_target)
 *
 * Output inductor:
 *   Lo = Vo * (1 - D_eff) / (Fs * dIo)
 *
 * Series inductance (ZVS assist):
 *   0.5 * Lr * Ip^2 > 0.5 * Coss * (Vin/2)^2
 *   Lr_min = Coss * (Vin/2)^2 / Ip^2
 *
 * Compared to PSFB:
 *   - Half the number of primary switches (2 vs 4)
 *   - Half the primary voltage swing (Vin/2 vs Vin)
 *   - For the same output, turns ratio n is ~half → higher primary currents
 *   - Simpler gate drive (only one leg to drive)
 *   - Typically suited for medium power (up to ~500 W)
 *   - Split capacitors must handle full primary current ripple
 */
class Pshb : public MAS::PhaseShiftedHalfBridge, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 5;

    double computedOutputInductance = 0;
    double computedSeriesInductance = 0;
    double computedMagnetizingInductance = 0;
    double computedDeadTime = 200e-9;
    double computedEffectiveDutyCycle = 0;
    double computedDiodeVoltageDrop = 0.6;
    double mosfetCoss = 200e-12;

    // Half-bridge voltage factor (split-cap → ±Vin/2)
    static constexpr double BRIDGE_VOLTAGE_FACTOR = 0.5;

    mutable double lastDutyCycleLoss = 0.0;
    mutable double lastEffectiveDutyCycle = 0.0;
    mutable double lastZvsMarginLagging = 0.0;
    mutable double lastZvsLoadThreshold = 0.0;
    mutable double lastResonantTransitionTime = 0.0;
    mutable double lastPrimaryPeakCurrent = 0.0;

    mutable std::vector<Waveform> extraLoVoltageWaveforms;
    mutable std::vector<Waveform> extraLoCurrentWaveforms;
    mutable std::vector<Waveform> extraLo2VoltageWaveforms;
    mutable std::vector<Waveform> extraLo2CurrentWaveforms;
    mutable std::vector<Waveform> extraLrVoltageWaveforms;
    mutable std::vector<Waveform> extraLrCurrentWaveforms;

public:
    bool _assertErrors = false;

    Pshb(const json& j);
    Pshb() {};

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { numPeriodsToExtract = value; }
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { numSteadyStatePeriods = value; }

    double get_computed_output_inductance() const { return computedOutputInductance; }
    void set_computed_output_inductance(double value) { computedOutputInductance = value; }
    double get_computed_series_inductance() const { return computedSeriesInductance; }
    void set_computed_series_inductance(double value) { computedSeriesInductance = value; }
    double get_computed_magnetizing_inductance() const { return computedMagnetizingInductance; }
    void set_computed_magnetizing_inductance(double value) { computedMagnetizingInductance = value; }
    double get_computed_dead_time() const { return computedDeadTime; }
    void set_computed_dead_time(double value) { computedDeadTime = value; }
    double get_computed_effective_duty_cycle() const { return computedEffectiveDutyCycle; }
    double get_bridge_voltage_factor() const { return BRIDGE_VOLTAGE_FACTOR; }

    /// Per-MOSFET output capacitance (F) used for ZVS energy and resonant
    /// transition-time predictions. Default 200 pF. Schema-less; set via setter.
    double get_mosfet_coss() const { return mosfetCoss; }
    void set_mosfet_coss(double value) { mosfetCoss = value; }

    /// Duty-cycle loss for the last solved OP, ΔD = 4·Lk·Io·Fs/(n·Vin/2).
    double get_last_duty_cycle_loss() const { return lastDutyCycleLoss; }
    /// Effective duty cycle actually delivered to the secondary.
    double get_last_effective_duty_cycle() const { return lastEffectiveDutyCycle; }
    /// Lagging-leg ZVS energy margin (J) using Vbus = Vin/2.
    double get_last_zvs_margin_lagging() const { return lastZvsMarginLagging; }
    /// Output-current ZVS load threshold (A).
    double get_last_zvs_load_threshold() const { return lastZvsLoadThreshold; }
    /// Quarter-period of the ZVS resonant tank (s).
    double get_last_resonant_transition_time() const { return lastResonantTransitionTime; }
    /// Primary peak current at the lagging-leg switching instant (A).
    double get_last_primary_peak_current() const { return lastPrimaryPeakCurrent; }

    // ---- Topology interface ----
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const PshbOperatingPoint& pshbOpPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // ---- PSHB-specific calculations ----
    static double compute_effective_duty_cycle(double phaseShiftDeg);
    static double compute_output_voltage(double Vin, double Deff, double n,
                                         double Vd, BRectifierType rectType);
    static double compute_turns_ratio(double Vin, double Vo, double Deff,
                                      double Vd, BRectifierType rectType);
    static double compute_output_inductance(double Vo, double Deff, double Fs,
                                            double Io, double rippleRatio);
    static double compute_primary_rms_current(double Io, double n, double Deff);

    // ---- Extra components ----
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // ---- SPICE simulation ----
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


class AdvancedPshb : public Pshb {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredMagnetizingInductance;
    std::optional<double> desiredSeriesInductance;
    std::optional<double> desiredOutputInductance;

public:
    AdvancedPshb() = default;
    ~AdvancedPshb() = default;
    AdvancedPshb(const json& j);

    Inputs process();

    const double& get_desired_magnetizing_inductance() const { return desiredMagnetizingInductance; }
    void set_desired_magnetizing_inductance(const double& value) { desiredMagnetizingInductance = value; }
    const std::vector<double>& get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double>& value) { desiredTurnsRatios = value; }
    std::optional<double> get_desired_series_inductance() const { return desiredSeriesInductance; }
    void set_desired_series_inductance(std::optional<double> value) { desiredSeriesInductance = value; }
    std::optional<double> get_desired_output_inductance() const { return desiredOutputInductance; }
    void set_desired_output_inductance(std::optional<double> value) { desiredOutputInductance = value; }
};

void from_json(const json& j, AdvancedPshb& x);
void to_json(json& j, const AdvancedPshb& x);

inline void from_json(const json& j, AdvancedPshb& x) {
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_phase_shift(get_stack_optional<double>(j, "maximumPhaseShift"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<PshbOperatingPoint>>());
    x.set_output_inductance(get_stack_optional<double>(j, "outputInductance"));
    x.set_rectifier_type(get_stack_optional<BRectifierType>(j, "rectifierType"));
    x.set_series_inductance(get_stack_optional<double>(j, "seriesInductance"));
    x.set_use_leakage_inductance(get_stack_optional<bool>(j, "useLeakageInductance"));
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_magnetizing_inductance(j.at("desiredMagnetizingInductance").get<double>());
    x.set_desired_series_inductance(get_stack_optional<double>(j, "desiredSeriesInductance"));
    x.set_desired_output_inductance(get_stack_optional<double>(j, "desiredOutputInductance"));
}

inline void to_json(json& j, const AdvancedPshb& x) {
    j = json::object();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumPhaseShift"] = x.get_maximum_phase_shift();
    j["operatingPoints"] = x.get_operating_points();
    j["outputInductance"] = x.get_output_inductance();
    j["rectifierType"] = x.get_rectifier_type();
    j["seriesInductance"] = x.get_series_inductance();
    j["useLeakageInductance"] = x.get_use_leakage_inductance();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredMagnetizingInductance"] = x.get_desired_magnetizing_inductance();
    j["desiredSeriesInductance"] = x.get_desired_series_inductance();
    j["desiredOutputInductance"] = x.get_desired_output_inductance();
}

} // namespace OpenMagnetics
