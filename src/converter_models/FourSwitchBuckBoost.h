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
 * @brief Four-Switch Buck-Boost (4SBB / non-inverting buck-boost / H-bridge
 *        buck-boost / FSBB) DC-DC converter.
 *
 * Single-inductor synchronous topology with two half-bridges. Used in laptops,
 * USB-PD chargers, drones, single-cell-Li chargers, and automotive 12 V↔48 V
 * mild-hybrid rails — every application where Vin can be **above OR below**
 * Vout. Examples: TI LM5175/LM5176/LM5177/TPS55288/BQ25710/BQ25713;
 * ADI LT8390/LT8391/LT8392; Renesas ISL81801/ISL78224/ISL78226;
 * MPS MPQ4214/MP4247.
 *
 * The controller picks one of three operating regions per Vin value
 * (FOUR_SWITCH_BUCK_BOOST_PLAN.md §A.5):
 *
 *   BUCK       (Vo/Vin < 1 - hysteresis): Q3 always-ON, Q4 always-OFF,
 *              Q1/Q2 PWM. Inductor sees +Vin-Vo / -Vo (classical buck).
 *   BOOST      (Vo/Vin > 1 + hysteresis): Q1 always-ON, Q2 always-OFF,
 *              Q3/Q4 PWM. Inductor sees +Vin / -(Vo-Vin) (classical boost).
 *   BUCK_BOOST (transition band): all four switches modulate. Two
 *              sub-modes per `transitionMode`:
 *                * SIMULTANEOUS (Mode 1): Q1+Q3 ON during D, Q2+Q4 ON
 *                  during 1-D. Flyback-like; M = D/(1-D).
 *                * SPLIT_PWM (Mode 2, LM5176/LT8390 default): four
 *                  sub-intervals per Tsw; M = D_buck / (1-D_boost).
 *
 * The single inductor is sized for the **worst-case region** across the
 * Vin sweep:  L = max(L_buck @ Vin_max, L_boost @ Vin_min). Saturation
 * is dominated by the BOOST region @ Vin_min (DC bias = Iout·Vo/Vin_min);
 * core loss is typically dominated by BUCK_BOOST mode.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW
 * =====================================================================
 *
 *   Vin ──┬── Q1 ──┬─── L ───┬── Q3 ──┬── Vout
 *         │  (HS)  │         │  (HS)  │
 *         Cin     Q2       [SW2]     Q4         Co ── Rload
 *         │  (LS)  │         │  (LS)  │         │
 *         │        │         │        │         │
 *        GND ─── GND ─── [SW1] ── GND ── GND ─── GND
 *
 *   BUCK (Vin >> Vo):    Q3 ON, Q4 OFF, Q1/Q2 PWM(D_buck)
 *   BOOST (Vin << Vo):   Q1 ON, Q2 OFF, Q3/Q4 PWM(D_boost)
 *   BUCK_BOOST (Vin≈Vo): all four PWM
 *
 * =====================================================================
 * MODULATION CONVENTION
 * =====================================================================
 *
 *   BUCK       : D_buck = Vo / (Vin · η),       applied to Q1
 *   BOOST      : D_boost = 1 - Vin · η / Vo,    applied to Q4
 *   BUCK_BOOST : Mode 1 (simultaneous):  D = Vo / (Vin + Vo)
 *                Mode 2 (split-PWM):     D_buck + D_boost = 1 + δ
 *                                        with δ ≈ 0.05 (LM5176)
 *
 * =====================================================================
 * KEY EQUATIONS
 * =====================================================================
 *
 *   Buck region:
 *     M = D_buck                                                [TI SLVA535B Eq. 1]
 *     I_L,avg = Iout
 *     ΔI_L = Vo·(1−D_buck)/(L·Fs) = (Vin−Vo)·D_buck/(L·Fs)      [Eq. 3]
 *     L_min = Vo·(Vin_max−Vo)/(K·Fs·Vin_max·Iout)               [Eq. 5]
 *
 *   Boost region:
 *     M = 1/(1−D_boost)                                          [Eq. 2]
 *     I_L,avg = Iout/(1−D_boost) = Iout·Vo/Vin   ← largest DC bias
 *     ΔI_L = Vin·D_boost/(L·Fs)                                  [Eq. 4]
 *     L_min = Vin_min²·(Vo−Vin_min)/(K·Fs·Iout·Vo²)              [Eq. 8]
 *
 *   Buck-boost region (Mode 1, SIMULTANEOUS):
 *     M = D/(1−D),  D = Vo/(Vin+Vo)                              [Restrepo IEEE 8986423]
 *     ΔI_L = Vin·D/(L·Fs) (large — flyback-like)
 *
 *   Buck-boost region (Mode 2, SPLIT_PWM):
 *     M = D_buck / (1−D_boost)                                   [LM5176 d.s. §8.3.4]
 *
 *   DCM critical inductance:
 *     K_crit (buck)  = 1 − D_buck
 *     K_crit (boost) = D_boost · (1 − D_boost)²                  [Erickson §5]
 *     K_crit (BB-1)  = (1 − D)²
 *
 *   Worst-case design (the headline):
 *     Saturation peak : boost @ Vin_min  (I_pk = Iout·Vo/Vin_min + ΔI_L/2)
 *     Copper / RMS    : boost @ Vin_min
 *     Core loss       : buck-boost @ Vin ≈ Vo
 *     L sizing        : L = max(L_buck@VinMax, L_boost@VinMin)
 *
 * =====================================================================
 * REFERENCES
 * =====================================================================
 *
 * [1] TI SLVA535B — Basic Calculations of a 4-Switch Buck-Boost Power Stage.
 * [2] TI LM5175 / LM5176 / LM5177 datasheets — synchronous 4SBB controllers.
 * [3] ADI LT8390 / LT8391 / LT8392 datasheets — peakBuckPeakBoost family.
 * [4] Restrepo et al., IEEE 8986423 — operation modes, control, experiments.
 * [5] Yan et al., IEEE 9219154 — MPC smooth mode transition.
 * [6] Renesas ISL81801 / ISL78224 / ISL78226 — automotive bidirectional.
 * [7] Erickson & Maksimović, "Fundamentals of Power Electronics", 3rd ed.
 *     Ch. 5 (DCM), Ch. 6 (cascade connection).
 *
 * =====================================================================
 * ACCURACY DISCLAIMER
 * =====================================================================
 *
 * v1 limitations (per FOUR_SWITCH_BUCK_BOOST_PLAN.md §A.4):
 *   - Bootstrap-capacitor dynamics on Q1/Q3 high-side gates not modelled;
 *     SPICE uses ideal SW elements with PULSE drives.
 *   - Coss(Vds) variation, body-diode reverse recovery, gate-driver
 *     propagation delay, dead-time blanking not represented.
 *   - Mode-transition glitches at the BUCK ↔ BUCK_BOOST and
 *     BUCK_BOOST ↔ BOOST boundaries are NOT simulated; the analytical
 *     solver picks one region per OP based on the controller hysteresis.
 *     Real LM5176/LT8390 controllers blend smoothly over a ~5 % Vin band,
 *     introducing additional ripple that is under-predicted here.
 *   - Inductor saturation: peak current at the BUCK ↔ BUCK_BOOST mode
 *     boundary during a load step can exceed steady-state predictions
 *     by ~20 % (per ADI app note). v1 reports the steady-state worst
 *     case only.
 *   - Quasi-resonant ZVS extensions (LT8390 boost leg) omitted.
 *   - MPC controller mode (Yan IEEE 9219154) out of scope.
 *
 * =====================================================================
 * DISAMBIGUATION
 * =====================================================================
 *
 *   * 4SBB vs 2-switch buck-boost (inverting) — the 2-switch version
 *     inverts the output and is 2nd-order; 4SBB is non-inverting and uses
 *     two half-bridges with one inductor.
 *   * 4SBB vs SEPIC / Cuk / Zeta — all four are non-inverting buck-boost.
 *     SEPIC/Cuk/Zeta are 4th-order (two L + two C); 4SBB is 2nd-order
 *     (one L) but uses four active switches. SEPIC/Cuk/Zeta have RHPZ
 *     in CCM; 4SBB has RHPZ only in BOOST and BUCK_BOOST regions.
 *   * 4SBB vs DAB — DAB is isolated (transformer); 4SBB is non-isolated.
 *     Both have 4 switches, but DAB modulates phase shift, 4SBB modulates
 *     duty.
 */
class FourSwitchBuckBoost : public MAS::FourSwitchBuckBoost, public Topology {
public:
    enum class Region { BUCK, BOOST, BUCK_BOOST };

private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 30;

    // Maximum operating duty cycle. In all three regions (buck-only,
    // boost-only, buck-boost) the static duty fns previously had no
    // configurable ceiling and the SPLIT_PWM buck-boost path silently
    // clamped to [0.05, 0.95]. Now configurable; clamps are replaced
    // with throws. Mirrors Cuk (757d50cb), Sepic (9725ba41),
    // Zeta (31fa57c7), Weinberg (a61872f8).
    std::optional<double> maximumDutyCycle = 0.95;

    // Internal sizing knobs.
    double inductorRippleRatioDefault = 0.4;
    double outputRippleRatioDefault   = 0.01;
    double transitionHysteresisDefault = 0.15;

    // Switch / driver schema-less knobs (not exposed in MAS schema yet).
    double minOnTime  = 80e-9;
    double minOffTime = 80e-9;

protected:
    // Per-OP analytical diagnostics (Guide §2.4).
    mutable Region lastRegion                   = Region::BUCK;
    mutable double lastDutyCycleBuck            = 0.0;
    mutable double lastDutyCycleBoost           = 0.0;
    mutable double lastConversionRatio          = 0.0;
    mutable bool   lastIsCcm                    = true;
    mutable double lastDcmK                     = 0.0;
    mutable double lastDcmKcrit                 = 0.0;
    mutable double lastInductorAverageCurrent   = 0.0;
    mutable double lastInductorPeakToPeak       = 0.0;
    mutable double lastInductorPeakCurrent      = 0.0;
    mutable double lastInductorRmsCurrent       = 0.0;
    mutable double lastInputCurrentAverage      = 0.0;
    mutable double lastInputCurrentRipple       = 0.0;
    mutable double lastQ1RmsCurrent             = 0.0;
    mutable double lastQ2RmsCurrent             = 0.0;
    mutable double lastQ3RmsCurrent             = 0.0;
    mutable double lastQ4RmsCurrent             = 0.0;
    mutable double lastSizedInductance          = 0.0;
    mutable double lastSizedInputCapacitance    = 0.0;
    mutable double lastSizedOutputCapacitance   = 0.0;
    mutable double lastOutputVoltageRipple      = 0.0;

    // Extra-component waveforms (cleared in process_operating_points).
    mutable std::vector<Waveform> extraInductorVoltageWaveforms;
    mutable std::vector<Waveform> extraInductorCurrentWaveforms;
    mutable std::vector<Waveform> extraInputCapVoltageWaveforms;
    mutable std::vector<Waveform> extraInputCapCurrentWaveforms;
    mutable std::vector<Waveform> extraOutputCapVoltageWaveforms;
    mutable std::vector<Waveform> extraOutputCapCurrentWaveforms;

public:
    bool _assertErrors = false;

    FourSwitchBuckBoost(const json& j);
    FourSwitchBuckBoost() {}

    MAS::Topologies topology_kind() const override {
        return MAS::Topologies::FOUR_SWITCH_BUCK_BOOST_CONVERTER;
    }
    bool is_bridge_topology() const override { return true; }

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int v) { numPeriodsToExtract = v; }
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int v) { numSteadyStatePeriods = v; }

    std::optional<double> get_maximum_duty_cycle() const { return maximumDutyCycle; }
    void set_maximum_duty_cycle(std::optional<double> value) { this->maximumDutyCycle = value; }

    // Diagnostic accessors.
    Region get_last_region()                   const { return lastRegion; }
    double get_last_duty_cycle_buck()          const { return lastDutyCycleBuck; }
    double get_last_duty_cycle_boost()         const { return lastDutyCycleBoost; }
    double get_last_conversion_ratio()         const { return lastConversionRatio; }
    bool   get_last_is_ccm()                   const { return lastIsCcm; }
    double get_last_dcm_k()                    const { return lastDcmK; }
    double get_last_dcm_kcrit()                const { return lastDcmKcrit; }
    double get_last_inductor_average_current() const { return lastInductorAverageCurrent; }
    double get_last_inductor_peak_to_peak()    const { return lastInductorPeakToPeak; }
    double get_last_inductor_peak_current()    const { return lastInductorPeakCurrent; }
    double get_last_inductor_rms_current()     const { return lastInductorRmsCurrent; }
    double get_last_input_current_average()    const { return lastInputCurrentAverage; }
    double get_last_input_current_ripple()     const { return lastInputCurrentRipple; }
    double get_last_q1_rms_current()           const { return lastQ1RmsCurrent; }
    double get_last_q2_rms_current()           const { return lastQ2RmsCurrent; }
    double get_last_q3_rms_current()           const { return lastQ3RmsCurrent; }
    double get_last_q4_rms_current()           const { return lastQ4RmsCurrent; }
    double get_last_sized_inductance()         const { return lastSizedInductance; }
    double get_last_sized_input_capacitance()  const { return lastSizedInputCapacitance; }
    double get_last_sized_output_capacitance() const { return lastSizedOutputCapacitance; }
    double get_last_output_voltage_ripple()    const { return lastOutputVoltageRipple; }

    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios,
                                                         double magnetizingInductance) override;
    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const TopologyExcitation& excitation,
        double inductance);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode = ExtraComponentsMode::IDEAL,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // Static analytical helpers.
    static Region classify_region(double Vin, double Vo, double hysteresisRatio,
                                  double minOnTime, double minOffTime, double Fs);
    static double compute_buck_duty(double Vin, double Vo, double efficiency, double maximumDutyCycle = 0.95);
    static double compute_boost_duty(double Vin, double Vo, double efficiency, double maximumDutyCycle = 0.95);
    static void   compute_buck_boost_duties(double Vin, double Vo, double efficiency,
                                            TransitionMode mode,
                                            double& D_buck_out, double& D_boost_out,
                                            double maximumDutyCycle = 0.95);
    static double compute_inductor_buck_region(double Vin, double Vo, double Iout,
                                               double Fs, double rippleRatio);
    static double compute_inductor_boost_region(double Vin, double Vo, double Iout,
                                                double Fs, double rippleRatio);
    static double compute_worst_case_inductance(double VinMin, double VinMax, double Vo,
                                                double Iout, double Fs, double rippleRatio);
    static double compute_inductor_dc_bias_boost(double Vin, double Vo, double Iout);

    // SPICE.
    std::string generate_ngspice_circuit(double inductance,
                                         size_t inputVoltageIndex = 0,
                                         size_t operatingPointIndex = 0);
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance);
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        double inductance, size_t numberOfPeriods = 2);
};

class AdvancedFourSwitchBuckBoost : public FourSwitchBuckBoost {
private:
    double desiredInductance = 0.0;

public:
    AdvancedFourSwitchBuckBoost() = default;
    ~AdvancedFourSwitchBuckBoost() = default;

    AdvancedFourSwitchBuckBoost(const json& j);

    Inputs process();

    // Issue M1: build DesignRequirements directly from desired* fields,
    // bypassing the parent's PDR which auto-sizes inductance from ripple.
    DesignRequirements process_design_requirements() override;

    const double& get_desired_inductance() const { return desiredInductance; }
    double& get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double& value) { desiredInductance = value; }
};

void from_json(const json& j, AdvancedFourSwitchBuckBoost& x);
void to_json(json& j, const AdvancedFourSwitchBuckBoost& x);

inline void from_json(const json& j, AdvancedFourSwitchBuckBoost& x) {
    x.set_bidirectional(get_stack_optional<bool>(j, "bidirectional"));
    x.set_control_mode(get_stack_optional<MAS::ControlMode>(j, "controlMode"));
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<TopologyExcitation>>());
    x.set_output_voltage_ripple_ratio(get_stack_optional<double>(j, "outputVoltageRippleRatio"));
    x.set_phase_count(get_stack_optional<int64_t>(j, "phaseCount"));
    x.set_transition_hysteresis_ratio(get_stack_optional<double>(j, "transitionHysteresisRatio"));
    x.set_transition_mode(get_stack_optional<MAS::TransitionMode>(j, "transitionMode"));
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json& j, const AdvancedFourSwitchBuckBoost& x) {
    j = json::object();
    j["bidirectional"] = x.get_bidirectional();
    j["controlMode"] = x.get_control_mode();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumSwitchCurrent"] = x.get_maximum_switch_current();
    j["operatingPoints"] = x.get_operating_points();
    j["outputVoltageRippleRatio"] = x.get_output_voltage_ripple_ratio();
    j["phaseCount"] = x.get_phase_count();
    j["transitionHysteresisRatio"] = x.get_transition_hysteresis_ratio();
    j["transitionMode"] = x.get_transition_mode();
    j["desiredInductance"] = x.get_desired_inductance();
}

} // namespace OpenMagnetics
