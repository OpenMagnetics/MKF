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
 * @brief Asymmetric Half-Bridge (AHB) DC-DC Converter
 *
 * Inherits from MAS::AsymmetricHalfBridge and the Topology interface.
 *
 * Reference: Imbertson & Mohan, "Asymmetrical Duty Cycle Permits Zero
 * Switching Loss in PWM Circuits with No Conduction Loss Penalty",
 * IEEE Trans. Industry Applications 29(1):121-125, 1993.
 * https://doi.org/10.1109/28.195897
 *
 * =====================================================================
 * TOPOLOGY DISAMBIGUATION (READ THIS FIRST)
 * =====================================================================
 *
 * The Asymmetric Half-Bridge (AHB) is a single-leg, two-switch isolated
 * PWM converter. The two MOSFETs (Q1 high-side, Q2 low-side) operate
 * with COMPLEMENTARY duty cycles D / (1-D) (with a small dead-time
 * between them). A DC-blocking capacitor Cb sits in series with the
 * primary winding; in steady state it charges to V_Cb = (1-D)*Vin and
 * makes the primary voltage SWING ASYMMETRICALLY between
 * +(1-D)*Vin (during D*Tsw) and -D*Vin (during (1-D)*Tsw).
 *
 * AHB is NOT:
 *   - The Pinheiro-Barbi 3-level "Phase-Shifted Half Bridge" (PSHB):
 *     PSHB uses a 3-level NPC leg (4 stacked switches + 2 clamp
 *     diodes) and delivers a 3-level primary voltage; AHB uses 2
 *     switches and delivers an asymmetric 2-level primary voltage.
 *     See PhaseShiftedHalfBridge.h for the PSHB documentation.
 *   - The split-capacitor Half-Bridge (HB): split-cap HB uses 2
 *     switches and 2 input caps in series, primary swings symmetrically
 *     between +Vin/2 and -Vin/2 (no DC-blocking cap, symmetric duty).
 *   - The Half-Bridge LLC resonant converter (HB-LLC): HB-LLC uses 2
 *     switches at 50% duty + a resonant tank Lr/Lm/Cr; controlled by
 *     frequency, not duty.
 *   - The Phase-Shifted Full Bridge (PSFB): PSFB uses 4 switches
 *     (2 legs) with phase shift between legs; primary swings
 *     symmetrically between +Vin and -Vin.
 *
 * =====================================================================
 * WHAT THIS MODEL IS  (V1 baseline, classic CT secondary)
 * =====================================================================
 *
 *   - Single leg: Q1 (HS) + Q2 (LS) with antiparallel body diodes
 *     between Vin_p and GND, midpoint = mid_A.
 *   - DC-blocking capacitor Cb in series between mid_A and the primary
 *     winding. In steady state V_Cb = (1-D)*Vin (Imbertson-Mohan eq. 3).
 *   - Primary winding (with leakage Llk and magnetizing Lm) carries
 *     asymmetric 2-level voltage:
 *         v_pri = +(1-D)*Vin  during Q1 ON (D*Tsw)
 *         v_pri = -D*Vin      during Q2 ON ((1-D)*Tsw)
 *     Volt-second balance is automatic: V_Cb tracks D continuously.
 *   - Secondary: rectifier-type aware (centerTapped / fullBridge /
 *     currentDoubler / ahbFlyback) + Lo + Co + Rload.
 *   - Conversion ratio:  Vo = 2*D*(1-D)*Vin/n   (CT/FB)
 *                        Vo =   D*(1-D)*Vin/n   (CD)
 *     The gain curve is NON-MONOTONIC; it peaks at D=0.5 with
 *     Vo,max = Vin/(2*n). Conventional control uses D in [0, 0.5].
 *
 * =====================================================================
 * IMPLEMENTATION STATUS  (P1-P12 complete — May 2026)
 * =====================================================================
 *
 * Full implementation of V1-V6 variants:
 *   V1 = CT secondary  (default)         V4 = synchronous rectifier (set_use_synchronous_rectifier(true))
 *   V2 = full-bridge secondary           V5 = AHB-Flyback (rectifierType = ahbFlyback)
 *   V3 = current-doubler secondary       V6 = multi-output (multiple outputVoltages)
 *
 * V5 (AHB-Flyback) is implemented as a two-switch / active-clamp
 * flyback (no series Cb; Lm carries the energy storage; mean(iLm) > 0);
 * the gain is Vo = D·Vin/((1-D)·n), singular at D=1.
 *
 * V6 multi-output: cross-regulation is APPROXIMATE — per-secondary
 * leakage inductance projection (DAB-style) is used but the result is
 * accurate to ~20-100% per secondary unless leakage is supplied.
 * A warn-once message is printed by process_operating_points.
 *
 * Comparison to PSHB (single-leg 3-level NPC implementation in MKF):
 *   - PSHB delivers a 3-level primary voltage (+Vin/2, 0, -Vin/2).
 *   - AHB delivers an asymmetric 2-level primary voltage
 *     (+(1-D)*Vin, -D*Vin).
 *   - PSHB uses 4 switches; AHB uses 2 switches.
 *   - PSHB conversion ratio is monotonic in Deff; AHB is non-monotonic
 *     in D and peaks at D=0.5.
 *
 * Comparison to ActiveClampForward:
 *   - Both share the DC-blocking-cap + asymmetric magnetizing ramp
 *     pattern.
 *   - ACF uses one main switch + one auxiliary clamp switch; AHB uses
 *     two equal-stress switches (HS + LS) on the same leg.
 *
 * =====================================================================
 * EQUATIONS CHEAT-SHEET (CCM, ideal)
 * =====================================================================
 *
 *   DC-blocking cap voltage : V_Cb = (1−D)·Vin                        [Imbertson eq. 3]
 *   Primary voltage         : v_pri = +(1−D)·Vin during D·Tsw,
 *                             v_pri = −D·Vin    during (1−D)·Tsw      ← asymmetric, automatic flux balance
 *   Conversion ratio (CT/FB): Vo = 2·D·(1−D)·Vin / n                  ← non-monotonic, peak at D=0.5
 *   Conversion ratio (CD)   : Vo = D·(1−D)·Vin / n
 *   Conversion ratio (V5)   : Vo = D·Vin / ((1−D)·n)                  (active-clamp flyback)
 *   Maximum gain at D=0.5   : Vo,max = Vin / (2·n)
 *   Switch peak voltage     : V_Q,pk = Vin (each switch in V1-V4);
 *                             Q2 sees Vin + Vo·n in V5 (clamp adder).
 *   Switch RMS              : I_Q1,rms = I_pri,Q1·√D ;
 *                             I_Q2,rms = I_pri,Q2·√(1−D)              ← unequal when D ≠ 0.5
 *   ZVS condition           : ½·(Llk+Lm,refl)·I_pri² ≥ 2·Coss·Vin²
 *   Dead-time               : t_d ≥ π·√(Llk·2·Coss)
 *   Steady-state ΔΦ         : ΔΦ_ss = D·(1−D)·Vin·Tsw / Np
 *   Bpk                     : D·(1−D)·Vin·Tsw / (2·Np·Ae)
 *   Magnetising DC bias (V1): mean(i_Lm) = (Io_per_inductor / n)·(1−2D)
 *
 *   Sizing — Cb (ripple ≤ 5 % V_Cb):  Cb ≥ I_pri,pk·D / (fsw·ΔV_Cb)   [ON AN-4153 eq. 16]
 *   Sizing — Lo (CCM)              :  Lo ≥ Vo·(1 − 2·D·(1−D)) / (ΔILo·fsw)
 *   Sizing — Lm (ZVS assist)       :  Lm ≤ ((1−D)·Vin·D·Tsw) / (2·I_m,target)
 *   Sizing — Lm (V5 30 % ripple)   :  Lm = Vin·D·n·(1−D)·Tsw / (0.30·Io)
 *   Sizing — Co (V5 flyback)       :  Co = Io·D / (fsw·ΔVo)
 *
 * =====================================================================
 * REFERENCES
 * =====================================================================
 *
 *   [1] P. Imbertson & N. Mohan, "Asymmetrical Duty Cycle Permits Zero
 *       Switching Loss in PWM Circuits with No Conduction Loss Penalty,"
 *       IEEE TIA, 29(1):121-125, 1993. https://doi.org/10.1109/28.195897
 *   [2] R. L. Steigerwald, "Designing High-Efficiency Asymmetric Half-Bridge
 *       PWM Converters," TI SLUP223 (Power Supply Design Seminar).
 *   [3] ON Semi AN-4153, "Designing Asymmetric PWM Half-Bridge DC/DC
 *       Converter with FPS."
 *   [4] STMicro AN2852, "EVL6591-90WADP - 90 W AC-DC AHB adapter."
 *   [5] Pressman, Billings & Morey, "Switching Power Supply Design" 3e Ch. 16.
 *   [6] Erickson & Maksimović, "Fundamentals of Power Electronics" 2e §6.3
 *       (DC analysis), §11.1.4 (small-signal).
 *   [7] TI SLUA121, "The Current-Doubler Rectifier: An Alternative
 *       Rectification Technique."
 *
 * =====================================================================
 * ACCURACY DISCLAIMER
 * =====================================================================
 *
 *   The analytical model assumes ideal switches/diodes, lossless reactives,
 *   ideal coupling K=1 (transformed to K=0.999 for SPICE numerical
 *   stability), and PWL transitions. The SPICE model omits MOSFET Coss
 *   variation with V_DS, body-diode reverse recovery, gate-driver
 *   propagation delay, transformer inter-winding capacitance, and
 *   core-loss resistance. Predicted η is an upper bound; expect 2-4 %
 *   loss vs. measured at full load.
 *
 *   V5 (AHB-Flyback) is a two-switch / active-clamp flyback approximation:
 *   the analytical model assumes the secondary perfectly clamps v_pri to
 *   −Vo·n during the OFF interval; in practice clamp-cap ringing and
 *   leakage-energy recycling will deviate. V5 gain is singular at D=1.
 *
 *   V6 (multi-output) cross-regulation is approximate: per-output current
 *   waveforms are projected via load-share (Vo_k·Io_k / Σ Vo_j·Io_j) of
 *   the primary current; without per-secondary leakage inductance the
 *   per-output magnetic loss may be off by 20-100 %. A warn-once message
 *   is printed by process_operating_points in this case.
 */
class AsymmetricHalfBridge : public MAS::AsymmetricHalfBridge, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 20;

    // Computed design values (populated by process_design_requirements once
    // P6 lands; default-zero in the stub).
    double computedTurnsRatio            = 0;
    double computedOutputInductance      = 0;   // Lo
    double computedMagnetizingInductance = 0;   // Lm
    double computedLeakageInductance     = 1e-6;// Llk
    double computedDcBlockingCapacitance = 0;   // Cb
    double computedOutputCapacitance     = 0;   // Co
    double computedDutyCycle             = 0;
    double computedDeadTime              = 100e-9;
    double computedDiodeVoltageDrop      = 0.6;

    // Schema-less device parameter (per-FET output capacitance).
    double mosfetCoss = 200e-12;

    // Schema-less V4 toggle. When true, the secondary diodes are replaced
    // by SR FETs in the SPICE netlist (driven by the complementary primary
    // PWM gate signals — see generate_ngspice_circuit) and the analytical
    // diode-voltage drop is reduced to ~Rds_on·I (modelled as a 50 mV
    // floor); see set_use_synchronous_rectifier docstring for caveats.
    bool useSynchronousRectifier = false;

    // Per-OP diagnostics (populated once P3 lands; default-zero in the stub).
    mutable double lastDutyCycle                       = 0;
    mutable double lastConversionRatio                 = 0;
    mutable double lastDcBlockingCapVoltage            = 0;  // V_Cb = (1-D)*Vin
    mutable double lastDcBlockingCapRipple             = 0;
    mutable double lastPrimaryPeakVoltagePositive      = 0;
    mutable double lastPrimaryPeakVoltageNegative      = 0;
    mutable double lastSwitchPeakVoltageQ1             = 0;
    mutable double lastSwitchPeakVoltageQ2             = 0;
    mutable double lastSwitchRmsCurrentQ1              = 0;
    mutable double lastSwitchRmsCurrentQ2              = 0;
    mutable double lastZvsMargin                       = 0;
    mutable double lastResonantTransitionTime          = 0;
    mutable double lastSteadyStateFluxExcursion        = 0;
    mutable double lastTransientFluxExcursionEstimate  = 0;
    mutable double lastMagnetizingCurrentRipple        = 0;
    mutable double lastOutputInductorRipple            = 0;
    mutable int    lastOperatingMode                   = 0; // 0=CCM, 1=DCM
    mutable int    lastRectifierType                   = 0; // 0=CT,1=CD,2=FB,3=Flyback

    // ---- Per-OP diagnostic vectors (one entry per V_in solved in
    // process_operating_points). Cleared at top, pushed at end of each iter.
    mutable std::vector<std::string> perOpName;
    mutable std::vector<double>      perOpDutyCycle;
    mutable std::vector<double>      perOpConversionRatio;
    mutable std::vector<double>      perOpDcBlockingCapVoltage;
    mutable std::vector<double>      perOpDcBlockingCapRipple;
    mutable std::vector<double>      perOpPrimaryPeakVoltagePositive;
    mutable std::vector<double>      perOpPrimaryPeakVoltageNegative;
    mutable std::vector<double>      perOpSwitchPeakVoltageQ1;
    mutable std::vector<double>      perOpSwitchPeakVoltageQ2;
    mutable std::vector<double>      perOpSwitchRmsCurrentQ1;
    mutable std::vector<double>      perOpSwitchRmsCurrentQ2;
    mutable std::vector<double>      perOpZvsMargin;
    mutable std::vector<double>      perOpResonantTransitionTime;
    mutable std::vector<double>      perOpSteadyStateFluxExcursion;
    mutable std::vector<double>      perOpTransientFluxExcursionEstimate;
    mutable std::vector<double>      perOpMagnetizingCurrentRipple;
    mutable std::vector<double>      perOpOutputInductorRipple;
    mutable std::vector<int>         perOpOperatingMode;
    mutable std::vector<int>         perOpRectifierType;

    // Extra-component waveforms (cleared in process_operating_points once P4 lands).
    mutable std::vector<Waveform> extraLoVoltageWaveforms;
    mutable std::vector<Waveform> extraLoCurrentWaveforms;
    mutable std::vector<Waveform> extraLo2VoltageWaveforms;   // CD only
    mutable std::vector<Waveform> extraLo2CurrentWaveforms;   // CD only
    mutable std::vector<Waveform> extraCbVoltageWaveforms;
    mutable std::vector<Waveform> extraCbCurrentWaveforms;


public:
    bool _assertErrors = false;

    AsymmetricHalfBridge(const json& j);
    AsymmetricHalfBridge() = default;

    MAS::Topology topology_kind() const override { return MAS::Topology::ASYMMETRIC_HALF_BRIDGE_CONVERTER; }
    bool is_bridge_topology() const override { return true; }

    // ---- Boilerplate accessors (Guide §2) ----
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { numPeriodsToExtract = value; }
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { numSteadyStatePeriods = value; }

    double get_computed_turns_ratio() const { return computedTurnsRatio; }
    double get_computed_output_inductance() const { return computedOutputInductance; }
    double get_computed_magnetizing_inductance() const { return computedMagnetizingInductance; }
    double get_computed_leakage_inductance() const { return computedLeakageInductance; }
    double get_computed_dc_blocking_capacitance() const { return computedDcBlockingCapacitance; }
    double get_computed_output_capacitance() const { return computedOutputCapacitance; }
    double get_computed_duty_cycle() const { return computedDutyCycle; }
    double get_computed_dead_time() const { return computedDeadTime; }
    double get_computed_diode_voltage_drop() const { return computedDiodeVoltageDrop; }

    /// Pin specific component values. When unset (==0), process_operating_*
    /// will size the corresponding element analytically (Lo: 30 % ripple
    /// rule). Used by AdvancedAHB and by tests that need to drive a known
    /// DCM design.
    void set_computed_output_inductance(double value) { computedOutputInductance = value; }
    void set_computed_leakage_inductance(double value) { computedLeakageInductance = value; }
    void set_computed_magnetizing_inductance(double value) { computedMagnetizingInductance = value; }
    void set_computed_dc_blocking_capacitance(double value) { computedDcBlockingCapacitance = value; }
    void set_computed_output_capacitance(double value) { computedOutputCapacitance = value; }

    /// Per-MOSFET output capacitance (F) used for ZVS energy and resonant
    /// transition-time predictions. Default 200 pF. Schema-less; set via setter.
    double get_mosfet_coss() const { return mosfetCoss; }
    void set_mosfet_coss(double value) { mosfetCoss = value; }

    /// V4 synchronous-rectifier toggle (schema-less). When true:
    ///   - generate_ngspice_circuit replaces every secondary diode with
    ///     an S_* voltage-controlled switch driven by the complementary
    ///     primary PWM gate (CT: D_ra→S_ra on g2, D_rb→S_rb on g1; FB:
    ///     four SR FETs paired with the four bridge diodes' phases; CD:
    ///     two SR FETs on g1/g2).
    ///   - The analytical diode-voltage drop used by process_design_*
    ///     and the per-OP analytical waveforms is reduced to 0.05 V
    ///     (Rds_on·I floor) instead of the default 0.6 V Si diode drop,
    ///     reflecting the lower SR conduction loss. This affects
    ///     compute_turns_ratio, the sized n, and the diagnostics that
    ///     embed Vd (lastConversionRatio, lastDiodeVoltageDrop).
    /// Caveat: SR phasing in V4 follows the primary PWM signals; explicit
    /// SR-driver delay / dead-time is not modeled in the analytical path.
    bool get_use_synchronous_rectifier() const { return useSynchronousRectifier; }
    void set_use_synchronous_rectifier(bool value) {
        useSynchronousRectifier = value;
        // Lower analytical Vd floor when SR active. The 50 mV figure is
        // representative of a 30 mΩ SR FET at 1.5 A average secondary
        // current; surface via lastDiodeVoltageDrop so users can adjust.
        computedDiodeVoltageDrop = value ? 0.05 : 0.6;
    }

    // ---- Per-OP diagnostics (Guide §2.4) ----
    double get_last_duty_cycle() const { return lastDutyCycle; }
    double get_last_conversion_ratio() const { return lastConversionRatio; }
    /// Steady-state DC-blocking cap voltage V_Cb = (1-D)*Vin.
    double get_last_dc_blocking_cap_voltage() const { return lastDcBlockingCapVoltage; }
    double get_last_dc_blocking_cap_ripple() const { return lastDcBlockingCapRipple; }
    /// Positive primary voltage during Q1 ON: +(1-D)*Vin.
    double get_last_primary_peak_voltage_positive() const { return lastPrimaryPeakVoltagePositive; }
    /// Negative primary voltage during Q2 ON: -D*Vin.
    double get_last_primary_peak_voltage_negative() const { return lastPrimaryPeakVoltageNegative; }
    double get_last_switch_peak_voltage_q1() const { return lastSwitchPeakVoltageQ1; }
    double get_last_switch_peak_voltage_q2() const { return lastSwitchPeakVoltageQ2; }
    /// Asymmetric: Q1 and Q2 carry different RMS currents when D != 0.5.
    double get_last_switch_rms_current_q1() const { return lastSwitchRmsCurrentQ1; }
    double get_last_switch_rms_current_q2() const { return lastSwitchRmsCurrentQ2; }
    /// ZVS energy margin: 0.5*(Llk+Lm,refl)*I_pri^2 - 2*Coss*Vin^2.
    double get_last_zvs_margin() const { return lastZvsMargin; }
    /// Quarter-period of the ZVS resonant tank: pi*sqrt(Llk * 2*Coss).
    double get_last_resonant_transition_time() const { return lastResonantTransitionTime; }
    /// Steady-state peak flux density-related quantity (delta-Phi).
    double get_last_steady_state_flux_excursion() const { return lastSteadyStateFluxExcursion; }
    /// Worst-case transient flux excursion under a Vin step (only computed
    /// if inputVoltageStepRange is supplied in the schema).
    double get_last_transient_flux_excursion_estimate() const { return lastTransientFluxExcursionEstimate; }
    double get_last_magnetizing_current_ripple() const { return lastMagnetizingCurrentRipple; }
    double get_last_output_inductor_ripple() const { return lastOutputInductorRipple; }
    int get_last_operating_mode() const { return lastOperatingMode; }
    int get_last_rectifier_type() const { return lastRectifierType; }

    const std::vector<std::string>& get_per_op_name()                              const { return perOpName; }
    const std::vector<double>&      get_per_op_duty_cycle()                        const { return perOpDutyCycle; }
    const std::vector<double>&      get_per_op_conversion_ratio()                  const { return perOpConversionRatio; }
    const std::vector<double>&      get_per_op_dc_blocking_cap_voltage()           const { return perOpDcBlockingCapVoltage; }
    const std::vector<double>&      get_per_op_dc_blocking_cap_ripple()            const { return perOpDcBlockingCapRipple; }
    const std::vector<double>&      get_per_op_primary_peak_voltage_positive()     const { return perOpPrimaryPeakVoltagePositive; }
    const std::vector<double>&      get_per_op_primary_peak_voltage_negative()     const { return perOpPrimaryPeakVoltageNegative; }
    const std::vector<double>&      get_per_op_switch_peak_voltage_q1()            const { return perOpSwitchPeakVoltageQ1; }
    const std::vector<double>&      get_per_op_switch_peak_voltage_q2()            const { return perOpSwitchPeakVoltageQ2; }
    const std::vector<double>&      get_per_op_switch_rms_current_q1()             const { return perOpSwitchRmsCurrentQ1; }
    const std::vector<double>&      get_per_op_switch_rms_current_q2()             const { return perOpSwitchRmsCurrentQ2; }
    const std::vector<double>&      get_per_op_zvs_margin()                        const { return perOpZvsMargin; }
    const std::vector<double>&      get_per_op_resonant_transition_time()          const { return perOpResonantTransitionTime; }
    const std::vector<double>&      get_per_op_steady_state_flux_excursion()       const { return perOpSteadyStateFluxExcursion; }
    const std::vector<double>&      get_per_op_transient_flux_excursion_estimate() const { return perOpTransientFluxExcursionEstimate; }
    const std::vector<double>&      get_per_op_magnetizing_current_ripple()        const { return perOpMagnetizingCurrentRipple; }
    const std::vector<double>&      get_per_op_output_inductor_ripple()            const { return perOpOutputInductorRipple; }
    const std::vector<int>&         get_per_op_operating_mode()                    const { return perOpOperatingMode; }
    const std::vector<int>&         get_per_op_rectifier_type()                    const { return perOpRectifierType; }

    // ---- Topology interface (P1 implements run_checks; P3-P7 implement the rest) ----
    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;

    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const MAS::AhbOperatingPoint& opPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // ---- Static analytical helpers (P2 deliverable; stubs throw) ----
    static double compute_duty_cycle(double Vin, double Vo, double n,
                                     double diodeDrop, double efficiency,
                                     AhbRectifierType rectType);
    static double compute_dc_blocking_cap_voltage(double Vin, double D);
    static double compute_conversion_ratio(double D, double n,
                                           AhbRectifierType rectType);
    static double compute_turns_ratio(double Vin_min, double Vo, double D_max,
                                      double diodeDrop,
                                      AhbRectifierType rectType);
    static double compute_lo_min(double Vo, double D, double dILo, double fsw);
    static double compute_lm_min_for_zvs(double Vin, double D, double Tsw,
                                         double Im_target);
    static double compute_cb_min(double Iprimary_pk, double D, double dVCb,
                                 double fsw);
    static double compute_zvs_energy_balance(double Llk, double Lm_refl,
                                             double Ipri, double Coss, double Vin);
    static double compute_dead_time(double Llk, double Coss);
    static double compute_steady_state_flux_excursion(double Vin, double D,
                                                      double Tsw, double Np, double Ae);
    static double compute_transient_flux_excursion(double Vin, double dVin, double D,
                                                   double Tsw, double Np, double Ae,
                                                   double Lm, double Rdcr);

    // ---- Extra-components factory (P9 deliverable; stub throws) ----
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // ---- SPICE (P7 deliverable; stubs throw) ----
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);

    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);
};


/**
 * @brief Advanced AHB sizing override: lets the user pin specific
 * turns-ratios, magnetizing inductance, leakage inductance, output
 * inductor, DC-blocking capacitance, and output capacitance instead of
 * letting process_design_requirements compute them from the operating
 * points. P12 deliverable per ASYMMETRIC_HALF_BRIDGE_PLAN.md.
 */
class AdvancedAsymmetricHalfBridge : public AsymmetricHalfBridge {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredMagnetizingInductance = 0;
    std::optional<double> desiredLeakageInductance;
    std::optional<double> desiredOutputInductance;
    std::optional<double> desiredDcBlockingCapacitance;
    std::optional<double> desiredOutputCapacitance;

public:
    AdvancedAsymmetricHalfBridge() = default;
    ~AdvancedAsymmetricHalfBridge() = default;
    AdvancedAsymmetricHalfBridge(const json& j);

    // Override parent's process_design_requirements() so the generic
    // construct → PDR flow returns DR populated from desired* fields
    // (where supplied) rather than only auto-sized values. Issue M1.
    DesignRequirements process_design_requirements() override;

    Inputs process();

    const std::vector<double>& get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double>& value) { desiredTurnsRatios = value; }
    double get_desired_magnetizing_inductance() const { return desiredMagnetizingInductance; }
    void set_desired_magnetizing_inductance(double value) { desiredMagnetizingInductance = value; }
    std::optional<double> get_desired_leakage_inductance() const { return desiredLeakageInductance; }
    void set_desired_leakage_inductance(std::optional<double> value) { desiredLeakageInductance = value; }
    std::optional<double> get_desired_output_inductance() const { return desiredOutputInductance; }
    void set_desired_output_inductance(std::optional<double> value) { desiredOutputInductance = value; }
    std::optional<double> get_desired_dc_blocking_capacitance() const { return desiredDcBlockingCapacitance; }
    void set_desired_dc_blocking_capacitance(std::optional<double> value) { desiredDcBlockingCapacitance = value; }
    std::optional<double> get_desired_output_capacitance() const { return desiredOutputCapacitance; }
    void set_desired_output_capacitance(std::optional<double> value) { desiredOutputCapacitance = value; }
};

void from_json(const json& j, AdvancedAsymmetricHalfBridge& x);
void to_json(json& j, const AdvancedAsymmetricHalfBridge& x);

inline void from_json(const json& j, AdvancedAsymmetricHalfBridge& x) {
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_dc_blocking_capacitance(get_stack_optional<double>(j, "dcBlockingCapacitance"));
    x.set_magnetizing_inductance(get_stack_optional<double>(j, "magnetizingInductance"));
    x.set_leakage_inductance(get_stack_optional<double>(j, "leakageInductance"));
    x.set_use_leakage_inductance(get_stack_optional<bool>(j, "useLeakageInductance"));
    x.set_output_inductance(get_stack_optional<double>(j, "outputInductance"));
    x.set_rectifier_type(get_stack_optional<AhbRectifierType>(j, "rectifierType"));
    x.set_maximum_duty_cycle(get_stack_optional<double>(j, "maximumDutyCycle"));
    x.set_input_voltage_step_range(get_stack_optional<double>(j, "inputVoltageStepRange"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<MAS::AhbOperatingPoint>>());

    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_magnetizing_inductance(j.at("desiredMagnetizingInductance").get<double>());
    x.set_desired_leakage_inductance(get_stack_optional<double>(j, "desiredLeakageInductance"));
    x.set_desired_output_inductance(get_stack_optional<double>(j, "desiredOutputInductance"));
    x.set_desired_dc_blocking_capacitance(get_stack_optional<double>(j, "desiredDcBlockingCapacitance"));
    x.set_desired_output_capacitance(get_stack_optional<double>(j, "desiredOutputCapacitance"));
}

inline void to_json(json& j, const AdvancedAsymmetricHalfBridge& x) {
    j = json::object();
    j["inputVoltage"] = x.get_input_voltage();
    j["efficiency"] = x.get_efficiency();
    j["dcBlockingCapacitance"] = x.get_dc_blocking_capacitance();
    j["magnetizingInductance"] = x.get_magnetizing_inductance();
    j["leakageInductance"] = x.get_leakage_inductance();
    j["useLeakageInductance"] = x.get_use_leakage_inductance();
    j["outputInductance"] = x.get_output_inductance();
    j["rectifierType"] = x.get_rectifier_type();
    j["maximumDutyCycle"] = x.get_maximum_duty_cycle();
    j["inputVoltageStepRange"] = x.get_input_voltage_step_range();
    j["operatingPoints"] = x.get_operating_points();

    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredMagnetizingInductance"] = x.get_desired_magnetizing_inductance();
    j["desiredLeakageInductance"] = x.get_desired_leakage_inductance();
    j["desiredOutputInductance"] = x.get_desired_output_inductance();
    j["desiredDcBlockingCapacitance"] = x.get_desired_dc_blocking_capacitance();
    j["desiredOutputCapacitance"] = x.get_desired_output_capacitance();
}

} // namespace OpenMagnetics
