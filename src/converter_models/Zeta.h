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
 * @brief Zeta — non-inverting fourth-order DC-DC converter, the dual of SEPIC.
 *
 * Inherits converter-level parameters from MAS::Zeta and the Topology
 * interface. Two wiring variants are supported via schema flags:
 *
 *   * V1 — uncoupled (default): L1 (magnetizing) and L2 (output filter)
 *          are independent inductors.
 *   * V2 — coupled-inductor (coupledInductor=true): L1 and L2 share one
 *          core (typically 1:1) with coupling coefficient
 *          `couplingCoefficient`. Ripple-steering per TI SLYT411 — the
 *          condition is identical to coupled-SEPIC because L1 and L2 see
 *          the same voltage in both subintervals.
 *
 * Optional `synchronousRectifier` replaces the catch diode with a
 * low-side MOSFET (LT8711-style). The high-side switch is modelled as
 * a PFET in v1 (NFET-bootstrap variant deferred — see the accuracy
 * disclaimer below).
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW
 * =====================================================================
 *
 *               Q1 (HIGH-SIDE)         Cc
 *   Vin ───┤S        D├────┬───────┤├───────┬─────▢▢▢───┬─────┬─── Vout > 0
 *               │              │             │   L2      │     │
 *               │              │             │           │     │
 *               │             ▽D1            │           Co   Rload
 *               │             (cathode       │           │     │
 *               │              to Cc/Q1      │           │     │
 *               │              drain)        │           │     │
 *              GND            GND            ▢▢▢ L1     GND   GND
 *                                             │
 *                                            GND
 *
 * Switch Q1 ON  : VL1 = +Vin,   VL2 = +Vin (both rising — Cc=Vin steady)
 *                 Diode D1 reverse-biased; iL1+iL2 flows in Q1.
 * Switch Q1 OFF : VL1 = -Vo,    VL2 = -Vo (falling)
 *                 Diode conducts iL1+iL2 to GND.
 *
 * Volt-second balance on either inductor → Vo/Vin = D/(1-D).
 * Charge-balance on Cc → I_Cc,avg = 0; VCc_steady = Vin.
 *
 * Distinctive vs SEPIC: the LC filter (L2 + Co) is on the OUTPUT side,
 * so output ripple is small (triangular, set by L2/Co); input current
 * is pulsed through Q1 (large ripple). DC bias allocation is OPPOSITE
 * of SEPIC: in Zeta L1 (magnetizing) carries the large
 * Iout·D/(1-D) bias and L2 (filter) carries Iout — note SEPIC has L1
 * carrying Iin = Iout·D/(1-D) too, so the *numbers* match SEPIC but
 * the *physical role* of L1 vs L2 is swapped.
 *
 * =====================================================================
 * KEY EQUATIONS (CCM, ideal)
 * =====================================================================
 *
 *   M(D)        = +D / (1 - D)                              [TI SLYT372, Erickson §6]
 *   D           = (Vo + Vd) / (Vin·η + Vo + Vd)             [TI SLVAFJ6 Eq. 6]
 *   VCc         = Vin                                        [TI SLYT372]
 *   VS,peak     = Vin + Vo                                   [TI SLVAFJ6 Eq. 7]
 *   VD,peak     = Vin + Vo                                   [TI SLVAFJ6 Eq. 8]
 *   IL1,avg     = Iout · D / (1 - D)                         [TI SLYT372]
 *   IL2,avg     = Iout                                       [TI SLYT372]
 *   ΔIL1       = Vin · D / (L1 · fsw)
 *   ΔIL2       = Vin · D / (L2 · fsw)
 *   IS,peak     = IL1,avg + IL2,avg + (ΔIL1 + ΔIL2)/2
 *   ICc,rms     ≈ Iout · √(D / (1 - D))                     [TI SLYT372]
 *   Cc ≥ Iout·D / (ΔVCc · fsw)                              [TI SLYT372]
 *   Co ≥ ΔIL2 / (8 · fsw · ΔVo)                              [Erickson §6, buck filter]
 *
 *   K(D)        = 2 · Le · fsw / R,    Le = L1·L2/(L1+L2)
 *   Kcrit(D)    = (1 - D)²
 *   CCM if K > Kcrit, else DCM (Erickson §5/6; MDPI Sensors 21/22/7434)
 *
 * Coupled-inductor (V2) zero-input-ripple condition (TI SLYT411):
 *   k · √(L2/L1) = 1.
 *
 * =====================================================================
 * REFERENCES
 * =====================================================================
 *
 * [1] TI SLYT372, "Designing DC/DC Converters Based on the ZETA Topology".
 * [2] TI SLVAFJ6, "How to Approach a Power-Supply Design Pt 4: Zeta Equations".
 * [3] TI SNVA608, "LM5085 Designing Non-Inverting Buck-Boost (Zeta) with a Buck P-FET".
 * [4] TI PMP9581 reference designs (REVB coupled, REVC uncoupled).
 * [5] TI SLYT411, "Benefits of a coupled-inductor SEPIC" (applies identically to coupled Zeta).
 * [6] ADI, "The Low Output Voltage Ripple Zeta DC/DC Converter Topology".
 * [7] ADI LT8711 datasheet (multi-topology controller).
 * [8] Erickson & Maksimović, "Fundamentals of Power Electronics" 3rd ed., Ch. 5/6/12.
 * [9] MDPI Sensors 21/22/7434, SEPIC/Cuk/Zeta DCM modelling.
 *
 * =====================================================================
 * ACCURACY DISCLAIMER
 * =====================================================================
 *
 * Idealised first-order model. Cc/Co ESR added as a small series R only.
 * MOSFET Coss(Vds), diode reverse-recovery, and core-loss resistance are
 * not modelled. DCM has multiple sub-modes (DCM1/DCM2 per MDPI Sensors
 * 21/22/7434); v1 detects DCM via K vs Kcrit and flags lastIsCcm=false
 * but emits CCM-shape waveforms (deep-DCM waveform reconstruction is a
 * known follow-up). High-side gate-driver bootstrap dynamics are NOT
 * represented — the SPICE netlist drives the switch with an ideal
 * voltage-controlled switch (parity with SEPIC's low-side switch).
 * NFET-bootstrap high-side variant (vs PFET) is deferred to v2.
 * ZVS / soft-switching N/A — Zeta is hard-switched.
 *
 * =====================================================================
 * DISAMBIGUATION
 * =====================================================================
 *
 *   * Zeta vs SEPIC — SEPIC has the LC filter on the INPUT side
 *                     (low-side switch, low input ripple, RHPZ in CCM);
 *                     Zeta has it on the OUTPUT side (high-side switch,
 *                     low output ripple, NO RHPZ in CCM). Same M=D/(1-D).
 *   * Zeta vs Cuk    — Cuk is inverting (Vo<0); Zeta is non-inverting.
 *   * Zeta vs Buck-Boost (non-isolated) — buck-boost inverts and is
 *                     2nd-order; Zeta is 4th-order, non-inverting.
 *   * Zeta vs Flyback — flyback uses an isolation transformer; Zeta is
 *                     non-isolated. The "isolated Zeta" variant
 *                     (out of scope) replaces L1 with a transformer.
 */
class Zeta : public MAS::Zeta, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 50;

    // Maximum operating duty cycle. Zeta CCM duty is identical to
    // SEPIC: D = (Vo+Vd)/(Vin*eta + Vo+Vd); for low Vin or high Vo
    // this approaches 1 and M = D/(1-D) blows up. Previously hard-coded
    // to 0.95 in calculate_duty_cycle; now configurable. Mirrors
    // Cuk (757d50cb), Sepic (9725ba41).
    std::optional<double> maximumDutyCycle = 0.95;

    // Internal sizing rules-of-thumb for L2 / Cc / Co (V1 only).
    double l2RipplePct = 0.30;        // ΔIL2 / IL2,avg target
    double ccRipplePct = 0.05;        // ΔVCc / VCc     target
    double coRipplePct = 0.01;        // ΔVo  / Vo      target

protected:
    // ---- Per-OP analytical diagnostics ----
    mutable double lastDutyCycle = 0.0;
    mutable double lastConversionRatio = 0.0;          // M(D) = +D/(1-D)
    mutable double lastCouplingCapVoltage = 0.0;       // VCc = Vin
    mutable double lastInputInductorAverage = 0.0;     // IL1,avg
    mutable double lastOutputInductorAverage = 0.0;    // IL2,avg = Iout
    mutable double lastInputInductorRipple = 0.0;      // ΔIL1
    mutable double lastOutputInductorRipple = 0.0;     // ΔIL2
    mutable double lastSwitchPeakVoltage = 0.0;        // VS,peak = Vin + Vo
    mutable double lastSwitchPeakCurrent = 0.0;
    mutable double lastDiodePeakReverseVoltage = 0.0;
    mutable double lastDiodePeakCurrent = 0.0;
    mutable double lastCouplingCapRmsCurrent = 0.0;
    mutable double lastDcmK = 0.0;
    mutable double lastDcmKcrit = 0.0;
    mutable bool   lastIsCcm = true;
    mutable double lastSizedL2 = 0.0;
    mutable double lastSizedCc = 0.0;
    mutable double lastSizedCo = 0.0;
    mutable double lastOutputVoltageRipple = 0.0;      // small (L2 in series with Co)
    mutable double lastInputCurrentRipple  = 0.0;      // large (Q1 pulsed input)

    // ---- Extra-component waveforms ----
    mutable std::vector<Waveform> extraL2VoltageWaveforms;
    mutable std::vector<Waveform> extraL2CurrentWaveforms;
    mutable std::vector<Waveform> extraCcVoltageWaveforms;
    mutable std::vector<Waveform> extraCcCurrentWaveforms;
    mutable std::vector<Waveform> extraCoVoltageWaveforms;
    mutable std::vector<Waveform> extraCoCurrentWaveforms;

public:
    bool _assertErrors = false;

    Zeta(const json& j);
    Zeta() {};

    MAS::Topologies topology_kind() const override { return MAS::Topologies::ZETA_CONVERTER; }

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }

    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    std::optional<double> get_maximum_duty_cycle() const { return maximumDutyCycle; }
    void set_maximum_duty_cycle(std::optional<double> value) { this->maximumDutyCycle = value; }

    // ---- Per-OP diagnostic accessors ----
    double get_last_duty_cycle()                  const { return lastDutyCycle; }
    double get_last_conversion_ratio()            const { return lastConversionRatio; }
    double get_last_coupling_cap_voltage()        const { return lastCouplingCapVoltage; }
    double get_last_input_inductor_average()      const { return lastInputInductorAverage; }
    double get_last_output_inductor_average()     const { return lastOutputInductorAverage; }
    double get_last_input_inductor_ripple()       const { return lastInputInductorRipple; }
    double get_last_output_inductor_ripple()      const { return lastOutputInductorRipple; }
    double get_last_switch_peak_voltage()         const { return lastSwitchPeakVoltage; }
    double get_last_switch_peak_current()         const { return lastSwitchPeakCurrent; }
    double get_last_diode_peak_reverse_voltage()  const { return lastDiodePeakReverseVoltage; }
    double get_last_diode_peak_current()          const { return lastDiodePeakCurrent; }
    double get_last_coupling_cap_rms_current()    const { return lastCouplingCapRmsCurrent; }
    double get_last_dcm_k()                       const { return lastDcmK; }
    double get_last_dcm_kcrit()                   const { return lastDcmKcrit; }
    bool   get_last_is_ccm()                      const { return lastIsCcm; }
    double get_last_sized_l2()                    const { return lastSizedL2; }
    double get_last_sized_cc()                    const { return lastSizedCc; }
    double get_last_sized_co()                    const { return lastSizedCo; }
    double get_last_output_voltage_ripple()       const { return lastOutputVoltageRipple; }
    double get_last_input_current_ripple()        const { return lastInputCurrentRipple; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, const TopologyExcitation& outputOperatingPoint, double inductanceL1);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    /**
     * @brief Emit ancillary components (L2 inductor, Cc coupling cap,
     *        Co output cap) as MAS::Inputs / CAS::Inputs entries.
     */
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode = ExtraComponentsMode::IDEAL,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // ---- Static analytical helpers ----
    static double calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency, double maximumDutyCycle = 0.95);
    static double calculate_conversion_ratio(double dutyCycle);                 // +D/(1-D)
    static double calculate_l1_min(double inputVoltage, double dutyCycle, double deltaIL1, double switchingFrequency);
    static double calculate_l2_min(double inputVoltage, double dutyCycle, double deltaIL2, double switchingFrequency);
    static double calculate_cc_min(double outputCurrent, double dutyCycle, double deltaVCc, double switchingFrequency);
    static double calculate_dcm_K(double L1, double L2, double switchingFrequency, double loadResistance);

    /**
     * @brief Generate an ngspice circuit for this Zeta (V1 uncoupled
     *        or V2 coupled, with optional synchronous rectifier).
     *
     * High-side switch implemented with an ideal SPICE switch driven
     * between Vin and Vsw nodes (PFET-equivalent). Includes IC
     * pre-charge for L1, L2, Cc (=Vin), Co (=Vo), GEAR integration,
     * R-snubber across the switch, and ESR on Cc/Co.
     */
    std::string generate_ngspice_circuit(
        double inductanceL1,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);

    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductanceL1);
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        double inductanceL1,
        size_t numberOfPeriods = 2);
};

class AdvancedZeta : public Zeta {
private:
    double desiredInductance;

public:
    AdvancedZeta() = default;
    ~AdvancedZeta() = default;

    AdvancedZeta(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }
};


void from_json(const json & j, AdvancedZeta & x);
void to_json(json & j, const AdvancedZeta & x);


inline void from_json(const json & j, AdvancedZeta& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<TopologyExcitation>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedZeta & x) {
    j = json::object();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumSwitchCurrent"] = x.get_maximum_switch_current();
    j["operatingPoints"] = x.get_operating_points();
    j["desiredInductance"] = x.get_desired_inductance();
}
} // namespace OpenMagnetics
