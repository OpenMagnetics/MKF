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
 * @brief SEPIC (Single-Ended Primary-Inductor Converter) — non-inverting,
 *        fourth-order, capacitive-coupling DC-DC converter.
 *
 * Inherits converter-level parameters from MAS::Sepic and the Topology
 * interface. Two wiring variants are supported via schema flags:
 *
 *   * V1 — uncoupled (default): L1 and L2 are independent inductors.
 *   * V2 — coupled-inductor (coupledInductor=true): L1 and L2 share one
 *          core (typically 1:1) with coupling coefficient
 *          `couplingCoefficient`. Ripple-steering per TI SLYT411 reduces
 *          input-current ripple.
 *
 * Optional `synchronousRectifier`: replaces the catch diode with a
 * low-side MOSFET (LTC1871-style) — common at low Vout.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW
 * =====================================================================
 *
 *           L1                 Cs
 *   Vin ───▢▢▢───┬─────────┤├─────────┬───▢▢▢───┬─────┬─── Vout > 0
 *                │                     │   L2    │     │
 *                │                     │         │     │
 *                S                     ▽D       Co   Rload
 *                │                     │         │     │
 *               GND                   GND       GND   GND
 *
 * Switch S ON  : VL1 = +Vin,   VL2 = +VCs ≈ +Vin (rising)
 *                Diode reverse-biased; iL1+iL2 flows in S.
 * Switch S OFF : VL1 = -Vo,    VL2 = -Vo (falling)
 *                Diode conducts iL1+iL2 to load+Co.
 *
 * Volt-second balance on either inductor → Vo/Vin = D/(1-D).
 * Charge-balance on Cs → I_Cs,avg = 0; VCs_steady = Vin.
 *
 * =====================================================================
 * KEY EQUATIONS (CCM, ideal)
 * =====================================================================
 *
 *   M(D)        = +D / (1 - D)                              [TI SNVA168E, Erickson §6.4]
 *   D           = (Vo + Vd) / (Vin·η + Vo + Vd)             [TI SNVA168E, with η, Vd]
 *   VCs         = Vin                                        [Wikipedia "SEPIC"]
 *   VS,peak     = Vin + Vo                                   [TI SNVA168E]
 *   VD,peak     = Vin + Vo                                   [TI SNVA168E]
 *   IL1,avg     = Iout · D / (1 - D)        (= Iin)         [TI SNVA168E]
 *   IL2,avg     = Iout                                       [TI SNVA168E]
 *   ΔIL1       = Vin · D / (L1 · fsw)
 *   ΔIL2       = Vin · D / (L2 · fsw)
 *   IS,peak     = IL1,avg + IL2,avg + (ΔIL1 + ΔIL2)/2
 *   ICs,rms     ≈ Iout · √(D / (1 - D))                     [TI SLYT309]
 *   Cs ≥ Iout·D / (ΔVCs · fsw)                              [TI SNVA168E]
 *
 *   K(D)        = 2 · Le · fsw / R,    Le = L1·L2/(L1+L2)
 *   Kcrit(D)    = (1 - D)²
 *   CCM if K > Kcrit, else DCM (Erickson §5.7)
 *
 * Coupled-inductor (V2) zero-input-ripple condition (TI SLYT411):
 *   k · √(L2/L1) = 1   (with k = couplingCoefficient).
 *
 * =====================================================================
 * REFERENCES
 * =====================================================================
 *
 * [1] TI SNVA168E (AN-1484), "Designing A SEPIC Converter".
 * [2] TI SLYT309, "Designing DC/DC converters based on SEPIC topology".
 * [3] TI SLYT411, "Benefits of a coupled-inductor SEPIC".
 * [4] TI SLUA158, "Versatile Low-Power SEPIC Converter".
 * [5] TI TIDA-00781, 12 W SEPIC reference design (12 V → 12 V @ 1 A).
 * [6] ADI LTC1871 datasheet (boost / SEPIC controller, synchronous
 *     low-side rectifier examples).
 * [7] Erickson & Maksimović, "Fundamentals of Power Electronics" 3rd
 *     ed. (2020), §5.7 (DCM SEPIC), §6.4 (gain derivation).
 * [8] Würth Elektronik ANP135, "SEPIC with coupled and uncoupled inductors".
 *
 * =====================================================================
 * ACCURACY DISCLAIMER
 * =====================================================================
 *
 * Idealised first-order model. Coupling-cap ESR is added as a small
 * series R only; Cs frequency-dependent loss is not modelled. MOSFET
 * Coss(Vds), diode reverse-recovery, gate-driver propagation delay,
 * and core-loss resistance are not modelled. DCM is detected and
 * flagged via lastIsCcm=false (CCM-shape waveforms are emitted in DCM
 * — full piecewise-linear three-sub-interval reconstruction is a
 * known follow-up).  ZVS / soft-switching N/A — SEPIC is hard-switched.
 *
 * =====================================================================
 * DISAMBIGUATION
 * =====================================================================
 *
 *   * SEPIC vs Cuk     — Cuk is inverting (Vo<0); SEPIC is non-inverting.
 *   * SEPIC vs Zeta    — Zeta has Cs in series with L1's output side
 *                        and needs a high-side switch driver.
 *   * SEPIC vs Buck-Boost (non-isolated) — buck-boost inverts and is
 *                        2nd-order; SEPIC is 4th-order, non-inverting.
 *   * SEPIC vs Flyback — flyback uses an isolation transformer with
 *                        primary-side energy storage; SEPIC stores in L1
 *                        and transfers via Cs. Often interchangeable
 *                        <50 W; SEPIC wins when isolation is unneeded
 *                        and input ripple matters.
 */
class Sepic : public MAS::Sepic, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 50;

    // Maximum operating duty cycle. SEPIC CCM duty is
    // D = (Vo+Vd)/(Vin*eta + Vo+Vd); for low Vin or high Vo this
    // approaches 1, where M = D/(1-D) blows up and regulation
    // collapses. Previously hard-coded to 0.95 in calculate_duty_cycle;
    // now configurable via this field. Mirrors Cuk (757d50cb).
    std::optional<double> maximumDutyCycle = 0.95;

    // Internal sizing rules-of-thumb for L2 / Cs / Co (V1 only).
    double l2RipplePct = 0.30;        // ΔIL2 / IL2,avg target
    double csRipplePct = 0.05;        // ΔVCs / VCs     target
    double coRipplePct = 0.01;        // ΔVo  / Vo      target

protected:
    // ---- Per-OP analytical diagnostics ----
    mutable double lastDutyCycle = 0.0;
    mutable double lastConversionRatio = 0.0;          // M(D) = +D/(1-D)
    mutable double lastCouplingCapVoltage = 0.0;       // VCs = Vin
    mutable double lastInputInductorAverage = 0.0;     // IL1,avg = Iin
    mutable double lastOutputInductorAverage = 0.0;    // IL2,avg = Iout
    mutable double lastInputInductorRipple = 0.0;      // ΔIL1
    mutable double lastOutputInductorRipple = 0.0;     // ΔIL2
    mutable double lastSwitchPeakVoltage = 0.0;        // VS,peak = Vin + Vo
    mutable double lastSwitchPeakCurrent = 0.0;        // ≈ IL1+IL2 + (ΔIL1+ΔIL2)/2
    mutable double lastDiodePeakReverseVoltage = 0.0;  // = VS,peak
    mutable double lastDiodePeakCurrent = 0.0;         // = IS,peak
    mutable double lastCouplingCapRmsCurrent = 0.0;    // Iout·√(D/(1-D))
    mutable double lastDcmK = 0.0;                     // K(D)
    mutable double lastDcmKcrit = 0.0;                 // (1-D)²
    mutable bool   lastIsCcm = true;                   // K > Kcrit
    mutable double lastSizedL2 = 0.0;                  // internally-sized L2
    mutable double lastSizedCs = 0.0;                  // internally-sized Cs
    mutable double lastSizedCo = 0.0;                  // internally-sized Co

    // ---- Per-OP diagnostic vectors (one entry per V_in × OP iteration) ----
    mutable std::vector<std::string> perOpName;
    mutable std::vector<double>      perOpDutyCycle;
    mutable std::vector<double>      perOpConversionRatio;
    mutable std::vector<double>      perOpCouplingCapVoltage;
    mutable std::vector<double>      perOpInputInductorAverage;
    mutable std::vector<double>      perOpOutputInductorAverage;
    mutable std::vector<double>      perOpInputInductorRipple;
    mutable std::vector<double>      perOpOutputInductorRipple;
    mutable std::vector<double>      perOpSwitchPeakVoltage;
    mutable std::vector<double>      perOpSwitchPeakCurrent;
    mutable std::vector<double>      perOpDiodePeakReverseVoltage;
    mutable std::vector<double>      perOpDiodePeakCurrent;
    mutable std::vector<double>      perOpCouplingCapRmsCurrent;
    mutable std::vector<bool>        perOpIsCcm;
    mutable std::vector<double>      perOpSizedCs;
    mutable std::vector<double>      perOpSizedCo;

    // ---- Extra-component waveforms (filled in process_operating_points_for_input_voltage,
    //      consumed by get_extra_components_inputs).  Cleared at the start of
    //      process_operating_points so a re-run does not accumulate. ----
    mutable std::vector<Waveform> extraL2VoltageWaveforms;
    mutable std::vector<Waveform> extraL2CurrentWaveforms;
    mutable std::vector<Waveform> extraCsVoltageWaveforms;
    mutable std::vector<Waveform> extraCsCurrentWaveforms;
    mutable std::vector<Waveform> extraCoVoltageWaveforms;
    mutable std::vector<Waveform> extraCoCurrentWaveforms;

public:
    bool _assertErrors = false;

    Sepic(const json& j);
    Sepic() {};

    MAS::Topology topology_kind() const override { return MAS::Topology::SEPIC_CONVERTER; }

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
    double get_last_sized_cs()                    const { return lastSizedCs; }
    double get_last_sized_co()                    const { return lastSizedCo; }

    // ---- Per-OP vector accessors ----
    const std::vector<std::string>& get_per_op_name()                       const { return perOpName; }
    const std::vector<double>&      get_per_op_duty_cycle()                 const { return perOpDutyCycle; }
    const std::vector<double>&      get_per_op_conversion_ratio()           const { return perOpConversionRatio; }
    const std::vector<double>&      get_per_op_coupling_cap_voltage()       const { return perOpCouplingCapVoltage; }
    const std::vector<double>&      get_per_op_input_inductor_average()     const { return perOpInputInductorAverage; }
    const std::vector<double>&      get_per_op_output_inductor_average()    const { return perOpOutputInductorAverage; }
    const std::vector<double>&      get_per_op_input_inductor_ripple()      const { return perOpInputInductorRipple; }
    const std::vector<double>&      get_per_op_output_inductor_ripple()     const { return perOpOutputInductorRipple; }
    const std::vector<double>&      get_per_op_switch_peak_voltage()        const { return perOpSwitchPeakVoltage; }
    const std::vector<double>&      get_per_op_switch_peak_current()        const { return perOpSwitchPeakCurrent; }
    const std::vector<double>&      get_per_op_diode_peak_reverse_voltage() const { return perOpDiodePeakReverseVoltage; }
    const std::vector<double>&      get_per_op_diode_peak_current()         const { return perOpDiodePeakCurrent; }
    const std::vector<double>&      get_per_op_coupling_cap_rms_current()   const { return perOpCouplingCapRmsCurrent; }
    const std::vector<bool>&        get_per_op_is_ccm()                     const { return perOpIsCcm; }
    const std::vector<double>&      get_per_op_sized_cs()                   const { return perOpSizedCs; }
    const std::vector<double>&      get_per_op_sized_co()                   const { return perOpSizedCo; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, const TopologyExcitation& outputOperatingPoint, double inductanceL1);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    /**
     * @brief Emit ancillary components (L2 inductor, Cs coupling cap,
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
    static double calculate_cs_min(double outputCurrent, double dutyCycle, double deltaVCs, double switchingFrequency);
    static double calculate_dcm_K(double L1, double L2, double switchingFrequency, double loadResistance);

    /**
     * @brief Generate an ngspice circuit for this SEPIC (V1 uncoupled
     *        or V2 coupled, with optional synchronous rectifier).
     *
     * Includes IC pre-charge for L1, L2, Cs (=Vin), Co (=Vo), GEAR
     * integration, R-snubber across the switch, and ESR on Cs/Co.
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

class AdvancedSepic : public Sepic {
private:
    double desiredInductance;

public:
    AdvancedSepic() = default;
    ~AdvancedSepic() = default;

    AdvancedSepic(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }
};


void from_json(const json & j, AdvancedSepic & x);
void to_json(json & j, const AdvancedSepic & x);


inline void from_json(const json & j, AdvancedSepic& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<TopologyExcitation>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedSepic & x) {
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
