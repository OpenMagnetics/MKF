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
 * @brief Cuk (inverting, fourth-order) DC-DC Converter.
 *
 * Inherits converter-level parameters from MAS::Cuk and the Topology
 * interface. This class implements all five Cuk topology variants from
 * CUK_PLAN.md:
 *
 *   * V1 — non-isolated, asynchronous (D1) freewheel, CCM (default)
 *   * V2 — coupled-inductor zero-ripple (Erickson §6); enable via
 *          schema flag `coupledInductor` + optional `couplingCoefficient`
 *   * V3 — isolated (Erickson §6, Cuk's 1977 patent); enable via
 *          schema flags `isolated`, `turnsRatio`, optional
 *          `couplingCapacitanceSecondary`
 *   * V4 — synchronous rectifier (S2 replaces D1, complementary PWM,
 *          no dead-time); enable via schema flag `synchronous`
 *   * V5 — bidirectional (S2 + S1 both 4-quadrant); implies V4
 *          synchronous regardless; enable via schema flag `bidirectional`
 *
 * V2 and V3 are mutually exclusive (V2 puts both windings on one core,
 * V3 separates them with galvanic isolation). All other combinations are
 * supported (V3+V4, V3+V5, V2+V4, V2+V5, V1 default).
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW (V1 non-isolated, single-cap coupling)
 * =====================================================================
 *
 *           L1               C1                 L2
 *   Vin ──┳──/\/\─────┳──┤├──┳──/\/\──┳── Vo (negative; load to GND)
 *         │           │      │        │
 *         │          S│     D│       Co
 *         │           │      │        │
 *        Cin         GND    GND      GND ── Rload
 *         │
 *        GND
 *
 * Switch S ON  : VL1 = +Vin,           VL2 = -VC1 + Vo (negative slope)
 * Switch S OFF : VL1 = +Vin - VC1,     VL2 = +Vo      (negative slope)
 *
 * Output Vo is *negative* (M(D) = -D/(1-D) for V1/V2; -D/((1-D)·n) for V3).
 * Iout reported as a positive scalar flowing OUT of the negative terminal
 * toward GND (passive sign convention for the load).
 *
 * V3 isolated extends the middle section to:
 *
 *   ... node_A ──┤Ca├── tx_pri_pos ─┤Lp ⟂ Ls├─ tx_sec_pos ──┤Cb├── node_B ...
 *
 * with K=0.999 mutual coupling between Lp and Ls (Lp = Lm, Ls = Lm/n²).
 *
 * The "primary" magnetic returned by process_design_requirements is:
 *   * V1/V2: input inductor L1 (single winding for V1, 2-winding 1:1 for V2)
 *   * V3:    isolation transformer Lp:Ls (turns_ratios={n}, sides={PRI,SEC})
 *
 * For V1/V2 the output inductor L2, coupling cap C1 and output cap Co
 * are sized internally and emitted via get_extra_components_inputs.
 * For V3 the input choke L1, output inductor L2, primary/secondary
 * coupling caps Ca/Cb and output cap Co are emitted as 5 ancillary
 * entries.
 *
 * =====================================================================
 * KEY EQUATIONS (CCM, ideal — Erickson-Maksimović 3rd ed. §2.4 / §6)
 * =====================================================================
 *
 *   M(D)        = -D / (1 - D)              (V1/V2)
 *   M(D,n)      = -D / ((1 - D) · n)        (V3, n = Np/Ns)
 *   VC1 / VCa   = Vin / (1 - D)
 *   VCb         = |Vo| / D                  (V3 only)
 *   VS,peak     = VCa
 *   VD,peak     = VCb (V3) | VCa (V1/V2)
 *   IL1,avg     = |Vo|·Iout/(Vin·η)         (= Iout·D/(1-D)/n in lossless)
 *   IL2,avg     = Iout
 *   IC1,rms     ≈ √(D·(1-D)) · (Iin + Iout)
 *   ΔIL1       = Vin · D / (L1 · fsw)
 *   ΔIL2       = |Vo| · (1-D) / (L2 · fsw)
 *
 *   K(D)        = 2 · Le · fsw / R,    Le = L1·L2/(L1+L2)
 *   Kcrit(D)    = (1 - D)²
 *   CCM if K > Kcrit, else DCM
 *   |M(D,K)|    = D / √K            (DCM gain magnitude)
 *   D2          = √K                 (second sub-interval, DCM)
 *
 *   ωRHP        ≈ R · (1 - D)² / (D · L2)  (dominant, Erickson §8)
 *   fc,max      ≤ ωRHP / (2π · 5)
 *
 * =====================================================================
 * REFERENCES
 * =====================================================================
 *
 * [1] Slobodan Cuk, "A New Optimum Topology Switching DC-to-DC
 *     Converter", IEEE PESC 1977.
 * [2] Cuk patents US4,257,087 and US4,355,352 (isolated variant).
 * [3] Erickson & Maksimović, "Fundamentals of Power Electronics"
 *     3rd ed. (2020), §2.4, §5, §6, §8, §11.
 * [4] TI LM2611 datasheet (SNOS965F).
 * [5] Simon Bramble, "Inverting DC-DC Converter Design",
 *     http://www.simonbramble.co.uk/
 *
 * =====================================================================
 * ACCURACY DISCLAIMER
 * =====================================================================
 *
 * SPICE model omits MOSFET Coss variation with Vds, body-diode reverse
 * recovery, gate-driver propagation delay, and core-loss resistance.
 * Predicted η is an upper bound; expect 2-5 % loss vs. measured at
 * full load. DCM operation is detected and flagged via lastIsCcm=false
 * and the |M(D,K)|=D/√K diagnostic; full DCM-specific *waveform*
 * synthesis (piecewise-linear three-sub-interval iL1/iL2 reconstruction)
 * remains a known follow-up — the analytical waveforms emitted in DCM
 * are CCM-shape approximations.
 *
 * =====================================================================
 * POLARITY CONVENTION
 * =====================================================================
 *
 * Vo < 0 in all variants. The user-facing operatingPoints input passes
 * |Vo| as a positive scalar (matching every other MKF non-isolated
 * converter); the internal calculate_duty_cycle and SPICE netlist
 * treat |Vo| as the magnitude. The emitted ConverterWaveforms
 * output_voltage is signed (negative).
 */
class Cuk : public MAS::Cuk, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 50;

    // Maximum operating duty cycle. Cuk CCM duty is
    // D = (|Vo|+Vd)*n / (Vin*eta + (|Vo|+Vd)*n); for low Vin or high
    // step-down/up ratio this approaches 1, where the conversion ratio
    // M = -D/(1-D) blows up and regulation collapses. Cuk previously
    // hard-coded a 0.95 ceiling in calculate_duty_cycle(); this field
    // makes it configurable so callers with stricter controller
    // constraints (e.g. maxD = 0.7) get a loud throw rather than
    // running well past their controller's range. Mirrors Flyback
    // (04272d7b), forward family (683e731c), Buck (2c9300c2), Boost
    // (96fdb52a), IsolatedBuck (703bc80e), IsolatedBuckBoost
    // (a158d548).
    std::optional<double> maximumDutyCycle = 0.95;

    // Internal sizing rules-of-thumb for L2 / C1 / Co (V1 only).
    // These are fixed defaults; future revisions will expose them via
    // schema fields and/or DesignRequirements (see CUK_PLAN.md §13).
    double l2RipplePct = 0.30;        // ΔIL2 / IL2,avg target
    double c1RipplePct = 0.05;        // ΔVC1 / VC1     target
    double coRipplePct = 0.01;        // ΔVo  / |Vo|    target

protected:
    // ---- Per-OP analytical diagnostics ----
    mutable double lastDutyCycle = 0.0;
    mutable double lastConversionRatio = 0.0;          // M(D) signed (negative)
    mutable double lastCouplingCapVoltage = 0.0;       // VC1
    mutable double lastInputInductorAverage = 0.0;     // IL1,avg = Iin
    mutable double lastOutputInductorAverage = 0.0;    // IL2,avg = Iout
    mutable double lastInputInductorRipple = 0.0;      // ΔIL1
    mutable double lastOutputInductorRipple = 0.0;     // ΔIL2
    mutable double lastSwitchPeakVoltage = 0.0;        // VS,peak = VC1
    mutable double lastSwitchPeakCurrent = 0.0;        // ≈ IL1,avg+IL2,avg + (ΔIL1+ΔIL2)/2
    mutable double lastDiodePeakReverseVoltage = 0.0;  // = VS,peak
    mutable double lastDiodePeakCurrent = 0.0;         // = IS,peak
    mutable double lastCouplingCapRmsCurrent = 0.0;    // √(D(1-D))·(Iin+Iout)
    mutable double lastRhpZeroFrequency = 0.0;         // fRHP
    mutable double lastRecommendedLoopBandwidth = 0.0; // fRHP/5
    mutable double lastDcmK = 0.0;                     // K(D)
    mutable double lastDcmKcrit = 0.0;                 // (1-D)²
    mutable bool   lastIsCcm = true;                   // K > Kcrit
    mutable double lastSizedL2 = 0.0;                  // internally-sized L2
    mutable double lastSizedC1 = 0.0;                  // internally-sized C1
    mutable double lastSizedCo = 0.0;                  // internally-sized Co
    // V3 isolated diagnostics. Populated only when isolated=true.
    mutable double lastSizedCa = 0.0;                  // primary-side coupling cap
    mutable double lastSizedCb = 0.0;                  // secondary-side coupling cap
    mutable double lastSizedLm = 0.0;                  // transformer magnetizing inductance
    mutable double lastTurnsRatio = 1.0;

    // ---- Per-OP diagnostic vectors ----
    mutable std::vector<std::string> perOpName;
    mutable std::vector<double>  perOpDutyCycle;
    mutable std::vector<double>  perOpConversionRatio;
    mutable std::vector<double>  perOpCouplingCapVoltage;
    mutable std::vector<double>  perOpInputInductorAverage;
    mutable std::vector<double>  perOpOutputInductorAverage;
    mutable std::vector<double>  perOpInputInductorRipple;
    mutable std::vector<double>  perOpOutputInductorRipple;
    mutable std::vector<double>  perOpSwitchPeakVoltage;
    mutable std::vector<double>  perOpSwitchPeakCurrent;
    mutable std::vector<double>  perOpDiodePeakReverseVoltage;
    mutable std::vector<double>  perOpDiodePeakCurrent;
    mutable std::vector<double>  perOpCouplingCapRmsCurrent;
    mutable std::vector<bool>    perOpIsCcm;
    mutable std::vector<double>  perOpSizedCa;
    mutable std::vector<double>  perOpSizedCb;
    mutable std::vector<double>  perOpSizedCo;
    mutable std::vector<double>  perOpRhpZeroFrequency;               // Np/Ns (1.0 in V1/V2)

    // ---- Extra-component waveforms (filled in process_operating_point_for_input_voltage,
    //      consumed by get_extra_components_inputs).  Cleared at the start of
    //      process_operating_points so a re-run does not accumulate. ----
    mutable std::vector<Waveform> extraL1VoltageWaveforms;   // V3 only (V1 reuses primary excitation)
    mutable std::vector<Waveform> extraL1CurrentWaveforms;
    mutable std::vector<Waveform> extraL2VoltageWaveforms;
    mutable std::vector<Waveform> extraL2CurrentWaveforms;
    mutable std::vector<Waveform> extraC1VoltageWaveforms;   // V1/V2: single C1; V3: unused
    mutable std::vector<Waveform> extraC1CurrentWaveforms;
    mutable std::vector<Waveform> extraCaVoltageWaveforms;   // V3 only: primary-side coupling cap
    mutable std::vector<Waveform> extraCaCurrentWaveforms;
    mutable std::vector<Waveform> extraCbVoltageWaveforms;   // V3 only: secondary-side coupling cap
    mutable std::vector<Waveform> extraCbCurrentWaveforms;
    mutable std::vector<Waveform> extraCoVoltageWaveforms;
    mutable std::vector<Waveform> extraCoCurrentWaveforms;

public:
    bool _assertErrors = false;

    Cuk(const json& j);
    Cuk() {};

    MAS::Topologies topology_kind() const override { return MAS::Topologies::CUK_CONVERTER; }

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
    double get_last_rhp_zero_frequency()          const { return lastRhpZeroFrequency; }
    double get_last_recommended_loop_bandwidth()  const { return lastRecommendedLoopBandwidth; }
    double get_last_dcm_k()                       const { return lastDcmK; }
    double get_last_dcm_kcrit()                   const { return lastDcmKcrit; }
    bool   get_last_is_ccm()                      const { return lastIsCcm; }
    double get_last_sized_l2()                    const { return lastSizedL2; }
    double get_last_sized_c1()                    const { return lastSizedC1; }
    double get_last_sized_co()                    const { return lastSizedCo; }
    double get_last_sized_ca()                    const { return lastSizedCa; }
    double get_last_sized_cb()                    const { return lastSizedCb; }
    double get_last_sized_lm()                    const { return lastSizedLm; }
    double get_last_turns_ratio()                 const { return lastTurnsRatio; }

    // ---- Per-OP vector accessors ----
    const std::vector<std::string>& get_per_op_name() const { return perOpName; }
    const std::vector<double>& get_per_op_duty_cycle() const { return perOpDutyCycle; }
    const std::vector<double>& get_per_op_conversion_ratio() const { return perOpConversionRatio; }
    const std::vector<double>& get_per_op_coupling_cap_voltage() const { return perOpCouplingCapVoltage; }
    const std::vector<double>& get_per_op_input_inductor_average() const { return perOpInputInductorAverage; }
    const std::vector<double>& get_per_op_output_inductor_average() const { return perOpOutputInductorAverage; }
    const std::vector<double>& get_per_op_input_inductor_ripple() const { return perOpInputInductorRipple; }
    const std::vector<double>& get_per_op_output_inductor_ripple() const { return perOpOutputInductorRipple; }
    const std::vector<double>& get_per_op_switch_peak_voltage() const { return perOpSwitchPeakVoltage; }
    const std::vector<double>& get_per_op_switch_peak_current() const { return perOpSwitchPeakCurrent; }
    const std::vector<double>& get_per_op_diode_peak_reverse_voltage() const { return perOpDiodePeakReverseVoltage; }
    const std::vector<double>& get_per_op_diode_peak_current() const { return perOpDiodePeakCurrent; }
    const std::vector<double>& get_per_op_coupling_cap_rms_current() const { return perOpCouplingCapRmsCurrent; }
    const std::vector<bool>& get_per_op_is_ccm() const { return perOpIsCcm; }
    const std::vector<double>& get_per_op_sized_ca() const { return perOpSizedCa; }
    const std::vector<double>& get_per_op_sized_cb() const { return perOpSizedCb; }
    const std::vector<double>& get_per_op_sized_co() const { return perOpSizedCo; }
    const std::vector<double>& get_per_op_rhp_zero_frequency() const { return perOpRhpZeroFrequency; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, const CukOperatingPoint& outputOperatingPoint, double inductanceL1);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    /**
     * @brief Emit ancillary components (L2 inductor, C1 coupling cap, Co
     *        output cap) as MAS::Inputs / CAS::Inputs entries, populated
     *        from the analytical waveforms stashed during
     *        process_operating_points_for_input_voltage.
     *
     * For Cuk V1, three ancillaries are returned in order: L2 (Inputs),
     * C1 (CAS::Inputs, role=DC_LINK as the closest semantic match for
     * an energy-transfer coupling cap), Co (CAS::Inputs, role=OUTPUT_FILTER).
     */
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode = ExtraComponentsMode::IDEAL,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // ---- Static analytical helpers (CUK_PLAN.md §2.13) ----
    //
    // For V3 isolated, turnsRatio = Np/Ns (Flyback convention). Conversion gain
    // becomes |Vo|/Vin = D / ((1-D) · turnsRatio); reduces to V1 when turnsRatio=1.
    static double calculate_duty_cycle(double inputVoltage, double outputVoltageMagnitude, double diodeVoltageDrop, double efficiency, double turnsRatio = 1.0, double maximumDutyCycle = 0.95);
    static double calculate_conversion_ratio(double dutyCycle);                 // -D/(1-D)
    static double calculate_coupling_cap_voltage(double inputVoltage, double dutyCycle);
    static double calculate_l1_min(double inputVoltage, double dutyCycle, double deltaIL1, double switchingFrequency);
    static double calculate_l2_min(double outputVoltageMagnitude, double dutyCycle, double deltaIL2, double switchingFrequency);
    static double calculate_c1_min(double outputCurrent, double dutyCycle, double deltaVC1, double switchingFrequency);
    static double calculate_dcm_K(double L1, double L2, double switchingFrequency, double loadResistance);
    static double calculate_rhp_zero_frequency(double loadResistance, double dutyCycle, double L2);

    /**
     * @brief Compute the input-choke inductance L1 for an isolated (V3) Cuk.
     *
     * In V1/V2 the magnetizing-inductance design requirement IS L1; in V3 the
     * design requirement is the transformer Lm and L1 is ancillary, sized
     * internally from the same ripple-ratio / max-switch-current rule used by
     * process_design_requirements for V1.  Throws if neither current_ripple_ratio
     * nor maximum_switch_current is provided (no silent fallbacks per
     * CLAUDE.md).
     */
    double compute_l1_for_isolated() const;

    /**
     * @brief Generate an ngspice circuit for this Cuk converter (V1).
     *
     * Uses the calculated design parameters (L1, sized L2, sized C1,
     * sized Co, duty cycle) to create a SPICE netlist that can be
     * simulated. Includes IC pre-charge for L1, L2, C1, Co (mandatory
     * per CUK_PLAN.md §6.2 — without IC, ngspice diverges on the
     * 4th-order LC tank cold-start), GEAR integration, snubbers, and
     * small parasitic R on inductors / ESR on caps.
     *
     * @param inductanceL1 L1 inductance in H
     * @param inputVoltageIndex Which input voltage to use (0=nom, 1=min, 2=max)
     * @param operatingPointIndex Which operating point to simulate
     * @return SPICE netlist string
     */
    std::string generate_ngspice_circuit(
        double inductanceL1,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);

    /**
     * @brief Simulate the Cuk converter and extract operating points.
     *
     * The "primary" winding emitted is the input inductor L1.
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductanceL1);

    /**
     * @brief Simulate and extract topology-level (Vin/Iin/Vout/Iout) waveforms.
     */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        double inductanceL1,
        size_t numberOfPeriods = 2);
};

class AdvancedCuk : public Cuk {
private:
    double desiredInductance;

public:
    AdvancedCuk() = default;
    ~AdvancedCuk() = default;

    AdvancedCuk(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }
};


void from_json(const json & j, AdvancedCuk & x);
void to_json(json & j, const AdvancedCuk & x);


inline void from_json(const json & j, AdvancedCuk& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<CukOperatingPoint>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedCuk & x) {
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
