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
 * @brief Cuk (inverting, fourth-order) DC-DC Converter — non-isolated V1.
 *
 * Inherits converter-level parameters from MAS::Cuk and the Topology
 * interface. This class implements the non-isolated, asynchronous-diode
 * (or synchronous), CCM Cuk topology. Coupled-inductor (V2) and isolated
 * (V3) variants are deferred per CUK_PLAN.md §11.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW (non-isolated Cuk, V1)
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
 * Output Vo is *negative* (M(D) = -D/(1-D)). Iout reported as a positive
 * scalar flowing OUT of the negative terminal toward GND (passive sign
 * convention for the load).
 *
 * The "primary" magnetic returned by process_design_requirements is
 * the input inductor L1 (single winding, one isolation side). The
 * output inductor L2 and coupling cap C1 are sized internally from
 * rules-of-thumb (ΔIL2_pct = 0.30, ΔVC1_pct = 0.05) and instantiated
 * inside the SPICE netlist; they are NOT exposed as DesignRequirements
 * in this initial revision.
 *
 * =====================================================================
 * KEY EQUATIONS (CCM, ideal — Erickson-Maksimović 3rd ed. §2.4)
 * =====================================================================
 *
 *   M(D)        = -D / (1 - D)                                 [E1]
 *   VC1         = Vin / (1 - D) = Vin + |Vo|                   [E2]
 *   VS,peak     = VD,peak = Vin + |Vo|                         [E3,E4]
 *   IL1,avg     = D · Iout / (1 - D)            (= Iin)        [E5a]
 *   IL2,avg     = Iout                                          [E5b]
 *   IC1,rms     ≈ √(D·(1-D)) · (Iin + Iout)                    [E10]
 *   ΔIL1       = Vin · D / (L1 · fsw)
 *   ΔIL2       = |Vo| · (1-D) / (L2 · fsw)
 *
 *   K(D)        = 2 · Le · fsw / R,    Le = L1·L2/(L1+L2)
 *   Kcrit(D)    = (1 - D)²
 *   CCM if K > Kcrit, else DCM
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
 * [2] Erickson & Maksimović, "Fundamentals of Power Electronics"
 *     3rd ed. (2020), §2.4, §5, §8, §11.
 * [3] TI LM2611 datasheet (SNOS965F).
 * [4] Simon Bramble, "Inverting DC-DC Converter Design",
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
 * but DCM-specific waveform construction is not yet implemented; the
 * CCM expressions are emitted regardless and the user is warned.
 *
 * =====================================================================
 * POLARITY CONVENTION
 * =====================================================================
 *
 * Vo < 0 in V1 (non-isolated). The user-facing operatingPoints input
 * passes |Vo| as a positive scalar (matching every other MKF
 * non-isolated converter); the internal calculate_duty_cycle and SPICE
 * netlist treat |Vo| as the magnitude. The emitted ConverterWaveforms
 * output_voltage is signed (negative).
 */
class Cuk : public MAS::Cuk, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 50;

    // Internal sizing rules-of-thumb for L2 / C1 / Co (V1 only).
    // These are fixed defaults; future revisions will expose them via
    // schema fields and/or DesignRequirements (see CUK_PLAN.md §13).
    double l2RipplePct = 0.30;        // ΔIL2 / IL2,avg target
    double c1RipplePct = 0.05;        // ΔVC1 / VC1     target
    double coRipplePct = 0.01;        // ΔVo  / |Vo|    target

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

public:
    bool _assertErrors = false;

    Cuk(const json& j);
    Cuk() {};

    MAS::Topologies topology_kind() const override { return MAS::Topologies::CUK_CONVERTER; }

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }

    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

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

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, const BaseOperatingPoint& outputOperatingPoint, double inductanceL1);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    // ---- Static analytical helpers (CUK_PLAN.md §2.13) ----
    static double calculate_duty_cycle(double inputVoltage, double outputVoltageMagnitude, double diodeVoltageDrop, double efficiency);
    static double calculate_conversion_ratio(double dutyCycle);                 // -D/(1-D)
    static double calculate_coupling_cap_voltage(double inputVoltage, double dutyCycle);
    static double calculate_l1_min(double inputVoltage, double dutyCycle, double deltaIL1, double switchingFrequency);
    static double calculate_l2_min(double outputVoltageMagnitude, double dutyCycle, double deltaIL2, double switchingFrequency);
    static double calculate_c1_min(double outputCurrent, double dutyCycle, double deltaVC1, double switchingFrequency);
    static double calculate_dcm_K(double L1, double L2, double switchingFrequency, double loadResistance);
    static double calculate_rhp_zero_frequency(double loadResistance, double dutyCycle, double L2);

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
    x.set_operating_points(j.at("operatingPoints").get<std::vector<BaseOperatingPoint>>());
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
