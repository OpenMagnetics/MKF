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
 * @brief Weinberg — current-fed, push-pull-derivative, isolated, boost-capable
 *        DC-DC converter (ESA Weinberg 1974; Weinberg & Schreuders PESC 1985).
 *
 * Inherits converter-level parameters from MAS::Weinberg and the Topology
 * interface. Two primary-side variants are supported via the `variant`
 * schema flag:
 *
 *   * V1 — `classic` (default): center-tapped push-pull primary (2
 *          switches Q1/Q2 with 180° phase + (2D−1)·Tsw/2 overlap when
 *          D > 0.5) plus per-switch energy-recovery diodes D3a/D3b
 *          clamping the leakage spike to Vin. Switch V_DS = 2·Vin/(1−D).
 *   * V2 — `bridge`: H-bridge primary (4 switches in two diagonal
 *          pairs, complementary-PWM driven). No D3 needed since the
 *          switches handle the leakage path. Switch V_DS = Vin (halved
 *          vs V1).
 *
 * Optional `synchronousRectifier` replaces the secondary CT-FW diodes
 * (D_pos / D_neg) with active SR MOSFETs (S_pos / S_neg) gated
 * complementary to the primary switches.
 *
 * Magnetic structure (Tier 3 simplified — mirrors SEPIC/Zeta clean
 * shape rather than emitting the full 4-winding xfmr): the main
 * magnetic returned by `process_design_requirements` is a TWO-WINDING
 * transformer (combined primary, combined secondary) with turns ratio
 * n = Np_total / Ns_total. The center-tapping (CT primary in V1, CT
 * secondary in both variants) is implied in the SPICE netlist only
 * (Lpri_a + Lpri_b on the primary in V1; Lsec_a + Lsec_b on the
 * secondary always; full 4-way K coupling at 0.99). The input coupled
 * inductor L1 (1:1, gapped, two windings on the same core) and the
 * output capacitor Co flow through `get_extra_components_inputs`.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW (V1 classic)
 * =====================================================================
 *
 *             L1a              Np_a          D_pos
 *  Vin ──┳──[~~~]───┳────────[==o==]────┬───┤▶├────┐
 *        │  k=0.999 │                   │          │
 *        │  L1b     │                   │ T 4-w    │
 *        │  [~~~]   │                   │ k=0.99   │
 *        └──────────┴──Np_b──[==o==]────┘          │
 *                        │     │                    │
 *                       Q1    Q2  ── center-tap ──→ vout ── Co ── Rload ── 0
 *                        │     │                    │
 *                     D3a/D3b  clamp to Vin         │
 *                        │     │                   D_neg
 *                       GND   GND  ◀── Ns_b ────────┘
 *
 * Modulation: Q1 phase = 0; Q2 phase = Tsw/2 (180° offset). Each
 * switch ON for D·Tsw within its own half-period. When D > 0.5 the
 * pulses overlap by (2D−1)·Tsw/2 per period.
 *
 * =====================================================================
 * KEY EQUATIONS (CCM, 1:1 input coupled-L, n = Np/Ns)
 * =====================================================================
 *
 *   Three operating regimes:
 *     Buck-like     (D < 0.5):    M = Vout/Vin = D / n
 *     Boundary      (D = 0.5):    M = 1 / n
 *     Boost-like    (D > 0.5):    M = 1 / (2 · n · (1 − D))           [Weinberg PESC 1985]
 *
 *   Overlap fraction              : max(0, 2D − 1)
 *   Switch peak voltage (V1)      : V_Q,pk = 2·Vin / (1 − D)          [clamped by D3]
 *   Switch peak voltage (V2)      : V_Q,pk = Vin                       [bridge]
 *   Iin_avg                        : Iout / (η · M)
 *   Switch peak current            : I_Q,pk ≈ Iin + ΔI_L1/2  (non-overlap)
 *                                            = Iin/2 + ΔI_L1/2  (overlap)
 *   Diode peak reverse voltage    : V_D,pk = 2·Vout                    [CT FW rectifier]
 *   Diode average current         : I_D,avg = Iout/2                   [per diode]
 *   Co RMS current                 : I_Co,rms ≈ Iout · √((1−D)/D)
 *   ω_RHPZ (boost regime)         : ≈ R · (1 − D)² / (n² · L1)         [Erickson §8]
 *
 *   Sizing (CCM, boost regime):
 *     L1   ≥ Vin · (2D − 1) · Tsw / (4 · ΔI_L1)
 *     Co   ≥ Iout · D / (ΔVo · fsw)
 *
 *   K_crit (boost) = 2 · L1 · fsw · (1 − D)² / R    ; CCM if K > K_crit
 *
 * =====================================================================
 * REFERENCES
 * =====================================================================
 *
 * [1] A. H. Weinberg, "A Boost Regulator with a New Energy-Transfer
 *     Principle", ESA/ESTEC Spacecraft Power Conditioning Electronics
 *     Seminar, 1974.
 * [2] A. H. Weinberg & L. Schreuders, "A High-Power High-Voltage DC-DC
 *     Converter for Space Applications", IEEE PESC ~1985.
 * [3] R. C. Larico & I. Barbi, "Three-Phase Weinberg Isolated DC-DC
 *     Converter", IEEE Trans. Ind. Electron., 59(2):888–896, 2012.
 * [4] Yadav & Gowrishankara, "Analysis and Development of Weinberg
 *     Converter with Fault-Tolerant Feature", IEEE 2016.
 * [5] Bhandari et al., "Low-output-current-ripple Weinberg converter
 *     for BDR in GEO satellites", J. Power Electron. 25:591, 2025.
 * [6] Marsala & Olivieri, "Comparison Between Two Current-Fed Push-Pull
 *     DC-DC Converters", IEEE 1996. (Watkins-Johnson disambiguation.)
 * [7] Erickson & Maksimović, "Fundamentals of Power Electronics" 3rd
 *     ed. (2020), §8 (RHP zero in boost-derived topologies).
 *
 * =====================================================================
 * ACCURACY DISCLAIMER
 * =====================================================================
 *
 * Idealised first-order CCM model. SPICE model omits MOSFET Coss(Vds),
 * body-diode reverse-recovery, gate-driver propagation delay,
 * transformer inter-winding capacitance, and core-loss resistance.
 * Predicted η is an upper bound; expect 2–4 % loss vs. measured at
 * full load. Energy-recovery diode D3 is modeled as ideal — actual
 * leakage-spike overshoot may exceed 2·Vin/(1−D) by 10–20 % unless an
 * explicit RCD snubber is added. DCM is detected via K vs Kcrit and
 * flagged (lastIsCcm=false) but emits CCM-shape waveforms — full-DCM
 * waveform reconstruction is a known follow-up. Three-phase Weinberg
 * (Larico-Barbi 2012), V3 bidirectional, V5 ZVS active-clamp, and
 * multi-output Weinberg are out of scope for this v1 add.
 *
 * =====================================================================
 * DISAMBIGUATION
 * =====================================================================
 *
 *   * Weinberg vs Watkins-Johnson — W-J uses a *single* input inductor
 *     and 1 switch; Weinberg uses a *coupled 1:1* input inductor with
 *     energy-recovery diode D3, which intrinsically handles the
 *     overlap interval without needing a freewheel switch. M(D) and
 *     V_Q,pk differ. See Marsala & Olivieri, IEEE 1996.
 *   * Weinberg vs plain current-fed push-pull — plain CFPP needs an
 *     explicit clamp/snubber for the leakage spike; Weinberg's D3 +
 *     1:1 coupled L is the recovery mechanism.
 *   * Weinberg vs PSFB / Phase-Shift Full-Bridge — PSFB is voltage-fed
 *     (capacitive input); Weinberg is current-fed (inductive input).
 *   * Weinberg vs Flyback — Flyback is a single-switch, single-magnetic
 *     boost-derived isolated topology suited to <100 W; Weinberg is
 *     push-pull-derived, two-magnetic, suited to 0.5–5 kW.
 */
class Weinberg : public MAS::Weinberg, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 50;     // Weinberg is slow to settle (large L1)

    // Internal sizing rules-of-thumb for L1 / Co.
    double l1RipplePct = 0.20;          // ΔIL1 / IL1_per_winding,avg target
    double coRipplePct = 0.01;          // ΔVo  / Vo                    target
    double leakageInductanceFraction = 0.001;  // Llk = 0.1 % of Lpri_half — D3 not modelled in V1 SPICE (see Weinberg.cpp), so keep leakage low.

protected:
    // ---- Per-OP analytical diagnostics (Guide §2.4) ----
    mutable double lastDutyCycle = 0.0;
    mutable double lastConversionRatio = 0.0;
    mutable int    lastOperatingRegime = 0;          // 0 = buck, 1 = boundary, 2 = boost
    mutable double lastOverlapFraction = 0.0;        // max(0, 2D-1)
    mutable double lastSwitchPeakVoltage = 0.0;      // V1: 2·Vin/(1-D); V2: Vin
    mutable double lastSwitchPeakCurrent = 0.0;
    mutable double lastDiodePeakReverseVoltage = 0.0; // 2·Vout
    mutable double lastDiodePeakCurrent = 0.0;
    mutable double lastEnergyRecoveryAvgCurrent = 0.0; // I_D3,avg (V1 only)
    mutable double lastInputInductorAverage = 0.0;    // Iin/2 per L1 winding
    mutable double lastInputInductorRipple = 0.0;     // ΔI_L1
    mutable double lastMagnetizingRipple = 0.0;
    mutable double lastFluxImbalanceMargin = 0.0;     // (D_Q1 − D_Q2)/D, surfaced
    mutable double lastRhpZeroFrequency = 0.0;
    mutable double lastDcmK = 0.0;
    mutable double lastDcmKcrit = 0.0;
    mutable bool   lastIsCcm = true;
    mutable double lastSizedL1 = 0.0;
    mutable double lastSizedCo = 0.0;
    mutable double lastOutputVoltageRipple = 0.0;

    // ---- Extra-component waveforms (cleared in process_operating_points) ----
    mutable std::vector<Waveform> extraL1VoltageWaveforms;   // per-winding voltage
    mutable std::vector<Waveform> extraL1CurrentWaveforms;   // per-winding current
    mutable std::vector<Waveform> extraCoVoltageWaveforms;
    mutable std::vector<Waveform> extraCoCurrentWaveforms;

public:
    bool _assertErrors = false;

    Weinberg(const json& j);
    Weinberg() {};

    MAS::Topologies topology_kind() const override { return MAS::Topologies::WEINBERG_CONVERTER; }

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }

    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    // ---- Per-OP diagnostic accessors ----
    double get_last_duty_cycle()                   const { return lastDutyCycle; }
    double get_last_conversion_ratio()             const { return lastConversionRatio; }
    int    get_last_operating_regime()             const { return lastOperatingRegime; }
    double get_last_overlap_fraction()             const { return lastOverlapFraction; }
    double get_last_switch_peak_voltage()          const { return lastSwitchPeakVoltage; }
    double get_last_switch_peak_current()          const { return lastSwitchPeakCurrent; }
    double get_last_diode_peak_reverse_voltage()   const { return lastDiodePeakReverseVoltage; }
    double get_last_diode_peak_current()           const { return lastDiodePeakCurrent; }
    double get_last_energy_recovery_avg_current()  const { return lastEnergyRecoveryAvgCurrent; }
    double get_last_input_inductor_average()       const { return lastInputInductorAverage; }
    double get_last_input_inductor_ripple()        const { return lastInputInductorRipple; }
    double get_last_magnetizing_ripple()           const { return lastMagnetizingRipple; }
    double get_last_flux_imbalance_margin()        const { return lastFluxImbalanceMargin; }
    double get_last_rhp_zero_frequency()           const { return lastRhpZeroFrequency; }
    double get_last_dcm_k()                        const { return lastDcmK; }
    double get_last_dcm_kcrit()                    const { return lastDcmKcrit; }
    bool   get_last_is_ccm()                       const { return lastIsCcm; }
    double get_last_sized_l1()                     const { return lastSizedL1; }
    double get_last_sized_co()                     const { return lastSizedCo; }
    double get_last_output_voltage_ripple()        const { return lastOutputVoltageRipple; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;

    OperatingPoint process_operating_points_for_input_voltage(
        double inputVoltage,
        const TopologyExcitation& outputOperatingPoint,
        double turnsRatio,
        double magnetizingInductance);

    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    /**
     * @brief Emit ancillary components (L1 input coupled inductor, Co
     *        output cap) as MAS::Inputs / CAS::Inputs entries.
     */
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode = ExtraComponentsMode::IDEAL,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // ---- Static analytical helpers ----
    static double calculate_duty_cycle(double inputVoltage, double outputVoltage,
                                       double turnsRatio, double diodeVoltageDrop,
                                       double efficiency);
    static double calculate_conversion_ratio_boost(double dutyCycle, double turnsRatio);  // 1/(2n(1-D))
    static double calculate_conversion_ratio_buck(double dutyCycle, double turnsRatio);   // D/n
    static double calculate_overlap_fraction(double dutyCycle);                            // max(0, 2D-1)
    static double calculate_switch_peak_voltage_classic(double inputVoltage, double dutyCycle); // 2*Vin/(1-D)
    static double calculate_switch_peak_voltage_bridge(double inputVoltage);               // = Vin
    static double calculate_l1_min(double inputVoltage, double dutyCycle, double deltaIL1, double switchingFrequency);
    static double calculate_dcm_K_boost(double L1, double switchingFrequency, double dutyCycle, double loadResistance);
    static double calculate_rhp_zero_freq(double loadResistance, double dutyCycle, double turnsRatio, double L1);
    static int    detect_operating_regime(double dutyCycle);   // 0 buck, 1 boundary, 2 boost

    /**
     * @brief Generate an ngspice circuit for this Weinberg (V1 classic
     *        push-pull or V2 H-bridge primary, with optional
     *        synchronous-rectifier secondary).
     *
     * Mandatory contents per WEINBERG_PLAN.md §6.2: K_in (L1a, L1b)
     * coupling at 0.999; K_xfmr (Lpri_a, Lpri_b, Lsec_a, Lsec_b) at
     * 0.99 (NEVER 1.0 — singular factorisation); explicit Llk on every
     * transformer winding (2 % of self); D3a/D3b energy-recovery
     * diodes per switch drain to vin (V1 only); Vab probe; snubbers on
     * every switch and rectifier diode; IC pre-charge for L1a/L1b/Co;
     * GEAR + ITL1=500 + UIC.
     */
    std::string generate_ngspice_circuit(
        double turnsRatio,
        double magnetizingInductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);

    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        double turnsRatio,
        double magnetizingInductance);

    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        double turnsRatio,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);
};

class AdvancedWeinberg : public Weinberg {
private:
    double desiredInductance;
    double desiredTurnsRatio;

public:
    AdvancedWeinberg() = default;
    ~AdvancedWeinberg() = default;

    AdvancedWeinberg(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    const double & get_desired_turns_ratio() const { return desiredTurnsRatio; }
    double & get_mutable_desired_turns_ratio() { return desiredTurnsRatio; }
    void set_desired_turns_ratio(const double & value) { this->desiredTurnsRatio = value; }
};


void from_json(const json & j, AdvancedWeinberg & x);
void to_json(json & j, const AdvancedWeinberg & x);


inline void from_json(const json & j, AdvancedWeinberg& x) {
    x.set_coupling_coefficient_input(get_stack_optional<double>(j, "couplingCoefficientInput"));
    x.set_coupling_coefficient_main(get_stack_optional<double>(j, "couplingCoefficientMain"));
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<TopologyExcitation>>());
    x.set_synchronous_rectifier(get_stack_optional<bool>(j, "synchronousRectifier"));
    x.set_variant(get_stack_optional<MAS::Variant>(j, "variant"));
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
    x.set_desired_turns_ratio(j.at("desiredTurnsRatio").get<double>());
}

inline void to_json(json & j, const AdvancedWeinberg & x) {
    j = json::object();
    j["couplingCoefficientInput"] = x.get_coupling_coefficient_input();
    j["couplingCoefficientMain"] = x.get_coupling_coefficient_main();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumSwitchCurrent"] = x.get_maximum_switch_current();
    j["operatingPoints"] = x.get_operating_points();
    j["synchronousRectifier"] = x.get_synchronous_rectifier();
    j["variant"] = x.get_variant();
    j["desiredInductance"] = x.get_desired_inductance();
    j["desiredTurnsRatio"] = x.get_desired_turns_ratio();
}
} // namespace OpenMagnetics
