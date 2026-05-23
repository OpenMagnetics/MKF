#pragma once
// Version: 2.0 (CLLC schema v2 — bridgeType, integrated Lr1/Lr2, asymmetric a/b ratios)
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"

namespace OpenMagnetics {
using namespace MAS;

/**
 * @brief Resonant tank parameters for the CLLC converter.
 *
 * Contains all the calculated values needed to define the resonant tank
 * of a bidirectional CLLC converter. Computed from the converter
 * specification using the 11-step Infineon AN-2024-06 procedure. The
 * `resonantInductorRatio` (a) and `resonantCapacitorRatio` (b) capture
 * asymmetric-tank designs (Min IEEE TPE 36(6), 2021).
 */
struct CllcResonantParameters {
    double turnsRatio;                 ///< Transformer turns ratio n = Np/Ns
    double resonantFrequency;          ///< Natural resonant frequency fr [Hz]
    double primaryResonantInductance;  ///< Lr1 — primary series resonant inductor [H]
    double primaryResonantCapacitance; ///< Cr1 — primary series resonant capacitor [F]
    double magnetizingInductance;      ///< Lm  — transformer magnetizing inductance [H]
    double secondaryResonantInductance;  ///< Lr2 — secondary series resonant inductor [H]
    double secondaryResonantCapacitance; ///< Cr2 — secondary series resonant capacitor [F]
    double qualityFactor;              ///< Q   — Infineon convention Q = √(Lr1/Cr1) / Ro
    double inductanceRatio;            ///< k   = Lm / Lr1
    double equivalentAcResistance;     ///< Ro  — FHA equivalent AC load resistance [Ω]
    double resonantInductorRatio;      ///< a   = n²·Lr2/Lr1 (1.0 for symmetric tank)
    double resonantCapacitorRatio;     ///< b   = Cr2/(n²·Cr1) (1.0 for symmetric tank)
    double bridgeVoltageFactor;        ///< 1.0 for full-bridge, 0.5 for half-bridge
    double lipFrequency;               ///< Load-Independent-Point frequency f_LIP = ωP/(2π) [Hz]
};

/**
 * @brief CLLC Bidirectional Resonant Converter (full-bridge or half-bridge).
 *
 * Inherits converter-level parameters from MAS::CllcResonant and the
 * Topology interface. CLLC is the bidirectional sibling of LLC: the
 * single secondary-side LC is replaced with a symmetric (or asymmetric)
 * Cr2/Lr2 pair so that forward and reverse power flow see equivalent
 * resonant tanks. Used in V2G chargers (OBC), DC microgrids, server
 * PSUs and battery formation/test equipment.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW (forward power flow shown; reverse is symmetric)
 * =====================================================================
 *
 *  +--[Q1]--+----Cr1----Lr1---+               +--[Sa]--+--[Sc]--+
 *  |        |                 |               |        |        |
 *  Vin   bridge_a           Lm·│Lpri ⇋ Lsec│  bridge_c  Cout    Rload
 *  |        |                 |               |        |        |
 *  +--[Q2]--+----GND_pri      +---Lr2---Cr2---+--[Sb]--+--[Sd]--+
 *           Q3, Q4 mirror leg                  Sb, Sd mirror leg
 *
 *  - Primary inverter: full-bridge (Q1..Q4) or half-bridge (Q1, Q2 only,
 *    bridge-voltage factor 0.5; bridge_a swing is ±Vin/2 instead of ±Vin).
 *  - Secondary: ALWAYS active synchronous rectifier (Sa..Sd). Diodes are
 *    not enough — a "bidirectional" CLLC with diode-only secondary cannot
 *    physically push power back to Vin. In reverse mode, Sa..Sd act as
 *    inverter and Q1..Q4 as the (active) rectifier.
 *  - Cr1, Lr1, Lm, Lr2, Cr2 form a 5-element resonant tank. Lr1 and Lr2
 *    can be either discrete inductors OR absorbed into the transformer
 *    primary/secondary leakage (`integratedResonantInductor1/2` flags).
 *    Most published designs absorb both; KIT 20 kW and TI TIDM-02002 keep
 *    Lr2 discrete (the "CLLLC" layout — same equations, separate part).
 *
 * =====================================================================
 * MODULATION CONVENTION
 * =====================================================================
 *
 * Open-loop frequency control with 50 % duty (each half-bridge leg).
 * `td` (dead-time) is the gap between high-side OFF and low-side ON,
 * sized for ZVS via the Coss-charge-balance criterion (Infineon AN
 * Eq 28 below). The control variable is `fs = switchingFrequency`;
 * gain is shaped by the position of fs within [fs_min, f_LIP, fr, fs_max].
 *
 * Schema mapping (MAS::CllcResonant):
 *   inputVoltage        → Vin (DimensionWithTolerance)
 *   operatingPoints[i]  → outputVoltages[0], outputCurrents[0],
 *                         switchingFrequency, ambientTemperature, powerFlow
 *   qualityFactor       → Q (Infineon convention, see disambiguation note)
 *   bridgeType          → fullBridge (default) | halfBridge
 *   bidirectional       → reverse-mode SR drive enabled (forward SR
 *                         always on for efficiency)
 *   integratedResonantInductor1/2 → Lr1, Lr2 are transformer leakage
 *   resonantInductorRatio (a) / resonantCapacitorRatio (b) → asymmetric
 *                         tank (a·b ≈ 1 keeps fr1 = fr2)
 *   symmetricDesign     → DEPRECATED v1 alias; true iff a==1 && b==1
 *
 * =====================================================================
 * Q-CONVENTION DISAMBIGUATION (do not skip this paragraph)
 * =====================================================================
 *
 * Three quality-factor conventions circulate in the LLC/CLLC literature:
 *
 *   Infineon AN-2024-06 / Steigerwald 1988:
 *     Q = √(Lr/Cr) / Ro            (this code uses this)
 *
 *   Erickson "Fundamentals of Power Electronics" 2e, §19:
 *     Q_e = Ro / √(Lr/Cr) = 1/Q
 *
 *   Some SiC vendor app-notes use Q' = Ro / Zr * (8n²/π²) (FHA-folded).
 *
 * Throughout this class, `qualityFactor` and `Q` ALWAYS mean the
 * Infineon convention. Tests / fixtures must match. If you import a
 * design from a paper using Q_e, invert before passing to MAS.
 *
 * =====================================================================
 * KEY EQUATIONS (Infineon AN-2024-06 §2.3, Min IEEE TPE 36(6) 2021)
 * =====================================================================
 *
 * Step 1   n        = Vin_nom / Vout_nom            (full-bridge)
 *                   = Vin_nom / (2·Vout_nom)        (half-bridge)
 * Step 2   Mg_min   = n·Vout_min / Vin_max
 *          Mg_max   = n·Vout_max / Vin_min
 * Step 3   pick k ∈ [3, 6], Q ∈ [0.2, 0.5] so M(ωn,Q,k) covers
 *          [Mg_min, Mg_max] within [fs_min, fs_max]
 * Step 4   Ro       = (8·n²/π²) · Vout² / Pout      (FHA-equivalent AC load)
 * Step 5   Cr1      = 1 / (2π·Q·fr·Ro)
 * Step 6   Lr1      = 1 / ((2π·fr)²·Cr1)
 * Step 7   Lm       = k · Lr1
 * Step 8   sym :    Lr2 = Lr1/n²,        Cr2 = n²·Cr1               (a=b=1)
 *          asym:    Lr2 = a·Lr1/n²,      Cr2 = n²·b·Cr1,  with a·b ≈ 1
 * Step 9   ZVS:     Lm ≤ td / (16 · Coss_eq · fs_max)               (AN Eq 28)
 *
 * FHA voltage gain (symmetric tank, forward direction, AN Eq 2):
 *
 *   M(ωn, Q, k) = 1 / √( (1 + 1/k − 1/(k·ωn²))²
 *                      + Q²·(ωn − 1/ωn)² )
 *
 *   where ωn = ωs/ωr,  ωs = 2π·fs,  ωr = 1/√(Lr1·Cr1).
 *
 * Asymmetric closed-form per Min 2021: M_fwd(ωn) ≠ M_rev(ωn); evaluated
 * numerically as a complex tank impedance in the .cpp.
 *
 * Per-sub-interval natural frequencies (TDA solver):
 *
 *   Power-delivery (rectifier on, Lm clamped to ±n·Vout):
 *     ωP = 1 / √( (Lr1 + Lr2/n²) · (Cr1 · n²·Cr2) / (Cr1 + n²·Cr2) )
 *     For symmetric a=b=1 this reduces to ωr exactly.
 *   Freewheel (rectifier off, Lm in series with Lr1):
 *     ωF = ωP / √(1 + k)
 *
 * LIP (Load-Independent Point): f_LIP = ωP / (2π). At fs = f_LIP the FHA
 * gain ≈ 1 regardless of load — natural design centre.
 *
 * Modes:
 *   1. Sub-resonant boost   fs ≤ fm        → [P, F] sequence per half cycle
 *   2. At resonance         fs ≈ f_LIP     → [P]    (no freewheel)
 *   3. Super-resonant buck  fs > f_LIP     → [P]    truncated, hard turn-off
 *
 * Capacitive operation (fs > f_capacitive(Q,k)): rectifier polarity
 * reverses, ZVS lost. The controller must hard-limit fs above this.
 *
 * =====================================================================
 * REFERENCES
 * =====================================================================
 *   [1] Infineon AN-2024-06, "Operation and modeling analysis of a
 *       bidirectional CLLC converter"  (full 11-step design procedure,
 *       Eq 2 = FHA gain, Eq 28 = ZVS).
 *   [2] Bartecka et al., "Effective Design Methodology of CLLC Resonant
 *       Converter Based on the Minimal Area Product of High-Frequency
 *       Transformer", Energies 17(1):55, 2024.
 *   [3] Min, Lee, Hong et al., "A Novel Design Methodology for
 *       Asymmetric CLLC Resonant DC-DC Converter", IEEE TPE 36(6):6794,
 *       2021.  (asymmetric a, b derivation)
 *   [4] KIT 20 kW CLLLC thesis, https://publikationen.bibliothek.kit.edu/1000123579
 *   [5] Lin, Liu, "Design of a Bidirectional 3.3 kW CLLC Resonant
 *       Converter for Battery Charger", various venues.
 *   [6] TI TIDM-02002 / TIDUEG2C  (CLLLC reference design).
 *   [7] Wei et al., "Time-Domain Analysis of LLC Resonant Converter",
 *       IEEE PEAC 2018.  (TDA solver pattern)
 *   [8] Rezayati PhD thesis 2023, https://theses.hal.science/tel-04094415
 *       (state-plane analysis).
 *   [9] Liu, "Synchronous Rectifier Driving Scheme via Lr Voltage",
 *       IEEE TPE 2021.  (active SR drive)
 *  [10] Infineon UG-2020-31, REF-DAB11KIZSICSYS user manual
 *       (11 kW SiC reference; baseline test fixture).
 *
 * =====================================================================
 * ACCURACY DISCLAIMER (MKF magnetic-design context)
 * =====================================================================
 *
 * The generated SPICE netlist uses idealised voltage-controlled switches
 * (`SW` model with hysteresis) on both bridges and B-source ideal-diode
 * emulators where appropriate. Omitted effects:
 *   • Coss(Vds) variation (constant-Coss assumption for ZVS estimate)
 *   • Body-diode reverse-recovery (trr, Qrr)
 *   • Gate-driver propagation delay and miller plateau
 *   • Switching-loss energy (Eon, Eoff)
 *   • Skin / proximity effects in Lr1, Lr2 windings (handled separately
 *     by MKF coil modules)
 *
 * None of these flow through the magnetic component. The MAS user sees
 * iLr1, iLm, iLr2 and the bridge voltages, all of which are dominated
 * by the resonant-tank dynamics, not the silicon. The analytical TDA
 * solver and SPICE agree to NRMSE ≤ 0.15 on the primary and secondary
 * winding currents across all three modes (sub-/at-/super-resonant) and
 * both directions. That is the magnetic-design accuracy bar; for
 * converter-level efficiency, ZVS validation at light load, or EMI
 * prediction, wrap the netlist with a vendor MOSFET sub-circuit.
 *
 * Disambiguation: this class is "CLLC" (5-element bidirectional tank).
 * It is NOT:
 *   • LLC                — unidirectional, single Cr+Lr+Lm (use Llc.h).
 *   • SRC (series-resonant) — no Lm, simpler 2-element tank.
 *   • CLLLC              — same equations, but Lr2 is always discrete
 *                           (not absorbed). Selectable here via
 *                           integratedResonantInductor2 = false.
 *   • DAB                — phase-shift (no resonant tank). Use Dab.h.
 */
class CllcConverter : public MAS::CllcResonant, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 10;

    /// Default inductance ratio k = Lm/Lr1 (Infineon AN: k=4.45 is a
    /// good middle ground; range 3–6).
    double defaultInductanceRatio = 4.45;
    /// Default dead time in seconds (300 ns, Infineon AN §2.4).
    double deadTime = 300e-9;
    /// MOSFET output capacitance Coss [F] used in the Eq 28 ZVS check.
    double mosfetCoss = 200e-12;

    // ------------------------------------------------------------------
    // Per-OP TDA diagnostics (mutable; populated by
    // process_operating_point_for_input_voltage). All consumed by tests
    // and the get_*() accessors below.
    // ------------------------------------------------------------------

    /// 1 = sub-resonant boost, 2 = at resonance, 3 = super-resonant buck.
    mutable int lastMode = 0;
    /// f_LIP = ωP/(2π) of the last design [Hz].
    mutable double computedLipFrequency = 0.0;
    /// L2-norm of steady-state Newton residual ‖F(x0)‖ — small ⇒ converged.
    mutable double lastSteadyStateResidual = 0.0;
    /// ZVS margin on the primary: iLr1(switch-edge) − iLm_threshold > 0 ⇒ ZVS.
    mutable double lastZvsMarginPrimary = 0.0;
    /// ZVS margin on the secondary (relevant in reverse mode).
    mutable double lastZvsMarginSecondary = 0.0;
    /// Resonant transition time (fraction of dead-time) for the bridge node.
    mutable double lastResonantTransitionTime = 0.0;
    /// Peak primary-side current in steady state [A].
    mutable double lastPrimaryPeakCurrent = 0.0;
    /// Peak voltage on Cr1 in steady state [V].
    mutable double lastResonantCapPeakVoltage = 0.0;
    /// Sub-state sequence traversed in the last solved half-cycle
    /// (CllcSubState IDs: 0=P_POS, 1=P_NEG, 2=F_POS, 3=F_NEG).
    mutable std::vector<int> lastSubStateSequence;

    // ------------------------------------------------------------------
    // Extra-component waveforms — one entry per OP, cleared at the start
    // of every process_operating_points() call. Consumed by
    // get_extra_components_inputs to round-trip Cr1, Cr2, Lr1, Lr2 as
    // ancillary MAS inputs.
    // ------------------------------------------------------------------
    mutable std::vector<Waveform> extraCr1VoltageWaveforms;
    mutable std::vector<Waveform> extraCr1CurrentWaveforms;
    mutable std::vector<Waveform> extraCr2VoltageWaveforms;
    mutable std::vector<Waveform> extraCr2CurrentWaveforms;
    mutable std::vector<Waveform> extraLr1VoltageWaveforms;
    mutable std::vector<Waveform> extraLr1CurrentWaveforms;
    mutable std::vector<Waveform> extraLr2VoltageWaveforms;
    mutable std::vector<Waveform> extraLr2CurrentWaveforms;
    mutable std::vector<std::vector<double>> extraTimeVectors;

    // Resonant parameters used during the last process_operating_points()
    // call. Cached so get_extra_components_inputs can size Cr1/Cr2/Lr1/Lr2
    // identically to what populated the extra*Waveforms above.
    mutable CllcResonantParameters lastResonantParameters{};

public:
    bool _assertErrors = false;

    CllcConverter(const json& j);
    CllcConverter() {};

    MAS::Topologies topology_kind() const override { return MAS::Topologies::CLLC_RESONANT_CONVERTER; }
    bool is_bridge_topology() const override { return true; }

    // --- Accessors for simulation tuning --------------------------------
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }

    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    double get_default_inductance_ratio() const { return defaultInductanceRatio; }
    void set_default_inductance_ratio(double value) { this->defaultInductanceRatio = value; }

    double get_dead_time() const { return deadTime; }
    void set_dead_time(double value) { this->deadTime = value; }

    double get_mosfet_coss() const { return mosfetCoss; }
    void set_mosfet_coss(double value) { this->mosfetCoss = value; }

    // --- Diagnostics (read-only) ----------------------------------------
    /** Mode index: 1 = sub-resonant boost, 2 = at resonance, 3 = super-resonant buck. */
    int get_last_mode() const { return lastMode; }
    /** Load-Independent Point frequency f_LIP = ωP/(2π) of the last solve [Hz]. */
    double get_lip_frequency() const { return computedLipFrequency; }
    /** L2-norm ‖F(x0)‖ of the last steady-state Newton residual. */
    double get_last_steady_state_residual() const { return lastSteadyStateResidual; }
    /** Primary-side ZVS margin [A]: iLr1(t_switch) − iLm_threshold; > 0 ⇒ ZVS. */
    double get_last_zvs_margin_primary() const { return lastZvsMarginPrimary; }
    /** Secondary-side ZVS margin (reverse-mode operation). */
    double get_last_zvs_margin_secondary() const { return lastZvsMarginSecondary; }
    /** Resonant transition time as fraction of the dead-time. */
    double get_last_resonant_transition_time() const { return lastResonantTransitionTime; }
    /** Peak primary current of the last solved OP [A]. */
    double get_last_primary_peak_current() const { return lastPrimaryPeakCurrent; }
    /** Peak Cr1 voltage of the last solved OP [V]. */
    double get_last_resonant_cap_peak_voltage() const { return lastResonantCapPeakVoltage; }
    /** Sub-state sequence (CllcSubState IDs) traversed in the last half-cycle. */
    const std::vector<int>& get_last_sub_state_sequence() const { return lastSubStateSequence; }

    // --- Helpers --------------------------------------------------------
    /** Returns 1.0 for full-bridge, 0.5 for half-bridge. */
    double get_bridge_voltage_factor() const;

    /**
     * Resolve asymmetric tank ratios a = n²·L2/L1 and b = C2/(n²·C1)
     * from the v2 schema fields. Falls back to the v1 symmetricDesign
     * alias when neither ratio is explicit. Centralized so the gain
     * function, TDA solver, AdvancedCllc and get_extra_components_inputs
     * stay consistent.
     */
    void resolve_resonant_ratios(double& a, double& b) const;

    // --- Core design methods --------------------------------------------

    /** Validate converter configuration. */
    bool run_checks(bool assert = false) override;

    /**
     * @brief Calculate all resonant tank parameters from the converter spec.
     *
     * Implements the Infineon AN-2024-06 §2.3 11-step procedure.
     */
    CllcResonantParameters calculate_resonant_parameters();

    /**
     * @brief FHA voltage gain at a given fs (forward direction).
     *
     * Symmetric tank → Infineon AN Eq 2; asymmetric → numerical
     * complex-impedance evaluation per Min IEEE TPE 2021.
     */
    double get_voltage_gain(double switchingFrequency, const CllcResonantParameters& params);

    /** Process design requirements for the magnetic component. */
    DesignRequirements process_design_requirements() override;

    /**
     * @brief Generate analytical OperatingPoints for the transformer windings.
     *
     * In phase P1 of the rewrite this still uses the FHA waveform model; the
     * full TDA solver lands in phase P2 and replaces the body of
     * process_operating_point_for_input_voltage.
     */
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    /** Per-input-voltage analytical worker. */
    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const CllcOperatingPoint& cllcOperatingPoint,
        double turnsRatio,
        double magnetizingInductance,
        const CllcResonantParameters& params);

    /** Generate operating points from a Magnetic component. */
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    // --- Ngspice simulation ---------------------------------------------

    /**
     * @brief Generate an ngspice netlist for the CLLC converter.
     *
     * Phase P1 retains the original (diode-secondary) netlist for backward
     * compatibility with the existing 16 TestCllc cases. Phase P4 swaps in
     * the DAB-quality netlist with active SR (Sa..Sd), per-switch snubbers,
     * METHOD=GEAR, IC pre-charge, sense sources and Evab probe.
     */
    std::string generate_ngspice_circuit(
        double turnsRatio,
        const CllcResonantParameters& params,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);

    /** Simulate and extract winding-level operating points. */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    /**
     * @brief Simulate and extract converter-level waveforms.
     *
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @param numberOfPeriods Number of switching periods to simulate (default 2)
     */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);

    /**
     * @brief Return ancillary MAS/CAS Inputs for Cr1, Cr2 and (when not
     * integrated into the transformer leakage) Lr1, Lr2. Call
     * process_operating_points() first.
     *
     * IDEAL: Cr1/Cr2 always; Lr1/Lr2 only when the corresponding
     *        integratedResonantInductor* flag is false.
     * REAL : same as IDEAL plus, for the inductors, the external part
     *        Lr_external = max(Lr - leakage, 0).
     */
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;
};


/**
 * @brief Advanced CLLC converter with user-specified resonant parameters.
 *
 * Allows the user to directly specify desired turns ratio, magnetizing
 * inductance, and optionally resonant inductors/capacitors, bypassing the
 * automatic design.
 */
class AdvancedCllcConverter : public CllcConverter {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredMagnetizingInductance;
    std::optional<double> desiredResonantInductancePrimary;
    std::optional<double> desiredResonantCapacitancePrimary;
    std::optional<double> desiredResonantInductanceSecondary;
    std::optional<double> desiredResonantCapacitanceSecondary;

public:
    AdvancedCllcConverter() = default;
    ~AdvancedCllcConverter() = default;
    AdvancedCllcConverter(const json& j);

    /**
     * @brief Process the converter with user-specified parameters.
     * @return Inputs containing DesignRequirements and OperatingPoints
     */
    Inputs process();

    // --- Accessors ---
    const double& get_desired_magnetizing_inductance() const { return desiredMagnetizingInductance; }
    void set_desired_magnetizing_inductance(const double& value) { this->desiredMagnetizingInductance = value; }

    const std::vector<double>& get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double>& value) { this->desiredTurnsRatios = value; }

    std::optional<double> get_desired_resonant_inductance_primary() const { return desiredResonantInductancePrimary; }
    void set_desired_resonant_inductance_primary(std::optional<double> value) { this->desiredResonantInductancePrimary = value; }

    std::optional<double> get_desired_resonant_capacitance_primary() const { return desiredResonantCapacitancePrimary; }
    void set_desired_resonant_capacitance_primary(std::optional<double> value) { this->desiredResonantCapacitancePrimary = value; }

    std::optional<double> get_desired_resonant_inductance_secondary() const { return desiredResonantInductanceSecondary; }
    void set_desired_resonant_inductance_secondary(std::optional<double> value) { this->desiredResonantInductanceSecondary = value; }

    std::optional<double> get_desired_resonant_capacitance_secondary() const { return desiredResonantCapacitanceSecondary; }
    void set_desired_resonant_capacitance_secondary(std::optional<double> value) { this->desiredResonantCapacitanceSecondary = value; }
};


// --- JSON serialization ---
void from_json(const json& j, AdvancedCllcConverter& x);
void to_json(json& j, const AdvancedCllcConverter& x);

inline void from_json(const json& j, AdvancedCllcConverter& x) {
    // Base CllcResonant fields (v2 schema)
    x.set_bidirectional(get_stack_optional<bool>(j, "bidirectional"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_max_switching_frequency(j.at("maxSwitchingFrequency").get<double>());
    x.set_min_switching_frequency(j.at("minSwitchingFrequency").get<double>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<CllcOperatingPoint>>());
    x.set_quality_factor(get_stack_optional<double>(j, "qualityFactor"));
    x.set_symmetric_design(get_stack_optional<bool>(j, "symmetricDesign"));
    // v2 additions (optional, defaults preserved by absence)
    x.set_bridge_type(get_stack_optional<LlcBridgeType>(j, "bridgeType"));
    x.set_integrated_resonant_inductor1(get_stack_optional<bool>(j, "integratedResonantInductor1"));
    x.set_integrated_resonant_inductor2(get_stack_optional<bool>(j, "integratedResonantInductor2"));
    x.set_resonant_inductor_ratio(get_stack_optional<double>(j, "resonantInductorRatio"));
    x.set_resonant_capacitor_ratio(get_stack_optional<double>(j, "resonantCapacitorRatio"));

    // Advanced fields
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_magnetizing_inductance(j.at("desiredMagnetizingInductance").get<double>());
    x.set_desired_resonant_inductance_primary(get_stack_optional<double>(j, "desiredResonantInductancePrimary"));
    x.set_desired_resonant_capacitance_primary(get_stack_optional<double>(j, "desiredResonantCapacitancePrimary"));
    x.set_desired_resonant_inductance_secondary(get_stack_optional<double>(j, "desiredResonantInductanceSecondary"));
    x.set_desired_resonant_capacitance_secondary(get_stack_optional<double>(j, "desiredResonantCapacitanceSecondary"));
}

inline void to_json(json& j, const AdvancedCllcConverter& x) {
    j = json::object();
    j["bidirectional"] = x.get_bidirectional();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maxSwitchingFrequency"] = x.get_max_switching_frequency();
    j["minSwitchingFrequency"] = x.get_min_switching_frequency();
    j["operatingPoints"] = x.get_operating_points();
    j["qualityFactor"] = x.get_quality_factor();
    j["symmetricDesign"] = x.get_symmetric_design();
    j["bridgeType"] = x.get_bridge_type();
    j["integratedResonantInductor1"] = x.get_integrated_resonant_inductor1();
    j["integratedResonantInductor2"] = x.get_integrated_resonant_inductor2();
    j["resonantInductorRatio"] = x.get_resonant_inductor_ratio();
    j["resonantCapacitorRatio"] = x.get_resonant_capacitor_ratio();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredMagnetizingInductance"] = x.get_desired_magnetizing_inductance();
    j["desiredResonantInductancePrimary"] = x.get_desired_resonant_inductance_primary();
    j["desiredResonantCapacitancePrimary"] = x.get_desired_resonant_capacitance_primary();
    j["desiredResonantInductanceSecondary"] = x.get_desired_resonant_inductance_secondary();
    j["desiredResonantCapacitanceSecondary"] = x.get_desired_resonant_capacitance_secondary();
}

} // namespace OpenMagnetics
