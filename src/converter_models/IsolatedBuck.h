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
 * @brief Isolated Buck Converter (also known as "Flybuck" / "Fly-Buck")
 *
 * A buck converter whose inductor has been replaced by a coupled
 * inductor (transformer). The primary winding still acts as the buck
 * filter inductor and produces the (non-isolated) primary output
 * V_pri = D · V_in (ideal). Each secondary winding rectifies the
 * mirror-imaged voltage during the switch OFF interval, giving
 * isolated outputs whose nominal value is set by the turns ratio:
 *     V_sec_k = V_pri / N_k − V_diode
 * Distinct from a true Flyback in that the primary side IS regulated
 * (it's a real synchronous buck loop) — the secondary outputs are
 * loosely cross-regulated and depend on the primary loop driving the
 * primary V_pri to its target. Hence: cheap isolated bias rails.
 *
 *      Vin ─[S1]─┬─ Lpri ─┬── V_pri  (regulated, non-isolated)
 *                │        │            ↓
 *               [S2/D]   Cpri        Rload_pri
 *                │        │
 *               GND      GND
 *                │
 *      ╞═════════╪═══════════╡  K=0.99 coupled
 *                │
 *                Lsec_k ── D_k ── Vout_k  (isolated)
 *                                  │
 *                                 Cout_k − Rload_k
 *
 *  KEY EQUATIONS (CCM, lossless):
 *    (1)  D       = V_pri / V_in                           (ideal)
 *         D       = V_pri / (V_in · η)                     (with η)
 *    (2)  ΔIL_mag = (V_in − V_pri) · t_on / Lmag
 *    (3)  V_sec_k = V_pri / N_k − V_diode                  (nominal)
 *    (4)  N_k     = V_pri / (V_sec_k + V_diode)            (turns ratio)
 *    (5)  IL_pri_avg = I_out_pri + Σ I_out_sec_k / N_k     (reflected)
 *    (6)  IL_pri_pk  = IL_pri_avg + ΔIL_mag / 2
 *
 *  DCM BOUNDARY (CCM ↔ DCM):
 *    Same condition as a regular buck on the primary inductor:
 *    K = 2·Lmag·Fs / R_load_pri,   K_crit = 1 − D.
 *    CCM if K ≥ K_crit, else DCM. In a properly-designed Flybuck the
 *    secondary is sized to keep the primary inductor in CCM under the
 *    expected load range; the secondaries themselves are always in
 *    DCM-like rectifier conduction (only conduct during S OFF).
 *
 *  REFERENCE DESIGNS (verified via vendor product pages):
 *    • Low corner  — TI LM5160 "Flybuck Quick Start" (24 V → ±12 V
 *      bias @ 100 mA, 350 kHz).
 *    • Mid corner  — TI LM5017 Flybuck app note SNVA674A (48 V → 5 V
 *      isolated bias @ 1 A, 200 kHz).
 *    • High corner — TI LM5160-Q1 automotive Flybuck reference
 *      (60 V → 12 V isolated bias @ 1 A, 350 kHz, AEC-Q100).
 *
 *  TEXTBOOK CROSS-CHECK:
 *    [1] Erickson & Maksimović, "Fundamentals of Power Electronics",
 *        3rd ed., Ch. 2.2 (buck CCM analysis), Ch. 5 (DCM boundary).
 *    [2] TI app note SNVA674A "Designing an Isolated Buck (Flybuck)
 *        Converter" (LM5017 / LM5160).
 */
class IsolatedBuck : public MAS::IsolatedBuck, public Topology {
private:
    // Simulation tuning
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 5;

    // Maximum operating duty cycle. Isolated-buck (Flybuck) primary is
    // a regulated synchronous buck whose CCM duty is
    // D = V_pri/(V_in*eta); at low Vin this approaches 1. Real
    // controllers cap D well below 1 to keep the secondary-rectifier
    // OFF interval long enough for energy transfer. Throw on violation
    // (no silent clamp). Mirrors Flyback (04272d7b), forward family
    // (683e731c), Buck (2c9300c2), Boost (96fdb52a).
    std::optional<double> maximumDutyCycle = 0.95;

    // Per-OP diagnostic outputs (mutable — populated by
    // process_operating_point_for_input_voltage()). Allow tests and
    // the Web frontend to inspect what the analytical model decided
    // without re-running it.
    mutable double lastDutyCycle = 0.0;                    // D = V_pri·η/V_in (CCM, asymptotic)
    mutable double lastMagnetizingCurrentRipple = 0.0;     // ΔIL_mag pk-pk [A]
    mutable double lastPrimaryAverageCurrent = 0.0;        // IL_pri_avg [A]
    mutable double lastPrimaryPeakCurrent = 0.0;           // IL_pri_pk = avg + ΔIL_mag/2 [A]
    mutable double lastSecondaryPeakCurrent = 0.0;         // worst-case secondary peak [A]
    mutable bool   lastIsCcm = true;                       // false → primary inductor entered DCM

    mutable std::vector<Waveform> extraLoVoltageWaveforms;
    mutable std::vector<Waveform> extraLoCurrentWaveforms;
    mutable std::vector<double>   extraLoInductances;

public:
    bool _assertErrors = false;

    IsolatedBuck(const json& j);
    IsolatedBuck() {
    };

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    std::optional<double> get_maximum_duty_cycle() const { return maximumDutyCycle; }
    void set_maximum_duty_cycle(std::optional<double> value) { this->maximumDutyCycle = value; }

    // Per-OP analytical diagnostics (read-only; populated by
    // process_operating_point_for_input_voltage()).
    double get_last_duty_cycle() const { return lastDutyCycle; }
    double get_last_magnetizing_current_ripple() const { return lastMagnetizingCurrentRipple; }
    double get_last_primary_average_current() const { return lastPrimaryAverageCurrent; }
    double get_last_primary_peak_current() const { return lastPrimaryPeakCurrent; }
    double get_last_secondary_peak_current() const { return lastSecondaryPeakCurrent; }
    bool   get_last_is_ccm() const { return lastIsCcm; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint process_operating_point_for_input_voltage(double inputVoltage, const IsolatedBuckOperatingPoint& outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance);
    double calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency);

    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    /**
     * @brief Generate an ngspice circuit for this Isolated Buck (Flybuck) converter
     * 
     * @param turnsRatios Turns ratios for secondary windings
     * @param magnetizingInductance Magnetizing inductance in H
     * @param inputVoltageIndex Which input voltage to use (0=nom, 1=min, 2=max)
     * @param operatingPointIndex Which operating point to simulate
     * @return SPICE netlist string
     */
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);
    
    /**
     * @brief Simulate the Isolated Buck converter and extract operating points
     * 
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);
    
    /**
     * @brief Simulate and extract topology-level waveforms for converter validation
     * 
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @param numberOfPeriods Number of switching periods to simulate (default 2)
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);

};

class AdvancedIsolatedBuck : public IsolatedBuck {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;

protected:
public:
    AdvancedIsolatedBuck() = default;
    ~AdvancedIsolatedBuck() = default;

    AdvancedIsolatedBuck(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    std::vector<double> & get_mutable_desired_turns_ratios() { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

};


void from_json(const json & j, AdvancedIsolatedBuck & x);
void to_json(json & j, const AdvancedIsolatedBuck & x);


inline void from_json(const json & j, AdvancedIsolatedBuck& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<IsolatedBuckOperatingPoint>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedIsolatedBuck & x) {
    j = json::object();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumSwitchCurrent"] = x.get_maximum_switch_current();
    j["operatingPoints"] = x.get_operating_points();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredInductance"] = x.get_desired_inductance();
}
} // namespace OpenMagnetics