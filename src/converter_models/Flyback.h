#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"


namespace OpenMagnetics {
using namespace MAS;

class FlybackOperatingPoint : public MAS::FlybackOperatingPoint{
public:
    double resolve_switching_frequency(double inputVoltage, double diodeVoltageDrop, std::optional<double> inductance = std::nullopt, std::optional<std::vector<double>> turnsRatios = std::nullopt, double efficiency = 0.85);
    FlybackModes resolve_mode(std::optional<double> currentRippleRatio = std::nullopt);
};


/**
 * @brief Flyback Converter (multi-output capable)
 *
 * Buck-boost-derived isolated converter with energy stored in a coupled
 * inductor (transformer with airgap). Supports multiple secondaries via
 * power-share allocation. Resolved per-OP via FlybackOperatingPoint::
 * resolve_mode() / resolve_switching_frequency() (CCM / DCM / QRM / BMO).
 *
 *  TOPOLOGY (single secondary, dot polarities matter — primary and
 *  secondary are wound with OPPOSITE polarity so energy is stored
 *  during t_on and released during t_off):
 *
 *      Vin ──[Sw]──┬─ Lp •          • Ls ─[D]──┬── Vout
 *                  │   │            │          │
 *                 GND  •────────────•         Cout
 *                                              │
 *                                            GND_sec
 *
 *  MODES:
 *    CCM (Continuous Conduction) — I_pri never reaches 0; smaller ΔIL,
 *      higher device stress at fixed power.
 *    DCM (Discontinuous Conduction) — I_pri = 0 for part of t_off;
 *      simpler control, larger ΔIL.
 *    QRM (Quasi-Resonant) — turn-on at valley of Lp·Cds resonance after
 *      demagnetization → soft-switching, variable Fs.
 *    BMO (Boundary / Critical Mode) — turn-on exactly when IL_pri = 0.
 *
 *  KEY EQUATIONS (single secondary, n = N_pri/N_sec):
 *    (1) CCM:  Vout = (Vin · D) / (n · (1 − D))                [Erickson 6.3]
 *    (2) DCM:  Vout = Vin · D · √(Tsw / (2 · Lp · Iout / Vin))
 *    (3) ΔIL_pri = Vin · D · Tsw / Lp
 *    (4) IL_pri,pk = IL_pri,avg + ΔIL_pri/2  (CCM)
 *                  = ΔIL_pri                  (DCM, IL_min = 0)
 *    (5) Vds_max  = Vin_max + n · (Vout + Vd) + V_leakage_spike
 *
 *  REFERENCE DESIGNS (verified via vendor product pages):
 *    See TestTopologyFlyback.cpp — three EVMs spanning low / mid / high
 *    power: TI **LM5180-Q1EVM** (low ~3 W, primary-side-regulated),
 *    TI **LM5021EVM** (mid, current-mode flyback controller), and
 *    TI **UCC28C44EVM-296** (high, current-mode controller demo).
 *
 *  TEXTBOOK CROSS-CHECK:
 *    [1] Erickson & Maksimović, "Fundamentals of Power Electronics",
 *        3rd ed., Chapter 6.3 (flyback CCM/DCM analysis).
 *    [2] Pressman, "Switching Power Supply Design", 3rd ed., Chapter 4.
 *    [3] TI SLUP254 — "Designing Flyback Converters Using Peak-Current-
 *        Mode Control" (worked design example, 48 V → 12 V CCM).
 *
 *  ACCURACY DISCLAIMER:
 *    The analytical model assumes ideal coupling (k = 1), zero leakage,
 *    instantaneous diode commutation, and constant-frequency operation
 *    (QRM/BMO are emitted with a representative Fs from
 *    resolve_switching_frequency()). Real converters show leakage
 *    spikes at switch turn-off, reverse-recovery on the secondary
 *    diode, and Fs hopping in QRM. Diagnostics fields (`last*`)
 *    capture the last call to
 *    `process_operating_points_for_input_voltage()`.
 */
class Flyback : public MAS::Flyback, public Topology {
private:
    std::optional<double> maximumDrainSourceVoltage = 600;
    std::vector<OpenMagnetics::FlybackOperatingPoint> operatingPoints;
    std::optional<double> maximumDutyCycle = 0.5;
    double efficiency = 1;
    int numPeriodsToExtract = 5;  // Number of periods to extract from simulation
    int numSteadyStatePeriods = 5;  // Number of steady-state cycles to skip

    // ---- Per-OP analytical diagnostics ----
    // Populated by process_operating_points_for_input_voltage() for the
    // most recent call. Used by golden-quality tests to assert mode and
    // operating-point numerics without reaching into OperatingPoint
    // waveform arrays.
    mutable FlybackModes lastMode = FlybackModes::CONTINUOUS_CONDUCTION_MODE;
    mutable double lastDutyCycle = 0.0;             // D
    mutable double lastSwitchingFrequency = 0.0;    // Fs [Hz]
    mutable double lastPrimaryAverageCurrent = 0.0; // IL_pri,avg [A]
    mutable double lastPrimaryPeakToPeak = 0.0;     // ΔIL_pri,pp [A]
    mutable double lastPrimaryPeakCurrent = 0.0;    // IL_pri,pk [A]
    mutable double lastSecondaryPeakCurrent = 0.0;  // IL_sec,pk lumped [A]
    mutable bool   lastIsCcm = true;                // false → DCM/QRM/BMO

    // ---- Per-OP diagnostic vectors (one entry per V_in × OP iteration)
    // populated by process_operating_points so the wizard can render a
    // Diagnostics table with one column per OP. Cleared at the start of
    // each process_operating_points call.
    mutable std::vector<std::string>  perOpName;
    mutable std::vector<FlybackModes> perOpMode;
    mutable std::vector<double>       perOpDutyCycle;
    mutable std::vector<double>       perOpSwitchingFrequency;
    mutable std::vector<double>       perOpPrimaryAverageCurrent;
    mutable std::vector<double>       perOpPrimaryPeakToPeak;
    mutable std::vector<double>       perOpPrimaryPeakCurrent;
    mutable std::vector<double>       perOpSecondaryPeakCurrent;
    mutable std::vector<bool>         perOpIsCcm;

public:
    // ---- Per-OP diagnostic accessors ----
    FlybackModes get_last_mode() const { return lastMode; }
    double get_last_duty_cycle() const { return lastDutyCycle; }
    double get_last_switching_frequency() const { return lastSwitchingFrequency; }
    double get_last_primary_average_current() const { return lastPrimaryAverageCurrent; }
    double get_last_primary_peak_to_peak() const { return lastPrimaryPeakToPeak; }
    double get_last_primary_peak_current() const { return lastPrimaryPeakCurrent; }
    double get_last_secondary_peak_current() const { return lastSecondaryPeakCurrent; }
    bool   get_last_is_ccm() const { return lastIsCcm; }

    // ---- Per-OP vector accessors ----
    const std::vector<std::string>&  get_per_op_name()                      const { return perOpName; }
    const std::vector<FlybackModes>& get_per_op_mode()                      const { return perOpMode; }
    const std::vector<double>&       get_per_op_duty_cycle()                const { return perOpDutyCycle; }
    const std::vector<double>&       get_per_op_switching_frequency()       const { return perOpSwitchingFrequency; }
    const std::vector<double>&       get_per_op_primary_average_current()   const { return perOpPrimaryAverageCurrent; }
    const std::vector<double>&       get_per_op_primary_peak_to_peak()      const { return perOpPrimaryPeakToPeak; }
    const std::vector<double>&       get_per_op_primary_peak_current()      const { return perOpPrimaryPeakCurrent; }
    const std::vector<double>&       get_per_op_secondary_peak_current()    const { return perOpSecondaryPeakCurrent; }
    const std::vector<bool>&         get_per_op_is_ccm()                    const { return perOpIsCcm; }

    MAS::Topology topology_kind() const override { return MAS::Topology::FLYBACK_CONVERTER; }


    const std::vector<OpenMagnetics::FlybackOperatingPoint> & get_operating_points() const { return operatingPoints; }
    std::vector<OpenMagnetics::FlybackOperatingPoint> & get_mutable_operating_points() { return operatingPoints; }
    void set_operating_points(const std::vector<OpenMagnetics::FlybackOperatingPoint> & value) { this->operatingPoints = value; }

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    bool _assertErrors = false;

    Flyback(const json& j);
    Flyback() {
        maximumDrainSourceVoltage = 600;
        maximumDutyCycle = 0.5;
        efficiency = 1;
        // Also set the parent's efficiency via setter (the local 'efficiency' shadows the parent's)
        set_efficiency(1);
    };

    bool run_checks(bool assert = false) override;

    // According to Worked Example (7), pages 135-144 — Designing the Flyback Transformer of Switching Power Supplies A - Z (Second Edition) by Sanjaya Maniktala
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, OpenMagnetics::FlybackOperatingPoint outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance, std::optional<FlybackModes> customMode=std::nullopt, std::optional<double> customDutyCycle=std::nullopt, std::optional<double> customDeadTime=std::nullopt);
    
    /**
     * @brief Generate an ngspice circuit for this flyback converter
     * 
     * Uses the calculated design parameters (turns ratio, inductance, duty cycle)
     * to create a SPICE netlist that can be simulated.
     * 
     * @param turnsRatios Turns ratios for each secondary winding
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
     * @brief Simulate the flyback converter and extract operating points from waveforms
     * 
     * This method extracts winding voltage/current waveforms suitable for magnetic
     * component analysis (core losses, winding losses, etc.)
     * 
     * 1. Generates ngspice circuits for each input voltage / operating point combination
     * 2. Runs ngspice simulations
     * 3. Extracts waveforms and converts them to MAS::OperatingPoint format
     * 
     * @param turnsRatios Turns ratios for each secondary
     * @param magnetizingInductance Magnetizing inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 0);
    
    /**
     * @brief Generate an ngspice circuit using a real Magnetic component model
     * 
     * Uses CircuitSimulatorExporterNgspiceModel::export_magnetic_as_subcircuit to
     * generate a realistic magnetic component with winding resistances, leakage
     * inductance, and AC resistance modeling via ladder networks.
     * 
     * @param magnetic The Magnetic component to use in the simulation
     * @param inputVoltageIndex Which input voltage to use (0=nom, 1=min, 2=max)
     * @param operatingPointIndex Which operating point to simulate
     * @return SPICE netlist string
     */
    std::string generate_ngspice_circuit_with_magnetic(
        const Magnetic& magnetic,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);
    
    /**
     * @brief Simulate the flyback converter using a real Magnetic component
     * 
     * This method uses the actual magnetic component model (with winding resistance,
     * leakage inductance, coupling coefficient) instead of ideal coupled inductors.
     * The magnetic model is generated using CircuitSimulatorExporterNgspiceModel.
     * 
     * @param magnetic The Magnetic component to simulate with
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_with_magnetic_and_extract_operating_points(
        const Magnetic& magnetic);
    
    /**
     * @brief Simulate and extract topology-level waveforms for converter validation
     * 
     * This method extracts all relevant converter signals (input/output voltages,
     * switch node, currents) to validate that the simulation matches expected
     * flyback converter behavior.
     * 
     * @param turnsRatios Turns ratios for each secondary
     * @param magnetizingInductance Magnetizing inductance in H
     * @param numberOfPeriods Number of switching periods to simulate (default 2)
     * @return Vector of OperatingPoints for each operating condition
     */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);
    
    static double get_total_input_power(std::vector<double> outputCurrents, std::vector<double> outputVoltages, double efficiency, double diodeVoltageDrop);
    static double get_total_input_power(double outputCurrent, double outputVoltage, double efficiency, double diodeVoltageDrop);
    static double get_total_input_current(std::vector<double> outputCurrents, double inputVoltage, std::vector<double> outputVoltages, double diodeVoltageDrop);
    static double get_minimum_output_reflected_voltage(double maximumDrainSourceVoltage, double maximumInputVoltage, double safetyMargin=0.85);

};

class AdvancedFlyback : public Flyback {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;
    std::vector<std::vector<double>> desiredDutyCycle;
    std::optional<std::vector<double>> desiredDeadTime;

protected:
public:
    AdvancedFlyback() = default;
    ~AdvancedFlyback() = default;

    AdvancedFlyback(const json& j);

    /**
     * @brief Override the parent Flyback's process_design_requirements()
     *        so that callers using the generic
     *        construct-from-json → process_design_requirements → ...
     *        flow succeed even when neither maximumDutyCycle nor
     *        maximumDrainSourceVoltage is supplied.
     *
     * In the Advanced workflow the design has already been picked by
     * the user — desiredInductance and desiredTurnsRatios fix L and n
     * directly, so the parent's iterative search over
     * (maximumDutyCycle, maximumDrainSourceVoltage) is moot.
     *
     * Mirrors the DesignRequirements-building portion of
     * AdvancedFlyback::process(); see Flyback.cpp for the override.
     * Resolves Issue M1 (Advanced* API contract) for AdvancedFlyback.
     */
    DesignRequirements process_design_requirements() override;

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    const std::vector<std::vector<double>> & get_desired_duty_cycle() const { return desiredDutyCycle; }
    std::vector<std::vector<double>> & get_mutable_desired_duty_cycle() { return desiredDutyCycle; }
    void set_desired_duty_cycle(const std::vector<std::vector<double>> & value) { this->desiredDutyCycle = value; }

    std::optional<std::vector<double>> get_desired_dead_time() const { return desiredDeadTime; }
    void set_desired_dead_time(std::optional<std::vector<double>> value) { this->desiredDeadTime = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    std::vector<double> & get_mutable_desired_turns_ratios() { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

};


void from_json(const json & j, FlybackOperatingPoint & x);
void to_json(json & j, const FlybackOperatingPoint & x);

inline void from_json(const json & j, FlybackOperatingPoint& x) {
    x.set_output_voltages(j.at("outputVoltages").get<std::vector<double>>());
    x.set_output_currents(j.at("outputCurrents").get<std::vector<double>>());
    x.set_switching_frequency(get_stack_optional<double>(j, "switchingFrequency"));
    x.set_mode(get_stack_optional<FlybackModes>(j, "mode"));
    x.set_ambient_temperature(j.at("ambientTemperature").get<double>());
}

inline void to_json(json & j, const FlybackOperatingPoint & x) {
    j = json::object();
    j["outputVoltages"] = x.get_output_voltages();
    j["outputCurrents"] = x.get_output_currents();
    j["switchingFrequency"] = x.get_switching_frequency();
    j["mode"] = x.get_mode();
    j["ambientTemperature"] = x.get_ambient_temperature();
}

void from_json(const json & j, Flyback & x);
void to_json(json & j, const Flyback & x);

inline void from_json(const json & j, Flyback& x) {
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_maximum_drain_source_voltage(get_stack_optional<double>(j, "maximumDrainSourceVoltage"));
    x.set_maximum_duty_cycle(get_stack_optional<double>(j, "maximumDutyCycle"));
    x.set_current_ripple_ratio(j.at("currentRippleRatio").get<double>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<FlybackOperatingPoint>>());
    x.set_efficiency(j.at("efficiency").get<double>());
}

inline void to_json(json & j, const Flyback & x) {
    j = json::object();
    j["inputVoltage"] = x.get_input_voltage();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["maximumDrainSourceVoltage"] = x.get_maximum_drain_source_voltage();
    j["maximumDutyCycle"] = x.get_maximum_duty_cycle();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["operatingPoints"] = x.get_operating_points();
    j["efficiency"] = x.get_efficiency();
}

void from_json(const json & j, AdvancedFlyback & x);
void to_json(json & j, const AdvancedFlyback & x);

inline void from_json(const json & j, AdvancedFlyback& x) {
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
    x.set_desired_dead_time(get_stack_optional<std::vector<double>>(j, "desiredDeadTime"));
    x.set_desired_duty_cycle(j.at("desiredDutyCycle").get<std::vector<std::vector<double>>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<FlybackOperatingPoint>>());
    x.set_efficiency(j.at("efficiency").get<double>());
    x.set_current_ripple_ratio(std::numeric_limits<double>::quiet_NaN());
}

inline void to_json(json & j, const AdvancedFlyback & x) {
    j = json::object();
    j["inputVoltage"] = x.get_input_voltage();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["desiredInductance"] = x.get_desired_inductance();
    j["desiredDutyCycle"] = x.get_desired_duty_cycle();
    j["desiredDeadTime"] = x.get_desired_dead_time();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["operatingPoints"] = x.get_operating_points();
    j["efficiency"] = x.get_efficiency();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
}
} // namespace OpenMagnetics