#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief Push-Pull Converter (center-tapped primary, two low-side switches)
 *
 * Transformer-isolated forward-class converter with two primary half-windings
 * sharing a center tap tied to Vin and two low-side switches S1/S2 driven
 * 180° out of phase. Each switch conducts at most half the period (D ≤ 0.5)
 * to avoid simultaneous conduction (which would short the primary). Output
 * is rectified by a center-tapped secondary feeding an LC filter.
 *
 *  TOPOLOGY (single secondary, dot polarities matter):
 *
 *               Vin (center tap)
 *                    │
 *           ┌────────┴────────┐
 *           │                 │
 *        Lpri_top •         • Lpri_bot
 *           │                 │
 *          [S1]              [S2]            (180° out of phase, D≤0.5)
 *           │                 │
 *          GND               GND
 *
 *           Lsec_top •───[D1]─┐
 *                              ├──[Lout]──┬── Vout
 *           Lsec_bot •───[D2]─┘           │
 *                    │                  Cout─── Rload
 *                   GND                   │
 *                                        GND
 *
 *  KEY EQUATIONS (η=1, Vd=0, single secondary, N = Np:Ns):
 *    (1) Vout = (Vin / N) · 2·D                                  [Pressman 5.4]
 *        ⇒ D_per_switch  = (Vout · N) / (2 · Vin)        (≤ 0.5)
 *    (2) On-time per switch:  t1 = D · Tsw = (Vout · N) / (2 · Vin · Fs)
 *    (3) Magnetizing peak:    ΔImag = Vin · t1 / Lm
 *    (4) Primary current during conduction (reflected secondary):
 *        IL_pri,avg(on) = Iout / N          (flat-top approximation)
 *        IL_pri,pk      = Iout/N + ΔImag/2
 *    (5) CCM iff (Iout · (1 − ripple/2)) / N > ΔImag/2
 *    (6) Vds_max ≈ 2 · Vin (worst case, ideal coupling)
 *
 *  REFERENCE DESIGNS (verified via vendor product pages):
 *    See TestTopologyPushPull.cpp — three TI push-pull driver designs:
 *    **SN6501** (low ~1.75 W, 5 V, 410 kHz internal osc),
 *    **SN6505B** (mid ~5 W, 5 V, 420 kHz, 1 A peak), and
 *    **SN6507** (high, 24 V wide-Vin, programmable 100 kHz–2 MHz,
 *    integrated 0.5-A NMOS push-pull driver).
 *
 *  TEXTBOOK CROSS-CHECK:
 *    [1] Pressman, "Switching Power Supply Design", 3rd ed., Ch. 5
 *        (push-pull center-tapped topology, D ≤ 0.5 constraint, flux walk).
 *    [2] Erickson & Maksimović, "Fundamentals of Power Electronics",
 *        3rd ed., §6.2 (transformer-isolated forward family).
 *    [3] TI SLLA284 — "Isolation in AMC1xxx Devices Using Push-Pull
 *        Transformer Drivers" (SN6501/SN6505 application note).
 *
 *  ACCURACY DISCLAIMER:
 *    The analytical model assumes ideal coupling (k = 1), zero leakage,
 *    instantaneous diode commutation, and constant-frequency operation.
 *    Real push-pulls show flux-walk asymmetry (drives need DC blocking or
 *    current-mode control), leakage spikes at switch turn-off, and
 *    rectifier reverse-recovery. Diagnostics fields (`last*`) capture the
 *    last call to `process_operating_points_for_input_voltage()`.
 */
class PushPull : public MAS::PushPull, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 500;

    mutable std::vector<Waveform> extraLoVoltageWaveforms;
    mutable std::vector<Waveform> extraLoCurrentWaveforms;
    mutable double                extraLoInductance = 0.0;

    // ---- Per-OP analytical diagnostics ----
    // Populated by process_operating_points_for_input_voltage() for the
    // most recent call. Used by golden-quality tests to assert mode and
    // operating-point numerics without reaching into OperatingPoint
    // waveform arrays.
    mutable double lastDutyCycle = 0.0;             // D per switch, 0..0.5
    mutable double lastSwitchingFrequency = 0.0;    // Fs [Hz]
    mutable double lastPrimaryAverageCurrent = 0.0; // Midpoint of model ramp = (1-r/2)·Iout/N [A]
    mutable double lastPrimaryPeakCurrent = 0.0;    // IL_pri,pk = (1-r/2)·Iout/N + ΔImag/2 [A]
    mutable double lastMagnetizingPeakCurrent = 0.0;// ΔImag = Vin·t1/Lm [A]
    mutable bool   lastIsCcm = true;                // true → CCM, false → DCM

public:
    // ---- Per-OP diagnostic accessors ----
    double get_last_duty_cycle() const { return lastDutyCycle; }
    double get_last_switching_frequency() const { return lastSwitchingFrequency; }
    double get_last_primary_average_current() const { return lastPrimaryAverageCurrent; }
    double get_last_primary_peak_current() const { return lastPrimaryPeakCurrent; }
    double get_last_magnetizing_peak_current() const { return lastMagnetizingPeakCurrent; }
    bool   get_last_is_ccm() const { return lastIsCcm; }

    MAS::Topologies topology_kind() const override { return MAS::Topologies::PUSH_PULL_CONVERTER; }
    bool is_bridge_topology() const override { return true; }

    bool _assertErrors = false;

    PushPull(const json& j);
    PushPull() {
    };
    
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    bool run_checks(bool assert = false) override;

    Inputs process();
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, const PushPullOperatingPoint& outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance, double outputInductance);
    double get_output_inductance(double mainSecondaryTurnsRatio);
    double get_maximum_duty_cycle();

    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    /**
     * @brief Generate an ngspice circuit for this Push-Pull converter
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
     * @brief Simulate the Push-Pull converter and extract operating points
     * 
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);
    
    /**
     * @brief Simulate and extract operating points from topology waveforms
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

class AdvancedPushPull : public PushPull {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;
    std::optional<double> desiredOutputInductance;

protected:
public:
    bool _assertErrors = false;

    AdvancedPushPull() = default;
    ~AdvancedPushPull() = default;

    AdvancedPushPull(const json& j);

    Inputs process();

    // Issue M1: build DesignRequirements directly from desired* fields,
    // bypassing the parent's PDR which auto-sizes turns ratios/inductance.
    DesignRequirements process_design_requirements() override;

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    std::optional<double> get_desired_output_inductance() const { return desiredOutputInductance; }
    void set_desired_output_inductance(std::optional<double> value) { this->desiredOutputInductance = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

};


void from_json(const json & j, AdvancedPushPull & x);
void to_json(json & j, const AdvancedPushPull & x);

inline void from_json(const json & j, AdvancedPushPull& x) {
    x.set_current_ripple_ratio(j.at("currentRippleRatio").get<double>());
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_duty_cycle(get_stack_optional<double>(j, "dutyCycle"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_drain_source_voltage(get_stack_optional<double>(j, "maximumDrainSourceVoltage"));
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<PushPullOperatingPoint>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
    x.set_desired_output_inductance(get_stack_optional<double>(j, "desiredOutputInductance"));
}

inline void to_json(json & j, const AdvancedPushPull & x) {
    j = json::object();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["dutyCycle"] = x.get_duty_cycle();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumDrainSourceVoltage"] = x.get_maximum_drain_source_voltage();
    j["maximumSwitchCurrent"] = x.get_maximum_switch_current();
    j["operatingPoints"] = x.get_operating_points();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredInductance"] = x.get_desired_inductance();
    j["desiredOutputInductance"] = x.get_desired_output_inductance();
}
} // namespace OpenMagnetics