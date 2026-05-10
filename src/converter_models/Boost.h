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
 * @brief Boost (step-up) DC-DC Converter
 *
 * Inherits converter-level parameters from MAS::Boost and the Topology
 * interface. Pure DC-DC step-up converter (PFC is a separate topology
 * implemented in PowerFactorCorrection.h/.cpp).
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW
 * =====================================================================
 *
 *           L                  D
 *    +---/\/\/\---+----|>|----+----+
 *    |            |           |    |
 *   Vin          [Q]         Co   R_load (Vout)
 *    |            |           |    |
 *    +------------+-----------+----+
 *                 |
 *               PWM (D, Fs)
 *
 * Single-switch hard-switched topology. Switch Q ON: inductor charges
 * from Vin (slope = Vin/L). Switch Q OFF: inductor freewheels through
 * D into Co, transferring energy to the output (slope = (Vin-Vo-Vd)/L,
 * which is negative for Vo > Vin - Vd).
 *
 * The implementation models the synchronous variant by setting the
 * diode-voltage-drop Vd = 0; both EVMs at the mid and high power corner
 * (TPS61178EVM-792, LM5122EVM-1PH) are synchronous boost designs.
 *
 * =====================================================================
 * KEY EQUATIONS (CCM, lossless except optional efficiency η)
 * =====================================================================
 *
 * Voltage conversion ratio (CCM, with diode drop Vd and efficiency η):
 *   D = 1 - (Vin · η) / (Vout + Vd)                            [Eq. 1]
 *
 * Inductor on-time and ripple:
 *   t_on   = D / Fs                                            [Eq. 2]
 *   ΔIL_pp = Vin · t_on / L                                    [Eq. 3]
 *
 * Inductor average current (= input current, equals diode-side average
 * current scaled by 1/(1-D)):
 *   IL_avg = Iout · (Vout + Vd) / Vin                          [Eq. 4]
 *
 * Inductor minimum (CCM check):
 *   IL_min = IL_avg − ΔIL_pp / 2                               [Eq. 5]
 *   If IL_min < 0  →  converter is in DCM, equations below replace
 *   the CCM ones (see Boost.cpp lines 64–72).
 *
 * Inductor sizing for a target ripple ΔIL_max at worst-case Vin:
 *   L_min = Vin · (Vout − Vin) / (ΔIL_max · Fs · Vout)         [Eq. 6]
 *
 * (Eq. 6 is the standard boost result D·(1−D)·Vout / (ΔIL·Fs)
 *  rewritten using D = 1 − Vin/Vout. Implemented at Boost.cpp:150.)
 *
 * =====================================================================
 * KEY EQUATIONS (DCM)
 * =====================================================================
 *
 * Triggered when IL_min < 0. Inductor current is a triangular pulse
 * starting and ending at zero each switching period, followed by a
 * dead-time interval where iL = 0:
 *
 *   t_on  = sqrt( 2 · Iout · L · (Vout + Vd − Vin) /
 *                ( Fs · Vin² ) )                               [Eq. 7]
 *   t_off = t_on · ( (Vout + Vd) / (Vout + Vd − Vin) − 1 )     [Eq. 8]
 *   t_dt  = 1/Fs − t_on − t_off                                [Eq. 9]
 *
 * (Implemented at Boost.cpp:65–66.)
 *
 * =====================================================================
 * MAS WAVEFORM CONVENTION
 * =====================================================================
 *
 * Single magnetic excitation (the inductor itself, not a transformer):
 *   Voltage across L:
 *     V(t) = +Vin              during t_on
 *          = +Vin − Vout − Vd  during t_off  (negative)
 *     → bipolar rectangular waveform, peak-to-peak = Vout + Vd
 *   Current through L:
 *     i(t) = piecewise-linear triangular ramp,
 *            DC offset = IL_avg, peak-to-peak = ΔIL_pp
 *     → TRIANGULAR (CCM) or TRIANGULAR_WITH_DEADTIME (DCM)
 *
 * Only one winding is created in process_operating_points_for_input_voltage()
 * (Boost.cpp:60–81); turns ratio is irrelevant for a single-inductor
 * topology and DesignRequirements.turns_ratios is left empty.
 *
 * =====================================================================
 * REFERENCES — INDUSTRY EVMs USED IN GOLDEN-QUALITY TESTS
 * =====================================================================
 *
 * [1] TI TPS61089 + EVM TPS61089EVM-742  — low corner (~25 W)
 *     User Guide: SLVUAM6 (07 Apr 2016).
 *     Vin 2.7–12 V, Vout 4.5–12.6 V, 7 A peak switch current,
 *     Fs 200 kHz–2.2 MHz adjustable, synchronous boost.
 *
 * [2] TI TPS61178 + EVM TPS61178EVM-792  — mid corner (~32 W)
 *     Datasheet rev. E (09 Aug 2019).
 *     Vin 6–12 V → Vout 16 V @ 2 A, 96 % peak efficiency,
 *     20 V / 10 A integrated synchronous boost IC.
 *
 * [3] TI LM5122 + EVM LM5122EVM-1PH      — high corner (~108 W)
 *     Datasheet rev. H (09 Jun 2017).
 *     Vin 9–20 V → Vout 24 V @ 4.5 A, 65 V wide-Vin sync-boost
 *     controller, peak current-mode control.
 *
 * Textbook cross-check:
 *   [4] Erickson & Maksimović, "Fundamentals of Power Electronics",
 *       3rd ed., Chapter 2.3 (boost CCM analysis), Chapter 5
 *       (DCM boundary condition K = 2L·Fs/R).
 */
class Boost : public MAS::Boost, public Topology {
private:
    int numPeriodsToExtract = 5;  // Number of periods to extract from simulation
    int numSteadyStatePeriods = 50;

    // ---- Per-OP analytical diagnostics ----
    // Populated by process_operating_points_for_input_voltage() for the
    // most recent call. Mirrors the Dab.cpp `lastXxx` pattern. Used by
    // golden-quality tests to assert CCM/DCM mode and operating-point
    // numerics without reaching into the OperatingPoint waveform arrays.
    mutable double lastDutyCycle = 0.0;             // D = 1 - Vin·η/(Vo+Vd)
    mutable double lastInductorAverageCurrent = 0;  // IL_avg [A]
    mutable double lastInductorPeakToPeak = 0;      // ΔIL_pp [A]
    mutable double lastPeakInductorCurrent = 0;     // IL_avg + ΔIL_pp/2 (CCM) or ΔIL_pp (DCM)
    mutable bool   lastIsCcm = true;                // false → DCM (IL_min < 0)
    mutable double lastConductionRatio = 1.0;       // (t_on + t_off) / T  (1.0 in CCM)

public:
    bool _assertErrors = false;

    Boost(const json& j);
    Boost() {
    };

    MAS::Topologies topology_kind() const override { return MAS::Topologies::BOOST_CONVERTER; }

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    // ---- Per-OP diagnostic accessors ----
    // Populated by process_operating_points_for_input_voltage().
    double get_last_duty_cycle() const { return lastDutyCycle; }
    double get_last_inductor_average_current() const { return lastInductorAverageCurrent; }
    double get_last_inductor_peak_to_peak() const { return lastInductorPeakToPeak; }
    double get_last_peak_inductor_current() const { return lastPeakInductorCurrent; }
    bool   get_last_is_ccm() const { return lastIsCcm; }
    double get_last_conduction_ratio() const { return lastConductionRatio; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) override;

    OperatingPoint process_operating_points_for_input_voltage(double inputVoltage, const BaseOperatingPoint& outputOperatingPoint, double inductance);
    double calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);
    double resolve_efficiency();

    /**
     * @brief Generate an ngspice circuit for this Boost converter
     * 
     * Uses the calculated design parameters (inductance, duty cycle)
     * to create a SPICE netlist that can be simulated.
     * 
     * @param inductance Inductance in H
     * @param inputVoltageIndex Which input voltage to use (0=nom, 1=min, 2=max)
     * @param operatingPointIndex Which operating point to simulate
     * @return SPICE netlist string
     */
    std::string generate_ngspice_circuit(
        double inductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);
    
    /**
     * @brief Simulate the Boost converter and extract operating points from waveforms
     * 
     * @param inductance Inductance in H
     * @return Vector of ConverterWaveforms extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance);
    
    /**
     * @brief Simulate the Boost converter and extract operating points from waveforms
     * 
     * @param inductance Inductance in H
     * @param numberOfPeriods Number of switching periods to simulate (default 2)
     * @return Vector of ConverterWaveforms extracted from simulation
     */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        double inductance,
        size_t numberOfPeriods = 2);
};

class AdvancedBoost : public Boost {
private:
    double desiredInductance;

protected:
public:
    AdvancedBoost() = default;
    ~AdvancedBoost() = default;

    AdvancedBoost(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

};


void from_json(const json & j, AdvancedBoost & x);
void to_json(json & j, const AdvancedBoost & x);


inline void from_json(const json & j, AdvancedBoost& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<BaseOperatingPoint>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedBoost & x) {
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