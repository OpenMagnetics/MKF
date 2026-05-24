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
 * @brief Buck (Step-Down) Converter — synchronous or diode rectified
 *
 * Single-switch hard-switched DC-DC converter. The high-side switch
 * connects Vin to the inductor input during t_on; the freewheeling
 * path (low-side diode in async, low-side MOSFET in synchronous EVMs)
 * carries the inductor current during t_off.
 *
 *      Vin ──[Sw]──┬── L ──┬── Vout
 *                  │       │
 *                 [D]     Cout
 *                  │       │
 *                 GND     GND
 *
 *  KEY EQUATIONS (CCM, lossless):
 *    (1)  D     = Vout / Vin                                (ideal)
 *         D     = (Vout + Vd) / ((Vin + Vd) · η)            (with Vd, η)
 *    (2)  ΔIL   = (Vin − Vout) · t_on / L
 *    (3)  IL_avg= Iout
 *    (4)  IL_pk = IL_avg + ΔIL/2
 *
 *  DCM BOUNDARY (CCM ↔ DCM):
 *    (5)  K = 2·L·Fs / R_load,        K_crit = 1 − D
 *         CCM if K ≥ K_crit, else DCM.
 *  In DCM, the duty cycle is solved from a quadratic so that
 *  area-balance Iout = (1/2)·ΔIL·(t_on + t_off)/T is satisfied;
 *  see process_operating_points_for_input_voltage().
 *
 *  REFERENCE DESIGNS (verified via vendor product pages):
 *    • TPS54202EVM-716    — 8–28 V → 5 V @ 2 A, 500 kHz sync buck
 *      (low corner ~10 W).
 *    • LMR33630ADDAEVM    — 12 V → 5 V @ 3 A, 400 kHz sync buck
 *      (mid corner ~15 W).
 *    • LM5146-Q1-EVM12V   — 15–85 V → 12 V @ 8 A, 400 kHz sync-buck
 *      controller (high corner ~96 W).
 *
 *  TEXTBOOK CROSS-CHECK:
 *    [1] Erickson & Maksimović, "Fundamentals of Power Electronics",
 *        3rd ed., Chapter 2.2 (buck CCM analysis), Chapter 5
 *        (DCM boundary K = 2L·Fs/R).
 */
class Buck : public MAS::Buck, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 50;  // Number of periods to extract from simulation

    // Maximum allowed operating duty cycle. The buck-family CCM duty is
    // D = (Vout + Vd) / ((Vin + Vd) · η); for low-line corners with
    // Vin_min close to Vout this approaches 1. We enforce an explicit
    // upper bound and throw when the required D exceeds it (rather than
    // silently clamping, which would make Vout collapse). Mirrors the
    // pattern in Flyback (see Flyback.h::maximumDutyCycle, commit
    // 04272d7b) and the forward family (commit 683e731c).
    std::optional<double> maximumDutyCycle = 0.95;

    // ---- Per-OP analytical diagnostics ----
    // Populated by process_operating_points_for_input_voltage() for the
    // most recent call. Used by golden-quality tests to assert CCM/DCM
    // mode and operating-point numerics without reaching into the
    // OperatingPoint waveform arrays.
    mutable double lastDutyCycle = 0.0;             // D = Vout·η/(Vin) (asymptotic)
    mutable double lastInductorAverageCurrent = 0;  // IL_avg = Iout (CCM)
    mutable double lastInductorPeakToPeak = 0;      // ΔIL_pp [A]
    mutable double lastPeakInductorCurrent = 0;     // IL_avg + ΔIL_pp/2 (CCM) or ΔIL_pp (DCM)
    mutable bool   lastIsCcm = true;                // false → DCM (IL_min < 0)
    mutable double lastConductionRatio = 1.0;       // (t_on + t_off) / T (1.0 in CCM)

public:
    bool _assertErrors = false;

    Buck(const json& j);
    Buck() {
    };

    MAS::Topologies topology_kind() const override { return MAS::Topologies::BUCK_CONVERTER; }

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    std::optional<double> get_maximum_duty_cycle() const { return maximumDutyCycle; }
    void set_maximum_duty_cycle(std::optional<double> value) { this->maximumDutyCycle = value; }

    // ---- Per-OP diagnostic accessors ----
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

    /**
     * @brief Generate an ngspice circuit for this Buck converter
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
     * @brief Simulate the Buck converter and extract operating points from waveforms
     * 
     * @param inductance Inductance in H
     * @return Vector of OperatingPoints extracted from simulation
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance);
    
    /**
     * @brief Simulate and extract operating points from topology waveforms
     * 
     * Runs simulation and extracts operating points from the resulting waveforms.
     * Similar to simulate_and_extract_operating_points but returns raw OperatingPoint
     * objects for further processing.
     * 
     * @param inductance Inductance in H
     * @param numberOfPeriods Number of switching periods to simulate (default 2)
     * @return Vector of OperatingPoints extracted from simulation waveforms
     */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        double inductance,
        size_t numberOfPeriods = 2);
};

class AdvancedBuck : public Buck {
private:
    double desiredInductance;

protected:
public:
    AdvancedBuck() = default;
    ~AdvancedBuck() = default;

    AdvancedBuck(const json& j);

    Inputs process();

    /**
     * @brief Override the parent Buck's process_design_requirements()
     *        so that downstream pipelines reading the
     *        DesignRequirements (CoreAdviser, MagneticAdviser, …)
     *        see the user's ``desiredInductance`` as the L target.
     *
     * Buck::process_design_requirements() sets only
     * ``magnetizingInductance.minimum`` (the ripple-floor). With that
     * the CoreAdviser shortlists the smallest core satisfying the
     * floor, ignoring the L the caller actually asked for. The
     * Advanced workflow says "I already know my L"; setting nominal
     * from desiredInductance lets MKF design for it.
     */
    DesignRequirements process_design_requirements() override;

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

};


void from_json(const json & j, AdvancedBuck & x);
void to_json(json & j, const AdvancedBuck & x);


inline void from_json(const json & j, AdvancedBuck& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<BaseOperatingPoint>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedBuck & x) {
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