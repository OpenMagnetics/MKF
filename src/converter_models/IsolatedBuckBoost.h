#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"

using namespace MAS;

namespace OpenMagnetics {


class IsolatedBuckBoost : public MAS::IsolatedBuckBoost, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 20;  // Increased from 5 to ensure steady state

    // Maximum operating duty cycle. Isolated buck-boost (flyback-style)
    // CCM duty is D = Vout/(Vin+Vout); at low Vin this approaches 1.
    // Real controllers cap D well below 1 (typ. 0.5-0.85) so the
    // secondary-rectifier conduction interval (1-D) stays long enough
    // for energy transfer. Throw on violation (no silent clamp).
    // Mirrors Flyback (04272d7b), forward family (683e731c), Buck
    // (2c9300c2), Boost (96fdb52a), IsolatedBuck (703bc80e).
    std::optional<double> maximumDutyCycle = 0.85;

    // Per-OP diagnostic outputs (mutable — populated by
    // processOperatingPointsForInputVoltage()). Same pattern as
    // IsolatedBuck — exposes the analytical model's intermediate
    // results so tests / Web frontend can verify without re-deriving
    // the equations.
    //
    // Topology refresher (η=1, Vd=0):
    //   D       = Vout / (Vin + Vout)              (flyback duty)
    //   ΔIm_pp  = Vin · D / (L · Fs)
    //           = Vin · Vout / ((Vin+Vout) · L · Fs)
    //   Im_avg  = Iout_pri + Σ Iout_sec / n        (sum of reflected loads)
    //   Im_pk   = Im_avg + ΔIm_pp / 2
    //   K       = 2 · L · Fs / R_eff_pri,   K_crit = (1 − D)²
    //   CCM iff K ≥ K_crit (Erickson §5.2 buck-boost criterion).
    mutable double lastDutyCycle = 0.0;
    mutable double lastMagnetizingCurrentRipple = 0.0;
    mutable double lastPrimaryAverageCurrent = 0.0;
    mutable double lastPrimaryPeakCurrent = 0.0;
    mutable double lastSecondaryPeakCurrent = 0.0;
    mutable bool   lastIsCcm = true;

    // ---- Per-OP diagnostic vectors ----
    mutable std::vector<std::string> perOpName;
    mutable std::vector<double>  perOpDutyCycle;
    mutable std::vector<double>  perOpMagnetizingCurrentRipple;
    mutable std::vector<double>  perOpPrimaryAverageCurrent;
    mutable std::vector<double>  perOpPrimaryPeakCurrent;
    mutable std::vector<double>  perOpSecondaryPeakCurrent;
    mutable std::vector<bool>    perOpIsCcm;

    mutable std::vector<Waveform> extraLoVoltageWaveforms;
    mutable std::vector<Waveform> extraLoCurrentWaveforms;
    mutable std::vector<double>   extraLoInductances;

public:
    bool _assertErrors = false;

    IsolatedBuckBoost(const json& j);
    IsolatedBuckBoost() {
    };

    MAS::Topologies topology_kind() const override { return MAS::Topologies::ISOLATED_BUCK_BOOST_CONVERTER; }

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    std::optional<double> get_maximum_duty_cycle() const { return maximumDutyCycle; }
    void set_maximum_duty_cycle(std::optional<double> value) { this->maximumDutyCycle = value; }

    // Per-OP analytical diagnostics (read-only; populated by
    // processOperatingPointsForInputVoltage()).
    double get_last_duty_cycle() const { return lastDutyCycle; }
    double get_last_magnetizing_current_ripple() const { return lastMagnetizingCurrentRipple; }
    double get_last_primary_average_current() const { return lastPrimaryAverageCurrent; }
    double get_last_primary_peak_current() const { return lastPrimaryPeakCurrent; }
    double get_last_secondary_peak_current() const { return lastSecondaryPeakCurrent; }
    bool   get_last_is_ccm() const { return lastIsCcm; }

    // ---- Per-OP vector accessors ----
    const std::vector<std::string>& get_per_op_name() const { return perOpName; }
    const std::vector<double>& get_per_op_duty_cycle() const { return perOpDutyCycle; }
    const std::vector<double>& get_per_op_magnetizing_current_ripple() const { return perOpMagnetizingCurrentRipple; }
    const std::vector<double>& get_per_op_primary_average_current() const { return perOpPrimaryAverageCurrent; }
    const std::vector<double>& get_per_op_primary_peak_current() const { return perOpPrimaryPeakCurrent; }
    const std::vector<double>& get_per_op_secondary_peak_current() const { return perOpSecondaryPeakCurrent; }
    const std::vector<bool>& get_per_op_is_ccm() const { return perOpIsCcm; }

    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint processOperatingPointsForInputVoltage(double inputVoltage, IsolatedBuckBoostOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance);
    double calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency);

    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    /**
     * @brief Generate an ngspice circuit for this Isolated Buck-Boost converter
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
     * @brief Simulate the Isolated Buck-Boost converter and extract operating points
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

class AdvancedIsolatedBuckBoost : public IsolatedBuckBoost {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;

protected:
public:


    AdvancedIsolatedBuckBoost() = default;
    ~AdvancedIsolatedBuckBoost() = default;

    AdvancedIsolatedBuckBoost(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    std::vector<double> & get_mutable_desired_turns_ratios() { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

};


void from_json(const json & j, AdvancedIsolatedBuckBoost & x);
void to_json(json & j, const AdvancedIsolatedBuckBoost & x);


inline void from_json(const json & j, AdvancedIsolatedBuckBoost& x) {
    x.set_current_ripple_ratio(get_stack_optional<double>(j, "currentRippleRatio"));
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_switch_current(get_stack_optional<double>(j, "maximumSwitchCurrent"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<IsolatedBuckBoostOperatingPoint>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json & j, const AdvancedIsolatedBuckBoost & x) {
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