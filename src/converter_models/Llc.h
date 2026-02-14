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
 * @brief LLC Resonant Converter (Half-Bridge or Full-Bridge)
 *
 * Inherits converter-level parameters from MAS::LlcResonant and
 * the Topology interface, exactly mirroring the TwoSwitchForward pattern.
 *
 * Analytical equations based on Runo Nielsen's Time Domain Approach (TDA)
 * (cf. www.runonielsen.dk/LLC_LCC.pdf, www.runonielsen.dk/llc.pdf)
 * and standard LLC design methodology (TI SLUA119, ST AN2450, Infineon AN 2012-09).
 *
 * Key design flow (Runo Nielsen / TDA):
 * 1. Fictitious output voltage: Vo ≈ ½·Vin_nom (HB) or Vo ≈ Vin_nom (FB)
 * 2. Turns ratio: n = Vo / Vout
 * 3. Rac = (8 * n^2) / pi^2 * Rload (FHA-referred AC load resistance)
 * 4. Characteristic impedance: Zr = Q * Rac
 * 5. Resonant components: Ls = Zr / (2*pi*fr), C = 1 / (2*pi*fr*Zr)
 * 6. Magnetizing inductance: L = Ln * Ls (Ln = inductance ratio)
 * 7. Two resonant frequencies:
 *      w1 = 1/sqrt(Ls*C)    — power delivery (diodes ON)
 *      w0 = 1/sqrt((Ls+L)*C) — freewheeling (diodes OFF)
 * 8. LIP frequency: f1 = w1/(2*pi)
 *
 * Waveform generation uses piecewise time-domain solution:
 *   Phase A (power delivery): sinusoidal at w1 + linear magnetizing ramp
 *   Phase B (freewheeling):   sinusoidal at w0, IL = ILs
 *   Steady-state via bisection on Vc0 (capacitor start voltage)
 */
class Llc : public MAS::LlcResonant, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 5;

    // Computed resonant-tank values (filled by process_design_requirements)
    double computedResonantInductance = 0;   // Ls
    double computedResonantCapacitance = 0;  // C = C1+C2
    double computedInductanceRatio = 5;      // Ln = L / Ls
    double computedDeadTime = 200e-9;        // Default 200 ns dead time

public:
    bool _assertErrors = false;

    Llc(const json& j);
    Llc() {};

    // Simulation tuning 
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    // Computed resonant-tank accessors 
    double get_computed_resonant_inductance() const { return computedResonantInductance; }
    double get_computed_resonant_capacitance() const { return computedResonantCapacitance; }
    double get_computed_inductance_ratio() const { return computedInductanceRatio; }
    double get_computed_dead_time() const { return computedDeadTime; }
    void set_computed_inductance_ratio(double value) { this->computedInductanceRatio = value; }

    // Topology interface 
    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;

    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    // Per-operating-point analytical waveforms (Runo Nielsen TDA)
    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const LlcOperatingPoint& llcOpPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // Helper: bridge voltage factor 
    /** Returns 0.5 for half-bridge, 1.0 for full-bridge */
    double get_bridge_voltage_factor() const;

    /** Returns the resonant frequency (user-provided or geometric mean of fsw range) */
    double get_effective_resonant_frequency() const;

    // SPICE simulation 
    std::string generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);

    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);
};


// 
// AdvancedLlc: user supplies desired turns ratios & inductances directly
// 
class AdvancedLlc : public Llc {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredMagnetizingInductance;
    std::optional<double> desiredResonantInductance;
    std::optional<double> desiredResonantCapacitance;

public:
    AdvancedLlc() = default;
    ~AdvancedLlc() = default;

    AdvancedLlc(const json& j);

    Inputs process();

    const double& get_desired_magnetizing_inductance() const { return desiredMagnetizingInductance; }
    void set_desired_magnetizing_inductance(const double& value) { this->desiredMagnetizingInductance = value; }

    const std::vector<double>& get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double>& value) { this->desiredTurnsRatios = value; }

    std::optional<double> get_desired_resonant_inductance() const { return desiredResonantInductance; }
    void set_desired_resonant_inductance(std::optional<double> value) { this->desiredResonantInductance = value; }

    std::optional<double> get_desired_resonant_capacitance() const { return desiredResonantCapacitance; }
    void set_desired_resonant_capacitance(std::optional<double> value) { this->desiredResonantCapacitance = value; }
};


// 
// JSON serialization
// 
void from_json(const json& j, AdvancedLlc& x);
void to_json(json& j, const AdvancedLlc& x);

inline void from_json(const json& j, AdvancedLlc& x) {
    // LlcResonant base fields
    x.set_bridge_type(get_stack_optional<LlcBridgeType>(j, "bridgeType"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_integrated_resonant_inductor(get_stack_optional<bool>(j, "integratedResonantInductor"));
    x.set_max_switching_frequency(j.at("maxSwitchingFrequency").get<double>());
    x.set_min_switching_frequency(j.at("minSwitchingFrequency").get<double>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<LlcOperatingPoint>>());
    x.set_quality_factor(get_stack_optional<double>(j, "qualityFactor"));
    x.set_resonant_frequency(get_stack_optional<double>(j, "resonantFrequency"));

    // AdvancedLlc extra fields - use get_stack_optional or check if keys exist
    auto desiredTurnsRatios = get_stack_optional<std::vector<double>>(j, "desiredTurnsRatios");
    if (desiredTurnsRatios.has_value()) {
        x.set_desired_turns_ratios(desiredTurnsRatios.value());
    }
    auto desiredMagnetizingInductance = get_stack_optional<double>(j, "desiredMagnetizingInductance");
    if (desiredMagnetizingInductance.has_value()) {
        x.set_desired_magnetizing_inductance(desiredMagnetizingInductance.value());
    }
    x.set_desired_resonant_inductance(get_stack_optional<double>(j, "desiredResonantInductance"));
    x.set_desired_resonant_capacitance(get_stack_optional<double>(j, "desiredResonantCapacitance"));
}

inline void to_json(json& j, const AdvancedLlc& x) {
    j = json::object();
    j["bridgeType"] = x.get_bridge_type();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["integratedResonantInductor"] = x.get_integrated_resonant_inductor();
    j["maxSwitchingFrequency"] = x.get_max_switching_frequency();
    j["minSwitchingFrequency"] = x.get_min_switching_frequency();
    j["operatingPoints"] = x.get_operating_points();
    j["qualityFactor"] = x.get_quality_factor();
    j["resonantFrequency"] = x.get_resonant_frequency();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredMagnetizingInductance"] = x.get_desired_magnetizing_inductance();
    j["desiredResonantInductance"] = x.get_desired_resonant_inductance();
    j["desiredResonantCapacitance"] = x.get_desired_resonant_capacitance();
}

} // namespace OpenMagnetics
