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
 * @brief Resonant tank parameters for the CLLC converter.
 *
 * Contains all the calculated values needed to define the resonant tank
 * of a bidirectional CLLC converter. These are computed from the converter
 * specifications using FHA (First Harmonic Approximation) analysis as
 * described in Infineon AN-2024-06 and Bartecka et al., Energies 2024.
 */
struct CllcResonantParameters {
    double turnsRatio;                ///< Transformer turns ratio n = Np/Ns
    double resonantFrequency;         ///< Natural resonant frequency fr [Hz]
    double primaryResonantInductance; ///< L1 - primary series resonant inductor [H]
    double primaryResonantCapacitance;///< C1 - primary series resonant capacitor [F]
    double magnetizingInductance;     ///< Lm - transformer magnetizing inductance [H]
    double secondaryResonantInductance;  ///< L2 - secondary series resonant inductor [H]
    double secondaryResonantCapacitance; ///< C2 - secondary series resonant capacitor [F]
    double qualityFactor;             ///< Q - quality factor of resonant tank
    double inductanceRatio;           ///< k = Lm / L1
    double equivalentAcResistance;    ///< Ro - FHA equivalent AC load resistance [Ohm]
    double resonantInductorRatio;     ///< a = n^2*L2/L1 (1.0 for symmetric)
    double resonantCapacitorRatio;    ///< b = C2/(n^2*C1) (1.0 for symmetric)
};

/**
 * @brief CLLC Bidirectional Resonant Converter model.
 *
 * Implements the design equations and waveform generation for a full-bridge
 * CLLC bidirectional resonant converter. The resonant tank consists of:
 *   Primary:   C1 -- L1 -- [Transformer with Lm] -- L2 -- C2   :Secondary
 *
 * Operating modes:
 *   - fs < fr : Boost mode (voltage gain > 1 at resonance point)
 *   - fs = fr : Resonant mode (highest efficiency)
 *   - fs > fr : Buck mode (voltage gain < 1)
 *
 * Supports both forward (primary→secondary) and reverse (secondary→primary) power flow.
 *
 * Design methodology follows:
 *   [1] Infineon AN-2024-06: "Operation and modeling analysis of a bidirectional CLLC converter"
 *   [2] Bartecka et al., "Effective Design Methodology of CLLC Resonant Converter Based on
 *       the Minimal Area Product of High-Frequency Transformer", Energies 2024, 17, 55.
 */
class CllcConverter : public MAS::CllcResonant, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 10;

    /// Default inductance ratio k = Lm/L1 (from Infineon AN Table: k=4.45 is a good middle ground)
    double defaultInductanceRatio = 4.45;
    /// Default dead time in seconds (300ns from Infineon AN Section 2.4)
    double deadTime = 300e-9;

public:
    bool _assertErrors = false;

    CllcConverter(const json& j);
    CllcConverter() {};

    // --- Accessors for simulation parameters ---
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { this->numPeriodsToExtract = value; }

    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { this->numSteadyStatePeriods = value; }

    double get_default_inductance_ratio() const { return defaultInductanceRatio; }
    void set_default_inductance_ratio(double value) { this->defaultInductanceRatio = value; }

    double get_dead_time() const { return deadTime; }
    void set_dead_time(double value) { this->deadTime = value; }

    // --- Core design methods ---

    /**
     * @brief Check converter configuration validity.
     */
    bool run_checks(bool assert = false) override;

    /**
     * @brief Calculate all resonant tank parameters from the converter specification.
     *
     * Implements the 11-step design procedure from Infineon AN-2024-06 Section 2.3:
     *   Step 1: n = Vin_nom / Vout_nom
     *   Step 2: Mg_min, Mg_max
     *   Step 3: Choose k, Q (defaults or from user)
     *   Step 4: Ro = 8n²/(π²) * Vout²/Pout
     *   Step 5: C1 = 1/(2π*Q*fr*Ro)
     *   Step 6: L1 = 1/((2πfr)²*C1)
     *   Step 7: Lm = k * L1
     *   Step 8-9: L2, C2 (from a, b ratios; symmetric: a=b=1)
     *
     * @return CllcResonantParameters with all computed values
     */
    CllcResonantParameters calculate_resonant_parameters();

    /**
     * @brief Compute FHA voltage gain at a given switching frequency.
     *
     * Uses the transfer function from Infineon AN Equation 2 (forward mode):
     *   H(jωs) = nVout/Vin = f(ωs, ωr, Q, k, a, b)
     *
     * @param switchingFrequency Switching frequency fs [Hz]
     * @param params Resonant tank parameters
     * @return Voltage gain |H(jω)|
     */
    double get_voltage_gain(double switchingFrequency, const CllcResonantParameters& params);

    /**
     * @brief Process design requirements from converter specification.
     *
     * Calculates turns ratio and magnetizing inductance suitable for MAS
     * DesignRequirements. The magnetizing inductance is set to the computed Lm,
     * and the turns ratio is the transformer ratio n.
     *
     * @return DesignRequirements for the magnetic component
     */
    DesignRequirements process_design_requirements() override;

    /**
     * @brief Generate analytical OperatingPoints for the transformer windings.
     *
     * Computes voltage and current waveforms for primary and secondary windings
     * using FHA analysis. The waveforms are:
     *   Primary voltage:  bipolar rectangular ±Vin with dead time
     *   Primary current:  sinusoidal resonant + triangular magnetizing
     *   Secondary voltage: bipolar rectangular ±Vout*n (referred)
     *   Secondary current: sinusoidal (Ip - Im) scaled by turns ratio
     *
     * @param turnsRatios Vector of turns ratios [n] (primary-to-secondary)
     * @param magnetizingInductance Lm in Henries
     * @return Vector of OperatingPoint for each (inputVoltage, operatingPoint) combination
     */
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    /**
     * @brief Generate analytical operating point for a specific input voltage.
     */
    OperatingPoint process_operating_points_for_input_voltage(
        double inputVoltage,
        const CllcOperatingPoint& cllcOperatingPoint,
        double turnsRatio,
        double magnetizingInductance,
        const CllcResonantParameters& params);

    /**
     * @brief Generate operating points from a Magnetic component.
     */
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    // --- Ngspice simulation methods ---

    /**
     * @brief Generate an ngspice netlist for the CLLC converter.
     *
     * Creates a full-bridge CLLC circuit with:
     *   - Primary full-bridge (4 voltage-controlled switches, complementary pairs)
     *   - Series resonant components C1, L1 on primary
     *   - Coupled inductors for transformer (Lpri with Lm, Lsec)
     *   - Series resonant components L2, C2 on secondary
     *   - Secondary full-bridge (diode rectifier)
     *   - Output filter capacitor + load resistor
     *
     * @param turnsRatio Transformer turns ratio n
     * @param params Resonant tank parameters
     * @param inputVoltageIndex Index into collected input voltages
     * @param operatingPointIndex Index into operating points
     * @return Netlist string for ngspice
     */
    std::string generate_ngspice_circuit(
        double turnsRatio,
        const CllcResonantParameters& params,
        size_t inputVoltageIndex = 0,
        size_t operatingPointIndex = 0);

    /**
     * @brief Simulate and extract winding-level operating points.
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    /**
     * @brief Simulate and extract converter-level waveforms.
     * 
     * @param turnsRatios Turns ratios for each winding
     * @param magnetizingInductance Magnetizing inductance in H
     * @param numberOfPeriods Number of switching periods to simulate (default 2)
     * @return Vector of ConverterWaveforms extracted from simulation
     */
    std::vector<ConverterWaveforms> simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods = 2);
};


/**
 * @brief Advanced CLLC converter with user-specified resonant parameters.
 *
 * Allows the user to directly specify desired turns ratio, magnetizing inductance,
 * and optionally resonant inductors/capacitors, bypassing the automatic design.
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
    // Base CllcResonant fields
    x.set_bidirectional(get_stack_optional<bool>(j, "bidirectional"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_max_switching_frequency(j.at("maxSwitchingFrequency").get<double>());
    x.set_min_switching_frequency(j.at("minSwitchingFrequency").get<double>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<CllcOperatingPoint>>());
    x.set_quality_factor(get_stack_optional<double>(j, "qualityFactor"));
    x.set_symmetric_design(get_stack_optional<bool>(j, "symmetricDesign"));

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
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredMagnetizingInductance"] = x.get_desired_magnetizing_inductance();
    j["desiredResonantInductancePrimary"] = x.get_desired_resonant_inductance_primary();
    j["desiredResonantCapacitancePrimary"] = x.get_desired_resonant_capacitance_primary();
    j["desiredResonantInductanceSecondary"] = x.get_desired_resonant_inductance_secondary();
    j["desiredResonantCapacitanceSecondary"] = x.get_desired_resonant_capacitance_secondary();
}

} // namespace OpenMagnetics
