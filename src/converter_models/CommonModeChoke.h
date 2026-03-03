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
 * @brief Structure holding CMC simulation waveforms for analysis
 */
struct CmcSimulationWaveforms {
    std::vector<double> time;
    double frequency;

    // Input signals
    std::vector<double> inputVoltage;              // Common mode noise source

    // CMC line currents
    std::vector<std::vector<double>> windingCurrents;  // Current through each winding

    // LISN output (measured noise)
    std::vector<double> lisnVoltage;              // Voltage at LISN measurement point

    // Metadata
    std::string operatingPointName;
    double commonModeAttenuation;   // Attenuation in dB at test frequency
    double commonModeImpedance;     // Impedance in Ohms at test frequency
    double theoreticalImpedance;    // Theoretical Z = 2*pi*f*L in Ohms
};

/**
 * @brief Common Mode Choke (CMC) for EMI filtering
 *
 * Physical principle:
 *   - CM currents (same phase on all lines) produce ADDITIVE flux → high Z
 *   - DM currents (opposite phase) produce CANCELLING flux → low Z
 *
 * Specification modes (all resolved to impedance points inside the constructor):
 *   1. minimumImpedance[]    : direct (frequency, impedance) pairs
 *   2. targetInsertionLoss[] : (frequency, insertionLoss[dB]) pairs
 *                              Z_cm = Z_line·(10^(IL/20)−1)
 *   3. Estimate from noise   : parasiticCap_pF + dvdt_V_ns + safetyMargin_dB
 *                              → CM current estimate → required Z at 150 kHz
 *
 * Output:
 *   designRequirements.magnetizingInductance = L_cm = max(Z/(2πf))
 *   operatingPoints: sinusoidal CM excitation at dominant frequency
 *
 * Supported configurations: 2 windings (single-phase), 3 (three-phase), 4 (three-phase+neutral)
 *
 * Also provides ngspice simulation methods for CISPR 16 LISN testing and
 * realistic line + switching noise simulation.
 */
class CommonModeChoke : public Topology {
private:
    int numPeriodsToExtract   = 5;
    int numSteadyStatePeriods = 3;

    // ── Parsed fields ──────────────────────────────────────────────
    DimensionWithTolerance operatingVoltage;
    double operatingCurrent         = 1.0;
    double lineFrequency            = 50.0;
    double lineImpedance            = 50.0;
    double ambientTemperature       = 25.0;
    double maximumDcResistance      = 0.0;
    double maximumLeakageInductance = 0.0;
    int    numberOfWindings         = 2;

    // Noise-estimation fields (specMode = "Estimate from noise")
    double parasiticCap_pF  = 0.0;
    double dvdt_V_ns        = 0.0;
    double safetyMargin_dB  = 6.0;

    struct ImpedancePoint     { double frequency; double impedance; };
    struct InsertionLossPoint { double frequency; double insertionLoss; };

    std::vector<ImpedancePoint>     minimumImpedance;
    std::vector<InsertionLossPoint> targetInsertionLoss;

    // ── Computed ───────────────────────────────────────────────────
    double computedInductance = 0.0;
    double dominantFrequency  = 0.0;
    double dominantImpedance  = 0.0;

    // ── Internal helpers ───────────────────────────────────────────
    void computeFromNoiseParams();
    void resolveInductance();

public:
    bool _assertErrors = false;

    CommonModeChoke(const json& j);
    CommonModeChoke() {}

    // ── Accessors ──────────────────────────────────────────────────
    int    get_num_periods_to_extract()   const { return numPeriodsToExtract; }
    void   set_num_periods_to_extract(int v)    { numPeriodsToExtract = v; }
    int    get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void   set_num_steady_state_periods(int v)  { numSteadyStatePeriods = v; }

    double get_computed_inductance() const { return computedInductance; }
    double get_dominant_frequency()  const { return dominantFrequency; }
    double get_dominant_impedance()  const { return dominantImpedance; }
    int    get_number_of_windings()  const { return numberOfWindings; }
    double get_operating_current()   const { return operatingCurrent; }
    double get_line_frequency()      const { return lineFrequency; }
    double get_line_impedance()      const { return lineImpedance; }
    double get_ambient_temperature() const { return ambientTemperature; }

    void   set_operating_current(double v)  { operatingCurrent = v; }
    void   set_line_frequency(double v)     { lineFrequency = v; }
    void   set_line_impedance(double v)     { lineImpedance = v; }
    void   set_ambient_temperature(double v){ ambientTemperature = v; }
    void   set_number_of_windings(int v)    { numberOfWindings = v; }

    const DimensionWithTolerance& get_operating_voltage() const { return operatingVoltage; }
    void set_operating_voltage(const DimensionWithTolerance& v) { operatingVoltage = v; }

    // ── Topology interface ─────────────────────────────────────────
    bool run_checks(bool assert_errors = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    // ── Static helpers (public for unit tests) ─────────────────────
    static double impedanceToInductance(double impedance_Ohms, double frequency_Hz);
    static double inductanceToImpedance(double inductance_H,   double frequency_Hz);
    static double insertionLossToImpedance(double insertionLoss_dB, double lineImpedance_Ohms);
    static double noiseParamsToImpedance(double parasiticCap_pF,
                                         double dvdt_V_ns,
                                         double lineImpedance_Ohms,
                                         double safetyMargin_dB,
                                         double testFrequency_Hz = 150e3);
    static std::vector<std::string> windingNames(int numWindings);

    // ── Ngspice simulation methods ────────────────────────────────
    std::string generate_ngspice_circuit(double inductance, double frequency = 150000);
    std::string generate_realistic_cmc_circuit(
        double inductance,
        double parasiticCap_pF,
        double dvdt_V_ns,
        double lineFrequency_Hz = 50.0);
    std::vector<CmcSimulationWaveforms> simulate_and_extract_waveforms(
        double inductance,
        const std::vector<double>& frequencies);
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance);
    std::vector<OperatingPoint> simulate_realistic_cmc(
        double inductance,
        double parasiticCap_pF,
        double dvdt_V_ns);
};

// Backward-compatible alias
using Cmc = CommonModeChoke;


// ─── AdvancedCommonModeChoke — "I know the design I want" variant ─────
class AdvancedCommonModeChoke : public CommonModeChoke {
private:
    double desiredInductance = 0.0;

public:
    AdvancedCommonModeChoke() = default;
    ~AdvancedCommonModeChoke() = default;
    AdvancedCommonModeChoke(const json& j);

    Inputs process();

    const double& get_desired_inductance() const  { return desiredInductance; }
    double& get_mutable_desired_inductance()       { return desiredInductance; }
    void set_desired_inductance(const double& v)   { desiredInductance = v; }
};

// Backward-compatible alias
using AdvancedCmc = AdvancedCommonModeChoke;

void from_json(const json& j, AdvancedCommonModeChoke& x);
void to_json(json& j, const AdvancedCommonModeChoke& x);

inline void from_json(const json& j, AdvancedCommonModeChoke& x) {
    CommonModeChoke base(j);
    static_cast<CommonModeChoke&>(x) = base;
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
}

inline void to_json(json& j, const AdvancedCommonModeChoke& x) {
    j = json::object();
    j["desiredInductance"] = x.get_desired_inductance();
}

} // namespace OpenMagnetics
