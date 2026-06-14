#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "processors/NgspiceRunner.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief DMC configuration alias.
 *
 * The MAS schema defines `Configuration` (singlePhase, singlePhaseBalanced,
 * threePhase, threePhaseWithNeutral) — that's the single source of truth.
 * `DmcConfiguration` is a backward-compatible alias so older MKF call sites
 * that use the long name keep compiling.
 */
using DmcConfiguration = MAS::Configuration;

/**
 * @brief Result of attenuation verification via ngspice simulation.
 */
struct DmcAttenuationResult {
    double frequency;           ///< Test frequency in Hz
    double requiredAttenuation; ///< Required attenuation in dB
    double measuredAttenuation; ///< Measured attenuation in dB (from ngspice)
    double theoreticalAttenuation; ///< Theoretical attenuation from LC filter model
    bool passed;                ///< Whether the requirement was met
    std::string message;        ///< Human-readable result message
};

/**
 * @brief Structure holding DMC simulation waveforms for analysis
 */
struct DmcSimulationWaveforms {
    std::vector<double> time;
    double frequency;

    // Input/Output signals
    std::vector<double> inputVoltage;    // Noise source voltage
    std::vector<double> outputVoltage;   // Filtered output voltage
    std::vector<double> inductorCurrent; // Current through DMC

    // Metadata
    std::string operatingPointName;
    double dmAttenuation;  // Attenuation in dB at test frequency
};

/**
 * @brief Differential Mode Choke (DMC) converter model for EMI filter applications.
 *
 * MAS::DifferentialModeChoke owns the spec fields. MKF inherits and adds the
 * MKF-only engineering layer (LC filter sizing, ngspice circuit generation,
 * design verification).
 */
class DifferentialModeChoke : public MAS::DifferentialModeChoke {
public:
    bool _assertErrors = false;

    DifferentialModeChoke() = default;
    DifferentialModeChoke(const json& j);

    /**
     * @brief Process the DMC specification into a complete Inputs structure.
     */
    Inputs process();

    /**
     * @brief Generate design requirements for DMC optimization.
     */
    DesignRequirements process_design_requirements();

    /**
     * @brief Generate operating points representing DMC excitation.
     */
    std::vector<OperatingPoint> process_operating_points();

    /**
     * @brief Resolve the peak winding current.
     *
     * Returns the explicitly provided peakCurrent when present. Otherwise the
     * DMC was sized in "help me with the design" mode (impedance + operating
     * current, no converter peak current), so derive the peak from the operating
     * current using the supplied, explicit ripple/margin fraction:
     *     peak = operatingCurrent * (1 + rippleFraction)
     * Throws [MISSING_DATA] only when neither a peakCurrent nor a positive
     * operatingCurrent is available — i.e. there is genuinely no current info.
     */
    double resolve_peak_current(double rippleFraction) const;

    /**
     * @brief Get the number of windings based on configuration.
     */
    int get_number_of_windings() const {
        auto cfg = MAS::DifferentialModeChoke::get_configuration().value_or(MAS::Configuration::SINGLE_PHASE);
        switch (cfg) {
            case MAS::Configuration::SINGLE_PHASE: return 1;
            case MAS::Configuration::SINGLE_PHASE_BALANCED: return 2;
            case MAS::Configuration::THREE_PHASE: return 3;
            case MAS::Configuration::THREE_PHASE_WITH_NEUTRAL: return 4;
        }
        return 1;
    }

    /**
     * @brief Convenience accessor — returns configuration with explicit default.
     * MAS getter returns std::optional<Configuration>; legacy MKF callers
     * expect a plain enum value.
     */
    MAS::Configuration get_configuration_or_default() const {
        return MAS::DifferentialModeChoke::get_configuration().value_or(MAS::Configuration::SINGLE_PHASE);
    }

    /**
     * @brief Generate an ngspice circuit for DMC LC filter testing.
     */
    std::string generate_ngspice_circuit(double inductance, double frequency = 150000);

    /**
     * @brief Simulate DMC and extract waveforms.
     *
     * @param numberOfPeriods  Number of steady-state periods to extract per
     *                         frequency (defaults to 1). Mirrors the
     *                         pattern used by every other converter
     *                         (Buck, Boost, AHB, Cllc, ...): the ngspice
     *                         runner is configured with extractOnePeriod
     *                         + numberOfPeriods so the caller never has
     *                         to tile waveforms post-hoc.
     */
    std::vector<DmcSimulationWaveforms> simulate_and_extract_waveforms(
        double inductance,
        const std::vector<double>& frequencies,
        size_t numberOfPeriods = 1);

    /**
     * @brief Simulate and extract operating points from simulation.
     */
    std::vector<OperatingPoint> simulate_and_extract_operating_points(double inductance);

    /**
     * @brief Verify that a design meets the attenuation requirements using ngspice.
     */
    std::vector<DmcAttenuationResult> verify_attenuation(
        double inductance,
        std::optional<double> capacitance = std::nullopt);

    /**
     * @brief Calculate the required inductance to meet attenuation requirements.
     */
    static double calculate_required_inductance(
        double targetAttenuation,
        double frequency,
        double capacitance);

    /**
     * @brief Propose a design based on the specifications.
     */
    json propose_design();

    // ---- Diagnostics (populated by process_design_requirements) ----
    double get_computed_inductance()         const { return computedInductance; }
    double get_computed_min_frequency()      const { return computedMinFrequency; }
    double get_computed_max_frequency()      const { return computedMaxFrequency; }
    double get_computed_impedance_at_min_freq() const { return computedImpedanceAtMinFreq; }
    int    get_computed_number_windings()    const { return computedNumberWindings; }

private:
    mutable double computedInductance         = 0.0;
    mutable double computedMinFrequency       = 0.0;  // lowest impedance spec frequency [Hz]
    mutable double computedMaxFrequency       = 0.0;  // highest impedance spec frequency [Hz]
    mutable double computedImpedanceAtMinFreq = 0.0;  // |Z| at the lowest spec freq [Ω]
    mutable int    computedNumberWindings     = 1;
};

} // namespace OpenMagnetics
