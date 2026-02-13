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
 * @brief Dual Active Bridge (DAB) DC-DC Converter
 *
 * Inherits converter-level parameters from MAS::DualActiveBridge and
 * the Topology interface, exactly mirroring the LLC pattern.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW
 * =====================================================================
 *
 *   +--[Q1]--+--[Q3]--+             +--[Q5]--+--[Q7]--+
 *   |        |        |             |        |        |
 *  V1     bridge_a  V1           V2     bridge_c  V2
 *   |        |        |             |        |        |
 *   +--[Q2]--+--[Q4]--+             +--[Q6]--+--[Q8]--+
 *            |     L_series              |
 *            +---/\/\/---[Transformer]---+
 *                         Np : Ns
 *
 * Two full H-bridges connected through a high-frequency transformer
 * and a series inductance L (leakage or external coupling inductor).
 * Power transfer controlled by phase shift phi between the two bridges.
 *
 * =====================================================================
 * KEY EQUATIONS (Single Phase Shift / SPS modulation)
 * =====================================================================
 *
 * References:
 *   [1] TI TIDA-010054 (tidues0e.pdf) - Section 2.3
 *   [2] Demetriades PhD thesis (FULLTEXT01.pdf) - Chapter 6
 *   [3] Shao et al. IEEE TPEL 2021 - DAB Modeling & Control Review
 *
 * Turns ratio:
 *   N = V1_nom / V2_nom  (primary-to-secondary)
 *
 * Voltage conversion ratio:
 *   d = N * V2 / V1   (d = 1 at nominal operating point)
 *
 * Power transfer (SPS):
 *   P = N * V1 * V2 * phi * (pi - |phi|) / (2 * pi^2 * Fs * L)
 *
 * Phase shift for desired power:
 *   phi = (pi/2) * (1 - sqrt(1 - 8*Fs*L*P / (N*V1*V2)))
 *
 * Inductor current (piecewise linear, referred to primary):
 *   Ibase = V1 / (2*pi*Fs*L)
 *   i1 = 0.5 * (2*phi - (1-d)*pi) * Ibase
 *   i2 = 0.5 * (2*d*phi + (1-d)*pi) * Ibase
 *
 * ZVS boundaries:
 *   Primary:   phi_zvs > (1 - 1/d) * pi/2
 *   Secondary: phi_zvs > (1 - d) * pi/2
 *
 * Waveform: 4 intervals per period, piecewise linear inductor current.
 *   Interval 1 (0 to t_phi):       slope = (V1 + N*V2)/L
 *   Interval 2 (t_phi to Ts/2):    slope = (V1 - N*V2)/L
 *   Interval 3 (Ts/2 to Ts/2+t_phi): slope = -(V1 + N*V2)/L
 *   Interval 4 (Ts/2+t_phi to Ts): slope = -(V1 - N*V2)/L
 *
 * Transformer excitation (for MAS OperatingPoint):
 *   Primary:   V(t) = bipolar rectangular +/- V1
 *              I(t) = piecewise linear iL(t)
 *   Secondary: V(t) = bipolar rectangular +/- V2
 *              I(t) = N * iL(t)
 *
 * Magnetizing inductance:
 *   Lm is the transformer magnetizing inductance.
 *   Lm >> L (typically 10x-100x the series inductance).
 *   Lm determines the magnetizing current ripple.
 *   The magnetizing current is a triangular wave superimposed on iL.
 *   For MAS design: Lm_min = V1 / (4 * Fs * Delta_Im_max)
 */
class Dab : public MAS::DualActiveBridge, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 5;

    // Computed design values (filled by process_design_requirements)
    double computedSeriesInductance = 0;     // L (leakage + external)
    double computedMagnetizingInductance = 0; // Lm of transformer
    double computedDeadTime = 200e-9;         // Default 200 ns dead time
    double computedPhaseShift = 0;            // phi in radians (computed if not given)

public:
    bool _assertErrors = false;

    Dab(const json& j);
    Dab() {};

    // ---- Simulation tuning ----
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { numPeriodsToExtract = value; }
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { numSteadyStatePeriods = value; }

    // ---- Computed value accessors ----
    double get_computed_series_inductance() const { return computedSeriesInductance; }
    void set_computed_series_inductance(double value) { computedSeriesInductance = value; }

    double get_computed_magnetizing_inductance() const { return computedMagnetizingInductance; }
    void set_computed_magnetizing_inductance(double value) { computedMagnetizingInductance = value; }

    double get_computed_dead_time() const { return computedDeadTime; }
    void set_computed_dead_time(double value) { computedDeadTime = value; }

    double get_computed_phase_shift() const { return computedPhaseShift; }

    // ---- Topology interface ----
    bool run_checks(bool assert = false) override;

    DesignRequirements process_design_requirements() override;

    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    // Per-operating-point analytical waveforms
    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const DabOperatingPoint& dabOpPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // ---- DAB-specific calculations ----

    /** Compute power transfer for given parameters */
    static double compute_power(double V1, double V2, double N,
                                double phi, double Fs, double L);

    /** Compute required series inductance for desired power at given phase shift */
    static double compute_series_inductance(double V1, double V2, double N,
                                            double phi, double Fs, double P);

    /** Compute phase shift for desired power with given inductance */
    static double compute_phase_shift(double V1, double V2, double N,
                                      double Fs, double L, double P);

    /** Compute voltage conversion ratio d = N*V2/V1 */
    static double compute_voltage_ratio(double V1, double V2, double N);

    /** Compute inductor current at switching instants i1 and i2 */
    static void compute_switching_currents(double V1, double V2, double N,
                                           double phi, double Fs, double L,
                                           double& i1, double& i2);

    /** Compute primary RMS current */
    static double compute_primary_rms_current(double i1, double i2, double phi);

    /** Check ZVS condition for both bridges */
    static bool check_zvs_primary(double phi, double d);
    static bool check_zvs_secondary(double phi, double d);

    // ---- SPICE simulation ----
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
        double magnetizingInductance);
};


/**
 * @brief AdvancedDab: user supplies desired turns ratios & inductances directly
 */
class AdvancedDab : public Dab {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredMagnetizingInductance;
    std::optional<double> desiredSeriesInductance;

public:
    AdvancedDab() = default;
    ~AdvancedDab() = default;

    AdvancedDab(const json& j);

    Inputs process();

    const double& get_desired_magnetizing_inductance() const { return desiredMagnetizingInductance; }
    void set_desired_magnetizing_inductance(const double& value) { desiredMagnetizingInductance = value; }

    const std::vector<double>& get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double>& value) { desiredTurnsRatios = value; }

    std::optional<double> get_desired_series_inductance() const { return desiredSeriesInductance; }
    void set_desired_series_inductance(std::optional<double> value) { desiredSeriesInductance = value; }
};


// ---- JSON serialization ----
void from_json(const json& j, AdvancedDab& x);
void to_json(json& j, const AdvancedDab& x);

inline void from_json(const json& j, AdvancedDab& x) {
    // DualActiveBridge base fields
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<DabOperatingPoint>>());
    x.set_series_inductance(get_stack_optional<double>(j, "seriesInductance"));
    x.set_use_leakage_inductance(get_stack_optional<bool>(j, "useLeakageInductance"));

    // AdvancedDab extra fields
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_magnetizing_inductance(j.at("desiredMagnetizingInductance").get<double>());
    x.set_desired_series_inductance(get_stack_optional<double>(j, "desiredSeriesInductance"));
}

inline void to_json(json& j, const AdvancedDab& x) {
    j = json::object();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["operatingPoints"] = x.get_operating_points();
    j["seriesInductance"] = x.get_series_inductance();
    j["useLeakageInductance"] = x.get_use_leakage_inductance();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredMagnetizingInductance"] = x.get_desired_magnetizing_inductance();
    j["desiredSeriesInductance"] = x.get_desired_series_inductance();
}

} // namespace OpenMagnetics
