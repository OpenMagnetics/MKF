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
 * @brief Phase-Shifted Full Bridge (PSFB) DC-DC Converter
 *
 * Inherits from MAS::PhaseShiftFullBridge and the Topology interface.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW
 * =====================================================================
 *
 *   +--[QA]---+---[QC]--+
 *   |         |         |
 *  Vin    bridge_A   Vin
 *   |         |         |
 *   +--[QB]---+---[QD]--+
 *             |   Lk/Lr
 *             +--/\/\/---[T1 Np:Ns]---[Rectifier]---[Lo]---Vo
 *
 * Full-bridge primary with phase-shift control between legs.
 * QA-QB form the "active" (leading) leg.
 * QC-QD form the "passive" (lagging) leg.
 * All switches operate at 50% duty, the phase shift between legs
 * determines the effective duty cycle and power transfer.
 *
 * Rectifier options:
 *   - CENTER_TAPPED: two diodes, center-tapped secondary
 *   - CURRENT_DOUBLER: two inductors on secondary
 *   - FULL_BRIDGE: four diodes
 *
 * =====================================================================
 * KEY EQUATIONS
 * =====================================================================
 *
 * References:
 *   [1] TI TIDU248 - Phase-Shifted Full Bridge Design Guide
 *   [2] Sabate et al., "Design Considerations for High-Voltage
 *       High-Power Full-Bridge Zero-Voltage-Switched PWM Converter"
 *
 * Effective duty cycle (ratio of overlap between diagonal switches):
 *   D_eff = phase_shift / 180  (phase_shift in degrees)
 *
 * Output voltage (center-tapped rectifier):
 *   Vo = Vin * D_eff / n - Vd
 *   where n = Np/Ns (per half-secondary)
 *
 * Output voltage (current doubler):
 *   Vo = Vin * D_eff / (2*n) - Vd
 *
 * Output voltage (full-bridge rectifier):
 *   Vo = Vin * D_eff / n - 2*Vd
 *
 * Turns ratio (from desired D_eff at nominal):
 *   n = Vin_nom * D_eff_nom / (Vo + Vd)           (center-tapped)
 *   n = Vin_nom * D_eff_nom / (2 * (Vo + Vd))     (current doubler)
 *   n = Vin_nom * D_eff_nom / (Vo + 2*Vd)         (full-bridge rect.)
 *
 * Primary voltage waveform (3-level):
 *   +Vin during power transfer (diagonal pair ON)
 *   0V during freewheeling (same-leg pair ON)
 *   -Vin during opposite power transfer
 *   0V during opposite freewheeling
 *
 * Primary current:
 *   During power transfer: i_pri rises with slope (Vin - n*Vo) / Lk
 *   During freewheeling: i_pri circulates, slow decay
 *   Transitions have resonant edges (ZVS)
 *
 * Output inductor (single inductor, center-tapped/full-bridge rect.):
 *   Lo = Vo * (1 - D_eff) / (Fs * delta_Io)
 *   where delta_Io = r * Io (r = current_ripple_ratio)
 *
 * Magnetizing inductance:
 *   Primary sees bipolar rectangular +/-Vin (each for D_eff*Ts/2)
 *   B_peak = Vin * D_eff / (2 * Fs * Np * Ae)
 *   Lm = Vin / (4 * Fs * Im_peak)
 *
 * Series inductance (leakage + external Lr):
 *   Provides ZVS energy: 0.5*Lr*Ip^2 > 0.5*Coss*Vin^2
 *   Lr_min = Coss * Vin^2 / Ip^2
 *   Duty cycle loss: t_loss = Lr * Io_ref / (n * Vin)
 */
class Psfb : public MAS::PhaseShiftFullBridge, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 5;

    double computedOutputInductance = 0;
    double computedSeriesInductance = 0;
    double computedMagnetizingInductance = 0;
    double computedDeadTime = 200e-9;
    double computedEffectiveDutyCycle = 0;
    double computedDiodeVoltageDrop = 0.6;  // Default

public:
    bool _assertErrors = false;

    Psfb(const json& j);
    Psfb() {};

    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int value) { numPeriodsToExtract = value; }
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int value) { numSteadyStatePeriods = value; }

    double get_computed_output_inductance() const { return computedOutputInductance; }
    void set_computed_output_inductance(double value) { computedOutputInductance = value; }
    double get_computed_series_inductance() const { return computedSeriesInductance; }
    void set_computed_series_inductance(double value) { computedSeriesInductance = value; }
    double get_computed_magnetizing_inductance() const { return computedMagnetizingInductance; }
    void set_computed_magnetizing_inductance(double value) { computedMagnetizingInductance = value; }
    double get_computed_dead_time() const { return computedDeadTime; }
    void set_computed_dead_time(double value) { computedDeadTime = value; }
    double get_computed_effective_duty_cycle() const { return computedEffectiveDutyCycle; }

    // ---- Topology interface ----
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const PsfbOperatingPoint& psfbOpPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // ---- PSFB-specific calculations ----
    static double compute_effective_duty_cycle(double phaseShiftDeg);
    static double compute_output_voltage(double Vin, double Deff, double n,
                                         double Vd, PsfbRectifierType rectType);
    static double compute_turns_ratio(double Vin, double Vo, double Deff,
                                      double Vd, PsfbRectifierType rectType);
    static double compute_output_inductance(double Vo, double Deff, double Fs,
                                            double Io, double rippleRatio);
    static double compute_primary_rms_current(double Io, double n, double Deff);

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


class AdvancedPsfb : public Psfb {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredMagnetizingInductance;
    std::optional<double> desiredSeriesInductance;
    std::optional<double> desiredOutputInductance;

public:
    AdvancedPsfb() = default;
    ~AdvancedPsfb() = default;
    AdvancedPsfb(const json& j);

    Inputs process();

    const double& get_desired_magnetizing_inductance() const { return desiredMagnetizingInductance; }
    void set_desired_magnetizing_inductance(const double& value) { desiredMagnetizingInductance = value; }
    const std::vector<double>& get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double>& value) { desiredTurnsRatios = value; }
    std::optional<double> get_desired_series_inductance() const { return desiredSeriesInductance; }
    void set_desired_series_inductance(std::optional<double> value) { desiredSeriesInductance = value; }
    std::optional<double> get_desired_output_inductance() const { return desiredOutputInductance; }
    void set_desired_output_inductance(std::optional<double> value) { desiredOutputInductance = value; }
};

void from_json(const json& j, AdvancedPsfb& x);
void to_json(json& j, const AdvancedPsfb& x);

inline void from_json(const json& j, AdvancedPsfb& x) {
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_maximum_phase_shift(get_stack_optional<double>(j, "maximumPhaseShift"));
    x.set_operating_points(j.at("operatingPoints").get<std::vector<PsfbOperatingPoint>>());
    x.set_output_inductance(get_stack_optional<double>(j, "outputInductance"));
    x.set_rectifier_type(get_stack_optional<PsfbRectifierType>(j, "rectifierType"));
    x.set_series_inductance(get_stack_optional<double>(j, "seriesInductance"));
    x.set_use_leakage_inductance(get_stack_optional<bool>(j, "useLeakageInductance"));
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_desired_magnetizing_inductance(j.at("desiredMagnetizingInductance").get<double>());
    x.set_desired_series_inductance(get_stack_optional<double>(j, "desiredSeriesInductance"));
    x.set_desired_output_inductance(get_stack_optional<double>(j, "desiredOutputInductance"));
}

inline void to_json(json& j, const AdvancedPsfb& x) {
    j = json::object();
    j["efficiency"] = x.get_efficiency();
    j["inputVoltage"] = x.get_input_voltage();
    j["maximumPhaseShift"] = x.get_maximum_phase_shift();
    j["operatingPoints"] = x.get_operating_points();
    j["outputInductance"] = x.get_output_inductance();
    j["rectifierType"] = x.get_rectifier_type();
    j["seriesInductance"] = x.get_series_inductance();
    j["useLeakageInductance"] = x.get_use_leakage_inductance();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredMagnetizingInductance"] = x.get_desired_magnetizing_inductance();
    j["desiredSeriesInductance"] = x.get_desired_series_inductance();
    j["desiredOutputInductance"] = x.get_desired_output_inductance();
}

} // namespace OpenMagnetics
