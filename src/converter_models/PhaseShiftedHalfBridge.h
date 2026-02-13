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
 * @brief Phase-Shifted Half Bridge (PSHB) DC-DC Converter
 *
 * Inherits from MAS::PhaseShiftFullBridge and the Topology interface.
 * Reuses the same JSON/MAS schema as the full-bridge variant but applies
 * a half-bridge voltage factor (Vin/2) throughout.
 *
 * =====================================================================
 * TOPOLOGY OVERVIEW
 * =====================================================================
 *
 * Half-bridge with split-capacitor input and phase-shift control:
 *
 *          +Vin ----+---- C1 ----+---- midCap
 *                   |            |
 *                  [QA]          |
 *                   |            |
 *                midSW ---[Lr]---+---[T1 Np:Ns]---[Rectifier]---[Lo]---Vo
 *                   |
 *                  [QB]
 *                   |
 *          GND -----+
 *
 * QA and QB switch complementarily at ~50% duty.  The transformer
 * primary is connected between the switch mid-point and the capacitor
 * mid-point.  The capacitor divider forces the primary voltage to
 * swing between +Vin/2 and -Vin/2.
 *
 * "Phase-shifted" refers to the timing/overlap control of the two
 * switches relative to each other (analogous to the leading/lagging
 * leg concept in the full-bridge, but here applied within a single
 * leg via asymmetric duty or dead-time modulation).  The effective
 * duty cycle D_eff controls the power transfer.
 *
 * Alternatively, in literature this is sometimes called an
 * "asymmetric half-bridge" or "complementary half-bridge."
 *
 * =====================================================================
 * KEY EQUATIONS
 * =====================================================================
 *
 * Effective primary voltage amplitude:
 *   Vpri_pk = Vin / 2   (half of full-bridge)
 *
 * Output voltage (center-tapped rectifier):
 *   Vo = (Vin/2) * D_eff / n - Vd
 *
 * Output voltage (current doubler):
 *   Vo = (Vin/2) * D_eff / (2*n) - Vd
 *
 * Output voltage (full-bridge rectifier):
 *   Vo = (Vin/2) * D_eff / n - 2*Vd
 *
 * Turns ratio (center-tapped, for target D_eff at nominal Vin):
 *   n = (Vin_nom/2) * D_eff_nom / (Vo + Vd)
 *
 * Effective duty cycle:
 *   D_eff = phase_shift / 180  (phase_shift in degrees)
 *
 * Primary voltage waveform (3-level, same shape as PSFB but half amplitude):
 *   +(Vin/2)  during power transfer
 *   0         during freewheeling
 *   -(Vin/2)  during opposite power transfer
 *   0         during opposite freewheeling
 *
 * Primary current: same shape as PSFB, but currents are higher for
 *   the same output power because reflected load current = Io/n and
 *   n is smaller (due to half the primary voltage).
 *
 * Magnetizing inductance:
 *   Im_peak = (Vin/2) * D_eff / (4 * Fs * Lm)
 *   Lm = (Vin/2) * D_eff / (4 * Fs * Im_target)
 *
 * Output inductor:
 *   Lo = Vo * (1 - D_eff) / (Fs * dIo)
 *
 * Series inductance (ZVS assist):
 *   0.5 * Lr * Ip^2 > 0.5 * Coss * (Vin/2)^2
 *   Lr_min = Coss * (Vin/2)^2 / Ip^2
 *
 * Compared to PSFB:
 *   - Half the number of primary switches (2 vs 4)
 *   - Half the primary voltage swing (Vin/2 vs Vin)
 *   - For the same output, turns ratio n is ~half → higher primary currents
 *   - Simpler gate drive (only one leg to drive)
 *   - Typically suited for medium power (up to ~500 W)
 *   - Split capacitors must handle full primary current ripple
 */
class Pshb : public MAS::PhaseShiftFullBridge, public Topology {
private:
    int numPeriodsToExtract = 5;
    int numSteadyStatePeriods = 5;

    double computedOutputInductance = 0;
    double computedSeriesInductance = 0;
    double computedMagnetizingInductance = 0;
    double computedDeadTime = 200e-9;
    double computedEffectiveDutyCycle = 0;
    double computedDiodeVoltageDrop = 0.6;

    // Half-bridge voltage factor
    static constexpr double BRIDGE_VOLTAGE_FACTOR = 0.5;

public:
    bool _assertErrors = false;

    Pshb(const json& j);
    Pshb() {};

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
    double get_bridge_voltage_factor() const { return BRIDGE_VOLTAGE_FACTOR; }

    // ---- Topology interface ----
    bool run_checks(bool assert = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const PsfbOperatingPoint& pshbOpPoint,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // ---- PSHB-specific calculations ----
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


class AdvancedPshb : public Pshb {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredMagnetizingInductance;
    std::optional<double> desiredSeriesInductance;
    std::optional<double> desiredOutputInductance;

public:
    AdvancedPshb() = default;
    ~AdvancedPshb() = default;
    AdvancedPshb(const json& j);

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

void from_json(const json& j, AdvancedPshb& x);
void to_json(json& j, const AdvancedPshb& x);

inline void from_json(const json& j, AdvancedPshb& x) {
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

inline void to_json(json& j, const AdvancedPshb& x) {
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
