#pragma once
// Version: 1.0 — initial CLLLC bidirectional symmetric resonant converter.
// See src/converter_models/CLLLC_PLAN.md for full design rationale.
//
// CLLLC vs neighbours:
//  - LLC:  one Lr (primary) + one Cr (primary), parallel Lm. Forward only.
//  - CLLC: one Lr (primary) + TWO Cr, parallel Lm. Bidirectional, asymmetric.
//  - CLLLC: TWO Lr (one each side) + TWO Cr (one each side) + parallel Lm.
//           Bidirectional and SYMMETRIC if Lr2=Lr1/n^2 and Cr2=Cr1*n^2.
//  - DAB:  no resonant tank, just transformer leakage. Bidirectional, 8 sw.
//  - SRC/PRC: series-resonant only (no Lm parallel branch).
//
// Headline diagnostics that justify a separate model from LLC and CLLC:
//  - lastZvsMarginPrimaryLagging AND lastZvsMarginSecondaryLagging
//    (two ZVS conditions, one per bridge).
//  - lastCurrentSharingRatio = max(i_Lr1) / (max(i_Lr2)/n) — should be
//    1.0 +/- 5% when tankSymmetryRatio == 1.0.
//
// References (full bibliography in CLLLC_PLAN.md sec. 9):
//  TI TIDM-02002 (6.6 kW OBC), TI TIDM-02013 (7.4 kW OBC),
//  MDPI Electronics 12(7):1605, IEEE 6972179 (time-domain), Huang SLUP263.
//
// Accuracy disclaimer (v1):
//  - FHA gain accuracy degrades farther from f_r1 than for LLC.
//  - Coss(Vds) variation not modelled.
//  - Body-diode reverse recovery omitted.
//  - Direction-reversal transient near gain=1 NOT modelled; v1 picks one
//    direction per OP from powerFlowDirection.
//  - Inter-winding capacitance not modelled.
//  - Leakage tolerances assumed ideal.

#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "converter_models/Topology.h"
#include "processors/NgspiceRunner.h"

namespace OpenMagnetics {
using namespace MAS;

class Clllc : public MAS::ClllcResonant, public Topology {
public:
    enum class PowerFlowDirection { FORWARD, REVERSE };
    enum class OperatingRegion { ABOVE_RESONANCE, AT_RESONANCE, BELOW_RESONANCE };

private:
    // Simulation tuning
    int numPeriodsToExtract   = 5;
    int numSteadyStatePeriods = 4;   // 5-state system needs more than LLC

    // Computed design values (filled by process_design_requirements)
    double computedPrimarySeriesInductance      = 0;   // Lr1
    double computedSecondarySeriesInductance    = 0;   // Lr2 = Lr1/n^2
    double computedPrimaryResonantCapacitance   = 0;   // Cr1
    double computedSecondaryResonantCapacitance = 0;   // Cr2 = Cr1*n^2
    double computedMagnetizingInductance        = 0;   // Lm
    double computedTurnsRatio                   = 0;   // n = Np/Ns
    double computedDeadTime                     = 200e-9;
    double computedInductanceRatioK             = 0;   // K = Lm/Lr1
    double computedQualityFactor                = 0;
    double computedPrimaryResonantFrequency     = 0;   // f_r1
    double computedMagnetizingResonantFreqFwd   = 0;   // f_m1
    double computedMagnetizingResonantFreqRev   = 0;   // f_m2 (= f_m1 if sym)

    // User overrides (skip MAS schema; settable via C++)
    double userPrimarySeriesInductance    = 0; // 0 = unset
    double userPrimaryResonantCapacitance = 0; // 0 = unset

    // Schema-less device parameters (Coss for ZVS sizing)
    double mosfetCossPrimary   = 200e-12;
    double mosfetCossSecondary = 400e-12;

    // Per-OP diagnostics
    mutable PowerFlowDirection lastPowerFlowDirection = PowerFlowDirection::FORWARD;
    mutable OperatingRegion lastOperatingRegion       = OperatingRegion::AT_RESONANCE;
    mutable int lastModeForward = 1;
    mutable int lastModeReverse = 1;

    // ZVS diagnostics — populated for BOTH bridges (headline)
    mutable double lastZvsMarginPrimaryLagging   = 0.0;
    mutable double lastZvsMarginSecondaryLagging = 0.0;
    mutable double lastZvsLoadThresholdPrimary   = 0.0;
    mutable double lastZvsLoadThresholdSecondary = 0.0;
    mutable double lastResonantTransitionTime    = 0.0;

    mutable double lastPrimaryPeakCurrent     = 0.0;
    mutable double lastSecondaryPeakCurrent   = 0.0;
    mutable double lastPrimaryRmsCurrent      = 0.0;
    mutable double lastSecondaryRmsCurrent    = 0.0;
    mutable double lastMagnetizingPeakCurrent = 0.0;
    mutable double lastCr1PeakVoltage         = 0.0;
    mutable double lastCr2PeakVoltage         = 0.0;
    mutable double lastCurrentSharingRatio    = 1.0;
    mutable double lastSteadyStateResidual    = 0.0;

    // Extra-component waveforms (cleared in process_operating_points)
    mutable std::vector<Waveform> extraLr1CurrentWaveforms;
    mutable std::vector<Waveform> extraLr1VoltageWaveforms;
    mutable std::vector<Waveform> extraLr2CurrentWaveforms;
    mutable std::vector<Waveform> extraLr2VoltageWaveforms;
    mutable std::vector<Waveform> extraCr1CurrentWaveforms;
    mutable std::vector<Waveform> extraCr1VoltageWaveforms;
    mutable std::vector<Waveform> extraCr2CurrentWaveforms;
    mutable std::vector<Waveform> extraCr2VoltageWaveforms;
    mutable std::vector<Waveform> extraLmCurrentWaveforms;
    mutable std::vector<std::vector<double>> extraTimeVectors;

public:
    bool _assertErrors = false;

    Clllc(const json& j);
    Clllc() {}

    // Tuning ─────────────────────────────────────────────────
    int get_num_periods_to_extract() const { return numPeriodsToExtract; }
    void set_num_periods_to_extract(int v) { this->numPeriodsToExtract = v; }
    int get_num_steady_state_periods() const { return numSteadyStatePeriods; }
    void set_num_steady_state_periods(int v) { this->numSteadyStatePeriods = v; }

    // Computed accessors ─────────────────────────────────────
    double get_computed_primary_series_inductance() const      { return computedPrimarySeriesInductance; }
    double get_computed_secondary_series_inductance() const    { return computedSecondarySeriesInductance; }
    double get_computed_primary_resonant_capacitance() const   { return computedPrimaryResonantCapacitance; }
    double get_computed_secondary_resonant_capacitance() const { return computedSecondaryResonantCapacitance; }
    double get_computed_magnetizing_inductance() const         { return computedMagnetizingInductance; }
    double get_computed_turns_ratio() const                    { return computedTurnsRatio; }
    double get_computed_dead_time() const                      { return computedDeadTime; }
    double get_computed_inductance_ratio_k() const             { return computedInductanceRatioK; }
    double get_computed_quality_factor() const                 { return computedQualityFactor; }
    double get_computed_primary_resonant_frequency() const     { return computedPrimaryResonantFrequency; }
    double get_computed_magnetizing_resonant_freq_fwd() const  { return computedMagnetizingResonantFreqFwd; }
    double get_computed_magnetizing_resonant_freq_rev() const  { return computedMagnetizingResonantFreqRev; }
    void set_computed_dead_time(double v) { this->computedDeadTime = v; }

    // User overrides for primary tank components
    void set_user_primary_series_inductance(double v)    { userPrimarySeriesInductance = v; }
    void set_user_primary_resonant_capacitance(double v) { userPrimaryResonantCapacitance = v; }
    double get_user_primary_series_inductance() const    { return userPrimarySeriesInductance; }
    double get_user_primary_resonant_capacitance() const { return userPrimaryResonantCapacitance; }

    // Coss (for ZVS sizing diagnostics)
    void set_mosfet_coss_primary(double v)   { mosfetCossPrimary = v; }
    void set_mosfet_coss_secondary(double v) { mosfetCossSecondary = v; }
    double get_mosfet_coss_primary() const   { return mosfetCossPrimary; }
    double get_mosfet_coss_secondary() const { return mosfetCossSecondary; }

    // Per-OP diagnostics ─────────────────────────────────────
    PowerFlowDirection get_last_power_flow_direction() const { return lastPowerFlowDirection; }
    OperatingRegion    get_last_operating_region() const     { return lastOperatingRegion; }
    int get_last_mode_forward() const { return lastModeForward; }
    int get_last_mode_reverse() const { return lastModeReverse; }
    double get_last_zvs_margin_primary_lagging() const   { return lastZvsMarginPrimaryLagging; }
    double get_last_zvs_margin_secondary_lagging() const { return lastZvsMarginSecondaryLagging; }
    double get_last_zvs_load_threshold_primary() const   { return lastZvsLoadThresholdPrimary; }
    double get_last_zvs_load_threshold_secondary() const { return lastZvsLoadThresholdSecondary; }
    double get_last_resonant_transition_time() const     { return lastResonantTransitionTime; }
    double get_last_primary_peak_current() const     { return lastPrimaryPeakCurrent; }
    double get_last_secondary_peak_current() const   { return lastSecondaryPeakCurrent; }
    double get_last_primary_rms_current() const      { return lastPrimaryRmsCurrent; }
    double get_last_secondary_rms_current() const    { return lastSecondaryRmsCurrent; }
    double get_last_magnetizing_peak_current() const { return lastMagnetizingPeakCurrent; }
    double get_last_cr1_peak_voltage() const { return lastCr1PeakVoltage; }
    double get_last_cr2_peak_voltage() const { return lastCr2PeakVoltage; }
    double get_last_current_sharing_ratio() const { return lastCurrentSharingRatio; }
    double get_last_steady_state_residual() const { return lastSteadyStateResidual; }

    // Bridge factor (0.5 HB / 1.0 FB) per side
    double get_primary_bridge_voltage_factor() const;
    double get_secondary_bridge_voltage_factor() const;

    // Effective resonant frequency (user value or geo-mean of fsw range)
    double get_effective_resonant_frequency() const;

    // Effective inductance ratio K (user override or computed)
    double get_effective_inductance_ratio_k() const;
    // Effective tank symmetry (defaults to 1.0)
    double get_effective_tank_symmetry_ratio() const;

    // Topology interface ─────────────────────────────────────
    bool run_checks(bool assertErrors = false) override;
    DesignRequirements process_design_requirements() override;
    std::vector<OperatingPoint> process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) override;

    OperatingPoint process_operating_point_for_input_voltage(
        double inputVoltage,
        const ClllcOperatingPoint& op,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance);

    // Static analytical helpers
    static double compute_primary_resonant_frequency(double Lr1, double Cr1);
    static double compute_magnetizing_resonant_frequency_fwd(double Lr1, double Lm, double Cr1);
    static double compute_magnetizing_resonant_frequency_rev(double Lr2, double Lm, double n, double Cr2);
    static double compute_fha_gain(double f_n, double K, double Q, double a, double b);
    static double compute_inductance_ratio_K(double Lm, double Lr1);
    static double compute_quality_factor(double Lr1, double Cr1, double R_load_reflected);
    static double compute_zvs_lm_max(double T_dead, double C_oss, double Fs);
    static double compute_turns_ratio(double V_HV, double V_LV);
    static double compute_symmetric_lr2(double Lr1, double n);
    static double compute_symmetric_cr2(double Cr1, double n);

    // Extra components factory
    std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt) override;

    // SPICE
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


// ─────────────────────────────────────────────────────────────
// AdvancedClllc: bypass design step; user supplies all components.
// ─────────────────────────────────────────────────────────────
class AdvancedClllc : public Clllc {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredMagnetizingInductance         = 0;
    std::optional<double> desiredPrimarySeriesInductance;
    std::optional<double> desiredPrimaryResonantCapacitance;

public:
    AdvancedClllc() = default;
    ~AdvancedClllc() = default;
    AdvancedClllc(const json& j);

    Inputs process();

    const double& get_desired_magnetizing_inductance() const { return desiredMagnetizingInductance; }
    void set_desired_magnetizing_inductance(const double& v) { this->desiredMagnetizingInductance = v; }
    const std::vector<double>& get_desired_turns_ratios() const { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double>& v) { this->desiredTurnsRatios = v; }
    std::optional<double> get_desired_primary_series_inductance() const { return desiredPrimarySeriesInductance; }
    void set_desired_primary_series_inductance(std::optional<double> v) { this->desiredPrimarySeriesInductance = v; }
    std::optional<double> get_desired_primary_resonant_capacitance() const { return desiredPrimaryResonantCapacitance; }
    void set_desired_primary_resonant_capacitance(std::optional<double> v) { this->desiredPrimaryResonantCapacitance = v; }
};


// ─────────────────────────────────────────────────────────────
// JSON serialization for AdvancedClllc
// ─────────────────────────────────────────────────────────────
void from_json(const json& j, AdvancedClllc& x);
void to_json(json& j, const AdvancedClllc& x);

inline void from_json(const json& j, AdvancedClllc& x) {
    // ClllcResonant base fields
    x.set_bridge_type_primary(get_stack_optional<LlcBridgeType>(j, "bridgeTypePrimary"));
    x.set_bridge_type_secondary(get_stack_optional<LlcBridgeType>(j, "bridgeTypeSecondary"));
    x.set_control_strategy(get_stack_optional<ClllcControlStrategy>(j, "controlStrategy"));
    x.set_efficiency(get_stack_optional<double>(j, "efficiency"));
    x.set_high_voltage_bus_voltage(j.at("highVoltageBusVoltage").get<DimensionWithTolerance>());
    x.set_low_voltage_bus_voltage(j.at("lowVoltageBusVoltage").get<DimensionWithTolerance>());
    x.set_inductance_ratio_k(get_stack_optional<double>(j, "inductanceRatioK"));
    x.set_integrated_resonant_inductors(get_stack_optional<bool>(j, "integratedResonantInductors"));
    x.set_max_switching_frequency(j.at("maxSwitchingFrequency").get<double>());
    x.set_min_switching_frequency(j.at("minSwitchingFrequency").get<double>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<ClllcOperatingPoint>>());
    x.set_primary_resonant_capacitance(get_stack_optional<double>(j, "primaryResonantCapacitance"));
    x.set_primary_resonant_frequency(get_stack_optional<double>(j, "primaryResonantFrequency"));
    x.set_primary_series_inductance(get_stack_optional<double>(j, "primarySeriesInductance"));
    x.set_quality_factor(get_stack_optional<double>(j, "qualityFactor"));
    x.set_tank_symmetry_ratio(get_stack_optional<double>(j, "tankSymmetryRatio"));

    // Advanced extra fields
    auto desiredTurnsRatios = get_stack_optional<std::vector<double>>(j, "desiredTurnsRatios");
    if (desiredTurnsRatios.has_value()) {
        x.set_desired_turns_ratios(desiredTurnsRatios.value());
    }
    auto desiredMagnetizingInductance = get_stack_optional<double>(j, "desiredMagnetizingInductance");
    if (desiredMagnetizingInductance.has_value()) {
        x.set_desired_magnetizing_inductance(desiredMagnetizingInductance.value());
    }
    x.set_desired_primary_series_inductance(get_stack_optional<double>(j, "desiredPrimarySeriesInductance"));
    x.set_desired_primary_resonant_capacitance(get_stack_optional<double>(j, "desiredPrimaryResonantCapacitance"));
}

inline void to_json(json& j, const AdvancedClllc& x) {
    j = json::object();
    j["bridgeTypePrimary"] = x.get_bridge_type_primary();
    j["bridgeTypeSecondary"] = x.get_bridge_type_secondary();
    j["controlStrategy"] = x.get_control_strategy();
    j["efficiency"] = x.get_efficiency();
    j["highVoltageBusVoltage"] = x.get_high_voltage_bus_voltage();
    j["lowVoltageBusVoltage"] = x.get_low_voltage_bus_voltage();
    j["inductanceRatioK"] = x.get_inductance_ratio_k();
    j["integratedResonantInductors"] = x.get_integrated_resonant_inductors();
    j["maxSwitchingFrequency"] = x.get_max_switching_frequency();
    j["minSwitchingFrequency"] = x.get_min_switching_frequency();
    j["operatingPoints"] = x.get_operating_points();
    j["primaryResonantCapacitance"] = x.get_primary_resonant_capacitance();
    j["primaryResonantFrequency"] = x.get_primary_resonant_frequency();
    j["primarySeriesInductance"] = x.get_primary_series_inductance();
    j["qualityFactor"] = x.get_quality_factor();
    j["tankSymmetryRatio"] = x.get_tank_symmetry_ratio();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["desiredMagnetizingInductance"] = x.get_desired_magnetizing_inductance();
    j["desiredPrimarySeriesInductance"] = x.get_desired_primary_series_inductance();
    j["desiredPrimaryResonantCapacitance"] = x.get_desired_primary_resonant_capacitance();
}

} // namespace OpenMagnetics
