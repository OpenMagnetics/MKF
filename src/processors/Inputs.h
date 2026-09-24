#pragma once

#include "json.hpp"
#include "constructive_models/MasMigration.h"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "Definitions.h"
#include "support/Exceptions.h"
using json = nlohmann::json;

using namespace MAS;

namespace OpenMagnetics {

class Inputs : public MAS::Inputs {
  public:
    Inputs(json j, bool processWaveform = true, std::optional<std::variant<double, std::vector<double>>> magnetizingInductance = std::nullopt) {
        OpenMagnetics::compat::migrate_pre_1_0(j);
        throw_if_json_lacks_required_fields(j);
        from_json(j, *this);
        auto check_passed = check_integrity();
        if (!check_passed.first) {
            throw std::runtime_error("Missing inputs");
        }
        if (processWaveform) {
            process(magnetizingInductance);
        }
    }
    Inputs(std::vector<OperatingPoint> operatingPoints, DesignRequirements designRequirements, bool processWaveform = true, std::optional<std::variant<double, std::vector<double>>> magnetizingInductance = std::nullopt) {
        set_operating_points(operatingPoints);
        set_design_requirements(designRequirements);
        auto check_passed = check_integrity();
        if (!check_passed.first) {
            throw std::runtime_error("Missing inputs: " + check_passed.second);
        }
        if (processWaveform) {
            process(magnetizingInductance);
        }
    }
    Inputs() = default;
    virtual ~Inputs() = default;

    // ABT #1329: MAS requires inputs.designRequirements and inputs.operatingPoints. Without this
    // check a missing one surfaced as a raw "[json.exception.out_of_range.403] key
    // 'designRequirements' not found" from the generated parser, naming no MAS field.
    static void throw_if_json_lacks_required_fields(const json& inputsJson);
    // ABT #1329: insulation coordination (clearance, creepage, distance through insulation, lead
    // sleeves) is defined by the standards the design names. designRequirements.insulation.standards
    // is optional in MAS, so an insulation block without it is valid and simply names no standard
    // to coordinate by: true only when the insulation block carries a standards list. (An EMPTY
    // list, which MKF's own quick inputs build, keeps its previous meaning: coordinate by no table.)
    bool has_insulation_coordination_requirements() const;
    std::pair<bool, std::string> check_integrity();
    void process(std::optional<std::variant<double, std::vector<double>>> magnetizingInductance = std::nullopt);
    static OperatingPoint process_operating_point(OperatingPoint operatingPoint, double magnetizingInductance, std::optional<std::vector<double>> turnsRatios = std::nullopt, bool isDmcTopology = false);

    static bool is_waveform_sampled(Waveform waveform);
    static bool is_waveform_imported(Waveform waveform);
    static bool is_multiport_inductor(OperatingPoint operatingPoint, std::optional<std::vector<IsolationSide>> isolationSides = std::nullopt);
    static bool can_be_common_mode_choke(OperatingPoint operatingPoint);
    static bool can_be_differential_mode_choke(OperatingPoint operatingPoint, bool isDmcTopology = false);

    static Waveform calculate_sampled_waveform(Waveform waveform, double frequency=0, std::optional<size_t> numberPoints=std::nullopt, std::optional<size_t> maximumNumberPoints=std::nullopt);
    static ProcessedWaveform calculate_processed_data(Waveform waveform, std::optional<double> frequency=std::nullopt, bool includeAdvancedData=true, std::optional<ProcessedWaveform> processed=std::nullopt);
    static ProcessedWaveform calculate_processed_data(SignalDescriptor excitation, Waveform sampledWaveform, bool includeAdvancedData=true, std::optional<ProcessedWaveform> processed=std::nullopt);
    static ProcessedWaveform calculate_processed_data(Harmonics harmonics, Waveform waveform, bool includeAdvancedData=true, std::optional<ProcessedWaveform> processed=std::nullopt);
    static Harmonics calculate_harmonics_data(Waveform waveform, double frequency);
    static OperatingPointExcitation prune_harmonics(OperatingPointExcitation excitation, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);
    static SignalDescriptor prune_harmonics(SignalDescriptor signalDescriptor, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);
    static OperatingPoint prune_harmonics(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);

    // static OperatingPointExcitation reflect_waveforms(OperatingPointExcitation excitation, double ratio);
    static SignalDescriptor reflect_waveform(SignalDescriptor excitation, double ratio);
    static SignalDescriptor reflect_waveform(SignalDescriptor signal, double ratio, WaveformLabel label);

    static bool is_standardized(SignalDescriptor signal);
    static SignalDescriptor standardize_waveform(SignalDescriptor parameter, double frequency);
    static Waveform reconstruct_signal(Harmonics harmonics, double frequency);
    OperatingPoint get_operating_point(size_t index);
    const OperatingPoint& get_operating_point_ref(size_t index) const { return get_operating_points()[index]; }
    OperatingPointExcitation get_winding_excitation(size_t operatingPointIndex=0, size_t windingIndex=0);
    OperatingPointExcitation get_primary_excitation(size_t operatingPointIndex=0);
    // Perf (Phase 6): const-ref overload avoids deep-copying the OperatingPoint
    // just to read a single excitation. Hot in AreaProduct and Loss filters.
    static OperatingPointExcitation get_primary_excitation(const OperatingPoint& operatingPoint);

    static SignalDescriptor calculate_induced_voltage(OperatingPointExcitation& excitation, double magnetizingInductance, bool compress=true);

    /**
     * @brief Peak excursion of ∫V dt across the primary voltage waveform.
     *
     * Returns max|∫V dt| over the waveform's samples. This is the V·s
     * quantity that pairs with Faraday's law (`B_peak = V·s / (N · A_e)`)
     * for non-sinusoidal voltages — DAB squares, CLLLC trapezoids,
     * push-pull rectangulars, etc. — where the sinusoidal V_peak/ω
     * shortcut underestimates B_peak by ~36 %.
     *
     * Fallback: if the operating point doesn't carry a sampled voltage
     * waveform but has a processed peak + frequency, return V_peak/ω
     * (the sinusoidal estimate). Returns 0 if neither is available.
     *
     * Iterates over ALL operating points and ALL excitations (or the
     * single OP / primary excitation overload below), returning the
     * worst-case excursion.
     */
    static double calculate_max_volt_seconds(const OperatingPointExcitation& excitation);
    static double calculate_max_volt_seconds(const OperatingPoint& operatingPoint);
    static bool include_dc_offset_into_magnetizing_current(OperatingPoint operatingPoint, std::vector<double> turnsRatios);
    // addOffset: anchor the magnetizing current's DC level to the winding current
    // (see include_dc_offset_into_magnetizing_current). alternatingConduction: true
    // when the operating point has several windings conducting alternately (flyback
    // style); the winding current then contains commutation steps, so the DC anchor
    // is peak-based on the volt-second integral instead of the measured midpoint
    // (ABT #907). Pass excitations_per_winding().size() > 1 at the call site.
    // The ampere-turn current of an operating point: i_m(t) = sum_k c_k (N_k/N_ref) i_k(t), the
    // single current that, in the reference winding, drives the same core flux as all windings
    // together. Built from the winding current WAVEFORMS on one common period grid, so phase and
    // dot direction count: c_k = +1 for a primary-side winding, -1 otherwise (MAS's excitation
    // sign convention: primary-side currents passive, the others source), N_k/N_ref from the
    // design's turnsRatios (Np/Nk against winding 0, re-referred to the reference) or 1:1 when
    // none are given; the reference is the first primary-side winding. Its processed values come
    // from the summed waveform (offset = signed DC, peak = magnitude peak). Throws when a winding
    // has no current waveform (processed values alone carry no phase), when the windings'
    // frequencies differ, or when the isolation sides are not given for every winding.
    static SignalDescriptor calculate_ampere_turn_current(const OperatingPoint& operatingPoint,
                                                          const std::vector<double>& turnsRatios,
                                                          const std::vector<IsolationSide>& isolationSides);

    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            double magnetizingInductance,
                                                            bool compress,
                                                            bool addOffset,
                                                            bool alternatingConduction);
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            Waveform voltageSampledWaveform,
                                                            double magnetizingInductance,
                                                            bool compress,
                                                            bool addOffset,
                                                            bool alternatingConduction);
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            double magnetizingInductance,
                                                            bool compress,
                                                            double dcCurrent=0);
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            Waveform voltageSampledWaveform,
                                                            double magnetizingInductance,
                                                            bool compress,
                                                            double dcCurrent=0);
    static SignalDescriptor add_offset_to_excitation(SignalDescriptor signalDescriptor,
                                                             double offset,
                                                             double frequency);
    static OperatingPointExcitation get_excitation_with_proportional_current(OperatingPointExcitation excitation, double proportion);
    static OperatingPointExcitation get_excitation_with_proportional_voltage(OperatingPointExcitation excitation, double proportion);

    static void make_waveform_size_power_of_two(OperatingPoint* operatingPoint);

    static Inputs create_quick_operating_point(double frequency,
                                                      double magnetizingInductance,
                                                      double temperature,
                                                      WaveformLabel waveShape,
                                                      double peakToPeak,
                                                      double dutyCycle,
                                                      double dcCurrent,
                                                      std::vector<double> turnsRatios={});
    static Inputs create_quick_operating_point_only_current(double frequency,
                                                                  double magnetizingInductance,
                                                                  double temperature,
                                                                  WaveformLabel waveShape,
                                                                  double peakToPeak,
                                                                  double dutyCycle,
                                                                  double dcCurrent,
                                                                  std::vector<double> turnsRatios={});
    static Inputs create_quick_operating_point_only_current(double frequency,
                                                                  double magnetizingInductance,
                                                                  double temperature,
                                                                  WaveformLabel waveShape,
                                                                  std::vector<double> peakToPeaks,
                                                                  double dutyCycle,
                                                                  double dcCurrent,
                                                                  std::vector<double> turnsRatios={});
    static OperatingPoint create_operating_point_with_sinusoidal_current_mask(double frequency,
                                                                              double magnetizingInductance,
                                                                              double temperature,
                                                                              std::vector<double> turnsRatios,
                                                                              std::vector<double> currentPeakMask,
                                                                              double currentOffset = 0);

    static WaveformLabel try_guess_waveform_label(Waveform waveform);
    // ABT #1330: whether a processed description defines a waveform on its own. A signal given only
    // by its processed parameters (label, peakToPeak, offset, dutyCycle...) is rebuilt with
    // create_waveform, which reads the duty cycle for every non-sinusoidal label and the dead time
    // for the *_WITH_DEADTIME ones; CUSTOM and RECTANGULAR_DCM have no parametric form at all.
    // Throws naming what is missing (`what` says which signal), never assumes a value.
    static void throw_if_processed_cannot_define_waveform(const ProcessedWaveform& processed, const std::string& what);
    static Waveform create_waveform(ProcessedWaveform processed, double frequency);
    // phase: radians, applied to SINUSOIDAL waveforms only (positive phase
    // advances the waveform — e.g. phase=pi/2 turns sin(wt) into sin(wt+pi/2)=cos(wt),
    // so the resulting signal leads the same-frequency sine by 90 degrees).
    // Ignored for non-sinusoidal labels (use `skew` for those).
    static Waveform create_waveform(WaveformLabel label, double peakToPeak, double frequency, double dutyCycle=0.5, double offset=0, double deadTime=0, double skew=0, double phase=0);
    static ProcessedWaveform calculate_basic_processed_data(Waveform waveform);
    static Waveform compress_waveform(const Waveform& waveform);

    static Waveform calculate_derivative_waveform(Waveform waveform);
    static Waveform calculate_integral_waveform(Waveform waveform, bool subtractAverage=false);
    static double calculate_waveform_average(Waveform waveform) ;
    static Waveform multiply_waveform(Waveform waveform, double scalarValue);
    static Waveform sum_waveform(Waveform waveform, double scalarValue);
    // Pointwise add two waveforms. Both must have the same number of data points.
    // Time vector is taken from `a` (must be aligned with `b`).
    static Waveform add_waveforms(Waveform a, const Waveform& b);
    // Pointwise subtract: returns a - b. Same length requirement.
    static Waveform subtract_waveforms(Waveform a, const Waveform& b);
    // Circularly shift waveform samples by `samples` positions to the right
    // (i.e. delay by samples * dt). Negative values shift left (advance).
    // Operates on the data vector only; time vector (if any) is preserved.
    static Waveform shift_waveform(Waveform waveform, int samples);
    // SignalDescriptor wrappers: perform the op on the underlying sampled
    // waveform and recompute harmonics and processed data at `frequency`.
    static SignalDescriptor add_signals(const SignalDescriptor& a, const SignalDescriptor& b, double frequency);
    static SignalDescriptor subtract_signals(const SignalDescriptor& a, const SignalDescriptor& b, double frequency);
    static SignalDescriptor shift_signal(const SignalDescriptor& s, int samples, double frequency);

    static double try_guess_duty_cycle(Waveform waveform, WaveformLabel label=WaveformLabel::CUSTOM, double frequency=0);
    static double calculate_instantaneous_power(OperatingPointExcitation excitation);

    static double calculate_waveform_coefficient(OperatingPoint* operatingPoint);

    static void scale_time_to_frequency(Inputs& inputs, double newFrequency, bool cleanFrequencyDependentFields=false, bool processSignals=false, bool useCurrentAsBase=true);
    static void scale_time_to_frequency(OperatingPoint& operatingPoint, double newFrequency, bool cleanFrequencyDependentFields=false, bool processSignals=false, bool useCurrentAsBase=true);
    static void scale_time_to_frequency(OperatingPointExcitation& excitation, double newFrequency, bool cleanFrequencyDependentFields=false, bool processSignals=false, bool useCurrentAsBase=true);
    static Waveform scale_time_to_frequency(Waveform waveform, double newFrequency);
    
    void set_operating_point_by_index(const OperatingPoint& value, size_t index) {
        get_mutable_operating_points()[index] = value;
    }
    static void set_current_as_magnetizing_current(OperatingPoint* operatingPoint);

    static double get_switching_frequency(const OperatingPointExcitation& excitation);
    static double get_magnetic_flux_density_peak(OperatingPointExcitation excitation, double switchingFrequency);
    static double get_magnetic_flux_density_peak_to_peak(OperatingPointExcitation excitation, double switchingFrequency);
    static SignalDescriptor get_multiport_inductor_magnetizing_current(OperatingPoint operatingPoint);
    static SignalDescriptor get_common_mode_choke_magnetizing_current(OperatingPoint operatingPoint);
    static SignalDescriptor get_differential_mode_choke_magnetizing_current(OperatingPoint operatingPoint);

    DimensionWithTolerance get_altitude();
    Cti get_cti();
    IsolationClass get_insulation_type();
    DimensionWithTolerance get_main_supply_voltage();
    OvervoltageCategory get_overvoltage_category();
    PollutionDegree get_pollution_degree();
    std::vector<InsulationStandards> get_standards();
    WiringTechnology get_wiring_technology() const;
    double get_maximum_voltage_peak();
    double get_maximum_voltage_rms();
    double get_maximum_current_effective_frequency();
    double get_maximum_current_peak();
    double get_maximum_current_rms();
    double get_maximum_frequency();
    double get_maximum_temperature() const;
    double get_maximum_voltage_peak(size_t windingIndex);
    double get_maximum_voltage_rms(size_t windingIndex);
    double get_maximum_current_peak(size_t windingIndex);
    double get_maximum_current_rms(size_t windingIndex);
    double get_maximum_current_effective_frequency(size_t windingIndex);
    double get_maximum_current_dc_bias(size_t windingIndex = 0);
    SignalDescriptor get_current_with_effective_maximum(size_t windingIndex = 0);
    std::vector<IsolationSide> get_isolation_sides_used();

};
void from_json(const json& j, Inputs& x);
void to_json(json& j, const Inputs& x);
void to_file(std::filesystem::path filepath, const Inputs & x);

inline void from_json(const json& j, Inputs& x) {
    x.set_design_requirements(j.at("designRequirements").get<DesignRequirements>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<OperatingPoint>>());
}
inline void to_json(json& j, const Inputs& x) {
    j = json::object();
    j["designRequirements"] = x.get_design_requirements();
    j["operatingPoints"] = x.get_operating_points();
}
inline void to_file(std::filesystem::path filepath, const Inputs & x) {
    json masJson;
    to_json(masJson, x);

    std::ofstream myfile;
    myfile.open(filepath);
    myfile << masJson;
}


bool operator==(Inputs lhs, Inputs rhs);

inline bool operator==(Inputs lhs, Inputs rhs) {
    // TODO: Add more comparisond
    bool isEqual = resolve_dimensional_values(lhs.get_design_requirements().get_magnetizing_inductance()) == resolve_dimensional_values(rhs.get_design_requirements().get_magnetizing_inductance()) &&
                   lhs.get_operating_points().size() == rhs.get_operating_points().size();

    if (isEqual) {
        // FIXED: BUG-09 — Check sizes match
        if (lhs.get_design_requirements().get_turns_ratios().size() != rhs.get_design_requirements().get_turns_ratios().size()) return false;
        for (size_t i = 0; i < lhs.get_design_requirements().get_turns_ratios().size(); ++i) {
            isEqual &= resolve_dimensional_values(lhs.get_design_requirements().get_turns_ratios()[i]) == resolve_dimensional_values(rhs.get_design_requirements().get_turns_ratios()[i]);
        }
        for (size_t i = 0; i < lhs.get_operating_points().size(); ++i) {
            isEqual &= lhs.get_operating_points()[i].get_excitations_per_winding().size() == rhs.get_operating_points()[i].get_excitations_per_winding().size();
        }
    }

    return isEqual;
}

} // namespace OpenMagnetics
