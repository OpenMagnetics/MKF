#pragma once

#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
#include "Definitions.h"
using json = nlohmann::json;

using namespace MAS;

namespace OpenMagnetics {

class Inputs : public MAS::Inputs {
  public:
    Inputs(const json& j, bool processWaveform = true, std::optional<std::variant<double, std::vector<double>>> magnetizingInductance = std::nullopt) {
        from_json(j, *this);
        auto check_passed = check_integrity();
        if (!check_passed.first) {
            throw std::runtime_error("Missing inputs");
        }
        if (processWaveform) {
            process(magnetizingInductance);
        }
    }
    Inputs() = default;
    virtual ~Inputs() = default;

    std::pair<bool, std::string> check_integrity();
    void process(std::optional<std::variant<double, std::vector<double>>> magnetizingInductance = std::nullopt);
    static OperatingPoint process_operating_point(OperatingPoint operatingPoint, double magnetizingInductance, std::optional<std::vector<double>> turnsRatios = std::nullopt);

    static bool is_waveform_sampled(Waveform waveform);
    static bool is_waveform_imported(Waveform waveform);
    static bool is_multiport_inductor(OperatingPoint operatingPoint, std::optional<std::vector<IsolationSide>> isolationSides = std::nullopt);
    static bool can_be_common_mode_choke(OperatingPoint operatingPoint);

    static Waveform calculate_sampled_waveform(Waveform waveform, double frequency=0, std::optional<size_t> numberPoints=std::nullopt);
    static Processed calculate_processed_data(Waveform waveform, std::optional<double> frequency=std::nullopt, bool includeAdvancedData=true, std::optional<Processed> processed=std::nullopt);
    static Processed calculate_processed_data(SignalDescriptor excitation, Waveform sampledWaveform, bool includeAdvancedData=true, std::optional<Processed> processed=std::nullopt);
    static Processed calculate_processed_data(Harmonics harmonics, Waveform waveform, bool includeAdvancedData=true, std::optional<Processed> processed=std::nullopt);
    static Harmonics calculate_harmonics_data(Waveform waveform, double frequency);
    static OperatingPointExcitation prune_harmonics(OperatingPointExcitation excitation, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);
    static SignalDescriptor prune_harmonics(SignalDescriptor signalDescriptor, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);

    // static OperatingPointExcitation reflect_waveforms(OperatingPointExcitation excitation, double ratio);
    static SignalDescriptor reflect_waveform(SignalDescriptor excitation, double ratio);
    static SignalDescriptor reflect_waveform(SignalDescriptor signal, double ratio, WaveformLabel label);

    static bool is_standardized(SignalDescriptor signal);
    static SignalDescriptor standardize_waveform(SignalDescriptor parameter, double frequency);
    static Waveform reconstruct_signal(Harmonics harmonics, double frequency);
    OperatingPoint get_operating_point(size_t index);
    OperatingPointExcitation get_winding_excitation(size_t operatingPointIndex=0, size_t windingIndex=0);
    OperatingPointExcitation get_primary_excitation(size_t operatingPointIndex=0);
    static OperatingPointExcitation get_primary_excitation(OperatingPoint operatingPoint);

    static SignalDescriptor calculate_induced_voltage(OperatingPointExcitation& excitation, double magnetizingInductance, bool compress=true);
    static bool include_dc_offset_into_magnetizing_current(OperatingPoint operatingPoint, std::vector<double> turnsRatios);
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            double magnetizingInductance,
                                                            bool compress,
                                                            bool addOffset);
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            Waveform voltageSampledWaveform,
                                                            double magnetizingInductance,
                                                            bool compress,
                                                            bool addOffset);
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
    static Waveform create_waveform(Processed processed, double frequency);
    static Waveform create_waveform(WaveformLabel label, double peakToPeak, double frequency, double dutyCycle=0.5, double offset=0, double deadTime=0, double skew=0);
    static Processed calculate_basic_processed_data(Waveform waveform);
    static Waveform compress_waveform(Waveform waveform);

    static Waveform calculate_derivative_waveform(Waveform waveform);
    static Waveform calculate_integral_waveform(Waveform waveform, bool subtractAverage=false);
    static double calculate_waveform_average(Waveform waveform) ;
    static Waveform multiply_waveform(Waveform waveform, double scalarValue);
    static Waveform sum_waveform(Waveform waveform, double scalarValue);

    static double try_guess_duty_cycle(Waveform waveform, WaveformLabel label=WaveformLabel::CUSTOM);
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

    static double get_switching_frequency(OperatingPointExcitation excitation);
    static double get_magnetic_flux_density_peak(OperatingPointExcitation excitation, double switchingFrequency);
    static double get_magnetic_flux_density_peak_to_peak(OperatingPointExcitation excitation, double switchingFrequency);
    static SignalDescriptor get_multiport_inductor_magnetizing_current(OperatingPoint operatingPoint);
    static SignalDescriptor get_common_mode_choke_magnetizing_current(OperatingPoint operatingPoint);

    DimensionWithTolerance get_altitude();
    Cti get_cti();
    InsulationType get_insulation_type();
    DimensionWithTolerance get_main_supply_voltage();
    OvervoltageCategory get_overvoltage_category();
    PollutionDegree get_pollution_degree();
    std::vector<InsulationStandards> get_standards();
    WiringTechnology get_wiring_technology();
    double get_maximum_voltage_peak();
    double get_maximum_voltage_rms();
    double get_maximum_current_effective_frequency();
    double get_maximum_current_peak();
    double get_maximum_current_rms();
    double get_maximum_frequency();
    double get_maximum_temperature();
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
