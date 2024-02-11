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
using json = nlohmann::json;

namespace OpenMagnetics {

class InputsWrapper : public Inputs {
  public:
    InputsWrapper(const json& j, bool processWaveform = true) {
        from_json(j, *this);
        auto check_passed = check_integrity();
        if (!check_passed.first) {
            throw std::runtime_error("Missing inputs");
        }
        if (processWaveform) {
            process_waveforms();
        }
    }
    InputsWrapper() = default;
    virtual ~InputsWrapper() = default;

    std::pair<bool, std::string> check_integrity();
    void process_waveforms();
    static OperatingPoint process_operating_point(OperatingPoint operatingPoint, double magnetizingInductance);

    static bool is_waveform_sampled(Waveform waveform);
    static Waveform calculate_sampled_waveform(Waveform waveform, double frequency = 0);
    static Processed calculate_processed_data(SignalDescriptor excitation, Waveform sampledWaveform, bool includeAdvancedData = true, std::optional<Processed> processed = std::nullopt);
    static Processed calculate_processed_data(Harmonics harmonics, Waveform waveform, bool includeAdvancedData = true, std::optional<Processed> processed = std::nullopt);
    static Harmonics calculate_harmonics_data(Waveform waveform, double frequency);
    static SignalDescriptor reflect_waveform(SignalDescriptor excitation, double ratio);
    static SignalDescriptor reflect_waveform(SignalDescriptor signal, double ratio, WaveformLabel label);

    static SignalDescriptor standarize_waveform(SignalDescriptor parameter, double frequency);
    OperatingPoint get_operating_point(size_t index);
    OperatingPointExcitation get_winding_excitation(size_t operatingPointIndex = 0, size_t windingIndex = 0);
    OperatingPointExcitation get_primary_excitation(size_t operatingPointIndex = 0);
    static OperatingPointExcitation get_primary_excitation(OperatingPoint operatingPoint);

    static SignalDescriptor calculate_induced_voltage(OperatingPointExcitation& excitation,
                                                        double magnetizingInductance);
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            double magnetizingInductance,
                                                            bool compress = true,
                                                            double offset=0);
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            Waveform sampledWaveform,
                                                            double magnetizingInductance,
                                                            bool compress = true,
                                                            double offset=0);
    static SignalDescriptor add_offset_to_excitation(SignalDescriptor signalDescriptor,
                                                             double offset,
                                                             double frequency);
    static void make_waveform_size_power_of_two(OperatingPoint* operatingPoint);

    static InputsWrapper create_quick_operating_point(double frequency,
                                                      double magnetizingInductance,
                                                      double temperature,
                                                      WaveformLabel waveShape,
                                                      double peakToPeak,
                                                      double dutyCycle,
                                                      double dcCurrent,
                                                      std::vector<double> turnsRatios = {});
    static InputsWrapper create_quick_operating_point_only_current(double frequency,
                                                                  double magnetizingInductance,
                                                                  double temperature,
                                                                  WaveformLabel waveShape,
                                                                  double peakToPeak,
                                                                  double dutyCycle,
                                                                  double dcCurrent,
                                                                  std::vector<double> turnsRatios = {});
    static InputsWrapper create_quick_operating_point_only_current(double frequency,
                                                                  double magnetizingInductance,
                                                                  double temperature,
                                                                  WaveformLabel waveShape,
                                                                  std::vector<double> peakToPeaks,
                                                                  double dutyCycle,
                                                                  double dcCurrent,
                                                                  std::vector<double> turnsRatios = {});

    static WaveformLabel try_guess_waveform_label(Waveform waveform);
    static Waveform create_waveform(Processed processed, double frequency);
    static Processed calculate_basic_processed_data(Waveform waveform);
    static Waveform compress_waveform(Waveform waveform);

    static Waveform calculate_derivative_waveform(Waveform waveform);
    static Waveform calculate_integral_waveform(Waveform waveform);

    static double try_guess_duty_cycle(Waveform waveform, WaveformLabel label = WaveformLabel::CUSTOM);
    static double calculate_instantaneous_power(OperatingPointExcitation excitation);

    static double calculate_waveform_coefficient(OperatingPoint* operatingPoint);

    static void scale_time_to_frequency(InputsWrapper& inputs, double newFrequency);
    static void scale_time_to_frequency(OperatingPoint& operatingPoint, double newFrequency);
    static void scale_time_to_frequency(OperatingPointExcitation& excitation, double newFrequency);
    static Waveform scale_time_to_frequency(Waveform waveform, double newFrequency);
    
    void set_operating_point_by_index(const OperatingPoint& value, size_t index) {
        get_mutable_operating_points()[index] = value;
    }

    DimensionWithTolerance get_altitude();
    Cti get_cti();
    InsulationType get_insulation_type();
    DimensionWithTolerance get_main_supply_voltage();
    OvervoltageCategory get_overvoltage_category();
    PollutionDegree get_pollution_degree();
    std::vector<InsulationStandards> get_standards();
    double get_maximum_voltage_peak();
    double get_maximum_voltage_rms();
    double get_maximum_current_effective_frequency();
    double get_maximum_current_peak();
    double get_maximum_current_rms();
    double get_maximum_frequency();
    double get_maximum_voltage_peak(size_t windingIndex);
    double get_maximum_voltage_rms(size_t windingIndex);
    double get_maximum_current_peak(size_t windingIndex);
    double get_maximum_current_rms(size_t windingIndex);

    void from_json(const json& j, Inputs& x);
    void to_json(json& j, const Inputs& x);

    inline void from_json(const json& j, InputsWrapper& x) {
        x.set_design_requirements(j.at("designRequirements").get<DesignRequirements>());
        x.set_operating_points(j.at("operatingPoints").get<std::vector<OperatingPoint>>());
    }
    inline void to_json(json& j, const InputsWrapper& x) {
        j = json::object();
        j["designRequirements"] = x.get_design_requirements();
        j["operatingPoints"] = x.get_operating_points();
    }
};
} // namespace OpenMagnetics
