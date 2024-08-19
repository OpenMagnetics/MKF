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
    static bool is_waveform_imported(Waveform waveform);
    static bool is_multiport_inductor(OperatingPoint operatingPoint);
    static Waveform calculate_sampled_waveform(Waveform waveform, double frequency = 0);
    static Processed calculate_processed_data(Waveform waveform, std::optional<double> frequency = std::nullopt, bool includeAdvancedData = true, std::optional<Processed> processed = std::nullopt);
    static Processed calculate_processed_data(SignalDescriptor excitation, Waveform sampledWaveform, bool includeAdvancedData = true, std::optional<Processed> processed = std::nullopt);
    static Processed calculate_processed_data(Harmonics harmonics, Waveform waveform, bool includeAdvancedData = true, std::optional<Processed> processed = std::nullopt);
    static Harmonics calculate_harmonics_data(Waveform waveform, double frequency);
    static SignalDescriptor reflect_waveform(SignalDescriptor excitation, double ratio);
    static SignalDescriptor reflect_waveform(SignalDescriptor signal, double ratio, WaveformLabel label);

    static SignalDescriptor standarize_waveform(SignalDescriptor parameter, double frequency);
    static Waveform reconstruct_signal(Harmonics harmonics, double frequency);
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
    static void set_current_as_magnetizing_current(OperatingPoint* operatingPoint);

    static double get_switching_frequency(OperatingPointExcitation excitation);
    static double get_magnetic_flux_density_peak(OperatingPointExcitation excitation, double switchingFrequency);
    static double get_magnetic_flux_density_peak_to_peak(OperatingPointExcitation excitation, double switchingFrequency);
    static SignalDescriptor get_multiport_inductor_magnetizing_current(OperatingPoint operatingPoint);

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
    double get_maximum_temperature();
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


    class CircuitSimulationReader {
      public:
        enum class DataType : int {
            TIME,
            VOLTAGE,
            CURRENT
        };

      private:
        class CircuitSimulationSignal {
            public:
                CircuitSimulationSignal() = default;
                virtual ~CircuitSimulationSignal() = default;
                std::string name;
                std::vector<double> data;
                DataType type;
                int windingIndex;
                int operatingPointIndex;

            public:

        };


        std::vector<CircuitSimulationSignal> _columns;
        std::vector<Waveform> _waveforms;
        CircuitSimulationSignal _time;

        std::vector<std::string> _timeAliases = {"TIME", "Time", "time", "[s]"};
        std::vector<std::string> _currentAliases = {"CURRENT", "CURR", "Current", "Curr", "I(", "current", "curr", "i(", "[A]", "Ip", "Is", "It"};
        std::vector<std::string> _voltageAliases = {"VOLTAGE", "VOLT", "Voltage", "Volt", "V(", "voltage", "volt", "v(", "[V]", "Vp", "Vs", "Vt"};

      public:

        CircuitSimulationReader() = default;
        virtual ~CircuitSimulationReader() = default;

        CircuitSimulationReader(std::string filePathOrFile);

        void process_line(std::string line, char separator);
        bool extract_winding_indexes(size_t numberWindings);
        bool extract_column_types(double frequency);
        OperatingPoint extract_operating_point(size_t numberWindings, double frequency, std::optional<std::vector<std::map<std::string, std::string>>> mapColumnNames = std::nullopt, double ambientTemperature = 25);
        std::vector<std::map<std::string, std::string>> extract_map_column_names(size_t numberWindings, double frequency);
        bool assign_column_names(std::vector<std::map<std::string, std::string>> columnNames);
        std::vector<std::string> extract_column_names();
        Waveform extract_waveform(CircuitSimulationSignal signal, double frequency, bool sample=true);
        static CircuitSimulationSignal find_time(std::vector<CircuitSimulationSignal> columns);
        static Waveform get_one_period(Waveform waveform, double frequency, bool sample=true);
        static char guess_separator(std::string line);
        static bool can_be_voltage(std::vector<double> data, double limit=0.05);
        static bool can_be_time(std::vector<double> data);
        static bool can_be_current(std::vector<double> data, double limit=0.05);
        std::optional<InputsWrapper::CircuitSimulationReader::DataType> guess_type_by_name(std::string name);
    };
};
} // namespace OpenMagnetics
