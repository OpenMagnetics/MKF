#pragma once
#include "json.hpp"
#include "MAS.hpp"
#include "MagneticWrapper.h"
#include <random>
#include <chrono>

using ordered_json = nlohmann::ordered_json;

namespace OpenMagnetics {


enum class CircuitSimulatorExporterModels : int {
    SIMBA,
    NGSPICE
};


class CircuitSimulatorExporterModel {
    private:
    public:
        std::string programName = "Default";
        virtual ordered_json export_magnetic_as_subcircuit(MagneticWrapper magnetic, std::optional<std::string> filePathOrFile = std::nullopt) = 0;
        static std::shared_ptr<CircuitSimulatorExporterModel> factory(CircuitSimulatorExporterModels programName);

};

class CircuitSimulatorExporter {
    private:
    protected:
        std::shared_ptr<CircuitSimulatorExporterModel> _model;
    public:
        CircuitSimulatorExporter(CircuitSimulatorExporterModels program);
        CircuitSimulatorExporter();
        void export_magnetic_as_subcircuit(MagneticWrapper magnetic, std::string outputFilename, std::optional<std::string> filePathOrFile = std::nullopt);
};

class CircuitSimulatorExporterSimbaModel : public CircuitSimulatorExporterModel {
    enum class SimbaSupportedDeviceTypes : int {
        AIR_GAP,
        LINEAR_CORE,
        WINDING,
        ELECTRICAL_PIN
    };


    public:
        std::string programName = "Simba";
        double _scale;
        double _modelSize = 40;
        std::string generate_id();
        std::mt19937 _gen;
        CircuitSimulatorExporterSimbaModel();
        ordered_json export_magnetic_as_subcircuit(MagneticWrapper magnetic, std::optional<std::string> filePathOrFile = std::nullopt);
        ordered_json create_device(std::string libraryName, std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_air_gap(std::vector<int> coordinates, double area, double length, int angle, std::string name);
        ordered_json create_core(double initialPermeability, std::vector<int> coordinates, double area, double length, int angle, std::string name);
        ordered_json create_winding(size_t numberTurns, std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_pin(std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_resistor(double resistance, std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_inductor(double inductance, std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_connector(std::vector<int> startingCoordinates, std::vector<int> endingCoordinates, std::string name);
        ordered_json merge_connectors(ordered_json connectors);
};

class CircuitSimulatorExporterNgspiceModel : public CircuitSimulatorExporterModel {
    public:
        std::string programName = "Ngspice";
        ordered_json export_magnetic_as_subcircuit(MagneticWrapper magnetic, std::optional<std::string> filePathOrFile = std::nullopt);
};

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
    std::optional<CircuitSimulationReader::DataType> guess_type_by_name(std::string name);
};


} // namespace OpenMagnetics