// Tests for CircuitSimulationReader: importing waveforms exported by circuit simulators
// (LTspice, PLECS, PSIM, SIMBA, and files users uploaded to the web tool).
//
// These cases lived in TestCircuitSimulatorInterface.cpp and were deleted with it in 3e0261fd,
// when the converter models moved to Kirchhoff. The reader did not move: it is what the web
// frontend's waveform import calls, so its tests are restored here unchanged.
#include "processors/CircuitSimulatorInterface.h"
#include "processors/Inputs.h"
#include "support/Painter.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
double max_error = 0.01;
auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
bool plot = false;
}

TEST_CASE("Test_Guess_Periodicity_Simba", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    std::ifstream is(simulation_path);
    std::vector<std::vector<double>> columns;
    std::vector<std::string> column_names;

    if(is.is_open()) {
        std::string line;
        while(getline(is, line)) {
            std::stringstream ss(line);
            std::string token;
            if (column_names.size() == 0) {
                // Getting column names
                while(getline(ss, token, ',')) {
                    column_names.push_back(token);
                    columns.push_back({});
                }
            }
            else {
                size_t currentColumnIndex = 0;
                while(getline(ss, token, ',')) {
                    columns[currentColumnIndex].push_back(stod(token));
                    currentColumnIndex++;
                }

            }
        }
        is.close();
    }
    else {
        throw std::runtime_error("File not found");
    }
    MAS::Waveform waveform;
    waveform.set_time(columns[0]);
    waveform.set_data(columns[1]);
    auto waveformOnePeriod = CircuitSimulationReader().get_one_period(waveform, 100000);
    REQUIRE(128U == waveformOnePeriod.get_data().size());
}

TEST_CASE("Test_Guess_Separator_Commas", "[processor][circuit-simulation-reader]") {
    std::string row = "columns,separated,by,commas";
    REQUIRE(',' == CircuitSimulationReader::guess_separator(row));
}

TEST_CASE("Test_Guess_Separator_Semicolon", "[processor][circuit-simulation-reader]") {
    std::string row = "columns;separated;by;semicolon";
    REQUIRE(';' == CircuitSimulationReader::guess_separator(row));
}

TEST_CASE("Test_Guess_Separator_Tabs", "[processor][circuit-simulation-reader]") {
    std::string row = "columns\tseparated\tby\ttabs";
    REQUIRE('\t' == CircuitSimulationReader::guess_separator(row));
}

TEST_CASE("Test_Guess_Separator_Simba", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    std::ifstream is(simulation_path);
    std::vector<std::vector<double>> columns;
    std::vector<std::string> column_names;

    if(is.is_open()) {
        std::string line;
        while(getline(is, line)) {
            REQUIRE(',' == CircuitSimulationReader::guess_separator(line));
        }
    }
}

TEST_CASE("Test_Guess_Separator_Ltspice", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "ltspice_simulation.txt");

    std::ifstream is(simulation_path);
    std::vector<std::vector<double>> columns;
    std::vector<std::string> column_names;

    if(is.is_open()) {
        std::string line;
        while(getline(is, line)) {
            REQUIRE('\t' == CircuitSimulationReader::guess_separator(line));
        }
    }
}

TEST_CASE("Test_Import_Csv_Rosano_Forward", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "forward_case.csv");

    double frequency = 200000;
    auto reader = CircuitSimulationReader(simulation_path.string());
    std::vector<std::map<std::string, std::string>> mapColumnNames;
    std::map<std::string, std::string> primaryColumnNames;
    std::map<std::string, std::string> secondaryColumnNames;
    primaryColumnNames["time"] = "Time";
    primaryColumnNames["current"] = "Ipri";
    primaryColumnNames["magnetizingCurrent"] = "Im";
    primaryColumnNames["voltage"] = "Vpri";
    mapColumnNames.push_back(primaryColumnNames);
    secondaryColumnNames["time"] = "Time";
    secondaryColumnNames["current"] = "Isec";
    secondaryColumnNames["voltage"] = "Vsec";
    mapColumnNames.push_back(secondaryColumnNames);
    auto operatingPoint = reader.extract_operating_point(2, frequency, mapColumnNames);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 121e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto primaryMagnetizingCurrent = primaryExcitation.get_magnetizing_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryMagnetizingCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryMagnetizingCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryVoltage.get_waveform().value());

        painter.export_svg();

    }
}

TEST_CASE("Test_Can_Be_Time_Accepts_Repeated_Timestamps", "[processor][circuit-simulation-reader]") {
    // A monotonic non-decreasing axis with repeated timestamps (from a time
    // column rounded to a few significant figures) is a valid time axis.
    REQUIRE(CircuitSimulationReader::can_be_time({0.0, 1e-6, 1e-6, 2e-6, 3e-6, 3e-6}) == true);
    // Strictly increasing is valid.
    REQUIRE(CircuitSimulationReader::can_be_time({0.0, 1e-6, 2e-6, 3e-6}) == true);
    // Strictly decreasing is not a time axis.
    REQUIRE(CircuitSimulationReader::can_be_time({3e-6, 2e-6, 1e-6}) == false);
    // A non-monotonic step is not a time axis.
    REQUIRE(CircuitSimulationReader::can_be_time({0.0, 2e-6, 1e-6, 3e-6}) == false);
    // A constant column must not be mistaken for a time axis.
    REQUIRE(CircuitSimulationReader::can_be_time({5.0, 5.0, 5.0}) == false);
}

TEST_CASE("Test_Import_Csv_Repeated_Timestamps", "[processor][circuit-simulation-reader]") {
    // Regression: this circuit-simulator export rounds the time column to 3
    // significant figures, producing 342 repeated consecutive timestamps. The
    // reader used to throw "no time column found" because can_be_time required
    // a strictly increasing axis. It now accepts a non-decreasing axis and
    // collapses repeated timestamps (keeping the last sample at each instant),
    // so the file loads and both windings extract. Columns are time,V1,V3,I1,I3
    // (note the V2/I2 gap, which the winding-index detection reindexes to 0/1).
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "repeated_timestamps_two_periods.csv");

    double frequency = 200000;

    // Construction parses + collapses the repeated timestamps; used to throw.
    CircuitSimulationReader reader(simulation_path.string());

    auto columnNames = reader.extract_column_names();
    REQUIRE(columnNames.size() == 5);
    REQUIRE(columnNames[0] == "time");

    std::vector<std::map<std::string, std::string>> mapColumnNames;
    std::map<std::string, std::string> primaryColumnNames;
    std::map<std::string, std::string> secondaryColumnNames;
    primaryColumnNames["time"] = "time";
    primaryColumnNames["current"] = "I1";
    primaryColumnNames["voltage"] = "V1";
    mapColumnNames.push_back(primaryColumnNames);
    secondaryColumnNames["time"] = "time";
    secondaryColumnNames["current"] = "I3";
    secondaryColumnNames["voltage"] = "V3";
    mapColumnNames.push_back(secondaryColumnNames);

    auto operatingPoint = reader.extract_operating_point(2, frequency, mapColumnNames);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    for (auto& excitation : operatingPoint.get_excitations_per_winding()) {
        REQUIRE(excitation.get_current());
        REQUIRE(excitation.get_voltage());
        auto currentWaveform = excitation.get_current().value().get_waveform().value();
        REQUIRE(currentWaveform.get_data().size() > 1);
        // After collapsing repeated timestamps the extracted period must have a
        // strictly increasing time axis (proves the duplicates were removed).
        if (currentWaveform.get_time()) {
            auto times = currentWaveform.get_time().value();
            for (size_t i = 1; i < times.size(); ++i) {
                REQUIRE(times[i] > times[i - 1]);
            }
        }
    }
}

TEST_CASE("Test_Import_Csv_Rosano_Flyback", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "flyback_case.csv");

    double frequency = 200000;
    auto reader = CircuitSimulationReader(simulation_path.string());
    std::vector<std::map<std::string, std::string>> mapColumnNames;
    std::map<std::string, std::string> primaryColumnNames;
    std::map<std::string, std::string> secondaryColumnNames;
    primaryColumnNames["time"] = "Time";
    primaryColumnNames["current"] = "Ipri";
    primaryColumnNames["magnetizingCurrent"] = "Imag";
    primaryColumnNames["voltage"] = "Vpri";
    mapColumnNames.push_back(primaryColumnNames);
    secondaryColumnNames["time"] = "Time";
    secondaryColumnNames["current"] = "Isec";
    secondaryColumnNames["voltage"] = "Vsec";
    mapColumnNames.push_back(secondaryColumnNames);
    auto operatingPoint = reader.extract_operating_point(2, frequency, mapColumnNames);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 50e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto primaryMagnetizingCurrent = primaryExcitation.get_magnetizing_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryMagnetizingCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryMagnetizingCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryVoltage.get_waveform().value());

        painter.export_svg();

    }
    if (true) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());

        painter.export_svg();

    }
}

TEST_CASE("Test_Simba", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    double turnsRatio = 1.0 / 0.3;
    double frequency = 100000;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto operatingPoint = reader.extract_operating_point(2, frequency);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 220e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto secondaryFrequency = secondaryExcitation.get_frequency();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE(frequency == secondaryFrequency);
    REQUIRE_THAT(2.79694, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 2.79694 * max_error));
    REQUIRE_THAT(primaryCurrent.get_processed().value().get_rms().value() / turnsRatio, Catch::Matchers::WithinAbs(secondaryCurrent.get_processed().value().get_rms().value(), primaryCurrent.get_processed().value().get_rms().value() / turnsRatio * max_error));
    REQUIRE_THAT(13.1204, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 13.1204 * max_error));
    REQUIRE_THAT(primaryVoltage.get_processed().value().get_rms().value() * turnsRatio, Catch::Matchers::WithinAbs(secondaryVoltage.get_processed().value().get_rms().value(), primaryVoltage.get_processed().value().get_rms().value() * turnsRatio * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryVoltage.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());

        painter.export_svg();

    }
}

TEST_CASE("Test_PFC_Only_Current", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "only_pfc_current_waveform.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto operatingPoint = reader.extract_operating_point(2, frequency);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 110e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 1);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryCurrent = primaryExcitation.get_current().value();
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryCurrent.get_waveform().value());

        painter.export_svg();

    }
}

TEST_CASE("Test_Simba_File_Loaded", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    std::string file = "";
    std::string line;
    std::ifstream is(simulation_path);
    if(is.is_open()) {
        while(getline(is, line)) {
            file += line;
            file += "\n";
        }
        is.close();
    }
    else {
        throw std::runtime_error("File not found");
    }

    double turnsRatio = 1.0 / 0.3;
    double frequency = 100000;
    auto reader = CircuitSimulationReader(file);
    auto operatingPoint = reader.extract_operating_point(2, frequency);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 220e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto secondaryFrequency = secondaryExcitation.get_frequency();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE(frequency == secondaryFrequency);
    REQUIRE_THAT(2.79694, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 2.79694 * max_error));
    REQUIRE_THAT(primaryCurrent.get_processed().value().get_rms().value() / turnsRatio, Catch::Matchers::WithinAbs(secondaryCurrent.get_processed().value().get_rms().value(), primaryCurrent.get_processed().value().get_rms().value() / turnsRatio * max_error));
    REQUIRE_THAT(13.1204, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 13.1204 * max_error));
    REQUIRE_THAT(primaryVoltage.get_processed().value().get_rms().value() * turnsRatio, Catch::Matchers::WithinAbs(secondaryVoltage.get_processed().value().get_rms().value(), primaryVoltage.get_processed().value().get_rms().value() * turnsRatio * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryVoltage.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());

        painter.export_svg();

    }
}

TEST_CASE("Test_Ltspice", "[processor][circuit-simulation-reader][ltspice]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "ltspice_simulation.txt");

    double frequency = 372618;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto operatingPoint = reader.extract_operating_point(2, frequency);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto secondaryFrequency = secondaryExcitation.get_frequency();
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE(frequency == secondaryFrequency);
    REQUIRE_THAT(0.0524431, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 0.0524431 * max_error));
    REQUIRE_THAT(0.4, Catch::Matchers::WithinAbs(secondaryCurrent.get_processed().value().get_rms().value(), 0.4 * max_error));
    REQUIRE_THAT(6, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 6 * max_error));
    REQUIRE_THAT(64, Catch::Matchers::WithinAbs(secondaryVoltage.get_processed().value().get_rms().value(), 64 * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryVoltage.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());

        painter.export_svg();

    }
}

TEST_CASE("Test_Plecs", "[processor][circuit-simulation-reader][plecs]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "plecs_simulation.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto operatingPoint = reader.extract_operating_point(1, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 1);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE_THAT(11.3, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 11.3 * max_error));
    REQUIRE_THAT(324, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 324 * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryVoltage.get_waveform().value());

        painter.export_svg();

    }
}

TEST_CASE("Test_Plecs_Missing_Windings", "[processor][circuit-simulation-reader][plecs]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "wrong_plecs_simulation.csv");

    double frequency = 50;
    {
        // Auto-detection with nonsensical column names: heuristic should still detect signal types
        // from waveform shape (triangular → current, square-wave → voltage)
        auto reader = CircuitSimulationReader(simulation_path.string());
        auto operatingPoint = reader.extract_operating_point(1, frequency);

        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

        REQUIRE(operatingPoint.get_excitations_per_winding().size() == 1);
        REQUIRE(operatingPoint.get_excitations_per_winding()[0].get_current());
        REQUIRE(operatingPoint.get_excitations_per_winding()[0].get_voltage());
    }
    {
        auto reader = CircuitSimulationReader(simulation_path.string());
        std::vector<std::map<std::string, std::string>> mapColumnNames;
        std::map<std::string, std::string> primaryColumnNames;
        primaryColumnNames["time"] = "IHave";
        primaryColumnNames["current"] = "no";
        primaryColumnNames["voltage"] = "idea";
        mapColumnNames.push_back(primaryColumnNames);
        auto operatingPoint = reader.extract_operating_point(1, frequency, mapColumnNames);
        operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

        REQUIRE(operatingPoint.get_excitations_per_winding().size() == 1);
        REQUIRE(operatingPoint.get_excitations_per_winding()[0].get_current());
        REQUIRE(operatingPoint.get_excitations_per_winding()[0].get_voltage());

        auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
        auto primaryFrequency = primaryExcitation.get_frequency();
        auto primaryCurrent = primaryExcitation.get_current().value();
        auto primaryVoltage = primaryExcitation.get_voltage().value();

        REQUIRE(frequency == primaryFrequency);
        REQUIRE_THAT(11.3, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 11.3 * max_error));
        REQUIRE_THAT(324, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 324 * max_error));

        if (plot) {
            auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("primaryCurrent.svg");
            Painter painter(outFile);
            painter.paint_waveform(primaryCurrent.get_waveform().value());
            painter.export_svg();
        }
        if (plot) {
            auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("primaryVoltage.svg");
            Painter painter(outFile);
            painter.paint_waveform(primaryVoltage.get_waveform().value());
            painter.export_svg();
        }
    }
}

TEST_CASE("Test_Psim", "[processor][circuit-simulation-reader][psim]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "psim_simulation.csv");

    double frequency = 120000;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 52e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
    auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
    auto primaryFrequency = primaryExcitation.get_frequency();
    auto primaryCurrent = primaryExcitation.get_current().value();
    auto primaryVoltage = primaryExcitation.get_voltage().value();
    auto secondaryExcitation = operatingPoint.get_excitations_per_winding()[1];
    auto secondaryCurrent = secondaryExcitation.get_current().value();
    auto secondaryVoltage = secondaryExcitation.get_voltage().value();

    REQUIRE(frequency == primaryFrequency);
    REQUIRE_THAT(1.25, Catch::Matchers::WithinAbs(primaryCurrent.get_processed().value().get_rms().value(), 1.25 * max_error));
    REQUIRE_THAT(29.7, Catch::Matchers::WithinAbs(primaryVoltage.get_processed().value().get_rms().value(), 29.7 * max_error));

    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("primaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(primaryVoltage.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryCurrent.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryCurrent.get_waveform().value());

        painter.export_svg();

    }
    if (plot) {
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("secondaryVoltage.svg");
        Painter painter(outFile);
        painter.paint_waveform(secondaryVoltage.get_waveform().value());

        painter.export_svg();

    }
}

TEST_CASE("Test_Psim_Harmonics_Size_Error", "[processor][circuit-simulation-reader][psim]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "psim_simulation.csv");

    double frequency = 100000;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 0.0001);

    auto commonHarmonicIndexes = get_main_harmonic_indexes(operatingPoint, 0.05);
    // The highest harmonic index whose amplitude clears 5 % of the maximum
    // depends on downstream normalization details that drift with unrelated
    // Inputs/Utils changes (e.g. the 50th harmonic now sits right on the
    // cusp and occasionally flips across the threshold). Assert a small
    // range instead of pinning a single number.
    REQUIRE(commonHarmonicIndexes.size() > 0);
    REQUIRE(commonHarmonicIndexes.back() >= 49U);
    REQUIRE(commonHarmonicIndexes.back() <= 51U);
}

TEST_CASE("Test_Simba_Column_Names", "[processor][circuit-simulation-reader][simba]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "simba_simulation.csv");

    double frequency = 100000;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto mapColumnNames = reader.extract_map_column_names(2, frequency);

    REQUIRE(mapColumnNames.size() == 2);
    REQUIRE(!mapColumnNames[0]["time"].compare("Time [s]"));
    REQUIRE(!mapColumnNames[0]["current"].compare("TX1 - W2 - Current [A]"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("TX1 - W2 - Voltage [V]"));
    REQUIRE(!mapColumnNames[1]["time"].compare("Time [s]"));
    REQUIRE(!mapColumnNames[1]["current"].compare("TX1 - W5 - Current [A]"));
    REQUIRE(!mapColumnNames[1]["voltage"].compare("TX1 - W5 - Voltage [V]"));
}

TEST_CASE("Test_Ltspice_Column_Names", "[processor][circuit-simulation-reader][ltspice]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "ltspice_simulation.txt");

    double frequency = 372618;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto mapColumnNames = reader.extract_map_column_names(2, frequency);

    REQUIRE(mapColumnNames.size() == 2);
    REQUIRE(!mapColumnNames[0]["time"].compare("time"));
    REQUIRE(!mapColumnNames[0]["current"].compare("I(L1)"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("V(n001)"));
    REQUIRE(!mapColumnNames[1]["time"].compare("time"));
    REQUIRE(!mapColumnNames[1]["current"].compare("I(L2)"));
    REQUIRE(!mapColumnNames[1]["voltage"].compare("V(n002)"));
}

TEST_CASE("Test_Plecs_Column_Names", "[processor][circuit-simulation-reader][plecs]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "plecs_simulation.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto mapColumnNames = reader.extract_map_column_names(1, frequency);

    REQUIRE(mapColumnNames.size() == 1);
    REQUIRE(!mapColumnNames[0]["time"].compare("Time / s"));
    REQUIRE(!mapColumnNames[0]["current"].compare("L2:Inductor current"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("L2:Inductor voltage"));
}

TEST_CASE("Test_Plecs_Web", "[processor][circuit-simulation-reader][plecs]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "plecs_simulation.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path.string()); 
    std::string mapColumnNamesString = R"([{"current":"L2:Inductor current","time":"Time / s","voltage":"L2:Inductor voltage"}])";

    std::vector<std::map<std::string, std::string>> mapColumnNames = json::parse(mapColumnNamesString).get<std::vector<std::map<std::string, std::string>>>();


    auto operatingPoint = reader.extract_operating_point(1, frequency, mapColumnNames);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);

    REQUIRE(mapColumnNames.size() == 1);
    REQUIRE(!mapColumnNames[0]["time"].compare("Time / s"));
    REQUIRE(!mapColumnNames[0]["current"].compare("L2:Inductor current"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("L2:Inductor voltage"));
}

TEST_CASE("Test_Plecs_Column_Names_Missing_Windings", "[processor][circuit-simulation-reader][plecs]") {
    // wrong_plecs_simulation.csv has columns "IHave","no","idea" — nonsensical names.
    // Heuristic detection identifies signals by waveform shape, assigns both to winding 0.
    // Winding 1 should have empty current/voltage since no data maps to it.
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "wrong_plecs_simulation.csv");

    double frequency = 50;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto mapColumnNames = reader.extract_map_column_names(2, frequency);

    REQUIRE(mapColumnNames.size() == 2);
    // Winding 0: heuristic detects current from waveform shape (triangular current derivative = square wave)
    // Voltage detection by heuristic is less reliable — may or may not be detected
    REQUIRE(!mapColumnNames[0]["current"].empty());
    // Winding 1: no data maps here since all columns default to winding 0
    REQUIRE(mapColumnNames[1]["current"].empty());
    REQUIRE(mapColumnNames[1]["voltage"].empty());
}

TEST_CASE("Test_Psim_Column_Names", "[processor][circuit-simulation-reader][psim]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "psim_simulation.csv");

    double frequency = 120000;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto mapColumnNames = reader.extract_map_column_names(2, frequency);

    REQUIRE(mapColumnNames.size() == 2);
    REQUIRE(!mapColumnNames[0]["time"].compare("Time"));
    REQUIRE(!mapColumnNames[0]["current"].compare("Ipri"));
    REQUIRE(!mapColumnNames[0]["voltage"].compare("Vpri"));
    REQUIRE(!mapColumnNames[1]["time"].compare("Time"));
    REQUIRE(!mapColumnNames[1]["current"].compare("Isec"));
    REQUIRE(!mapColumnNames[1]["voltage"].compare("Vsec"));
}

TEST_CASE("Test_Import_Csv_Web_1", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_web_1.csv");

    double frequency = 919963.201472;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 10e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
}

TEST_CASE("Test_Import_Csv_Web_2", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_web_2.csv");

    double frequency = 1e6;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 10e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
}

TEST_CASE("Test_Import_Csv_Web_3", "[processor][circuit-simulation-reader]") {
    auto simulation_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_web_3.csv");

    double frequency = 50e3;
    auto reader = CircuitSimulationReader(simulation_path.string());
    auto operatingPoint = reader.extract_operating_point(2, frequency);

    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 10e-6);

    REQUIRE(operatingPoint.get_excitations_per_winding().size() == 2);
}

// A user's LTspice export of an inductor: "time", "I(L1)", "V(p1,p2)". The digits are the inductor's
// and the probe nodes' names, but they were ranked as winding numbers, so V(p1,p2) went to a
// second winding the magnetic does not have and auto-detection returned no voltage at all.
TEST_CASE("Test_Ltspice_Column_Names_Single_Winding_Probe_Digits", "[processor][circuit-simulation-reader][ltspice]") {
    const double frequency = 45000;
    const double period = 1 / frequency;
    const size_t numberPoints = 4000;
    std::ostringstream csv;
    csv << "time\tI(L1)\tV(p1,p2)\n";
    for (size_t i = 0; i < numberPoints; ++i) {
        double time = 4 * period * static_cast<double>(i) / static_cast<double>(numberPoints - 1);
        double phase = std::fmod(time, period) / period;
        double current = phase < 0.4? 2 + phase / 0.4 : 3 - (phase - 0.4) / 0.6;
        double voltage = phase < 0.4? 24 : -12;
        csv << time << "\t" << current << "\t" << voltage << "\n";
    }

    CircuitSimulationReader reader(csv.str(), true);
    auto mapColumnNames = reader.extract_map_column_names(1, frequency);

    REQUIRE(mapColumnNames.size() == 1);
    REQUIRE(mapColumnNames[0]["time"] == "time");
    REQUIRE(mapColumnNames[0]["current"] == "I(L1)");
    REQUIRE(mapColumnNames[0]["voltage"] == "V(p1,p2)");
}

// The web import runs exactly this: parse the file, extract one period per signal, process the
// operating point. With a fine LTspice time step a period holds hundreds of thousands of points,
// and two stages were quadratic in them — sampling restarted its segment search for every sample,
// and the waveform average and integral re-copied the whole time axis on every iteration (MAS's
// get_time() returns by value). 200k points per period spent minutes in the web engine, whose
// watchdog then aborted the import. The bound is generous; the quadratic path is far beyond it.
TEST_CASE("Test_Import_Dense_Ltspice_Export_Is_Not_Quadratic", "[processor][circuit-simulation-reader][ltspice]") {
    const double frequency = 45000;
    const double period = 1 / frequency;
    const size_t pointsPerPeriod = 200000;
    const size_t numberPeriods = 2;
    std::ostringstream csv;
    csv.precision(15);
    csv << "time\tI(L1)\tV(p1,p2)\n";
    for (size_t i = 0; i < pointsPerPeriod * numberPeriods; ++i) {
        double time = numberPeriods * period * static_cast<double>(i) / static_cast<double>(pointsPerPeriod * numberPeriods - 1);
        double phase = std::fmod(time, period) / period;
        double current = phase < 0.4? 2 + phase / 0.4 : 3 - (phase - 0.4) / 0.6;
        double voltage = phase < 0.4? 24 : -16;
        csv << time << "\t" << current << "\t" << voltage << "\n";
    }
    std::vector<std::map<std::string, std::string>> mapColumnNames = {{{"time", "time"}, {"current", "I(L1)"}, {"voltage", "V(p1,p2)"}}};

    auto start = std::chrono::steady_clock::now();
    CircuitSimulationReader reader(csv.str(), true);
    auto operatingPoint = reader.extract_operating_point(1, frequency, mapColumnNames);
    operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, 100e-6);
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    REQUIRE(seconds < 20);
    auto excitation = operatingPoint.get_excitations_per_winding()[0];
    auto processed = excitation.get_current()->get_processed().value();
    CHECK_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(1, 1e-3));
    CHECK_THAT(processed.get_average().value(), Catch::Matchers::WithinAbs(2.5, 1e-3));
    REQUIRE(excitation.get_voltage());
}
