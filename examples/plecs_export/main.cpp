// PLECS export example for MKF.
//
// Usage:
//   plecs_export_example <input.json> <output.plecs> <subcircuit|symbol> [frequency_Hz] [temperature_C]
//
// Reads a MAS Magnetic JSON file and writes a .plecs snippet that can be
// copy-pasted into a larger PLECS schematic file.

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#include <nlohmann/json.hpp>

#include "processors/CircuitSimulatorInterface.h"

namespace {

void print_usage(std::ostream& os) {
    os << "Usage: plecs_export_example <input.json> <output.plecs> "
       << "<subcircuit|symbol> [frequency_Hz] [temperature_C]\n"
       << "  subcircuit mode requires frequency_Hz and temperature_C.\n"
       << "  symbol mode ignores frequency_Hz and temperature_C.\n";
}

int fail(const std::string& message) {
    std::cerr << "error: " << message << "\n";
    print_usage(std::cerr);
    return 1;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        return fail("not enough arguments");
    }

    const std::string input_path = argv[1];
    const std::string output_path = argv[2];
    const std::string mode = argv[3];

    if (mode != "subcircuit" && mode != "symbol") {
        return fail("mode must be 'subcircuit' or 'symbol' (got '" + mode + "')");
    }

    double frequency = 0.0;
    double temperature = 0.0;
    if (mode == "subcircuit") {
        if (argc < 6) {
            return fail("subcircuit mode requires frequency_Hz and temperature_C");
        }
        try {
            frequency = std::stod(argv[4]);
            temperature = std::stod(argv[5]);
        } catch (const std::exception& e) {
            return fail(std::string{"could not parse frequency/temperature: "} + e.what());
        }
    }

    std::ifstream input_file(input_path);
    if (!input_file.is_open()) {
        return fail("could not open input file: " + input_path);
    }

    nlohmann::json j;
    try {
        input_file >> j;
    } catch (const std::exception& e) {
        return fail(std::string{"failed to parse JSON: "} + e.what());
    }

    std::string content;
    try {
        OpenMagnetics::Magnetic magnetic(j);
        OpenMagnetics::CircuitSimulatorExporterPlecsModel exporter;
        if (mode == "subcircuit") {
            content = exporter.export_magnetic_as_subcircuit(magnetic, frequency, temperature);
        } else {
            content = exporter.export_magnetic_as_symbol(magnetic);
        }
    } catch (const std::exception& e) {
        return fail(std::string{"PLECS export failed: "} + e.what());
    }

    std::ofstream output_file(output_path);
    if (!output_file.is_open()) {
        return fail("could not open output file for writing: " + output_path);
    }
    output_file << content;
    if (!output_file) {
        return fail("failed while writing output file: " + output_path);
    }

    std::cout << "Wrote " << output_path << "\n";
    return 0;
}
