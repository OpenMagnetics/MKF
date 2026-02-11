#pragma once
#include "json.hpp"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "processors/CircuitSimulatorInterface.h"
#include "support/Utils.h"
#include <string>
#include <vector>
#include <map>
#include <optional>
#include <filesystem>

#ifdef ENABLE_NGSPICE
#include "sharedspice.h"
#endif

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief Result of a circuit simulation containing extracted waveforms
 */
struct SimulationResult {
    bool success = false;
    std::string errorMessage;
    std::vector<Waveform> waveforms;
    std::vector<std::string> waveformNames;
    std::optional<OperatingPoint> operatingPoint;
    double simulationTime = 0.0;  // Wall-clock time in seconds
};

/**
 * @brief Configuration for running a simulation
 */
struct SimulationConfig {
    double stopTime = 0.0;              // Simulation stop time (0 = auto from netlist)
    double stepSize = 0.0;              // Step size (0 = auto)
    size_t steadyStateCycles = 5;       // Number of cycles to skip for steady state
    double frequency = 0.0;             // Operating frequency (for waveform extraction)
    bool extractOnePeriod = true;       // Extract only one period of steady-state
    size_t numberOfPeriods = 2;         // Number of periods to extract (when extractOnePeriod = false)
    std::string workingDirectory = "";  // Working directory for simulation files
    bool keepTempFiles = false;         // Keep temporary files after simulation
    double timeout = 60.0;              // Timeout in seconds
};

/**
 * @brief NgspiceRunner - Runs ngspice simulations and extracts waveforms
 * 
 * This class provides two modes of operation:
 * 1. Command-line mode: Runs ngspice as an external process (always available)
 * 2. Shared library mode: Uses ngspice shared library API (requires ENABLE_NGSPICE)
 * 
 * Usage:
 * @code
 *   NgspiceRunner runner;
 *   
 *   // Run from netlist string
 *   SimulationConfig config;
 *   config.frequency = 100e3;
 *   auto result = runner.run_simulation(netlistString, config);
 *   
 *   // Run from file
 *   auto result = runner.run_simulation_file("/path/to/circuit.cir", config);
 *   
 *   // Run with magnetic model
 *   auto result = runner.simulate_magnetic_circuit(magnetic, converterNetlist, operatingFrequency);
 * @endcode
 */
class NgspiceRunner {
public:
    enum class ExecutionMode {
        COMMAND_LINE,       // Run ngspice as external process
        SHARED_LIBRARY      // Use ngspice shared library (requires ENABLE_NGSPICE)
    };

private:
    ExecutionMode _mode;
    std::string _ngspicePath;  // Path to ngspice executable (for COMMAND_LINE mode)
    bool _verbose = false;
    
#ifdef ENABLE_NGSPICE
    // Shared library mode data
    static NgspiceRunner* _instance;  // For callback routing
    std::vector<std::string> _capturedOutput;
    std::vector<double> _timeData;
    std::map<std::string, std::vector<double>> _vectorData;
    bool _simulationComplete = false;
    bool _simulationError = false;
    std::string _errorMessage;
    
    // Callback functions for ngspice shared library
    static int ng_getchar(char* outputreturn, int ident, void* userdata);
    static int ng_getstat(char* outputreturn, int ident, void* userdata);
    static int ng_exit(int exitstatus, bool immediate, bool quitexit, int ident, void* userdata);
    static int ng_data(vecvaluesall* vecvals, int numvecs, int ident, void* userdata);
    static int ng_initdata(vecinfoall* vecinfo, int ident, void* userdata);
    static int ng_thread_runs(bool running, int ident, void* userdata);
#endif

    // Helper methods
    std::string create_temp_directory();
    std::string write_netlist_to_file(const std::string& netlist, const std::string& directory);
    std::string find_ngspice_executable();
    SimulationResult run_command_line(const std::string& netlistPath, const SimulationConfig& config);
    SimulationResult run_shared_library(const std::string& netlist, const SimulationConfig& config);
    SimulationResult parse_raw_file(const std::string& rawFilePath, const SimulationConfig& config);
    SimulationResult parse_csv_output(const std::string& csvFilePath, const SimulationConfig& config);

public:
    /**
     * @brief Construct NgspiceRunner with automatic mode detection
     * 
     * Will use shared library mode if ENABLE_NGSPICE is defined and ngspice
     * is properly initialized, otherwise falls back to command-line mode.
     */
    NgspiceRunner();
    
    /**
     * @brief Construct NgspiceRunner with specific execution mode
     * @param mode Execution mode to use
     * @param ngspicePath Path to ngspice executable (for COMMAND_LINE mode)
     */
    NgspiceRunner(ExecutionMode mode, const std::string& ngspicePath = "");
    
    ~NgspiceRunner();
    
    /**
     * @brief Set verbose mode for debugging output
     */
    void set_verbose(bool verbose) { _verbose = verbose; }
    
    /**
     * @brief Get current execution mode
     */
    ExecutionMode get_mode() const { return _mode; }
    
    /**
     * @brief Check if ngspice is available
     * @return true if ngspice can be executed
     */
    bool is_available();
    
    /**
     * @brief Run simulation from netlist string
     * @param netlist SPICE netlist as string
     * @param config Simulation configuration
     * @return SimulationResult containing waveforms or error information
     */
    SimulationResult run_simulation(const std::string& netlist, const SimulationConfig& config);
    
    /**
     * @brief Run simulation from netlist file
     * @param netlistPath Path to SPICE netlist file
     * @param config Simulation configuration
     * @return SimulationResult containing waveforms or error information
     */
    SimulationResult run_simulation_file(const std::string& netlistPath, const SimulationConfig& config);
    
    /**
     * @brief Simulate a circuit containing a magnetic component
     * 
     * This method:
     * 1. Exports the magnetic component as a SPICE subcircuit
     * 2. Combines it with the provided converter circuit
     * 3. Runs the simulation
     * 4. Extracts waveforms for each winding
     * 
     * @param magnetic The magnetic component to simulate
     * @param converterCircuit The converter circuit netlist (should reference the magnetic subcircuit)
     * @param frequency Operating frequency for waveform extraction
     * @param config Additional simulation configuration
     * @return SimulationResult with operating point containing per-winding excitations
     */
    SimulationResult simulate_magnetic_circuit(
        const Magnetic& magnetic,
        const std::string& converterCircuit,
        double frequency,
        const SimulationConfig& config = SimulationConfig());
    
    /**
     * @brief Extract operating point from simulation result
     * 
     * Uses CircuitSimulationReader to parse waveforms and create an OperatingPoint
     * 
     * @param result The simulation result containing waveforms
     * @param numberWindings Number of windings to extract
     * @param frequency Operating frequency
     * @param ambientTemperature Ambient temperature for the operating point
     * @return OperatingPoint with per-winding excitations
     */
    static OperatingPoint extract_operating_point(
        const SimulationResult& result,
        size_t numberWindings,
        double frequency,
        double ambientTemperature = 25.0);
    
    /**
     * @brief Mapping for waveform names per winding
     * 
     * Each entry in the vector corresponds to a winding (index 0 = primary, etc.)
     * The map keys can be "voltage" and "current", values are the waveform names
     * to look for in the SimulationResult.
     */
    using WaveformNameMapping = std::vector<std::map<std::string, std::string>>;
    
    /**
     * @brief Extract operating point from simulation result with explicit waveform mapping
     * 
     * This overload allows specifying exactly which waveform names correspond to
     * voltage and current for each winding.
     * 
     * @param result The simulation result containing waveforms
     * @param waveformMapping Per-winding mapping of signal type to waveform name
     *        Example: { {{"voltage", "v(pri_in)"}, {"current", "i(vpri_sense)"}},
     *                   {{"voltage", "v(vout)"}, {"current", "i(vsec_sense)"}} }
     * @param frequency Operating frequency
     * @param windingNames Optional names for each winding (e.g., "Primary", "Secondary 0")
     * @param ambientTemperature Ambient temperature for the operating point
     * @param flipCurrentSign Per-winding flag to flip current sign (ngspice convention)
     * @return OperatingPoint with per-winding excitations
     */
    static OperatingPoint extract_operating_point(
        const SimulationResult& result,
        const WaveformNameMapping& waveformMapping,
        double frequency,
        const std::vector<std::string>& windingNames = {},
        double ambientTemperature = 25.0,
        const std::vector<bool>& flipCurrentSign = {});
        
    /**
     * @brief Generate a simple test circuit for a magnetic component
     * 
     * Creates a basic test circuit that applies voltage to the primary
     * and loads the secondaries with resistors.
     * 
     * @param magnetic The magnetic component
     * @param frequency Operating frequency
     * @param primaryVoltage Primary winding voltage (peak)
     * @param loadResistances Load resistance for each secondary (in ohms)
     * @return Complete netlist string ready for simulation
     */
    static std::string generate_test_circuit(
        const Magnetic& magnetic,
        double frequency,
        double primaryVoltage,
        const std::vector<double>& loadResistances);
};

} // namespace OpenMagnetics
