#include "processors/NgspiceRunner.h"
#include "processors/CircuitSimulatorInterface.h"
#include "support/Exceptions.h"
#include "Defaults.h"

#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>
#include <regex>
#include <cstdlib>
#include <array>
#include <memory>
#include <iostream>
#include <algorithm>
#include <cctype>

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#define popen _popen
#define pclose _pclose
#else
#include <unistd.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <signal.h>
#endif

namespace OpenMagnetics {

#ifdef ENABLE_NGSPICE
// Static instance for callback routing
NgspiceRunner* NgspiceRunner::_instance = nullptr;

// Callback implementations for ngspice shared library
int NgspiceRunner::ng_getchar(char* outputreturn, int ident, void* userdata) {
    if (_instance) {
        _instance->_capturedOutput.push_back(std::string(outputreturn));
    }
    return 0;
}

int NgspiceRunner::ng_getstat(char* outputreturn, int ident, void* userdata) {
    // Detect simulation completion via status message
    // In WASM, ng_thread_runs may not work properly, so we detect "--ready--"
    if (_instance && outputreturn) {
        std::string status(outputreturn);
        if (status.find("--ready--") != std::string::npos) {
            _instance->_simulationComplete = true;
        }
    }
    return 0;
}

int NgspiceRunner::ng_exit(int exitstatus, bool immediate, bool quitexit, int ident, void* userdata) {
    if (_instance) {
        _instance->_simulationComplete = true;
        if (exitstatus != 0) {
            _instance->_simulationError = true;
            _instance->_errorMessage = "ngspice exited with status " + std::to_string(exitstatus);
        }
    }
    return exitstatus;
}

int NgspiceRunner::ng_data(vecvaluesall* vecvals, int numvecs, int ident, void* userdata) {
    if (!_instance) return 0;
    
    // Store vector data as it comes in
    for (int i = 0; i < numvecs; i++) {
        std::string name = vecvals->vecsa[i]->name;
        double value = vecvals->vecsa[i]->creal;  // Real part
        
        // Check for time vector - ngspice may prefix with plot name like "tran1.time"
        // Use case-insensitive check for "time" at the end of the name
        std::string lowerName = name;
        std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(), ::tolower);
        bool isTime = (lowerName == "time") || 
                      (lowerName.size() > 5 && lowerName.substr(lowerName.size() - 5) == ".time");
        
        if (isTime) {
            _instance->_timeData.push_back(value);
        } else {
            _instance->_vectorData[name].push_back(value);
        }
    }
    return 0;
}

int NgspiceRunner::ng_initdata(vecinfoall* vecinfo, int ident, void* userdata) {
    if (_instance && _instance->_verbose) {
        std::cout << "ngspice init: " << vecinfo->veccount << " vectors" << std::endl;
    }
    return 0;
}

int NgspiceRunner::ng_thread_runs(bool running, int ident, void* userdata) {
    if (_instance) {
        if (!running) {
            _instance->_simulationComplete = true;
        }
    }
    return 0;
}
#endif

NgspiceRunner::NgspiceRunner() {
#ifdef ENABLE_NGSPICE
    _mode = ExecutionMode::SHARED_LIBRARY;
    _instance = this;
    
    // Initialize ngspice shared library
    int ret = ngSpice_Init(ng_getchar, ng_getstat, ng_exit, ng_data, ng_initdata, ng_thread_runs, nullptr);
    if (ret != 0) {
        std::cout << "Warning: Failed to initialize ngspice shared library, falling back to command-line mode" << std::endl;
        _mode = ExecutionMode::COMMAND_LINE;
    }
#else
    _mode = ExecutionMode::COMMAND_LINE;
#endif
    _ngspicePath = find_ngspice_executable();
}

NgspiceRunner::NgspiceRunner(ExecutionMode mode, const std::string& ngspicePath) 
    : _mode(mode), _ngspicePath(ngspicePath) {
    
#ifdef ENABLE_NGSPICE
    if (_mode == ExecutionMode::SHARED_LIBRARY) {
        _instance = this;
        int ret = ngSpice_Init(ng_getchar, ng_getstat, ng_exit, ng_data, ng_initdata, ng_thread_runs, nullptr);
        if (ret != 0) {
            throw std::runtime_error("Failed to initialize ngspice shared library");
        }
    }
#else
    if (_mode == ExecutionMode::SHARED_LIBRARY) {
        throw std::runtime_error("ENABLE_NGSPICE not defined - shared library mode unavailable");
    }
#endif
    
    if (_mode == ExecutionMode::COMMAND_LINE && _ngspicePath.empty()) {
        _ngspicePath = find_ngspice_executable();
    }
}

NgspiceRunner::~NgspiceRunner() {
#ifdef ENABLE_NGSPICE
    if (_mode == ExecutionMode::SHARED_LIBRARY && _instance == this) {
        _instance = nullptr;
    }
#endif
}

std::string NgspiceRunner::find_ngspice_executable() {
    // Common locations for ngspice
    std::vector<std::string> candidates = {
        "ngspice",                          // In PATH
        "/usr/bin/ngspice",
        "/usr/local/bin/ngspice",
        "/opt/ngspice/bin/ngspice",
#ifdef _WIN32
        "C:\\Program Files\\ngspice\\bin\\ngspice.exe",
        "C:\\ngspice\\bin\\ngspice.exe",
#endif
    };
    
    for (const auto& path : candidates) {
        // Try to execute with --version to check if it works
#ifdef _WIN32
        std::string cmd = "\"" + path + "\" --version >nul 2>&1";
#else
        std::string cmd = path + " --version >/dev/null 2>&1";
#endif
        if (std::system(cmd.c_str()) == 0) {
            return path;
        }
    }
    
    return "ngspice";  // Default, hope it's in PATH
}

bool NgspiceRunner::is_available() {
#ifdef ENABLE_NGSPICE
    if (_mode == ExecutionMode::SHARED_LIBRARY) {
        return true;  // Already initialized in constructor
    }
#endif
    
    // Check command-line availability
#ifdef _WIN32
    std::string cmd = "\"" + _ngspicePath + "\" --version >nul 2>&1";
#else
    std::string cmd = _ngspicePath + " --version >/dev/null 2>&1";
#endif
    return std::system(cmd.c_str()) == 0;
}

std::string NgspiceRunner::create_temp_directory() {
    std::string tempDir;
    
#ifdef _WIN32
    char tempPath[MAX_PATH];
    GetTempPathA(MAX_PATH, tempPath);
    tempDir = std::string(tempPath) + "ngspice_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    _mkdir(tempDir.c_str());
#else
    tempDir = "/tmp/ngspice_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    mkdir(tempDir.c_str(), 0755);
#endif
    
    return tempDir;
}

std::string NgspiceRunner::write_netlist_to_file(const std::string& netlist, const std::string& directory) {
    std::string filePath = directory + "/circuit.cir";
    std::ofstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to create netlist file: " + filePath);
    }
    file << netlist;
    file.close();
    return filePath;
}

SimulationResult NgspiceRunner::run_simulation(const std::string& netlist, const SimulationConfig& config) {
    std::string workDir = config.workingDirectory;
    if (workDir.empty()) {
        workDir = create_temp_directory();
    }
    
    std::string netlistPath = write_netlist_to_file(netlist, workDir);
    
    SimulationResult result;
    
#ifdef ENABLE_NGSPICE
    if (_mode == ExecutionMode::SHARED_LIBRARY) {
        result = run_shared_library(netlist, config);
    } else {
        result = run_command_line(netlistPath, config);
    }
#else
    result = run_command_line(netlistPath, config);
#endif
    
    // Cleanup temp files if requested
    if (!config.keepTempFiles && config.workingDirectory.empty()) {
        std::filesystem::remove_all(workDir);
    }
    
    return result;
}

SimulationResult NgspiceRunner::run_simulation_file(const std::string& netlistPath, const SimulationConfig& config) {
#ifdef ENABLE_NGSPICE
    if (_mode == ExecutionMode::SHARED_LIBRARY) {
        // Read file and run via shared library
        std::ifstream file(netlistPath);
        if (!file.is_open()) {
            SimulationResult result;
            result.success = false;
            result.errorMessage = "Failed to open netlist file: " + netlistPath;
            return result;
        }
        std::stringstream buffer;
        buffer << file.rdbuf();
        return run_shared_library(buffer.str(), config);
    }
#endif
    return run_command_line(netlistPath, config);
}

SimulationResult NgspiceRunner::run_command_line(const std::string& netlistPath, const SimulationConfig& config) {
    SimulationResult result;
    auto startTime = std::chrono::steady_clock::now();
    
    // Get directory of netlist
    std::filesystem::path netlistDir = std::filesystem::path(netlistPath).parent_path();
    std::string outputFile = (netlistDir / "output.csv").string();
    std::string rawFile = (netlistDir / "output.raw").string();
    
    // Build ngspice command
    // We'll add commands to save data and write to CSV
    std::string controlScript = netlistDir.string() + "/control.sp";
    {
        std::ofstream ctrlFile(controlScript);
        ctrlFile << "* Control script for batch simulation\n";
        ctrlFile << ".include " << netlistPath << "\n\n";
        ctrlFile << ".control\n";
        ctrlFile << "set filetype=ascii\n";  // ASCII raw file for easier parsing
        ctrlFile << "run\n";
        ctrlFile << "write " << rawFile << " all\n";
        ctrlFile << "wrdata " << outputFile << " all\n";
        ctrlFile << "quit\n";
        ctrlFile << ".endc\n";
        ctrlFile << ".end\n";
        ctrlFile.close();
    }
    
    // Run ngspice in batch mode
#ifdef _WIN32
    std::string cmd = "\"" + _ngspicePath + "\" -b \"" + controlScript + "\" 2>&1";
#else
    std::string cmd = _ngspicePath + " -b \"" + controlScript + "\" 2>&1";
#endif
    
    if (_verbose) {
        std::cerr << "Running: " << cmd << std::endl;
    }
    
    // Execute and capture output
    std::array<char, 4096> buffer;
    std::string cmdOutput;
    
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        result.success = false;
        result.errorMessage = "Failed to execute ngspice";
        return result;
    }
    
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        cmdOutput += buffer.data();
    }
    
    int exitCode = pclose(pipe);
    
    auto endTime = std::chrono::steady_clock::now();
    result.simulationTime = std::chrono::duration<double>(endTime - startTime).count();
    
    if (_verbose) {
        std::cerr << "ngspice output:\n" << cmdOutput << std::endl;
    }
    
    // Check for errors in output
    if (cmdOutput.find("Error") != std::string::npos || 
        cmdOutput.find("error") != std::string::npos ||
        cmdOutput.find("ERROR") != std::string::npos) {
        // Check if it's a fatal error
        if (cmdOutput.find("Fatal") != std::string::npos || exitCode != 0) {
            result.success = false;
            result.errorMessage = "ngspice error: " + cmdOutput;
            return result;
        }
    }
    
    // Try to parse raw file first (more complete data)
    if (std::filesystem::exists(rawFile)) {
        result = parse_raw_file(rawFile, config);
        result.simulationTime = std::chrono::duration<double>(endTime - startTime).count();
    } else if (std::filesystem::exists(outputFile)) {
        result = parse_csv_output(outputFile, config);
        result.simulationTime = std::chrono::duration<double>(endTime - startTime).count();
    } else {
        result.success = false;
        result.errorMessage = "No output file generated by ngspice";
    }
    
    return result;
}

#ifdef ENABLE_NGSPICE
SimulationResult NgspiceRunner::run_shared_library(const std::string& netlist, const SimulationConfig& config) {
    SimulationResult result;
    auto startTime = std::chrono::steady_clock::now();
    
    // Clear previous data
    _capturedOutput.clear();
    _timeData.clear();
    _vectorData.clear();
    _simulationComplete = false;
    _simulationError = false;
    _errorMessage.clear();
    
    // Send netlist to ngspice
    // IMPORTANT: Build lineStorage FIRST completely, THEN build pointers
    // Otherwise vector reallocation invalidates c_str() pointers
    std::vector<std::string> lineStorage;
    std::istringstream stream(netlist);
    std::string line;
    
    int lineNum = 0;
    while (std::getline(stream, line)) {
        lineStorage.push_back(line);
    }
    
    // Now build pointer array - lineStorage won't reallocate anymore
    std::vector<char*> lines;
    lines.reserve(lineStorage.size() + 1);
    for (auto& s : lineStorage) {
        lines.push_back(const_cast<char*>(s.c_str()));
    }
    lines.push_back(nullptr);
    
    int ret = ngSpice_Circ(lines.data());
    
    // Print any captured output from ngspice (errors should appear here)
    if (!_capturedOutput.empty()) {
    
    if (ret != 0) {
        result.success = false;
        result.errorMessage = "Failed to load circuit into ngspice (ret=" + std::to_string(ret) + ")";
        if (!_capturedOutput.empty()) {
            result.errorMessage += ". Last ngspice output: " + _capturedOutput.back();
        }
        return result;
    }
    
    // Run simulation
    ret = ngSpice_Command(const_cast<char*>("run"));
    if (ret != 0) {
        result.success = false;
        result.errorMessage = "Failed to run simulation";
        return result;
    }
    
    // Wait for completion with timeout
    auto timeoutEnd = startTime + std::chrono::duration<double>(config.timeout);
    int timeoutCheckCount = 0;
    while (!_simulationComplete) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        timeoutCheckCount++;
        if (std::chrono::steady_clock::now() > timeoutEnd) {
            auto elapsedTime = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
            ngSpice_Command(const_cast<char*>("stop"));
            result.success = false;
            result.errorMessage = "Simulation timeout after " + std::to_string(elapsedTime) + "s (limit: " + std::to_string(config.timeout) + "s)";
            return result;
        }
    }
    
    auto endTime = std::chrono::steady_clock::now();
    result.simulationTime = std::chrono::duration<double>(endTime - startTime).count();
    
    if (_simulationError) {
        result.success = false;
        result.errorMessage = _errorMessage;
        return result;
    }
    
    // Convert stored data to waveforms
    if (_timeData.empty()) {
        result.success = false;
        result.errorMessage = "No time data captured from simulation";
        return result;
    }
    
    result.success = true;
    
    // Create time waveform
    Waveform timeWaveform;
    timeWaveform.set_time(_timeData);
    timeWaveform.set_data(_timeData);
    result.waveforms.push_back(timeWaveform);
    result.waveformNames.push_back("time");
    
    // Create waveforms for each vector
    for (const auto& [name, data] : _vectorData) {
        if (data.size() == _timeData.size()) {
            Waveform waveform;
            waveform.set_time(_timeData);
            waveform.set_data(data);
            result.waveforms.push_back(waveform);
            result.waveformNames.push_back(name);
        }
    }
    
    // Extract periods if requested
    // For converter waveforms, find the switch-ON edge and extract N periods for better visualization
    // We search backwards to use the last (settled) periods after initial transients have decayed
    if (config.extractOnePeriod && config.frequency > 0) {
        double period = 1.0 / config.frequency;
        const size_t numPeriodsToExtract = config.numberOfPeriods;  // Use configured number of periods
        
        // Get time data from first waveform (index 1, since 0 is time itself)
        if (result.waveforms.size() > 1 && result.waveforms[1].get_time()) {
            auto time = result.waveforms[1].get_time().value();
            auto voltageData = result.waveforms[1].get_data();  // First signal is typically voltage
            
            // Find voltage range using percentiles to ignore spikes
            // Sort a copy to find percentiles
            std::vector<double> sortedVoltage = voltageData;
            std::sort(sortedVoltage.begin(), sortedVoltage.end());
            size_t p5_idx = sortedVoltage.size() * 5 / 100;
            size_t p95_idx = sortedVoltage.size() * 95 / 100;
            double vMin = sortedVoltage[p5_idx];
            double vMax = sortedVoltage[p95_idx];
            double vRange = vMax - vMin;
            double threshold = vMin + vRange * 0.5;  // 50% of range as threshold
            
            // Search BACKWARDS to find the LAST rising edge that allows full period extraction
            // This uses settled waveforms after initial transients have decayed
            double minEdgeTime = time.back() - numPeriodsToExtract * period;
            size_t edgeIndex = 0;
            
            // Find the last rising edge before minEdgeTime
            for (size_t i = time.size() - 1; i > 0; --i) {
                if (time[i] <= minEdgeTime) {
                    // Look for rising edge crossing the threshold
                    if (voltageData[i] > threshold && voltageData[i-1] <= threshold) {
                        edgeIndex = i;
                        break;
                    }
                }
            }
            
            // If no edge found searching backwards, try forwards as fallback
            if (edgeIndex == 0) {
                for (size_t i = 1; i < time.size(); ++i) {
                    if (time[i] <= minEdgeTime) {
                        if (voltageData[i] > threshold && voltageData[i-1] <= threshold) {
                            edgeIndex = i;
                            break;
                        }
                    }
                }
            }
            
            // If still no edge found, start from where we have numPeriodsToExtract periods left
            if (edgeIndex == 0) {
                for (size_t i = 0; i < time.size(); ++i) {
                    if (time[i] >= time.back() - numPeriodsToExtract * period) {
                        edgeIndex = i;
                        break;
                    }
                }
            }
            
            // Find period end (numPeriodsToExtract periods after edge, or end of data)
            size_t periodEndIndex = time.size();  // Default to end of data
            double extractStartTime = time[edgeIndex];
            double targetEndTime = extractStartTime + numPeriodsToExtract * period;
            double tolerance = period * 0.001;  // 0.1% tolerance
            for (size_t i = edgeIndex; i < time.size(); ++i) {
                if (time[i] >= targetEndTime - tolerance) {
                    periodEndIndex = i + 1;
                    break;
                }
            }
            
            // Ensure we have at least some data
            if (periodEndIndex <= edgeIndex) {
                periodEndIndex = time.size();
            }
            
            // Extract periods for all waveforms using the found indices
            for (size_t i = 1; i < result.waveforms.size(); ++i) {
                auto wfTime = result.waveforms[i].get_time().value();
                auto wfData = result.waveforms[i].get_data();
                
                // Clamp indices to valid range
                size_t startIdx = std::min(edgeIndex, wfData.size() - 1);
                size_t endIdx = std::min(periodEndIndex, wfData.size());
                
                if (endIdx > startIdx) {
                    auto periodTime = std::vector<double>(wfTime.begin() + startIdx, wfTime.begin() + endIdx);
                    auto periodData = std::vector<double>(wfData.begin() + startIdx, wfData.begin() + endIdx);
                    
                    // Check if this is a voltage waveform (not a current measurement)
                    // Current measurements in ngspice have "#branch" in their name (e.g., "vpri_sense#branch")
                    std::string wfName = result.waveformNames[i];
                    std::transform(wfName.begin(), wfName.end(), wfName.begin(), ::tolower);
                    bool isCurrent = (wfName.find("#branch") != std::string::npos);
                    bool isVoltage = !isCurrent && (wfName != "time");
                    
                    // Note: Voltage clipping disabled - flyback waveforms have legitimate wide voltage swings
                    // The secondary winding voltage swings negative during ON and positive during OFF
                    // Clipping would distort the actual waveform shape
                    
                    // Offset time to start at 0
                    double offset = periodTime[0];
                    for (auto& t : periodTime) {
                        t -= offset;
                    }
                    
                    Waveform newWaveform;
                    newWaveform.set_time(periodTime);
                    newWaveform.set_data(periodData);
                    
                    // Store the waveform directly without resampling (preserve all data points)
                    result.waveforms[i] = newWaveform;
                }
            }
            
            // Update the "time" waveform at index 0 to match the processed waveforms
            // This is important because some code extracts "time" separately for the x-axis
            if (result.waveforms.size() > 1 && result.waveforms[1].get_time()) {
                auto processedTime = result.waveforms[1].get_time().value();
                Waveform updatedTimeWaveform;
                updatedTimeWaveform.set_time(processedTime);
                updatedTimeWaveform.set_data(processedTime);
                result.waveforms[0] = updatedTimeWaveform;
            }
        }
    }
    
    return result;
}
#endif

SimulationResult NgspiceRunner::parse_raw_file(const std::string& rawFilePath, const SimulationConfig& config) {
    SimulationResult result;
    std::ifstream file(rawFilePath);
    
    if (!file.is_open()) {
        result.success = false;
        result.errorMessage = "Failed to open raw file: " + rawFilePath;
        return result;
    }
    
    std::string line;
    std::vector<std::string> variableNames;
    size_t numPoints = 0;
    size_t numVariables = 0;
    bool inValues = false;
    bool isReal = false;
    bool isComplex = false;
    
    // Parse header
    while (std::getline(file, line)) {
        if (line.find("No. Variables:") != std::string::npos) {
            numVariables = std::stoul(line.substr(line.find(":") + 1));
        } else if (line.find("No. Points:") != std::string::npos) {
            numPoints = std::stoul(line.substr(line.find(":") + 1));
        } else if (line.find("Flags: real") != std::string::npos) {
            isReal = true;
        } else if (line.find("Flags: complex") != std::string::npos) {
            isComplex = true;
        } else if (line.find("Variables:") != std::string::npos) {
            // Read variable names
            for (size_t i = 0; i < numVariables; ++i) {
                std::getline(file, line);
                std::istringstream iss(line);
                size_t idx;
                std::string name, type;
                iss >> idx >> name >> type;
                variableNames.push_back(name);
            }
        } else if (line.find("Values:") != std::string::npos || line.find("Binary:") != std::string::npos) {
            inValues = true;
            break;
        }
    }
    
    if (!inValues || variableNames.empty()) {
        result.success = false;
        result.errorMessage = "Invalid raw file format";
        return result;
    }
    
    // Initialize data vectors
    std::vector<std::vector<double>> data(numVariables);
    
    // Parse values (ASCII format)
    // ngspice raw format: 
    //   <index>\t<first_value>
    //   \t<second_value>
    //   \t<third_value>
    //   ... (one value per line, indented with tabs or spaces)
    // Note: index lines may start with spaces before the index number
    size_t currentVar = 0;
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;
        
        // Trim leading spaces and check if this is an index line
        std::string trimmedLine = line;
        size_t firstNonSpace = line.find_first_not_of(" ");
        if (firstNonSpace != std::string::npos) {
            trimmedLine = line.substr(firstNonSpace);
        }
        
        // Check if this is an index line (trimmed line starts with digit and has tab separator)
        if (!trimmedLine.empty() && std::isdigit(static_cast<unsigned char>(trimmedLine[0])) 
            && trimmedLine.find('\t') != std::string::npos) {
            // This is a new data point - reset variable counter
            currentVar = 0;
            
            // Parse the first value after the tab
            size_t tabPos = trimmedLine.find('\t');
            if (tabPos != std::string::npos) {
                std::string valueStr = trimmedLine.substr(tabPos + 1);
                // Handle complex numbers (comma-separated real,imag)
                size_t commaPos = valueStr.find(',');
                if (commaPos != std::string::npos) {
                    valueStr = valueStr.substr(0, commaPos);
                }
                try {
                    double value = std::stod(valueStr);
                    if (currentVar < numVariables) {
                        data[currentVar].push_back(value);
                        currentVar++;
                    }
                } catch (...) {
                    // Skip invalid values
                }
            }
        } else if (!trimmedLine.empty() && (line[0] == '\t' || line[0] == ' ')) {
            // This is a continuation line with next variable value
            std::string valueStr = trimmedLine;
            // Handle complex numbers
            size_t commaPos = valueStr.find(',');
            if (commaPos != std::string::npos) {
                valueStr = valueStr.substr(0, commaPos);
            }
            try {
                double value = std::stod(valueStr);
                if (currentVar < numVariables) {
                    data[currentVar].push_back(value);
                    currentVar++;
                }
            } catch (...) {
                // Skip invalid values
            }
        }
    }
    
    file.close();
    
    if (data.empty() || data[0].empty()) {
        result.success = false;
        result.errorMessage = "No data parsed from raw file";
        return result;
    }
    
    result.success = true;
    
    // Convert to waveforms
    std::vector<double> timeData = data[0];  // First variable is typically time
    
    for (size_t i = 0; i < variableNames.size(); ++i) {
        Waveform waveform;
        waveform.set_time(timeData);
        waveform.set_data(data[i]);
        result.waveforms.push_back(waveform);
        
        // Normalize signal names from ngspice raw file format:
        // - v(node) -> node (voltage signals)
        // - i(source) -> source#branch (current signals)
        std::string signalName = variableNames[i];
        if (signalName.size() > 3 && signalName[0] == 'v' && signalName[1] == '(' && signalName.back() == ')') {
            // v(node) -> node
            signalName = signalName.substr(2, signalName.size() - 3);
        } else if (signalName.size() > 3 && signalName[0] == 'i' && signalName[1] == '(' && signalName.back() == ')') {
            // i(source) -> source#branch
            signalName = signalName.substr(2, signalName.size() - 3) + "#branch";
        }
        result.waveformNames.push_back(signalName);
    }
    
    // Extract one period if requested
    if (config.extractOnePeriod && config.frequency > 0 && result.waveforms.size() > 1) {
        CircuitSimulationReader reader;
        for (size_t i = 1; i < result.waveforms.size(); ++i) {
            result.waveforms[i] = reader.get_one_period(result.waveforms[i], config.frequency, true);
        }
        // Update the time waveform at index 0 to match the extracted period
        if (result.waveforms[1].get_time()) {
            auto processedTime = result.waveforms[1].get_time().value();
            Waveform updatedTimeWaveform;
            updatedTimeWaveform.set_time(processedTime);
            updatedTimeWaveform.set_data(processedTime);
            result.waveforms[0] = updatedTimeWaveform;
        }
    }
    
    return result;
}

SimulationResult NgspiceRunner::parse_csv_output(const std::string& csvFilePath, const SimulationConfig& config) {
    SimulationResult result;
    
    try {
        CircuitSimulationReader reader(csvFilePath);
        
        // Use CircuitSimulationReader's functionality to parse
        auto columnNames = reader.extract_column_names();
        
        result.success = true;
        result.waveformNames = columnNames;
        
        // The reader stores waveforms internally; we need to extract them
        // For now, we'll use the existing infrastructure
        
    } catch (const std::exception& e) {
        result.success = false;
        result.errorMessage = std::string("Failed to parse CSV: ") + e.what();
    }
    
    return result;
}

SimulationResult NgspiceRunner::simulate_magnetic_circuit(
    const Magnetic& magnetic,
    const std::string& converterCircuit,
    double frequency,
    const SimulationConfig& config) {
    
    // Generate the magnetic subcircuit
    CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::NGSPICE);
    std::string magneticSubcircuit = exporter.export_magnetic_as_subcircuit(
        const_cast<Magnetic&>(magnetic), 
        frequency, 
        defaults.ambientTemperature);
    
    // Combine subcircuit with converter circuit
    std::string fullNetlist = magneticSubcircuit + "\n\n" + converterCircuit;
    
    // Update config with frequency
    SimulationConfig runConfig = config;
    runConfig.frequency = frequency;
    
    // Run simulation
    SimulationResult result = run_simulation(fullNetlist, runConfig);
    
    if (result.success) {
        // Extract operating point for the magnetic component
        size_t numWindings = magnetic.get_coil().get_functional_description().size();
        result.operatingPoint = extract_operating_point(result, numWindings, frequency);
    }
    
    return result;
}

OperatingPoint NgspiceRunner::extract_operating_point(
    const SimulationResult& result,
    size_t numberWindings,
    double frequency,
    double ambientTemperature) {
    
    // Build a temporary CSV-like string from the waveforms
    std::stringstream csvData;
    
    // Header
    for (size_t i = 0; i < result.waveformNames.size(); ++i) {
        csvData << result.waveformNames[i];
        if (i < result.waveformNames.size() - 1) csvData << ",";
    }
    csvData << "\n";
    
    // Data rows
    if (!result.waveforms.empty() && result.waveforms[0].get_data().size() > 0) {
        size_t numPoints = result.waveforms[0].get_data().size();
        for (size_t row = 0; row < numPoints; ++row) {
            for (size_t col = 0; col < result.waveforms.size(); ++col) {
                auto& data = result.waveforms[col].get_data();
                if (row < data.size()) {
                    csvData << data[row];
                }
                if (col < result.waveforms.size() - 1) csvData << ",";
            }
            csvData << "\n";
        }
    }
    
    // Use CircuitSimulationReader to extract operating point
    CircuitSimulationReader reader(csvData.str());
    return reader.extract_operating_point(numberWindings, frequency, std::nullopt, ambientTemperature);
}

OperatingPoint NgspiceRunner::extract_operating_point(
    const SimulationResult& result,
    const WaveformNameMapping& waveformMapping,
    double frequency,
    const std::vector<std::string>& windingNames,
    double ambientTemperature,
    const std::vector<bool>& flipCurrentSign) {
    
    OperatingPoint operatingPoint;
    
    // Build a map from waveform names to indices for quick lookup
    // Handle ngspice naming: v(node) for voltage, i(vsource) for current
    std::map<std::string, size_t> nameToIndex;
    for (size_t i = 0; i < result.waveformNames.size(); ++i) {
        std::string name = result.waveformNames[i];
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        
        // Store original and lowercase
        nameToIndex[name] = i;
        nameToIndex[lower] = i;
        
        // Strip v() or i() prefix if present and store both versions
        // ngspice returns: v(node), i(vsource) - we need to match both with and without prefix
        if (lower.size() > 2 && (lower[0] == 'v' || lower[0] == 'i') && lower[1] == '(') {
            size_t closeParenPos = lower.rfind(')');
            if (closeParenPos != std::string::npos && closeParenPos > 2) {
                std::string stripped = lower.substr(2, closeParenPos - 2);
                nameToIndex[stripped] = i;
                // For current, ngspice uses i(vsource) but code may look for vsource#branch
                if (lower[0] == 'i') {
                    nameToIndex[stripped + "#branch"] = i;
                }
            }
        }
    }
    
    // Process each winding
    for (size_t windingIndex = 0; windingIndex < waveformMapping.size(); ++windingIndex) {
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        
        // Set winding name if provided
        if (windingIndex < windingNames.size()) {
            excitation.set_name(windingNames[windingIndex]);
        } else {
            excitation.set_name(windingIndex == 0 ? "Primary" : "Secondary " + std::to_string(windingIndex - 1));
        }
        
        const auto& mapping = waveformMapping[windingIndex];
        
        // Look for voltage waveform
        auto voltageIt = mapping.find("voltage");
        if (voltageIt != mapping.end()) {
            std::string voltageName = voltageIt->second;
            std::string voltageNameLower = voltageName;
            std::transform(voltageNameLower.begin(), voltageNameLower.end(), voltageNameLower.begin(), ::tolower);
            
            auto indexIt = nameToIndex.find(voltageNameLower);
            if (indexIt != nameToIndex.end()) {
                SignalDescriptor voltage;
                voltage.set_waveform(result.waveforms[indexIt->second]);
                excitation.set_voltage(voltage);
            }
        }
        
        // Look for current waveform
        auto currentIt = mapping.find("current");
        if (currentIt != mapping.end()) {
            std::string currentName = currentIt->second;
            std::string currentNameLower = currentName;
            std::transform(currentNameLower.begin(), currentNameLower.end(), currentNameLower.begin(), ::tolower);
            
            auto indexIt = nameToIndex.find(currentNameLower);
            if (indexIt != nameToIndex.end()) {
                Waveform currentWaveform = result.waveforms[indexIt->second];
                
                // Flip current sign if requested (common for ngspice sense resistors)
                bool shouldFlip = windingIndex < flipCurrentSign.size() ? flipCurrentSign[windingIndex] : false;
                if (shouldFlip) {
                    auto data = currentWaveform.get_data();
                    for (auto& d : data) d = -d;
                    currentWaveform.set_data(data);
                }
                
                SignalDescriptor current;
                current.set_waveform(currentWaveform);
                excitation.set_current(current);
            }
        }
        
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }
    
    // Set operating conditions
    OperatingConditions conditions;
    conditions.set_ambient_temperature(ambientTemperature);
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);
    
    return operatingPoint;
}

std::string NgspiceRunner::generate_test_circuit(
    const Magnetic& magnetic,
    double frequency,
    double primaryVoltage,
    const std::vector<double>& loadResistances) {
    
    std::stringstream circuit;
    auto& coil = magnetic.get_coil();
    size_t numWindings = coil.get_functional_description().size();
    std::string magneticRef = const_cast<Magnetic&>(magnetic).get_reference();
    
    // Replace spaces and special chars in reference
    std::replace(magneticRef.begin(), magneticRef.end(), ' ', '_');
    
    circuit << "* Test circuit for " << magneticRef << "\n\n";
    
    // Voltage source on primary
    double period = 1.0 / frequency;
    circuit << "* Primary voltage source (rectangular wave)\n";
    circuit << "Vpri P1+ 0 PULSE(0 " << primaryVoltage << " 0 1n 1n " << (period/2) << " " << period << ")\n";
    
    // Ground connection for primary negative
    circuit << "V_gnd_pri P1- 0 0\n";
    
    // Load resistors on secondaries
    for (size_t i = 1; i < numWindings; ++i) {
        double loadR = (i-1 < loadResistances.size()) ? loadResistances[i-1] : 100.0;
        circuit << "Rload" << i << " P" << (i+1) << "+ P" << (i+1) << "- " << loadR << "\n";
    }
    
    // Instantiate the magnetic subcircuit
    circuit << "\n* Magnetic component instance\n";
    circuit << "X1 ";
    for (size_t i = 0; i < numWindings; ++i) {
        circuit << "P" << (i+1) << "+ P" << (i+1) << "- ";
    }
    circuit << magneticRef << "\n";
    
    // Analysis commands
    size_t numCycles = 20;
    double stopTime = numCycles * period;
    double stepTime = period / 1000;  // 1000 points per cycle
    
    circuit << "\n* Transient analysis\n";
    circuit << ".tran " << stepTime << " " << stopTime << " " << (stopTime - 2*period) << "\n";
    
    // Probe all winding currents and voltages
    circuit << "\n* Probes\n";
    for (size_t i = 0; i < numWindings; ++i) {
        circuit << ".save v(P" << (i+1) << "+) v(P" << (i+1) << "-)\n";
    }
    
    circuit << "\n.end\n";
    
    return circuit.str();
}

} // namespace OpenMagnetics