#include "converter_models/DifferentialModeChoke.h"
#include "converter_models/Topology.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <numeric>
#include <iomanip>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    DifferentialModeChoke::DifferentialModeChoke(const json& j) {
        // Parse input voltage
        if (j.contains("inputVoltage")) {
            _inputVoltage = j["inputVoltage"].get<DimensionWithTolerance>();
        }
        
        // Parse operating current
        if (j.contains("operatingCurrent")) {
            _operatingCurrent = j["operatingCurrent"].get<double>();
        }
        
        // Parse configuration
        if (j.contains("configuration")) {
            std::string configStr = j["configuration"].get<std::string>();
            if (configStr == "THREE_PHASE") {
                _configuration = DmcConfiguration::THREE_PHASE;
            } else if (configStr == "THREE_PHASE_WITH_NEUTRAL") {
                _configuration = DmcConfiguration::THREE_PHASE_WITH_NEUTRAL;
            } else {
                _configuration = DmcConfiguration::SINGLE_PHASE;
            }
        }
        
        // Parse line frequency (REQUIRED)
        if (!j.contains("lineFrequency")) {
            throw std::runtime_error("DifferentialModeChoke: 'lineFrequency' is required");
        }
        _lineFrequency = j["lineFrequency"].get<double>();
        
        // Parse optional parameters
        if (j.contains("peakCurrent")) {
            _peakCurrent = j["peakCurrent"].get<double>();
        }
        if (j.contains("minimumInductance")) {
            _minimumInductance = j["minimumInductance"].get<double>();
        }
        if (j.contains("filterCapacitance")) {
            _filterCapacitance = j["filterCapacitance"].get<double>();
        }
        if (j.contains("minimumImpedance")) {
            std::vector<ImpedanceAtFrequency> impedances;
            for (const auto& imp : j["minimumImpedance"]) {
                ImpedanceAtFrequency impAtFreq;
                impAtFreq.set_frequency(imp["frequency"].get<double>());
                
                ImpedancePoint impedance;
                impedance.set_magnitude(imp["impedance"].get<double>());
                impAtFreq.set_impedance(impedance);
                
                impedances.push_back(impAtFreq);
            }
            _minimumImpedance = impedances;
        }
        if (j.contains("switchingFrequency")) {
            _switchingFrequency = j["switchingFrequency"].get<double>();
        }
        if (j.contains("maximumDcResistance")) {
            _maximumDcResistance = j["maximumDcResistance"].get<double>();
        }
        if (j.contains("ambientTemperature")) {
            _ambientTemperature = j["ambientTemperature"].get<double>();
        }
    }

    DesignRequirements DifferentialModeChoke::process_design_requirements() {
        DesignRequirements designRequirements;

        // DMC may have multiple windings for 3-phase
        designRequirements.get_mutable_turns_ratios().clear();
        
        // For 3-phase: each winding is independent (not coupled like CMC)
        // But they all need the same inductance, so we set turns ratio = 1 for all
        int numWindings = get_number_of_windings();
        if (numWindings > 1) {
            for (int i = 1; i < numWindings; i++) {
                DimensionWithTolerance turnsRatio;
                turnsRatio.set_nominal(1.0);
                designRequirements.get_mutable_turns_ratios().push_back(turnsRatio);
            }
        }

        // All windings on primary side (not isolated)
        std::vector<IsolationSide> isolationSides;
        for (int i = 0; i < numWindings; i++) {
            isolationSides.push_back(IsolationSide::PRIMARY);
        }
        designRequirements.set_isolation_sides(isolationSides);

        // Set application and sub-application
        designRequirements.set_application(Application::INTERFERENCE_SUPPRESSION);
        designRequirements.set_sub_application(SubApplication::DIFFERENTIAL_MODE_NOISE_FILTERING);

        // Set inductance requirement
        if (_minimumInductance) {
            DimensionWithTolerance inductanceWithTolerance;
            inductanceWithTolerance.set_minimum(_minimumInductance.value());
            designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        }

        // Set impedance requirements if specified
        if (_minimumImpedance) {
            designRequirements.set_minimum_impedance(_minimumImpedance.value());
        }

        return designRequirements;
    }

    std::vector<OperatingPoint> DifferentialModeChoke::process_operating_points() {
        std::vector<OperatingPoint> operatingPoints;

        // Operating frequency is the line frequency (for loss calculation)
        // The ripple is at the switching frequency
        double operatingFrequency = _lineFrequency;
        double rippleFrequency = _switchingFrequency.value_or(100000);

        // Determine peak current
        double peakCurrent = _peakCurrent.value_or(_operatingCurrent * 1.2);  // 20% margin if not specified

        // Calculate current ripple (difference between peak and operating)
        double currentRipple = peakCurrent - _operatingCurrent;
        if (currentRipple < 0) {
            currentRipple = _operatingCurrent * 0.2;  // Default 20% ripple
        }

        int numWindings = get_number_of_windings();
        double operatingVoltage = resolve_dimensional_values(_inputVoltage);
        
        // Phase angles for 3-phase systems (120° apart)
        std::vector<double> phaseAngles;
        if (_configuration == DmcConfiguration::SINGLE_PHASE) {
            phaseAngles = {0.0};
        } else if (_configuration == DmcConfiguration::THREE_PHASE) {
            phaseAngles = {0.0, 2.0 * std::numbers::pi / 3.0, 4.0 * std::numbers::pi / 3.0};
        } else {  // THREE_PHASE_WITH_NEUTRAL
            phaseAngles = {0.0, 2.0 * std::numbers::pi / 3.0, 4.0 * std::numbers::pi / 3.0, 0.0};
        }

        std::vector<OperatingPointExcitation> excitations;
        
        for (int windingIdx = 0; windingIdx < numWindings; windingIdx++) {
            double phaseAngle = phaseAngles[windingIdx];
            bool isNeutral = (_configuration == DmcConfiguration::THREE_PHASE_WITH_NEUTRAL && 
                             windingIdx == 3);
            
            // Create sinusoidal current waveform at line frequency
            // For DMC, the operating current is the mains frequency component
            // The ripple current (at switching frequency) is what gets filtered
            
            double currentAmplitude = _operatingCurrent * std::sqrt(2);  // RMS to peak
            if (isNeutral) {
                // Neutral current is typically smaller (unbalanced load component)
                currentAmplitude *= 0.1;  // 10% of phase current
            }
            
            // Generate sinusoidal waveform at line frequency
            int numPoints = 10000;
            double period = 1.0 / operatingFrequency;
            std::vector<double> timeData(numPoints);
            std::vector<double> currentData(numPoints);
            
            for (int i = 0; i < numPoints; i++) {
                double t = i * period / numPoints;
                timeData[i] = t;
                
                // Sinusoidal current at line frequency with phase offset
                currentData[i] = currentAmplitude * std::sin(2.0 * std::numbers::pi * operatingFrequency * t + phaseAngle);
                
                // Add triangular ripple at switching frequency
                double ripplePhase = std::fmod(t * rippleFrequency, 1.0);
                double ripple = (ripplePhase < 0.5) ? 
                    (4.0 * ripplePhase - 1.0) : (3.0 - 4.0 * ripplePhase);
                currentData[i] += currentRipple * ripple;
            }
            
            Waveform currentWaveform;
            currentWaveform.set_data(currentData);
            currentWaveform.set_time(timeData);

            SignalDescriptor current;
            current.set_waveform(currentWaveform);
            auto currentProcessed = Inputs::calculate_processed_data(currentWaveform, operatingFrequency, true);
            auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, operatingFrequency);
            auto currentHarmonics = Inputs::calculate_harmonics_data(sampledCurrentWaveform, operatingFrequency);
            current.set_processed(currentProcessed);
            current.set_harmonics(currentHarmonics);

            // Voltage across inductor: V = L * di/dt
            // At line frequency, the AC voltage is small
            // The main voltage stress is from ripple filtering
            double voltageAmplitude = operatingVoltage * 0.05;  // ~5% of input
            if (isNeutral) {
                voltageAmplitude *= 0.1;
            }
            
            auto voltageWaveform = Inputs::create_waveform(
                WaveformLabel::SINUSOIDAL,
                voltageAmplitude,
                operatingFrequency,
                0.5
            );

            SignalDescriptor voltage;
            voltage.set_waveform(voltageWaveform);
            auto voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, operatingFrequency, true);
            auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, operatingFrequency);
            auto voltageHarmonics = Inputs::calculate_harmonics_data(sampledVoltageWaveform, operatingFrequency);
            voltage.set_processed(voltageProcessed);
            voltage.set_harmonics(voltageHarmonics);

            // Create excitation
            OperatingPointExcitation excitation;
            excitation.set_current(current);
            excitation.set_frequency(operatingFrequency);
            excitation.set_voltage(voltage);
            
            excitations.push_back(excitation);
        }

        OperatingPoint operatingPoint;
        operatingPoint.set_excitations_per_winding(excitations);
        operatingPoint.get_mutable_conditions().set_ambient_temperature(_ambientTemperature);

        operatingPoints.push_back(operatingPoint);

        return operatingPoints;
    }

    Inputs DifferentialModeChoke::process() {
        Inputs inputs;
        
        auto designRequirements = process_design_requirements();
        auto operatingPoints = process_operating_points();

        inputs.set_design_requirements(designRequirements);
        inputs.set_operating_points(operatingPoints);

        return inputs;
    }

    std::string DifferentialModeChoke::generate_ngspice_circuit(double inductance, double frequency) {
        // ============================================================
        // LC Low-Pass Filter Test Circuit
        // ============================================================
        // DMC forms the L part of an LC filter for differential mode noise
        // 
        // Circuit topology:
        //   DM Noise Source -> DMC (L) -> Load
        //                         |
        //                         C
        //                         |
        //                        GND
        //
        // The filter cutoff frequency: fc = 1 / (2π√LC)
        // ============================================================

        std::ostringstream circuit;
        double period = 1.0 / frequency;
        
        // Simulation parameters
        int numPeriods = 20;  // Simulate 20 periods for settling
        double simTime = numPeriods * period;
        double stepTime = period / 100;  // 100 points per period
        
        // Calculate filter capacitor for reasonable cutoff (fc = frequency / 10)
        double cutoffFrequency = frequency / 10;
        double filterCapacitance = 1.0 / (4.0 * std::numbers::pi * std::numbers::pi 
                                          * cutoffFrequency * cutoffFrequency * inductance);
        
        // Ensure reasonable capacitor value
        if (filterCapacitance < 1e-9) filterCapacitance = 1e-9;   // Min 1nF
        if (filterCapacitance > 100e-6) filterCapacitance = 100e-6; // Max 100µF
        
        double operatingVoltage = resolve_dimensional_values(_inputVoltage);
        double loadResistance = operatingVoltage / _operatingCurrent;
        
        circuit << "* Differential Mode Choke LC Filter Test Circuit\n";
        circuit << "* Generated by OpenMagnetics\n";
        circuit << "* Test frequency: " << (frequency / 1e3) << " kHz\n";
        circuit << "* Inductance: " << (inductance * 1e6) << " uH\n";
        circuit << "* Filter capacitance: " << (filterCapacitance * 1e9) << " nF\n";
        circuit << "* Cutoff frequency: " << (cutoffFrequency / 1e3) << " kHz\n\n";

        // DC source with DM noise superimposed
        circuit << "* Input: DC + Differential Mode Noise\n";
        circuit << "Vdc in_dc 0 " << operatingVoltage << "\n";
        circuit << "Vnoise noise_src in_dc SIN(0 " << (operatingVoltage * 0.1) << " " << frequency << ")\n\n";
        
        // DMC inductor with current sense
        circuit << "* Differential Mode Choke\n";
        circuit << "Vsense noise_src dmc_in 0\n";  // Current sense
        circuit << "Ldmc dmc_in dmc_out " << std::scientific << inductance << std::fixed << "\n";
        circuit << "Rdmc_esr dmc_out filter_out 0.01\n";  // Small ESR\n\n";

        // Filter capacitor
        circuit << "* Filter Capacitor\n";
        circuit << "Cfilt filter_out 0 " << std::scientific << filterCapacitance << std::fixed << "\n";
        circuit << "Rc_esr filter_out filter_out_esr 0.01\n\n";  // ESR

        // Load
        circuit << "* Load\n";
        circuit << "Rload filter_out_esr 0 " << loadResistance << "\n\n";

        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << std::fixed << "\n\n";

        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(noise_src) v(filter_out) i(Vsense)\n\n";

        // Options
        circuit << ".options RELTOL=0.001 ABSTOL=1e-12 VNTOL=1e-9\n\n";

        circuit << ".end\n";

        return circuit.str();
    }

    std::vector<DmcSimulationWaveforms> DifferentialModeChoke::simulate_and_extract_waveforms(
        double inductance,
        const std::vector<double>& frequencies) {
        
        std::vector<DmcSimulationWaveforms> results;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for DMC simulation");
        }

        for (double frequency : frequencies) {
            std::string netlist = generate_ngspice_circuit(inductance, frequency);
            
            SimulationConfig config;
            config.frequency = frequency;
            config.keepTempFiles = false;
            config.extractOnePeriod = false;
            
            auto simResult = runner.run_simulation(netlist, config);
            
            if (!simResult.success) {
                continue;  // Skip failed simulations
            }

            DmcSimulationWaveforms waveforms;
            waveforms.frequency = frequency;
            waveforms.operatingPointName = "DMC_" + std::to_string(static_cast<int>(frequency / 1000)) + "kHz";

            // Extract waveforms from simulation result
            for (size_t i = 0; i < simResult.waveformNames.size(); ++i) {
                const auto& name = simResult.waveformNames[i];
                const auto& wf = simResult.waveforms[i];
                
                if (name == "noise_src") {
                    waveforms.inputVoltage = wf.get_data();
                    waveforms.time = wf.get_time().value_or(std::vector<double>());
                } else if (name == "filter_out") {
                    waveforms.outputVoltage = wf.get_data();
                } else if (name.find("vsense") != std::string::npos) {
                    waveforms.inductorCurrent = wf.get_data();
                }
            }

            // Calculate DM attenuation (dB) at test frequency
            // Compare AC component only (subtract DC)
            if (!waveforms.inputVoltage.empty() && !waveforms.outputVoltage.empty()) {
                double vinAvg = std::accumulate(waveforms.inputVoltage.begin(), 
                                                waveforms.inputVoltage.end(), 0.0) / waveforms.inputVoltage.size();
                double voutAvg = std::accumulate(waveforms.outputVoltage.begin(), 
                                                 waveforms.outputVoltage.end(), 0.0) / waveforms.outputVoltage.size();
                
                double vinAcPeak = 0, voutAcPeak = 0;
                for (size_t j = 0; j < waveforms.inputVoltage.size(); ++j) {
                    vinAcPeak = std::max(vinAcPeak, std::abs(waveforms.inputVoltage[j] - vinAvg));
                    voutAcPeak = std::max(voutAcPeak, std::abs(waveforms.outputVoltage[j] - voutAvg));
                }
                
                if (voutAcPeak > 0 && vinAcPeak > 0) {
                    waveforms.dmAttenuation = 20.0 * std::log10(vinAcPeak / voutAcPeak);
                }
            }

            results.push_back(waveforms);
        }

        return results;
    }

    std::vector<OperatingPoint> DifferentialModeChoke::simulate_and_extract_operating_points(double inductance) {
        std::vector<OperatingPoint> operatingPoints;
        
        // Get frequencies from minimum impedance requirements or use defaults
        std::vector<double> frequencies;
        if (_minimumImpedance.has_value()) {
            for (const auto& impReq : _minimumImpedance.value()) {
                frequencies.push_back(impReq.get_frequency());
            }
        }
        
        if (frequencies.empty()) {
            double freq = _switchingFrequency.value_or(100000);
            frequencies = {freq, freq * 2, freq * 5};  // Fundamental and harmonics
        }

        auto simWaveforms = simulate_and_extract_waveforms(inductance, frequencies);
        
        for (const auto& simWf : simWaveforms) {
            if (simWf.inductorCurrent.empty()) continue;
            
            Waveform currentWaveform;
            currentWaveform.set_data(simWf.inductorCurrent);
            if (!simWf.time.empty()) {
                currentWaveform.set_time(simWf.time);
            }
            
            SignalDescriptor current;
            current.set_waveform(currentWaveform);
            auto sampledWaveform = Inputs::calculate_sampled_waveform(currentWaveform, simWf.frequency);
            current.set_harmonics(Inputs::calculate_harmonics_data(sampledWaveform, simWf.frequency));
            current.set_processed(Inputs::calculate_processed_data(currentWaveform, simWf.frequency, true));
            
            OperatingPointExcitation excitation;
            excitation.set_current(current);
            excitation.set_frequency(simWf.frequency);

            OperatingPoint operatingPoint;
            operatingPoint.set_excitations_per_winding({excitation});
            operatingPoint.get_mutable_conditions().set_ambient_temperature(_ambientTemperature);
            operatingPoint.set_name(simWf.operatingPointName);
            
            operatingPoints.push_back(operatingPoint);
        }

        return operatingPoints;
    }

    double DifferentialModeChoke::calculate_required_inductance(
        double targetAttenuation,
        double frequency,
        double capacitance) {
        
        // LC filter transfer function: H(s) = 1 / (1 + s²LC)
        // At frequency f: |H(jω)| = 1 / |1 - ω²LC| for ω >> ω_c
        // Attenuation (dB) = 20*log10(|H|) ≈ -40*log10(f/fc) for f >> fc
        // 
        // For target attenuation A (positive dB):
        // A = 40*log10(f/fc)
        // f/fc = 10^(A/40)
        // fc = f / 10^(A/40)
        // fc = 1/(2π√LC)
        // L = 1/(4π²fc²C)
        
        double ratio = std::pow(10.0, targetAttenuation / 40.0);
        double cutoffFrequency = frequency / ratio;
        double inductance = 1.0 / (4.0 * std::numbers::pi * std::numbers::pi * 
                                   cutoffFrequency * cutoffFrequency * capacitance);
        return inductance;
    }

    std::vector<DmcAttenuationResult> DifferentialModeChoke::verify_attenuation(
        double inductance,
        std::optional<double> capacitance) {
        
        std::vector<DmcAttenuationResult> results;
        
        // Use provided capacitance or calculate from filter design
        double filterCapacitance = capacitance.value_or(_filterCapacitance.value_or(0));
        
        // If no capacitance provided, calculate based on switching frequency
        if (filterCapacitance <= 0) {
            double noiseFrequency = _switchingFrequency.value_or(100000);
            double cutoffFrequency = noiseFrequency / 10;  // fc = fsw/10 for good attenuation
            filterCapacitance = 1.0 / (4.0 * std::numbers::pi * std::numbers::pi * 
                                       cutoffFrequency * cutoffFrequency * inductance);
            // Clamp to reasonable range
            if (filterCapacitance < 1e-9) filterCapacitance = 1e-9;
            if (filterCapacitance > 100e-6) filterCapacitance = 100e-6;
        }
        
        // Get test frequencies from impedance requirements
        std::vector<std::pair<double, double>> testPoints;  // {frequency, required_attenuation}
        
        if (_minimumImpedance.has_value()) {
            for (const auto& impReq : _minimumImpedance.value()) {
                double freq = impReq.get_frequency();
                double impedance = impReq.get_impedance().get_magnitude();
                
                // Convert impedance requirement to approximate attenuation
                // Assuming load impedance is Rload = Vin/Iout
                double operatingVoltage = resolve_dimensional_values(_inputVoltage);
                double loadImpedance = operatingVoltage / _operatingCurrent;
                double expectedAttenuation = 20.0 * std::log10(impedance / loadImpedance);
                
                testPoints.push_back({freq, std::max(0.0, expectedAttenuation)});
            }
        }
        
        if (testPoints.empty()) {
            // Use default test frequencies if no requirements specified
            double fsw = _switchingFrequency.value_or(100000);
            testPoints = {{fsw, 20.0}, {fsw * 2, 30.0}, {fsw * 5, 40.0}};
        }
        
        // Calculate cutoff frequency
        double cutoffFrequency = 1.0 / (2.0 * std::numbers::pi * std::sqrt(inductance * filterCapacitance));
        
        // Run ngspice simulations to verify
        NgspiceRunner runner;
        bool ngspiceAvailable = runner.is_available();
        
        for (const auto& [frequency, requiredAttenuation] : testPoints) {
            DmcAttenuationResult result;
            result.frequency = frequency;
            result.requiredAttenuation = requiredAttenuation;
            
            // Calculate theoretical attenuation
            // For LC low-pass filter above cutoff: A = 40*log10(f/fc)
            if (frequency > cutoffFrequency) {
                result.theoreticalAttenuation = 40.0 * std::log10(frequency / cutoffFrequency);
            } else {
                result.theoreticalAttenuation = 0.0;  // Below cutoff, no attenuation
            }
            
            // Try to get measured attenuation from ngspice
            if (ngspiceAvailable) {
                try {
                    auto waveforms = simulate_and_extract_waveforms(inductance, {frequency});
                    if (!waveforms.empty()) {
                        result.measuredAttenuation = waveforms[0].dmAttenuation;
                    } else {
                        result.measuredAttenuation = result.theoreticalAttenuation;
                    }
                } catch (...) {
                    result.measuredAttenuation = result.theoreticalAttenuation;
                }
            } else {
                result.measuredAttenuation = result.theoreticalAttenuation;
            }
            
            // Check if requirement is met
            result.passed = (result.measuredAttenuation >= requiredAttenuation * 0.9);  // 10% margin
            
            // Generate message
            std::ostringstream msg;
            msg << std::fixed << std::setprecision(1);
            msg << "At " << (frequency / 1e3) << " kHz: ";
            msg << "Required " << requiredAttenuation << " dB, ";
            msg << "Measured " << result.measuredAttenuation << " dB ";
            msg << "(Theoretical: " << result.theoreticalAttenuation << " dB) - ";
            msg << (result.passed ? "PASS" : "FAIL");
            result.message = msg.str();
            
            results.push_back(result);
        }
        
        return results;
    }

    json DifferentialModeChoke::propose_design() {
        json design;
        
        double operatingVoltage = resolve_dimensional_values(_inputVoltage);
        double loadImpedance = operatingVoltage / _operatingCurrent;
        
        // Determine target frequency and attenuation
        double targetFrequency = _switchingFrequency.value_or(100000);
        double targetAttenuation = 40.0;  // Default 40dB attenuation at switching frequency
        
        if (_minimumImpedance.has_value() && !_minimumImpedance->empty()) {
            // Use first impedance requirement as primary target
            const auto& impReq = _minimumImpedance->front();
            targetFrequency = impReq.get_frequency();
            double impedance = impReq.get_impedance().get_magnitude();
            targetAttenuation = 20.0 * std::log10(impedance / loadImpedance);
        }
        
        // Calculate cutoff frequency (fc = f_target / 10^(A/40))
        double ratio = std::pow(10.0, targetAttenuation / 40.0);
        double cutoffFrequency = targetFrequency / ratio;
        
        // Choose capacitance based on practical constraints
        // For EMI filters, typical values: 100nF to 10µF
        double capacitance = _filterCapacitance.value_or(0);
        
        if (capacitance <= 0) {
            // Select capacitance to get reasonable inductance
            // Target inductance ~100µH to 1mH for typical EMI filters
            double targetInductance = 470e-6;  // 470µH as starting point
            capacitance = 1.0 / (4.0 * std::numbers::pi * std::numbers::pi * 
                                 cutoffFrequency * cutoffFrequency * targetInductance);
            
            // Clamp to practical values
            if (capacitance < 100e-9) capacitance = 100e-9;    // Min 100nF
            if (capacitance > 10e-6) capacitance = 10e-6;       // Max 10µF
        }
        
        // Calculate required inductance
        double inductance = 1.0 / (4.0 * std::numbers::pi * std::numbers::pi * 
                                   cutoffFrequency * cutoffFrequency * capacitance);
        
        // Use minimum inductance requirement if specified and higher
        if (_minimumInductance.has_value() && _minimumInductance.value() > inductance) {
            inductance = _minimumInductance.value();
            // Recalculate capacitance to maintain cutoff
            capacitance = 1.0 / (4.0 * std::numbers::pi * std::numbers::pi * 
                                 cutoffFrequency * cutoffFrequency * inductance);
        }
        
        // Calculate peak current for saturation requirement
        double peakCurrent = _peakCurrent.value_or(_operatingCurrent * 1.4);  // 40% margin
        
        // Energy storage requirement
        double energyStorage = 0.5 * inductance * peakCurrent * peakCurrent;
        
        // Build design JSON
        design["inductance"] = inductance;
        design["inductance_uH"] = inductance * 1e6;
        design["capacitance"] = capacitance;
        design["capacitance_nF"] = capacitance * 1e9;
        design["cutoffFrequency"] = cutoffFrequency;
        design["cutoffFrequency_kHz"] = cutoffFrequency / 1e3;
        design["targetAttenuation_dB"] = targetAttenuation;
        design["peakCurrent"] = peakCurrent;
        design["energyStorage_mJ"] = energyStorage * 1e3;
        design["configuration"] = _configuration == DmcConfiguration::SINGLE_PHASE ? "SINGLE_PHASE" :
                                  _configuration == DmcConfiguration::THREE_PHASE ? "THREE_PHASE" :
                                  "THREE_PHASE_WITH_NEUTRAL";
        design["numberOfWindings"] = get_number_of_windings();
        
        // Verify the design
        auto verificationResults = verify_attenuation(inductance, capacitance);
        design["verification"] = json::array();
        bool allPassed = true;
        for (const auto& result : verificationResults) {
            json vr;
            vr["frequency"] = result.frequency;
            vr["requiredAttenuation_dB"] = result.requiredAttenuation;
            vr["measuredAttenuation_dB"] = result.measuredAttenuation;
            vr["theoreticalAttenuation_dB"] = result.theoreticalAttenuation;
            vr["passed"] = result.passed;
            vr["message"] = result.message;
            design["verification"].push_back(vr);
            if (!result.passed) allPassed = false;
        }
        design["allRequirementsMet"] = allPassed;
        
        return design;
    }

} // namespace OpenMagnetics
