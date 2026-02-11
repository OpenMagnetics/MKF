#include "converter_models/CommonModeChoke.h"
#include "converter_models/Topology.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <numbers>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    CommonModeChoke::CommonModeChoke(const json& j) {
        // Parse configuration (number of phases)
        if (j.contains("configuration")) {
            std::string configStr = j["configuration"].get<std::string>();
            if (configStr == "Single Phase" || configStr == "singlePhase" || configStr == "SINGLE_PHASE") {
                _configuration = CmcConfiguration::SINGLE_PHASE;
            } else if (configStr == "Three Phase" || configStr == "threePhase" || configStr == "THREE_PHASE") {
                _configuration = CmcConfiguration::THREE_PHASE;
            } else if (configStr == "Three Phase With Neutral" || configStr == "threePhaseWithNeutral" || 
                       configStr == "THREE_PHASE_WITH_NEUTRAL" || configStr == "Three Phase + Neutral") {
                _configuration = CmcConfiguration::THREE_PHASE_WITH_NEUTRAL;
            }
        } else if (j.contains("numberOfWindings")) {
            // Alternative: specify number of windings directly
            int numWindings = j["numberOfWindings"].get<int>();
            if (numWindings == 2) {
                _configuration = CmcConfiguration::SINGLE_PHASE;
            } else if (numWindings == 3) {
                _configuration = CmcConfiguration::THREE_PHASE;
            } else if (numWindings == 4) {
                _configuration = CmcConfiguration::THREE_PHASE_WITH_NEUTRAL;
            }
        }
        
        // Parse operating voltage
        if (j.contains("operatingVoltage")) {
            _operatingVoltage = j["operatingVoltage"].get<DimensionWithTolerance>();
        }
        
        // Parse operating current (line current)
        if (j.contains("operatingCurrent")) {
            _operatingCurrent = j["operatingCurrent"].get<double>();
        }
        
        // Parse neutral current (for 4-winding config)
        if (j.contains("neutralCurrent")) {
            _neutralCurrent = j["neutralCurrent"].get<double>();
        }
        
        // Parse minimum impedance requirements
        if (j.contains("minimumImpedance")) {
            for (const auto& imp : j["minimumImpedance"]) {
                ImpedanceAtFrequency impAtFreq;
                impAtFreq.set_frequency(imp["frequency"].get<double>());
                
                ImpedancePoint impedance;
                impedance.set_magnitude(imp["impedance"].get<double>());
                impAtFreq.set_impedance(impedance);
                
                _minimumImpedance.push_back(impAtFreq);
            }
        }
        
        // Parse line frequency (mains frequency for differential mode current) - REQUIRED
        if (!j.contains("lineFrequency")) {
            throw std::runtime_error("CommonModeChoke: 'lineFrequency' is required");
        }
        _lineFrequency = j["lineFrequency"].get<double>();
        
        // Parse optional parameters
        if (j.contains("lineImpedance")) {
            _lineImpedance = j["lineImpedance"].get<double>();
        }
        if (j.contains("maximumDcResistance")) {
            _maximumDcResistance = j["maximumDcResistance"].get<double>();
        }
        if (j.contains("maximumLeakageInductance")) {
            _maximumLeakageInductance = j["maximumLeakageInductance"].get<double>();
        }
        if (j.contains("ambientTemperature")) {
            _ambientTemperature = j["ambientTemperature"].get<double>();
        }
    }

    size_t CommonModeChoke::get_number_of_windings() const {
        switch (_configuration) {
            case CmcConfiguration::SINGLE_PHASE:
                return 2;
            case CmcConfiguration::THREE_PHASE:
                return 3;
            case CmcConfiguration::THREE_PHASE_WITH_NEUTRAL:
                return 4;
            default:
                return 2;
        }
    }

    DesignRequirements CommonModeChoke::process_design_requirements() {
        DesignRequirements designRequirements;

        size_t numWindings = get_number_of_windings();

        // CMC has N identical windings (all 1:1 turns ratios)
        // For N windings, we need N-1 turns ratios (all 1:1)
        for (size_t i = 1; i < numWindings; ++i) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(1.0);
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }

        // Set minimum impedance requirements - crucial for CMC performance
        designRequirements.set_minimum_impedance(_minimumImpedance);

        // All windings are on the primary side (line side)
        std::vector<IsolationSide> isolationSides;
        for (size_t i = 0; i < numWindings; ++i) {
            isolationSides.push_back(IsolationSide::PRIMARY);
        }
        designRequirements.set_isolation_sides(isolationSides);

        // Set application and sub-application for proper adviser configuration
        designRequirements.set_application(Application::INTERFERENCE_SUPPRESSION);
        designRequirements.set_sub_application(SubApplication::COMMON_MODE_NOISE_FILTERING);

        // Magnetizing inductance is not the primary concern for CMCs,
        // but set a minimum to ensure the core provides adequate impedance
        // The actual inductance will be determined by impedance requirements
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(1e-6);  // 1 µH minimum
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        return designRequirements;
    }

    std::vector<OperatingPoint> CommonModeChoke::process_operating_points() {
        std::vector<OperatingPoint> operatingPoints;

        // The operating frequency is the LINE frequency (mains), NOT the noise frequency
        // The differential mode current flows at mains frequency (50/60 Hz)
        // This is what causes magnetic excitation and determines core/copper losses
        // The noise frequencies in minimumImpedance are filtering requirements, not excitation
        double operatingFrequency = _lineFrequency;  // 50 or 60 Hz typically

        double operatingVoltage = resolve_dimensional_values(_operatingVoltage);
        
        // Create voltage waveform (same for all windings - minimal in differential mode)
        auto voltageWaveform = Inputs::create_waveform(
            WaveformLabel::SINUSOIDAL,
            operatingVoltage * 0.01,  // Small voltage (leakage drop only)
            operatingFrequency,
            0.5
        );

        SignalDescriptor voltage;
        voltage.set_waveform(voltageWaveform);
        auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, operatingFrequency);
        voltage.set_harmonics(Inputs::calculate_harmonics_data(sampledVoltageWaveform, operatingFrequency));
        voltage.set_processed(Inputs::calculate_processed_data(voltageWaveform, operatingFrequency, true));

        std::vector<OperatingPointExcitation> excitations;

        if (_configuration == CmcConfiguration::SINGLE_PHASE) {
            // Single-phase: 2 windings with opposite phase (180° shift)
            // Line winding
            auto currentWaveform1 = Inputs::create_waveform(
                WaveformLabel::SINUSOIDAL, 
                _operatingCurrent * std::sqrt(2) * 2,  // Peak-to-peak
                operatingFrequency, 
                0.5
            );

            SignalDescriptor winding1Current;
            winding1Current.set_waveform(currentWaveform1);
            auto sampledCurrentWaveform1 = Inputs::calculate_sampled_waveform(currentWaveform1, operatingFrequency);
            winding1Current.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform1, operatingFrequency));
            winding1Current.set_processed(Inputs::calculate_processed_data(currentWaveform1, operatingFrequency, true));

            OperatingPointExcitation winding1Excitation;
            winding1Excitation.set_current(winding1Current);
            winding1Excitation.set_frequency(operatingFrequency);
            winding1Excitation.set_voltage(voltage);
            winding1Excitation.set_name("Line");
            excitations.push_back(winding1Excitation);

            // Neutral winding (180° phase shift - opposite)
            auto currentWaveform2 = Inputs::create_waveform(
                WaveformLabel::SINUSOIDAL, 
                _operatingCurrent * std::sqrt(2) * 2,
                operatingFrequency, 
                0.5
            );
            auto waveformData = currentWaveform2.get_data();
            for (auto& point : waveformData) {
                point = -point;  // 180° shift
            }
            currentWaveform2.set_data(waveformData);

            SignalDescriptor winding2Current;
            winding2Current.set_waveform(currentWaveform2);
            auto sampledCurrentWaveform2 = Inputs::calculate_sampled_waveform(currentWaveform2, operatingFrequency);
            winding2Current.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform2, operatingFrequency));
            winding2Current.set_processed(Inputs::calculate_processed_data(currentWaveform2, operatingFrequency, true));

            OperatingPointExcitation winding2Excitation;
            winding2Excitation.set_current(winding2Current);
            winding2Excitation.set_frequency(operatingFrequency);
            winding2Excitation.set_voltage(voltage);
            winding2Excitation.set_name("Neutral");
            excitations.push_back(winding2Excitation);

        } else {
            // Three-phase (with or without neutral): currents 120° apart
            // Phase angles: L1=0°, L2=120°, L3=240°
            std::vector<double> phaseAngles = {0.0, 2.0 * std::numbers::pi / 3.0, 4.0 * std::numbers::pi / 3.0};
            std::vector<std::string> phaseNames = {"L1", "L2", "L3"};

            for (size_t phase = 0; phase < 3; ++phase) {
                // Create sinusoidal current with phase shift
                size_t numPoints = 128;  // Power of 2 for FFT
                double period = 1.0 / operatingFrequency;
                std::vector<double> data;
                std::vector<double> time;
                
                double peakCurrent = _operatingCurrent * std::sqrt(2);
                for (size_t i = 0; i < numPoints; ++i) {
                    double t = i * period / (numPoints - 1);
                    double angle = 2.0 * std::numbers::pi * operatingFrequency * t + phaseAngles[phase];
                    data.push_back(peakCurrent * std::sin(angle));
                    time.push_back(t);
                }

                Waveform currentWaveform;
                currentWaveform.set_data(data);
                currentWaveform.set_time(time);

                SignalDescriptor phaseCurrent;
                phaseCurrent.set_waveform(currentWaveform);
                auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, operatingFrequency);
                phaseCurrent.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform, operatingFrequency));
                phaseCurrent.set_processed(Inputs::calculate_processed_data(currentWaveform, operatingFrequency, true));

                OperatingPointExcitation phaseExcitation;
                phaseExcitation.set_current(phaseCurrent);
                phaseExcitation.set_frequency(operatingFrequency);
                phaseExcitation.set_voltage(voltage);
                phaseExcitation.set_name(phaseNames[phase]);
                excitations.push_back(phaseExcitation);
            }

            // Add neutral winding if 4-winding configuration
            if (_configuration == CmcConfiguration::THREE_PHASE_WITH_NEUTRAL) {
                // Neutral current: In a balanced system, neutral current is zero
                // In unbalanced systems or with harmonics, neutral carries the sum of phase currents
                // For modeling, we use a smaller current representing unbalance/harmonics
                double neutralCurrentValue = _neutralCurrent.value_or(_operatingCurrent * 0.1);  // Default 10% of line
                
                auto neutralWaveform = Inputs::create_waveform(
                    WaveformLabel::SINUSOIDAL, 
                    neutralCurrentValue * std::sqrt(2) * 2,  // Peak-to-peak
                    operatingFrequency, 
                    0.5
                );

                SignalDescriptor neutralCurrent;
                neutralCurrent.set_waveform(neutralWaveform);
                auto sampledNeutralWaveform = Inputs::calculate_sampled_waveform(neutralWaveform, operatingFrequency);
                neutralCurrent.set_harmonics(Inputs::calculate_harmonics_data(sampledNeutralWaveform, operatingFrequency));
                neutralCurrent.set_processed(Inputs::calculate_processed_data(neutralWaveform, operatingFrequency, true));

                OperatingPointExcitation neutralExcitation;
                neutralExcitation.set_current(neutralCurrent);
                neutralExcitation.set_frequency(operatingFrequency);
                neutralExcitation.set_voltage(voltage);
                neutralExcitation.set_name("Neutral");
                excitations.push_back(neutralExcitation);
            }
        }

        OperatingPoint operatingPoint;
        operatingPoint.set_excitations_per_winding(excitations);
        operatingPoint.get_mutable_conditions().set_ambient_temperature(_ambientTemperature);

        operatingPoints.push_back(operatingPoint);

        return operatingPoints;
    }

    Inputs CommonModeChoke::process() {
        Inputs inputs;
        
        auto designRequirements = process_design_requirements();
        auto operatingPoints = process_operating_points();

        inputs.set_design_requirements(designRequirements);
        inputs.set_operating_points(operatingPoints);

        return inputs;
    }

    std::string CommonModeChoke::generate_ngspice_circuit(double inductance, double frequency) {
        // ============================================================
        // LISN (Line Impedance Stabilization Network) Test Circuit
        // ============================================================
        // Based on CISPR 16-1-2 standard for conducted emissions testing
        // 
        // The LISN:
        // - Provides a known impedance (50Ω) to the EUT
        // - Isolates the AC mains from measurement equipment
        // - Extracts common mode noise for measurement
        //
        // Circuit topology:
        //   CM Noise Source -> CMC -> LISN -> 50Ω Measurement
        //
        // For single-phase: Line and Neutral through CMC
        // For three-phase: L1, L2, L3 through CMC (and N if 4-winding)
        // ============================================================

        size_t numWindings = get_number_of_windings();
        
        std::ostringstream circuit;
        double period = 1.0 / frequency;
        
        // Simulation parameters
        int numPeriods = 20;  // Simulate 20 periods for settling
        double simTime = numPeriods * period;
        double stepTime = period / 100;  // 100 points per period
        
        circuit << "* Common Mode Choke EMI Test Circuit - LISN Configuration\n";
        circuit << "* Generated by OpenMagnetics\n";
        circuit << "* Configuration: " << (numWindings == 2 ? "Single Phase" : 
                                           (numWindings == 3 ? "Three Phase" : "Three Phase + Neutral")) << "\n";
        circuit << "* Test frequency: " << (frequency / 1e3) << " kHz\n";
        circuit << "* Inductance per winding: " << (inductance * 1e6) << " uH\n\n";

        // Common mode noise source (represents switching converter noise)
        // Connected to ground, couples into both lines equally
        circuit << "* Common Mode Noise Source (switching noise)\n";
        circuit << "Vcm_noise cm_src 0 SIN(0 1 " << frequency << ")\n\n";
        
        // Noise coupling capacitors (CM noise couples to both lines)
        circuit << "* CM Noise Coupling (parasitic capacitance to ground)\n";
        for (size_t w = 0; w < numWindings; ++w) {
            circuit << "Ccm" << w << " cm_src cmc_in" << w << " 100p\n";
        }
        circuit << "\n";

        // CMC coupled inductors
        // Windings are wound in the same direction - CM currents add, DM currents cancel
        circuit << "* Common Mode Choke\n";
        for (size_t w = 0; w < numWindings; ++w) {
            circuit << "Lcmc" << w << " cmc_in" << w << " cmc_out" << w << " " 
                    << std::scientific << inductance << std::fixed << "\n";
        }
        
        // Couple all windings together (perfect coupling for ideal CMC)
        circuit << "\n* CMC Coupling (k ~ 1 for tight coupling)\n";
        for (size_t i = 0; i < numWindings; ++i) {
            for (size_t j = i + 1; j < numWindings; ++j) {
                circuit << "K" << i << "_" << j << " Lcmc" << i << " Lcmc" << j << " 0.99\n";
            }
        }
        circuit << "\n";

        // Simplified LISN per CISPR 16
        // Each line has: Series L (50µH) + parallel RC to ground
        circuit << "* LISN Network (simplified CISPR 16)\n";
        for (size_t w = 0; w < numWindings; ++w) {
            circuit << "* Line " << w << " LISN\n";
            circuit << "Llisn" << w << " cmc_out" << w << " lisn_mid" << w << " 50u\n";
            circuit << "Clisn" << w << " lisn_mid" << w << " 0 1u\n";
            circuit << "Rlisn" << w << " lisn_mid" << w << " lisn_out" << w << " 5\n";
        }
        circuit << "\n";

        // 50Ω measurement resistors (standard EMI receiver input impedance)
        circuit << "* 50Ohm Measurement Point\n";
        for (size_t w = 0; w < numWindings; ++w) {
            circuit << "Rmeas" << w << " lisn_out" << w << " 0 50\n";
        }
        circuit << "\n";

        // Current sense resistors for waveform extraction
        circuit << "* Current Sense\n";
        for (size_t w = 0; w < numWindings; ++w) {
            circuit << "Vsense" << w << " cmc_in" << w << " cmc_in" << w << "_sense 0\n";
        }
        circuit << "\n";

        // AC load (represents the EUT - Equipment Under Test)
        circuit << "* AC Load (EUT)\n";
        double loadResistance = resolve_dimensional_values(_operatingVoltage) / _operatingCurrent;
        if (_configuration == CmcConfiguration::SINGLE_PHASE) {
            circuit << "Rload cmc_in0_sense cmc_in1_sense " << loadResistance << "\n";
        } else {
            // Delta-connected load for three-phase
            for (size_t w = 0; w < 3; ++w) {
                size_t next = (w + 1) % 3;
                circuit << "Rload" << w << " cmc_in" << w << "_sense cmc_in" << next << "_sense " 
                        << (loadResistance * 3) << "\n";
            }
        }
        circuit << "\n";

        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << std::fixed << "\n\n";

        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(cm_src)";
        for (size_t w = 0; w < numWindings; ++w) {
            circuit << " v(lisn_out" << w << ") i(Vsense" << w << ")";
        }
        circuit << "\n\n";

        // Options
        circuit << ".options RELTOL=0.001 ABSTOL=1e-12 VNTOL=1e-9\n\n";

        circuit << ".end\n";

        return circuit.str();
    }

    std::vector<CmcSimulationWaveforms> CommonModeChoke::simulate_and_extract_waveforms(
        double inductance,
        const std::vector<double>& frequencies) {
        
        std::vector<CmcSimulationWaveforms> results;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for CMC simulation");
        }

        size_t numWindings = get_number_of_windings();

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

            CmcSimulationWaveforms waveforms;
            waveforms.frequency = frequency;
            waveforms.operatingPointName = "CMC_" + std::to_string(static_cast<int>(frequency / 1000)) + "kHz";

            // Extract waveforms from simulation result
            for (size_t i = 0; i < simResult.waveformNames.size(); ++i) {
                const auto& name = simResult.waveformNames[i];
                const auto& wf = simResult.waveforms[i];
                
                if (name == "cm_src") {
                    waveforms.inputVoltage = wf.get_data();
                    waveforms.time = wf.get_time().value_or(std::vector<double>());
                } else if (name.find("lisn_out") != std::string::npos) {
                    waveforms.lisnVoltage = wf.get_data();
                } else if (name.find("vsense") != std::string::npos) {
                    waveforms.windingCurrents.push_back(wf.get_data());
                }
            }

            // Calculate theoretical impedance: Z = 2*pi*f*L
            waveforms.theoreticalImpedance = 2.0 * M_PI * frequency * inductance;

            // Calculate CM attenuation (dB) = 20 * log10(Vin / Vout)
            if (!waveforms.inputVoltage.empty() && !waveforms.lisnVoltage.empty()) {
                double vinPeak = *std::max_element(waveforms.inputVoltage.begin(), waveforms.inputVoltage.end());
                double voutPeak = *std::max_element(waveforms.lisnVoltage.begin(), waveforms.lisnVoltage.end());
                if (voutPeak > 0 && vinPeak > 0) {
                    waveforms.commonModeAttenuation = 20.0 * std::log10(vinPeak / voutPeak);
                }
            }

            // Calculate actual impedance from voltage/current (Z = V/I)
            // The CM current is the sum of all winding currents
            if (!waveforms.inputVoltage.empty() && !waveforms.windingCurrents.empty()) {
                double vinPeak = *std::max_element(waveforms.inputVoltage.begin(), waveforms.inputVoltage.end());
                
                // Sum CM currents from all windings
                double totalCmCurrentPeak = 0;
                for (const auto& windingCurrent : waveforms.windingCurrents) {
                    if (!windingCurrent.empty()) {
                        double iPeak = *std::max_element(windingCurrent.begin(), windingCurrent.end());
                        totalCmCurrentPeak += iPeak;
                    }
                }
                
                // Z = V / I (for CM path, all currents flow in same direction)
                if (totalCmCurrentPeak > 1e-12) {
                    waveforms.commonModeImpedance = vinPeak / totalCmCurrentPeak;
                } else {
                    // If current is negligible, use theoretical value
                    waveforms.commonModeImpedance = waveforms.theoreticalImpedance;
                }
            } else {
                waveforms.commonModeImpedance = waveforms.theoreticalImpedance;
            }

            results.push_back(waveforms);
        }

        return results;
    }

    std::vector<OperatingPoint> CommonModeChoke::simulate_and_extract_operating_points(double inductance) {
        std::vector<OperatingPoint> operatingPoints;
        
        // Get frequencies from minimum impedance requirements
        std::vector<double> frequencies;
        for (const auto& impReq : _minimumImpedance) {
            frequencies.push_back(impReq.get_frequency());
        }
        
        if (frequencies.empty()) {
            frequencies = {100000, 1000000};  // Default: 100kHz and 1MHz
        }

        auto simWaveforms = simulate_and_extract_waveforms(inductance, frequencies);
        
        for (const auto& simWf : simWaveforms) {
            std::vector<OperatingPointExcitation> excitations;
            
            for (size_t w = 0; w < simWf.windingCurrents.size(); ++w) {
                Waveform currentWaveform;
                currentWaveform.set_data(simWf.windingCurrents[w]);
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
                
                excitations.push_back(excitation);
            }

            OperatingPoint operatingPoint;
            operatingPoint.set_excitations_per_winding(excitations);
            operatingPoint.get_mutable_conditions().set_ambient_temperature(_ambientTemperature);
            operatingPoint.set_name(simWf.operatingPointName);
            
            operatingPoints.push_back(operatingPoint);
        }

        return operatingPoints;
    }

} // namespace OpenMagnetics
