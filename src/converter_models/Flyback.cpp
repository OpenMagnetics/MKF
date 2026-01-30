#include "converter_models/Flyback.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/CircuitSimulatorInterface.h"
#include "support/Utils.h"
#include <cfloat>
#include <algorithm>
#include <map>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    Flyback::Flyback(const json& j) {
        from_json(j, *this);
    }

    AdvancedFlyback::AdvancedFlyback(const json& j) {
        from_json(j, *this);
    }

    double calculate_BMO_duty_cycle(double outputVoltage, double inputVoltage, double turnsRatio){
        return (turnsRatio * outputVoltage) / (inputVoltage + turnsRatio * outputVoltage);
    }

    double calculate_BMO_primary_current_peak(double outputCurrent, double efficiency, double dutyCycle, double turnsRatio){
        return (2 * outputCurrent) / (efficiency * (1 - dutyCycle) * turnsRatio);
    }


    double calculate_QRM_frequency(double magnetizingInductance, double outputPower, double outputVoltage, double minimumInputVoltage, double turnsRatio, 
                                   double diodeVoltageDrop, double efficiency, double drainSourceCapacitance = 100e-12) {

        double dt = std::numbers::pi * sqrt(magnetizingInductance * drainSourceCapacitance);
        double a = pow((outputVoltage + diodeVoltageDrop + 1.0 / turnsRatio * minimumInputVoltage),2);
        double b = efficiency * minimumInputVoltage * minimumInputVoltage * pow((outputVoltage + diodeVoltageDrop),2);
        double c = outputVoltage + diodeVoltageDrop + 1.0 / turnsRatio * minimumInputVoltage;
        double d = sqrt(outputPower/(efficiency * magnetizingInductance));
        double e = minimumInputVoltage * (outputVoltage + diodeVoltageDrop);
        double f = sqrt(4 * dt + (2 * magnetizingInductance * outputPower* a) / b);
        double g = (1.414 * magnetizingInductance * c * d) / e;
        double h = 4.0 / pow((f + g),2);

        return h;
    }


    FlybackModes FlybackOperatingPoint::resolve_mode(std::optional<double> currentRippleRatio) {
        if (get_mode()) {
            return get_mode().value();
        }
        else {
            if (!currentRippleRatio) {
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Either current ripple ratio or mode is needed for the Flyback OperatingPoint Mode");
            }
            auto mode = currentRippleRatio.value() < 1? FlybackModes::CONTINUOUS_CONDUCTION_MODE : FlybackModes::DISCONTINUOUS_CONDUCTION_MODE;
            return mode;
        }
    }
    double FlybackOperatingPoint::resolve_switching_frequency(double inputVoltage, double diodeVoltageDrop, std::optional<double> inductance, std::optional<std::vector<double>> turnsRatios, double efficiency) {
        if (get_switching_frequency()) {
            return get_switching_frequency().value();
        }
        else {
            if (!get_mode()) {
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Either switching frequency or mode is needed for the Flyback OperatingPoint");
            }
            auto mode = get_mode().value();
            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "Switching Frequency is needed for CCM");
                }
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "Switching Frequency is needed for DCM");
                }
                case FlybackModes::QUASI_RESONANT_MODE: {
                    if (!inductance) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "Inductance in missing for switching frequency calculation");
                    }
                    if (!turnsRatios) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "TurnsRatios in missing for switching frequency calculation");
                    }
                    double totalOutputVoltageReflectedPrimaryMinusDiode = 0;

                    for (size_t secondaryIndex = 0; secondaryIndex < get_output_voltages().size(); ++secondaryIndex) {
                        auto outputVoltage = get_output_voltages()[secondaryIndex];
                        auto turnsRatio = turnsRatios.value()[secondaryIndex];
                        totalOutputVoltageReflectedPrimaryMinusDiode += outputVoltage * turnsRatio;
                    }
                    
                    double totalOutputPower = Flyback::get_total_input_power(get_output_currents(), get_output_voltages(), 1, 0);

                    double switchingFrequency = calculate_QRM_frequency(inductance.value(), totalOutputPower, totalOutputVoltageReflectedPrimaryMinusDiode / turnsRatios.value()[0], inputVoltage, turnsRatios.value()[0], diodeVoltageDrop, efficiency);
                    return switchingFrequency;
                }
                case FlybackModes::BOUNDARY_MODE_OPERATION: {
                    if (!inductance) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "Inductance in missing for switching frequency calculation");
                    }
                    if (!turnsRatios) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "TurnsRatios in missing for switching frequency calculation");
                    }
                    double currentPeak = 0;
                    double switchingFrequency = 0;
                    for (size_t secondaryIndex = 0; secondaryIndex < get_output_voltages().size(); ++secondaryIndex) {
                        auto outputCurrent = get_output_currents()[secondaryIndex];
                        auto outputVoltage = get_output_voltages()[secondaryIndex];
                        auto turnsRatio = turnsRatios.value()[secondaryIndex];
                        auto dutyCycleMaximum = calculate_BMO_duty_cycle((outputVoltage + diodeVoltageDrop), outputVoltage, turnsRatio);
                        currentPeak = std::max(currentPeak, calculate_BMO_primary_current_peak(outputCurrent, efficiency, dutyCycleMaximum, turnsRatio)); // hardcoded
                        double tOn = (currentPeak * inductance.value()) / inputVoltage;
                        double tOff = (currentPeak * inductance.value()) / (turnsRatio * outputVoltage);

                        switchingFrequency = std::max(switchingFrequency, 1.0 / (tOn + tOff));
                    }
                    return switchingFrequency;
                }
                default:
                  throw InvalidInputException(ErrorCode::INVALID_INPUT, "Unknown mode in Flyback");
            }
        }
    }

    OperatingPoint Flyback::process_operating_points_for_input_voltage(double inputVoltage, FlybackOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, std::optional<FlybackModes> customMode, std::optional<double> customDutyCycle, std::optional<double> customDeadTime) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.resolve_switching_frequency(inputVoltage, get_diode_voltage_drop(), inductance, turnsRatios, get_efficiency());

        double deadTime = 0;
        double maximumReflectedOutputVoltage = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
            double outputVoltage = outputOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop();
            maximumReflectedOutputVoltage = std::max(maximumReflectedOutputVoltage, outputVoltage * turnsRatios[secondaryIndex]);
        }

        double primaryVoltavePeaktoPeak = inputVoltage + maximumReflectedOutputVoltage;

        double totalOutputPower = get_total_input_power(outputOperatingPoint.get_output_currents(), outputOperatingPoint.get_output_voltages(), 1, 0);
        double maximumEffectiveLoadCurrent = totalOutputPower / outputOperatingPoint.get_output_voltages()[0];
        double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
        double totalInputPower = get_total_input_power(outputOperatingPoint.get_output_currents(), outputOperatingPoint.get_output_voltages(), get_efficiency(), 0);
        double averageInputCurrent = totalInputPower / inputVoltage;

        double dutyCycle;
        if (customDutyCycle) {
            dutyCycle = customDutyCycle.value();
        }
        else {
            dutyCycle = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
        }

        double centerSecondaryCurrentRampLumped = maximumEffectiveLoadCurrent / (1 - dutyCycle);
        double centerPrimaryCurrentRamp = centerSecondaryCurrentRampLumped / turnsRatios[0];


        if (customDeadTime) {
            deadTime = customDeadTime.value();
        }

        if (dutyCycle > 1 ) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "dutyCycle cannot be larger than one: " + std::to_string(dutyCycle));
        }

        double primaryCurrentAverage = centerPrimaryCurrentRamp;
        double currentRippleRatio;
        if (std::isnan(get_current_ripple_ratio())) {
            double primaryCurrentPeakToPeak = inputVoltage * dutyCycle / switchingFrequency / inductance;
            currentRippleRatio = primaryCurrentPeakToPeak / centerPrimaryCurrentRamp;
        }
        else {
            currentRippleRatio = get_current_ripple_ratio();
        }
        double primaryCurrentPeakToPeak = centerPrimaryCurrentRamp * currentRippleRatio * 2;
        double primaryCurrentOffset = primaryCurrentAverage - primaryCurrentPeakToPeak / 2;
        primaryCurrentOffset = std::max(0.0, primaryCurrentOffset);

        FlybackModes mode;
        if (customMode) {
            mode = customMode.value();
        }
        else {
            if (primaryCurrentOffset > 0) {
                mode = FlybackModes::CONTINUOUS_CONDUCTION_MODE;
            }
            else {
                mode = FlybackModes::DISCONTINUOUS_CONDUCTION_MODE;
            }
        }

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, primaryCurrentOffset, deadTime);

            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    break;
                }
                case FlybackModes::QUASI_RESONANT_MODE:
                case FlybackModes::BOUNDARY_MODE_OPERATION:
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    break;
                }
            }

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // Secondaries
        for (size_t secondaryIndex = 0; secondaryIndex < turnsRatios.size(); ++secondaryIndex) {

            Waveform currentWaveform;
            Waveform voltageWaveform;

            double secondaryPower = get_total_input_power(outputOperatingPoint.get_output_currents()[secondaryIndex], outputOperatingPoint.get_output_voltages()[secondaryIndex], 1, 0);
            double powerDivider = secondaryPower / totalOutputPower;

            double minimumSecondaryVoltage = -inputVoltage / turnsRatios[secondaryIndex];
            double maximumSecondaryVoltage = outputOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop();
            double secondaryVoltagePeaktoPeak = maximumSecondaryVoltage - minimumSecondaryVoltage;
            double secondaryCurrentAverage = centerPrimaryCurrentRamp * turnsRatios[secondaryIndex] * powerDivider;
            double secondaryCurrentPeaktoPeak = secondaryCurrentAverage * currentRippleRatio * 2;
            double secondaryCurrentOffset = std::max(0.0, secondaryCurrentAverage - secondaryCurrentPeaktoPeak / 2);

            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::SECONDARY_RECTANGULAR, secondaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    break;
                }
                case FlybackModes::QUASI_RESONANT_MODE:
                case FlybackModes::BOUNDARY_MODE_OPERATION:
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME, secondaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    break;
                }
            }

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Secondary " + std::to_string(secondaryIndex));
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    std::string Flyback::generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {
        
        // Get input voltages
        std::vector<double> inputVoltages;
        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
        }
        
        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("operatingPointIndex out of range");
        }
        
        double inputVoltage = inputVoltages[inputVoltageIndex];
        auto opPoint = get_mutable_operating_points()[operatingPointIndex];  // Copy to allow modification
        
        // Calculate switching frequency and duty cycle
        double switchingFrequency = opPoint.resolve_switching_frequency(inputVoltage, get_diode_voltage_drop(), magnetizingInductance, turnsRatios, get_efficiency());
        
        // Calculate duty cycle
        double totalOutputPower = get_total_input_power(opPoint.get_output_currents(), opPoint.get_output_voltages(), 1, 0);
        double maximumEffectiveLoadCurrent = totalOutputPower / opPoint.get_output_voltages()[0];
        double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
        double totalInputPower = get_total_input_power(opPoint.get_output_currents(), opPoint.get_output_voltages(), get_efficiency(), 0);
        double averageInputCurrent = totalInputPower / inputVoltage;
        double dutyCycle = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
        
        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        double secondaryInductance = magnetizingInductance / (turnsRatios[0] * turnsRatios[0]);
        
        // Simulation: 10 periods, skip first 5 for steady state
        double simTime = 10 * period;
        double startTime = 5 * period;
        double stepTime = period / 200;
        
        circuit << "* Flyback Converter - Generated from MAS::Flyback\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << "%\n";
        circuit << "* Lp=" << (magnetizingInductance*1e6) << "uH, N=" << turnsRatios[0] << "\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // PWM Switch (ideal: zero on-resistance, infinite off-resistance)
        circuit << "* PWM Switch (ideal)\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 5 0 1p 1p " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW(Vt=2.5 Vh=0.5 Ron=1u Roff=1G)\n";
        circuit << "S1 vin_dc pri_p pwm_ctrl 0 SW1\n\n";
        
        // Primary current sense (for waveform extraction)
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense pri_p pri_in 0\n\n";
        
        // Flyback Transformer (ideal coupling = 1)
        circuit << "* Flyback Transformer (ideal coupled inductors, k=1)\n";
        circuit << "Lp pri_in 0 " << magnetizingInductance << " IC=0\n";
        circuit << "Ls sec_in 0 " << secondaryInductance << " IC=0\n";
        circuit << "Kfly Lp Ls -1\n\n";
        
        // Output Rectifier (ideal diode)
        circuit << "* Output Rectifier (ideal diode)\n";
        circuit << ".model DIDEAL D(Is=1e-14 Rs=1u N=0.001)\n";
        circuit << "Dout sec_in sec_p DIDEAL\n\n";
        
        // Secondary current sense
        circuit << "* Secondary current sense\n";
        circuit << "Vsec_sense sec_p vout 0\n\n";
        
        // Output capacitor and load
        double loadResistance = opPoint.get_output_voltages()[0] / opPoint.get_output_currents()[0];
        circuit << "* Output Filter and Load\n";
        circuit << "Cout vout 0 10u IC=" << opPoint.get_output_voltages()[0] << "\n";
        circuit << "Rload vout 0 " << loadResistance << "\n\n";
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << stepTime << " " << simTime << " " << startTime << "\n\n";
        
        // Save primary and secondary voltages and currents
        // For secondary, save both sec_in (winding voltage) and vout (output voltage)
        circuit << "* Output signals (voltage and current for each winding)\n";
        circuit << ".save v(pri_in) v(sec_in) v(vout) i(Vpri_sense) i(Vsec_sense)\n\n";
        
        // Options for convergence
        circuit << ".ic v(vout)=" << opPoint.get_output_voltages()[0] << "\n";
        circuit << ".options RELTOL=0.001 ABSTOL=1n VNTOL=1u\n\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> Flyback::simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
            inputVoltagesNames.push_back("Nom.");
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
            inputVoltagesNames.push_back("Min.");
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
            inputVoltagesNames.push_back("Max.");
        }
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            double inputVoltage = inputVoltages[inputVoltageIndex];
            
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto flybackOpPoint = get_mutable_operating_points()[opIndex];  // Copy to allow modification
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                // Calculate switching frequency for waveform extraction
                double switchingFrequency = flybackOpPoint.resolve_switching_frequency(
                    inputVoltage, get_diode_voltage_drop(), magnetizingInductance, turnsRatios, get_efficiency());
                
                // Run simulation
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.keepTempFiles = false;
                
                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Define waveform name mapping for flyback circuit
                // Primary: v(pri_in) for voltage, i(vpri_sense) for current
                // Secondary: v(sec_in) for WINDING voltage (not v(vout) which is DC output), i(vsec_sense) for current
                NgspiceRunner::WaveformNameMapping waveformMapping = {
                    {{"voltage", "v(pri_in)"}, {"current", "i(vpri_sense)"}},
                    {{"voltage", "v(sec_in)"}, {"current", "i(vsec_sense)"}}
                };
                
                std::vector<std::string> windingNames = {"Primary", "Secondary 0"};
                std::vector<bool> flipCurrentSign = {false, false};  // Keep current signs as-is
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    flybackOpPoint.get_ambient_temperature(),
                    flipCurrentSign);
                
                // Set name
                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt. (simulated)";
                if (get_operating_points().size() > 1) {
                    name += " op. point " + std::to_string(opIndex);
                }
                operatingPoint.set_name(name);
                
                operatingPoints.push_back(operatingPoint);
            }
        }
        
        return operatingPoints;
    }

    std::string Flyback::generate_ngspice_circuit_with_magnetic(
        const Magnetic& magneticConst,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {
        
        // Make a copy since some methods are non-const
        Magnetic magnetic = magneticConst;
        
        // Get input voltages
        std::vector<double> inputVoltages;
        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
        }
        
        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("operatingPointIndex out of range");
        }
        
        double inputVoltage = inputVoltages[inputVoltageIndex];
        auto opPoint = get_mutable_operating_points()[operatingPointIndex];
        
        // Extract turns ratios and magnetizing inductance from the Magnetic
        auto coil = magnetic.get_coil();
        size_t numWindings = coil.get_functional_description().size();
        if (numWindings < 2) {
            throw std::invalid_argument("Magnetic must have at least 2 windings for flyback");
        }
        
        double primaryTurns = coil.get_functional_description()[0].get_number_turns();
        std::vector<double> turnsRatios;
        for (size_t i = 1; i < numWindings; ++i) {
            double secondaryTurns = coil.get_functional_description()[i].get_number_turns();
            turnsRatios.push_back(primaryTurns / secondaryTurns);
        }
        
        double magnetizingInductance = resolve_dimensional_values(
            MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
        
        // Calculate switching frequency and duty cycle
        double switchingFrequency = opPoint.resolve_switching_frequency(
            inputVoltage, get_diode_voltage_drop(), magnetizingInductance, turnsRatios, get_efficiency());
        
        double totalOutputPower = get_total_input_power(opPoint.get_output_currents(), opPoint.get_output_voltages(), 1, 0);
        double maximumEffectiveLoadCurrent = totalOutputPower / opPoint.get_output_voltages()[0];
        double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
        double totalInputPower = get_total_input_power(opPoint.get_output_currents(), opPoint.get_output_voltages(), get_efficiency(), 0);
        double averageInputCurrent = totalInputPower / inputVoltage;
        double dutyCycle = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
        
        // Generate the magnetic subcircuit using CircuitSimulatorExporter
        CircuitSimulatorExporterNgspiceModel ngspiceExporter;
        std::string magneticSubcircuit = ngspiceExporter.export_magnetic_as_subcircuit(
            magnetic, switchingFrequency, opPoint.get_ambient_temperature());
        
        // Build the main circuit netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        
        // Simulation: 10 periods, skip first 5 for steady state
        double simTime = 10 * period;
        double startTime = 5 * period;
        double stepTime = period / 200;
        
        circuit << "* Flyback Converter with Real Magnetic Component - Generated from MAS::Flyback\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << "%\n";
        circuit << "* Magnetic: " << magnetic.get_reference() << "\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, N=" << turnsRatios[0] << "\n\n";
        
        // Include the magnetic subcircuit definition
        circuit << "* Magnetic Component Subcircuit\n";
        circuit << magneticSubcircuit << "\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // PWM Switch (ideal: zero on-resistance, infinite off-resistance)
        circuit << "* PWM Switch (ideal)\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 5 0 1p 1p " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW(Vt=2.5 Vh=0.5 Ron=1u Roff=1G)\n";
        circuit << "S1 vin_dc pri_p pwm_ctrl 0 SW1\n\n";
        
        // Primary current sense (for waveform extraction)
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense pri_p pri_in 0\n\n";
        
        // Instantiate the magnetic component subcircuit
        // Subcircuit pins are: P1+ P1- P2+ P2- ... for each winding
        // We connect: pri_in to P1+, 0 to P1-, sec_in to P2+, 0 to P2-
        std::string subcktName = fix_filename(magnetic.get_reference());
        circuit << "* Magnetic component instance\n";
        circuit << "X1 pri_in 0 sec_in 0 " << subcktName << "\n\n";
        
        // Output Rectifier (ideal diode)
        circuit << "* Output Rectifier (ideal diode)\n";
        circuit << ".model DIDEAL D(Is=1e-14 Rs=1u N=0.001)\n";
        circuit << "Dout sec_in sec_p DIDEAL\n\n";
        
        // Secondary current sense
        circuit << "* Secondary current sense\n";
        circuit << "Vsec_sense sec_p vout 0\n\n";
        
        // Output capacitor and load
        double loadResistance = opPoint.get_output_voltages()[0] / opPoint.get_output_currents()[0];
        circuit << "* Output Filter and Load\n";
        circuit << "Cout vout 0 10u IC=" << opPoint.get_output_voltages()[0] << "\n";
        circuit << "Rload vout 0 " << loadResistance << "\n\n";
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << stepTime << " " << simTime << " " << startTime << "\n\n";
        
        // Save primary and secondary voltages and currents
        circuit << "* Output signals (voltage and current for each winding)\n";
        circuit << ".save v(pri_in) v(sec_in) v(vout) i(Vpri_sense) i(Vsec_sense)\n\n";
        
        // Options for convergence
        circuit << ".ic v(vout)=" << opPoint.get_output_voltages()[0] << "\n";
        circuit << ".options RELTOL=0.001 ABSTOL=1n VNTOL=1u\n\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }

    std::vector<OperatingPoint> Flyback::simulate_with_magnetic_and_extract_operating_points(
        const Magnetic& magneticConst) {
        
        // Make a copy since some methods are non-const
        Magnetic magnetic = magneticConst;
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Extract turns ratios and magnetizing inductance from the Magnetic
        auto coil = magnetic.get_coil();
        size_t numWindings = coil.get_functional_description().size();
        if (numWindings < 2) {
            throw std::invalid_argument("Magnetic must have at least 2 windings for flyback");
        }
        
        double primaryTurns = coil.get_functional_description()[0].get_number_turns();
        std::vector<double> turnsRatios;
        for (size_t i = 1; i < numWindings; ++i) {
            double secondaryTurns = coil.get_functional_description()[i].get_number_turns();
            turnsRatios.push_back(primaryTurns / secondaryTurns);
        }
        
        double magnetizingInductance = resolve_dimensional_values(
            MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
            inputVoltagesNames.push_back("Nom.");
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
            inputVoltagesNames.push_back("Min.");
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
            inputVoltagesNames.push_back("Max.");
        }
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            double inputVoltage = inputVoltages[inputVoltageIndex];
            
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto flybackOpPoint = get_mutable_operating_points()[opIndex];
                
                // Generate circuit with real magnetic
                std::string netlist = generate_ngspice_circuit_with_magnetic(magnetic, inputVoltageIndex, opIndex);
                
                // Calculate switching frequency for waveform extraction
                double switchingFrequency = flybackOpPoint.resolve_switching_frequency(
                    inputVoltage, get_diode_voltage_drop(), magnetizingInductance, turnsRatios, get_efficiency());
                
                // Run simulation
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.keepTempFiles = false;
                
                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Define waveform name mapping for flyback circuit
                NgspiceRunner::WaveformNameMapping waveformMapping = {
                    {{"voltage", "v(pri_in)"}, {"current", "i(vpri_sense)"}},
                    {{"voltage", "v(sec_in)"}, {"current", "i(vsec_sense)"}}
                };
                
                std::vector<std::string> windingNames = {"Primary", "Secondary 0"};
                std::vector<bool> flipCurrentSign = {false, false};
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    flybackOpPoint.get_ambient_temperature(),
                    flipCurrentSign);
                
                // Set name
                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt. (simulated with " + magnetic.get_reference() + ")";
                if (get_operating_points().size() > 1) {
                    name += " op. point " + std::to_string(opIndex);
                }
                operatingPoint.set_name(name);
                
                operatingPoints.push_back(operatingPoint);
            }
        }
        
        return operatingPoints;
    }

    std::vector<FlybackTopologyWaveforms> Flyback::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
        
        std::vector<FlybackTopologyWaveforms> topologyWaveforms;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Collect input voltages to simulate
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        inputVoltages.push_back(get_input_voltage().get_nominal().value());
        inputVoltagesNames.push_back("Nom.");
        
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
            inputVoltagesNames.push_back("Min.");
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
            inputVoltagesNames.push_back("Max.");
        }
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            double inputVoltage = inputVoltages[inputVoltageIndex];
            
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto flybackOpPoint = get_mutable_operating_points()[opIndex];
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                // Calculate switching frequency
                double switchingFrequency = flybackOpPoint.resolve_switching_frequency(
                    inputVoltage, get_diode_voltage_drop(), magnetizingInductance, turnsRatios, get_efficiency());
                
                // Calculate duty cycle (same formula as in generate_ngspice_circuit)
                double totalOutputPower = get_total_input_power(flybackOpPoint.get_output_currents(), flybackOpPoint.get_output_voltages(), 1, 0);
                double maximumEffectiveLoadCurrent = totalOutputPower / flybackOpPoint.get_output_voltages()[0];
                double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
                double totalInputPower = get_total_input_power(flybackOpPoint.get_output_currents(), flybackOpPoint.get_output_voltages(), get_efficiency(), 0);
                double averageInputCurrent = totalInputPower / inputVoltage;
                double dutyCycle = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
                
                // Run simulation
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.keepTempFiles = false;
                
                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Build name-to-index map for waveform lookup (case-insensitive)
                std::map<std::string, size_t> nameToIndex;
                for (size_t i = 0; i < simResult.waveformNames.size(); ++i) {
                    std::string lower = simResult.waveformNames[i];
                    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                    nameToIndex[lower] = i;
                }
                
                // Helper lambda to get waveform data by name
                auto getWaveformData = [&](const std::string& name) -> std::vector<double> {
                    std::string lower = name;
                    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                    auto it = nameToIndex.find(lower);
                    if (it != nameToIndex.end()) {
                        return simResult.waveforms[it->second].get_data();
                    }
                    return {};
                };
                
                // Extract topology waveforms
                FlybackTopologyWaveforms waveforms;
                waveforms.frequency = switchingFrequency;
                waveforms.inputVoltageValue = inputVoltage;
                waveforms.outputVoltageValue = flybackOpPoint.get_output_voltages()[0];
                waveforms.dutyCycle = dutyCycle;
                
                // Set name
                waveforms.operatingPointName = inputVoltagesNames[inputVoltageIndex] + " input";
                if (get_operating_points().size() > 1) {
                    waveforms.operatingPointName += " op. point " + std::to_string(opIndex);
                }
                
                // Extract time vector
                waveforms.time = getWaveformData("time");
                
                // Extract voltage signals
                waveforms.switchNodeVoltage = getWaveformData("v(pri_in)");
                waveforms.secondaryWindingVoltage = getWaveformData("v(sec_in)");
                waveforms.outputVoltage = getWaveformData("v(vout)");
                
                // Extract current signals
                waveforms.primaryCurrent = getWaveformData("i(vpri_sense)");
                waveforms.secondaryCurrent = getWaveformData("i(vsec_sense)");
                
                topologyWaveforms.push_back(waveforms);
            }
        }
        
        return topologyWaveforms;
    }

    double Flyback::get_total_input_power(std::vector<double> outputCurrents, std::vector<double> outputVoltages, double efficiency, double diodeVoltageDrop) {
        double totalPower = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputCurrents.size(); ++secondaryIndex) {
            totalPower += outputCurrents[secondaryIndex] * (outputVoltages[secondaryIndex] + diodeVoltageDrop);
        }

        return totalPower / efficiency;
    }


    double Flyback::get_total_input_current(std::vector<double> outputCurrents, double inputVoltage, std::vector<double> outputVoltages, double diodeVoltageDrop) {
        double totalCurrent = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputCurrents.size(); ++secondaryIndex) {
            totalCurrent += outputCurrents[secondaryIndex] * (outputVoltages[secondaryIndex] + diodeVoltageDrop) / inputVoltage;
        }

        return totalCurrent;
    }


    double Flyback::get_total_input_power(double outputCurrent, double outputVoltage, double efficiency, double diodeVoltageDrop) {
        double totalPower = outputCurrent * (outputVoltage + diodeVoltageDrop);

        return totalPower / efficiency;
    }

    double Flyback::get_minimum_output_reflected_voltage(double maximumDrainSourceVoltage, double maximumInputVoltage, double safetyMargin) {
        return maximumDrainSourceVoltage * safetyMargin - maximumInputVoltage;
    }

    bool Flyback::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        for (size_t flybackOperatingPointIndex = 1; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            if (get_operating_points()[flybackOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[flybackOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
                if (!assert) {
                    return false;
                }
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Different operating points cannot have different number of output currents");
            }
        }
        if (!get_input_voltage().get_nominal() && !get_input_voltage().get_maximum() && !get_input_voltage().get_minimum()) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "No input voltage introduced");
        }

        return true;
    }

    DesignRequirements Flyback::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_maximum_drain_source_voltage() && !get_maximum_duty_cycle()) {
            throw std::invalid_argument("Missing both maximum duty cycle and maximum drain source voltage");
        }
        double globalNeededInductance = 0;
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size(), 0);

        if (get_maximum_duty_cycle()) {
            double maximumDutyCycle = get_maximum_duty_cycle().value();
            if (maximumDutyCycle > 1 || maximumDutyCycle < 0) {
                throw std::invalid_argument("maximumDutyCycle must be between 0 and 1");
            }
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];

                double totalOutputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), 1, 0);
                double totalInputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), get_efficiency(), 0);
                double maximumEffectiveLoadCurrent = totalOutputPower / flybackOperatingPoint.get_output_voltages()[0];
                double averageInputCurrent = totalInputPower / minimumInputVoltage;
                double maximumEffectiveLoadCurrentReflected = averageInputCurrent * (1 - maximumDutyCycle) / maximumDutyCycle;

                auto turnsRatioFirstOutput = maximumEffectiveLoadCurrent / maximumEffectiveLoadCurrentReflected;
                turnsRatios[0] = std::max(turnsRatios[0], turnsRatioFirstOutput);

                for (size_t secondaryIndex = 1; secondaryIndex < flybackOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                    auto turnsRatio = turnsRatioFirstOutput * (flybackOperatingPoint.get_output_voltages()[0] + get_diode_voltage_drop()) / (flybackOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
                    turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatio);
                }
            }
        }

        if (get_maximum_drain_source_voltage()) {
            std::vector<double> turnsRatiosFromMaximumDrainSourceVoltage(get_operating_points()[0].get_output_voltages().size(), 0);
            double maximumDrainSourceVoltage = get_maximum_drain_source_voltage().value();
            auto minimumOutputReflectedVoltage = get_minimum_output_reflected_voltage(maximumDrainSourceVoltage, maximumInputVoltage);
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
                for (size_t secondaryIndex = 0; secondaryIndex < flybackOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                    auto turnsRatio = minimumOutputReflectedVoltage / (flybackOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
                    turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex] = std::max(turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex], turnsRatio);
                }
            }

            for (size_t secondaryIndex = 0; secondaryIndex < get_operating_points()[0].get_output_voltages().size(); ++secondaryIndex) {
                if (turnsRatios[secondaryIndex] > 1) {
                    turnsRatios[secondaryIndex] = std::min(turnsRatios[secondaryIndex], turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex]);
                }
                else {
                    turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex]);
                }
            }
        }

        for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
            double switchingFrequency = flybackOperatingPoint.resolve_switching_frequency(minimumInputVoltage, get_diode_voltage_drop());
            double totalOutputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), 1, 0);
            double maximumEffectiveLoadCurrent = totalOutputPower / flybackOperatingPoint.get_output_voltages()[0];
            double dutyCycle = 0;
            if (get_maximum_duty_cycle()) {
                dutyCycle = get_maximum_duty_cycle().value();
            }
            else {
                double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
                double totalInputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), get_efficiency(), 0);
                double averageInputCurrent = totalInputPower / minimumInputVoltage;
                dutyCycle = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
            }

            double centerSecondaryCurrentRampLumped = maximumEffectiveLoadCurrent / (1 - dutyCycle);
            double centerPrimaryCurrentRamp = centerSecondaryCurrentRampLumped / turnsRatios[0];
            double tOn = dutyCycle / switchingFrequency;
            double voltsSeconds = minimumInputVoltage * tOn;
            double neededInductance = voltsSeconds / get_current_ripple_ratio() / centerPrimaryCurrentRamp;
            globalNeededInductance = std::max(globalNeededInductance, neededInductance);
        }

        // Add maximum inductance to stay in DCM
        double maximumInductance = 0;
        for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
            if (flybackOperatingPoint.get_mode()) {
                if (flybackOperatingPoint.get_mode().value() != MAS::FlybackModes::CONTINUOUS_CONDUCTION_MODE) {
                    double totalOutputPower = Flyback::get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), 1, get_diode_voltage_drop());
                    double switchingFrequency = flybackOperatingPoint.resolve_switching_frequency(minimumInputVoltage, get_diode_voltage_drop());
                    double mainOutputVoltage = flybackOperatingPoint.get_output_voltages()[0];
                    double aux = (mainOutputVoltage + get_diode_voltage_drop()) * turnsRatios[0];
                    // Accoding to Switch-Mode Power Supplies, 2nd ed; Christophe Basso; page 747.
                    double maximumInductanceThisPoint = get_efficiency() * pow(minimumInputVoltage, 2) * pow(aux, 2) / (2 * totalOutputPower * switchingFrequency * (minimumInputVoltage + aux) * (aux + get_efficiency() * minimumInputVoltage));
                    maximumInductance = std::max(maximumInductance, maximumInductanceThisPoint);
                }
            }
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        for (auto turnsRatio : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(globalNeededInductance, 10));

        if (maximumInductance > 0) {
            inductanceWithTolerance.set_maximum(roundFloat(maximumInductance, 10));
        }

        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::FLYBACK_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Flyback::process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
            inputVoltagesNames.push_back("Nom.");
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
            inputVoltagesNames.push_back("Min.");
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
            inputVoltagesNames.push_back("Max.");
        }

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto mode = get_mutable_operating_points()[flybackOperatingPointIndex].resolve_mode(get_current_ripple_ratio());
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, magnetizingInductance, mode);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(flybackOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> Flyback::process_operating_points(OpenMagnetics::Magnetic magnetic) {
        run_checks(_assertErrors);

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedFlyback::process() {
        Flyback::run_checks(_assertErrors);

        Inputs inputs;

        double globalNeededInductance = get_desired_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;


        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
            inputVoltagesNames.push_back("Nom.");
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
            inputVoltagesNames.push_back("Max.");
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
            inputVoltagesNames.push_back("Min.");
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        for (auto turnsRatio : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(globalNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::FLYBACK_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                std::optional<double> customDeadTime = std::nullopt;
                if (get_desired_duty_cycle().size() <= flybackOperatingPointIndex) {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing duty cycle for flybackOperatingPointIndex: " + std::to_string(flybackOperatingPointIndex));
                }
                double customDutyCycle = get_desired_duty_cycle()[flybackOperatingPointIndex][inputVoltageIndex];

                if (get_desired_dead_time()) {
                    if (get_desired_dead_time()->size() <= flybackOperatingPointIndex) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing dead time for flybackOperatingPointIndex: " + std::to_string(flybackOperatingPointIndex));
                    }
                    customDeadTime = get_desired_dead_time().value()[flybackOperatingPointIndex];
                }

                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, globalNeededInductance, std::nullopt, customDutyCycle, customDeadTime);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(flybackOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }
} // namespace OpenMagnetics
