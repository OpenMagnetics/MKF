#include "converter_models/Buck.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/NgspiceRunner.h"
#include "support/Utils.h"
#include <cfloat>
#include <algorithm>
#include <map>
#include <sstream>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double Buck::calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency) {
        auto dutyCycle = (outputVoltage + diodeVoltageDrop) / ((inputVoltage + diodeVoltageDrop) * efficiency);
        if (dutyCycle >= 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Duty cycle must be smaller than 1");
        }
        return dutyCycle;
    }

    Buck::Buck(const json& j) {
        from_json(j, *this);
    }

    AdvancedBuck::AdvancedBuck(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint Buck::process_operating_points_for_input_voltage(double inputVoltage, const BuckOperatingPoint& outputOperatingPoint, double inductance) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double outputVoltage = outputOperatingPoint.get_output_voltage();
        double outputCurrent = outputOperatingPoint.get_output_current();
        double diodeVoltageDrop = get_diode_voltage_drop();
        double efficiency = 1;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }

        auto dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, diodeVoltageDrop, efficiency);

        auto tOn = dutyCycle / switchingFrequency;
        auto primaryCurrentPeakToPeak = (inputVoltage - outputVoltage) * tOn / inductance;
        auto minimumCurrent = outputCurrent - primaryCurrentPeakToPeak / 2;
        auto primaryVoltageMinimum = -outputVoltage -diodeVoltageDrop;
        auto primaryVoltageMaximum = inputVoltage - outputVoltage;
        auto primaryVoltagePeaktoPeak = primaryVoltageMaximum - primaryVoltageMinimum;

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            if (minimumCurrent < 0) {
                tOn = sqrt(2 * outputCurrent * inductance * (outputVoltage + diodeVoltageDrop) / (switchingFrequency * (inputVoltage - outputVoltage) * (inputVoltage + diodeVoltageDrop)));
                auto tOff = tOn * ((inputVoltage + diodeVoltageDrop) / (outputVoltage + diodeVoltageDrop) - 1);
                auto deadTime = 1.0 / switchingFrequency - tOn - tOff;
                primaryCurrentPeakToPeak = (inputVoltage - outputVoltage) * tOn / inductance;
                outputCurrent = primaryCurrentPeakToPeak / 2;

                currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR_WITH_DEADTIME, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, outputCurrent, deadTime);
                voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);

            }
            else {
                currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, outputCurrent, 0);
                voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, 0);
            }

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    bool Buck::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        if (!get_input_voltage().get_nominal() && !get_input_voltage().get_maximum() && !get_input_voltage().get_minimum()) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "No input voltage introduced");
        }

        return true;
    }

    DesignRequirements Buck::process_design_requirements() {
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Missing both current ripple ratio and maximum switch current");
        }

        double maximumCurrentRiple = 0;
        if (get_current_ripple_ratio()) {
            double currentRippleRatio = get_current_ripple_ratio().value();
            double maximumOutputCurrent = 0;

            for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
                auto buckOperatingPoint = get_operating_points()[buckOperatingPointIndex];
                maximumOutputCurrent = std::max(maximumOutputCurrent, buckOperatingPoint.get_output_current());
            }

            maximumCurrentRiple = currentRippleRatio * maximumOutputCurrent;
        }

        if (get_maximum_switch_current()) {
            double maximumSwitchCurrent = get_maximum_switch_current().value();
            double maximumOutputCurrent = 0;

            for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
                auto buckOperatingPoint = get_operating_points()[buckOperatingPointIndex];
                maximumOutputCurrent = std::max(maximumOutputCurrent, buckOperatingPoint.get_output_current());
            }

            maximumCurrentRiple = (maximumSwitchCurrent - maximumOutputCurrent) * 2;
        }

        double maximumNeededInductance = 0;
        for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
            auto buckOperatingPoint = get_operating_points()[buckOperatingPointIndex];
            auto switchingFrequency = buckOperatingPoint.get_switching_frequency();
            auto outputVoltage = buckOperatingPoint.get_output_voltage();
            auto desiredInductance = outputVoltage * (maximumInputVoltage - outputVoltage) / (maximumCurrentRiple * switchingFrequency * maximumInputVoltage);
            maximumNeededInductance = std::max(maximumNeededInductance, desiredInductance);
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::BUCK_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Buck::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[buckOperatingPointIndex], magnetizingInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(buckOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> Buck::process_operating_points(OpenMagnetics::Magnetic magnetic) {
        run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedBuck::process() {
        Buck::run_checks(_assertErrors);

        Inputs inputs;

        double maximumNeededInductance = get_desired_inductance();

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;


        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        DesignRequirements designRequirements;

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::BUCK_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[buckOperatingPointIndex], maximumNeededInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(buckOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

    std::string Buck::generate_ngspice_circuit(
        double inductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames_;
    collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames_);
        
        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("operatingPointIndex out of range");
        }
        
        double inputVoltage = inputVoltages[inputVoltageIndex];
        auto opPoint = get_operating_points()[operatingPointIndex];
        
        double outputVoltage = opPoint.get_output_voltage();
        double outputCurrent = opPoint.get_output_current();
        double switchingFrequency = opPoint.get_switching_frequency();
        double diodeVoltageDrop = get_diode_voltage_drop();
        double efficiency = 1.0;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }
        
        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, diodeVoltageDrop, efficiency);
        
        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        
        // Simulation: run 10x the extraction periods for settling, then extract the last N periods
        int periodsToExtract = get_num_periods_to_extract();
    int numSteadyStatePeriods = get_num_steady_state_periods();
    const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;
    double simTime = numPeriodsTotal * period;
    double startTime = numSteadyStatePeriods * period;
        double stepTime = period / 200;
        
        circuit << "* Buck Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vout=" << outputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* L=" << (inductance*1e6) << "uH, Iout=" << outputCurrent << "A\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // PWM High-side Switch (ideal)
        circuit << "* PWM High-side Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 5 0 10n 10n " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=2.5 VH=0.5\n";
        circuit << "S1 vin_dc sw pwm_ctrl 0 SW1\n\n";
        
        // Freewheeling Diode (low-side)
        circuit << "* Freewheeling Diode\n";
        circuit << ".model DIDEAL D(IS=1e-14 RS=1e-6)\n";
        circuit << "D1 0 sw DIDEAL\n\n";
        
        // Inductor with current sense
        circuit << "* Inductor with current sense\n";
        circuit << "Vl_sense sw l_in 0\n";
        circuit << "L1 l_in vout " << std::scientific << inductance << std::fixed << "\n\n";
        
        // Output capacitor and load
        double loadResistance = outputVoltage / outputCurrent;
        circuit << "* Output Filter and Load\n";
        circuit << "Cout vout 0 100u IC=" << outputVoltage << "\n";
        circuit << "Rload vout 0 " << loadResistance << "\n\n";
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << "\n\n";
        
        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(sw) v(l_in) v(vout) i(Vl_sense)\n\n";
        
        // Options for convergence
        circuit << ".options RELTOL=0.001 ABSTOL=1e-9 VNTOL=1e-6 ITL1=1000 ITL4=1000\n";
        circuit << ".ic v(vout)=" << outputVoltage << "\n\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> Buck::simulate_and_extract_operating_points(double inductance) {
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            double inputVoltage = inputVoltages[inputVoltageIndex];
            
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto buckOpPoint = get_operating_points()[opIndex];
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(inductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = buckOpPoint.get_switching_frequency();
                
                // Run simulation
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = 1;  // Only one period for operating points
                config.keepTempFiles = false;
                
                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Buck has only one winding (the inductor)
                // The inductor voltage is the difference between switch node and output
                NgspiceRunner::WaveformNameMapping waveformMapping = {
                    {{"voltage", "sw"}, {"current", "vl_sense#branch"}}
                };
                
                std::vector<std::string> windingNames = {"Primary"};
                std::vector<bool> flipCurrentSign = {false};
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    buckOpPoint.get_ambient_temperature(),
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

    std::vector<ConverterWaveforms> Buck::simulate_and_extract_topology_waveforms(double inductance) {
    
    std::vector<ConverterWaveforms> results;
    
    NgspiceRunner runner;
    if (!runner.is_available()) {
        throw std::runtime_error("ngspice is not available for simulation");
    }
    
    std::vector<double> inputVoltages;
    std::vector<std::string> inputVoltagesNames;
    collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
    
    for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
        for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
            auto opPoint = get_operating_points()[opIndex];
            
            std::string netlist = generate_ngspice_circuit(inductance, inputVoltageIndex, opIndex);
            double switchingFrequency = opPoint.get_switching_frequency();
            
            SimulationConfig config;
            config.frequency = switchingFrequency;
            config.extractOnePeriod = true;
            config.numberOfPeriods = 2;
            config.keepTempFiles = false;
            
            auto simResult = runner.run_simulation(netlist, config);
            if (!simResult.success) {
                throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
            }
            
            std::map<std::string, size_t> nameToIndex;
            for (size_t i = 0; i < simResult.waveformNames.size(); ++i) {
                std::string lower = simResult.waveformNames[i];
                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                nameToIndex[lower] = i;
            }
            auto getWaveform = [&](const std::string& name) -> Waveform {
                std::string lower = name;
                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                auto it = nameToIndex.find(lower);
                if (it != nameToIndex.end()) return simResult.waveforms[it->second];
                return Waveform();
            };
            
            ConverterWaveforms wf;
            wf.set_switching_frequency(switchingFrequency);
            std::string name = inputVoltagesNames[inputVoltageIndex] + " input";
            if (get_operating_points().size() > 1) {
                name += " op. point " + std::to_string(opIndex);
            }
            wf.set_operating_point_name(name);
            
            wf.set_input_voltage(getWaveform("sw"));
            wf.set_input_current(getWaveform("vl_sense#branch"));
            
            wf.get_mutable_output_voltages().push_back(getWaveform("vout"));
            wf.get_mutable_output_currents().push_back(getWaveform("vl_sense#branch"));
            
            results.push_back(wf);
        }
    }
    
    return results;
}

} // namespace OpenMagnetics
