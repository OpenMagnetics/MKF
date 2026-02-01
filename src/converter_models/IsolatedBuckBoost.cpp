#include "converter_models/IsolatedBuckBoost.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include <sstream>
#include <algorithm>
#include <map>
#include "support/Exceptions.h"
#include "processors/NgspiceRunner.h"

namespace OpenMagnetics {

    double IsolatedBuckBoost::calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency) {
        auto dutyCycle = outputVoltage / (inputVoltage + outputVoltage) * efficiency;
        if (dutyCycle >= 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Duty cycle must be smaller than 1");
        }
        return dutyCycle;
    }

    IsolatedBuckBoost::IsolatedBuckBoost(const json& j) {
        from_json(j, *this);
    }

    AdvancedIsolatedBuckBoost::AdvancedIsolatedBuckBoost(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint IsolatedBuckBoost::processOperatingPointsForInputVoltage(double inputVoltage, IsolatedBuckBoostOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double primaryOutputVoltage = outputOperatingPoint.get_output_voltages()[0];
        double primaryOutputCurrent = outputOperatingPoint.get_output_currents()[0];
        double diodeVoltageDrop = get_diode_voltage_drop();
        double totalReflectedSecondaryCurrent = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size() - 1; ++secondaryIndex) {
            totalReflectedSecondaryCurrent += outputOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex];
        }

        double efficiency = 1;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }


        auto dutyCycle = calculate_duty_cycle(inputVoltage, primaryOutputVoltage, efficiency);

        auto tOn = dutyCycle / switchingFrequency;

        auto primaryCurrentPeakToPeak = (inputVoltage * primaryOutputVoltage) / (inputVoltage + primaryOutputVoltage) / (switchingFrequency * inductance);

        auto primaryVoltaveMaximum = inputVoltage;
        auto primaryVoltaveMinimum = primaryOutputVoltage - diodeVoltageDrop;
        auto primaryVoltavePeaktoPeak = primaryVoltaveMaximum - primaryVoltaveMinimum;

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, primaryOutputCurrent, 0);
            voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, 0);

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // Secondaries
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size() - 1; ++secondaryIndex) {
            Waveform currentWaveform;
            Waveform voltageWaveform;
            double secondaryOutputCurrent = outputOperatingPoint.get_output_currents()[secondaryIndex + 1];

            auto secondaryCurrentMaximum = (1 + dutyCycle) / (1 - dutyCycle) * secondaryOutputCurrent - secondaryOutputCurrent;
            auto secondaryCurrentMinimum = 0;
            auto secondaryCurrentPeakToPeak = secondaryCurrentMaximum - secondaryCurrentMinimum;

            auto secondaryVoltaveMaximum = inputVoltage / turnsRatios[secondaryIndex] - diodeVoltageDrop;
            auto secondaryVoltaveMinimum = (primaryOutputVoltage - diodeVoltageDrop) / turnsRatios[secondaryIndex] + diodeVoltageDrop;
            auto secondaryVoltavePeaktoPeak = secondaryVoltaveMaximum - secondaryVoltaveMinimum;

            currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, secondaryCurrentPeakToPeak, switchingFrequency, 1.0 - dutyCycle, secondaryOutputCurrent, 0, tOn);
            voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, secondaryVoltavePeaktoPeak, switchingFrequency, 1.0 - dutyCycle, 0, 0, tOn);

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Secondary " + std::to_string(secondaryIndex));
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    bool IsolatedBuckBoost::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        for (size_t isolatedbuckBoostOperatingPointIndex = 1; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
            if (get_operating_points()[isolatedbuckBoostOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[isolatedbuckBoostOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
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

    DesignRequirements IsolatedBuckBoost::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double efficiency = 1;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Missing both current ripple ratio and maximum swtich current");
        }

        // Turns ratio calculation
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size() - 1, 0);
        for (auto isolatedbuckBoostOperatingPoint : get_operating_points()) {
            double primaryVoltage = isolatedbuckBoostOperatingPoint.get_output_voltages()[0];
            for (size_t secondaryIndex = 0; secondaryIndex < isolatedbuckBoostOperatingPoint.get_output_voltages().size() - 1; ++secondaryIndex) {
                auto turnsRatio = primaryVoltage / (isolatedbuckBoostOperatingPoint.get_output_voltages()[secondaryIndex + 1] + get_diode_voltage_drop());
                turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatio);
            }
        }

        // Inductance calculation

        double maximumCurrentRiple = 0;
        if (get_current_ripple_ratio()) {
            double currentRippleRatio = get_current_ripple_ratio().value();
            double maximumOutputCurrent = 0;

            for (size_t isolatedbuckOperatingPointIndex = 0; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
                auto isolatedbuckOperatingPoint = get_operating_points()[isolatedbuckOperatingPointIndex];
                auto totalCurrentInPoint = isolatedbuckOperatingPoint.get_output_currents()[0];
                double totalReflectedSecondaryCurrent = 0;
                for (size_t secondaryIndex = 0; secondaryIndex < isolatedbuckOperatingPoint.get_output_currents().size() - 1; ++secondaryIndex) {
                    totalReflectedSecondaryCurrent += isolatedbuckOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex];
                }
                totalCurrentInPoint += totalReflectedSecondaryCurrent;
                maximumOutputCurrent = std::max(maximumOutputCurrent, totalCurrentInPoint);
            }

            maximumCurrentRiple = currentRippleRatio * maximumOutputCurrent;
        }

        if (get_maximum_switch_current()) {
            for (auto isolatedbuckBoostOperatingPoint : get_operating_points()) {
                auto primaryCurrent = isolatedbuckBoostOperatingPoint.get_output_currents()[0];
                double totalReflectedSecondaryCurrent = 0;
                for (size_t secondaryIndex = 0; secondaryIndex < isolatedbuckBoostOperatingPoint.get_output_currents().size() - 1; ++secondaryIndex) {
                    totalReflectedSecondaryCurrent += isolatedbuckBoostOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex];
                }
                double primaryOutputVoltage = isolatedbuckBoostOperatingPoint.get_output_voltages()[0];
                auto dutyCycle = calculate_duty_cycle(minimumInputVoltage, primaryOutputVoltage, efficiency);
                maximumCurrentRiple = get_maximum_switch_current().value() - (primaryCurrent + totalReflectedSecondaryCurrent) / (1 - dutyCycle);
            }
        }

        double maximumNeededInductance = 0;
        for (size_t isolatedbuckBoostOperatingPointIndex = 0; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
            auto isolatedbuckBoostOperatingPoint = get_operating_points()[isolatedbuckBoostOperatingPointIndex];
            double primaryVoltage = isolatedbuckBoostOperatingPoint.get_output_voltages()[0];
            auto switchingFrequency = isolatedbuckBoostOperatingPoint.get_switching_frequency();
            auto desiredInductance = primaryVoltage * maximumInputVoltage / (primaryVoltage + maximumInputVoltage) / (2 * maximumCurrentRiple* switchingFrequency);

            maximumNeededInductance = std::max(maximumNeededInductance, desiredInductance);
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        for (auto turnsRatio : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::ISOLATED_BUCK_BOOST_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> IsolatedBuckBoost::process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) {
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
            for (size_t isolatedbuckBoostOperatingPointIndex = 0; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[isolatedbuckBoostOperatingPointIndex], turnsRatios, magnetizingInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(isolatedbuckBoostOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> IsolatedBuckBoost::process_operating_points(Magnetic magnetic) {
        IsolatedBuckBoost::run_checks(_assertErrors);

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedIsolatedBuckBoost::process() {
        IsolatedBuckBoost::run_checks(_assertErrors);

        Inputs inputs;

        double maximumNeededInductance = get_desired_inductance();
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
        inductanceWithTolerance.set_nominal(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::ISOLATED_BUCK_BOOST_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t isolatedbuckBoostOperatingPointIndex = 0; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[isolatedbuckBoostOperatingPointIndex], turnsRatios, maximumNeededInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(isolatedbuckBoostOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

    std::string IsolatedBuckBoost::generate_ngspice_circuit(
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
        auto opPoint = get_operating_points()[operatingPointIndex];
        
        double switchingFrequency = opPoint.get_switching_frequency();
        double efficiency = 1.0;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }
        double dutyCycle = calculate_duty_cycle(inputVoltage, opPoint.get_output_voltages()[0], efficiency);
        
        // Number of secondaries
        size_t numSecondaries = turnsRatios.size();
        
        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        
        // Simulation: run 10x the extraction periods for settling
        int periodsToExtract = get_num_periods_to_extract();
        const int numPeriodsTotal = 10 * periodsToExtract;
        double simTime = numPeriodsTotal * period;
        double startTime = (numPeriodsTotal - periodsToExtract) * period;
        double stepTime = period / 200;
        
        circuit << "* Isolated Buck-Boost Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, " << numSecondaries << " secondaries\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // PWM Switch
        circuit << "* PWM Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 5 0 10n 10n " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=2.5 VH=0.5\n";
        circuit << "S1 vin_dc pri_p pwm_ctrl 0 SW1\n\n";
        
        // Primary current sense
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense pri_p pri_in 0\n\n";
        
        // Isolated Buck-Boost Transformer (similar to flyback, stores energy)
        circuit << "* Transformer - Primary and " << numSecondaries << " secondaries\n";
        circuit << "Lpri pri_in 0 " << std::scientific << magnetizingInductance << std::fixed << "\n";
        
        // Secondary windings (opposite polarity for buck-boost action)
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double secondaryInductance = magnetizingInductance / (turnsRatios[secIdx] * turnsRatios[secIdx]);
            circuit << "Lsec" << secIdx << " 0 sec" << secIdx << "_in " << std::scientific << secondaryInductance << std::fixed << "\n";
        }
        circuit << "\n";
        
        // Coupling coefficients
        circuit << "* Transformer coupling\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "Kps" << secIdx << " Lpri Lsec" << secIdx << " 1\n";
        }
        // Couple secondaries to each other
        for (size_t i = 0; i < numSecondaries; ++i) {
            for (size_t j = i + 1; j < numSecondaries; ++j) {
                circuit << "Kss" << i << "_" << j << " Lsec" << i << " Lsec" << j << " 1\n";
            }
        }
        circuit << "\n";
        
        // Diode model
        circuit << "* Diode model\n";
        circuit << ".model DIDEAL D(IS=1e-14 RS=1e-6)\n\n";
        
        // Output stages for each secondary
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double outputVoltage = opPoint.get_output_voltages()[secIdx];
            double outputCurrent = opPoint.get_output_currents()[secIdx];
            double loadResistance = outputVoltage / outputCurrent;
            
            circuit << "* Secondary " << secIdx << " output stage\n";
            circuit << "Dout" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_p DIDEAL\n";
            circuit << "Vsec_sense" << secIdx << " sec" << secIdx << "_p vout" << secIdx << " 0\n";
            circuit << "Cout" << secIdx << " vout" << secIdx << " 0 10u IC=" << outputVoltage << "\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << "\n\n";
        
        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(pri_in) i(Vpri_sense)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) v(vout" << secIdx << ") i(Vsec_sense" << secIdx << ")";
        }
        circuit << "\n\n";
        
        // Options for convergence
        circuit << ".options RELTOL=0.001 ABSTOL=1e-9 VNTOL=1e-6 ITL1=1000 ITL4=1000\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << ".ic v(vout" << secIdx << ")=" << opPoint.get_output_voltages()[secIdx] << "\n";
        }
        circuit << "\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> IsolatedBuckBoost::simulate_and_extract_operating_points(
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
        
        size_t numSecondaries = turnsRatios.size();
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto ibbOpPoint = get_operating_points()[opIndex];
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = ibbOpPoint.get_switching_frequency();
                
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
                
                // Build waveform mapping
                NgspiceRunner::WaveformNameMapping waveformMapping;
                
                // Primary winding
                waveformMapping.push_back({{"voltage", "pri_in"}, {"current", "vpri_sense#branch"}});
                
                // Secondary windings
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    std::string secVoltage = "sec" + std::to_string(secIdx) + "_in";
                    std::string secCurrent = "vsec_sense" + std::to_string(secIdx) + "#branch";
                    waveformMapping.push_back({{"voltage", secVoltage}, {"current", secCurrent}});
                }
                
                std::vector<std::string> windingNames;
                windingNames.push_back("Primary");
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    windingNames.push_back("Secondary " + std::to_string(secIdx));
                }
                
                std::vector<bool> flipCurrentSign(windingNames.size(), false);
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    ibbOpPoint.get_ambient_temperature(),
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

    std::vector<IsolatedBuckBoostTopologyWaveforms> IsolatedBuckBoost::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
        
        std::vector<IsolatedBuckBoostTopologyWaveforms> topologyWaveforms;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Collect input voltages to simulate
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
        
        size_t numSecondaries = turnsRatios.size();
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            double inputVoltage = inputVoltages[inputVoltageIndex];
            
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto ibbOpPoint = get_operating_points()[opIndex];
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = ibbOpPoint.get_switching_frequency();
                double efficiency = 1.0;
                if (get_efficiency()) {
                    efficiency = get_efficiency().value();
                }
                double dutyCycle = calculate_duty_cycle(inputVoltage, ibbOpPoint.get_output_voltages()[0], efficiency);
                
                // Run simulation
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = 2;  // Two periods for topology waveform visualization
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
                IsolatedBuckBoostTopologyWaveforms waveforms;
                waveforms.frequency = switchingFrequency;
                waveforms.inputVoltageValue = inputVoltage;
                waveforms.dutyCycle = dutyCycle;
                
                // Set output voltage values
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    waveforms.outputVoltageValues.push_back(ibbOpPoint.get_output_voltages()[secIdx]);
                }
                
                // Set name
                waveforms.operatingPointName = inputVoltagesNames[inputVoltageIndex] + " input";
                if (get_operating_points().size() > 1) {
                    waveforms.operatingPointName += " op. point " + std::to_string(opIndex);
                }
                
                // Extract time vector
                waveforms.time = getWaveformData("time");
                
                // Extract primary winding signals
                waveforms.primaryVoltage = getWaveformData("pri_in");
                waveforms.primaryCurrent = getWaveformData("vpri_sense#branch");
                
                // Extract secondary signals
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    waveforms.secondaryWindingVoltages.push_back(getWaveformData("sec" + std::to_string(secIdx) + "_in"));
                    waveforms.outputVoltages.push_back(getWaveformData("vout" + std::to_string(secIdx)));
                    waveforms.secondaryCurrents.push_back(getWaveformData("vsec_sense" + std::to_string(secIdx) + "#branch"));
                }
                
                topologyWaveforms.push_back(waveforms);
            }
        }
        
        return topologyWaveforms;
    }

} // namespace OpenMagnetics
