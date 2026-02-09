#include "converter_models/SingleSwitchForward.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/CircuitSimulatorInterface.h"
#include "support/Utils.h"
#include <cfloat>
#include <sstream>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double SingleSwitchForward::get_total_reflected_secondary_current(const ForwardOperatingPoint& forwardOperatingPoint, const std::vector<double>& turnsRatios, double rippleRatio) {
        double totalReflectedSecondaryCurrent = 0;

        if (turnsRatios.size() != forwardOperatingPoint.get_output_currents().size() + 1) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Turns ratios must have one more position than outputs for the demagnetization winding");
        }

        for (size_t secondaryIndex = 0; secondaryIndex < forwardOperatingPoint.get_output_currents().size(); ++secondaryIndex) {
            totalReflectedSecondaryCurrent += forwardOperatingPoint.get_output_currents()[secondaryIndex] / turnsRatios[secondaryIndex + 1] * rippleRatio;
        }

        return totalReflectedSecondaryCurrent;
    }

    double SingleSwitchForward::get_maximum_duty_cycle() {
        if (get_duty_cycle()) {
            return get_duty_cycle().value();
        }
        else {
            return 0.45;
        }
    }

    SingleSwitchForward::SingleSwitchForward(const json& j) {
        from_json(j, *this);
    }

    AdvancedSingleSwitchForward::AdvancedSingleSwitchForward(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint SingleSwitchForward::process_operating_points_for_input_voltage(double inputVoltage, const ForwardOperatingPoint& outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance, double mainOutputInductance) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double mainOutputCurrent = outputOperatingPoint.get_output_currents()[0];
        double mainOutputVoltage = outputOperatingPoint.get_output_voltages()[0];
        double mainSecondaryTurnsRatio = turnsRatios[1];
        double diodeVoltageDrop = get_diode_voltage_drop();

        auto dutyCycle = get_maximum_duty_cycle();

        // Assume CCM
        auto period = 1.0 / switchingFrequency;
        auto t1 = period / 2 * (mainOutputVoltage + diodeVoltageDrop) / (inputVoltage / mainSecondaryTurnsRatio);
        if (t1 > period / 2) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "T1 cannot be larger than period/2, wrong topology configuration");
        }

        auto magnetizationCurrent = inputVoltage * t1 / inductance;
        auto minimumPrimaryCurrent = -magnetizationCurrent / 2;
        auto maximumPrimaryCurrent = magnetizationCurrent / 2;

        std::vector<double> minimumSecondaryCurrents;
        std::vector<double> maximumSecondaryCurrents;

        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
            auto outputCurrentRipple = get_current_ripple_ratio() * outputOperatingPoint.get_output_currents()[secondaryIndex];
            auto minimumSecondaryCurrent = outputOperatingPoint.get_output_currents()[secondaryIndex] - outputCurrentRipple / 2;
            minimumSecondaryCurrents.push_back(minimumSecondaryCurrent);

            auto maximumSecondaryCurrent = outputOperatingPoint.get_output_currents()[secondaryIndex] + outputCurrentRipple / 2;
            maximumSecondaryCurrents.push_back(maximumSecondaryCurrent);
            size_t turnsRatioSecondaryIndex = 1 + secondaryIndex;  // To avoid demagnetization winding
            minimumPrimaryCurrent += minimumSecondaryCurrent / turnsRatios[turnsRatioSecondaryIndex];
            maximumPrimaryCurrent += maximumSecondaryCurrent / turnsRatios[turnsRatioSecondaryIndex];
        }


        if (minimumPrimaryCurrent < 0) { // DCM
            t1 = sqrt(2 * mainOutputCurrent * mainOutputInductance * (mainOutputVoltage + diodeVoltageDrop) / (switchingFrequency * (inputVoltage / mainSecondaryTurnsRatio - diodeVoltageDrop - mainOutputVoltage) * (inputVoltage / mainSecondaryTurnsRatio)));
            if (t1 > period / 2) {
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "T1 cannot be larger than period/2, wrong topology configuration");
            }
            minimumPrimaryCurrent = 0;
            maximumPrimaryCurrent = magnetizationCurrent;

            for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                auto outputCurrentRipple = get_current_ripple_ratio() * outputOperatingPoint.get_output_currents()[secondaryIndex];
                minimumSecondaryCurrents[secondaryIndex] = 0;
                maximumSecondaryCurrents[secondaryIndex] = outputCurrentRipple;

                size_t turnsRatioSecondaryIndex = 1 + secondaryIndex;  // To avoid demagnetization winding
                minimumPrimaryCurrent += minimumSecondaryCurrents[secondaryIndex] / turnsRatios[turnsRatioSecondaryIndex];
                maximumPrimaryCurrent += maximumSecondaryCurrents[secondaryIndex] / turnsRatios[turnsRatioSecondaryIndex];
            }
        }

        double deadTime = period - t1 * 2;
        // Primary
        {
            double primaryCurrentPeakToPeak = maximumPrimaryCurrent - minimumPrimaryCurrent;
            double primaryVoltagePeaktoPeak = 2 * inputVoltage;
            auto currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, minimumPrimaryCurrent, deadTime);
            auto voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }
        // Demagnetization winding
        {
            double primaryVoltagePeaktoPeak = 2 * inputVoltage;
            auto currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME, magnetizationCurrent, switchingFrequency, dutyCycle, minimumPrimaryCurrent, deadTime);
            auto voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Demagnetization winding");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {

            double secondaryCurrentPeakToPeak = maximumSecondaryCurrents[secondaryIndex] - minimumSecondaryCurrents[secondaryIndex];

            size_t turnsRatioSecondaryIndex = 1 + secondaryIndex;  // To avoid demagnetization winding
            double minimumSecondaryVoltage = -(inputVoltage + diodeVoltageDrop) / turnsRatios[turnsRatioSecondaryIndex];
            double maximumSecondaryVoltage = inputVoltage / turnsRatios[turnsRatioSecondaryIndex];
            double secondaryVoltagePeakToPeak = maximumSecondaryVoltage - minimumSecondaryVoltage;
            double secondaryVoltageOffset = maximumSecondaryVoltage + minimumSecondaryVoltage;

            auto currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, secondaryCurrentPeakToPeak, switchingFrequency, dutyCycle, minimumSecondaryCurrents[secondaryIndex], 0);
            auto voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, secondaryVoltagePeakToPeak, switchingFrequency, dutyCycle, secondaryVoltageOffset, deadTime);
            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Secondary " + std::to_string(secondaryIndex));
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    bool SingleSwitchForward::run_checks(bool assert) {
        return ForwardConverterUtils::run_checks_common(this, assert);
    }

    DesignRequirements SingleSwitchForward::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double diodeVoltageDrop = get_diode_voltage_drop();
        double dutyCycle = get_maximum_duty_cycle();

        // Turns ratio calculation
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size() + 1, 0);

        // Demagnetization winding has the same number of turns as primary.
        turnsRatios[0] = 1;

        for (auto forwardOperatingPoint : get_operating_points()) {
            for (size_t secondaryIndex = 0; secondaryIndex < forwardOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                double turnsRatio = maximumInputVoltage * dutyCycle / (forwardOperatingPoint.get_output_voltages()[secondaryIndex] + diodeVoltageDrop);
                turnsRatios[secondaryIndex + 1] = std::max(turnsRatios[secondaryIndex + 1], turnsRatio);
            }
        }

        // Inductance calculation
        double minimumNeededInductance = 0;
        for (auto forwardOperatingPoint : get_operating_points()) {
            double switchingFrequency = forwardOperatingPoint.get_switching_frequency();
            double totalReflectedSecondaryCurrent = get_total_reflected_secondary_current(forwardOperatingPoint, turnsRatios);

            double neededInductance = minimumInputVoltage / (switchingFrequency * totalReflectedSecondaryCurrent);
            minimumNeededInductance = std::max(minimumNeededInductance, neededInductance);
        }

        if (get_maximum_switch_current()) {
            double maximumSwitchCurrent = get_maximum_switch_current().value();
            // According to https://www.analog.com/cn/resources/technical-articles/high-frequency-forward-pull-dc-dc-converter.html
            for (auto forwardOperatingPoint : get_operating_points()) {
                double switchingFrequency = forwardOperatingPoint.get_switching_frequency();
                double totalReflectedSecondaryCurrent = get_total_reflected_secondary_current(forwardOperatingPoint, turnsRatios, 1 + get_current_ripple_ratio());

                double minimumInductance = maximumInputVoltage * dutyCycle / switchingFrequency / (maximumSwitchCurrent - totalReflectedSecondaryCurrent);
                minimumNeededInductance = std::max(minimumNeededInductance, minimumInductance);
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
        inductanceWithTolerance.set_minimum(roundFloat(minimumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        designRequirements.set_isolation_sides(
            ForwardConverterUtils::create_isolation_sides(get_operating_points()[0].get_output_currents().size(), true));
        designRequirements.set_topology(Topologies::SINGLE_SWITCH_FORWARD_CONVERTER);
        return designRequirements;
    }

    double SingleSwitchForward::get_output_inductance(double secondaryTurnsRatio, size_t outputIndex) {
        double minimumOutputInductance = 0;
        auto dutyCycle = get_maximum_duty_cycle();
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        for (auto outputOperatingPoint : get_operating_points()) {
            double outputVoltage = outputOperatingPoint.get_output_voltages()[outputIndex];
            double switchingFrequency = outputOperatingPoint.get_switching_frequency();
            auto tOn = dutyCycle / switchingFrequency;
            auto outputInductance = (maximumInputVoltage / secondaryTurnsRatio - get_diode_voltage_drop() - outputVoltage) * tOn / get_current_ripple_ratio();
            minimumOutputInductance = std::max(minimumOutputInductance, outputInductance);
        }
        return minimumOutputInductance;
    }

    std::vector<OperatingPoint> SingleSwitchForward::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        ForwardConverterUtils::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        std::vector<double> outputInductancePerSecondary;

    // Iterate over secondaries (turnsRatios[0] is demagnetization winding)
    for (size_t secondaryIndex = 0; secondaryIndex < turnsRatios.size() - 1; ++secondaryIndex) {
        outputInductancePerSecondary.push_back(get_output_inductance(turnsRatios[secondaryIndex + 1], secondaryIndex));
    }


        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];

            for (size_t forwardOperatingPointIndex = 0; forwardOperatingPointIndex < get_operating_points().size(); ++forwardOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[forwardOperatingPointIndex], turnsRatios, magnetizingInductance, outputInductancePerSecondary[0]);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(forwardOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> SingleSwitchForward::process_operating_points(Magnetic magnetic) {
        SingleSwitchForward::run_checks(_assertErrors);

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(_magnetizingInductanceModel);;  // hardcoded
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedSingleSwitchForward::process() {
        SingleSwitchForward::run_checks(_assertErrors);

        Inputs inputs;

        double minimumNeededInductance = get_desired_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();
        std::vector<double> outputInductancePerSecondary;

        if (get_desired_output_inductances()) {
            outputInductancePerSecondary = get_desired_output_inductances().value();
        }
        else {
            // For SingleSwitchForward, turnsRatios[0] is demagnetization winding, so secondary indices start at 1
            for (size_t secondaryIndex = 1; secondaryIndex < turnsRatios.size(); ++secondaryIndex) {
                // secondaryIndex - 1 because output_voltages/output_currents don't include demag winding
                auto minimumOutputInductance = get_output_inductance(turnsRatios[secondaryIndex], secondaryIndex - 1);
                outputInductancePerSecondary.push_back(minimumOutputInductance);
            }
        }

        if (turnsRatios.size() != get_operating_points()[0].get_output_currents().size() + 1) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Turns ratios must have one more position than outputs for the demagnetization winding");
        }

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        ForwardConverterUtils::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        DesignRequirements designRequirements;

        designRequirements.get_mutable_turns_ratios().clear();
        for (auto turnsRatio : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(minimumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        designRequirements.set_isolation_sides(
            ForwardConverterUtils::create_isolation_sides(get_operating_points()[0].get_output_currents().size(), true));
        designRequirements.set_topology(Topologies::SINGLE_SWITCH_FORWARD_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t forwardOperatingPointIndex = 0; forwardOperatingPointIndex < get_operating_points().size(); ++forwardOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[forwardOperatingPointIndex], turnsRatios, minimumNeededInductance, outputInductancePerSecondary[0]);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(forwardOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

    std::string SingleSwitchForward::generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {
        
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        ForwardConverterUtils::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("operatingPointIndex out of range");
        }
        
        double inputVoltage = inputVoltages[inputVoltageIndex];
        auto opPoint = get_operating_points()[operatingPointIndex];
        
        double switchingFrequency = opPoint.get_switching_frequency();
        double dutyCycle = get_maximum_duty_cycle();
        double diodeVoltageDrop = get_diode_voltage_drop();
        
        // Number of secondaries (turns ratios[0] is demagnetization winding)
        size_t numSecondaries = turnsRatios.size() - 1;
        
        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        
        // Simulation: run steady-state periods for settling, then extract the last N periods
        int periodsToExtract = get_num_periods_to_extract();
        int numSteadyStatePeriods = get_num_steady_state_periods();
        const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;  // Steady state + extraction
        double simTime = numPeriodsTotal * period;
        double startTime = numSteadyStatePeriods * period;  // Start extracting after steady state
        double stepTime = period / 200;
        
        circuit << "* Single-Switch Forward Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, " << numSecondaries << " secondaries\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // PWM Main Switch
        circuit << "* PWM Main Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 5 0 10n 10n " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=2.5 VH=0.5\n";
        circuit << "S1 vin_dc pri_p pwm_ctrl 0 SW1\n\n";
        
        // Primary current sense
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense pri_p pri_in 0\n\n";
        
        // Forward Transformer with demagnetization winding
        // Primary winding
        circuit << "* Forward Transformer\n";
        circuit << "Lpri pri_in 0 " << std::scientific << magnetizingInductance << std::fixed << "\n";
        
        // Demagnetization winding (same turns as primary for single-switch forward = 1:1)
        circuit << "Ldemag demag_in 0 " << std::scientific << magnetizingInductance << std::fixed << "\n";
        
        // Secondary windings
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double secondaryInductance = magnetizingInductance / (turnsRatios[secIdx + 1] * turnsRatios[secIdx + 1]);
            circuit << "Lsec" << secIdx << " sec" << secIdx << "_in 0 " << std::scientific << secondaryInductance << std::fixed << "\n";
        }
        
        // Couple all windings together using pairwise K statements (ngspice requires this)
        circuit << "* Coupling: All windings coupled pairwise\n";
        circuit << "Kpri_demag Lpri Ldemag 0.9999\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "Kpri_sec" << secIdx << " Lpri Lsec" << secIdx << " 0.9999\n";
            circuit << "Kdemag_sec" << secIdx << " Ldemag Lsec" << secIdx << " 0.9999\n";
        }
        // Couple secondaries to each other if there are multiple
        for (size_t i = 0; i < numSecondaries; ++i) {
            for (size_t j = i + 1; j < numSecondaries; ++j) {
                circuit << "Ksec" << i << "_sec" << j << " Lsec" << i << " Lsec" << j << " 0.9999\n";
            }
        }
        circuit << "\n";
        
        // Demagnetization diode - conducts when main switch is off to reset the core
        // The demag winding has opposite polarity through the coupling definition
        circuit << "* Demagnetization Diode and current sense\n";
        circuit << ".model DIDEAL D(IS=1e-14 RS=1e-6)\n";
        circuit << "Vdemag_sense demag_in demag_sense 0\n";
        circuit << "Ddemag demag_sense vin_dc DIDEAL\n\n";
        
        // Output stages for each secondary
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "* Secondary " << secIdx << " output stage\n";
            
            // Forward rectifier - conducts when main switch is on
            circuit << "Dfwd" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_rect DIDEAL\n";
            
            // Freewheeling diode - conducts when main switch is off
            circuit << "Dfw" << secIdx << " 0 sec" << secIdx << "_rect DIDEAL\n";
            
            // Snubber resistors for convergence
            circuit << "Rsnub_fwd" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_rect 1MEG\n";
            circuit << "Rsnub_fw" << secIdx << " 0 sec" << secIdx << "_rect 1MEG\n";
            
            // Output inductor (LC filter)
            double outputVoltage = opPoint.get_output_voltages()[secIdx];
            double outputCurrent = opPoint.get_output_currents()[secIdx];
            double outputInductance = get_output_inductance(turnsRatios[secIdx + 1], secIdx);
            
            circuit << "Vsec_sense" << secIdx << " sec" << secIdx << "_rect sec" << secIdx << "_l_in 0\n";
            circuit << "Lout" << secIdx << " sec" << secIdx << "_l_in vout" << secIdx << " " << std::scientific << outputInductance << std::fixed << "\n";
            
            double loadResistance = outputVoltage / outputCurrent;
            circuit << "Cout" << secIdx << " vout" << secIdx << " 0 100u IC=" << outputVoltage << "\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << "\n\n";
        
        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(pri_in) i(Vpri_sense) v(demag_in) i(Vdemag_sense)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) i(Vsec_sense" << secIdx << ") v(vout" << secIdx << ")";
        }
        circuit << "\n\n";
        
        // Options
        circuit << ".options RELTOL=0.001 ABSTOL=1e-9 VNTOL=1e-6 ITL1=1000 ITL4=1000\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << ".ic v(vout" << secIdx << ")=" << opPoint.get_output_voltages()[secIdx] << "\n";
        }
        circuit << "\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> SingleSwitchForward::simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        ForwardConverterUtils::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        size_t numSecondaries = turnsRatios.size() - 1;
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto forwardOpPoint = get_operating_points()[opIndex];
                
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = forwardOpPoint.get_switching_frequency();
                
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = 1;
                config.keepTempFiles = false;
                
                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Define waveform name mapping
                NgspiceRunner::WaveformNameMapping waveformMapping;
                
                // Primary winding
                waveformMapping.push_back({{"voltage", "pri_in"}, {"current", "vpri_sense#branch"}});
                
                // Demagnetization winding
                waveformMapping.push_back({{"voltage", "demag_in"}, {"current", "vdemag_sense#branch"}});
                
                // Secondary windings
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    std::string voltageName = "sec" + std::to_string(secIdx) + "_in";
                    std::string currentName = "vsec_sense" + std::to_string(secIdx) + "#branch";
                    waveformMapping.push_back({{"voltage", voltageName}, {"current", currentName}});
                }
                
                std::vector<std::string> windingNames;
                windingNames.push_back("Primary");
                windingNames.push_back("Demagnetization winding");
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    windingNames.push_back("Secondary " + std::to_string(secIdx));
                }
                
                std::vector<bool> flipCurrentSign(2 + numSecondaries, false);
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    forwardOpPoint.get_ambient_temperature(),
                    flipCurrentSign);
                
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
    
    std::vector<OperatingPoint> SingleSwitchForward::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
        // For Single Switch Forward converter, topology waveforms are the same as operating points
        // The operating point already contains all winding voltages and currents
        return simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
    }
} // namespace OpenMagnetics
