#include "converter_models/Forward.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>

namespace OpenMagnetics {

    double get_total_reflected_secondary_current(ForwardOperatingPoint forwardOperatingPoint, std::vector<double> turnsRatios, double rippleRatio=1) {
        double totalReflectedSecondaryCurrent = 0;

        if (turnsRatios.size() != forwardOperatingPoint.get_output_currents().size() + 1) {
            throw std::runtime_error("Turns ratios must have one more position than outputs for the demagnetization winding");
        }

        for (size_t secondaryIndex = 0; secondaryIndex < forwardOperatingPoint.get_output_currents().size(); ++secondaryIndex) {
            totalReflectedSecondaryCurrent += forwardOperatingPoint.get_output_currents()[secondaryIndex] / turnsRatios[secondaryIndex + 1] * rippleRatio;
        }

        return totalReflectedSecondaryCurrent;
    }

    double Forward::get_maximum_duty_cycle() {
        if (get_duty_cycle()) {
            return get_duty_cycle().value();
        }
        else {
            return 0.45;
        }
    }

    Forward::Forward(const json& j) {
        from_json(j, *this);
    }

    AdvancedForward::AdvancedForward(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint Forward::process_operating_points_for_input_voltage(double inputVoltage, ForwardOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, double mainOutputInductance) {

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
            throw std::runtime_error("T1 cannot be larger than period/2, wrong topology configuration");
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
                throw std::runtime_error("T1 cannot be larger than period/2, wrong topology configuration");
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
            double primaryVoltavePeaktoPeak = 2 * inputVoltage;
            auto currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, minimumPrimaryCurrent, deadTime);
            auto voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }
        // Demagnetization winding
        {
            double primaryVoltavePeaktoPeak = 2 * inputVoltage;
            auto currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME, magnetizationCurrent, switchingFrequency, dutyCycle, minimumPrimaryCurrent, deadTime);
            auto voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
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

    bool Forward::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw std::runtime_error("At least one operating point is needed");
        }
        for (size_t forwardOperatingPointIndex = 1; forwardOperatingPointIndex < get_operating_points().size(); ++forwardOperatingPointIndex) {
            if (get_operating_points()[forwardOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw std::runtime_error("Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[forwardOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
                if (!assert) {
                    return false;
                }
                throw std::runtime_error("Different operating points cannot have different number of output currents");
            }
        }
        if (!get_input_voltage().get_nominal() && !get_input_voltage().get_maximum() && !get_input_voltage().get_minimum()) {
            if (!assert) {
                return false;
            }
            throw std::runtime_error("No input voltage introduced");
        }

        return true;
    }

    DesignRequirements Forward::process_design_requirements() {
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
                turnsRatios[secondaryIndex + 1] = std::max(turnsRatios[secondaryIndex], turnsRatio);
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
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0)); // For primary
        isolationSides.push_back(get_isolation_side_from_index(0)); // For demagnetization winding
        for (size_t windingIndex = 0; windingIndex < get_operating_points()[0].get_output_currents().size(); ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex + 1));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::SINGLE_SWITCH_FORWARD_CONVERTER);
        return designRequirements;
    }

    double Forward::get_output_inductance(double secondaryTurnsRatio, size_t outputIndex) {
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

    std::vector<OperatingPoint> Forward::process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) {
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

        std::vector<double> outputInductancePerSecondary;

        for (size_t forwardOperatingPointIndex = 0; forwardOperatingPointIndex < get_operating_points().size(); ++forwardOperatingPointIndex) {
            outputInductancePerSecondary.push_back(get_output_inductance(turnsRatios[forwardOperatingPointIndex + 1], forwardOperatingPointIndex));
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

    std::vector<OperatingPoint> Forward::process_operating_points(Magnetic magnetic) {
        Forward::run_checks(_assertErrors);

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedForward::process() {
        Forward::run_checks(_assertErrors);

        Inputs inputs;

        double minimumNeededInductance = get_desired_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();
        std::vector<double> outputInductancePerSecondary;

        if (get_desired_output_inductances()) {
            outputInductancePerSecondary = get_desired_output_inductances().value();
        }
        else {
            for (size_t secondaryIndex = 0; secondaryIndex < turnsRatios.size(); ++secondaryIndex) {
                auto minimumOutputInductance = get_output_inductance(turnsRatios[secondaryIndex], secondaryIndex);
                outputInductancePerSecondary.push_back(minimumOutputInductance);
            }
        }

        if (turnsRatios.size() != get_operating_points()[0].get_output_currents().size() + 1) {
            throw std::runtime_error("Turns ratios must have one more position than outputs for the demagnetization winding");
        }

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
        inductanceWithTolerance.set_nominal(roundFloat(minimumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0)); // For primary
        isolationSides.push_back(get_isolation_side_from_index(0)); // For demagnetization winding
        for (size_t windingIndex = 0; windingIndex < get_operating_points()[0].get_output_currents().size(); ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex + 1));
        }
        designRequirements.set_isolation_sides(isolationSides);
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
} // namespace OpenMagnetics
