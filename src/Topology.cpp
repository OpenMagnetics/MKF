#include "Topology.h"
#include "Utils.h"
#include <cfloat>

namespace OpenMagnetics {

    Flyback::Flyback(const json& j) {
        from_json(j, *this);
    }

    AdvancedFlyback::AdvancedFlyback(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint Flyback::processOperatingPointsForInputVoltage(double inputVoltage, Flyback::FlybackOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, std::optional<Flyback::Modes> customMode, std::optional<double> customDutyCycle, std::optional<double> customDeadTime) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();

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
            throw std::runtime_error("dutyCycle cannot be larger than one: " + std::to_string(dutyCycle));
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

        Flyback::Modes mode;
        if (customMode) {
            mode = customMode.value();
        }
        else {
            if (primaryCurrentOffset > 0) {
                mode = Flyback::Modes::ContinuousCurrentMode;
            }
            else {
                mode = Flyback::Modes::DiscontinuousCurrentMode;
            }
        }

        // Primary
        {
            OperatingPointExcitation excitation;
            Waveform currentWaveform;
            Waveform voltageWaveform;
            Processed currentProcessed;
            Processed voltageProcessed;

            currentWaveform = InputsWrapper::create_waveform(WaveformLabel::FLYBACK_PRIMARY, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, primaryCurrentOffset, deadTime);
            currentProcessed.set_label(WaveformLabel::FLYBACK_PRIMARY);
            currentProcessed.set_peak_to_peak(primaryCurrentPeakToPeak);
            currentProcessed.set_peak(primaryCurrentOffset + primaryCurrentPeakToPeak / 2);
            currentProcessed.set_duty_cycle(dutyCycle);
            currentProcessed.set_offset(primaryCurrentOffset);
            currentProcessed.set_dead_time(deadTime);

            voltageProcessed.set_peak_to_peak(primaryVoltavePeaktoPeak);
            voltageProcessed.set_peak(inputVoltage);
            voltageProcessed.set_duty_cycle(dutyCycle);
            voltageProcessed.set_offset(0);
            voltageProcessed.set_dead_time(deadTime);
            switch (mode) {
                case Flyback::Modes::ContinuousCurrentMode: {
                    voltageWaveform = InputsWrapper::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    voltageProcessed.set_label(WaveformLabel::RECTANGULAR);
                    break;
                }
                case Flyback::Modes::DiscontinuousCurrentMode: {
                    voltageWaveform = InputsWrapper::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    voltageProcessed.set_label(WaveformLabel::RECTANGULAR_WITH_DEADTIME);
                    break;
                }
            }

            excitation.set_frequency(switchingFrequency);
            SignalDescriptor current;
            current.set_waveform(currentWaveform);
            currentProcessed = InputsWrapper::calculate_processed_data(currentWaveform, switchingFrequency, true, currentProcessed);
            auto sampledCurrentWaveform = OpenMagnetics::InputsWrapper::calculate_sampled_waveform(currentWaveform, switchingFrequency);
            auto currentHarmonics = OpenMagnetics::InputsWrapper::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency);
            current.set_processed(currentProcessed);
            current.set_harmonics(currentHarmonics);
            excitation.set_current(current);
            SignalDescriptor voltage;
            voltage.set_waveform(voltageWaveform);
            voltageProcessed = InputsWrapper::calculate_processed_data(voltageWaveform, switchingFrequency, true, voltageProcessed);
            auto sampledVoltageWaveform = OpenMagnetics::InputsWrapper::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
            auto voltageHarmonics = OpenMagnetics::InputsWrapper::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency);
            voltage.set_processed(voltageProcessed);
            voltage.set_harmonics(voltageHarmonics);
            excitation.set_voltage(voltage);
            json isolationSideJson;
            to_json(isolationSideJson, get_isolation_side_from_index(0));
            excitation.set_name(isolationSideJson);
            excitation = InputsWrapper::prune_harmonics(excitation, Defaults().harmonicAmplitudeThreshold, 1);

            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // Secondaries
        for (size_t secondaryIndex = 0; secondaryIndex < turnsRatios.size(); ++secondaryIndex) {

            OperatingPointExcitation excitation;
            Waveform currentWaveform;
            Waveform voltageWaveform;
            Processed currentProcessed;
            Processed voltageProcessed;

            double secondaryPower = get_total_input_power(outputOperatingPoint.get_output_currents()[secondaryIndex], outputOperatingPoint.get_output_voltages()[secondaryIndex], 1, 0);
            double powerDivider = secondaryPower / totalOutputPower;

            double secondaryVoltagePeaktoPeak = inputVoltage / turnsRatios[secondaryIndex] + get_diode_voltage_drop() + outputOperatingPoint.get_output_voltages()[secondaryIndex];
            double secondaryCurrentAverage = centerPrimaryCurrentRamp * turnsRatios[secondaryIndex] * powerDivider;
            double secondaryCurrentPeaktoPeak = secondaryCurrentAverage * currentRippleRatio * 2;
            double secondaryCurrentOffset = std::max(0.0, secondaryCurrentAverage - secondaryCurrentPeaktoPeak / 2);


            currentProcessed.set_peak_to_peak(secondaryCurrentPeaktoPeak);
            currentProcessed.set_peak(secondaryCurrentOffset + secondaryCurrentPeaktoPeak / 2);
            currentProcessed.set_duty_cycle(dutyCycle);
            currentProcessed.set_offset(secondaryCurrentOffset);
            currentProcessed.set_dead_time(deadTime);

            voltageProcessed.set_peak_to_peak(secondaryVoltagePeaktoPeak);
            voltageProcessed.set_peak(outputOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
            voltageProcessed.set_duty_cycle(dutyCycle);
            voltageProcessed.set_offset(0);
            voltageProcessed.set_dead_time(deadTime);

            switch (mode) {
                case Flyback::Modes::ContinuousCurrentMode: {
                    voltageWaveform = InputsWrapper::create_waveform(WaveformLabel::RECTANGULAR, secondaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    currentWaveform = InputsWrapper::create_waveform(WaveformLabel::FLYBACK_SECONDARY, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    voltageProcessed.set_label(WaveformLabel::SECONDARY_RECTANGULAR);
                    currentProcessed.set_label(WaveformLabel::FLYBACK_SECONDARY);
                    break;
                }
                case Flyback::Modes::DiscontinuousCurrentMode: {
                    voltageWaveform = InputsWrapper::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, secondaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    currentWaveform = InputsWrapper::create_waveform(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    voltageProcessed.set_label(WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
                    currentProcessed.set_label(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
                    break;
                }
            }

            excitation.set_frequency(switchingFrequency);
            SignalDescriptor current;
            current.set_waveform(currentWaveform);
            currentProcessed = InputsWrapper::calculate_processed_data(currentWaveform, switchingFrequency, true, currentProcessed);
            auto sampledCurrentWaveform = OpenMagnetics::InputsWrapper::calculate_sampled_waveform(currentWaveform, switchingFrequency);
            auto currentHarmonics = OpenMagnetics::InputsWrapper::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency);
            current.set_processed(currentProcessed);
            current.set_harmonics(currentHarmonics);
            excitation.set_current(current);
            SignalDescriptor voltage;
            voltage.set_waveform(voltageWaveform);
            voltageProcessed = InputsWrapper::calculate_processed_data(voltageWaveform, switchingFrequency, true, voltageProcessed);
            auto sampledVoltageWaveform = OpenMagnetics::InputsWrapper::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
            auto voltageHarmonics = OpenMagnetics::InputsWrapper::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency);
            voltage.set_processed(voltageProcessed);
            voltage.set_harmonics(voltageHarmonics);
            excitation.set_voltage(voltage);
            json isolationSideJson;
            to_json(isolationSideJson, get_isolation_side_from_index(secondaryIndex + 1));
            excitation.set_name(isolationSideJson);
            excitation = InputsWrapper::prune_harmonics(excitation, Defaults().harmonicAmplitudeThreshold, 1);
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    double Flyback::calculate_maximum_duty_cycle(double minimumInputVoltage, double outputReflectedVoltage, Flyback::Modes mode) {
        switch (mode) {
            case Flyback::Modes::ContinuousCurrentMode:
                return outputReflectedVoltage / (outputReflectedVoltage + minimumInputVoltage);
            case Flyback::Modes::DiscontinuousCurrentMode:
                return outputReflectedVoltage / (outputReflectedVoltage + minimumInputVoltage);
        }

        throw std::runtime_error("Unknown flyback mode");
    }

    double Flyback::get_total_input_power(std::vector<double> outputCurrents, std::vector<double> outputVoltages, double efficiency, double diodeVoltageDrop) {
        double totalPower = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputCurrents.size(); ++secondaryIndex) {
            totalPower += outputCurrents[secondaryIndex] * (outputVoltages[secondaryIndex] + diodeVoltageDrop);
        }

        return totalPower / efficiency;
    }


    double Flyback::get_total_input_power(double outputCurrent, double outputVoltage, double efficiency, double diodeVoltageDrop) {
        double totalPower = outputCurrent * (outputVoltage + diodeVoltageDrop);

        return totalPower / efficiency;
    }

    double Flyback::get_needed_inductance(double inputVoltage, double inputPower, double dutyCycle, double frequency, double currentRippleRatio) {
        double averageInputCurrent = inputPower / (inputVoltage * dutyCycle);
        double inputCurrentRipple = currentRippleRatio * averageInputCurrent;
        return inputVoltage * dutyCycle / (frequency * inputCurrentRipple);
    }

    double Flyback::get_minimum_output_reflected_voltage(double maximumDrainSourceVoltage, double maximumInputVoltage, double safetyMargin) {
        return maximumDrainSourceVoltage * safetyMargin - maximumInputVoltage;
    }

    bool Flyback::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw std::runtime_error("At least one operating point is needed");
        }
        for (size_t flybackOperatingPointIndex = 1; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            if (get_operating_points()[flybackOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw std::runtime_error("Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[flybackOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
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

    InputsWrapper Flyback::process() {
        Flyback::run_checks(_assertErrors);

        InputsWrapper inputs;

        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        auto mode = get_current_ripple_ratio() < 1? Flyback::Modes::ContinuousCurrentMode : Flyback::Modes::DiscontinuousCurrentMode;

        if (!get_maximum_drain_source_voltage() && !get_maximum_duty_cycle()) {
            throw std::invalid_argument("Missing both maximum duty cycle and maximum drain source voltage");
        }
        double maximumDutyCycle;
        double maximumNeededInductance = 0;
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size(), 0);

        if (get_maximum_duty_cycle()) {
            maximumDutyCycle = get_maximum_duty_cycle().value();
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
        else {
            double maximumDrainSourceVoltage = get_maximum_drain_source_voltage().value();
            auto minimumOutputReflectedVoltage = get_minimum_output_reflected_voltage(maximumDrainSourceVoltage, maximumInputVoltage);
            maximumDutyCycle = calculate_maximum_duty_cycle(minimumInputVoltage, minimumOutputReflectedVoltage, mode);
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
                for (size_t secondaryIndex = 0; secondaryIndex < flybackOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                    auto turnsRatio = minimumOutputReflectedVoltage / (flybackOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
                    turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatio);
                }
            }
        }

        for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
            double switchingFrequency = flybackOperatingPoint.get_switching_frequency();
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
            maximumNeededInductance = std::max(maximumNeededInductance, neededInductance);
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
        inductanceWithTolerance.set_nominal(roundFloat(maximumNeededInductance, 10));
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
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, maximumNeededInductance, mode);
                json mierda;
                to_json(mierda, operatingPoint);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(flybackOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        for (auto operatingPoint : inputs.get_mutable_operating_points()) {
            auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
            auto primarycurrentProcessed = primaryExcitation.get_current()->get_processed().value();
        }

        return inputs;
    }

    InputsWrapper AdvancedFlyback::process() {
        Flyback::run_checks(_assertErrors);

        InputsWrapper inputs;

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
        designRequirements.set_topology(Topologies::FLYBACK_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                std::optional<double> customDeadTime = std::nullopt;
                if (get_desired_duty_cycle().size() <= flybackOperatingPointIndex) {
                    throw std::runtime_error("Missing duty cycle for flybackOperatingPointIndex: " + std::to_string(flybackOperatingPointIndex));
                }
                double customDutyCycle = get_desired_duty_cycle()[flybackOperatingPointIndex][inputVoltageIndex];

                if (get_desired_dead_time()) {
                    if (get_desired_dead_time()->size() <= flybackOperatingPointIndex) {
                        throw std::runtime_error("Missing dead time for flybackOperatingPointIndex: " + std::to_string(flybackOperatingPointIndex));
                    }
                    customDeadTime = get_desired_dead_time().value()[flybackOperatingPointIndex];
                }

                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, maximumNeededInductance, std::nullopt, customDutyCycle, customDeadTime);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(flybackOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        for (auto operatingPoint : inputs.get_mutable_operating_points()) {
            auto primaryExcitation = operatingPoint.get_excitations_per_winding()[0];
            auto primarycurrentProcessed = primaryExcitation.get_current()->get_processed().value();
        }

        return inputs;
    }

} // namespace OpenMagnetics
