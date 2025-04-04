#include "Topology.h"
#include "Utils.h"
#include <cfloat>

namespace OpenMagnetics {

    Flyback::Flyback(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint Flyback::processOperatingPointsForInputVoltage(double inputVoltage, Flyback::FlybackOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, Flyback::Modes mode) {
        double totalInputPower = get_total_input_power(outputOperatingPoint.get_output_currents(), outputOperatingPoint.get_output_voltages(), get_efficiency());

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double inputCurrent = totalInputPower / inputVoltage;

        double deadTime = 0;
        double maximumReflectedOutputVoltage = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
            double outputVoltage = outputOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop();
            maximumReflectedOutputVoltage = std::max(maximumReflectedOutputVoltage, outputVoltage * turnsRatios[secondaryIndex]);
        }

        double primaryVoltavePeaktoPeak = inputVoltage + maximumReflectedOutputVoltage;
        double dutyCycle = sqrt(get_current_ripple_ratio() * 2 * totalInputPower / pow(inputVoltage, 2) * inductance * switchingFrequency);

        if (dutyCycle > 1 ) {
            throw std::runtime_error("dutyCycle cannot be larger than one: " + std::to_string(dutyCycle));
        }

        double primaryCurrentAverage = totalInputPower / (inputVoltage * dutyCycle);
        double primaryCurrentPeakToPeak = inputVoltage * dutyCycle / (inductance * switchingFrequency);
        double primaryCurrentOffset = primaryCurrentAverage - primaryCurrentPeakToPeak / 2;
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
                    voltageWaveform = InputsWrapper::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle);
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

            double secondaryPower = get_total_input_power(outputOperatingPoint.get_output_currents()[secondaryIndex], outputOperatingPoint.get_output_voltages()[secondaryIndex], get_efficiency());
            double powerDivider = secondaryPower / totalInputPower;

            double secondaryVoltagePeaktoPeak = inputVoltage / turnsRatios[secondaryIndex] + get_diode_voltage_drop() + outputOperatingPoint.get_output_voltages()[secondaryIndex];
            double secondaryCurrentPeaktoPeak = primaryCurrentPeakToPeak * turnsRatios[secondaryIndex] * powerDivider;
            double secondaryCurrentOffset = primaryCurrentOffset * turnsRatios[secondaryIndex] * powerDivider;

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
                    voltageWaveform = InputsWrapper::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, secondaryVoltagePeaktoPeak, switchingFrequency, dutyCycle);
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

    double Flyback::get_maximum_duty_cycle(double minimumInputVoltage, double outputReflectedVoltage, Flyback::Modes mode) {
        switch (mode) {
            case Flyback::Modes::ContinuousCurrentMode:
                return outputReflectedVoltage / (outputReflectedVoltage + minimumInputVoltage);
            case Flyback::Modes::DiscontinuousCurrentMode:
                return outputReflectedVoltage / (outputReflectedVoltage + minimumInputVoltage);
        }

        throw std::runtime_error("Unknown flyback mode");
    }

    double Flyback::get_total_input_power(std::vector<double> outputCurrents, std::vector<double> outputVoltages, double efficiency) {
        double totalPower = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputCurrents.size(); ++secondaryIndex) {
            totalPower += outputCurrents[secondaryIndex] * (outputVoltages[secondaryIndex] + get_diode_voltage_drop());
        }

        return totalPower / efficiency;
    }


    double Flyback::get_total_input_power(double outputCurrent, double outputVoltage, double efficiency) {
        double totalPower = outputCurrent * (outputVoltage + get_diode_voltage_drop());

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
        auto checksPassed = Flyback::run_checks(_assertErrors);

        InputsWrapper inputs;

        auto maximumDrainSourceVoltage = get_maximum_drain_source_voltage();
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        auto mode = get_current_ripple_ratio() < 1? Flyback::Modes::ContinuousCurrentMode : Flyback::Modes::DiscontinuousCurrentMode;

        auto minimumOutputReflectedVoltage = get_minimum_output_reflected_voltage(maximumDrainSourceVoltage, maximumInputVoltage);
        double maximumDutyCycle = get_maximum_duty_cycle(minimumInputVoltage, minimumOutputReflectedVoltage, mode);
        double maximumNeededInductance = 0;
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size(), 0);
        for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
            double switchingFrequency = flybackOperatingPoint.get_switching_frequency();
            double totalInputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), get_efficiency());
            double neededInductance = get_needed_inductance(minimumInputVoltage, totalInputPower, maximumDutyCycle, switchingFrequency, get_current_ripple_ratio());
            maximumNeededInductance = std::max(maximumNeededInductance, neededInductance);
            for (size_t secondaryIndex = 0; secondaryIndex < flybackOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                auto turnsRatio = minimumOutputReflectedVoltage / (flybackOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
                turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatio);
            }
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

} // namespace OpenMagnetics
