#include "converter_models/Topology.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/Utils.h"
#include <cfloat>

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
                throw std::runtime_error("Either current ripple ratio or mode is needed for the Flyback OperatingPoint Mode");
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
                throw std::runtime_error("Either switching frequency or mode is needed for the Flyback OperatingPoint");
            }
            auto mode = get_mode().value();
            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    throw std::runtime_error("Switching Frequency is needed for CCM");
                }
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    throw std::runtime_error("Switching Frequency is needed for DCM");
                }
                case FlybackModes::QUASI_RESONANT_MODE: {
                    if (!inductance) {
                        throw std::runtime_error("Inductance in missing for switching frequency calculation");
                    }
                    if (!turnsRatios) {
                        throw std::runtime_error("TurnsRatios in missing for switching frequency calculation");
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
                        throw std::runtime_error("Inductance in missing for switching frequency calculation");
                    }
                    if (!turnsRatios) {
                        throw std::runtime_error("TurnsRatios in missing for switching frequency calculation");
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
                  throw std::runtime_error("Unknown mode in Flyback");
            }
        }
    }

    OperatingPoint Flyback::processOperatingPointsForInputVoltage(double inputVoltage, FlybackOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, std::optional<FlybackModes> customMode, std::optional<double> customDutyCycle, std::optional<double> customDeadTime) {

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
            OperatingPointExcitation excitation;
            Waveform currentWaveform;
            Waveform voltageWaveform;
            Processed currentProcessed;
            Processed voltageProcessed;

            currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, primaryCurrentOffset, deadTime);
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
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    voltageProcessed.set_label(WaveformLabel::RECTANGULAR);
                    break;
                }
                case FlybackModes::QUASI_RESONANT_MODE:
                case FlybackModes::BOUNDARY_MODE_OPERATION:
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    voltageProcessed.set_label(WaveformLabel::RECTANGULAR_WITH_DEADTIME);
                    break;
                }
            }

            excitation.set_frequency(switchingFrequency);
            SignalDescriptor current;
            current.set_waveform(currentWaveform);
            currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true, currentProcessed);
            auto sampledCurrentWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
            auto currentHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency);
            current.set_processed(currentProcessed);
            current.set_harmonics(currentHarmonics);
            excitation.set_current(current);
            SignalDescriptor voltage;
            voltage.set_waveform(voltageWaveform);
            voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true, voltageProcessed);
            auto sampledVoltageWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
            auto voltageHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency);
            voltage.set_processed(voltageProcessed);
            voltage.set_harmonics(voltageHarmonics);
            excitation.set_voltage(voltage);
            json isolationSideJson;
            to_json(isolationSideJson, get_isolation_side_from_index(0));
            excitation.set_name(isolationSideJson);
            excitation = Inputs::prune_harmonics(excitation, Defaults().harmonicAmplitudeThreshold, 1);

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
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, secondaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    voltageProcessed.set_label(WaveformLabel::SECONDARY_RECTANGULAR);
                    currentProcessed.set_label(WaveformLabel::FLYBACK_SECONDARY);
                    break;
                }
                case FlybackModes::QUASI_RESONANT_MODE:
                case FlybackModes::BOUNDARY_MODE_OPERATION:
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, secondaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    voltageProcessed.set_label(WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
                    currentProcessed.set_label(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
                    break;
                }
            }

            excitation.set_frequency(switchingFrequency);
            SignalDescriptor current;
            current.set_waveform(currentWaveform);
            currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true, currentProcessed);
            auto sampledCurrentWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
            auto currentHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency);
            current.set_processed(currentProcessed);
            current.set_harmonics(currentHarmonics);
            excitation.set_current(current);
            SignalDescriptor voltage;
            voltage.set_waveform(voltageWaveform);
            voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true, voltageProcessed);
            auto sampledVoltageWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
            auto voltageHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency);
            voltage.set_processed(voltageProcessed);
            voltage.set_harmonics(voltageHarmonics);
            excitation.set_voltage(voltage);
            json isolationSideJson;
            to_json(isolationSideJson, get_isolation_side_from_index(secondaryIndex + 1));
            excitation.set_name(isolationSideJson);
            excitation = Inputs::prune_harmonics(excitation, Defaults().harmonicAmplitudeThreshold, 1);
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
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

    DesignRequirements Flyback::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_maximum_drain_source_voltage() && !get_maximum_duty_cycle()) {
            throw std::invalid_argument("Missing both maximum duty cycle and maximum drain source voltage");
        }
        double maximumNeededInductance = 0;
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
            maximumNeededInductance = std::max(maximumNeededInductance, neededInductance);
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
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, magnetizingInductance, mode);

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

    Inputs Flyback::process() {
        Flyback::run_checks(_assertErrors);

        Inputs inputs;
        auto designRequirements = process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto turnsRatio : designRequirements.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(turnsRatio));
        }
        auto desiredMagnetizingInductance = resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        auto operatingPoints = process_operating_points(turnsRatios, desiredMagnetizingInductance);

        inputs.set_design_requirements(designRequirements);
        inputs.set_operating_points(operatingPoints);

        return inputs;
    }

    std::vector<OperatingPoint> Flyback::process_operating_points(Magnetic magnetic) {
        Flyback::run_checks(_assertErrors);

        std::vector<OperatingPoint> operatingPoints;

        if (!get_maximum_drain_source_voltage() && !get_maximum_duty_cycle()) {
            throw std::invalid_argument("Missing both maximum duty cycle and maximum drain source voltage");
        }
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();

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
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, magnetizingInductance, mode);

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

    Inputs AdvancedFlyback::process() {
        Flyback::run_checks(_assertErrors);

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

        return inputs;
    }

/** ------------------
 *  Inverter topology
 *  ------------------
 * */

    MyInverter::MyInverter(const json& j) {
        from_json(j, *this);
    }

    bool MyInverter::run_checks(bool assert) {
        (void) assert;
        return true;
    }

    DesignRequirements MyInverter::process_design_requirements() {
        DesignRequirements designRequirements;
        designRequirements.set_magnetizing_inductance();
        designRequirements.set_operating_temperature();
        designRequirements.set_application(1);      //POWER
        designRequirements.set_sub_application(3); //POWER_FILTERING
        return designRequirements;
    }


    std::complex<double> compute_load_impedance(const InverterLoad& load, double omega) {
        using namespace std::complex_literals;

        auto type = load.get_load_type();

        switch (type) {
            case LoadType::GRID: {
                // Grid modeled as R + jωL
                double R = load.get_grid_resistance().value_or(0.0);
                double L = load.get_grid_inductance().value_or(0.0);
                return std::complex<double>(R, omega * L);
            }

            case LoadType::R_L: {
                // General load R + jωL
                double R = load.get_resistance().value_or(0.0);
                double L = load.get_inductance().value_or(0.0);
                return std::complex<double>(R, omega * L);
            }

            default:
                throw std::runtime_error("Unknown load type");
        }
    }


    std::complex<double> compute_filter_impedance(const InverterDownstreamFilter& filter, double omega) {
        using namespace std::complex_literals;

        auto topology = filter.get_filter_topology();

        // === Build element impedances with ESR ===

        // Inductor 1
        double L1 = filter.get_inductor().get_inductance();
        double ESR_L1 = filter.get_inductor().get_esr().value_or(0.0);
        std::complex<double> ZL1(ESR_L1, omega * L1);

        // Capacitor (if present)
        std::complex<double> ZC;
        if (filter.get_capacitor()) {
            double C = filter.get_capacitor()->get_capacitance();
            double ESR_C = filter.get_capacitor()->get_esr().value_or(0.0);
            // Capacitor ESR is in series with reactive part
            ZC = std::complex<double>(ESR_C, -1.0 / (omega * C));
        }

        // Inductor 2 (if present)
        std::complex<double> ZL2;
        if (filter.get_inductor2()) {
            double L2 = filter.get_inductor2()->get_inductance();
            double ESR_L2 = filter.get_inductor2()->get_esr().value_or(0.0);
            ZL2 = std::complex<double>(ESR_L2, omega * L2);
        }

        // === Switch based on topology ===
        switch (topology) {
            case FilterTopologies::L:
                // Just series L1
                return ZL1;

            case FilterTopologies::LC:
                if (!filter.get_capacitor()) {
                    throw std::runtime_error("LC topology requires a capacitor");
                }
                // L1 in series with capacitor to ground (parallel branch)
                return (ZL1 * ZC) / (ZL1 + ZC);

            case FilterTopologies::LCL:
                if (!filter.get_capacitor() || !filter.get_inductor2()) {
                    throw std::runtime_error("LCL topology requires capacitor and second inductor");
                }
                // L1 in series with (ZC || ZL2)
                return ZL1 + (ZC * ZL2) / (ZC + ZL2);

            default:
                throw std::runtime_error("Unknown filter topology");
        }
    }

    /// --- FFT (simple DFT for clarity, can be optimized) ---
    std::vector<std::complex<double>> compute_fft(const std::vector<double>& signal) {
        int N = signal.size();
        std::vector<std::complex<double>> spectrum(N);

        for (int k = 0; k < N; ++k) {
            std::complex<double> sum = 0.0;
            for (int n = 0; n < N; ++n) {
                double angle = -2.0 * M_PI * k * n / N;
                sum += signal[n] * std::exp(std::complex<double>(0, angle));
            }
            spectrum[k] = sum / (double)N;
        }
        return spectrum;
    }

    /// Helper: dq -> abc transformation
    ABCVoltages dq_to_abc(const std::complex<double>& Vdq, double theta) {
        // Vdq = Vd + jVq
        double Vd = Vdq.real();
        double Vq = Vdq.imag();

        // Inverse Park transform
        double Va = Vd * cos(theta) - Vq * sin(theta);
        double Vb = Vd * cos(theta - 2.0*M_PI/3.0) - Vq * sin(theta - 2.0*M_PI/3.0);
        double Vc = Vd * cos(theta + 2.0*M_PI/3.0) - Vq * sin(theta + 2.0*M_PI/3.0);

        return {Va, Vb, Vc};
    }

        /// Clarke transform: abc -> alpha-beta
    std::pair<double,double> abc_to_alphabeta(const ABCVoltages& v) {
        double v_alpha = (2.0/3.0) * (v.Va - 0.5*v.Vb - 0.5*v.Vc);
        double v_beta  = (2.0/3.0) * ((sqrt(3)/2.0) * (v.Vb - v.Vc));
        return {v_alpha, v_beta};
    }

    /// SVPWM modulation: compute duty cycles for abc legs
    ABCVoltages svpwm_modulation(const ABCVoltages& Vabc,
                                double ma,
                                double Vdc,
                                double fsw) {
        auto [v_alpha, v_beta] = abc_to_alphabeta(Vabc);

        // Scale with modulation index
        v_alpha *= ma;
        v_beta  *= ma;

        double Ts = 1.0 / fsw;

        // Reference angle
        double theta = atan2(v_beta, v_alpha);
        if (theta < 0) theta += 2*M_PI;

        // Sector detection (1–6, each 60°)
        int sector = int(theta / (M_PI/3.0)) + 1;
        if (sector > 6) sector = 6;

        // Magnitude
        double Vref = sqrt(v_alpha*v_alpha + v_beta*v_beta);

        // Compute T1, T2 using standard SVPWM geometry
        double T1, T2;
        double angle = theta - (sector-1)*(M_PI/3.0);

        T1 = (Ts * sqrt(3) * Vref / Vdc) * sin(M_PI/3.0 - angle);
        T2 = (Ts * sqrt(3) * Vref / Vdc) * sin(angle);
        double T0 = Ts - T1 - T2;

        if (T0 < 0) T0 = 0; // clamp

        // Duty cycles depend on sector
        double Ta, Tb, Tc;
        switch (sector) {
            case 1: Ta = (T1+T2+T0/2)/Ts; Tb = (T2+T0/2)/Ts; Tc = T0/2/Ts; break;
            case 2: Ta = (T1+T0/2)/Ts; Tb = (T1+T2+T0/2)/Ts; Tc = T0/2/Ts; break;
            case 3: Ta = T0/2/Ts; Tb = (T1+T2+T0/2)/Ts; Tc = (T2+T0/2)/Ts; break;
            case 4: Ta = T0/2/Ts; Tb = (T1+T0/2)/Ts; Tc = (T1+T2+T0/2)/Ts; break;
            case 5: Ta = (T2+T0/2)/Ts; Tb = T0/2/Ts; Tc = (T1+T2+T0/2)/Ts; break;
            case 6: Ta = (T1+T2+T0/2)/Ts; Tb = T0/2/Ts; Tc = (T1+T0/2)/Ts; break;
            default: Ta=Tb=Tc=0.5; break;
        }

        // Return equivalent duty values as "abc voltages"
        return {Ta, Tb, Tc};
    }

    ABCVoltages compute_voltage_references(const Inverter& inverter,
                                        const InverterOperatingPoint& op_point,
                                        const Modulation& modulation,
                                        double grid_angle_rad) {
        using namespace std::complex_literals;

        const InverterLoad& load = op_point.get_load();
        double omega = 2.0 * M_PI * op_point.get_fundamental_frequency();

        std::complex<double> Vref_dq;

        // --- Step 1: dq reference (same as before) ---
        switch (load.get_load_type()) {
            case LoadType::GRID: {
                double Vg_rms = load.get_phase_voltage().value_or(230.0);
                std::complex<double> Vg = std::complex<double>(Vg_rms, 0.0);
                std::complex<double> Zg = std::complex<double>(
                    load.get_grid_resistance().value_or(0.0),
                    omega * load.get_grid_inductance().value_or(0.0)
                );
                double P = op_point.get_output_power().value_or(0.0);
                double pf = op_point.get_power_factor().value_or(1.0);
                double phi = acos(pf);
                std::complex<double> Iph = (P / (Vg_rms * pf)) * (cos(-phi) + 1i*sin(-phi));
                Vref_dq = Vg - Iph * Zg;
                break;
            }
            case LoadType::R_L: {
                std::complex<double> Zload = std::complex<double>(
                    load.get_resistance().value_or(0.0),
                    omega * load.get_inductance().value_or(0.0)
                );
                double P = op_point.get_output_power().value_or(0.0);
                double pf = op_point.get_power_factor().value_or(1.0);
                double phi = acos(pf);
                std::complex<double> Iph = (P / pf) / inverter.get_line_rms_current().get_nominal()
                                        * (cos(-phi) + 1i*sin(-phi));
                Vref_dq = Iph * Zload;
                break;
            }
            default:
                throw std::runtime_error("Unknown load type");
        }

        // --- Step 2: dq -> abc ---
        ABCVoltages Vabc = dq_to_abc(Vref_dq, grid_angle_rad);

        // --- Step 3: Modulation strategy ---
        double ma = modulation.get_modulation_depth();
        switch (modulation.get_modulation_strategy()) {
            case ModulationStrategy::SPWM: {
                Vabc.Va *= ma;
                Vabc.Vb *= ma;
                Vabc.Vc *= ma;
                break;
            }
            case ModulationStrategy::THIPWM: {
                double k = modulation.get_third_harmonic_injection_coefficient().value_or(1.0/6.0);
                Vabc.Va = ma * (Vabc.Va + k * sin(3.0 * grid_angle_rad));
                Vabc.Vb = ma * (Vabc.Vb + k * sin(3.0 * grid_angle_rad - 2.0*M_PI/3.0));
                Vabc.Vc = ma * (Vabc.Vc + k * sin(3.0 * grid_angle_rad + 2.0*M_PI/3.0));
                break;
            }
            case ModulationStrategy::SVPWM: {
                Vabc = svpwm_modulation(Vabc, ma,
                                        inverter.get_dc_bus_voltage().get_nominal(),
                                        modulation.get_switching_frequency());
                break;
            }
            default:
                throw std::runtime_error("Unknown modulation strategy");
        }

        return Vabc;
    }

    double compute_carrier(const Modulation& modulation, double t) {
        double fsw = modulation.get_switching_frequency();
        double Ts = 1.0 / fsw;

        // Normalize time within switching period
        double tau = fmod(t, Ts) / Ts;  // [0,1)

        switch (modulation.get_pwm_type()) {
            case PwmType::SAWTOOTH: {
                // Linear ramp: -1 → +1
                return 2.0 * tau - 1.0;
            }

            case PwmType::TRIANGULAR: {
                // Symmetric triangular: -1 → +1 → -1
                if (tau < 0.5) {
                    return 4.0 * tau - 1.0;    // ramp -1 → +1
                } else {
                    return -4.0 * tau + 3.0;   // ramp +1 → -1
                }
            }

            default:
                throw std::runtime_error("Unknown PWM carrier type");
        }
    }

    PwmSignals compare_with_carrier(const ABCVoltages& Vabc,
                                    double carrier,
                                    double Vdc,
                                    const Modulation& modulation) {
        // Normalize into duty cycles [0,1]
        double d_a = 0.5 * (Vabc.Va / (Vdc / 2.0) + 1.0);
        double d_b = 0.5 * (Vabc.Vb / (Vdc / 2.0) + 1.0);
        double d_c = 0.5 * (Vabc.Vc / (Vdc / 2.0) + 1.0);

        // Clamp to [0,1]
        d_a = std::clamp(d_a, 0.0, 1.0);
        d_b = std::clamp(d_b, 0.0, 1.0);
        d_c = std::clamp(d_c, 0.0, 1.0);

        // Apply deadtime correction (if set)
        if (modulation.get_deadtime()) {
            double t_dead = *modulation.get_deadtime();
            double Ts = 1.0 / modulation.get_switching_frequency();
            double delta_d = t_dead / Ts;  // fraction of duty cycle lost to deadtime

            d_a = std::clamp(d_a - delta_d, 0.0, 1.0);
            d_b = std::clamp(d_b - delta_d, 0.0, 1.0);
            d_c = std::clamp(d_c - delta_d, 0.0, 1.0);
        }

        // Apply rise time correction (if set)
        if (modulation.get_rise_time()) {
            double t_rise = *modulation.get_rise_time();
            double Ts = 1.0 / modulation.get_switching_frequency();
            double delta_d = t_rise / Ts;

            d_a = std::clamp(d_a - delta_d/2, 0.0, 1.0);
            d_b = std::clamp(d_b - delta_d/2, 0.0, 1.0);
            d_c = std::clamp(d_c - delta_d/2, 0.0, 1.0);
        }

        // Carrier comparison (carrier is [-1,1], convert duty back to same scale)
        double Va_norm = 2.0*d_a - 1.0;
        double Vb_norm = 2.0*d_b - 1.0;
        double Vc_norm = 2.0*d_c - 1.0;

        bool Sa = (Va_norm > carrier);
        bool Sb = (Vb_norm > carrier);
        bool Sc = (Vc_norm > carrier);

        return {Sa, Sb, Sc};
    }

    NodeResult solve_filter_topology(const InverterDownstreamFilter& filter,
                                    const InverterLoad& load,
                                    double omega,
                                    std::complex<double> Vinv)
    {
        using namespace std::complex_literals;

        // --- Build element impedances ---
        std::complex<double> ZL1(filter.get_inductor().get_esr().value_or(0.0),
                                omega * filter.get_inductor().get_inductance());

        std::complex<double> ZC(1e9, 0.0); // open if no C
        if (filter.get_capacitor()) {
            double C = filter.get_capacitor()->get_capacitance();
            double ESRc = filter.get_capacitor()->get_esr().value_or(0.0);
            ZC = std::complex<double>(ESRc, -1.0/(omega*C));
        }

        std::complex<double> ZL2(1e9, 0.0); // open if no L2
        if (filter.get_inductor2()) {
            double L2 = filter.get_inductor2()->get_inductance();
            double ESR2 = filter.get_inductor2()->get_esr().value_or(0.0);
            ZL2 = std::complex<double>(ESR2, omega * L2);
        }

        std::complex<double> Zload = compute_load_impedance(load, omega);

        // --- Solve depending on topology ---
        std::complex<double> vNode;
        switch (filter.get_filter_topology()) {
            case FilterTopologies::L: {
                vNode = Vinv * (Zload / (ZL1 + Zload));
                break;
            }
            case FilterTopologies::LC: {
                std::complex<double> Zeq = (ZC * Zload) / (ZC + Zload);
                vNode = Vinv * (Zeq / (ZL1 + Zeq));
                break;
            }
            case FilterTopologies::LCL: {
                std::complex<double> Zeq = (ZC * (ZL2 + Zload)) / (ZC + ZL2 + Zload);
                vNode = Vinv * (Zeq / (ZL1 + Zeq));
                break;
            }
            default:
                throw std::runtime_error("Unknown filter topology");
        }

        std::complex<double> vL1 = Vinv - vNode;
        std::complex<double> iL1 = vL1 / ZL1;

        return {vNode, vL1, iL1};
    }

    HarmonicsBundle compute_harmonics(const Modulation& modulation,
                                    const ABCVoltages& Vabc,
                                    double Vdc,
                                    std::complex<double> Vfund,
                                    std::complex<double> Ifund,
                                    double f1,
                                    const InverterDownstreamFilter& filter,
                                    const InverterLoad& load,
                                    int Nperiods = 5,
                                    int samplesPerPeriod = 200) {
        double fsw = modulation.get_switching_frequency();
        double Ts = 1.0 / fsw;
        double fs = fsw * samplesPerPeriod;

        // --- 1. Generate switching waveform (leg A as example) ---
        std::vector<double> waveform;
        int Nsamples = Nperiods * samplesPerPeriod;
        for (int n = 0; n < Nsamples; ++n) {
            double t = n * Ts / samplesPerPeriod;
            double carrier = compute_carrier(modulation, t);
            PwmSignals gates = compare_with_carrier(Vabc, carrier, Vdc, modulation);

            double vA = gates.Sa ? +Vdc/2.0 : -Vdc/2.0;
            waveform.push_back(vA);
        }

        // --- 2. FFT ---
        auto spectrum = compute_fft(waveform);
        int N = spectrum.size();

        // --- 3. Build harmonics up to 5*fsw ---
        HarmonicsBundle bundle;
        double fmax = 5.0 * fsw;

        for (int k = 0; k < N/2; ++k) {
            double f = k * fs / N;
            if (f <= fmax) {
                double omega_k = 2.0 * M_PI * f;

                std::complex<double> Vph = spectrum[k];

                // Use the *same filter equations* as fundamental, but per harmonic
                NodeResult node = solve_filter_topology(filter, load, omega_k, Vph);

                std::complex<double> Iph = node.iL1;

                bundle.Vharm.push_back({f, Vph});
                bundle.Iharm.push_back({f, Iph});
            }
        }

        // --- 4. Override fundamental bins with clean phasors ---
        auto vIt = std::min_element(bundle.Vharm.begin(), bundle.Vharm.end(),
            [f1](const HarmonicComponent& a, const HarmonicComponent& b) {
                return std::abs(a.frequency - f1) < std::abs(b.frequency - f1);
            });
        if (vIt != bundle.Vharm.end()) {
            vIt->phasor = Vfund;
            vIt->frequency = f1;
        }

        auto iIt = std::min_element(bundle.Iharm.begin(), bundle.Iharm.end(),
            [f1](const HarmonicComponent& a, const HarmonicComponent& b) {
                return std::abs(a.frequency - f1) < std::abs(b.frequency - f1);
            });
        if (iIt != bundle.Iharm.end()) {
            iIt->phasor = Ifund;
            iIt->frequency = f1;
        }

        return bundle;
    }

    std::vector<OperatingPoint> MyInverter::process_operating_points() {
        std::vector<OperatingPoint> operatingPointsResult;

        for (auto& op_point : this->operatingPoints) {
            double f1 = op_point.get_fundamental_frequency();
            double omega = 2.0 * M_PI * f1;

            // --- Step 1: Compute dq reference → abc
            ABCVoltages Vabc = compute_voltage_references(
                *this,
                op_point,
                op_point.get_modulation(),
                op_point.get_grid_angle()
            );

            // --- Step 2: Fundamental phasors (clean values)
            std::complex<double> Vfund = std::polar(
                op_point.get_phase_voltage().value_or(230.0),
                op_point.get_grid_angle()
            );

            // Solve full filter/load topology at fundamental
            NodeResult node = solve_filter_topology(
                op_point.get_filter(),
                op_point.get_load(),
                omega,
                Vfund
            );

            std::complex<double> Ifund = node.iL1;

            // --- Step 3: Harmonic analysis (both V and I at once)
            HarmonicsBundle bundle = compute_harmonics(
                op_point.get_modulation(),
                Vabc,
                this->get_dc_bus_voltage().get_nominal(),
                Vfund,
                Ifund,
                f1,
                op_point.get_filter(),
                op_point.get_load()
            );

            // --- Step 4: Signal descriptors
            SignalDescriptor voltageSig;
            voltageSig.set_fundamental_frequency(f1);
            voltageSig.set_harmonics(bundle.Vharm);
            voltageSig.set_rms(compute_rms_from_harmonics(bundle.Vharm));

            SignalDescriptor currentSig;
            currentSig.set_fundamental_frequency(f1);
            currentSig.set_harmonics(bundle.Iharm);
            currentSig.set_rms(compute_rms_from_harmonics(bundle.Iharm));

            // --- Step 5: Assemble excitation
            OperatingPointExcitation excitation;
            excitation.set_voltage(voltageSig);
            excitation.set_current(currentSig);

            // --- Step 6: Build operating point result
            OperatingPoint result;
            result.set_excitation(excitation);

            operatingPointsResult.push_back(result);
        }

        return operatingPointsResult;
    }

    Inputs MyInverter::process() {
        if (!run_checks(_assertErrors)) {
            throw std::runtime_error("Configuration checks failed");
        }
        Inputs inputs;
        auto designRequirements = process_design_requirements();
        auto operatingPoints = process_operating_points();
        inputs.set_design_requirements(designRequirements);
        inputs.set_operating_points(operatingPoints);
        return inputs;
    }

} // namespace OpenMagnetics
