#include "converter_models/Topology.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/Utils.h"
#include <cfloat>

#ifdef DEBUG_PLOTS
#include "DebugPlots.hpp"
namespace DebugPlots = debug_plots; // alias
#endif

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

    using ABCVoltages     = MyInverter::ABCVoltages;
    using PwmSignals      = MyInverter::PwmSignals;
    using NodeResult      = MyInverter::NodeResult;
    using HarmonicsBundle = MyInverter::HarmonicsBundle;

    MyInverter::MyInverter(const json& j) {
        from_json(j, *this);
    }

    bool MyInverter::run_checks(bool assert) {
        (void) assert;
        return true;
    }

    DesignRequirements MyInverter::process_design_requirements() {
        DesignRequirements designRequirements;
        DimensionWithTolerance inductanceWithTolerance;
        DimensionWithTolerance operatingTemp;

        for (size_t inverterOperatingPointsIndex = 0; inverterOperatingPointsIndex < get_operating_points().size(); ++inverterOperatingPointsIndex) {
            auto inverterOperatingPoint = get_operating_points()[inverterOperatingPointsIndex];
            operatingTemp.set_nominal(inverterOperatingPoint.get_operating_temperature());
        }

        if (auto filter = this->get_downstream_filter()) {
            inductanceWithTolerance.set_nominal(filter->get_inductor().get_desired_inductance().get_nominal());
        }
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        designRequirements.set_operating_temperature(operatingTemp);
        designRequirements.set_application(Application::POWER);
        designRequirements.set_sub_application(SubApplication::POWER_FILTERING);
        return designRequirements;
    }

    inline std::complex<double> Zdc(const DcBusCapacitor& Capacitor, double omega,
                                    double Rs_supply = 0.0, double ESL = 0.0) {
        
        const std::complex<double> jω(0.0, omega);
        const double ESR = Capacitor.get_resistance().value_or(0.0);
        const double C   = Capacitor.get_capacitance();
        const std::complex<double> zC = (omega == 0.0) ? std::complex<double>(1e9, 0.0) : 1.0 / (jω * C);
        const std::complex<double> zL = jω * ESL;
        return std::complex<double>(Rs_supply + ESR, 0.0) + zL + zC; // series model
    }

    inline std::size_t next_pow2(std::size_t n) {
        std::size_t p = 1;
        while (p < n) p <<= 1;
        return p;
    }

    inline void bit_reverse(std::vector<std::complex<double>>& a) {
        const std::size_t n = a.size();
        std::size_t j = 0;
        for (std::size_t i = 1; i < n; ++i) {
            std::size_t bit = n >> 1;
            for (; j & bit; bit >>= 1) {
                j ^= bit;
            }
            j ^= bit;
            if (i < j) {
                std::swap(a[i], a[j]);
            }
        }
    }

    // In-place radix-2 Cooley–Tukey FFT. No scaling.
    inline void fft_pow2(std::vector<std::complex<double>>& a, bool inverse) {
        const std::size_t n = a.size();
        bit_reverse(a);
        for (std::size_t len = 2; len <= n; len <<= 1) {
            const double ang = 2.0 * M_PI / static_cast<double>(len) * (inverse ? +1.0 : -1.0);
            const std::complex<double> wlen(std::cos(ang), std::sin(ang));
            for (std::size_t i = 0; i < n; i += len) {
                std::complex<double> w(1.0, 0.0);
                const std::size_t half = len >> 1;
                for (std::size_t j = 0; j < half; ++j) {
                    const std::complex<double> u = a[i + j];
                    const std::complex<double> v = a[i + j + half] * w;
                    a[i + j]         = u + v;
                    a[i + j + half]  = u - v;
                    w *= wlen;
                }
            }
        }
    }

    // Bluestein’s algorithm for arbitrary N
    inline std::vector<std::complex<double>> bluestein_dft(const std::vector<std::complex<double>>& x, int sign) {
        const std::size_t N = x.size();
        if (N == 0) return {};
        if (N == 1) return x;

        const std::size_t M = next_pow2(2 * N - 1);

        // Chirps
        std::vector<std::complex<double>> a(N), b(M);
        const double s = static_cast<double>(sign);

        for (std::size_t n = 0; n < N; ++n) {
            const double ang = -s * M_PI * (static_cast<double>(n) * static_cast<double>(n)) / static_cast<double>(N);
            const std::complex<double> a_chirp(std::cos(ang), std::sin(ang));
            a[n] = x[n] * a_chirp;
        }

        b[0] = std::complex<double>(1.0, 0.0);
        for (std::size_t n = 1; n < N; ++n) {
            const double ang = +s * M_PI * (static_cast<double>(n) * static_cast<double>(n)) / static_cast<double>(N);
            const std::complex<double> val(std::cos(ang), std::sin(ang));
            b[n]      = val;
            b[M - n]  = val; // symmetry
        }

        // Zero-pad A to length M
        std::vector<std::complex<double>> A(M, std::complex<double>(0.0, 0.0));
        for (std::size_t n = 0; n < N; ++n) {
            A[n] = a[n];
        }

        // FFTs
        fft_pow2(A, false);
        fft_pow2(b, false);

        // Point-wise multiply
        for (std::size_t i = 0; i < M; ++i) {
            A[i] *= b[i];
        }

        // Inverse FFT (unnormalized)
        fft_pow2(A, true);

        // Scale by 1/M
        const double invM = 1.0 / static_cast<double>(M);
        std::vector<std::complex<double>> Y(N);
        for (std::size_t k = 0; k < N; ++k) {
            const double ang = -s * M_PI * (static_cast<double>(k) * static_cast<double>(k)) / static_cast<double>(N);
            const std::complex<double> out_chirp(std::cos(ang), std::sin(ang));
            Y[k] = A[k] * invM * out_chirp;
        }

        return Y;
    }

// ---- Drop-in replacements ---------------------------------------------------

    std::vector<std::complex<double>> MyInverter::compute_fft(const std::vector<double>& signal) {
        const std::size_t N = signal.size();
        std::vector<std::complex<double>> x(N);
        for (std::size_t n = 0; n < N; ++n) {
            x[n] = std::complex<double>(signal[n], 0.0);
        }

        // Forward DFT (sign=+1), then apply 1/N scaling to match your original
        std::vector<std::complex<double>> Y = bluestein_dft(x, +1);
        const double invN = (N > 0) ? (1.0 / static_cast<double>(N)) : 0.0;
        for (auto& v : Y) {
            v *= invN;
        }
        return Y;
    }

    inline std::vector<std::complex<double>> compute_inverse_fft(const std::vector<std::complex<double>>& X) {
        // Inverse DFT (sign=-1), leave unscaled (same as your old code)
        return bluestein_dft(X, -1);
    }

    std::complex<double> MyInverter::compute_load_impedance(const InverterLoad& load, double omega) {
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


    std::complex<double> MyInverter::compute_filter_impedance(const InverterDownstreamFilter& filter, double omega) {
        using namespace std::complex_literals;

        auto topology = filter.get_filter_topology();

        // === Build element impedances with ESR ===

        // Inductor 1
        double L1 = filter.get_inductor().get_desired_inductance().get_nominal().value();
        double ESR_L1 = filter.get_inductor().get_resistance().value_or(0.0);
        std::complex<double> ZL1(ESR_L1, omega * L1);

        // Capacitor (if present)
        std::complex<double> ZC;
        if (filter.get_capacitor()) {
            double C = filter.get_capacitor()->get_desired_capacitance();
            double ESR_C = filter.get_capacitor()->get_resistance().value_or(0.0);
            // Capacitor ESR is in series with reactive part
            ZC = std::complex<double>(ESR_C, -1.0 / (omega * C));
        }

        // Inductor 2 (if present)
        std::complex<double> ZL2;
        if (filter.get_inductor2()) {
            double L2 = filter.get_inductor2()->get_desired_inductance2().get_nominal().value();
            double ESR_L2 = filter.get_inductor2()->get_resistance().value_or(0.0);
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

    /// Helper: dq -> abc transformation
    ABCVoltages MyInverter::dq_to_abc(const std::complex<double>& Vdq, double theta) {
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
    std::pair<double,double> MyInverter::abc_to_alphabeta(const MyInverter::ABCVoltages& v) {
        double v_alpha = (2.0/3.0) * (v.Va - 0.5*v.Vb - 0.5*v.Vc);
        double v_beta  = (2.0/3.0) * ((sqrt(3)/2.0) * (v.Vb - v.Vc));
        return {v_alpha, v_beta};
    }

    ABCVoltages MyInverter::svpwm_modulation(const MyInverter::ABCVoltages& Vabc,
                                            double ma,
                                            double Vdc,
                                            double fsw) {
        auto [alphaRef, betaRef] = abc_to_alphabeta(Vabc);

        alphaRef *= ma;
        betaRef  *= ma;

        const double switchingPeriod = 1.0 / fsw;

        // Reference vector angle and magnitude
        double refAngleRad = std::atan2(betaRef, alphaRef);
        if (refAngleRad < 0) refAngleRad += 2.0 * M_PI;

        int sector = int(refAngleRad / (M_PI / 3.0)) + 1;
        if (sector > 6) sector = 6;

        const double refMagnitude = std::hypot(alphaRef, betaRef);
        const double angleInSector = refAngleRad - (sector - 1) * (M_PI / 3.0);

        // SVPWM dwell times
        double tVec1 = (switchingPeriod * std::sqrt(3.0) * refMagnitude / Vdc) * std::sin(M_PI/3.0 - angleInSector);
        double tVec2 = (switchingPeriod * std::sqrt(3.0) * refMagnitude / Vdc) * std::sin(angleInSector);
        double tZero = switchingPeriod - tVec1 - tVec2;
        if (tZero < 0) tZero = 0;

        // Leg duties
        double dutyA, dutyB, dutyC;
        switch (sector) {
            case 1: dutyA = (tVec1+tVec2+tZero/2)/switchingPeriod; dutyB = (tVec2+tZero/2)/switchingPeriod;       dutyC = (tZero/2)/switchingPeriod; break;
            case 2: dutyA = (tVec1+tZero/2)/switchingPeriod;       dutyB = (tVec1+tVec2+tZero/2)/switchingPeriod; dutyC = (tZero/2)/switchingPeriod; break;
            case 3: dutyA = (tZero/2)/switchingPeriod;             dutyB = (tVec1+tVec2+tZero/2)/switchingPeriod; dutyC = (tVec2+tZero/2)/switchingPeriod; break;
            case 4: dutyA = (tZero/2)/switchingPeriod;             dutyB = (tVec1+tZero/2)/switchingPeriod;       dutyC = (tVec1+tVec2+tZero/2)/switchingPeriod; break;
            case 5: dutyA = (tVec2+tZero/2)/switchingPeriod;       dutyB = (tZero/2)/switchingPeriod;             dutyC = (tVec1+tVec2+tZero/2)/switchingPeriod; break;
            case 6: dutyA = (tVec1+tVec2+tZero/2)/switchingPeriod; dutyB = (tZero/2)/switchingPeriod;             dutyC = (tVec1+tZero/2)/switchingPeriod; break;
            default: dutyA = dutyB = dutyC = 0.5; break;
        }

        // ✅ Convert duties back to voltages in [−Vdc/2, +Vdc/2]
        return {
            (dutyA - 0.5) * Vdc,
            (dutyB - 0.5) * Vdc,
            (dutyC - 0.5) * Vdc
        };
    }

    ABCVoltages MyInverter::compute_voltage_references(
            const TwoLevelInverter& inverter,
            const InverterOperatingPoint& op_point,
            const Modulation& modulation,
            double grid_angle_rad) 
    {
        using namespace std::complex_literals;

        const InverterLoad& load = op_point.get_load();
        double omega = 2.0 * M_PI * op_point.get_fundamental_frequency();

        std::complex<double> Vref_dq;

        switch (load.get_load_type()) {
            case LoadType::GRID: {
                double Vg_rms = load.get_phase_voltage().value_or(230.0);
                std::complex<double> Vg = std::complex<double>(Vg_rms, 0.0);
                std::complex<double> Zg(load.get_grid_resistance().value_or(0.0),
                                        omega * load.get_grid_inductance().value_or(0.0));
                double P = op_point.get_output_power().value_or(0.0);
                double pf = op_point.get_power_factor().value_or(1.0);
                double phi = acos(pf);
                std::complex<double> Iph = (P / (Vg_rms * pf)) * (cos(-phi) + 1i*sin(-phi));
                Vref_dq = Vg - Iph * Zg;
                break;
            }
            case LoadType::R_L: {
                std::complex<double> Zload(load.get_resistance().value_or(0.0),
                                        omega * load.get_inductance().value_or(0.0));
                double P = op_point.get_output_power().value_or(0.0);
                double pf = op_point.get_power_factor().value_or(1.0);
                double phi = acos(pf);
                std::complex<double> Iph = (P / pf) / inverter.get_line_rms_current().get_nominal().value()
                                        * (cos(-phi) + 1i*sin(-phi));
                Vref_dq = Iph * Zload;
                break;
            }
            default:
                throw std::runtime_error("Unknown load type");
        }

        // dq -> abc
        ABCVoltages Vabc = dq_to_abc(Vref_dq, grid_angle_rad);

        // apply modulation
        double ma = modulation.get_modulation_depth();
        switch (modulation.get_modulation_strategy()) {
            case ModulationStrategy::SPWM:
                Vabc.Va *= ma;
                Vabc.Vb *= ma;
                Vabc.Vc *= ma;
                break;
            case ModulationStrategy::THIPWM: {
                double k = modulation.get_third_harmonic_injection_coefficient().value_or(1.0/6.0);
                double sinA = std::sin(grid_angle_rad);
                double sinB = std::sin(grid_angle_rad - 2.0*M_PI/3.0);
                double sinC = std::sin(grid_angle_rad + 2.0*M_PI/3.0);
                double third = std::sin(3.0 * grid_angle_rad);
                double Vdc_nom = inverter.get_dc_bus_voltage().get_nominal().value();
                Vabc.Va = (Vdc_nom/2.0) * ma * (sinA + k * third);
                Vabc.Vb = (Vdc_nom/2.0) * ma * (sinB + k * third);
                Vabc.Vc = (Vdc_nom/2.0) * ma * (sinC + k * third);
                break;
            }
            case ModulationStrategy::SVPWM:
                Vabc = svpwm_modulation(Vabc, ma,
                                        inverter.get_dc_bus_voltage().get_nominal().value(),
                                        modulation.get_switching_frequency());
                break;
            default:
                throw std::runtime_error("Unknown modulation strategy");
        }

        // enforce single-phase after modulation
        if (inverter.get_number_of_legs() == 2) {
            Vabc.Vb = -Vabc.Va;
            Vabc.Vc = 0.0;
        }

        return Vabc;
    }

    double MyInverter::compute_carrier(const Modulation& modulation, double t) {
        const double switchingFrequency = modulation.get_switching_frequency();
        const double switchingPeriod    = 1.0 / switchingFrequency;

        // phase inside the present PWM period, in [0,1)
        const double phaseInPeriod = std::fmod(t, switchingPeriod) / switchingPeriod;

        switch (modulation.get_pwm_type()) {
            case PwmType::SAWTOOTH:
                return 2.0 * phaseInPeriod - 1.0;   // -1 → +1
            case PwmType::TRIANGULAR:
                return (phaseInPeriod < 0.5)
                    ?  4.0 * phaseInPeriod - 1.0   // -1 → +1
                    : -4.0 * phaseInPeriod + 3.0;  // +1 → -1
            default:
                throw std::runtime_error("Unknown PWM carrier type");
        }
    }

    PwmSignals MyInverter::compare_with_carrier(const MyInverter::ABCVoltages& Vabc,
                                                double carrier,
                                                double Vdc,
                                                const Modulation& modulation) {
        // ---- Duty cycles (D) ----
        auto toDuty = [&](double v_leg)->double {
            // map leg ref to [0,1] duty wrt ±Vdc/2
            double duty = 0.5 * (v_leg / (Vdc / 2.0) + 1.0);
            return std::clamp(duty, 0.0, 1.0);
        };

        double dutyA = toDuty(Vabc.Va);
        double dutyB = toDuty(Vabc.Vb);
        double dutyC = toDuty(Vabc.Vc);

        // Deadtime & rise-time corrections expressed as duty fractions
        const double switchingPeriod = 1.0 / modulation.get_switching_frequency();

        if (modulation.get_deadtime()) {
            const double deadtimeFraction = *modulation.get_deadtime() / switchingPeriod;
            dutyA = std::clamp(dutyA - deadtimeFraction, 0.0, 1.0);
            dutyB = std::clamp(dutyB - deadtimeFraction, 0.0, 1.0);
            dutyC = std::clamp(dutyC - deadtimeFraction, 0.0, 1.0);
        }
        if (modulation.get_rise_time()) {
            const double riseFraction = *modulation.get_rise_time() / switchingPeriod;
            dutyA = std::clamp(dutyA - 0.5*riseFraction, 0.0, 1.0);
            dutyB = std::clamp(dutyB - 0.5*riseFraction, 0.0, 1.0);
            dutyC = std::clamp(dutyC - 0.5*riseFraction, 0.0, 1.0);
        }

        // ---- Comparator inputs (normalized back to [-1, +1]) ----
        const double compA = 2.0 * dutyA - 1.0;
        const double compB = 2.0 * dutyB - 1.0;
        const double compC = 2.0 * dutyC - 1.0;

        // ---- Switch states (S) for *upper* devices ----
        const bool gateUpperAOn = (compA > carrier);
        const bool gateUpperBOn = (compB > carrier);
        const bool gateUpperCOn = (compC > carrier);

        return { gateUpperAOn, gateUpperBOn, gateUpperCOn };
    }

    NodeResult MyInverter::solve_filter_topology(const InverterDownstreamFilter& filter,
                                    const InverterLoad& load,
                                    double omega,
                                    std::complex<double> Vinv)
    {
        using namespace std::complex_literals;

        // --- Build element impedances ---
        std::complex<double> ZL1(filter.get_inductor().get_resistance().value_or(0.0),
                                omega * filter.get_inductor().get_desired_inductance().get_nominal().value());

        std::complex<double> ZC(1e9, 0.0); // open if no C
        if (filter.get_capacitor()) {
            double C = filter.get_capacitor()->get_desired_capacitance();
            double ESRc = filter.get_capacitor()->get_resistance().value_or(0.0);
            ZC = std::complex<double>(ESRc, -1.0/(omega*C));
        }

        std::complex<double> ZL2(1e9, 0.0); // open if no L2
        if (filter.get_inductor2()) {
            double L2 = filter.get_inductor2()->get_desired_inductance2().get_nominal().value();
            double ESR2 = filter.get_inductor2()->get_resistance().value_or(0.0);
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

    HarmonicsBundle MyInverter::compute_harmonics(const Modulation& modulation,
                                                const InverterOperatingPoint& op_point,
                                                const ABCVoltages& Vabc_ref,
                                                double Vdc_nom,
                                                std::complex<double> Vfund,
                                                std::complex<double> Ifund,
                                                double f1,
                                                const InverterDownstreamFilter& filter,
                                                const InverterLoad& load,
                                                int Nperiods,
                                                int samplesPerPeriod)
    {
    const double fsw = modulation.get_switching_frequency();
    const double fs  = fsw * samplesPerPeriod;                 // Hz
    const int    N   = static_cast<int>(std::llround(Nperiods * fs / f1)); // simulate Nperiods fundamentals
    const double Ts  = 1.0 / fsw;                              // keep if you use Ts elsewhere

    std::vector<int>    sa(N), sb(N), sc(N);
    std::vector<double> va(N), vb(N), vc(N);
    std::vector<double> referenceA(N), referenceB(N), referenceC(N), carrierSig(N);

    ABCVoltages Vabc_use = Vabc_ref;

    #ifdef DEBUG_PLOTS
    std::vector<int> gateA(N), gateB(N), gateC(N);
    #endif

    for (int n = 0; n < N; ++n) {
        double t = static_cast<double>(n) / fs;
        double theta = 2.0 * M_PI * f1 * t + op_point.get_current_phase_angle().value_or(0.0);

        // compute instantaneous references
        ABCVoltages Vabc_t = compute_voltage_references(
            *this, op_point, modulation, theta);

        double carrier = compute_carrier(modulation, t);

        referenceA[n] = Vabc_t.Va / (Vdc_nom/2.0);
        referenceB[n] = Vabc_t.Vb / (Vdc_nom/2.0);
        referenceC[n] = Vabc_t.Vc / (Vdc_nom/2.0);
        carrierSig[n] = carrier;
        PwmSignals gates = compare_with_carrier(Vabc_t, carrier, Vdc_nom, modulation);

    #ifdef DEBUG_PLOTS
        gateA[n] = gates.gateUpperAOn ? 1 : 0;
        gateB[n] = gates.gateUpperBOn ? 1 : 0;
        gateC[n] = gates.gateUpperCOn ? 1 : 0;
    #endif

        sa[n] = gates.gateUpperAOn ? +1 : -1;
        sb[n] = gates.gateUpperBOn ? +1 : -1;
        sc[n] = gates.gateUpperCOn ? +1 : -1;

        va[n] = 0.5 * Vdc_nom * sa[n];
        vb[n] = 0.5 * Vdc_nom * sb[n];
        vc[n] = 0.5 * Vdc_nom * sc[n];
    }

    std::vector<std::complex<double>> Vin_a(N), Vin_b(N), Vin_c(N);

    if (get_number_of_legs() == 3) {
        // --- 1) FFT of pole voltages ---
        auto Va = compute_fft(va);
        auto Vb = compute_fft(vb);
        auto Vc = compute_fft(vc);

        for (int k = 0; k < N; ++k) {
            auto V0 = (Va[k] + Vb[k] + Vc[k]) / 3.0;
            Vin_a[k] = Va[k] - V0;
            Vin_b[k] = Vb[k] - V0;
            Vin_c[k] = Vc[k] - V0;
        }
    } else {
        // ✅ single-phase: FFT directly on line-to-line time waveform
        std::vector<double> v_ab_time(N);
        for (int n = 0; n < N; ++n) {
            v_ab_time[n] = va[n] - vb[n];
        }
        auto Vab = compute_fft(v_ab_time);

        for (int k = 0; k < N; ++k) {
            Vin_a[k] =  0.5 * Vab[k];
            Vin_b[k] = -0.5 * Vab[k];
            Vin_c[k] = std::complex<double>(0.0, 0.0);
        }
    }

        // --- 3) Solve filter per harmonic to get L1 current spectra ---
        std::vector<std::complex<double>> IL1_a(N), IL1_b(N), IL1_c(N);
        for (int k = 0; k < N; ++k) {
            double f     = k * fs / N;
            double omega = 2.0 * M_PI * f;

            NodeResult na = solve_filter_topology(filter, load, omega, Vin_a[k]);
            IL1_a[k] = na.iL1;

            if (get_number_of_legs() == 3) {
                NodeResult nb = solve_filter_topology(filter, load, omega, Vin_b[k]);
                NodeResult nc = solve_filter_topology(filter, load, omega, Vin_c[k]);
                IL1_b[k] = nb.iL1;
                IL1_c[k] = nc.iL1;
            } else {
                IL1_b[k] = -IL1_a[k];
                IL1_c[k] = std::complex<double>(0.0, 0.0);
            }
        }

        // --- 4) Inverse FFT currents to time domain, build instantaneous input power p(t) ---
        auto iLa_t_c = compute_inverse_fft(IL1_a);
        auto iLb_t_c = compute_inverse_fft(IL1_b);
        auto iLc_t_c = compute_inverse_fft(IL1_c);

        std::vector<double> p_t(N);
        if (get_number_of_legs() == 3) {
            for (int n = 0; n < N; ++n) {
                double pin = va[n] * iLa_t_c[n].real()
                        + vb[n] * iLb_t_c[n].real()
                        + vc[n] * iLc_t_c[n].real();
                p_t[n] = pin;
            }
        } else {
            for (int n=0; n<N; ++n) {
                double v_ab = va[n] - vb[n];
                double i_line = iLa_t_c[n].real();
                p_t[n] = v_ab * i_line;
            }
        }
        double p_avg = std::accumulate(p_t.begin(), p_t.end(), 0.0) / N;
        for (auto& x : p_t) x -= p_avg;

        // --- 5) DC ripple ---
        auto Pfft = compute_fft(p_t);
        std::vector<std::complex<double>> Idc(N), Vdc_fft(N);
        for (int k = 0; k < N; ++k) {
            double f     = k * fs / N;
            double omega = 2.0 * M_PI * f;
            if (k == 0) { Idc[k] = Vdc_fft[k] = 0.0; continue; }
            Idc[k]     = Pfft[k] / std::complex<double>(Vdc_nom, 0.0);
            Vdc_fft[k] = Zdc(get_dc_bus_capacitor(), omega) * Idc[k];
        }
        auto vdc_t_c = compute_inverse_fft(Vdc_fft);
        std::vector<double> vdc_t(N);
        for (int n = 0; n < N; ++n) vdc_t[n] = vdc_t_c[n].real();

        // --- 6) Remodulate pole voltages with ripple ---
        for (int n = 0; n < N; ++n) {
            va[n] = 0.5 * (Vdc_nom + vdc_t[n]) * sa[n];
            vb[n] = 0.5 * (Vdc_nom + vdc_t[n]) * sb[n];
            vc[n] = 0.5 * (Vdc_nom + vdc_t[n]) * sc[n];
        }

        // --- 7) Final harmonic bundle ---
        HarmonicsBundle bundle;
        const double fmax = 5.0 * fsw;
        for (int k = 0; k < N/2; ++k) {
            double f     = k * fs / N;
            if (f > fmax) break;
            double omega = 2.0 * M_PI * f;

            NodeResult na = solve_filter_topology(filter, load, omega, Vin_a[k]);
            bundle.Vharm.get_mutable_frequencies().push_back(f);
            bundle.Vharm.get_mutable_amplitudes().push_back(std::abs(na.vNode));

            bundle.Iharm.get_mutable_frequencies().push_back(f);
            bundle.Iharm.get_mutable_amplitudes().push_back(std::abs(na.iL1));
        }

        // Override fundamental bins
        auto& freqs = bundle.Vharm.get_mutable_frequencies();
        auto& vamps = bundle.Vharm.get_mutable_amplitudes();
        if (!freqs.empty()) {
            auto vIt = std::min_element(freqs.begin(), freqs.end(),
                [f1](double a, double b) { return std::abs(a - f1) < std::abs(b - f1); });
            if (vIt != freqs.end()) vamps[std::distance(freqs.begin(), vIt)] = std::abs(Vfund);
        }
        auto& ifreqs = bundle.Iharm.get_mutable_frequencies();
        auto& iamps  = bundle.Iharm.get_mutable_amplitudes();
        if (!ifreqs.empty()) {
            auto iIt = std::min_element(ifreqs.begin(), ifreqs.end(),
                [f1](double a, double b) { return std::abs(a - f1) < std::abs(b - f1); });
            if (iIt != ifreqs.end()) iamps[std::distance(ifreqs.begin(), iIt)] = std::abs(Ifund);
        }

#ifdef DEBUG_PLOTS
    debug_plots::init_folder();

    // Carrier vs Reference (1 fundamental period)
    debug_plots::plot_carrier_vs_refs(
    carrierSig, referenceA, referenceB, referenceC, f1);

    debug_plots::plot_pwm_signals(gateA, gateB, gateC, f1, fs);

    // Comparator outputs va,vb,vc over 2 Fsw
    debug_plots::plot_va_vb_vc_short(va, vb, vc, fsw, fs);

    // Comparator outputs va,vb,vc over full fundamental
    debug_plots::plot_va_vb_vc_fundamental(va, vb, vc, f1, fs);

    // FFT of vL1 & iL1 (initial)
    debug_plots::plot_fft_vl1_il1(bundle.Vharm.get_frequencies(),
                                  bundle.Vharm.get_amplitudes(),
                                  bundle.Iharm.get_amplitudes());

    // Instantaneous power
    debug_plots::plot_power(p_t, f1);

    // DC link ripple
    debug_plots::plot_vdc_ripple(vdc_t, f1);

    // Final FFT of vL1 & iL1 (after remodulation)
    debug_plots::plot_final_fft_vl1_il1(bundle.Vharm.get_frequencies(),
                                        bundle.Vharm.get_amplitudes(),
                                        bundle.Iharm.get_amplitudes());
#endif

        return bundle;
    }

    std::vector<OperatingPoint> MyInverter::process_operating_points() {
        std::vector<OperatingPoint> operatingPointsResult;

        for (auto& op_point : this->inverterOperatingPoints) {
            double f1 = op_point.get_fundamental_frequency();
            double omega = 2.0 * M_PI * f1;

            // --- Step 1: Compute dq reference → abc
            ABCVoltages Vabc = compute_voltage_references(
                *this,
                op_point,
                this->get_modulation().value(),
                op_point.get_current_phase_angle().value()
            );

            // --- Step 2: Fundamental phasors (clean values)
            std::complex<double> Vfund = std::polar(
                op_point.get_load().get_phase_voltage().value_or(230.0),
                op_point.get_current_phase_angle().value_or(0.0)
            );

            // Solve full filter/load topology at fundamental
            NodeResult node = solve_filter_topology(
                this->get_downstream_filter().value(),
                op_point.get_load(),
                omega,
                Vfund
            );

            std::complex<double> Ifund = node.iL1;

            // --- Step 3: Harmonic analysis (both V and I at once)
            HarmonicsBundle bundle = compute_harmonics(
                this->get_modulation().value(),
                op_point,
                Vabc,
                this->get_dc_bus_voltage().get_nominal().value(),
                Vfund,
                Ifund,
                f1,
                this->get_downstream_filter().value(),  // filter
                op_point.get_load()
            );

            // --- Step 4: Signal descriptors
            SignalDescriptor voltageSig;
            voltageSig.set_harmonics(bundle.Vharm);

            SignalDescriptor currentSig;
            currentSig.set_harmonics(bundle.Iharm);

            // --- Step 5: Assemble excitation
            OperatingPointExcitation excitation;
            excitation.set_voltage(voltageSig);
            excitation.set_current(currentSig);

            // --- Step 6: Build operating point result
            OperatingPoint result;
            result.set_excitations_per_winding({excitation});
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
