#include "converter_models/Flyback.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
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
            double secondaryVoltageOffset = maximumSecondaryVoltage + minimumSecondaryVoltage;
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

                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, maximumNeededInductance, std::nullopt, customDutyCycle, customDeadTime);

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
