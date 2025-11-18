#include "converter_models/Buck.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>

namespace OpenMagnetics {

    double Buck::calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency) {
        auto dutyCycle = (outputVoltage + diodeVoltageDrop) / ((inputVoltage + diodeVoltageDrop) * efficiency);
        if (dutyCycle >= 1) {
            throw std::runtime_error("Duty cycle must be smaller than 1");
        }
        return dutyCycle;
    }

    Buck::Buck(const json& j) {
        from_json(j, *this);
    }

    AdvancedBuck::AdvancedBuck(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint Buck::process_operating_points_for_input_voltage(double inputVoltage, BuckOperatingPoint outputOperatingPoint, double inductance) {

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
        auto primaryVoltaveMinimum = -outputVoltage -diodeVoltageDrop;
        auto primaryVoltaveMaximum = inputVoltage - outputVoltage;
        auto primaryVoltavePeaktoPeak = primaryVoltaveMaximum - primaryVoltaveMinimum;

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
                voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);

            }
            else {
                currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, outputCurrent, 0);
                voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, 0);
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
            throw std::runtime_error("At least one operating point is needed");
        }
        if (!get_input_voltage().get_nominal() && !get_input_voltage().get_maximum() && !get_input_voltage().get_minimum()) {
            if (!assert) {
                return false;
            }
            throw std::runtime_error("No input voltage introduced");
        }

        return true;
    }

    DesignRequirements Buck::process_design_requirements() {
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Missing both current ripple ratio and maximum swtich current");
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

    std::vector<OperatingPoint> Buck::process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) {
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

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
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
} // namespace OpenMagnetics
