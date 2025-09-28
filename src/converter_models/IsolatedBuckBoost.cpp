#include "converter_models/IsolatedBuckBoost.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>

namespace OpenMagnetics {

    double IsolatedBuckBoost::calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency) {
        auto dutyCycle = outputVoltage / (inputVoltage + outputVoltage) * efficiency;
        if (dutyCycle >= 1) {
            throw std::runtime_error("Duty cycle must be smaller than 1");
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
            throw std::runtime_error("At least one operating point is needed");
        }
        for (size_t isolatedbuckBoostOperatingPointIndex = 1; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
            if (get_operating_points()[isolatedbuckBoostOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw std::runtime_error("Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[isolatedbuckBoostOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
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

            for (size_t isolatedbuckBoostOperatingPointIndex = 0; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
                auto isolatedbuckBoostOperatingPoint = get_operating_points()[isolatedbuckBoostOperatingPointIndex];
                maximumOutputCurrent = std::max(maximumOutputCurrent, isolatedbuckBoostOperatingPoint.get_output_currents()[0]);
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
} // namespace OpenMagnetics
