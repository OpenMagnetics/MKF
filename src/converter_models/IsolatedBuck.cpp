#include "converter_models/IsolatedBuck.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>

namespace OpenMagnetics {

    double IsolatedBuck::calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency) {
        auto dutyCycle = outputVoltage / inputVoltage * efficiency;
        if (dutyCycle >= 1) {
            throw std::runtime_error("Duty cycle must be smaller than 1");
        }
        return dutyCycle;
    }

    IsolatedBuck::IsolatedBuck(const json& j) {
        from_json(j, *this);
    }

    AdvancedIsolatedBuck::AdvancedIsolatedBuck(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint IsolatedBuck::processOperatingPointsForInputVoltage(double inputVoltage, IsolatedBuckOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance) {

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
        auto period = 1.0 / switchingFrequency;

        auto tOn = dutyCycle * period;
        auto magnetizingCurrentRipple = (inputVoltage - primaryOutputVoltage) * tOn / inductance;
        auto primaryCurrentMaximum = primaryOutputCurrent + totalReflectedSecondaryCurrent + magnetizingCurrentRipple / 2;
        auto primaryCurrentMinimum = primaryOutputCurrent - totalReflectedSecondaryCurrent * (2 * dutyCycle) / (1 - dutyCycle) - magnetizingCurrentRipple / 2;
        auto primaryCurrentPeakToPeak = primaryCurrentMaximum - primaryCurrentMinimum;

        auto primaryVoltaveMaximum = inputVoltage - primaryOutputVoltage;
        auto primaryVoltaveMinimum = -primaryOutputVoltage;
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

            auto secondaryVoltaveMaximum = (inputVoltage - primaryOutputVoltage) / turnsRatios[secondaryIndex] - diodeVoltageDrop;
            auto secondaryVoltaveMinimum = -primaryOutputVoltage / turnsRatios[secondaryIndex] + diodeVoltageDrop;

            // Current
            {
                std::vector<double> data = {
                    0,
                    0,
                    secondaryOutputCurrent + secondaryCurrentMinimum,
                    secondaryOutputCurrent + secondaryCurrentMaximum
                };
                std::vector<double> time = {
                    0,
                    tOn,
                    tOn,
                    period
                };
                currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                currentWaveform.set_data(data);
                currentWaveform.set_time(time);
            }
            // Voltage
            {
                std::vector<double> data = {
                    secondaryVoltaveMinimum,
                    secondaryVoltaveMinimum,
                    secondaryVoltaveMaximum,
                    secondaryVoltaveMaximum
                };
                std::vector<double> time = {
                    0,
                    tOn,
                    tOn,
                    period
                };
                voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                voltageWaveform.set_data(data);
                voltageWaveform.set_time(time);
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

    bool IsolatedBuck::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw std::runtime_error("At least one operating point is needed");
        }
        for (size_t isolatedbuckOperatingPointIndex = 1; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
            if (get_operating_points()[isolatedbuckOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw std::runtime_error("Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[isolatedbuckOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
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

    DesignRequirements IsolatedBuck::process_design_requirements() {
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Missing both current ripple ratio and maximum swtich current");
        }

        // Turns ratio calculation
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size() - 1, 0);
        for (auto isolatedbuckOperatingPoint : get_operating_points()) {
            double primaryVoltage = isolatedbuckOperatingPoint.get_output_voltages()[0];
            for (size_t secondaryIndex = 0; secondaryIndex < isolatedbuckOperatingPoint.get_output_voltages().size() - 1; ++secondaryIndex) {
                auto turnsRatio = primaryVoltage / (isolatedbuckOperatingPoint.get_output_voltages()[secondaryIndex + 1] + get_diode_voltage_drop());
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
                maximumOutputCurrent = std::max(maximumOutputCurrent, isolatedbuckOperatingPoint.get_output_currents()[0]);
            }

            maximumCurrentRiple = currentRippleRatio * maximumOutputCurrent;
        }

        if (get_maximum_switch_current()) {
            for (auto isolatedbuckOperatingPoint : get_operating_points()) {
                auto primaryCurrent = isolatedbuckOperatingPoint.get_output_currents()[0];
                double totalReflectedSecondaryCurrent = 0;
                for (size_t secondaryIndex = 0; secondaryIndex < isolatedbuckOperatingPoint.get_output_currents().size() - 1; ++secondaryIndex) {
                    totalReflectedSecondaryCurrent += isolatedbuckOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex];
                }

                maximumCurrentRiple = (get_maximum_switch_current().value() - primaryCurrent - totalReflectedSecondaryCurrent) * 2;
            }
        }

        double maximumNeededInductance = 0;
        for (size_t isolatedbuckOperatingPointIndex = 0; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
            auto isolatedbuckOperatingPoint = get_operating_points()[isolatedbuckOperatingPointIndex];
            double primaryVoltage = isolatedbuckOperatingPoint.get_output_voltages()[0];
            auto switchingFrequency = isolatedbuckOperatingPoint.get_switching_frequency();
            auto desiredInductance =  (maximumInputVoltage - primaryVoltage) * primaryVoltage / (maximumInputVoltage * switchingFrequency * maximumCurrentRiple);
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
        designRequirements.set_topology(Topologies::ISOLATED_BUCK_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> IsolatedBuck::process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) {
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
            for (size_t isolatedbuckOperatingPointIndex = 0; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[isolatedbuckOperatingPointIndex], turnsRatios, magnetizingInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(isolatedbuckOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> IsolatedBuck::process_operating_points(Magnetic magnetic) {
        IsolatedBuck::run_checks(_assertErrors);

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedIsolatedBuck::process() {
        IsolatedBuck::run_checks(_assertErrors);

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
        designRequirements.set_topology(Topologies::ISOLATED_BUCK_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t isolatedbuckOperatingPointIndex = 0; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[isolatedbuckOperatingPointIndex], turnsRatios, maximumNeededInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(isolatedbuckOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }
} // namespace OpenMagnetics
