#include "converter_models/PushPull.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/CircuitSimulatorInterface.h"
#include "support/Utils.h"
#include <cfloat>
#include <sstream>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double get_total_reflected_secondary_current(const PushPullOperatingPoint& pushPullOperatingPoint, const std::vector<double>& turnsRatios) {
        double totalReflectedSecondaryCurrent = 0;
        // Main secondary
        totalReflectedSecondaryCurrent += pushPullOperatingPoint.get_output_currents()[0] / turnsRatios[1];
        for (size_t secondaryIndex = 0; secondaryIndex < pushPullOperatingPoint.get_output_voltages().size() - 1; ++secondaryIndex) {
            totalReflectedSecondaryCurrent += pushPullOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex + 3];
        }

        return totalReflectedSecondaryCurrent;
    }

    std::vector<double> convert_turns_ratios(const std::vector<double>& turnsRatios) {
        std::vector<double> newTurnsRatios;

        // Second Primary
        {
            newTurnsRatios.push_back(1);
        }
        // First Secondary
        {
            newTurnsRatios.push_back(turnsRatios[0]);
        }
        // Second Secondary
        {
            newTurnsRatios.push_back(turnsRatios[0]);
        }
        for (size_t auxiliarySecondaryIndex = 1; auxiliarySecondaryIndex < turnsRatios.size(); ++auxiliarySecondaryIndex) {
            auto turnsRatio = turnsRatios[auxiliarySecondaryIndex];
            newTurnsRatios.push_back(turnsRatio);
        }

        return newTurnsRatios;
    }

    double PushPull::get_maximum_duty_cycle() {
        if (get_duty_cycle()) {
            return get_duty_cycle().value();
        }
        else {
            return 0.5;
        }
    }

    PushPull::PushPull(const json& j) {
        from_json(j, *this);
    }

    AdvancedPushPull::AdvancedPushPull(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint PushPull::process_operating_points_for_input_voltage(double inputVoltage, const PushPullOperatingPoint& outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance, double outputInductance) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double mainOutputVoltage = outputOperatingPoint.get_output_voltages()[0];
        double mainOutputCurrent = outputOperatingPoint.get_output_currents()[0];
        double diodeVoltageDrop = get_diode_voltage_drop();

        double mainSecondaryTurnsRatio = turnsRatios[1];

        auto inductorCurrentRipple = get_current_ripple_ratio() * mainOutputCurrent;
        auto period = 1.0 / switchingFrequency;
        auto t1 = period / 2 * (mainOutputVoltage + diodeVoltageDrop) / (inputVoltage / mainSecondaryTurnsRatio);
        if (t1 > period / 2) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "T1 cannot be larger than period/2, wrong topology configuration");
        }

        auto magnetizationCurrent = inputVoltage * t1 / inductance;
        auto minimumSecondaryCurrent = mainOutputCurrent - inductorCurrentRipple / 2;
        auto maximumSecondaryCurrent = mainOutputCurrent + inductorCurrentRipple / 2;
        auto minimumPrimaryCurrent = minimumSecondaryCurrent / mainSecondaryTurnsRatio - magnetizationCurrent / 2;
        auto maximumPrimaryCurrent = minimumSecondaryCurrent / mainSecondaryTurnsRatio + magnetizationCurrent / 2;

        for (size_t auxiliarySecondaryIndex = 1; auxiliarySecondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++auxiliarySecondaryIndex) {

            auto auxiliaryInductorCurrentRipple = get_current_ripple_ratio() * outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex];
            auto minimumAuxiliarySecondaryCurrent = outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex] - auxiliaryInductorCurrentRipple / 2;
            auto maximumAuxiliarySecondaryCurrent = outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex] + auxiliaryInductorCurrentRipple / 2;
            size_t turnsRatioAuxiliarySecondaryIndex = 2 + auxiliarySecondaryIndex;
            minimumPrimaryCurrent += minimumAuxiliarySecondaryCurrent / turnsRatios[turnsRatioAuxiliarySecondaryIndex];
            maximumPrimaryCurrent += maximumAuxiliarySecondaryCurrent / turnsRatios[turnsRatioAuxiliarySecondaryIndex];
        }


        if (minimumPrimaryCurrent > 0) {  // CCM
            double minimumPrimarySideTransformerCurrent = minimumPrimaryCurrent;
            double maximumPrimarySideTransformerCurrent = maximumPrimaryCurrent;
            double minimumPrimarySideTransformerVoltage = -inputVoltage;
            double maximumPrimarySideTransformerVoltage = inputVoltage;

            double minimumSecondarySideTransformerCurrentT1OfFET = minimumSecondaryCurrent;
            double maximumSecondarySideTransformerCurrentT1OfFET = maximumSecondaryCurrent;
            double minimumSecondarySideTransformerCurrentT2OtherFET = (minimumSecondaryCurrent / mainSecondaryTurnsRatio + magnetizationCurrent / 2) / 2 * mainSecondaryTurnsRatio - inductorCurrentRipple / 2;
            double maximumSecondarySideTransformerCurrentT2OtherFET = (minimumSecondaryCurrent / mainSecondaryTurnsRatio + magnetizationCurrent / 2) / 2 * mainSecondaryTurnsRatio;
            double minimumSecondarySideTransformerCurrentT2OfFET = minimumSecondaryCurrent - minimumSecondarySideTransformerCurrentT2OtherFET;
            double maximumSecondarySideTransformerCurrentT2OfFET = maximumSecondaryCurrent - maximumSecondarySideTransformerCurrentT2OtherFET;
            double minimumSecondarySideTransformerVoltage = -inputVoltage / mainSecondaryTurnsRatio;
            double maximumSecondarySideTransformerVoltage = inputVoltage / mainSecondaryTurnsRatio;
            // First Primary
            {
                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        minimumPrimarySideTransformerCurrent,
                        maximumPrimarySideTransformerCurrent,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        maximumPrimarySideTransformerVoltage,
                        maximumPrimarySideTransformerVoltage,
                        0,
                        0,
                        minimumPrimarySideTransformerVoltage,
                        minimumPrimarySideTransformerVoltage,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "First primary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }

            // Second Primary
            {
                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        0,
                        0,
                        minimumPrimarySideTransformerCurrent,
                        maximumPrimarySideTransformerCurrent,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        minimumPrimarySideTransformerVoltage,
                        minimumPrimarySideTransformerVoltage,
                        0,
                        0,
                        maximumPrimarySideTransformerVoltage,
                        maximumPrimarySideTransformerVoltage,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Second primary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }

            // First Main Secondary
            {
                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        0,
                        0,
                        maximumSecondarySideTransformerCurrentT2OtherFET,
                        minimumSecondarySideTransformerCurrentT2OtherFET,
                        minimumSecondarySideTransformerCurrentT1OfFET,
                        maximumSecondarySideTransformerCurrentT1OfFET,
                        maximumSecondarySideTransformerCurrentT2OfFET,
                        minimumSecondarySideTransformerCurrentT2OfFET
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        minimumSecondarySideTransformerVoltage,
                        minimumSecondarySideTransformerVoltage,
                        0,
                        0,
                        maximumSecondarySideTransformerVoltage,
                        maximumSecondarySideTransformerVoltage,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "First secondary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }

            // Second Main Secondary
            {
                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        minimumSecondarySideTransformerCurrentT1OfFET,
                        maximumSecondarySideTransformerCurrentT1OfFET,
                        maximumSecondarySideTransformerCurrentT2OfFET,
                        minimumSecondarySideTransformerCurrentT2OfFET,
                        0,
                        0,
                        maximumSecondarySideTransformerCurrentT2OtherFET,
                        minimumSecondarySideTransformerCurrentT2OtherFET
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        maximumSecondarySideTransformerVoltage,
                        maximumSecondarySideTransformerVoltage,
                        0,
                        0,
                        minimumSecondarySideTransformerVoltage,
                        minimumSecondarySideTransformerVoltage,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Second secondary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }

            for (size_t auxiliarySecondaryIndex = 1; auxiliarySecondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++auxiliarySecondaryIndex) {

                auto auxiliaryInductorCurrentRipple = get_current_ripple_ratio() * outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex];
                auto minimumAuxiliarySecondaryCurrent = outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex] - auxiliaryInductorCurrentRipple / 2;
                auto maximumAuxiliarySecondaryCurrent = outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex] + auxiliaryInductorCurrentRipple / 2;
                size_t turnsRatioAuxiliarySecondaryIndex = 2 + auxiliarySecondaryIndex;
                double turnsRatioAuxiliarySecondary = turnsRatios[turnsRatioAuxiliarySecondaryIndex];
                minimumPrimaryCurrent += minimumAuxiliarySecondaryCurrent / turnsRatioAuxiliarySecondary;
                maximumPrimaryCurrent += maximumAuxiliarySecondaryCurrent / turnsRatioAuxiliarySecondary;

                double minimumAuxiliarySecondarySideTransformerCurrentT1OfFET = minimumAuxiliarySecondaryCurrent;
                double maximumAuxiliarySecondarySideTransformerCurrentT1OfFET = maximumAuxiliarySecondaryCurrent;
                double minimumAuxiliarySecondarySideTransformerCurrentT2OtherFET = (minimumAuxiliarySecondaryCurrent / turnsRatioAuxiliarySecondary + magnetizationCurrent / 2) / 2 * turnsRatioAuxiliarySecondary - inductorCurrentRipple / 2;
                double maximumAuxiliarySecondarySideTransformerCurrentT2OtherFET = (minimumAuxiliarySecondaryCurrent / turnsRatioAuxiliarySecondary + magnetizationCurrent / 2) / 2 * turnsRatioAuxiliarySecondary;
                double minimumAuxiliarySecondarySideTransformerCurrentT2OfFET = minimumAuxiliarySecondaryCurrent - minimumAuxiliarySecondarySideTransformerCurrentT2OtherFET;
                double maximumAuxiliarySecondarySideTransformerCurrentT2OfFET = maximumAuxiliarySecondaryCurrent - maximumAuxiliarySecondarySideTransformerCurrentT2OtherFET;
                double minimumAuxiliarySecondarySideTransformerVoltage = -inputVoltage / turnsRatioAuxiliarySecondary;
                double maximumAuxiliarySecondarySideTransformerVoltage = inputVoltage / turnsRatioAuxiliarySecondary;

                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        minimumAuxiliarySecondarySideTransformerCurrentT1OfFET,
                        maximumAuxiliarySecondarySideTransformerCurrentT1OfFET,
                        maximumAuxiliarySecondarySideTransformerCurrentT2OfFET,
                        minimumAuxiliarySecondarySideTransformerCurrentT2OfFET,
                        0,
                        0,
                        maximumAuxiliarySecondarySideTransformerCurrentT2OtherFET,
                        minimumAuxiliarySecondarySideTransformerCurrentT2OtherFET
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        maximumAuxiliarySecondarySideTransformerVoltage,
                        maximumAuxiliarySecondarySideTransformerVoltage,
                        0,
                        0,
                        minimumAuxiliarySecondarySideTransformerVoltage,
                        minimumAuxiliarySecondarySideTransformerVoltage,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Auxiliary " + std::to_string(auxiliarySecondaryIndex));
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);

            }

        }
        else {  // DCM
            auto t1 = sqrt(2 * mainOutputCurrent * outputInductance * (mainOutputVoltage + diodeVoltageDrop) / (2 * switchingFrequency * (inputVoltage / mainSecondaryTurnsRatio - diodeVoltageDrop - mainOutputVoltage) * (inputVoltage / mainSecondaryTurnsRatio)));
            auto t2 = t1 * (inputVoltage / mainSecondaryTurnsRatio) / (mainOutputVoltage + diodeVoltageDrop) - t1;
            if (t1 + t2 > period / 2) {
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "T1 + T2 cannot be larger than period/2, wrong topology configuration");
            }

            auto minimumSecondaryCurrent = 0;
            auto maximumSecondaryCurrent = inductorCurrentRipple;
            auto minimumPrimaryCurrent = 0;
            auto maximumPrimaryCurrent = inductorCurrentRipple / mainSecondaryTurnsRatio + magnetizationCurrent / 2;

            for (size_t auxiliarySecondaryIndex = 1; auxiliarySecondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++auxiliarySecondaryIndex) {
                auto auxiliaryInductorCurrentRipple = get_current_ripple_ratio() * outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex];
                size_t turnsRatioAuxiliarySecondaryIndex = 2 + auxiliarySecondaryIndex;
                maximumPrimaryCurrent += auxiliaryInductorCurrentRipple / turnsRatios[turnsRatioAuxiliarySecondaryIndex] + magnetizationCurrent / 2;
            }

            double minimumPrimarySideTransformerCurrent = minimumPrimaryCurrent;
            double maximumPrimarySideTransformerCurrent = maximumPrimaryCurrent;
            double minimumPrimarySideTransformerVoltage = -inputVoltage;
            double maximumPrimarySideTransformerVoltage = inputVoltage;

            double maximumSecondarySideTransformerCurrentT1OfFET = maximumSecondaryCurrent;
            double minimumSecondarySideTransformerCurrentT2OtherFET = (minimumSecondaryCurrent / mainSecondaryTurnsRatio + magnetizationCurrent / 2) / 2 * mainSecondaryTurnsRatio - inductorCurrentRipple / 2;
            double maximumSecondarySideTransformerCurrentT2OtherFET = (minimumSecondaryCurrent / mainSecondaryTurnsRatio + magnetizationCurrent / 2) / 2 * mainSecondaryTurnsRatio;
            double minimumSecondarySideTransformerCurrentT2OfFET = 0;
            double maximumSecondarySideTransformerCurrentT2OfFET = maximumSecondaryCurrent - maximumSecondarySideTransformerCurrentT2OtherFET;
            double minimumSecondarySideTransformerVoltage = -inputVoltage / mainSecondaryTurnsRatio;
            double maximumSecondarySideTransformerVoltage = inputVoltage / mainSecondaryTurnsRatio;

            double minimumPrimarySideTransformerVoltageT3 = -(mainOutputVoltage + diodeVoltageDrop) * mainSecondaryTurnsRatio;
            double maximumPrimarySideTransformerVoltageT3 = (mainOutputVoltage + diodeVoltageDrop) * mainSecondaryTurnsRatio;
            double minimumSecondarySideTransformerCurrentT3 = 0;
            double maximumSecondarySideTransformerCurrentT3 = maximumSecondarySideTransformerCurrentT2OtherFET - maximumSecondarySideTransformerCurrentT2OfFET;
            double minimumSecondarySideTransformerVoltageT3 = -mainOutputVoltage - diodeVoltageDrop;
            double maximumSecondarySideTransformerVoltageT3 = mainOutputVoltage + diodeVoltageDrop;

            // First Primary
            {
                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        minimumPrimarySideTransformerCurrent,
                        maximumPrimarySideTransformerCurrent,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        maximumPrimarySideTransformerVoltage,
                        maximumPrimarySideTransformerVoltage,
                        0,
                        0,
                        minimumPrimarySideTransformerVoltageT3,
                        minimumPrimarySideTransformerVoltageT3,
                        minimumPrimarySideTransformerVoltage,
                        minimumPrimarySideTransformerVoltage,
                        0,
                        0,
                        maximumPrimarySideTransformerVoltageT3,
                        maximumPrimarySideTransformerVoltageT3
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        t1 + t2,
                        t1 + t2,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period / 2 + t1 + t2,
                        period / 2 + t1 + t2,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "First primary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }

            // Second Primary
            {
                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        0,
                        minimumPrimarySideTransformerCurrent,
                        maximumPrimarySideTransformerCurrent,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        minimumPrimarySideTransformerVoltage,
                        minimumPrimarySideTransformerVoltage,
                        0,
                        0,
                        maximumPrimarySideTransformerVoltageT3,
                        maximumPrimarySideTransformerVoltageT3,
                        maximumPrimarySideTransformerVoltage,
                        maximumPrimarySideTransformerVoltage,
                        0,
                        0,
                        minimumPrimarySideTransformerVoltageT3,
                        minimumPrimarySideTransformerVoltageT3
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        t1 + t2,
                        t1 + t2,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period / 2 + t1 + t2,
                        period / 2 + t1 + t2,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Second primary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }

            // First Main Secondary
            {
                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        0,
                        0,
                        maximumSecondarySideTransformerCurrentT2OtherFET,
                        minimumSecondarySideTransformerCurrentT2OtherFET,
                        maximumSecondarySideTransformerCurrentT3,
                        minimumSecondarySideTransformerCurrentT3,
                        maximumSecondarySideTransformerCurrentT1OfFET,
                        maximumSecondarySideTransformerCurrentT2OfFET,
                        minimumSecondarySideTransformerCurrentT2OfFET,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        t1 + t2,
                        t1 + t2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period / 2 + t1 + t2,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        minimumSecondarySideTransformerVoltage,
                        minimumSecondarySideTransformerVoltage,
                        0,
                        0,
                        maximumSecondarySideTransformerVoltageT3,
                        maximumSecondarySideTransformerVoltageT3,
                        maximumSecondarySideTransformerVoltage,
                        maximumSecondarySideTransformerVoltage,
                        0,
                        0,
                        minimumSecondarySideTransformerVoltageT3,
                        minimumSecondarySideTransformerVoltageT3
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        t1 + t2,
                        t1 + t2,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period / 2 + t1 + t2,
                        period / 2 + t1 + t2,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "First secondary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }

            // Second Main Secondary
            {
                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        0,
                        maximumSecondarySideTransformerCurrentT1OfFET,
                        maximumSecondarySideTransformerCurrentT2OfFET,
                        minimumSecondarySideTransformerCurrentT2OfFET,
                        0,
                        maximumSecondarySideTransformerCurrentT2OtherFET,
                        minimumSecondarySideTransformerCurrentT2OtherFET,
                        maximumSecondarySideTransformerCurrentT3,
                        minimumSecondarySideTransformerCurrentT3
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        t1 + t2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period / 2 + t1 + t2,
                        period / 2 + t1 + t2,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        maximumSecondarySideTransformerVoltage,
                        maximumSecondarySideTransformerVoltage,
                        0,
                        0,
                        minimumSecondarySideTransformerVoltageT3,
                        minimumSecondarySideTransformerVoltageT3,
                        minimumSecondarySideTransformerVoltage,
                        minimumSecondarySideTransformerVoltage,
                        0,
                        0,
                        maximumSecondarySideTransformerVoltageT3,
                        maximumSecondarySideTransformerVoltageT3
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        t1 + t2,
                        t1 + t2,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period / 2 + t1 + t2,
                        period / 2 + t1 + t2,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Second secondary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }

            for (size_t auxiliarySecondaryIndex = 1; auxiliarySecondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++auxiliarySecondaryIndex) {

                auto auxiliaryInductorCurrentRipple = get_current_ripple_ratio() * outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex];
                auto auxiliaryOutputVoltage = outputOperatingPoint.get_output_voltages()[auxiliarySecondaryIndex];
                auto minimumAuxiliarySecondaryCurrent = outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex] - auxiliaryInductorCurrentRipple / 2;
                auto maximumAuxiliarySecondaryCurrent = outputOperatingPoint.get_output_currents()[auxiliarySecondaryIndex] + auxiliaryInductorCurrentRipple / 2;
                size_t turnsRatioAuxiliarySecondaryIndex = 2 + auxiliarySecondaryIndex;
                double turnsRatioAuxiliarySecondary = turnsRatios[turnsRatioAuxiliarySecondaryIndex];
                minimumPrimaryCurrent += minimumAuxiliarySecondaryCurrent / turnsRatioAuxiliarySecondary;
                maximumPrimaryCurrent += maximumAuxiliarySecondaryCurrent / turnsRatioAuxiliarySecondary;

                double maximumAuxiliarySecondarySideTransformerCurrentT1OfFET = maximumAuxiliarySecondaryCurrent;
                double minimumAuxiliarySecondarySideTransformerCurrentT2OtherFET = (minimumAuxiliarySecondaryCurrent / turnsRatioAuxiliarySecondary + magnetizationCurrent / 2) / 2 * turnsRatioAuxiliarySecondary - inductorCurrentRipple / 2;
                double maximumAuxiliarySecondarySideTransformerCurrentT2OtherFET = (minimumAuxiliarySecondaryCurrent / turnsRatioAuxiliarySecondary + magnetizationCurrent / 2) / 2 * turnsRatioAuxiliarySecondary;
                double minimumAuxiliarySecondarySideTransformerCurrentT2OfFET = 0;
                double maximumAuxiliarySecondarySideTransformerCurrentT2OfFET = maximumAuxiliarySecondaryCurrent - maximumAuxiliarySecondarySideTransformerCurrentT2OtherFET;
                double minimumAuxiliarySecondarySideTransformerVoltage = -inputVoltage / turnsRatioAuxiliarySecondary;
                double maximumAuxiliarySecondarySideTransformerVoltage = inputVoltage / turnsRatioAuxiliarySecondary;

                double minimumAuxiliarySecondarySideTransformerCurrentT3 = 0;
                double maximumAuxiliarySecondarySideTransformerCurrentT3 = maximumAuxiliarySecondarySideTransformerCurrentT2OtherFET - maximumAuxiliarySecondarySideTransformerCurrentT2OfFET;
                double minimumAuxiliarySecondarySideTransformerVoltageT3 = -auxiliaryOutputVoltage - diodeVoltageDrop;
                double maximumAuxiliarySecondarySideTransformerVoltageT3 = auxiliaryOutputVoltage + diodeVoltageDrop;

                Waveform currentWaveform;
                Waveform voltageWaveform;
                // Current
                {
                    std::vector<double> data = {
                        0,
                        0,
                        maximumAuxiliarySecondarySideTransformerCurrentT2OtherFET,
                        minimumAuxiliarySecondarySideTransformerCurrentT2OtherFET,
                        maximumAuxiliarySecondarySideTransformerCurrentT3,
                        minimumAuxiliarySecondarySideTransformerCurrentT3,
                        maximumAuxiliarySecondarySideTransformerCurrentT1OfFET,
                        maximumAuxiliarySecondarySideTransformerCurrentT2OfFET,
                        minimumAuxiliarySecondarySideTransformerCurrentT2OfFET,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        t1 + t2,
                        t1 + t2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period / 2 + t1 + t2,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                // Voltage
                {
                    std::vector<double> data = {
                        minimumAuxiliarySecondarySideTransformerVoltage,
                        minimumAuxiliarySecondarySideTransformerVoltage,
                        0,
                        0,
                        maximumAuxiliarySecondarySideTransformerVoltageT3,
                        maximumAuxiliarySecondarySideTransformerVoltageT3,
                        maximumAuxiliarySecondarySideTransformerVoltage,
                        maximumAuxiliarySecondarySideTransformerVoltage,
                        0,
                        0,
                        minimumAuxiliarySecondarySideTransformerVoltageT3,
                        minimumAuxiliarySecondarySideTransformerVoltageT3
                    };
                    std::vector<double> time = {
                        0,
                        t1,
                        t1,
                        t1 + t2,
                        t1 + t2,
                        period / 2,
                        period / 2,
                        period / 2 + t1,
                        period / 2 + t1,
                        period / 2 + t1 + t2,
                        period / 2 + t1 + t2,
                        period
                    };
                    voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    voltageWaveform.set_data(data);
                    voltageWaveform.set_time(time);
                }

                auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Auxiliary " + std::to_string(auxiliarySecondaryIndex));
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);

            }

        }
        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    bool PushPull::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        for (size_t pushPullOperatingPointIndex = 1; pushPullOperatingPointIndex < get_operating_points().size(); ++pushPullOperatingPointIndex) {
            if (get_operating_points()[pushPullOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[pushPullOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
                if (!assert) {
                    return false;
                }
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Different operating points cannot have different number of output currents");
            }
        }
        if (!get_input_voltage().get_nominal() && !get_input_voltage().get_maximum() && !get_input_voltage().get_minimum()) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "No input voltage introduced");
        }

        return true;
    }

    DesignRequirements PushPull::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double diodeVoltageDrop = get_diode_voltage_drop();
        double dutyCycle = get_maximum_duty_cycle();

        double efficiency = 1;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }

        // Turns ratio calculation
        std::vector<double> turnsRatios;
        turnsRatios.push_back(1);  // Second primary
        turnsRatios.push_back(0);  // First secondary
        turnsRatios.push_back(0);  // Second secondary

        // Remaining auxiliary secondaries
        for (size_t secondaryIndex = 0; secondaryIndex < get_operating_points()[0].get_output_voltages().size() - 1; ++secondaryIndex) {
            turnsRatios.push_back(0);
        }

        for (auto pushPullOperatingPoint : get_operating_points()) {
            double mainSecondaryVoltage = pushPullOperatingPoint.get_output_voltages()[0];
            double mainSecondaryTurnsRatio = dutyCycle * 2 * minimumInputVoltage / (mainSecondaryVoltage + diodeVoltageDrop);
            turnsRatios[1] = std::max(turnsRatios[1], mainSecondaryTurnsRatio);
            turnsRatios[2] = std::max(turnsRatios[2], mainSecondaryTurnsRatio);
            for (size_t secondaryIndex = 0; secondaryIndex < pushPullOperatingPoint.get_output_voltages().size() - 1; ++secondaryIndex) {
                double auxiliarySecondaryVoltage = pushPullOperatingPoint.get_output_voltages()[secondaryIndex + 1];
                auto turnsRatio = dutyCycle * 2 * minimumInputVoltage / (auxiliarySecondaryVoltage + diodeVoltageDrop);
                turnsRatios[secondaryIndex + 3] = std::max(turnsRatios[secondaryIndex + 3], turnsRatio);
            }
        }

        // Inductance calculation
        double minimumNeededInductance = 0;
        for (auto pushPullOperatingPoint : get_operating_points()) {
            double switchingFrequency = pushPullOperatingPoint.get_switching_frequency();
            double totalOutputPower = 0;
            for (size_t secondaryIndex = 0; secondaryIndex < pushPullOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                double secondaryVoltage = pushPullOperatingPoint.get_output_voltages()[secondaryIndex];
                double secondaryCurrent = pushPullOperatingPoint.get_output_currents()[secondaryIndex];
                totalOutputPower += (secondaryVoltage + diodeVoltageDrop) * secondaryCurrent;
            }
            double period = 1.0 / switchingFrequency;
            double tOn = period * dutyCycle;
            double primaryCurrent = totalOutputPower / minimumInputVoltage / efficiency;
            double neededInductance = minimumInputVoltage * tOn / primaryCurrent;
            minimumNeededInductance = std::max(minimumNeededInductance, neededInductance);
        }

        if (get_maximum_switch_current()) {
            double maximumDrainSourceVoltage = 0;
            if (get_maximum_drain_source_voltage()) {
                maximumDrainSourceVoltage = get_maximum_drain_source_voltage().value();
            }
            double maximumSwitchCurrent = get_maximum_switch_current().value();
            // According to https://www.analog.com/cn/resources/technical-articles/high-frequency-push-pull-dc-dc-converter.html
            for (auto pushPullOperatingPoint : get_operating_points()) {
                double switchingFrequency = pushPullOperatingPoint.get_switching_frequency();
                double totalReflectedSecondaryCurrent = get_total_reflected_secondary_current(pushPullOperatingPoint, turnsRatios);

                double minimumInductance = 1.0 / switchingFrequency / 4 * (maximumInputVoltage - maximumDrainSourceVoltage) / (maximumSwitchCurrent - totalReflectedSecondaryCurrent);
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
        isolationSides.push_back(get_isolation_side_from_index(0));
        isolationSides.push_back(get_isolation_side_from_index(0));
        isolationSides.push_back(get_isolation_side_from_index(1));
        isolationSides.push_back(get_isolation_side_from_index(1));
        for (size_t windingIndex = 4; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex - 2));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::PUSH_PULL_CONVERTER);
        return designRequirements;
    }

    double PushPull::get_output_inductance(double mainSecondaryTurnsRatio) {
        double minimumOutputInductance = 0;
        auto dutyCycle = get_maximum_duty_cycle();
        
        for (auto outputOperatingPoint : get_operating_points()) {
            double mainOutputVoltage = outputOperatingPoint.get_output_voltages()[0];
            double mainOutputCurrent = outputOperatingPoint.get_output_currents()[0];
            double switchingFrequency = outputOperatingPoint.get_switching_frequency();
            
            // In Push-Pull, output sees pulsating DC at 2x switching frequency
            // During the off-time (1 - 2*D) portion, inductor discharges
            // L = V * t / ΔI
            // Where ΔI = ripple_ratio * I_out
            double offTime = (1.0 - 2.0 * dutyCycle) / (2.0 * switchingFrequency);  // Half period off-time
            double deltaI = get_current_ripple_ratio() * mainOutputCurrent;
            
            // Prevent division by zero if duty cycle is exactly 0.5
            if (offTime > 0 && deltaI > 0) {
                auto outputInductance = mainOutputVoltage * offTime / deltaI;
                minimumOutputInductance = std::max(minimumOutputInductance, outputInductance);
            }
        }
        
        // If calculation gives zero (e.g., D=0.5), use a reasonable default
        if (minimumOutputInductance < 1e-9) {
            minimumOutputInductance = 10e-6;  // 10 µH default
        }
        
        return minimumOutputInductance;
    }

    std::vector<OperatingPoint> PushPull::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        auto minimumOutputInductance = get_output_inductance(turnsRatios[1]);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t pushPullOperatingPointIndex = 0; pushPullOperatingPointIndex < get_operating_points().size(); ++pushPullOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[pushPullOperatingPointIndex], turnsRatios, magnetizingInductance, minimumOutputInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(pushPullOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    Inputs PushPull::process() {
        PushPull::run_checks(_assertErrors);

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

    std::vector<OperatingPoint> PushPull::process_operating_points(Magnetic magnetic) {
        PushPull::run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedPushPull::process() {
        PushPull::run_checks(_assertErrors);

        Inputs inputs;

        double minimumNeededInductance = get_desired_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();
        double minimumOutputInductance = 0;
        if (get_desired_output_inductance()) {
            minimumOutputInductance = get_desired_output_inductance().value();
        }
        else {
            minimumOutputInductance = get_output_inductance(turnsRatios[0]);
        }

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;


        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        DesignRequirements designRequirements;

        auto convertedTurnsRatios = convert_turns_ratios(turnsRatios);
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
        isolationSides.push_back(get_isolation_side_from_index(0));
        isolationSides.push_back(get_isolation_side_from_index(0));
        isolationSides.push_back(get_isolation_side_from_index(1));
        isolationSides.push_back(get_isolation_side_from_index(1));
        for (size_t windingIndex = 4; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex - 2));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::PUSH_PULL_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t pushPullOperatingPointIndex = 0; pushPullOperatingPointIndex < get_operating_points().size(); ++pushPullOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[pushPullOperatingPointIndex], convertedTurnsRatios, minimumNeededInductance, minimumOutputInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(pushPullOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

    std::string PushPull::generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames_;
    collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames_);
        
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
        double outputVoltage = opPoint.get_output_voltages()[0];
        double outputCurrent = opPoint.get_output_currents()[0];
        
        // turnsRatios[0] is second primary (always 1)
        // turnsRatios[1] and turnsRatios[2] are the main secondary halves
        // In push-pull, we use the secondary turns ratio for voltage transformation
        double mainTurnsRatio = turnsRatios[1];
        
        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period / 2 * dutyCycle;  // Each switch conducts for half period
        
        // Simulation: run steady-state periods for settling, then extract the last N periods
        int periodsToExtract = get_num_periods_to_extract();
        int numSteadyStatePeriods = get_num_steady_state_periods();
        const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;  // Steady state + extraction
        double simTime = numPeriodsTotal * period;
        double startTime = numSteadyStatePeriods * period;  // Start extracting after steady state
        double stepTime = period / 400;  // 400 points per period for smooth waveforms
        
        circuit << "* Push-Pull Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vout=" << outputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, N=" << mainTurnsRatio << "\n\n";
        
        // DC Input (center tap)
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // PWM Switches - alternating conduction
        circuit << "* PWM Switches (alternating)\n";
        circuit << "Vpwm1 pwm_ctrl1 0 PULSE(0 5 0 10n 10n " << tOn << " " << period << ")\n";
        circuit << "Vpwm2 pwm_ctrl2 0 PULSE(0 5 " << (period/2) << " 10n 10n " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=2.5 VH=0.5 RON=0.01\n\n";
        
        // Push-Pull Transformer - Two independent forward converter halves
        // Model as two separate forward converters operating alternately
        // CRITICAL: In Push-Pull, when one primary half is energized, corresponding secondary produces positive voltage
        circuit << "* Push-Pull as Two Independent Forward Converters\n";
        
        // Top half - operates when S1 is ON
        circuit << "* Top Half\n";
        circuit << "S1 vin_dc sw1_node pwm_ctrl1 0 SW1\n";
        circuit << "Vpri_top_sense sw1_node pri_top 0\n";
        circuit << "Lpri_top pri_top 0 " << std::scientific << magnetizingInductance << std::fixed << "\n";
        
        double secondaryInductance = magnetizingInductance / (mainTurnsRatio * mainTurnsRatio);
        circuit << "Lsec_top 0 sec_top " << std::scientific << secondaryInductance << std::fixed << "\n";  // Reversed: ground to sec_top
        circuit << "K_top Lpri_top Lsec_top 0.99\n";
        
        // Snubber for top half (RCD snubber to clamp voltage spikes)
        circuit << "Rsnub_top pri_top snub_top_c 10\n";
        circuit << "Csnub_top snub_top_c 0 220p\n";
        circuit << "Dsnub_top snub_top_c vin_dc DSNUB\n\n";
        
        // Bottom half - operates when S2 is ON
        circuit << "* Bottom Half\n";
        circuit << "S2 vin_dc sw2_node pwm_ctrl2 0 SW1\n";
        circuit << "Vpri_bot_sense sw2_node pri_bot 0\n";
        circuit << "Lpri_bot pri_bot 0 " << std::scientific << magnetizingInductance << std::fixed << "\n";
        circuit << "Lsec_bot 0 sec_bot " << std::scientific << secondaryInductance << std::fixed << "\n";  // Reversed: ground to sec_bot
        circuit << "K_bot Lpri_bot Lsec_bot 0.99\n";
        
        // Snubber for bottom half (RCD snubber to clamp voltage spikes)
        circuit << "Rsnub_bot pri_bot snub_bot_c 10\n";
        circuit << "Csnub_bot snub_bot_c 0 220p\n";
        circuit << "Dsnub_bot snub_bot_c vin_dc DSNUB\n\n";
        
        // Output rectifiers and filter  
        circuit << "* Output Rectifiers and Filter\n";
        circuit << ".model DIDEAL D(IS=1e-14 RS=0.01 CJO=1e-12)\n";
        circuit << ".model DSNUB D(IS=1e-14 RS=0.1 CJO=1e-11)\n";
        
        // Output rectifiers
        circuit << "Dsec_top sec_top sec_rect DIDEAL\n";
        circuit << "Dsec_bot sec_bot sec_rect DIDEAL\n";
        
        // Secondary current sense
        circuit << "Vsec_sense sec_rect sec_l_in 0\n";
        
        // Output inductor and capacitor
        double outputInductance = get_output_inductance(mainTurnsRatio);
        double loadResistance = outputVoltage / outputCurrent;
        circuit << "Lout sec_l_in vout " << std::scientific << outputInductance << std::fixed << "\n";
        circuit << "Cout vout 0 100u IC=" << outputVoltage << "\n";
        circuit << "Rload vout 0 " << loadResistance << "\n\n";
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << " UIC\n\n";  // Add UIC flag
        
        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(pri_top) v(pri_bot) i(Vpri_top_sense) i(Vpri_bot_sense)";
        circuit << " v(sec_top) v(sec_bot) i(Vsec_sense) v(vout)\n\n";
        
        // Options - tighter tolerances for cleaner waveforms
        circuit << ".options RELTOL=0.001 ABSTOL=1e-9 VNTOL=1e-6 TRTOL=7 ITL1=500 ITL4=100\n";
        circuit << ".ic v(vout)=" << outputVoltage << "\n\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> PushPull::simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto ppOpPoint = get_operating_points()[opIndex];
                
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = ppOpPoint.get_switching_frequency();
                
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = 1;
                config.keepTempFiles = false;
                
                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Define waveform name mapping for push-pull (4 windings)
                NgspiceRunner::WaveformNameMapping waveformMapping;
                
                // First primary (top)
                waveformMapping.push_back({{"voltage", "pri_top"}, {"current", "vpri_top_sense#branch"}});
                
                // Second primary (bottom)
                waveformMapping.push_back({{"voltage", "pri_bot"}, {"current", "vpri_bot_sense#branch"}});
                
                // First secondary (top)
                waveformMapping.push_back({{"voltage", "sec_top"}, {"current", "vsec_sense#branch"}});
                
                // Second secondary (bottom)
                waveformMapping.push_back({{"voltage", "sec_bot"}, {"current", "vsec_sense#branch"}});
                
                std::vector<std::string> windingNames = {"First primary", "Second primary", "First secondary", "Second secondary"};
                std::vector<bool> flipCurrentSign(4, false);
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    ppOpPoint.get_ambient_temperature(),
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
    
    std::vector<ConverterWaveforms> PushPull::simulate_and_extract_topology_waveforms(const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
    
    std::vector<ConverterWaveforms> results;
    
    NgspiceRunner runner;
    if (!runner.is_available()) {
        throw std::runtime_error("ngspice is not available for simulation");
    }
    
    std::vector<double> inputVoltages;
    std::vector<std::string> inputVoltagesNames;
    collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
    
    for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
        for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
            auto opPoint = get_operating_points()[opIndex];
            
            std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
            double switchingFrequency = opPoint.get_switching_frequency();
            
            SimulationConfig config;
            config.frequency = switchingFrequency;
            config.extractOnePeriod = true;
            config.numberOfPeriods = 2;
            config.keepTempFiles = false;
            
            auto simResult = runner.run_simulation(netlist, config);
            if (!simResult.success) {
                throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
            }
            
            std::map<std::string, size_t> nameToIndex;
            for (size_t i = 0; i < simResult.waveformNames.size(); ++i) {
                std::string lower = simResult.waveformNames[i];
                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                nameToIndex[lower] = i;
            }
            auto getWaveform = [&](const std::string& name) -> Waveform {
                std::string lower = name;
                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                auto it = nameToIndex.find(lower);
                if (it != nameToIndex.end()) return simResult.waveforms[it->second];
                return Waveform();
            };
            
            ConverterWaveforms wf;
            wf.set_switching_frequency(switchingFrequency);
            std::string name = inputVoltagesNames[inputVoltageIndex] + " input";
            if (get_operating_points().size() > 1) {
                name += " op. point " + std::to_string(opIndex);
            }
            wf.set_operating_point_name(name);
            
            wf.set_input_voltage(getWaveform("pri_top"));
            wf.set_input_current(getWaveform("vpri_top_sense#branch"));
            
            wf.get_mutable_output_voltages().push_back(getWaveform("vout"));
            wf.get_mutable_output_currents().push_back(getWaveform("vsec_sense#branch"));
            
            results.push_back(wf);
        }
    }
    
    return results;
}
} // namespace OpenMagnetics
