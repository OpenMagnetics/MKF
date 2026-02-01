#include "converter_models/PushPull.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include <sstream>
#include <algorithm>
#include <map>
#include "support/Exceptions.h"
#include "processors/NgspiceRunner.h"

namespace OpenMagnetics {

    double get_total_reflected_secondary_current(PushPullOperatingPoint pushPullOperatingPoint, std::vector<double> turnsRatios) {
        double totalReflectedSecondaryCurrent = 0;
        // Main secondary
        totalReflectedSecondaryCurrent += pushPullOperatingPoint.get_output_currents()[0] / turnsRatios[1];
        for (size_t secondaryIndex = 0; secondaryIndex < pushPullOperatingPoint.get_output_voltages().size() - 1; ++secondaryIndex) {
            totalReflectedSecondaryCurrent += pushPullOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex + 3];
        }

        return totalReflectedSecondaryCurrent;
    }

    std::vector<double> convert_turns_ratios(std::vector<double> turnsRatios) {
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

    OperatingPoint PushPull::process_operating_points_for_input_voltage(double inputVoltage, PushPullOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, double outputInductance) {

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
                Processed currentProcessed;
                Processed voltageProcessed;

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
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        for (auto outputOperatingPoint : get_operating_points()) {
            double mainOutputVoltage = outputOperatingPoint.get_output_voltages()[0];
            double switchingFrequency = outputOperatingPoint.get_switching_frequency();
            auto tOn = dutyCycle / switchingFrequency;
            auto outputInductance = (maximumInputVoltage / mainSecondaryTurnsRatio - get_diode_voltage_drop() - mainOutputVoltage) * tOn / get_current_ripple_ratio();
            minimumOutputInductance = std::max(minimumOutputInductance, outputInductance);
        }
        return minimumOutputInductance;
    }

    std::vector<OperatingPoint> PushPull::process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) {
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

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
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
        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
        }
        
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
        
        // turnsRatios[0] is the main secondary turns ratio
        double secondaryTurnsRatio = turnsRatios[0];
        double outputInductance = get_output_inductance(secondaryTurnsRatio);
        
        double outputVoltage = opPoint.get_output_voltages()[0];
        double outputCurrent = opPoint.get_output_currents()[0];
        double loadResistance = outputVoltage / outputCurrent;
        
        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double halfPeriod = period / 2.0;
        double tOn = halfPeriod * dutyCycle;
        
        // Simulation: run 10x the extraction periods for settling
        int periodsToExtract = get_num_periods_to_extract();
        const int numPeriodsTotal = 10 * periodsToExtract;
        double simTime = numPeriodsTotal * period;
        double startTime = (numPeriodsTotal - periodsToExtract) * period;
        double stepTime = period / 200;
        
        circuit << "* Push-Pull Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, N=" << secondaryTurnsRatio << "\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // PWM Switches - two switches alternating at half-period
        circuit << "* PWM Switches - alternating\n";
        circuit << "Vpwm1 pwm1_ctrl 0 PULSE(0 5 0 10n 10n " << tOn << " " << period << ")\n";
        circuit << "Vpwm2 pwm2_ctrl 0 PULSE(0 5 " << halfPeriod << " 10n 10n " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=2.5 VH=0.5\n";
        circuit << "S1 vin_dc pri1_p pwm1_ctrl 0 SW1\n";
        circuit << "S2 vin_dc pri2_p pwm2_ctrl 0 SW1\n\n";
        
        // Primary current sense
        circuit << "* Primary current sense\n";
        circuit << "Vpri1_sense pri1_p pri1_in 0\n";
        circuit << "Vpri2_sense pri2_p pri2_in 0\n\n";
        
        // Push-Pull Transformer (center-tapped primary and secondary)
        // Primary 1 and Primary 2 have same inductance, opposite polarity
        // Secondary 1 and Secondary 2 have same inductance, opposite polarity
        circuit << "* Push-Pull Transformer\n";
        circuit << "Lpri1 pri1_in ct_pri " << std::scientific << magnetizingInductance << std::fixed << "\n";
        circuit << "Lpri2 pri2_in ct_pri " << std::scientific << magnetizingInductance << std::fixed << "\n";
        circuit << "Rct ct_pri 0 1m\n";  // Center tap connected to ground via small resistance
        
        double secondaryInductance = magnetizingInductance / (secondaryTurnsRatio * secondaryTurnsRatio);
        circuit << "Lsec1 sec1_in ct_sec " << std::scientific << secondaryInductance << std::fixed << "\n";
        circuit << "Lsec2 sec2_in ct_sec " << std::scientific << secondaryInductance << std::fixed << "\n";
        circuit << "\n";
        
        // Coupling - all windings coupled
        circuit << "* Transformer coupling\n";
        circuit << "K12 Lpri1 Lpri2 1\n";
        circuit << "Ks12 Lsec1 Lsec2 1\n";
        circuit << "K1s1 Lpri1 Lsec1 1\n";
        circuit << "K1s2 Lpri1 Lsec2 1\n";
        circuit << "K2s1 Lpri2 Lsec1 1\n";
        circuit << "K2s2 Lpri2 Lsec2 1\n\n";
        
        // Diode model
        circuit << "* Diode model\n";
        circuit << ".model DIDEAL D(IS=1e-14 RS=1e-6)\n\n";
        
        // Output stage - full-wave rectifier with LC filter
        circuit << "* Output stage - full-wave rectifier\n";
        circuit << "D1 sec1_in rect_out DIDEAL\n";
        circuit << "D2 sec2_in rect_out DIDEAL\n";
        circuit << "Vsec_sense ct_sec sec_sense 0\n";
        
        // LC filter
        circuit << "Lout rect_out L_out " << std::scientific << outputInductance << std::fixed << "\n";
        circuit << "VL_sense L_out vout 0\n";
        circuit << "Cout vout 0 100u IC=" << outputVoltage << "\n";
        circuit << "Rload vout 0 " << loadResistance << "\n\n";
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << "\n\n";
        
        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(pri1_in) v(pri2_in) v(sec1_in) v(sec2_in) v(vout)";
        circuit << " i(Vpri1_sense) i(Vpri2_sense) i(Vsec_sense) i(VL_sense)\n\n";
        
        // Options for convergence
        circuit << ".options RELTOL=0.001 ABSTOL=1e-9 VNTOL=1e-6 ITL1=1000 ITL4=1000\n";
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
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto pushPullOpPoint = get_operating_points()[opIndex];
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = pushPullOpPoint.get_switching_frequency();
                
                // Run simulation
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = 1;  // Only one period for operating points
                config.keepTempFiles = false;
                
                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Build waveform mapping for Push-Pull converter:
                // Primary 1, Primary 2, Secondary 1, Secondary 2
                NgspiceRunner::WaveformNameMapping waveformMapping;
                
                // Primary 1 winding
                waveformMapping.push_back({{"voltage", "pri1_in"}, {"current", "vpri1_sense#branch"}});
                
                // Primary 2 winding
                waveformMapping.push_back({{"voltage", "pri2_in"}, {"current", "vpri2_sense#branch"}});
                
                // Secondary 1 winding
                waveformMapping.push_back({{"voltage", "sec1_in"}, {"current", "vl_sense#branch"}});
                
                // Secondary 2 winding
                waveformMapping.push_back({{"voltage", "sec2_in"}, {"current", "vl_sense#branch"}});
                
                std::vector<std::string> windingNames = {"Primary 1", "Primary 2", "Secondary 1", "Secondary 2"};
                std::vector<bool> flipCurrentSign = {false, false, false, false};
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    pushPullOpPoint.get_ambient_temperature(),
                    flipCurrentSign);
                
                // Set name
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

    std::vector<PushPullTopologyWaveforms> PushPull::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
        
        std::vector<PushPullTopologyWaveforms> topologyWaveforms;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Collect input voltages to simulate
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
            double inputVoltage = inputVoltages[inputVoltageIndex];
            
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto pushPullOpPoint = get_operating_points()[opIndex];
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = pushPullOpPoint.get_switching_frequency();
                double dutyCycle = get_maximum_duty_cycle();
                
                // Run simulation
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = 2;  // Two periods for topology waveform visualization
                config.keepTempFiles = false;
                
                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Build name-to-index map for waveform lookup (case-insensitive)
                std::map<std::string, size_t> nameToIndex;
                for (size_t i = 0; i < simResult.waveformNames.size(); ++i) {
                    std::string lower = simResult.waveformNames[i];
                    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                    nameToIndex[lower] = i;
                }
                
                // Helper lambda to get waveform data by name
                auto getWaveformData = [&](const std::string& name) -> std::vector<double> {
                    std::string lower = name;
                    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                    auto it = nameToIndex.find(lower);
                    if (it != nameToIndex.end()) {
                        return simResult.waveforms[it->second].get_data();
                    }
                    return {};
                };
                
                // Extract topology waveforms
                PushPullTopologyWaveforms waveforms;
                waveforms.frequency = switchingFrequency;
                waveforms.inputVoltageValue = inputVoltage;
                waveforms.outputVoltageValue = pushPullOpPoint.get_output_voltages()[0];
                waveforms.dutyCycle = dutyCycle;
                
                // Set name
                waveforms.operatingPointName = inputVoltagesNames[inputVoltageIndex] + " input";
                if (get_operating_points().size() > 1) {
                    waveforms.operatingPointName += " op. point " + std::to_string(opIndex);
                }
                
                // Extract time vector
                waveforms.time = getWaveformData("time");
                
                // Extract primary winding signals
                waveforms.primary1Voltage = getWaveformData("pri1_in");
                waveforms.primary2Voltage = getWaveformData("pri2_in");
                waveforms.primary1Current = getWaveformData("vpri1_sense#branch");
                waveforms.primary2Current = getWaveformData("vpri2_sense#branch");
                
                // Extract secondary signals
                waveforms.secondaryVoltage = getWaveformData("sec1_in");  // Use sec1 as reference
                waveforms.secondaryCurrent = getWaveformData("vsec_sense#branch");
                waveforms.outputVoltage = getWaveformData("vout");
                waveforms.outputInductorCurrent = getWaveformData("vl_sense#branch");
                
                topologyWaveforms.push_back(waveforms);
            }
        }
        
        return topologyWaveforms;
    }

} // namespace OpenMagnetics
