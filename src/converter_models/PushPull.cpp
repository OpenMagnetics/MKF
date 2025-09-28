#include "converter_models/PushPull.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>

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
            throw std::runtime_error("T1 cannot be larger than period/2, wrong topology configuration");
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
                throw std::runtime_error("T1 + T2 cannot be larger than period/2, wrong topology configuration");
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
            throw std::runtime_error("At least one operating point is needed");
        }
        for (size_t pushPullOperatingPointIndex = 1; pushPullOperatingPointIndex < get_operating_points().size(); ++pushPullOperatingPointIndex) {
            if (get_operating_points()[pushPullOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw std::runtime_error("Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[pushPullOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
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
} // namespace OpenMagnetics
