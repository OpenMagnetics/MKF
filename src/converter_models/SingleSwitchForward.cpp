#include "converter_models/SingleSwitchForward.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/CircuitSimulatorInterface.h"
#include "support/Utils.h"
#include <cfloat>
#include <sstream>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double SingleSwitchForward::get_total_reflected_secondary_current(ForwardOperatingPoint forwardOperatingPoint, std::vector<double> turnsRatios, double rippleRatio) {
        double totalReflectedSecondaryCurrent = 0;

        if (turnsRatios.size() != forwardOperatingPoint.get_output_currents().size() + 1) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Turns ratios must have one more position than outputs for the demagnetization winding");
        }

        for (size_t secondaryIndex = 0; secondaryIndex < forwardOperatingPoint.get_output_currents().size(); ++secondaryIndex) {
            totalReflectedSecondaryCurrent += forwardOperatingPoint.get_output_currents()[secondaryIndex] / turnsRatios[secondaryIndex + 1] * rippleRatio;
        }

        return totalReflectedSecondaryCurrent;
    }

    double SingleSwitchForward::get_maximum_duty_cycle() {
        if (get_duty_cycle()) {
            return get_duty_cycle().value();
        }
        else {
            return 0.45;
        }
    }

    SingleSwitchForward::SingleSwitchForward(const json& j) {
        from_json(j, *this);
    }

    AdvancedSingleSwitchForward::AdvancedSingleSwitchForward(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint SingleSwitchForward::process_operating_points_for_input_voltage(double inputVoltage, ForwardOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, double mainOutputInductance) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double mainOutputCurrent = outputOperatingPoint.get_output_currents()[0];
        double mainOutputVoltage = outputOperatingPoint.get_output_voltages()[0];
        double mainSecondaryTurnsRatio = turnsRatios[1];
        double diodeVoltageDrop = get_diode_voltage_drop();

        // Assume CCM
        auto period = 1.0 / switchingFrequency;
        // Switch ON time from volt-second balance on the output inductor:
        //   Vout = Vin/n · D − Vd  ⇒  D = (Vout+Vd)/(Vin/n)
        // SSF has NO inherent /2 (unlike push-pull/half-bridge where D is
        // per-half-cycle).  Matches the SPICE netlist tOn computation
        // (generate_ngspice_circuit line ~519) so SPICE and analytical use
        // the same duty cycle and waveform timing.
        auto t1 = period * (mainOutputVoltage + diodeVoltageDrop) / (inputVoltage / mainSecondaryTurnsRatio);
        if (t1 > period / 2) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "T1 cannot be larger than period/2, wrong topology configuration");
        }

        auto magnetizationCurrent = inputVoltage * t1 / inductance;
        auto minimumPrimaryCurrent = -magnetizationCurrent / 2;
        auto maximumPrimaryCurrent = magnetizationCurrent / 2;

        std::vector<double> minimumSecondaryCurrents;
        std::vector<double> maximumSecondaryCurrents;

        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
            auto outputCurrentRipple = get_current_ripple_ratio() * outputOperatingPoint.get_output_currents()[secondaryIndex];
            auto minimumSecondaryCurrent = outputOperatingPoint.get_output_currents()[secondaryIndex] - outputCurrentRipple / 2;
            minimumSecondaryCurrents.push_back(minimumSecondaryCurrent);

            auto maximumSecondaryCurrent = outputOperatingPoint.get_output_currents()[secondaryIndex] + outputCurrentRipple / 2;
            maximumSecondaryCurrents.push_back(maximumSecondaryCurrent);
            size_t turnsRatioSecondaryIndex = 1 + secondaryIndex;  // To avoid demagnetization winding
            minimumPrimaryCurrent += minimumSecondaryCurrent / turnsRatios[turnsRatioSecondaryIndex];
            maximumPrimaryCurrent += maximumSecondaryCurrent / turnsRatios[turnsRatioSecondaryIndex];
        }


        if (minimumPrimaryCurrent < 0) { // DCM
            double sqrtArg = 2 * mainOutputCurrent * mainOutputInductance * (mainOutputVoltage + diodeVoltageDrop) / (switchingFrequency * (inputVoltage / mainSecondaryTurnsRatio - diodeVoltageDrop - mainOutputVoltage) * (inputVoltage / mainSecondaryTurnsRatio));
            if (sqrtArg < 0) {
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "SingleSwitchForward: negative value under sqrt in DCM t1 calculation");
            }
            t1 = sqrt(sqrtArg);
            if (t1 > period / 2) {
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "T1 cannot be larger than period/2, wrong topology configuration");
            }
            minimumPrimaryCurrent = 0;
            maximumPrimaryCurrent = magnetizationCurrent;

            for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                auto outputCurrentRipple = get_current_ripple_ratio() * outputOperatingPoint.get_output_currents()[secondaryIndex];
                minimumSecondaryCurrents[secondaryIndex] = 0;
                maximumSecondaryCurrents[secondaryIndex] = outputCurrentRipple;

                size_t turnsRatioSecondaryIndex = 1 + secondaryIndex;  // To avoid demagnetization winding
                minimumPrimaryCurrent += minimumSecondaryCurrents[secondaryIndex] / turnsRatios[turnsRatioSecondaryIndex];
                maximumPrimaryCurrent += maximumSecondaryCurrents[secondaryIndex] / turnsRatios[turnsRatioSecondaryIndex];
            }
        }

        double td = t1;  // Demagnetization time equals on-time for Nt=Np
        double deadTime = period - t1 - td;
        // Use actual duty cycle (t1/period) for secondary waveforms, not the fixed maximum
        double actualDutyCycle = t1 / period;
        // Primary — use custom waveform with explicit voltage levels (+Vin, -Vin, 0)
        {
            double primaryCurrentPeakToPeak = maximumPrimaryCurrent - minimumPrimaryCurrent;
            auto currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, primaryCurrentPeakToPeak, switchingFrequency, actualDutyCycle, minimumPrimaryCurrent, deadTime);
            Waveform voltageWaveform;
            voltageWaveform.set_data(std::vector<double>{0, inputVoltage, inputVoltage, -inputVoltage, -inputVoltage, 0, 0});
            voltageWaveform.set_time(std::vector<double>{0, 0, t1, t1, t1 + td, t1 + td, period});
            voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }
        // Demagnetization winding — same voltage shape but inverted polarity
        {
            auto currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME, magnetizationCurrent, switchingFrequency, actualDutyCycle, minimumPrimaryCurrent, deadTime);
            Waveform voltageWaveform;
            voltageWaveform.set_data(std::vector<double>{0, -inputVoltage, -inputVoltage, inputVoltage, inputVoltage, 0, 0});
            voltageWaveform.set_time(std::vector<double>{0, 0, t1, t1, t1 + td, t1 + td, period});
            voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Demagnetization winding");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {

            double secondaryCurrentPeakToPeak = maximumSecondaryCurrents[secondaryIndex] - minimumSecondaryCurrents[secondaryIndex];

            size_t turnsRatioSecondaryIndex = 1 + secondaryIndex;  // To avoid demagnetization winding
            double minimumSecondaryVoltage = -(inputVoltage + diodeVoltageDrop) / turnsRatios[turnsRatioSecondaryIndex];
            double maximumSecondaryVoltage = inputVoltage / turnsRatios[turnsRatioSecondaryIndex];
            double secondaryVoltagePeakToPeak = maximumSecondaryVoltage - minimumSecondaryVoltage;
            double secondaryVoltageOffset = maximumSecondaryVoltage + minimumSecondaryVoltage;

            auto currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, secondaryCurrentPeakToPeak, switchingFrequency, actualDutyCycle, minimumSecondaryCurrents[secondaryIndex], 0);
            auto voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, secondaryVoltagePeakToPeak, switchingFrequency, actualDutyCycle, secondaryVoltageOffset, deadTime);
            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Secondary " + std::to_string(secondaryIndex));
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    bool SingleSwitchForward::run_checks(bool assert) {
        return Topology::run_checks_common(this, assert);
    }

    DesignRequirements SingleSwitchForward::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double diodeVoltageDrop = get_diode_voltage_drop();
        double dutyCycle = get_maximum_duty_cycle();

        // Turns ratio calculation
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size() + 1, 0);

        // Demagnetization winding has the same number of turns as primary.
        turnsRatios[0] = 1;

        for (auto forwardOperatingPoint : get_operating_points()) {
            for (size_t secondaryIndex = 0; secondaryIndex < forwardOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                // Size turns ratio at Vin_MIN: D is highest at the low-line
                // corner. n = D_max * Vin_min / (Vout + Vd) so D(Vin_min) = D_max
                // and D(Vin > Vin_min) < D_max. The previous formula sized at
                // Vin_max, which understated n and produced infeasible designs
                // that violated maximumDutyCycle at low Vin (silently masked
                // by a downstream std::min clamp until that clamp was replaced
                // with a throw).
                double turnsRatio = minimumInputVoltage * dutyCycle / (forwardOperatingPoint.get_output_voltages()[secondaryIndex] + diodeVoltageDrop);
                turnsRatios[secondaryIndex + 1] = std::max(turnsRatios[secondaryIndex + 1], turnsRatio);
            }
        }

        // Inductance calculation
        double minimumNeededInductance = 0;
        for (auto forwardOperatingPoint : get_operating_points()) {
            double switchingFrequency = forwardOperatingPoint.get_switching_frequency();
            double totalReflectedSecondaryCurrent = get_total_reflected_secondary_current(forwardOperatingPoint, turnsRatios);

            double neededInductance = minimumInputVoltage / (switchingFrequency * totalReflectedSecondaryCurrent);
            minimumNeededInductance = std::max(minimumNeededInductance, neededInductance);
        }

        if (get_maximum_switch_current()) {
            double maximumSwitchCurrent = get_maximum_switch_current().value();
            // According to https://www.analog.com/cn/resources/technical-articles/high-frequency-forward-pull-dc-dc-converter.html
            for (auto forwardOperatingPoint : get_operating_points()) {
                double switchingFrequency = forwardOperatingPoint.get_switching_frequency();
                double totalReflectedSecondaryCurrent = get_total_reflected_secondary_current(forwardOperatingPoint, turnsRatios, 1 + get_current_ripple_ratio());

                double minimumInductance = maximumInputVoltage * dutyCycle / switchingFrequency / (maximumSwitchCurrent - totalReflectedSecondaryCurrent);
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
        designRequirements.set_isolation_sides(
            Topology::create_isolation_sides(get_operating_points()[0].get_output_currents().size(), true));
        designRequirements.set_topology(Topologies::SINGLE_SWITCH_FORWARD_CONVERTER);
        return designRequirements;
    }

    double SingleSwitchForward::get_output_inductance(double secondaryTurnsRatio, size_t outputIndex) {
        double minimumOutputInductance = 0;
        auto dutyCycle = get_maximum_duty_cycle();
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        for (auto outputOperatingPoint : get_operating_points()) {
            double outputVoltage = outputOperatingPoint.get_output_voltages()[outputIndex];
            double switchingFrequency = outputOperatingPoint.get_switching_frequency();
            auto tOn = dutyCycle / switchingFrequency;
            auto outputInductance = (maximumInputVoltage / secondaryTurnsRatio - get_diode_voltage_drop() - outputVoltage) * tOn / get_current_ripple_ratio();
            minimumOutputInductance = std::max(minimumOutputInductance, outputInductance);
        }
        return minimumOutputInductance;
    }

    std::vector<OperatingPoint> SingleSwitchForward::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        extraLoVoltageWaveforms.clear();
        extraLoCurrentWaveforms.clear();
        extraLoInductances.clear();

        std::vector<double> outputInductancePerSecondary;

        for (size_t forwardOperatingPointIndex = 0; forwardOperatingPointIndex < get_operating_points().size(); ++forwardOperatingPointIndex) {
            outputInductancePerSecondary.push_back(get_output_inductance(turnsRatios[forwardOperatingPointIndex + 1], forwardOperatingPointIndex));
        }

        // Populate extraLoInductances (one per secondary winding)
        size_t numSecondaries = turnsRatios.size() - 1;  // turnsRatios[0] is demagnetization
        for (size_t s = 0; s < numSecondaries; ++s) {
            double lo = get_output_inductance(turnsRatios[s + 1], s);
            extraLoInductances.push_back(lo);
        }

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];

            for (size_t forwardOperatingPointIndex = 0; forwardOperatingPointIndex < get_operating_points().size(); ++forwardOperatingPointIndex) {
                auto& fop = get_operating_points()[forwardOperatingPointIndex];
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, fop, turnsRatios, magnetizingInductance, outputInductancePerSecondary[0]);

                // Capture Lo waveforms for each secondary
                double switchingFrequency = fop.get_switching_frequency();
                double period = 1.0 / switchingFrequency;
                double diodeVoltageDrop = get_diode_voltage_drop();
                double mainSecondaryTurnsRatio = turnsRatios[1];
                double mainOutputVoltage = fop.get_output_voltages()[0];
                double actualDutyCycle = (mainOutputVoltage + diodeVoltageDrop) / (inputVoltage / mainSecondaryTurnsRatio);
                // No silent clamp: violating maximumDutyCycle means the design
                // cannot deliver Vout at this Vin. Throw with actionable info.
                // Tolerance allows for rounding in design_requirements (turns
                // ratio is rounded to 2 decimals) without firing on noise.
                {
                    const double dMax = get_maximum_duty_cycle();
                    constexpr double dutyTolerance = 0.01;  // 1 % over maxD
                    if (actualDutyCycle > dMax * (1.0 + dutyTolerance)) {
                        throw InvalidInputException(ErrorCode::INVALID_INPUT,
                            "SingleSwitchForward: required dutyCycle " + std::to_string(actualDutyCycle) +
                            " exceeds maximumDutyCycle " + std::to_string(dMax) +
                            " at Vin=" + std::to_string(inputVoltage) +
                            ". Increase turns ratio (lower mainSecondaryTurnsRatio), raise Vin_min, or relax maximumDutyCycle.");
                    }
                }
                double tOn = actualDutyCycle * period;
                double tOff = period - tOn;

                for (size_t s = 0; s < numSecondaries; ++s) {
                    double turnsRatioSec = turnsRatios[s + 1];
                    double outputVoltage = fop.get_output_voltages()[s];
                    double outputCurrent = fop.get_output_currents()[s];
                    double currentRipple = get_current_ripple_ratio() * outputCurrent;

                    // Lo current: triangular in CCM, DC output current with ripple
                    Waveform loCurrentWfm;
                    {
                        double iMin = outputCurrent - currentRipple / 2;
                        double iMax = outputCurrent + currentRipple / 2;
                        loCurrentWfm.set_data(std::vector<double>{iMin, iMax, iMax, iMin});
                        loCurrentWfm.set_time(std::vector<double>{0, tOn, tOn, period});
                        loCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
                    }

                    // Lo voltage: (Vsec - Vout) during on-time, (-Vout) during off-time
                    double vSecOnTime = inputVoltage / turnsRatioSec - diodeVoltageDrop - outputVoltage;
                    double vSecOffTime = -outputVoltage;
                    Waveform loVoltageWfm;
                    {
                        loVoltageWfm.set_data(std::vector<double>{vSecOnTime, vSecOnTime, vSecOffTime, vSecOffTime});
                        loVoltageWfm.set_time(std::vector<double>{0, tOn, tOn, period});
                        loVoltageWfm.set_ancillary_label(WaveformLabel::CUSTOM);
                    }

                    extraLoCurrentWaveforms.push_back(loCurrentWfm);
                    extraLoVoltageWaveforms.push_back(loVoltageWfm);
                }

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(forwardOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<std::variant<Inputs, CAS::Inputs>> SingleSwitchForward::get_extra_components_inputs(
        ExtraComponentsMode /*mode*/,
        std::optional<Magnetic> /*magnetic*/)
    {
        if (extraLoInductances.empty() || extraLoVoltageWaveforms.empty())
            throw std::runtime_error("SingleSwitchForward::get_extra_components_inputs: call process_operating_points() first");

        size_t numSecondaries = extraLoInductances.size();
        size_t nOps = extraLoVoltageWaveforms.size() / numSecondaries;

        std::vector<std::variant<Inputs, CAS::Inputs>> result;

        for (size_t s = 0; s < numSecondaries; ++s) {
            double lo = extraLoInductances[s];
            if (lo <= 0)
                throw std::runtime_error("SingleSwitchForward::get_extra_components_inputs: computedOutputInductance <= 0 for secondary " + std::to_string(s));

            Inputs masInputs;
            DesignRequirements dr;

            DimensionWithTolerance inductance;
            inductance.set_nominal(lo);
            dr.set_magnetizing_inductance(inductance);

            std::string inductorName = (numSecondaries == 1) ? "outputInductor" : ("outputInductor_" + std::to_string(s + 1));
            dr.set_name(inductorName);
            dr.set_topology(Topologies::SINGLE_SWITCH_FORWARD_CONVERTER);
            dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
            dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::SECONDARY});
            masInputs.set_design_requirements(dr);

            std::vector<OperatingPoint> masOps;
            for (size_t opIdx = 0; opIdx < nOps; ++opIdx) {
                size_t wfmIdx = opIdx * numSecondaries + s;
                OperatingPoint op;
                double fsw = get_operating_points()[opIdx % get_operating_points().size()].get_switching_frequency();
                auto excitation = complete_excitation(
                    extraLoCurrentWaveforms[wfmIdx],
                    extraLoVoltageWaveforms[wfmIdx],
                    fsw, "Primary");
                op.get_mutable_excitations_per_winding().push_back(excitation);
                masOps.push_back(op);
            }
            masInputs.set_operating_points(masOps);
            result.emplace_back(std::move(masInputs));
        }

        return result;
    }

    std::vector<OperatingPoint> SingleSwitchForward::process_operating_points(Magnetic magnetic) {
        SingleSwitchForward::run_checks(_assertErrors);

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedSingleSwitchForward::process() {
        SingleSwitchForward::run_checks(_assertErrors);

        Inputs inputs;

        double minimumNeededInductance = get_desired_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();
        std::vector<double> outputInductancePerSecondary;

        if (get_desired_output_inductances()) {
            outputInductancePerSecondary = get_desired_output_inductances().value();
        }
        else {
            // For SingleSwitchForward, turnsRatios[0] is demagnetization winding, so secondary indices start at 1
            for (size_t secondaryIndex = 1; secondaryIndex < turnsRatios.size(); ++secondaryIndex) {
                // secondaryIndex - 1 because output_voltages/output_currents don't include demag winding
                auto minimumOutputInductance = get_output_inductance(turnsRatios[secondaryIndex], secondaryIndex - 1);
                outputInductancePerSecondary.push_back(minimumOutputInductance);
            }
        }

        if (turnsRatios.size() != get_operating_points()[0].get_output_currents().size() + 1) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Turns ratios must have one more position than outputs for the demagnetization winding");
        }

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        inputs.set_design_requirements(process_design_requirements());

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t forwardOperatingPointIndex = 0; forwardOperatingPointIndex < get_operating_points().size(); ++forwardOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[forwardOperatingPointIndex], turnsRatios, minimumNeededInductance, outputInductancePerSecondary[0]);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(forwardOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

    DesignRequirements AdvancedSingleSwitchForward::process_design_requirements() {
        // Issue M1: build DR directly from desired* fields. Do NOT chain to
        // SingleSwitchForward::process_design_requirements() — the parent
        // sizes turns ratios from voltages and inductance from current
        // ripple, which would override the user-provided desired* values.
        std::vector<double> turnsRatios = get_desired_turns_ratios();
        double minimumNeededInductance = get_desired_inductance();

        if (turnsRatios.size() != get_operating_points()[0].get_output_currents().size() + 1) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Turns ratios must have one more position than outputs for the demagnetization winding");
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        for (auto turnsRatio : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(minimumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        designRequirements.set_isolation_sides(
            Topology::create_isolation_sides(get_operating_points()[0].get_output_currents().size(), true));
        designRequirements.set_topology(Topologies::SINGLE_SWITCH_FORWARD_CONVERTER);
        return designRequirements;
    }

    std::string SingleSwitchForward::generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {
        
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("operatingPointIndex out of range");
        }
        
        double inputVoltage = inputVoltages[inputVoltageIndex];
        auto opPoint = get_operating_points()[operatingPointIndex];
        
        double switchingFrequency = opPoint.get_switching_frequency();
        double diodeVoltageDrop = get_diode_voltage_drop();

        // Number of secondaries (turns ratios[0] is demagnetization winding)
        if (turnsRatios.empty()) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "SingleSwitchForward: turnsRatios must not be empty");
        }
        size_t numSecondaries = turnsRatios.size() - 1;

        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        // Compute actual duty cycle from voltage balance: D = (Vout + Vd) * n / Vin
        // turnsRatios[0] is demagnetization winding, turnsRatios[1] is first secondary
        double mainOutputVoltage = opPoint.get_output_voltages()[0];
        double mainSecondaryTurnsRatio = (numSecondaries > 0) ? turnsRatios[1] : 1.0;
        double dutyCycle = (mainOutputVoltage + diodeVoltageDrop) / (inputVoltage / mainSecondaryTurnsRatio);
        // No silent clamp: throw if the operating-point duty exceeds maxD,
        // with a 1 % tolerance for design-requirements rounding (n is
        // rounded to 2 decimals).
        {
            const double dMax = get_maximum_duty_cycle();
            constexpr double dutyTolerance = 0.01;
            if (dutyCycle > dMax * (1.0 + dutyTolerance)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "SingleSwitchForward::generate_ngspice_circuit: required dutyCycle " +
                    std::to_string(dutyCycle) + " exceeds maximumDutyCycle " + std::to_string(dMax) +
                    " at Vin=" + std::to_string(inputVoltage) +
                    ". Increase turns ratio, raise Vin_min, or relax maximumDutyCycle.");
            }
        }
        double tOn = period * dutyCycle;
        
        // Simulation: run steady-state periods for settling, then extract the last N periods
        int periodsToExtract = get_num_periods_to_extract();
        int numSteadyStatePeriods = get_num_steady_state_periods();
        const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;  // Steady state + extraction
        double simTime = numPeriodsTotal * period;
        double startTime = numSteadyStatePeriods * period;  // Start extracting after steady state
        double stepTime = period / 200;
        
        circuit << "* Single-Switch Forward Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, " << numSecondaries << " secondaries\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";

        // §8a.5 input-current sense: a zero-V source on the high-side
        // of S1's drain isolates the true switch-channel current from
        // the demag return current (which flows back into vin_dc via
        // Ddemag) and any future snubber currents. i(Vin) would
        // sum all of those, so it is NOT the converter's input draw.
        // Sign convention: + terminal is vin_dc, - is q1_drain, so
        // positive i(Vq1_sense) means current flowing OUT of vin_dc
        // INTO S1 (the converter's on-time draw).
        circuit << "* §8a.5 input-current sense\n";
        circuit << "Vq1_sense vin_dc q1_drain 0\n\n";

        // PWM Main Switch
        circuit << "* PWM Main Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 5 0 10n 10n " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=2.5 VH=0.5\n";
        circuit << "S1 q1_drain pri_p pwm_ctrl 0 SW1\n\n";
        
        // Primary current sense
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense pri_p pri_in 0\n\n";
        
        // Forward Transformer with demagnetization winding
        // Primary winding
        circuit << "* Forward Transformer\n";
        circuit << "Lpri pri_in 0 " << std::scientific << magnetizingInductance << std::fixed << "\n";
        
        // Demagnetization winding (same turns as primary for single-switch forward = 1:1)
        // CRITICAL: Terminals reversed to create opposite dot polarity for proper demagnetization
        // When switch turns off, primary voltage goes negative, demag voltage goes positive
        circuit << "Ldemag 0 demag_in " << std::scientific << magnetizingInductance << std::fixed << "\n";
        
        // Secondary windings
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double secondaryInductance = magnetizingInductance / (turnsRatios[secIdx + 1] * turnsRatios[secIdx + 1]);
            circuit << "Lsec" << secIdx << " sec" << secIdx << "_in 0 " << std::scientific << secondaryInductance << std::fixed << "\n";
        }
        
        // Couple all windings together using pairwise K statements (ngspice requires this)
        circuit << "* Coupling: All windings coupled pairwise\n";
        circuit << "Kpri_demag Lpri Ldemag 0.9999\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "Kpri_sec" << secIdx << " Lpri Lsec" << secIdx << " 0.9999\n";
            circuit << "Kdemag_sec" << secIdx << " Ldemag Lsec" << secIdx << " 0.9999\n";
        }
        // Couple secondaries to each other if there are multiple
        for (size_t i = 0; i < numSecondaries; ++i) {
            for (size_t j = i + 1; j < numSecondaries; ++j) {
                circuit << "Ksec" << i << "_sec" << j << " Lsec" << i << " Lsec" << j << " 0.9999\n";
            }
        }
        circuit << "\n";
        
        // Demagnetization diode - conducts when main switch is off to reset the core
        // The demag winding has opposite polarity through the coupling definition
        circuit << "* Demagnetization Diode and current sense\n";
        circuit << ".model DIDEAL D(IS=1e-14 RS=1e-6)\n";
        circuit << "Vdemag_sense demag_in demag_sense 0\n";
        circuit << "Ddemag demag_sense vin_dc DIDEAL\n\n";
        
        // Output stages for each secondary
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "* Secondary " << secIdx << " output stage\n";
            
            // Forward rectifier - conducts when main switch is on
            circuit << "Dfwd" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_rect DIDEAL\n";
            
            // Freewheeling diode - conducts when main switch is off
            circuit << "Dfw" << secIdx << " 0 sec" << secIdx << "_rect DIDEAL\n";
            
            // Snubber resistors for convergence
            circuit << "Rsnub_fwd" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_rect 1MEG\n";
            circuit << "Rsnub_fw" << secIdx << " 0 sec" << secIdx << "_rect 1MEG\n";
            
            // Output inductor (LC filter) with ESR for damping
            double outputVoltage = opPoint.get_output_voltages()[secIdx];
            double outputCurrent = opPoint.get_output_currents()[secIdx];
            double outputInductance = get_output_inductance(turnsRatios[secIdx + 1], secIdx);
            
            circuit << "Vsec_sense" << secIdx << " sec" << secIdx << "_rect sec" << secIdx << "_l_in 0\n";
            // Add series resistance to output inductor (10mΩ) for damping
            circuit << "Rlout" << secIdx << " sec" << secIdx << "_l_in lout" << secIdx << " 0.01\n";
            circuit << "Lout" << secIdx << " lout" << secIdx << " vout_c" << secIdx << " " << std::scientific << outputInductance << std::fixed << "\n";
            
            double loadResistance = outputVoltage / outputCurrent;
            // Add ESR to output capacitor (100mΩ) for damping
            circuit << "Rcout" << secIdx << " vout_c" << secIdx << " vout" << secIdx << " 0.1\n";
            circuit << "Cout" << secIdx << " vout" << secIdx << " 0 100u IC=" << outputVoltage << "\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << "\n\n";
        
        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(pri_in) i(Vpri_sense) v(demag_in) i(Vdemag_sense) v(vin_dc) i(Vin) i(Vq1_sense)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) i(Vsec_sense" << secIdx << ") v(vout" << secIdx << ")";
        }
        circuit << "\n\n";
        
        // Options
        circuit << ".options RELTOL=0.001 ABSTOL=1e-9 VNTOL=1e-6 ITL1=1000 ITL4=1000\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << ".ic v(vout" << secIdx << ")=" << opPoint.get_output_voltages()[secIdx] << "\n";
        }
        circuit << "\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> SingleSwitchForward::simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        if (turnsRatios.empty()) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "SingleSwitchForward: turnsRatios must not be empty");
        }
        size_t numSecondaries = turnsRatios.size() - 1;

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto forwardOpPoint = get_operating_points()[opIndex];
                
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = forwardOpPoint.get_switching_frequency();
                
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);

                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Define waveform name mapping
                NgspiceRunner::WaveformNameMapping waveformMapping;
                
                // Primary winding
                waveformMapping.push_back({{"voltage", "pri_in"}, {"current", "vpri_sense#branch"}});
                
                // Demagnetization winding
                waveformMapping.push_back({{"voltage", "demag_in"}, {"current", "vdemag_sense#branch"}});
                
                // Secondary windings
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    std::string voltageName = "sec" + std::to_string(secIdx) + "_in";
                    std::string currentName = "vsec_sense" + std::to_string(secIdx) + "#branch";
                    waveformMapping.push_back({{"voltage", voltageName}, {"current", currentName}});
                }
                
                std::vector<std::string> windingNames;
                windingNames.push_back("Primary");
                windingNames.push_back("Demagnetization winding");
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    windingNames.push_back("Secondary " + std::to_string(secIdx));
                }
                
                std::vector<bool> flipCurrentSign(2 + numSecondaries, false);
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    forwardOpPoint.get_ambient_temperature(),
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
    
    std::vector<ConverterWaveforms> SingleSwitchForward::simulate_and_extract_topology_waveforms(const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods) {
    
    std::vector<ConverterWaveforms> results;
    
    NgspiceRunner runner;
    if (!runner.is_available()) {
        throw std::runtime_error("ngspice is not available for simulation");
    }
    
    // Save original value and set the requested number of periods
    int originalNumPeriodsToExtract = get_num_periods_to_extract();
    set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
    
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
            config.numberOfPeriods = numberOfPeriods;
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
            
            // §5.1 / §8a.5 converter-port stream:
            //   input_voltage = v(vin_dc) (real DC rail, constant Vin).
            //   input_current = i(Vq1_sense) — zero-V sense source on
            //     S1's drain. Captures only the on-time switch-channel
            //     current. i(Vin) (the prior probe) sums Vq1_sense PLUS
            //     the demag-diode return current (Ddemag dumps energy
            //     back into vin_dc during off-time) PLUS any future
            //     snubber currents — not the converter's input draw.
            //     Sign: + terminal of Vq1_sense is vin_dc, so positive
            //     i(Vq1_sense) is current flowing OUT of vin_dc INTO
            //     the switch (non-negative in CCM/DCM).
            //   Bug C: the bare-node v(pri_in) probe used by the
            //     per-winding extractor is already a true winding-only
            //     voltage (Lpri is wired pri_in -> 0 — winding bottom
            //     is the netlist reference), so no Evpri_w E-source
            //     differential probe is needed. Likewise each secondary
            //     Lsec<i> is wired sec<i>_in -> 0.
            wf.set_input_voltage(getWaveform("vin_dc"));

            Waveform iIn = getWaveform("vq1_sense#branch");
            if (iIn.get_data().empty()) {
                throw std::runtime_error(
                    "SingleSwitchForward simulate_and_extract_topology_waveforms: "
                    "missing i(Vq1_sense) — netlist/save mismatch");
            }

            // ngspice SW1 model produces narrow di/dt spikes (~10^5 A)
            // around switching edges that are a numerical artefact of
            // the discontinuous RON/ROFF transition, not a physical
            // bound. Clamp samples to ±2·max(|i(Vpri_sense)|) — the
            // primary current sets the actual conduction-current
            // scale; this preserves all real samples while discarding
            // the artefact spikes that would poison the DC mean used
            // downstream by §5.1.
            //
            // numerical guard against ngspice SW-model di/dt spikes
            // (~10^5 A), not a physical bound.
            {
                auto& d = iIn.get_mutable_data();
                Waveform iPri = getWaveform("vpri_sense#branch");
                double iMax = 0.0;
                for (double v : iPri.get_data()) iMax = std::max(iMax, std::abs(v));
                if (iMax > 0.0) {
                    const double clamp = 2.0 * iMax;
                    for (auto& v : d) {
                        if (v >  clamp) v =  clamp;
                        if (v < -clamp) v = -clamp;
                    }
                }
            }
            wf.set_input_current(iIn);

            for (size_t secIdx = 0; secIdx + 1 < turnsRatios.size(); ++secIdx) {
                // turnsRatios[0] is the demagnetization-winding ratio;
                // physical secondaries are indexed 0..numSec-1 in the
                // netlist (vout0, vout1, ...) and start at
                // turnsRatios[1]. Iterating turnsRatios.size() (the old
                // bound) produced an empty extra entry per OP.
                wf.get_mutable_output_voltages().push_back(getWaveform("vout" + std::to_string(secIdx)));
                // Reconstruct Iout(t) = Vout(t)/Rload (DC by design).
                // Mirrors Buck.cpp / IsolatedBuck.cpp / AHB §5.1. Avoids
                // the i(Vsec_sense#branch) approach because that branch
                // is sensitive to inductor-startup transients within the
                // extraction window for some operating points.
                const double Vo_i = opPoint.get_output_voltages()[secIdx];
                const double Io_i = opPoint.get_output_currents()[secIdx];
                const double Rload_i = std::max(Vo_i / Io_i, 1e-3);
                Waveform ioutWf = getWaveform("vout" + std::to_string(secIdx));
                auto& ioutData = ioutWf.get_mutable_data();
                for (auto& v : ioutData) v = v / Rload_i;
                wf.get_mutable_output_currents().push_back(ioutWf);
            }
            
            results.push_back(wf);
        }
    }
    
    // Restore original value
    set_num_periods_to_extract(originalNumPeriodsToExtract);
    
    return results;
}
} // namespace OpenMagnetics
