#include "converter_models/TwoSwitchForward.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double TwoSwitchForward::get_total_reflected_secondary_current(const TopologyExcitation& forwardOperatingPoint, const std::vector<double>& turnsRatios, double rippleRatio) {
        double totalReflectedSecondaryCurrent = 0;

        if (turnsRatios.size() != forwardOperatingPoint.get_output_currents().size()) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Turns ratios must have same positions as outputs");
        }

        for (size_t secondaryIndex = 0; secondaryIndex < forwardOperatingPoint.get_output_currents().size(); ++secondaryIndex) {
            totalReflectedSecondaryCurrent += forwardOperatingPoint.get_output_currents()[secondaryIndex] / turnsRatios[secondaryIndex] * rippleRatio;
        }

        return totalReflectedSecondaryCurrent;
    }

    double TwoSwitchForward::get_maximum_duty_cycle() {
        if (get_duty_cycle()) {
            return get_duty_cycle().value();
        }
        else {
            return 0.45;
        }
    }

    TwoSwitchForward::TwoSwitchForward(const json& j) {
        from_json(j, *this);
    }

    AdvancedTwoSwitchForward::AdvancedTwoSwitchForward(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint TwoSwitchForward::process_operating_points_for_input_voltage(double inputVoltage, const TopologyExcitation& outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance, double mainOutputInductance) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double mainOutputCurrent = outputOperatingPoint.get_output_currents()[0];
        double mainOutputVoltage = outputOperatingPoint.get_output_voltages()[0];
        double mainSecondaryTurnsRatio = turnsRatios[0];
        double diodeVoltageDrop = get_diode_voltage_drop();

        // Assume CCM
        auto period = 1.0 / switchingFrequency;
        // Switch ON time from volt-second balance on the output inductor:
        //   Vout = Vin/n · D − Vd  ⇒  D = (Vout+Vd)/(Vin/n)
        // Two-Switch Forward has NO inherent /2 (unlike push-pull/half-bridge
        // where D is per-half-cycle).  Matches the SPICE netlist tOn
        // computation (generate_ngspice_circuit) so SPICE and analytical
        // use the same duty cycle and waveform timing.
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
            minimumPrimaryCurrent += minimumSecondaryCurrent / turnsRatios[secondaryIndex];
            maximumPrimaryCurrent += maximumSecondaryCurrent / turnsRatios[secondaryIndex];
        }

        if (minimumPrimaryCurrent < 0) { // DCM
            double sqrtArg = 2 * mainOutputCurrent * mainOutputInductance * (mainOutputVoltage + diodeVoltageDrop) / (switchingFrequency * (inputVoltage / mainSecondaryTurnsRatio - diodeVoltageDrop - mainOutputVoltage) * (inputVoltage / mainSecondaryTurnsRatio));
            if (sqrtArg < 0) {
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "TwoSwitchForward: negative value under sqrt in DCM t1 calculation");
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

                minimumPrimaryCurrent += minimumSecondaryCurrents[secondaryIndex] / turnsRatios[secondaryIndex];
                maximumPrimaryCurrent += maximumSecondaryCurrents[secondaryIndex] / turnsRatios[secondaryIndex];
            }
        }

        double minimumPrimarySideTransformerCurrentT1 = minimumPrimaryCurrent;
        double maximumPrimarySideTransformerCurrentT1 = maximumPrimaryCurrent;
        double minimumPrimarySideTransformerVoltage = -inputVoltage - 2 * diodeVoltageDrop;
        double maximumPrimarySideTransformerVoltage = inputVoltage;

        double minimumPrimarySideTransformerCurrentTd = 0;
        double maximumPrimarySideTransformerCurrentTd = magnetizationCurrent;

        double td = t1;
        double deadTime = period - t1 - td;

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;
            // Current
            {
                if (minimumPrimaryCurrent > 0) { // CCM
                    std::vector<double> data = {
                        0,
                        minimumPrimarySideTransformerCurrentT1,
                        maximumPrimarySideTransformerCurrentT1,
                        maximumPrimarySideTransformerCurrentTd,
                        minimumPrimarySideTransformerCurrentTd,
                        0,
                        0
                    };
                    std::vector<double> time = {
                        0,
                        0,
                        t1,
                        t1,
                        t1 + td,
                        period,
                        period
                    };
                    currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                    currentWaveform.set_data(data);
                    currentWaveform.set_time(time);
                }
                else  { // DCM
                    std::vector<double> data = {
                        minimumPrimarySideTransformerCurrentT1,
                        maximumPrimarySideTransformerCurrentT1,
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
            }
            // Voltage
            {
                std::vector<double> data = {
                    0,
                    maximumPrimarySideTransformerVoltage,
                    maximumPrimarySideTransformerVoltage,
                    minimumPrimarySideTransformerVoltage,
                    minimumPrimarySideTransformerVoltage,
                    0,
                    0
                };
                std::vector<double> time = {
                    0,
                    0,
                    t1,
                    t1,
                    t1 + td,
                    t1 + td,
                    period
                };
                voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                voltageWaveform.set_data(data);
                voltageWaveform.set_time(time);
            }

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "First primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }
        // Use actual duty cycle (t1/period) for secondary waveforms, not the fixed maximum
        double actualDutyCycle = t1 / period;
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {

            double secondaryCurrentPeakToPeak = maximumSecondaryCurrents[secondaryIndex] - minimumSecondaryCurrents[secondaryIndex];

            double minimumSecondaryVoltage = -(inputVoltage + 2 * diodeVoltageDrop) / turnsRatios[secondaryIndex];
            double maximumSecondaryVoltage = inputVoltage / turnsRatios[secondaryIndex];
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

        // ---- Diagnostics (per-OP, last-call-wins across multiple OPs) ----
        lastPrimaryPeakCurrent     = maximumPrimaryCurrent;
        lastMagnetizingPeakCurrent = magnetizationCurrent;
        if (!maximumSecondaryCurrents.empty()) {
            lastSecondaryPeakCurrent = maximumSecondaryCurrents[0];
        }
        lastIsCcm = (minimumPrimaryCurrent > 0);

        // Per-OP diagnostic snapshot.
        perOpMaximumDutyCycle.push_back(lastMaximumDutyCycle);
        perOpComputedMagnetizingInductance.push_back(lastComputedMagnetizingInductance);
        perOpComputedSecondaryTurnsRatio.push_back(lastComputedSecondaryTurnsRatio);
        perOpPrimaryPeakCurrent.push_back(lastPrimaryPeakCurrent);
        perOpSecondaryPeakCurrent.push_back(lastSecondaryPeakCurrent);
        perOpMagnetizingPeakCurrent.push_back(lastMagnetizingPeakCurrent);
        perOpIsCcm.push_back(lastIsCcm);

        return operatingPoint;
    }

    bool TwoSwitchForward::run_checks(bool assert) {
        return Topology::run_checks_common(this, assert);
    }

    DesignRequirements TwoSwitchForward::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double diodeVoltageDrop = get_diode_voltage_drop();
        double dutyCycle = get_maximum_duty_cycle();

        // Turns ratio calculation
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size(), 0);

        for (auto forwardOperatingPoint : get_operating_points()) {
            for (size_t secondaryIndex = 0; secondaryIndex < forwardOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                // Size turns ratio at Vin_MIN (worst case for required D).
                // See SingleSwitchForward.cpp for full rationale.
                double turnsRatio = minimumInputVoltage * dutyCycle / (forwardOperatingPoint.get_output_voltages()[secondaryIndex] + diodeVoltageDrop);
                turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatio);
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
            Topology::create_isolation_sides(get_operating_points()[0].get_output_currents().size(), false));
        designRequirements.set_topology(MAS::Topology::TWO_SWITCH_FORWARD_CONVERTER);

        // ---- Diagnostics ----
        lastMaximumDutyCycle              = dutyCycle;
        lastComputedMagnetizingInductance = minimumNeededInductance;
        lastComputedSecondaryTurnsRatio   = turnsRatios.empty() ? 0.0 : turnsRatios[0];

        return designRequirements;
    }

    double TwoSwitchForward::get_output_inductance(double secondaryTurnsRatio, size_t outputIndex) {
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

    std::vector<OperatingPoint> TwoSwitchForward::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        // Clear per-OP diagnostic vectors so the wizard table reflects this run only.
        perOpName.clear();
        perOpMaximumDutyCycle.clear();
        perOpComputedMagnetizingInductance.clear();
        perOpComputedSecondaryTurnsRatio.clear();
        perOpPrimaryPeakCurrent.clear();
        perOpSecondaryPeakCurrent.clear();
        perOpMagnetizingPeakCurrent.clear();
        perOpIsCcm.clear();


        extraLoVoltageWaveforms.clear();
        extraLoCurrentWaveforms.clear();
        extraLoInductances.clear();

        std::vector<double> outputInductancePerSecondary;

        size_t numSecondaries = turnsRatios.size();  // TwoSwitchForward has no demag winding
        for (size_t secondaryIndex = 0; secondaryIndex < numSecondaries; ++secondaryIndex) {
            outputInductancePerSecondary.push_back(get_output_inductance(turnsRatios[secondaryIndex], secondaryIndex));
            extraLoInductances.push_back(outputInductancePerSecondary.back());
        }

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];

            for (size_t forwardOperatingPointIndex = 0; forwardOperatingPointIndex < get_operating_points().size(); ++forwardOperatingPointIndex) {
                auto& fop = get_operating_points()[forwardOperatingPointIndex];
                std::string opName = inputVoltagesNames[inputVoltageIndex];
                if (get_operating_points().size() > 1) {
                    opName += " · OP" + std::to_string(forwardOperatingPointIndex);
                }
                perOpName.push_back(opName);
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, fop, turnsRatios, magnetizingInductance, outputInductancePerSecondary[0]);

                // Capture Lo waveforms for each secondary
                double switchingFrequency = fop.get_switching_frequency();
                double period = 1.0 / switchingFrequency;
                double diodeVoltageDrop = get_diode_voltage_drop();
                double mainSecondaryTurnsRatio = turnsRatios[0];
                double mainOutputVoltage = fop.get_output_voltages()[0];
                double actualDutyCycle = (mainOutputVoltage + diodeVoltageDrop) / (inputVoltage / mainSecondaryTurnsRatio);
                // No silent clamp: violating maximumDutyCycle means the design
                // cannot deliver Vout at this Vin. 1 % tolerance for the
                // design-requirements rounding-to-2-decimals.
                {
                    const double dMax = get_maximum_duty_cycle();
                    constexpr double dutyTolerance = 0.01;
                    if (actualDutyCycle > dMax * (1.0 + dutyTolerance)) {
                        throw InvalidInputException(ErrorCode::INVALID_INPUT,
                            "TwoSwitchForward: required dutyCycle " + std::to_string(actualDutyCycle) +
                            " exceeds maximumDutyCycle " + std::to_string(dMax) +
                            " at Vin=" + std::to_string(inputVoltage) +
                            ". Increase turns ratio, raise Vin_min, or relax maximumDutyCycle.");
                    }
                }
                double tOn = actualDutyCycle * period;

                for (size_t s = 0; s < numSecondaries; ++s) {
                    double turnsRatioSec = turnsRatios[s];
                    double outputVoltage = fop.get_output_voltages()[s];
                    double outputCurrent = fop.get_output_currents()[s];
                    double currentRipple = get_current_ripple_ratio() * outputCurrent;

                    Waveform loCurrentWfm;
                    {
                        double iMin = outputCurrent - currentRipple / 2;
                        double iMax = outputCurrent + currentRipple / 2;
                        loCurrentWfm.set_data(std::vector<double>{iMin, iMax, iMax, iMin});
                        loCurrentWfm.set_time(std::vector<double>{0, tOn, tOn, period});
                        loCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
                    }

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

    std::vector<std::variant<Inputs, CAS::Inputs>> TwoSwitchForward::get_extra_components_inputs(
        ExtraComponentsMode /*mode*/,
        std::optional<Magnetic> /*magnetic*/)
    {
        if (extraLoInductances.empty() || extraLoVoltageWaveforms.empty())
            throw std::runtime_error("TwoSwitchForward::get_extra_components_inputs: call process_operating_points() first");

        size_t numSecondaries = extraLoInductances.size();
        size_t nOps = extraLoVoltageWaveforms.size() / numSecondaries;

        std::vector<std::variant<Inputs, CAS::Inputs>> result;

        for (size_t s = 0; s < numSecondaries; ++s) {
            double lo = extraLoInductances[s];
            if (lo <= 0)
                throw std::runtime_error("TwoSwitchForward::get_extra_components_inputs: computedOutputInductance <= 0 for secondary " + std::to_string(s));

            Inputs masInputs;
            DesignRequirements dr;

            DimensionWithTolerance inductance;
            inductance.set_nominal(lo);
            dr.set_magnetizing_inductance(inductance);

            std::string inductorName = (numSecondaries == 1) ? "outputInductor" : ("outputInductor_" + std::to_string(s + 1));
            dr.set_name(inductorName);
            dr.set_topology(MAS::Topology::TWO_SWITCH_FORWARD_CONVERTER);
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

    std::vector<OperatingPoint> TwoSwitchForward::process_operating_points(Magnetic magnetic) {
        TwoSwitchForward::run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedTwoSwitchForward::process() {
        TwoSwitchForward::run_checks(_assertErrors);

        Inputs inputs;

        double minimumNeededInductance = get_desired_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();
        std::vector<double> outputInductancePerSecondary;

        if (get_desired_output_inductances()) {
            outputInductancePerSecondary = get_desired_output_inductances().value();
        }
        else {
            for (size_t secondaryIndex = 0; secondaryIndex < turnsRatios.size(); ++secondaryIndex) {
                auto minimumOutputInductance = get_output_inductance(turnsRatios[secondaryIndex], secondaryIndex);
                outputInductancePerSecondary.push_back(minimumOutputInductance);
            }
        }

        if (turnsRatios.size() != get_operating_points()[0].get_output_currents().size()) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Turns ratios must have same positions as outputs");
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

    DesignRequirements AdvancedTwoSwitchForward::process_design_requirements() {
        // Issue M1: build DR directly from desired* fields. Do NOT chain to
        // TwoSwitchForward::process_design_requirements() — the parent
        // auto-sizes turns ratios and inductance, defeating the "Advanced"
        // override semantics.
        std::vector<double> turnsRatios = get_desired_turns_ratios();
        double minimumNeededInductance = get_desired_inductance();

        if (turnsRatios.size() != get_operating_points()[0].get_output_currents().size()) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Turns ratios must have same positions as outputs");
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
            Topology::create_isolation_sides(get_operating_points()[0].get_output_currents().size(), false));
        designRequirements.set_topology(MAS::Topology::TWO_SWITCH_FORWARD_CONVERTER);
        return designRequirements;
    }
    
    std::string TwoSwitchForward::generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {

        // Minimal spice_config() consolidation — PULSE rise/fall only.
        const auto cfg = spice_config();

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
        double diodeVoltageDrop = get_diode_voltage_drop();

        size_t numSecondaries = turnsRatios.size();

        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        // Compute actual duty cycle from voltage balance: D = (Vout + Vd) / (Vin / n)
        double mainOutputVoltage = opPoint.get_output_voltages()[0];
        double mainSecondaryTurnsRatio = turnsRatios[0];
        double dutyCycle = (mainOutputVoltage + diodeVoltageDrop) / (inputVoltage / mainSecondaryTurnsRatio);
        // No silent clamp: throw if the operating-point duty exceeds maxD,
        // with 1 % tolerance for design-requirements rounding.
        {
            const double dMax = get_maximum_duty_cycle();
            constexpr double dutyTolerance = 0.01;
            if (dutyCycle > dMax * (1.0 + dutyTolerance)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "TwoSwitchForward::generate_ngspice_circuit: required dutyCycle " +
                    std::to_string(dutyCycle) + " exceeds maximumDutyCycle " + std::to_string(dMax) +
                    " at Vin=" + std::to_string(inputVoltage) +
                    ". Increase turns ratio, raise Vin_min, or relax maximumDutyCycle.");
            }
        }
        double tOn = period * dutyCycle;
        
        // Simulation timing
        int periodsToExtract = get_num_periods_to_extract();
        int numSteadyStatePeriods = get_num_steady_state_periods();
        const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;  // Steady state + extraction
        double simTime = numPeriodsTotal * period;
        double startTime = numSteadyStatePeriods * period;  // Start extracting after steady state
        double stepTime = period / 200;  // Finer time resolution for convergence
        
        circuit << "* Two-Switch Forward Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, " << numSecondaries << " secondaries\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";

        // §8a.5 input-current sense: a zero-V source upstream of S1's
        // drain isolates the true switch-channel current. The
        // Two-Switch Forward has S1 high-side and S2 low-side switching
        // together; the entire on-time bus current flows through S1, so
        // a single ammeter there captures the converter's draw. D2
        // (low-side clamp) returns demag current back to vin_dc during
        // off-time and is intentionally left on vin_dc so it does NOT
        // appear in i(Vq1_sense). D1 (high-side clamp 0 -> sw1_out)
        // is across S1's output, so it also does not flow through
        // Vq1_sense.
        // Sign convention: + terminal is vin_dc, - is q1_drain, so
        // positive i(Vq1_sense) means current flowing OUT of vin_dc
        // INTO S1 (non-negative in CCM/DCM).
        circuit << "* §8a.5 input-current sense\n";
        circuit << "Vq1_sense vin_dc q1_drain 0\n\n";

        // Two-Switch Forward uses two switches: S1 (high-side) and S2 (low-side)
        // Both switches turn on and off together
        circuit << "* PWM Switches (both controlled together)\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << cfg.pwmRise << " " << cfg.pwmFall
                << " " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                << " RON=" << cfg.swModelRON
                << " ROFF=" << std::scientific << cfg.swModelROFF << std::defaultfloat << "\n";
        // Use simple diode model without junction capacitance for better convergence
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::defaultfloat;
        if (!cfg.diodeExtra.empty()) circuit << " " << cfg.diodeExtra;
        circuit << ")\n\n";

        // High-side switch S1 with clamping diode D1
        circuit << "* High-side switch S1 and clamp diode D1\n";
        circuit << "S1 q1_drain sw1_out pwm_ctrl 0 SW1\n";
        circuit << "D1 0 sw1_out DIDEAL\n\n";
        
        // Primary current sense
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense sw1_out pri_in 0\n\n";
        
        // Transformer primary
        circuit << "* Forward Transformer\n";
        circuit << "Lpri pri_in pri_gnd " << std::scientific << magnetizingInductance << std::fixed << "\n";
        
        // Secondary windings
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double secondaryInductance = magnetizingInductance / (turnsRatios[secIdx] * turnsRatios[secIdx]);
            circuit << "Lsec" << secIdx << " sec" << secIdx << "_in 0 " << std::scientific << secondaryInductance << std::fixed << "\n";
        }
        
        // Couple all windings together using pairwise K statements
        circuit << "* Coupling: All windings coupled pairwise\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "Kpri_sec" << secIdx << " Lpri Lsec" << secIdx << " 0.9999\n";
        }
        for (size_t i = 0; i < numSecondaries; ++i) {
            for (size_t j = i + 1; j < numSecondaries; ++j) {
                circuit << "Ksec" << i << "_sec" << j << " Lsec" << i << " Lsec" << j << " 0.9999\n";
            }
        }
        circuit << "\n";
        
        // Low-side switch S2 with clamping diode D2
        circuit << "* Low-side switch S2 and clamp diode D2\n";
        circuit << "S2 pri_gnd 0 pwm_ctrl 0 SW1\n";
        circuit << "D2 pri_gnd vin_dc DIDEAL\n\n";

        // §8a.5 differential winding-voltage probes. The Two-Switch
        // primary winding sits between pri_in (top of S1/Vpri_sense
        // stack) and pri_gnd (top of S2) — NEITHER endpoint is the
        // netlist reference (ground), so the bare-node v(pri_in) probe
        // used previously lumped the winding voltage with the variable
        // pri_gnd offset (= 0 during on-time when S2 conducts, = Vin
        // during off-time when D2 clamps). The E-source Evpri_w exposes
        // v(pri_in) - v(pri_gnd) as a pure winding-only differential.
        // The secondary windings Lsec<i> are wired sec<i>_in -> 0 so
        // their bare-node probes are already winding-only — Bug C does
        // not require secondary E-sources here.
        circuit << "* §8a.5 differential primary winding probe\n";
        circuit << "Evpri_w vpri_w 0 pri_in pri_gnd 1\n\n";
        
        // Output stages for each secondary
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "* Secondary " << secIdx << " output stage\n";
            
            // Forward rectifier
            circuit << "Dfwd" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_rect DIDEAL\n";
            
            // Freewheeling diode
            circuit << "Dfw" << secIdx << " 0 sec" << secIdx << "_rect DIDEAL\n";
            
            // Snubber resistors for convergence
            circuit << "Rsnub_fwd" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_rect 1MEG\n";
            circuit << "Rsnub_fw" << secIdx << " 0 sec" << secIdx << "_rect 1MEG\n";
            
            double outputVoltage = opPoint.get_output_voltages()[secIdx];
            double outputCurrent = opPoint.get_output_currents()[secIdx];
            double outputInductance = get_output_inductance(turnsRatios[secIdx], secIdx);
            
            circuit << "Vsec_sense" << secIdx << " sec" << secIdx << "_rect sec" << secIdx << "_l_in 0\n";
            // Add small series resistance to output inductor for damping (helps convergence)
            circuit << "Rlout" << secIdx << " sec" << secIdx << "_l_in lout" << secIdx << "_node 0.01\n";
            circuit << "Lout" << secIdx << " lout" << secIdx << "_node vout" << secIdx << " " << std::scientific << outputInductance << std::fixed << "\n";
            
            double loadResistance = outputVoltage / outputCurrent;
            circuit << "Cout" << secIdx << " vout" << secIdx << " 0 100u IC=" << outputVoltage << "\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << " UIC\n\n";
        
        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(pri_in) v(pri_gnd) v(vpri_w) i(Vpri_sense) v(vin_dc) i(Vin) i(Vq1_sense)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) i(Vsec_sense" << secIdx << ") v(vout" << secIdx << ")";
        }
        circuit << "\n\n";
        
        // Options - tighter tolerances for better convergence with high voltage ratios
        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::defaultfloat
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << ".ic v(vout" << secIdx << ")=" << opPoint.get_output_voltages()[secIdx] << "\n";
        }
        circuit << "\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> TwoSwitchForward::simulate_and_extract_operating_points(
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
        
        size_t numSecondaries = turnsRatios.size();
        
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
                
                NgspiceRunner::WaveformNameMapping waveformMapping;

                // Primary winding — §8a.5 Bug C: read differential
                // v(vpri_w) (Evpri_w exposes pri_in - pri_gnd) instead
                // of bare v(pri_in), because pri_gnd is at a non-zero
                // reference (0 during S2 on-time, Vin during off-time
                // via D2 clamp), so v(pri_in) lumps the winding with
                // that offset.
                waveformMapping.push_back({{"voltage", "vpri_w"}, {"current", "vpri_sense#branch"}});
                
                // Secondary windings
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    std::string voltageName = "sec" + std::to_string(secIdx) + "_in";
                    std::string currentName = "vsec_sense" + std::to_string(secIdx) + "#branch";
                    waveformMapping.push_back({{"voltage", voltageName}, {"current", currentName}});
                }
                
                std::vector<std::string> windingNames;
                windingNames.push_back("Primary");
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    windingNames.push_back("Secondary " + std::to_string(secIdx));
                }
                
                std::vector<bool> flipCurrentSign(1 + numSecondaries, false);
                
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
    
    std::vector<ConverterWaveforms> TwoSwitchForward::simulate_and_extract_topology_waveforms(const std::vector<double>& turnsRatios,
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
            //   input_current = i(Vq1_sense) — zero-V sense source
            //     upstream of S1's drain. Two-Switch Forward has S1
            //     high-side + S2 low-side switching together; the full
            //     on-time bus current flows through S1, so a single
            //     ammeter captures the converter's draw. i(Vin) (the
            //     prior probe) sums Vq1_sense PLUS the D2-clamp return
            //     current (D2 dumps demag energy back into vin_dc
            //     during off-time) — not the converter's input draw.
            //     Sign: + terminal of Vq1_sense is vin_dc, so positive
            //     i(Vq1_sense) is current flowing OUT of vin_dc INTO
            //     S1 (non-negative in CCM/DCM).
            wf.set_input_voltage(getWaveform("vin_dc"));

            Waveform iIn = getWaveform("vq1_sense#branch");
            if (iIn.get_data().empty()) {
                throw std::runtime_error(
                    "TwoSwitchForward simulate_and_extract_topology_waveforms: "
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

            // TwoSwitchForward has no demag winding — every entry in
            // turnsRatios is a physical secondary (numSec = size()).
            for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
                wf.get_mutable_output_voltages().push_back(getWaveform("vout" + std::to_string(secIdx)));
                // Reconstruct Iout(t) = Vout(t)/Rload (DC by design).
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
