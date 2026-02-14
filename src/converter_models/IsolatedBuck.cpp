#include "converter_models/IsolatedBuck.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double IsolatedBuck::calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency) {
        auto dutyCycle = outputVoltage / inputVoltage * efficiency;
        if (dutyCycle >= 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Duty cycle must be smaller than 1");
        }
        return dutyCycle;
    }

    IsolatedBuck::IsolatedBuck(const json& j) {
        from_json(j, *this);
    }

    AdvancedIsolatedBuck::AdvancedIsolatedBuck(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint IsolatedBuck::processOperatingPointsForInputVoltage(double inputVoltage, const IsolatedBuckOperatingPoint& outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance) {

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

        auto primaryVoltageMaximum = inputVoltage - primaryOutputVoltage;
        auto primaryVoltageMinimum = -primaryOutputVoltage;
        auto primaryVoltagePeaktoPeak = primaryVoltageMaximum - primaryVoltageMinimum;

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, primaryOutputCurrent, 0);
            voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, 0);

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

            auto secondaryVoltageMaximum = (inputVoltage - primaryOutputVoltage) / turnsRatios[secondaryIndex] - diodeVoltageDrop;
            auto secondaryVoltageMinimum = -primaryOutputVoltage / turnsRatios[secondaryIndex] + diodeVoltageDrop;

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
                    secondaryVoltageMinimum,
                    secondaryVoltageMinimum,
                    secondaryVoltageMaximum,
                    secondaryVoltageMaximum
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
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        for (size_t isolatedbuckOperatingPointIndex = 1; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
            if (get_operating_points()[isolatedbuckOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[isolatedbuckOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
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

    DesignRequirements IsolatedBuck::process_design_requirements() {
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Missing both current ripple ratio and maximum switch current");
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
                auto totalCurrentInPoint = isolatedbuckOperatingPoint.get_output_currents()[0];
                double totalReflectedSecondaryCurrent = 0;
                for (size_t secondaryIndex = 0; secondaryIndex < isolatedbuckOperatingPoint.get_output_currents().size() - 1; ++secondaryIndex) {
                    totalReflectedSecondaryCurrent += isolatedbuckOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex];
                }
                totalCurrentInPoint += totalReflectedSecondaryCurrent;
                maximumOutputCurrent = std::max(maximumOutputCurrent, totalCurrentInPoint);
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

    std::vector<OperatingPoint> IsolatedBuck::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

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

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
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


        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

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

    std::vector<ConverterWaveforms> IsolatedBuck::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods) {

    // Save original value and set the requested number of periods
    int originalNumPeriodsToExtract = get_num_periods_to_extract();
    set_num_periods_to_extract(static_cast<int>(numberOfPeriods));

    auto operatingPoints = process_operating_points(turnsRatios, magnetizingInductance);
    std::vector<ConverterWaveforms> results;

    for (auto& op : operatingPoints) {
        ConverterWaveforms wf;
        wf.set_switching_frequency(op.get_excitations_per_winding()[0].get_frequency());
        if (op.get_name()) {
            wf.set_operating_point_name(op.get_name().value());
        }
        auto& priExc = op.get_excitations_per_winding()[0];
        if (priExc.get_voltage() && priExc.get_voltage()->get_waveform()) {
            wf.set_input_voltage(priExc.get_voltage()->get_waveform().value());
        }
        if (priExc.get_current() && priExc.get_current()->get_waveform()) {
            wf.set_input_current(priExc.get_current()->get_waveform().value());
        }
        for (size_t i = 1; i < op.get_excitations_per_winding().size(); ++i) {
            auto& secExc = op.get_excitations_per_winding()[i];
            if (secExc.get_voltage() && secExc.get_voltage()->get_waveform()) {
                wf.get_mutable_output_voltages().push_back(secExc.get_voltage()->get_waveform().value());
            }
            if (secExc.get_current() && secExc.get_current()->get_waveform()) {
                wf.get_mutable_output_currents().push_back(secExc.get_current()->get_waveform().value());
            }
        }
        results.push_back(wf);
    }

    // Restore original value
    set_num_periods_to_extract(originalNumPeriodsToExtract);

    return results;
}

    std::string IsolatedBuck::generate_ngspice_circuit(
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
        double primaryOutputVoltage = opPoint.get_output_voltages()[0];
        double primaryOutputCurrent = opPoint.get_output_currents()[0];

        double efficiency = 1;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }

        double dutyCycle = calculate_duty_cycle(inputVoltage, primaryOutputVoltage, efficiency);

        size_t numSecondaries = turnsRatios.size();

        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;

        // Simulation: run steady-state periods for settling, then extract the last N periods
        int periodsToExtract = get_num_periods_to_extract();
        int numSteadyStatePeriods = get_num_steady_state_periods();
        const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;
        double simTime = numPeriodsTotal * period;
        double startTime = numSteadyStatePeriods * period;
        double stepTime = period / 200;

        circuit << "* Isolated Buck (Flybuck) Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vout_pri=" << primaryOutputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, " << numSecondaries << " secondaries\n\n";

        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";

        // High-side switch (synchronous buck)
        circuit << "* High-side Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 5 0 10n 10n " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=2.5 VH=0.5 RON=0.01 ROFF=1e6\n";
        circuit << "S1 vin_dc sw_node pwm_ctrl 0 SW1\n\n";

        // Low-side switch (synchronous rectification)
        circuit << "* Low-side Switch (Synchronous Rectifier)\n";
        circuit << "Vpwm_inv pwm_inv 0 PULSE(5 0 0 10n 10n " << tOn << " " << period << ")\n";
        circuit << ".model SW2 SW VT=2.5 VH=0.5 RON=0.01 ROFF=1e6\n";
        circuit << "S2 sw_node 0 pwm_inv 0 SW2\n\n";

        // Coupled inductor (Primary = buck inductor)
        // Using proper buck topology: series inductor with parallel DC path resistor
        circuit << "* Coupled Inductor (Primary = Buck Inductor)\n";
        
        // Primary current sense - measure current entering the inductor
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense sw_node pri_in 0\n";
        
        // Primary inductor: from switch node to output (series inductor like buck converter)
        // This is the proper topology - inductor in series with output, not to ground
        circuit << "Lpri pri_in vpri_out " << std::scientific << magnetizingInductance << std::fixed << "\n\n";

        // Secondary windings
        // For flyback action, secondary windings referenced to ground
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double secondaryInductance = magnetizingInductance / (turnsRatios[secIdx] * turnsRatios[secIdx]);
            // Secondary winding: ground to sec_in
            circuit << "Lsec" << secIdx << " 0 sec" << secIdx << "_in " << std::scientific << secondaryInductance << std::fixed << "\n";
        }

        // Couple primary to each secondary with high coupling
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "Kpri_sec" << secIdx << " Lpri Lsec" << secIdx << " 0.99\n";
        }
        // Couple secondaries to each other
        for (size_t i = 0; i < numSecondaries; ++i) {
            for (size_t j = i + 1; j < numSecondaries; ++j) {
                circuit << "Ksec" << i << "_" << j << " Lsec" << i << " Lsec" << j << " 0.99\n";
            }
        }
        circuit << "\n";

        // Diode model
        circuit << "* Diode model\n";
        circuit << ".model DIDEAL D(IS=1e-14 RS=0.01 N=1.0)\n\n";

        // Primary output (buck output) - directly from inductor output
        // In this topology, vpri_out is the output node (after the inductor)
        circuit << "* Primary Output Stage (Buck - non-isolated)\n";
        double primaryLoadResistance = primaryOutputVoltage / primaryOutputCurrent;
        circuit << "Cpri vpri_out 0 100u IC=" << primaryOutputVoltage << "\n";
        circuit << "Rload_pri vpri_out 0 " << primaryLoadResistance << "\n\n";

        // Secondary output stages (isolated outputs)
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "* Secondary " << secIdx << " output stage (isolated)\n";
            // Add series resistance between winding and diode (models winding resistance)
            circuit << "Rsec" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_node 0.01\n";
            // Secondary rectifier: diode from secondary winding node to rectified node
            // Anode at secX_node (winding side), cathode at secX_rect (output side)
            circuit << "Dsec" << secIdx << " sec" << secIdx << "_node sec" << secIdx << "_rect DIDEAL\n";
            circuit << "Vsec_sense" << secIdx << " sec" << secIdx << "_rect vout" << secIdx << " 0\n";

            double outputVoltage = opPoint.get_output_voltages()[secIdx + 1];
            double outputCurrent = opPoint.get_output_currents()[secIdx + 1];
            double loadResistance = outputVoltage / outputCurrent;
            circuit << "Cout" << secIdx << " vout" << secIdx << " 0 100u IC=" << outputVoltage << "\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }

        // Transient Analysis with UIC (Use Initial Conditions)
        // UIC skips the DC operating point calculation and uses the specified initial conditions
        // This is necessary for circuits with inductors that don't have DC paths to ground
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << " UIC\n\n";

        // Save signals
        // For Flybuck topology:
        // - pri_in: inductor input (after current sense)
        // - vpri_out: primary output (inductor output / non-isolated output)
        // - voutX: secondary outputs (isolated outputs)
        // - secX_in: secondary winding inputs
        // - secX_rect: secondary rectified node (before output cap)
        circuit << "* Output signals\n";
        circuit << ".save v(sw_node) v(pri_in) v(vpri_out) i(Vpri_sense)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) v(sec" << secIdx << "_rect) i(Vsec_sense" << secIdx << ") v(vout" << secIdx << ")";
        }
        circuit << "\n\n";

        // Options (matching other converters for convergence)
        circuit << ".options RELTOL=0.001 ABSTOL=1e-9 VNTOL=1e-6 ITL1=1000 ITL4=1000\n";
        circuit << ".ic v(vpri_out)=" << primaryOutputVoltage << "\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << ".ic v(vout" << secIdx << ")=" << opPoint.get_output_voltages()[secIdx + 1] << "\n";
        }
        circuit << "\n";

        circuit << ".end\n";

        return circuit.str();
    }

    std::vector<OperatingPoint> IsolatedBuck::simulate_and_extract_operating_points(
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

        size_t numSecondaries = turnsRatios.size();

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto ibOpPoint = get_operating_points()[opIndex];

                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);

                double switchingFrequency = ibOpPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);

                if (!simResult.success) {
                    throw std::runtime_error("ngspice simulation failed: " + simResult.errorMessage);
                }

                // Define waveform name mapping
                NgspiceRunner::WaveformNameMapping waveformMapping;
                
                // Primary winding - voltage at inductor node (switching node), current through sense resistor
                // For Flybuck: pri_in is the inductor/switching node
                waveformMapping.push_back({{"voltage", "pri_in"}, {"current", "vpri_sense#branch"}});
                
                // Secondary windings
                for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                    // Use secX_rect for secondary voltage (rectified node before output cap)
                    std::string voltageName = "sec" + std::to_string(secIdx) + "_rect";
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
                    ibOpPoint.get_ambient_temperature(),
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

} // namespace OpenMagnetics
