#include "converter_models/TwoSwitchForward.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double TwoSwitchForward::get_total_reflected_secondary_current(const ForwardOperatingPoint& forwardOperatingPoint, const std::vector<double>& turnsRatios, double rippleRatio) {
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

    OperatingPoint TwoSwitchForward::process_operating_points_for_input_voltage(double inputVoltage, const ForwardOperatingPoint& outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance, double mainOutputInductance) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double mainOutputCurrent = outputOperatingPoint.get_output_currents()[0];
        double mainOutputVoltage = outputOperatingPoint.get_output_voltages()[0];
        double mainSecondaryTurnsRatio = turnsRatios[0];
        double diodeVoltageDrop = get_diode_voltage_drop();

        auto dutyCycle = get_maximum_duty_cycle();

        // Assume CCM
        auto period = 1.0 / switchingFrequency;
        auto t1 = period / 2 * (mainOutputVoltage + diodeVoltageDrop) / (inputVoltage / mainSecondaryTurnsRatio);

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
            t1 = sqrt(2 * mainOutputCurrent * mainOutputInductance * (mainOutputVoltage + diodeVoltageDrop) / (switchingFrequency * (inputVoltage / mainSecondaryTurnsRatio - diodeVoltageDrop - mainOutputVoltage) * (inputVoltage / mainSecondaryTurnsRatio)));
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
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {

            double secondaryCurrentPeakToPeak = maximumSecondaryCurrents[secondaryIndex] - minimumSecondaryCurrents[secondaryIndex];

            double minimumSecondaryVoltage = -(inputVoltage + 2 * diodeVoltageDrop) / turnsRatios[secondaryIndex];
            double maximumSecondaryVoltage = inputVoltage / turnsRatios[secondaryIndex];
            double secondaryVoltagePeakToPeak = maximumSecondaryVoltage - minimumSecondaryVoltage;
            double secondaryVoltageOffset = maximumSecondaryVoltage + minimumSecondaryVoltage;

            auto currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, secondaryCurrentPeakToPeak, switchingFrequency, dutyCycle, minimumSecondaryCurrents[secondaryIndex], 0);
            auto voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, secondaryVoltagePeakToPeak, switchingFrequency, dutyCycle, secondaryVoltageOffset, deadTime);
            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Secondary " + std::to_string(secondaryIndex));
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    bool TwoSwitchForward::run_checks(bool assert) {
        return ForwardConverterUtils::run_checks_common(this, assert);
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
                double turnsRatio = maximumInputVoltage * dutyCycle / (forwardOperatingPoint.get_output_voltages()[secondaryIndex] + diodeVoltageDrop);
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
            ForwardConverterUtils::create_isolation_sides(get_operating_points()[0].get_output_currents().size(), false));
        designRequirements.set_topology(Topologies::TWO_SWITCH_FORWARD_CONVERTER);
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

        ForwardConverterUtils::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        std::vector<double> outputInductancePerSecondary;

    for (size_t secondaryIndex = 0; secondaryIndex < turnsRatios.size(); ++secondaryIndex) {
        outputInductancePerSecondary.push_back(get_output_inductance(turnsRatios[secondaryIndex], secondaryIndex));
    }


        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];

            for (size_t forwardOperatingPointIndex = 0; forwardOperatingPointIndex < get_operating_points().size(); ++forwardOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[forwardOperatingPointIndex], turnsRatios, magnetizingInductance, outputInductancePerSecondary[0]);

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

    std::vector<OperatingPoint> TwoSwitchForward::process_operating_points(Magnetic magnetic) {
        TwoSwitchForward::run_checks(_assertErrors);

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(_magnetizingInductanceModel);;  // hardcoded
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


        ForwardConverterUtils::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

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
            ForwardConverterUtils::create_isolation_sides(get_operating_points()[0].get_output_currents().size(), false));
        designRequirements.set_topology(Topologies::TWO_SWITCH_FORWARD_CONVERTER);

        inputs.set_design_requirements(designRequirements);

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
    
    std::string TwoSwitchForward::generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames_;
    ForwardConverterUtils::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames_);
        
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
        
        size_t numSecondaries = turnsRatios.size();
        
        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        
        // Simulation timing
        int periodsToExtract = get_num_periods_to_extract();
        int numSteadyStatePeriods = get_num_steady_state_periods();
        const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;  // Steady state + extraction
        double simTime = numPeriodsTotal * period;
        double startTime = numSteadyStatePeriods * period;  // Start extracting after steady state
        double stepTime = period / 100;  // Reduced from 200 for faster simulation
        
        circuit << "* Two-Switch Forward Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, " << numSecondaries << " secondaries\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // Two-Switch Forward uses two switches: S1 (high-side) and S2 (low-side)
        // Both switches turn on and off together
        circuit << "* PWM Switches (both controlled together)\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 5 0 10n 10n " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=2.5 VH=0.5\n";
        circuit << ".model DIDEAL D(IS=1e-14 RS=0.01 CJO=1e-12)\n\n";
        
        // High-side switch S1 with clamping diode D1
        circuit << "* High-side switch S1 and clamp diode D1\n";
        circuit << "S1 vin_dc sw1_out pwm_ctrl 0 SW1\n";
        circuit << "D1 sw1_out vin_dc DIDEAL\n\n";
        
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
        circuit << "D2 0 pri_gnd DIDEAL\n\n";
        
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
            circuit << "Lout" << secIdx << " sec" << secIdx << "_l_in vout" << secIdx << " " << std::scientific << outputInductance << std::fixed << "\n";
            
            double loadResistance = outputVoltage / outputCurrent;
            circuit << "Cout" << secIdx << " vout" << secIdx << " 0 100u IC=" << outputVoltage << "\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << " UIC\n\n";
        
        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(pri_in) i(Vpri_sense)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) i(Vsec_sense" << secIdx << ") v(vout" << secIdx << ")";
        }
        circuit << "\n\n";
        
        // Options - relaxed tolerances for faster convergence
        circuit << ".options RELTOL=0.003 ABSTOL=1e-8 VNTOL=1e-5 TRTOL=10 ITL1=500 ITL4=100\n";
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
        ForwardConverterUtils::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        size_t numSecondaries = turnsRatios.size();
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto forwardOpPoint = get_operating_points()[opIndex];
                
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = forwardOpPoint.get_switching_frequency();
                
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = 1;
                config.keepTempFiles = false;
                
                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                NgspiceRunner::WaveformNameMapping waveformMapping;
                
                // Primary winding
                waveformMapping.push_back({{"voltage", "pri_in"}, {"current", "vpri_sense#branch"}});
                
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
    
    std::vector<OperatingPoint> TwoSwitchForward::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
        // For Two Switch Forward converter, topology waveforms are the same as operating points
        // The operating point already contains all winding voltages and currents
        return simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
    }
} // namespace OpenMagnetics
