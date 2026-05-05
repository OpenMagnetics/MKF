#include "converter_models/PowerFactorCorrection.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <numbers>
#include <numeric>
#include <algorithm>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    PowerFactorCorrection::PowerFactorCorrection(const json& j) : Topology(j) {
        // Required fields per the MAS schema. Guard up front so users see
        // PFC-specific messages instead of a less specific j.at() throw.
        if (!j.contains("inputVoltage")) {
            throw std::runtime_error("PowerFactorCorrection: 'inputVoltage' is required");
        }
        if (!j.contains("outputVoltage")) {
            throw std::runtime_error("PowerFactorCorrection: 'outputVoltage' is required");
        }
        if (!j.contains("outputPower")) {
            throw std::runtime_error("PowerFactorCorrection: 'outputPower' is required");
        }
        if (!j.contains("switchingFrequency")) {
            throw std::runtime_error("PowerFactorCorrection: 'switchingFrequency' is required");
        }
        if (!j.contains("ambientTemperature")) {
            throw std::runtime_error("PowerFactorCorrection: 'ambientTemperature' is required");
        }

        // Delegate to the MAS-generated parser. It maps the optional
        // 'mode' field through the PfcModes enum (continuousConductionMode,
        // discontinuousConductionMode, criticalConductionMode, transitionMode).
        MAS::from_json(j, static_cast<MAS::PowerFactorCorrection&>(*this));

        // MKF-only knob (not in the MAS schema).
        if (j.contains("numberOfPeriods")) {
            _numberOfPeriods = j["numberOfPeriods"].get<int>();
        }
    }

    bool PowerFactorCorrection::run_checks(bool assert) {
        bool valid = true;

        auto inputVoltage = get_input_voltage();
        double outputVoltage = get_output_voltage();
        double efficiency   = get_efficiency_or_default();

        // Check that output voltage is greater than peak input voltage
        double vinPeakMax = resolve_dimensional_values(inputVoltage, DimensionalValues::MAXIMUM) * std::sqrt(2);
        if (outputVoltage <= vinPeakMax) {
            if (assert) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "PFC output voltage must be greater than peak input voltage. Vout: " +
                    std::to_string(outputVoltage) + " <= Vin_peak_max: " + std::to_string(vinPeakMax));
            }
            valid = false;
        }

        // Check efficiency is reasonable
        if (efficiency <= 0 || efficiency > 1) {
            if (assert) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Efficiency must be between 0 and 1. Got: " + std::to_string(efficiency));
            }
            valid = false;
        }

        return valid;
    }

    double PowerFactorCorrection::calculate_duty_cycle(double vinPeak, double /*vout*/) {
        // For boost PFC: D = 1 - Vin/Vout
        // Accounting for diode drop: D = 1 - Vin/(Vout + Vd)
        return 1.0 - vinPeak / (get_output_voltage() + get_diode_voltage_drop_or_default());
    }

    double PowerFactorCorrection::calculate_inductance_ccm() {
        // For CCM PFC, worst case (maximum inductance requirement) is at minimum input voltage
        // at the peak of the AC line (where current is maximum)

        auto inputVoltage = get_input_voltage();
        double vinRmsMin  = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);
        double outputVoltage = get_output_voltage();

        // Calculate duty cycle at minimum input voltage peak
        double D = calculate_duty_cycle(vinPeakMin, outputVoltage);

        // Calculate input power (accounting for efficiency)
        double pinAvg = get_output_power() / get_efficiency_or_default();

        // Calculate average input current at minimum voltage
        if (vinRmsMin <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "PFC: minimum RMS input voltage must be positive");
        }
        double iinAvg = pinAvg / vinRmsMin;

        // Peak inductor current at line peak
        double iLpeak = iinAvg * std::sqrt(2);

        // Calculate ripple current based on ratio
        double deltaI = iLpeak * get_current_ripple_ratio_or_default();

        // CCM inductance formula
        double L = vinPeakMin * D / (deltaI * get_switching_frequency());

        return L;
    }

    double PowerFactorCorrection::calculate_inductance_dcm() {
        auto inputVoltage = get_input_voltage();
        double vinRmsMin  = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);

        double D      = calculate_duty_cycle(vinPeakMin, get_output_voltage());
        double pinAvg = get_output_power() / get_efficiency_or_default();

        double L = pow(vinPeakMin, 2) * pow(D, 2) / (2 * pinAvg * get_switching_frequency());
        return L;
    }

    double PowerFactorCorrection::calculate_inductance_crcm() {
        auto inputVoltage = get_input_voltage();
        double vinRmsMin  = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);

        double D      = calculate_duty_cycle(vinPeakMin, get_output_voltage());
        double pinAvg = get_output_power() / get_efficiency_or_default();
        double iinAvg = pinAvg / vinRmsMin;
        double iLpeak = iinAvg * std::sqrt(2) * 2;  // Peak current is 2x average in CrCM

        double L = vinPeakMin * D / (iLpeak * get_switching_frequency());
        return L;
    }

    double PowerFactorCorrection::calculate_peak_current(double vinPeak, double inductance) {
        double D = calculate_duty_cycle(vinPeak, get_output_voltage());

        double pinAvg = get_output_power() / get_efficiency_or_default();
        double iinRms = pinAvg / (vinPeak / std::sqrt(2));
        double iAvg = iinRms * std::sqrt(2);

        double deltaI = vinPeak * D / (inductance * get_switching_frequency());
        return iAvg + deltaI / 2;
    }

    std::string PowerFactorCorrection::determine_actual_mode(double inductance) {
        double L_critical = calculate_inductance_crcm();
        const double tolerance = 0.05;

        if (inductance > L_critical * (1.0 + tolerance)) {
            return "Continuous Conduction Mode";
        } else if (inductance < L_critical * (1.0 - tolerance)) {
            return "Discontinuous Conduction Mode";
        } else {
            return "Critical Conduction Mode";
        }
    }

    DesignRequirements PowerFactorCorrection::process_design_requirements() {
        DesignRequirements designRequirements;

        // PFC inductor is a single winding - no turns ratio
        designRequirements.get_mutable_turns_ratios().clear();

        // Single winding
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(IsolationSide::PRIMARY);
        designRequirements.set_isolation_sides(isolationSides);

        // Set application - PFC uses POWER application for low-loss materials
        designRequirements.set_application(Application::POWER);
        // Tag the topology so downstream advisers (CoreAdviser, MagneticFilter,
        // CoilAdviser) can distinguish PFC from generic boost/buck inductors and
        // apply the correct routing — same as CMC/DMC.
        designRequirements.set_topology(Topologies::POWER_FACTOR_CORRECTION);

        // Calculate required inductance based on mode
        std::string mode = get_mode_string();
        double inductance;
        if (mode == "Continuous Conduction Mode") {
            inductance = calculate_inductance_ccm();
        } else if (mode == "Discontinuous Conduction Mode") {
            inductance = calculate_inductance_dcm();
        } else if (mode == "Critical Conduction Mode" || mode == "Transition Mode") {
            inductance = calculate_inductance_crcm();
        } else {
            inductance = calculate_inductance_ccm();
        }

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(inductance);
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        return designRequirements;
    }

    std::vector<OperatingPoint> PowerFactorCorrection::process_operating_points(const std::vector<double>& /*turnsRatios*/, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;

        std::string mode = get_mode_string();

        // Calculate required inductance if not provided
        double L = magnetizingInductance;
        if (L <= 0) {
            if (mode == "Continuous Conduction Mode") {
                L = calculate_inductance_ccm();
            } else if (mode == "Discontinuous Conduction Mode") {
                L = calculate_inductance_dcm();
            } else {
                L = calculate_inductance_crcm();
            }
        }

        auto inputVoltage = get_input_voltage();
        double vinRmsMin  = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);

        double outputVoltage     = get_output_voltage();
        double outputPower       = get_output_power();
        double switchingFrequency= get_switching_frequency();
        double lineFrequency     = get_line_frequency_or_default();
        double diodeVoltageDrop  = get_diode_voltage_drop_or_default();
        double efficiency        = get_efficiency_or_default();
        double ambientTemperature= get_ambient_temperature();

        double pinAvg     = outputPower / efficiency;
        double iinRmsAvg  = pinAvg / vinRmsMin;
        double iLinePeak  = iinRmsAvg * std::sqrt(2);

        double mainsPeriod     = 1.0 / lineFrequency;
        double switchingPeriod = 1.0 / switchingFrequency;

        int actualPeriods = (_numberOfPeriods > 0) ? _numberOfPeriods : 2;
        double totalTime  = mainsPeriod * actualPeriods;
        double timeStep   = switchingPeriod / 4.0;
        size_t numPoints  = static_cast<size_t>(totalTime / timeStep) + 1;

        std::vector<double> currentData;
        std::vector<double> voltageData;
        std::vector<double> timeData;
        currentData.reserve(numPoints);
        voltageData.reserve(numPoints);
        timeData.reserve(numPoints);

        for (size_t i = 0; i < numPoints; ++i) {
            double t = i * timeStep;
            timeData.push_back(t);

            double theta = 2.0 * std::numbers::pi * t / mainsPeriod;

            double vinInst = vinPeakMin * std::abs(std::sin(theta));
            if (vinInst < vinPeakMin * 0.05) {
                vinInst = vinPeakMin * 0.05;
            }

            double D = 1.0 - vinInst / (outputVoltage + diodeVoltageDrop);
            if (D < 0.01) D = 0.01;
            if (D > 0.95) D = 0.95;

            double iAvgInst = iLinePeak * std::abs(std::sin(theta));
            double deltaI   = vinInst * D / (L * switchingFrequency);

            double tInSwitchCycle = std::fmod(t, switchingPeriod);
            double switchPhase    = tInSwitchCycle / switchingPeriod;

            double ripple;
            if (switchPhase < D) {
                ripple = -deltaI/2 + deltaI * (switchPhase / D);
            } else {
                ripple = deltaI/2 - deltaI * ((switchPhase - D) / (1 - D));
            }

            currentData.push_back(iAvgInst + ripple);

            if (switchPhase < D) {
                voltageData.push_back(vinInst);
            } else {
                voltageData.push_back(vinInst - outputVoltage - diodeVoltageDrop);
            }
        }

        Waveform currentWaveform;
        currentWaveform.set_data(currentData);
        currentWaveform.set_time(timeData);

        SignalDescriptor current;
        current.set_waveform(currentWaveform);
        auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, lineFrequency);
        current.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform, lineFrequency));
        auto currentProcessed = Inputs::calculate_processed_data(currentWaveform, lineFrequency, true);
        currentProcessed.set_label(WaveformLabel::CUSTOM);
        current.set_processed(currentProcessed);

        Waveform voltageWaveform;
        voltageWaveform.set_data(voltageData);
        voltageWaveform.set_time(timeData);

        SignalDescriptor voltage;
        voltage.set_waveform(voltageWaveform);
        auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, lineFrequency);
        voltage.set_harmonics(Inputs::calculate_harmonics_data(sampledVoltageWaveform, lineFrequency));
        auto voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, lineFrequency, true);
        voltageProcessed.set_label(WaveformLabel::CUSTOM);
        voltage.set_processed(voltageProcessed);

        OperatingPointExcitation excitation;
        excitation.set_current(current);
        excitation.set_frequency(lineFrequency);
        excitation.set_voltage(voltage);

        OperatingPoint operatingPoint;
        operatingPoint.set_excitations_per_winding({excitation});
        operatingPoint.get_mutable_conditions().set_ambient_temperature(ambientTemperature);
        operatingPoint.set_name("PFC_Line_Period_" + std::to_string(actualPeriods) + "_Periods");

        operatingPoints.push_back(operatingPoint);

        return operatingPoints;
    }

    std::string PowerFactorCorrection::generate_ngspice_circuit(double inductance,
                                                                 double /*dcResistance*/,
                                                                 double simulationTime,
                                                                 double timeStep) {
        std::ostringstream netlist;

        netlist << "* PFC Boost Converter - Ideal Behavioral Model\n";
        netlist << "* Generated by OpenMagnetics\n";
        netlist << "* Models ideal PFC with sinusoidal current envelope + switching ripple\n\n";

        auto inputVoltage = get_input_voltage();
        double vinRms;
        if (inputVoltage.get_nominal().has_value()) {
            vinRms = inputVoltage.get_nominal().value();
        } else {
            vinRms = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        }
        double vinPeak = vinRms * std::sqrt(2);
        double vout    = get_output_voltage();
        double pout    = get_output_power();
        double switchingFrequency = get_switching_frequency();
        double lineFrequency      = get_line_frequency_or_default();

        double iPeak = std::sqrt(2) * pout / vinRms;
        double tSw   = 1.0 / switchingFrequency;

        netlist << ".param vin_peak=" << vinPeak << "\n";
        netlist << ".param vout=" << vout << "\n";
        netlist << ".param fline=" << lineFrequency << "\n";
        netlist << ".param fsw=" << switchingFrequency << "\n";
        netlist << ".param L=" << inductance << "\n";
        netlist << ".param i_peak=" << iPeak << "\n\n";

        netlist << "* Rectified AC Input\n";
        netlist << "B_vin vin_rect 0 V=vin_peak*abs(sin(2*3.14159265*fline*time))\n\n";

        netlist << "* Ideal current envelope (sinusoidal, in phase with voltage)\n";
        netlist << "B_ienv i_env 0 V=i_peak*abs(sin(2*3.14159265*fline*time))\n\n";

        netlist << "* Instantaneous duty cycle\n";
        netlist << "B_duty duty 0 V=1-V(vin_rect)/vout\n\n";

        netlist << "* Current ripple amplitude\n";
        netlist << "B_rip ripple 0 V=V(vin_rect)*V(duty)/(L*fsw)/2\n\n";

        netlist << "* Sawtooth for triangular switching ripple\n";
        netlist << "V_saw saw 0 PULSE(-1 1 0 " << tSw/2 << " " << tSw/2 << " 1n " << tSw << ")\n\n";

        netlist << "* Total inductor current (envelope + ripple)\n";
        netlist << "B_iL i_L 0 V=V(i_env)+V(ripple)*V(saw)\n\n";

        netlist << "* Analysis\n";
        netlist << ".tran " << timeStep << " " << simulationTime << " 0 " << timeStep << "\n";
        netlist << ".save v(vin_rect) v(i_env) v(i_L) v(ripple)\n";
        netlist << ".end\n";

        return netlist.str();
    }

    PfcSimulationWaveforms PowerFactorCorrection::simulate_and_extract_waveforms(double inductance,
                                                                                   double /*dcResistance*/,
                                                                                   int numberOfCycles) {
        PfcSimulationWaveforms waveforms;
        double switchingFrequency = get_switching_frequency();
        double lineFrequency      = get_line_frequency_or_default();
        waveforms.switchingFrequency = switchingFrequency;
        waveforms.lineFrequency      = lineFrequency;

        auto inputVoltage = get_input_voltage();
        double vinRms;
        if (inputVoltage.get_nominal().has_value()) {
            vinRms = inputVoltage.get_nominal().value();
        } else {
            vinRms = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        }
        double vinPeak = vinRms * std::sqrt(2);
        double vout    = get_output_voltage();
        double pout    = get_output_power();

        double iPeak = std::sqrt(2) * pout / vinRms;

        double linePeriod      = 1.0 / lineFrequency;
        double switchingPeriod = 1.0 / switchingFrequency;
        double simulationTime  = numberOfCycles * linePeriod;

        double timeStep   = switchingPeriod / 4.0;
        size_t numPoints  = static_cast<size_t>(simulationTime / timeStep) + 1;

        waveforms.time.resize(numPoints);
        waveforms.inputVoltage.resize(numPoints);
        waveforms.inductorCurrent.resize(numPoints);
        waveforms.inductorVoltage.resize(numPoints);
        waveforms.currentEnvelope.resize(numPoints);
        waveforms.currentRipple.resize(numPoints);

        double omega_line = 2.0 * M_PI * lineFrequency;

        for (size_t i = 0; i < numPoints; i++) {
            double t = i * timeStep;
            waveforms.time[i] = t;

            double vin = vinPeak * std::abs(std::sin(omega_line * t));
            waveforms.inputVoltage[i] = vin;

            double iEnv = iPeak * std::abs(std::sin(omega_line * t));
            waveforms.currentEnvelope[i] = iEnv;

            double duty = 1.0 - vin / vout;
            if (duty < 0) duty = 0;
            if (duty > 1) duty = 1;

            double rippleAmplitude = vin * duty / (2.0 * inductance * switchingFrequency);
            waveforms.currentRipple[i] = rippleAmplitude;

            double tInSwitchingCycle = std::fmod(t, switchingPeriod);
            double sawtoothPhase = tInSwitchingCycle / switchingPeriod;
            double triangular;
            if (sawtoothPhase < 0.5) {
                triangular = -1.0 + 4.0 * sawtoothPhase;
            } else {
                triangular = 3.0 - 4.0 * sawtoothPhase;
            }

            waveforms.inductorCurrent[i] = iEnv + rippleAmplitude * triangular;

            if (sawtoothPhase < 0.5) {
                waveforms.inductorVoltage[i] = vin;
            } else {
                waveforms.inductorVoltage[i] = vin - vout;
            }
        }

        waveforms.inputCurrent = waveforms.inductorCurrent;

        double realPower = 0, vRms = 0, iRms = 0;
        for (size_t i = 0; i < numPoints; i++) {
            realPower += waveforms.inputVoltage[i] * waveforms.inputCurrent[i];
            vRms      += waveforms.inputVoltage[i] * waveforms.inputVoltage[i];
            iRms      += waveforms.inputCurrent[i] * waveforms.inputCurrent[i];
        }
        realPower /= numPoints;
        vRms = std::sqrt(vRms / numPoints);
        iRms = std::sqrt(iRms / numPoints);

        double apparentPower = vRms * iRms;
        waveforms.powerFactor = (apparentPower > 0) ? realPower / apparentPower : 1.0;
        waveforms.efficiency  = 1.0;
        waveforms.operatingPointName = "PFC_analytical";

        return waveforms;
    }

    std::vector<OperatingPoint> PowerFactorCorrection::simulate_and_extract_operating_points(double inductance,
                                                                                               double dcResistance) {
        std::vector<OperatingPoint> operatingPoints;

        int numPeriods = (_numberOfPeriods > 0) ? _numberOfPeriods : 2;
        PfcSimulationWaveforms waveforms = simulate_and_extract_waveforms(inductance, dcResistance, numPeriods);

        if (waveforms.inductorCurrent.empty() || waveforms.time.empty()) {
            return process_operating_points({}, inductance);
        }

        double switchingFrequency = get_switching_frequency();
        double lineFrequency      = get_line_frequency_or_default();
        double ambientTemperature = get_ambient_temperature();

        Waveform currentWaveform;
        currentWaveform.set_time(waveforms.time);
        currentWaveform.set_data(waveforms.inductorCurrent);

        SignalDescriptor current;
        current.set_waveform(currentWaveform);

        try {
            auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
            current.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency));
            auto currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true);
            currentProcessed.set_label(WaveformLabel::CUSTOM);
            current.set_processed(currentProcessed);
        } catch (const std::exception&) {
            Processed processed;
            if (!waveforms.inductorCurrent.empty()) {
                double iMax = *std::max_element(waveforms.inductorCurrent.begin(), waveforms.inductorCurrent.end());
                double iMin = *std::min_element(waveforms.inductorCurrent.begin(), waveforms.inductorCurrent.end());
                double iAvg = std::accumulate(waveforms.inductorCurrent.begin(), waveforms.inductorCurrent.end(), 0.0) /
                              waveforms.inductorCurrent.size();
                processed.set_peak(iMax);
                processed.set_peak_to_peak(iMax - iMin);
                processed.set_average(iAvg);

                double sumSq = 0;
                for (double val : waveforms.inductorCurrent) {
                    sumSq += val * val;
                }
                processed.set_rms(std::sqrt(sumSq / waveforms.inductorCurrent.size()));
            }
            current.set_processed(processed);
        }

        SignalDescriptor voltage;
        if (!waveforms.inductorVoltage.empty() && waveforms.inductorVoltage.size() == waveforms.time.size()) {
            Waveform voltageWaveform;
            voltageWaveform.set_time(waveforms.time);
            voltageWaveform.set_data(waveforms.inductorVoltage);
            voltage.set_waveform(voltageWaveform);

            try {
                auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
                voltage.set_harmonics(Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency));
                auto voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true);
                voltageProcessed.set_label(WaveformLabel::CUSTOM);
                voltage.set_processed(voltageProcessed);
            } catch (const std::exception&) {
                // Just use raw waveform without harmonics
            }
        }

        OperatingPointExcitation excitation;
        excitation.set_current(current);
        excitation.set_frequency(lineFrequency);
        excitation.set_voltage(voltage);

        OperatingPoint operatingPoint;
        operatingPoint.set_excitations_per_winding({excitation});
        operatingPoint.get_mutable_conditions().set_ambient_temperature(ambientTemperature);
        operatingPoint.set_name("PFC_simulated_" + std::to_string(static_cast<int>(lineFrequency)) + "Hz");

        operatingPoints.push_back(operatingPoint);

        return operatingPoints;
    }

} // namespace OpenMagnetics
