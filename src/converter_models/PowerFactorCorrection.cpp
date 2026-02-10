#include "converter_models/PowerFactorCorrection.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <numbers>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    PowerFactorCorrection::PowerFactorCorrection(const json& j) : Topology(j) {
        // Parse input voltage (AC RMS)
        if (j.contains("inputVoltage")) {
            _inputVoltage = j["inputVoltage"].get<DimensionWithTolerance>();
        }
        
        // Parse output voltage (DC)
        if (j.contains("outputVoltage")) {
            _outputVoltage = j["outputVoltage"].get<double>();
        }
        
        // Parse output power
        if (j.contains("outputPower")) {
            _outputPower = j["outputPower"].get<double>();
        }
        
        // Parse switching frequency
        if (j.contains("switchingFrequency")) {
            _switchingFrequency = j["switchingFrequency"].get<double>();
        }
        
        // Parse optional parameters with defaults
        if (j.contains("lineFrequency")) {
            _lineFrequency = j["lineFrequency"].get<double>();
        }
        if (j.contains("currentRippleRatio")) {
            _currentRippleRatio = j["currentRippleRatio"].get<double>();
        }
        if (j.contains("efficiency")) {
            _efficiency = j["efficiency"].get<double>();
        }
        if (j.contains("mode")) {
            _mode = j["mode"].get<std::string>();
        }
        if (j.contains("diodeVoltageDrop")) {
            _diodeVoltageDrop = j["diodeVoltageDrop"].get<double>();
        }
        if (j.contains("ambientTemperature")) {
            _ambientTemperature = j["ambientTemperature"].get<double>();
        }
    }

    bool PowerFactorCorrection::run_checks(bool assert) {
        bool valid = true;

        // Check that output voltage is greater than peak input voltage
        double vinPeakMax = resolve_dimensional_values(_inputVoltage, DimensionalValues::MAXIMUM) * std::sqrt(2);
        if (_outputVoltage <= vinPeakMax) {
            if (assert) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, 
                    "PFC output voltage must be greater than peak input voltage. Vout: " + 
                    std::to_string(_outputVoltage) + " <= Vin_peak_max: " + std::to_string(vinPeakMax));
            }
            valid = false;
        }

        // Check efficiency is reasonable
        if (_efficiency <= 0 || _efficiency > 1) {
            if (assert) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, 
                    "Efficiency must be between 0 and 1. Got: " + std::to_string(_efficiency));
            }
            valid = false;
        }

        return valid;
    }

    double PowerFactorCorrection::calculate_duty_cycle(double vinPeak, double vout) {
        // For boost PFC: D = 1 - Vin/Vout
        // Accounting for diode drop: D = 1 - Vin/(Vout + Vd)
        return 1.0 - vinPeak / (_outputVoltage + _diodeVoltageDrop);
    }

    double PowerFactorCorrection::calculate_inductance_ccm() {
        // For CCM PFC, worst case (maximum inductance requirement) is at minimum input voltage
        // at the peak of the AC line (where current is maximum)
        
        double vinRmsMin = resolve_dimensional_values(_inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);
        
        // Calculate duty cycle at minimum input voltage peak
        double D = calculate_duty_cycle(vinPeakMin, _outputVoltage);
        
        // Calculate input power (accounting for efficiency)
        double pinAvg = _outputPower / _efficiency;
        
        // Calculate average input current at minimum voltage
        double iinAvg = pinAvg / vinRmsMin;
        
        // Peak inductor current at line peak
        double iLpeak = iinAvg * std::sqrt(2);
        
        // Calculate ripple current based on ratio
        double deltaI = iLpeak * _currentRippleRatio;
        
        // CCM inductance formula: L = Vin * D / (deltaI * fsw)
        // At the peak of the line: Vin_peak * D / (deltaI * fsw)
        double L = vinPeakMin * D / (deltaI * _switchingFrequency);
        
        return L;
    }

    double PowerFactorCorrection::calculate_inductance_dcm() {
        // For DCM, inductance determines the power throughput
        // L = Vin^2 * D^2 / (2 * P * fsw)
        
        double vinRmsMin = resolve_dimensional_values(_inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);
        
        double D = calculate_duty_cycle(vinPeakMin, _outputVoltage);
        double pinAvg = _outputPower / _efficiency;
        
        double L = pow(vinPeakMin, 2) * pow(D, 2) / (2 * pinAvg * _switchingFrequency);
        
        return L;
    }

    double PowerFactorCorrection::calculate_inductance_crcm() {
        // For CrCM/TCM, inductance at boundary between CCM and DCM
        // This results in the current just reaching zero at the end of each cycle
        
        double vinRmsMin = resolve_dimensional_values(_inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);
        
        double D = calculate_duty_cycle(vinPeakMin, _outputVoltage);
        double pinAvg = _outputPower / _efficiency;
        double iinAvg = pinAvg / vinRmsMin;
        double iLpeak = iinAvg * std::sqrt(2) * 2;  // Peak current is 2x average in CrCM
        
        // At boundary: deltaI = 2 * I_avg, so ripple ratio = 2
        double L = vinPeakMin * D / (iLpeak * _switchingFrequency);
        
        return L;
    }

    double PowerFactorCorrection::calculate_peak_current(double vinPeak, double inductance) {
        // Calculate peak current at a given input voltage
        double D = calculate_duty_cycle(vinPeak, _outputVoltage);
        
        // Average input power
        double pinAvg = _outputPower / _efficiency;
        
        // Average current over the AC line cycle at this instantaneous voltage
        // For sinusoidal input: i(θ) = Ipk * sin(θ)
        // At line peak (θ = 90°): i = Ipk
        double iinRms = pinAvg / (vinPeak / std::sqrt(2));
        double iAvg = iinRms * std::sqrt(2);  // Peak of average current envelope
        
        // Ripple current: deltaI = Vin * D / (L * fsw)
        double deltaI = vinPeak * D / (inductance * _switchingFrequency);
        
        // Peak current = average + ripple/2
        return iAvg + deltaI / 2;
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
        // No specific sub-application for PFC boost inductor

        // Calculate required inductance based on mode
        double inductance;
        if (_mode == "Continuous Conduction Mode") {
            inductance = calculate_inductance_ccm();
        } else if (_mode == "Discontinuous Conduction Mode") {
            inductance = calculate_inductance_dcm();
        } else if (_mode == "Critical Conduction Mode" || _mode == "Transition Mode") {
            inductance = calculate_inductance_crcm();
        } else {
            inductance = calculate_inductance_ccm();  // Default to CCM
        }

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(inductance);
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        return designRequirements;
    }

    std::vector<OperatingPoint> PowerFactorCorrection::process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;

        // Calculate required inductance if not provided
        double L = magnetizingInductance;
        if (L <= 0) {
            if (_mode == "Continuous Conduction Mode") {
                L = calculate_inductance_ccm();
            } else if (_mode == "Discontinuous Conduction Mode") {
                L = calculate_inductance_dcm();
            } else {
                L = calculate_inductance_crcm();
            }
        }

        // Use minimum input voltage as worst case for current
        double vinRmsMin = resolve_dimensional_values(_inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);

        // ============================================================
        // PFC Waveform Generation: One complete mains cycle
        // ============================================================
        // The PFC inductor current follows a rectified sinusoidal envelope
        // (the line current shape) with high-frequency triangular ripple
        // superimposed at the switching frequency.
        //
        // For one half-cycle of the mains (e.g., 10ms for 50Hz):
        // - Envelope: I_avg(θ) = I_peak * |sin(θ)|
        // - Ripple: ΔI(θ) = V_in(θ) * D(θ) / (L * f_sw)
        // - D(θ) = 1 - V_in(θ) / V_out for boost PFC
        //
        // We simulate one complete mains cycle (full rectified sine)
        // which contains: f_sw / f_line switching cycles
        // ============================================================

        double pinAvg = _outputPower / _efficiency;
        double iinRmsAvg = pinAvg / vinRmsMin;
        double iLinePeak = iinRmsAvg * std::sqrt(2);  // Peak of the line current envelope

        // Number of switching cycles in one mains half-period
        double mainsHalfPeriod = 1.0 / (2.0 * _lineFrequency);
        int switchingCyclesPerHalfPeriod = static_cast<int>(std::round(_switchingFrequency * mainsHalfPeriod));
        
        // Use power-of-2 points per switching cycle for clean FFT
        int pointsPerSwitchingCycle = 32;
        int totalPoints = switchingCyclesPerHalfPeriod * pointsPerSwitchingCycle;
        
        // Ensure power of 2 total points for FFT
        int powerOf2 = 1;
        while (powerOf2 < totalPoints) powerOf2 *= 2;
        totalPoints = powerOf2;

        double switchingPeriod = 1.0 / _switchingFrequency;
        double totalTime = mainsHalfPeriod;  // One half of AC line cycle
        double dt = totalTime / totalPoints;

        std::vector<double> currentData(totalPoints);
        std::vector<double> voltageData(totalPoints);
        std::vector<double> timeData(totalPoints);

        for (int i = 0; i < totalPoints; ++i) {
            double t = i * dt;
            timeData[i] = t;

            // Line phase angle (0 to π for half cycle)
            double theta = std::numbers::pi * t / mainsHalfPeriod;
            
            // Instantaneous rectified input voltage
            double vinInst = vinPeakMin * std::abs(std::sin(theta));
            
            // Protect against division by zero near zero crossings
            if (vinInst < vinPeakMin * 0.05) {
                vinInst = vinPeakMin * 0.05;
            }
            
            // Duty cycle varies with instantaneous input voltage
            // D = 1 - Vin/Vout for boost
            double D = 1.0 - vinInst / (_outputVoltage + _diodeVoltageDrop);
            if (D < 0.01) D = 0.01;
            if (D > 0.95) D = 0.95;
            
            // Average current at this phase angle follows sinusoidal envelope
            double iAvgInst = iLinePeak * std::abs(std::sin(theta));
            
            // Ripple current at this point: ΔI = Vin * D / (L * fsw)
            double deltaI = vinInst * D / (L * _switchingFrequency);
            
            // Position within current switching cycle
            double tInSwitchCycle = std::fmod(t, switchingPeriod);
            double switchPhase = tInSwitchCycle / switchingPeriod;
            
            // Triangular ripple waveform
            // During on-time (0 to D): current ramps up
            // During off-time (D to 1): current ramps down
            double ripple;
            if (switchPhase < D) {
                // On-time: ramp up from -deltaI/2 to +deltaI/2
                ripple = -deltaI/2 + deltaI * (switchPhase / D);
            } else {
                // Off-time: ramp down from +deltaI/2 to -deltaI/2
                ripple = deltaI/2 - deltaI * ((switchPhase - D) / (1 - D));
            }
            
            // Total current = envelope + ripple
            currentData[i] = iAvgInst + ripple;
            
            // Voltage across inductor:
            // On-time: V_L = Vin (switch closed, diode reverse biased)
            // Off-time: V_L = Vin - Vout (switch open, diode conducting)
            if (switchPhase < D) {
                voltageData[i] = vinInst;  // On-time
            } else {
                voltageData[i] = vinInst - _outputVoltage - _diodeVoltageDrop;  // Off-time (negative)
            }
        }

        // Create current waveform
        Waveform currentWaveform;
        currentWaveform.set_data(currentData);
        currentWaveform.set_time(timeData);

        SignalDescriptor current;
        current.set_waveform(currentWaveform);
        auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, _switchingFrequency);
        current.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform, _switchingFrequency));
        current.set_processed(Inputs::calculate_processed_data(currentWaveform, _switchingFrequency, true));

        // Create voltage waveform
        Waveform voltageWaveform;
        voltageWaveform.set_data(voltageData);
        voltageWaveform.set_time(timeData);

        SignalDescriptor voltage;
        voltage.set_waveform(voltageWaveform);
        auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, _switchingFrequency);
        voltage.set_harmonics(Inputs::calculate_harmonics_data(sampledVoltageWaveform, _switchingFrequency));
        voltage.set_processed(Inputs::calculate_processed_data(voltageWaveform, _switchingFrequency, true));

        // Create operating point
        OperatingPointExcitation excitation;
        excitation.set_current(current);
        excitation.set_frequency(_switchingFrequency);
        excitation.set_voltage(voltage);

        OperatingPoint operatingPoint;
        operatingPoint.set_excitations_per_winding({excitation});
        operatingPoint.get_mutable_conditions().set_ambient_temperature(_ambientTemperature);
        operatingPoint.set_name("Vin_min_half_cycle");

        operatingPoints.push_back(operatingPoint);

        return operatingPoints;
    }

    std::string PowerFactorCorrection::generate_ngspice_circuit(double inductance,
                                                                 double dcResistance,
                                                                 double simulationTime,
                                                                 double timeStep) {
        std::ostringstream netlist;
        
        // Header
        netlist << "* PFC Boost Converter - Ideal Behavioral Model\n";
        netlist << "* Generated by OpenMagnetics\n";
        netlist << "* Models ideal PFC with sinusoidal current envelope + switching ripple\n\n";
        
        // Get operating parameters
        double vinPeak = _inputVoltage.get_nominal().value() * std::sqrt(2);
        double vinRms = _inputVoltage.get_nominal().value();
        double vout = _outputVoltage;
        double pout = _outputPower;
        
        // Peak current for unity power factor: Ipeak = sqrt(2) * Pout / Vin_rms
        double iPeak = std::sqrt(2) * pout / vinRms;
        
        // Switching period
        double tSw = 1.0 / _switchingFrequency;
        
        netlist << ".param vin_peak=" << vinPeak << "\n";
        netlist << ".param vout=" << vout << "\n";
        netlist << ".param fline=" << _lineFrequency << "\n";
        netlist << ".param fsw=" << _switchingFrequency << "\n";
        netlist << ".param L=" << inductance << "\n";
        netlist << ".param i_peak=" << iPeak << "\n\n";
        
        // Rectified AC input voltage (full-wave rectified sine)
        netlist << "* Rectified AC Input\n";
        netlist << "B_vin vin_rect 0 V=vin_peak*abs(sin(2*3.14159265*fline*time))\n\n";
        
        // Ideal sinusoidal current envelope (in phase with voltage for unity PF)
        netlist << "* Ideal current envelope (sinusoidal, in phase with voltage)\n";
        netlist << "B_ienv i_env 0 V=i_peak*abs(sin(2*3.14159265*fline*time))\n\n";
        
        // Duty cycle: D = 1 - Vin/Vout (for CCM boost)
        netlist << "* Instantaneous duty cycle\n";
        netlist << "B_duty duty 0 V=1-V(vin_rect)/vout\n\n";
        
        // Ripple amplitude: dI = Vin * D / (L * fsw)
        netlist << "* Current ripple amplitude\n";
        netlist << "B_rip ripple 0 V=V(vin_rect)*V(duty)/(L*fsw)/2\n\n";
        
        // Sawtooth for triangular ripple
        netlist << "* Sawtooth for triangular switching ripple\n";
        netlist << "V_saw saw 0 PULSE(-1 1 0 " << tSw/2 << " " << tSw/2 << " 1n " << tSw << ")\n\n";
        
        // Total inductor current = sinusoidal envelope + triangular ripple
        netlist << "* Total inductor current (envelope + ripple)\n";
        netlist << "B_iL i_L 0 V=V(i_env)+V(ripple)*V(saw)\n\n";
        
        // Simulation commands
        netlist << "* Analysis\n";
        netlist << ".tran " << timeStep << " " << simulationTime << " 0 " << timeStep << "\n";
        netlist << ".save v(vin_rect) v(i_env) v(i_L) v(ripple)\n";
        netlist << ".end\n";
        
        return netlist.str();
    }

    PfcSimulationWaveforms PowerFactorCorrection::simulate_and_extract_waveforms(double inductance,
                                                                                   double dcResistance,
                                                                                   int numberOfCycles) {
        PfcSimulationWaveforms waveforms;
        waveforms.switchingFrequency = _switchingFrequency;
        waveforms.lineFrequency = _lineFrequency;
        
        // Operating parameters
        double vinPeak = _inputVoltage.get_nominal().value() * std::sqrt(2);
        double vinRms = _inputVoltage.get_nominal().value();
        double vout = _outputVoltage;
        double pout = _outputPower;
        
        // Peak current for unity power factor: Ipeak = sqrt(2) * Pout / Vin_rms
        double iPeak = std::sqrt(2) * pout / vinRms;
        
        // Time parameters
        double linePeriod = 1.0 / _lineFrequency;
        double switchingPeriod = 1.0 / _switchingFrequency;
        double simulationTime = numberOfCycles * linePeriod;
        
        // Use 100 points per switching period for good resolution
        double timeStep = switchingPeriod / 100.0;
        size_t numPoints = static_cast<size_t>(simulationTime / timeStep) + 1;
        
        // Pre-allocate vectors
        waveforms.time.resize(numPoints);
        waveforms.inputVoltage.resize(numPoints);
        waveforms.inductorCurrent.resize(numPoints);
        waveforms.currentEnvelope.resize(numPoints);
        waveforms.currentRipple.resize(numPoints);
        
        // Generate waveforms analytically
        double omega_line = 2.0 * M_PI * _lineFrequency;
        
        for (size_t i = 0; i < numPoints; i++) {
            double t = i * timeStep;
            waveforms.time[i] = t;
            
            // Rectified input voltage: |Vin_peak * sin(ωt)|
            double vin = vinPeak * std::abs(std::sin(omega_line * t));
            waveforms.inputVoltage[i] = vin;
            
            // Sinusoidal current envelope (in phase with voltage for unity PF)
            double iEnv = iPeak * std::abs(std::sin(omega_line * t));
            waveforms.currentEnvelope[i] = iEnv;
            
            // Duty cycle: D = 1 - Vin/Vout (for CCM boost)
            double duty = 1.0 - vin / vout;
            if (duty < 0) duty = 0;
            if (duty > 1) duty = 1;
            
            // Current ripple amplitude: ΔI/2 = Vin * D / (2 * L * fsw)
            double rippleAmplitude = vin * duty / (2.0 * inductance * _switchingFrequency);
            waveforms.currentRipple[i] = rippleAmplitude;
            
            // Triangular ripple: sawtooth from -1 to +1 over switching period
            double tInSwitchingCycle = std::fmod(t, switchingPeriod);
            double sawtoothPhase = tInSwitchingCycle / switchingPeriod;  // 0 to 1
            double triangular;
            if (sawtoothPhase < 0.5) {
                triangular = -1.0 + 4.0 * sawtoothPhase;  // -1 to +1
            } else {
                triangular = 3.0 - 4.0 * sawtoothPhase;   // +1 to -1
            }
            
            // Total inductor current = envelope + ripple
            waveforms.inductorCurrent[i] = iEnv + rippleAmplitude * triangular;
        }
        
        // Input current is same as inductor current for boost PFC
        waveforms.inputCurrent = waveforms.inductorCurrent;
        
        // Calculate power factor (should be ~1.0 for ideal sinusoidal current)
        double realPower = 0, vRms = 0, iRms = 0;
        for (size_t i = 0; i < numPoints; i++) {
            realPower += waveforms.inputVoltage[i] * waveforms.inputCurrent[i];
            vRms += waveforms.inputVoltage[i] * waveforms.inputVoltage[i];
            iRms += waveforms.inputCurrent[i] * waveforms.inputCurrent[i];
        }
        realPower /= numPoints;
        vRms = std::sqrt(vRms / numPoints);
        iRms = std::sqrt(iRms / numPoints);
        
        double apparentPower = vRms * iRms;
        waveforms.powerFactor = (apparentPower > 0) ? realPower / apparentPower : 1.0;
        
        // Ideal model: 100% efficiency
        waveforms.efficiency = 1.0;
        
        waveforms.operatingPointName = "PFC_analytical";
        
        return waveforms;
    }

    std::vector<OperatingPoint> PowerFactorCorrection::simulate_and_extract_operating_points(double inductance,
                                                                                               double dcResistance) {
        std::vector<OperatingPoint> operatingPoints;
        
        // Run simulation for one complete line cycle
        PfcSimulationWaveforms waveforms = simulate_and_extract_waveforms(inductance, dcResistance, 1);
        
        // If simulation didn't return valid waveform data, fall back to analytical operating points
        if (waveforms.inductorCurrent.empty() || waveforms.time.empty()) {
            // Fall back to analytical method
            return process_operating_points({}, inductance);
        }
        
        // Create Waveform objects from simulation data
        Waveform currentWaveform;
        currentWaveform.set_time(waveforms.time);
        currentWaveform.set_data(waveforms.inductorCurrent);
        
        SignalDescriptor current;
        current.set_waveform(currentWaveform);
        
        // Only try to sample if we have valid data
        try {
            auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, _switchingFrequency);
            current.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform, _switchingFrequency));
            current.set_processed(Inputs::calculate_processed_data(currentWaveform, _switchingFrequency, true));
        } catch (const std::exception& e) {
            // If sampling fails, just use the raw waveform without harmonics
            Processed processed;
            if (!waveforms.inductorCurrent.empty()) {
                double iMax = *std::max_element(waveforms.inductorCurrent.begin(), waveforms.inductorCurrent.end());
                double iMin = *std::min_element(waveforms.inductorCurrent.begin(), waveforms.inductorCurrent.end());
                double iAvg = std::accumulate(waveforms.inductorCurrent.begin(), waveforms.inductorCurrent.end(), 0.0) / 
                              waveforms.inductorCurrent.size();
                processed.set_peak(iMax);
                processed.set_peak_to_peak(iMax - iMin);
                processed.set_average(iAvg);
                
                // Calculate RMS
                double sumSq = 0;
                for (double val : waveforms.inductorCurrent) {
                    sumSq += val * val;
                }
                processed.set_rms(std::sqrt(sumSq / waveforms.inductorCurrent.size()));
            }
            current.set_processed(processed);
        }
        
        // Create voltage waveform
        SignalDescriptor voltage;
        if (!waveforms.inductorVoltage.empty() && waveforms.inductorVoltage.size() == waveforms.time.size()) {
            Waveform voltageWaveform;
            voltageWaveform.set_time(waveforms.time);
            voltageWaveform.set_data(waveforms.inductorVoltage);
            voltage.set_waveform(voltageWaveform);
            
            try {
                auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, _switchingFrequency);
                voltage.set_harmonics(Inputs::calculate_harmonics_data(sampledVoltageWaveform, _switchingFrequency));
                voltage.set_processed(Inputs::calculate_processed_data(voltageWaveform, _switchingFrequency, true));
            } catch (const std::exception& e) {
                // Just use raw waveform without harmonics
            }
        }
        
        // Create operating point
        OperatingPointExcitation excitation;
        excitation.set_current(current);
        excitation.set_frequency(_switchingFrequency);
        excitation.set_voltage(voltage);
        
        OperatingPoint operatingPoint;
        operatingPoint.set_excitations_per_winding({excitation});
        operatingPoint.get_mutable_conditions().set_ambient_temperature(_ambientTemperature);
        operatingPoint.set_name("PFC_simulated_" + std::to_string(static_cast<int>(_lineFrequency)) + "Hz");
        
        operatingPoints.push_back(operatingPoint);
        
        return operatingPoints;
    }

} // namespace OpenMagnetics
