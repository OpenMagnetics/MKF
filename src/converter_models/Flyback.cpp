#include "converter_models/Flyback.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/CircuitSimulatorInterface.h"
#include "support/Utils.h"
#include <cfloat>
#include <algorithm>
#include <map>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    Flyback::Flyback(const json& j) {
        from_json(j, *this);
    }

    AdvancedFlyback::AdvancedFlyback(const json& j) {
        from_json(j, *this);
    }

    double calculate_BMO_duty_cycle(double outputVoltage, double inputVoltage, double turnsRatio){
        return (turnsRatio * outputVoltage) / (inputVoltage + turnsRatio * outputVoltage);
    }

    double calculate_BMO_primary_current_peak(double outputCurrent, double efficiency, double dutyCycle, double turnsRatio){
        return (2 * outputCurrent) / (efficiency * (1 - dutyCycle) * turnsRatio);
    }


    double calculate_QRM_frequency(double magnetizingInductance, double outputPower, double outputVoltage, double minimumInputVoltage, double turnsRatio, 
                                   double diodeVoltageDrop, double efficiency, double drainSourceCapacitance = 100e-12) {

        double dt = std::numbers::pi * sqrt(magnetizingInductance * drainSourceCapacitance);
        double a = pow((outputVoltage + diodeVoltageDrop + 1.0 / turnsRatio * minimumInputVoltage),2);
        double b = efficiency * minimumInputVoltage * minimumInputVoltage * pow((outputVoltage + diodeVoltageDrop),2);
        double c = outputVoltage + diodeVoltageDrop + 1.0 / turnsRatio * minimumInputVoltage;
        double d = sqrt(outputPower/(efficiency * magnetizingInductance));
        double e = minimumInputVoltage * (outputVoltage + diodeVoltageDrop);
        double f = sqrt(4 * dt + (2 * magnetizingInductance * outputPower* a) / b);
        double g = (1.414 * magnetizingInductance * c * d) / e;
        double h = 4.0 / pow((f + g),2);

        return h;
    }


    FlybackModes FlybackOperatingPoint::resolve_mode(std::optional<double> currentRippleRatio) {
        if (get_mode()) {
            return get_mode().value();
        }
        else {
            if (!currentRippleRatio) {
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Either current ripple ratio or mode is needed for the Flyback OperatingPoint Mode");
            }
            auto mode = currentRippleRatio.value() < 1? FlybackModes::CONTINUOUS_CONDUCTION_MODE : FlybackModes::DISCONTINUOUS_CONDUCTION_MODE;
            return mode;
        }
    }
    double FlybackOperatingPoint::resolve_switching_frequency(double inputVoltage, double diodeVoltageDrop, std::optional<double> inductance, std::optional<std::vector<double>> turnsRatios, double efficiency) {
        if (get_switching_frequency()) {
            return get_switching_frequency().value();
        }
        else {
            if (!get_mode()) {
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Either switching frequency or mode is needed for the Flyback OperatingPoint");
            }
            auto mode = get_mode().value();
            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "Switching Frequency is needed for CCM");
                }
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "Switching Frequency is needed for DCM");
                }
                case FlybackModes::QUASI_RESONANT_MODE: {
                    if (!inductance) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "Inductance in missing for switching frequency calculation");
                    }
                    if (!turnsRatios) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "TurnsRatios in missing for switching frequency calculation");
                    }
                    double totalOutputVoltageReflectedPrimaryMinusDiode = 0;

                    for (size_t secondaryIndex = 0; secondaryIndex < get_output_voltages().size(); ++secondaryIndex) {
                        auto outputVoltage = get_output_voltages()[secondaryIndex];
                        auto turnsRatio = turnsRatios.value()[secondaryIndex];
                        totalOutputVoltageReflectedPrimaryMinusDiode += outputVoltage * turnsRatio;
                    }
                    
                    double totalOutputPower = Flyback::get_total_input_power(get_output_currents(), get_output_voltages(), 1, 0);

                    double switchingFrequency = calculate_QRM_frequency(inductance.value(), totalOutputPower, totalOutputVoltageReflectedPrimaryMinusDiode / turnsRatios.value()[0], inputVoltage, turnsRatios.value()[0], diodeVoltageDrop, efficiency);
                    return switchingFrequency;
                }
                case FlybackModes::BOUNDARY_MODE_OPERATION: {
                    if (!inductance) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "Inductance in missing for switching frequency calculation");
                    }
                    if (!turnsRatios) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "TurnsRatios in missing for switching frequency calculation");
                    }
                    double currentPeak = 0;
                    double switchingFrequency = 0;
                    for (size_t secondaryIndex = 0; secondaryIndex < get_output_voltages().size(); ++secondaryIndex) {
                        auto outputCurrent = get_output_currents()[secondaryIndex];
                        auto outputVoltage = get_output_voltages()[secondaryIndex];
                        auto turnsRatio = turnsRatios.value()[secondaryIndex];
                        auto dutyCycleMaximum = calculate_BMO_duty_cycle((outputVoltage + diodeVoltageDrop), inputVoltage, turnsRatio);
                        currentPeak = std::max(currentPeak, calculate_BMO_primary_current_peak(outputCurrent, efficiency, dutyCycleMaximum, turnsRatio)); // hardcoded
                        double tOn = (currentPeak * inductance.value()) / inputVoltage;
                        double tOff = (currentPeak * inductance.value()) / (turnsRatio * outputVoltage);

                        switchingFrequency = std::max(switchingFrequency, 1.0 / (tOn + tOff));
                    }
                    return switchingFrequency;
                }
                default:
                  throw InvalidInputException(ErrorCode::INVALID_INPUT, "Unknown mode in Flyback");
            }
        }
    }

    OperatingPoint Flyback::process_operating_points_for_input_voltage(double inputVoltage, FlybackOperatingPoint outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance, std::optional<FlybackModes> customMode, std::optional<double> customDutyCycle, std::optional<double> customDeadTime) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.resolve_switching_frequency(inputVoltage, get_diode_voltage_drop(), inductance, turnsRatios, get_efficiency());
        
        if (switchingFrequency <= 0 || !std::isfinite(switchingFrequency)) {
            throw std::invalid_argument("Invalid switchingFrequency calculated: " + std::to_string(switchingFrequency));
        }

        double deadTime = 0;
        double maximumReflectedOutputVoltage = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
            double outputVoltage = outputOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop();
            maximumReflectedOutputVoltage = std::max(maximumReflectedOutputVoltage, outputVoltage * turnsRatios[secondaryIndex]);
        }

        double primaryVoltagePeaktoPeak = inputVoltage + maximumReflectedOutputVoltage;

        double totalOutputPower = get_total_input_power(outputOperatingPoint.get_output_currents(), outputOperatingPoint.get_output_voltages(), 1, 0);
        double maximumEffectiveLoadCurrent = totalOutputPower / outputOperatingPoint.get_output_voltages()[0];
        double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
        double totalInputPower = get_total_input_power(outputOperatingPoint.get_output_currents(), outputOperatingPoint.get_output_voltages(), get_efficiency(), 0);
        double averageInputCurrent = totalInputPower / inputVoltage;

        // Steady-state CCM duty (existing formula).
        double dCcm = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);

        // CCM/DCM boundary: critical primary inductance (Basso, "Switch-Mode Power Supplies",
        // 2nd ed., p. 747 — same expression already used in process_design_requirements below).
        // Below L_crit the converter is in DCM and the CCM duty formula is physically wrong.
        double mainOutputVoltageWithDiode = outputOperatingPoint.get_output_voltages()[0] + get_diode_voltage_drop();
        double aux = mainOutputVoltageWithDiode * turnsRatios[0];
        double criticalInductance = (totalOutputPower > 0)
            ? get_efficiency() * inputVoltage * inputVoltage * aux * aux
                / (2.0 * totalOutputPower * switchingFrequency
                   * (inputVoltage + aux) * (aux + get_efficiency() * inputVoltage))
            : 0.0;
        bool isDcm = inductance < criticalInductance;

        double dutyCycle;
        if (customDutyCycle) {
            dutyCycle = customDutyCycle.value();
        }
        else if (isDcm) {
            // DCM energy balance: Pout = ½·Lp·I_pk²·fsw with I_pk = Vin·D·T/Lp
            //   ⇒ D = sqrt(2 · Pin · Lp · fsw) / Vin    (Pin = Pout / efficiency)
            dutyCycle = std::sqrt(2.0 * totalInputPower * inductance * switchingFrequency) / inputVoltage;
        }
        else {
            dutyCycle = dCcm;
        }

        // Honour user-supplied maximumDutyCycle. Per project policy: throw, do not clamp.
        if (auto dMax = get_maximum_duty_cycle(); dMax && dutyCycle > dMax.value()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Required dutyCycle " + std::to_string(dutyCycle) +
                " exceeds maximumDutyCycle " + std::to_string(dMax.value()) +
                " at Vin=" + std::to_string(inputVoltage) + "V (mode=" +
                (isDcm ? "DCM" : "CCM") +
                "). Increase magnetizingInductance, lower switchingFrequency, or relax maximumDutyCycle.");
        }

        // If the caller forced a mode that disagrees with the physics, refuse — no silent fallback.
        // Only check user-explicit modes (set on the FlybackOperatingPoint), not ripple-heuristic
        // modes that the analytical pipeline derives and threads through customMode purely to pick
        // waveform labels.
        if (outputOperatingPoint.get_mode() && !customDutyCycle) {
            auto userMode = outputOperatingPoint.get_mode().value();
            bool forcedCcm = (userMode == FlybackModes::CONTINUOUS_CONDUCTION_MODE);
            bool forcedDcm = (userMode == FlybackModes::DISCONTINUOUS_CONDUCTION_MODE);
            if ((forcedCcm && isDcm) || (forcedDcm && !isDcm)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    std::string("Forced FlybackMode does not match the physics at Vin=") +
                    std::to_string(inputVoltage) + "V: requested " +
                    (forcedCcm ? "CCM" : "DCM") + " but Lp=" + std::to_string(inductance) +
                    " vs critical Lp=" + std::to_string(criticalInductance) + " puts it in " +
                    (isDcm ? "DCM" : "CCM") + ".");
            }
        }

        double centerSecondaryCurrentRampLumped = maximumEffectiveLoadCurrent / (1 - dutyCycle);
        double centerPrimaryCurrentRamp = centerSecondaryCurrentRampLumped / turnsRatios[0];


        if (customDeadTime) {
            deadTime = customDeadTime.value();
        }

        if (dutyCycle > 1 ) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "dutyCycle cannot be larger than one: " + std::to_string(dutyCycle));
        }

        double primaryCurrentAverage = centerPrimaryCurrentRamp;
        double currentRippleRatio;
        if (std::isnan(get_current_ripple_ratio())) {
            double primaryCurrentPeakToPeak = inputVoltage * dutyCycle / switchingFrequency / inductance;
            currentRippleRatio = primaryCurrentPeakToPeak / centerPrimaryCurrentRamp;
        }
        else {
            currentRippleRatio = get_current_ripple_ratio();
        }
        double primaryCurrentPeakToPeak = centerPrimaryCurrentRamp * currentRippleRatio * 2;
        double primaryCurrentOffset = primaryCurrentAverage - primaryCurrentPeakToPeak / 2;
        primaryCurrentOffset = std::max(0.0, primaryCurrentOffset);

        FlybackModes mode;
        if (customMode) {
            mode = customMode.value();
        }
        else {
            if (primaryCurrentOffset > 0) {
                mode = FlybackModes::CONTINUOUS_CONDUCTION_MODE;
            }
            else {
                mode = FlybackModes::DISCONTINUOUS_CONDUCTION_MODE;
            }
        }

        // ---- Populate per-OP diagnostics (golden-quality contract) ----
        // Scalar last_* fields reflect the MOST RECENT call. Vectors below
        // accumulate one entry per call so the wizard's Diagnostics table
        // can show every operating point in its own column.
        lastMode = mode;
        lastDutyCycle = dutyCycle;
        lastSwitchingFrequency = switchingFrequency;
        lastPrimaryAverageCurrent = primaryCurrentAverage;
        lastPrimaryPeakToPeak = primaryCurrentPeakToPeak;
        lastPrimaryPeakCurrent = primaryCurrentOffset + primaryCurrentPeakToPeak;
        lastSecondaryPeakCurrent = lastPrimaryPeakCurrent * turnsRatios[0];
        lastIsCcm = (mode == FlybackModes::CONTINUOUS_CONDUCTION_MODE);

        perOpMode.push_back(mode);
        perOpDutyCycle.push_back(dutyCycle);
        perOpSwitchingFrequency.push_back(switchingFrequency);
        perOpPrimaryAverageCurrent.push_back(primaryCurrentAverage);
        perOpPrimaryPeakToPeak.push_back(primaryCurrentPeakToPeak);
        perOpPrimaryPeakCurrent.push_back(lastPrimaryPeakCurrent);
        perOpSecondaryPeakCurrent.push_back(lastSecondaryPeakCurrent);
        perOpIsCcm.push_back(lastIsCcm);

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, primaryCurrentOffset, deadTime);

            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    break;
                }
                case FlybackModes::QUASI_RESONANT_MODE:
                case FlybackModes::BOUNDARY_MODE_OPERATION:
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    break;
                }
            }

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // Secondaries
        for (size_t secondaryIndex = 0; secondaryIndex < turnsRatios.size(); ++secondaryIndex) {

            Waveform currentWaveform;
            Waveform voltageWaveform;

            double secondaryPower = get_total_input_power(outputOperatingPoint.get_output_currents()[secondaryIndex], outputOperatingPoint.get_output_voltages()[secondaryIndex], 1, 0);
            double powerDivider = secondaryPower / totalOutputPower;

            double minimumSecondaryVoltage = -inputVoltage / turnsRatios[secondaryIndex];
            double maximumSecondaryVoltage = outputOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop();
            double secondaryVoltagePeaktoPeak = maximumSecondaryVoltage - minimumSecondaryVoltage;
            double secondaryCurrentAverage = centerPrimaryCurrentRamp * turnsRatios[secondaryIndex] * powerDivider;
            double secondaryCurrentPeaktoPeak = secondaryCurrentAverage * currentRippleRatio * 2;
            double secondaryCurrentOffset = std::max(0.0, secondaryCurrentAverage - secondaryCurrentPeaktoPeak / 2);

            // The SECONDARY_RECTANGULAR[_WITH_DEADTIME] *voltage* waveform builder ties
            // on/off amplitudes to the duty-cycle parameter (max = -peakToPeak·(1-D),
            // min = +peakToPeak·D). For the positive level to come out to (Vout+Vd) the
            // duty must be the volt-second-balance value D_sec = aux/(Vin+aux),
            // with aux = (Vout+Vd)·N_sec. In CCM the primary D already satisfies this
            // (modulo small η/Vd offsets the existing pipeline expects), so we keep
            // primary D for CCM to preserve historical behaviour. In DCM the primary D
            // is smaller than the volt-second value, which would shrink the secondary
            // rectangular's positive level below Vout+Vd; use D_sec for the secondary
            // *voltage* in that case so the reported peak is physically right.
            //
            // The FLYBACK_SECONDARY[_WITH_DEADTIME] *current* waveform parameter
            // `dutyCycle` has different semantics — it controls the ramp on/off split
            // such that the period-integrated average matches the requested load
            // current. Forcing D_sec there breaks I_avg = I_load. Always pass the
            // primary `dutyCycle` to the current waveforms.
            double secondaryVoltageDuty = dutyCycle;
            if (isDcm) {
                double secondaryAux = maximumSecondaryVoltage * turnsRatios[secondaryIndex];
                secondaryVoltageDuty = secondaryAux / (inputVoltage + secondaryAux);
            }

            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::SECONDARY_RECTANGULAR, secondaryVoltagePeaktoPeak, switchingFrequency, secondaryVoltageDuty, 0, deadTime);
                    currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    break;
                }
                case FlybackModes::QUASI_RESONANT_MODE:
                case FlybackModes::BOUNDARY_MODE_OPERATION:
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME, secondaryVoltagePeaktoPeak, switchingFrequency, secondaryVoltageDuty, 0, deadTime);
                    currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    break;
                }
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

    std::string Flyback::generate_ngspice_circuit(
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
        auto opPoint = get_mutable_operating_points()[operatingPointIndex];  // Copy to allow modification
        
        // Calculate switching frequency and duty cycle
        double switchingFrequency = opPoint.resolve_switching_frequency(inputVoltage, get_diode_voltage_drop(), magnetizingInductance, turnsRatios, get_efficiency());
        
        // Calculate duty cycle (mirrors process_operating_points_for_input_voltage):
        // CCM-then-DCM physics + maxD enforcement. See that routine for derivation.
        double totalOutputPower = get_total_input_power(opPoint.get_output_currents(), opPoint.get_output_voltages(), 1, 0);
        double maximumEffectiveLoadCurrent = totalOutputPower / opPoint.get_output_voltages()[0];
        double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
        double totalInputPower = get_total_input_power(opPoint.get_output_currents(), opPoint.get_output_voltages(), get_efficiency(), 0);
        double averageInputCurrent = totalInputPower / inputVoltage;
        double dCcm = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
        double mainOutputVoltageWithDiode = opPoint.get_output_voltages()[0] + get_diode_voltage_drop();
        double aux = mainOutputVoltageWithDiode * turnsRatios[0];
        double criticalInductance = (totalOutputPower > 0)
            ? get_efficiency() * inputVoltage * inputVoltage * aux * aux
                / (2.0 * totalOutputPower * switchingFrequency
                   * (inputVoltage + aux) * (aux + get_efficiency() * inputVoltage))
            : 0.0;
        bool isDcm = magnetizingInductance < criticalInductance;
        double dutyCycle = isDcm
            ? std::sqrt(2.0 * totalInputPower * magnetizingInductance * switchingFrequency) / inputVoltage
            : dCcm;
        if (auto dMax = get_maximum_duty_cycle(); dMax && dutyCycle > dMax.value()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Required dutyCycle " + std::to_string(dutyCycle) +
                " exceeds maximumDutyCycle " + std::to_string(dMax.value()) +
                " at Vin=" + std::to_string(inputVoltage) + "V (mode=" +
                (isDcm ? "DCM" : "CCM") + ") in generate_ngspice_circuit.");
        }

        // Capture the SPICE-input duty cycle (and frequency) so libMKF.cpp
        // can emit it as a "what the simulator actually ran" diagnostic
        // without re-deriving it from the trace (where ramp+commutation
        // would otherwise need careful zero-crossing detection).
        lastDutyCycle = dutyCycle;
        lastSwitchingFrequency = switchingFrequency;

        // Number of secondaries
        size_t numSecondaries = turnsRatios.size();

        // SPICE configuration knobs (registered defaults, see
        // Topology.cpp::spice_simulation_defaults()). Per-instance
        // overrides apply via Topology::set_spice_config().
        const SpiceSimulationConfig cfg = spice_config();

        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        
        // Simulation: run steady-state periods for settling, then extract the last N periods
        int periodsToExtract = get_num_periods_to_extract();
        int numSteadyStatePeriods = get_num_steady_state_periods();
        const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;  // Steady state + extraction
        double simTime = numPeriodsTotal * period;
        double startTime = numSteadyStatePeriods * period;  // Start extracting after steady state
        double stepTime = period / cfg.samplesPerPeriod;
        
        circuit << "* Flyback Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lp=" << (magnetizingInductance*1e6) << "uH, N1=" << turnsRatios[0] << ", " << numSecondaries << " secondaries\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // PWM Switch (ideal). PULSE timing emitted in scientific notation
        // so sub-µs values aren't silently rounded by std::fixed
        // (precision 6) to the wrong duty cycle / Fs.
        circuit << "* PWM Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                << tOn << " " << period << std::fixed
                << ")\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH << "\n";
        circuit << "S1 vin_dc pri_p pwm_ctrl 0 SW1\n";
        // Snubber RC across the switch — absorbs leakage-driven ringing
        // at turn-off and improves convergence near CCM/DCM boundary.
        circuit << "Rsnub_s1 vin_dc pri_p " << cfg.snubR << "\n"
                << "Csnub_s1 vin_dc pri_p " << std::scientific << cfg.snubC << std::fixed << "\n\n";
        
        // Primary current sense (for waveform extraction)
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense pri_p pri_in 0\n\n";
        
        // Flyback Transformer (ideal coupling = 1 for all windings)
        // NOTE: Secondary inductors have terminals swapped (0 to sec_N_in) to create
        // opposite dot polarity needed for flyback operation.
        circuit << "* Flyback Transformer - Primary and " << numSecondaries << " secondaries\n";
        circuit << "Lpri pri_in 0 " << std::scientific << magnetizingInductance << std::fixed << "\n";
        
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double secondaryInductance = magnetizingInductance / (turnsRatios[secIdx] * turnsRatios[secIdx]);
            circuit << "Lsec" << secIdx << " 0 sec" << secIdx << "_in " << std::scientific << secondaryInductance << std::fixed << "\n";
        }
        
        // ngspice requires separate K statements for each inductor pair
        // Couple primary to each secondary
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "K" << secIdx << " Lpri Lsec" << secIdx << " 1\n";
        }
        // Couple secondaries to each other (for proper cross-regulation)
        for (size_t i = 0; i < numSecondaries; ++i) {
            for (size_t j = i + 1; j < numSecondaries; ++j) {
                circuit << "K" << numSecondaries << "_" << i << "_" << j << " Lsec" << i << " Lsec" << j << " 1\n";
            }
        }
        circuit << "\n";
        
        // Output Rectifiers, current sense, and loads for each secondary
        circuit << "* Output Rectifier model\n";
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::fixed << ")\n\n";
        
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "* Secondary " << secIdx << " output stage\n";
            circuit << "Dout" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_p DIDEAL\n";
            circuit << "Vsec_sense" << secIdx << " sec" << secIdx << "_p vout" << secIdx << " 0\n";
            
            double outputCurrent = opPoint.get_output_currents()[secIdx];
            double outputVoltage = opPoint.get_output_voltages()[secIdx];
            if (outputCurrent <= 0.0) {
                throw std::invalid_argument(
                    "Flyback::generate_ngspice_circuit: outputCurrent must be > 0 "
                    "for secondary " + std::to_string(secIdx) +
                    " (received " + std::to_string(outputCurrent) +
                    "); cannot derive load resistance R_load = Vout / Iout.");
            }
            double loadResistance = outputVoltage / outputCurrent;
            circuit << "Cout" << secIdx << " vout" << secIdx << " 0 "
                    << std::scientific << cfg.outputCapacitance << std::fixed
                    << " IC=" << (outputVoltage * cfg.outputCapInitialChargeFraction) << "\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << "\n\n";
        
        // Save primary and all secondary voltages and currents
        circuit << "* Output signals\n";
        circuit << ".save v(pri_in) v(vin_dc)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) v(vout" << secIdx << ") i(Vsec_sense" << secIdx << ")";
        }
        circuit << " i(Vpri_sense)\n\n";
        
        // Solver options for convergence in switching circuits.
        // METHOD=GEAR + larger TRTOL handles the hard-switching event at
        // each PWM edge robustly.
        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::fixed
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        circuit << ".options METHOD=" << cfg.method << " TRTOL=" << cfg.trTol << "\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << ".ic v(vout" << secIdx << ")=" << opPoint.get_output_voltages()[secIdx] << "\n";
        }
        circuit << "\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> Flyback::simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods) {
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // If numberOfPeriods is specified, use it; otherwise use the member variable
        int numPeriods = (numberOfPeriods > 0) ? static_cast<int>(numberOfPeriods) : get_num_periods_to_extract();
        
        // Save original value and set the requested number of periods
        int originalNumPeriodsToExtract = get_num_periods_to_extract();
        set_num_periods_to_extract(numPeriods);
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            double inputVoltage = inputVoltages[inputVoltageIndex];
            
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto flybackOpPoint = get_mutable_operating_points()[opIndex];  // Copy to allow modification
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                // Calculate switching frequency for waveform extraction
                double switchingFrequency = flybackOpPoint.resolve_switching_frequency(
                    inputVoltage, get_diode_voltage_drop(), magnetizingInductance, turnsRatios, get_efficiency());
                
                // Run simulation
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);
                
                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Define waveform name mapping for flyback circuit with multiple secondaries
                // ngspice .save outputs voltage nodes as just the node name (e.g., "pri_in")
                // and branch currents as "vsourcename#branch" (e.g., "vpri_sense#branch")
                NgspiceRunner::WaveformNameMapping waveformMapping;
                
                // Primary winding
                waveformMapping.push_back({{"voltage", "pri_in"}, {"current", "vpri_sense#branch"}});
                
                // Secondary windings (one per turns ratio)
                for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
                    std::string voltageName = "sec" + std::to_string(secIdx) + "_in";
                    std::string currentName = "vsec_sense" + std::to_string(secIdx) + "#branch";
                    waveformMapping.push_back({{"voltage", voltageName}, {"current", currentName}});
                }
                
                std::vector<std::string> windingNames;
                windingNames.push_back("Primary");
                for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
                    windingNames.push_back("Secondary " + std::to_string(secIdx));
                }
                
                std::vector<bool> flipCurrentSign(1 + turnsRatios.size(), false);  // Keep current signs as-is
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    flybackOpPoint.get_ambient_temperature(),
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
        
        // Restore original value
        set_num_periods_to_extract(originalNumPeriodsToExtract);
        
        return operatingPoints;
    }

    std::string Flyback::generate_ngspice_circuit_with_magnetic(
        const Magnetic& magneticConst,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {
        
        // Make a copy since some methods are non-const
        Magnetic magnetic = magneticConst;
        
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
        auto opPoint = get_mutable_operating_points()[operatingPointIndex];
        
        // Extract turns ratios and magnetizing inductance from the Magnetic
        auto coil = magnetic.get_coil();
        size_t numWindings = coil.get_functional_description().size();
        if (numWindings < 2) {
            throw std::invalid_argument("Magnetic must have at least 2 windings for flyback");
        }
        
        double primaryTurns = coil.get_functional_description()[0].get_number_turns();
        std::vector<double> turnsRatios;
        for (size_t i = 1; i < numWindings; ++i) {
            double secondaryTurns = coil.get_functional_description()[i].get_number_turns();
            turnsRatios.push_back(primaryTurns / secondaryTurns);
        }
        
        double magnetizingInductance = resolve_dimensional_values(
            MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
        
        // Calculate switching frequency and duty cycle
        double switchingFrequency = opPoint.resolve_switching_frequency(
            inputVoltage, get_diode_voltage_drop(), magnetizingInductance, turnsRatios, get_efficiency());
        
        double totalOutputPower = get_total_input_power(opPoint.get_output_currents(), opPoint.get_output_voltages(), 1, 0);
        double maximumEffectiveLoadCurrent = totalOutputPower / opPoint.get_output_voltages()[0];
        double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
        double totalInputPower = get_total_input_power(opPoint.get_output_currents(), opPoint.get_output_voltages(), get_efficiency(), 0);
        double averageInputCurrent = totalInputPower / inputVoltage;
        // CCM-then-DCM physics + maxD enforcement (mirrors process_operating_points_for_input_voltage).
        double dCcm = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
        double mainOutputVoltageWithDiode = opPoint.get_output_voltages()[0] + get_diode_voltage_drop();
        double aux = mainOutputVoltageWithDiode * turnsRatios[0];
        double criticalInductance = (totalOutputPower > 0)
            ? get_efficiency() * inputVoltage * inputVoltage * aux * aux
                / (2.0 * totalOutputPower * switchingFrequency
                   * (inputVoltage + aux) * (aux + get_efficiency() * inputVoltage))
            : 0.0;
        bool isDcm = magnetizingInductance < criticalInductance;
        double dutyCycle = isDcm
            ? std::sqrt(2.0 * totalInputPower * magnetizingInductance * switchingFrequency) / inputVoltage
            : dCcm;
        if (auto dMax = get_maximum_duty_cycle(); dMax && dutyCycle > dMax.value()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Required dutyCycle " + std::to_string(dutyCycle) +
                " exceeds maximumDutyCycle " + std::to_string(dMax.value()) +
                " at Vin=" + std::to_string(inputVoltage) + "V (mode=" +
                (isDcm ? "DCM" : "CCM") + ") in generate_ngspice_circuit_with_magnetic.");
        }
        
        // Generate the magnetic subcircuit using CircuitSimulatorExporter
        CircuitSimulatorExporterNgspiceModel ngspiceExporter;
        std::string magneticSubcircuit = ngspiceExporter.export_magnetic_as_subcircuit(
            magnetic, switchingFrequency, opPoint.get_ambient_temperature());
        
        // Number of secondaries
        size_t numSecondaries = turnsRatios.size();

        // SPICE configuration knobs (registered defaults).
        const SpiceSimulationConfig cfg = spice_config();

        // Build the main circuit netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        
        // Simulation: run 10x the extraction periods for settling, then extract the last N periods
        int periodsToExtract = get_num_periods_to_extract();
        const int numPeriodsTotal = 10 * periodsToExtract;  // Run 10x more for settling
        double simTime = numPeriodsTotal * period;
        double startTime = (numPeriodsTotal - periodsToExtract) * period;
        double stepTime = period / cfg.samplesPerPeriod;
        
        circuit << "* Flyback Converter with Real Magnetic Component - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Magnetic: " << magnetic.get_reference() << "\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, N1=" << turnsRatios[0] << ", " << numSecondaries << " secondaries\n\n";
        
        // Include the magnetic subcircuit definition
        circuit << "* Magnetic Component Subcircuit\n";
        circuit << magneticSubcircuit << "\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // PWM Switch (ideal). Scientific PULSE timing — see notes in
        // generate_ngspice_circuit() above.
        circuit << "* PWM Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                << tOn << " " << period << std::fixed
                << ")\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH << "\n";
        circuit << "S1 vin_dc pri_p pwm_ctrl 0 SW1\n";
        // Real-magnetic path: use snubRReal if set (handles real Lk·di/dt
        // turn-off spike clamping), else fall back to snubR.
        double snubR_eff = cfg.snubRReal.value_or(cfg.snubR);
        circuit << "Rsnub_s1 vin_dc pri_p " << snubR_eff << "\n"
                << "Csnub_s1 vin_dc pri_p " << std::scientific << cfg.snubC << std::fixed << "\n\n";

        // Primary current sense (for waveform extraction)
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense pri_p pri_in 0\n\n";

        // Instantiate the magnetic component subcircuit
        // Subcircuit pins are: P1+ P1- P2+ P2- ... for each winding
        // We connect: pri_in to P1+, 0 to P1-, sec_N_in to PN+, 0 to PN-
        std::string subcktName = fix_filename(magnetic.get_reference());
        circuit << "* Magnetic component instance\n";
        circuit << "X1 pri_in 0";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " sec" << secIdx << "_in 0";
        }
        circuit << " " << subcktName << "\n\n";
        
        // Output Rectifier model
        circuit << "* Output Rectifier model\n";
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::fixed << ")\n\n";
        
        // Output stage for each secondary
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "* Secondary " << secIdx << " output stage\n";
            circuit << "Dout" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_p DIDEAL\n";
            circuit << "Vsec_sense" << secIdx << " sec" << secIdx << "_p vout" << secIdx << " 0\n";
            
            double outputCurrent = opPoint.get_output_currents()[secIdx];
            double outputVoltage = opPoint.get_output_voltages()[secIdx];
            if (outputCurrent <= 0.0) {
                throw std::invalid_argument(
                    "Flyback::generate_ngspice_circuit_with_magnetic: outputCurrent must be > 0 "
                    "for secondary " + std::to_string(secIdx));
            }
            double loadResistance = outputVoltage / outputCurrent;
            circuit << "Cout" << secIdx << " vout" << secIdx << " 0 "
                    << std::scientific << cfg.outputCapacitance << std::fixed
                    << " IC=" << (outputVoltage * cfg.outputCapInitialChargeFraction) << "\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << "\n\n";
        
        // Save primary and all secondary voltages and currents
        circuit << "* Output signals\n";
        circuit << ".save v(pri_in) v(vin_dc)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) v(vout" << secIdx << ") i(Vsec_sense" << secIdx << ")";
        }
        circuit << " i(Vpri_sense)\n\n";
        
        // Solver options (GEAR for hard-switching robustness)
        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::fixed
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        circuit << ".options METHOD=" << cfg.method << " TRTOL=" << cfg.trTol << "\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << ".ic v(vout" << secIdx << ")=" << opPoint.get_output_voltages()[secIdx] << "\n";
        }
        circuit << "\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }

    std::vector<OperatingPoint> Flyback::simulate_with_magnetic_and_extract_operating_points(
        const Magnetic& magneticConst) {
        
        // Make a copy since some methods are non-const
        Magnetic magnetic = magneticConst;
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Extract turns ratios and magnetizing inductance from the Magnetic
        auto coil = magnetic.get_coil();
        size_t numWindings = coil.get_functional_description().size();
        if (numWindings < 2) {
            throw std::invalid_argument("Magnetic must have at least 2 windings for flyback");
        }
        
        double primaryTurns = coil.get_functional_description()[0].get_number_turns();
        std::vector<double> turnsRatios;
        for (size_t i = 1; i < numWindings; ++i) {
            double secondaryTurns = coil.get_functional_description()[i].get_number_turns();
            turnsRatios.push_back(primaryTurns / secondaryTurns);
        }
        
        double magnetizingInductance = resolve_dimensional_values(
            MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            double inputVoltage = inputVoltages[inputVoltageIndex];
            
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto flybackOpPoint = get_mutable_operating_points()[opIndex];
                
                // Generate circuit with real magnetic
                std::string netlist = generate_ngspice_circuit_with_magnetic(magnetic, inputVoltageIndex, opIndex);
                
                // Calculate switching frequency for waveform extraction
                double switchingFrequency = flybackOpPoint.resolve_switching_frequency(
                    inputVoltage, get_diode_voltage_drop(), magnetizingInductance, turnsRatios, get_efficiency());
                
                // Run simulation
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);

                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }

                // Define waveform name mapping for flyback circuit with multiple secondaries
                // ngspice .save outputs voltage nodes as just the node name (e.g., "pri_in")
                // and branch currents as "vsourcename#branch" (e.g., "vpri_sense#branch")
                NgspiceRunner::WaveformNameMapping waveformMapping;

                // Primary winding
                waveformMapping.push_back({{"voltage", "pri_in"}, {"current", "vpri_sense#branch"}});

                // Secondary windings (one per turns ratio)
                for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
                    std::string voltageName = "sec" + std::to_string(secIdx) + "_in";
                    std::string currentName = "vsec_sense" + std::to_string(secIdx) + "#branch";
                    waveformMapping.push_back({{"voltage", voltageName}, {"current", currentName}});
                }

                std::vector<std::string> windingNames;
                windingNames.push_back("Primary");
                for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
                    windingNames.push_back("Secondary " + std::to_string(secIdx));
                }

                std::vector<bool> flipCurrentSign(1 + turnsRatios.size(), false);

                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    flybackOpPoint.get_ambient_temperature(),
                    flipCurrentSign);

                // Set name
                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt. (simulated with " + magnetic.get_reference() + ")";
                if (get_operating_points().size() > 1) {
                    name += " op. point " + std::to_string(opIndex);
                }
                operatingPoint.set_name(name);

                operatingPoints.push_back(operatingPoint);
            }
        }
        
        return operatingPoints;
    }

    std::vector<ConverterWaveforms> Flyback::simulate_and_extract_topology_waveforms(const std::vector<double>& turnsRatios,
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
            double switchingFrequency = opPoint.get_switching_frequency().value();
            
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
            
            // Per CONVERTER_MODELS_GOLDEN_GUIDE.md §5.0/§5.1: ConverterWaveforms is
            // the converter-port stream. input_voltage MUST be the DC source
            // (v(vin_dc)), NOT the bipolar across-Lpri swing v(pri_in) — the
            // latter belongs in excitations_per_winding[primary] (and is
            // already exposed via simulate_and_extract_operating_points).
            wf.set_input_voltage(getWaveform("vin_dc"));
            wf.set_input_current(getWaveform("vpri_sense#branch"));

            for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
                wf.get_mutable_output_voltages().push_back(getWaveform("vout" + std::to_string(secIdx)));
                wf.get_mutable_output_currents().push_back(getWaveform("vsec_sense" + std::to_string(secIdx) + "#branch"));
            }
            
            results.push_back(wf);
        }
    }
    
    // Restore original value
    set_num_periods_to_extract(originalNumPeriodsToExtract);
    
    return results;
}

    double Flyback::get_total_input_power(std::vector<double> outputCurrents, std::vector<double> outputVoltages, double efficiency, double diodeVoltageDrop) {
        double totalPower = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputCurrents.size(); ++secondaryIndex) {
            totalPower += outputCurrents[secondaryIndex] * (outputVoltages[secondaryIndex] + diodeVoltageDrop);
        }

        return totalPower / efficiency;
    }


    double Flyback::get_total_input_current(std::vector<double> outputCurrents, double inputVoltage, std::vector<double> outputVoltages, double diodeVoltageDrop) {
        double totalCurrent = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputCurrents.size(); ++secondaryIndex) {
            totalCurrent += outputCurrents[secondaryIndex] * (outputVoltages[secondaryIndex] + diodeVoltageDrop) / inputVoltage;
        }

        return totalCurrent;
    }


    double Flyback::get_total_input_power(double outputCurrent, double outputVoltage, double efficiency, double diodeVoltageDrop) {
        double totalPower = outputCurrent * (outputVoltage + diodeVoltageDrop);

        return totalPower / efficiency;
    }

    double Flyback::get_minimum_output_reflected_voltage(double maximumDrainSourceVoltage, double maximumInputVoltage, double safetyMargin) {
        return maximumDrainSourceVoltage * safetyMargin - maximumInputVoltage;
    }

    bool Flyback::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        for (size_t flybackOperatingPointIndex = 1; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            if (get_operating_points()[flybackOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[flybackOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
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

    DesignRequirements Flyback::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_maximum_drain_source_voltage() && !get_maximum_duty_cycle()) {
            throw std::invalid_argument("Missing both maximum duty cycle and maximum drain source voltage");
        }
        if (get_operating_points().empty()) {
            throw std::invalid_argument("Flyback requires at least one operating point");
        }
        double globalNeededInductance = 0;
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size(), 0);

        if (get_maximum_duty_cycle()) {
            double maximumDutyCycle = get_maximum_duty_cycle().value();
            if (maximumDutyCycle > 1 || maximumDutyCycle < 0) {
                throw std::invalid_argument("maximumDutyCycle must be between 0 and 1");
            }
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];

                double totalOutputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), 1, 0);
                double totalInputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), get_efficiency(), 0);
                double maximumEffectiveLoadCurrent = totalOutputPower / flybackOperatingPoint.get_output_voltages()[0];
                double averageInputCurrent = totalInputPower / minimumInputVoltage;
                double maximumEffectiveLoadCurrentReflected = averageInputCurrent * (1 - maximumDutyCycle) / maximumDutyCycle;

                auto turnsRatioFirstOutput = maximumEffectiveLoadCurrent / maximumEffectiveLoadCurrentReflected;
                turnsRatios[0] = std::max(turnsRatios[0], turnsRatioFirstOutput);

                for (size_t secondaryIndex = 1; secondaryIndex < flybackOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                    auto turnsRatio = turnsRatioFirstOutput * (flybackOperatingPoint.get_output_voltages()[0] + get_diode_voltage_drop()) / (flybackOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
                    turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatio);
                }
            }
        }

        if (get_maximum_drain_source_voltage()) {
            std::vector<double> turnsRatiosFromMaximumDrainSourceVoltage(get_operating_points()[0].get_output_voltages().size(), 0);
            double maximumDrainSourceVoltage = get_maximum_drain_source_voltage().value();
            auto minimumOutputReflectedVoltage = get_minimum_output_reflected_voltage(maximumDrainSourceVoltage, maximumInputVoltage);
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
                for (size_t secondaryIndex = 0; secondaryIndex < flybackOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                    auto turnsRatio = minimumOutputReflectedVoltage / (flybackOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
                    turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex] = std::max(turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex], turnsRatio);
                }
            }

            for (size_t secondaryIndex = 0; secondaryIndex < get_operating_points()[0].get_output_voltages().size(); ++secondaryIndex) {
                if (turnsRatios[secondaryIndex] > 1) {
                    turnsRatios[secondaryIndex] = std::min(turnsRatios[secondaryIndex], turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex]);
                }
                else {
                    turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex]);
                }
            }
        }

        for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
            double switchingFrequency = flybackOperatingPoint.resolve_switching_frequency(minimumInputVoltage, get_diode_voltage_drop());
            double totalOutputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), 1, 0);
            double maximumEffectiveLoadCurrent = totalOutputPower / flybackOperatingPoint.get_output_voltages()[0];
            double dutyCycle = 0;
            if (get_maximum_duty_cycle()) {
                dutyCycle = get_maximum_duty_cycle().value();
            }
            else {
                double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
                double totalInputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), get_efficiency(), 0);
                double averageInputCurrent = totalInputPower / minimumInputVoltage;
                dutyCycle = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
            }

            double centerSecondaryCurrentRampLumped = maximumEffectiveLoadCurrent / (1 - dutyCycle);
            double centerPrimaryCurrentRamp = centerSecondaryCurrentRampLumped / turnsRatios[0];
            double tOn = dutyCycle / switchingFrequency;
            double voltsSeconds = minimumInputVoltage * tOn;
            double neededInductance = voltsSeconds / get_current_ripple_ratio() / centerPrimaryCurrentRamp;
            globalNeededInductance = std::max(globalNeededInductance, neededInductance);
        }

        // Add maximum inductance to stay in DCM
        double maximumInductance = 0;
        for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
            if (flybackOperatingPoint.get_mode()) {
                if (flybackOperatingPoint.get_mode().value() != MAS::FlybackModes::CONTINUOUS_CONDUCTION_MODE) {
                    double totalOutputPower = Flyback::get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), 1, get_diode_voltage_drop());
                    double switchingFrequency = flybackOperatingPoint.resolve_switching_frequency(minimumInputVoltage, get_diode_voltage_drop());
                    double mainOutputVoltage = flybackOperatingPoint.get_output_voltages()[0];
                    double aux = (mainOutputVoltage + get_diode_voltage_drop()) * turnsRatios[0];
                    // Accoding to Switch-Mode Power Supplies, 2nd ed; Christophe Basso; page 747.
                    double maximumInductanceThisPoint = get_efficiency() * pow(minimumInputVoltage, 2) * pow(aux, 2) / (2 * totalOutputPower * switchingFrequency * (minimumInputVoltage + aux) * (aux + get_efficiency() * minimumInputVoltage));
                    maximumInductance = std::max(maximumInductance, maximumInductanceThisPoint);
                }
            }
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        bool dutyConstrained = get_maximum_duty_cycle().has_value();
        for (auto turnsRatio : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            double realizedTurnsRatio = roundFloat(turnsRatio, 2);
            // The maximum-duty-cycle turns ratio is an UPPER bound: the reflected secondary voltage
            // (and hence the required duty) RISES with the turns ratio, so the derivation above sets it
            // to the largest ratio that still meets D <= maximumDutyCycle at minimum input voltage.
            // Rounding that ratio UP to the realizable 2-decimal precision pushes the recomputed duty
            // just past the ceiling and trips the duty check in process_operating_points_for_input_voltage
            // (e.g. 2.2091 -> 2.21 takes D from 0.4500 to 0.4501 against a 0.45 ceiling). When the rounding
            // would exceed the derived ratio, floor it instead so the realized design honours the ceiling.
            if (dutyConstrained && realizedTurnsRatio > turnsRatio) {
                realizedTurnsRatio = floorFloat(turnsRatio, 2);
            }
            turnsRatioWithTolerance.set_nominal(realizedTurnsRatio);
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }
        DimensionWithTolerance inductanceWithTolerance;
        // If any operating point forces DCM/QRM/BMO and the Basso critical inductance
        // is below the ripple-sized "needed" inductance, the latter is a CCM-regime
        // value that contradicts the forced mode. Snap nominal to Lp_crit so the
        // design is self-consistent (otherwise the operating-point pipeline would
        // throw "Forced FlybackMode does not match the physics").
        double targetInductance = globalNeededInductance;
        if (maximumInductance > 0 && maximumInductance < globalNeededInductance) {
            targetInductance = maximumInductance;
        }
        // Set nominal to the calculated inductance - this is the target value
        inductanceWithTolerance.set_nominal(roundFloat(targetInductance, 10));
        // Set minimum slightly below nominal (allow 10% tolerance)
        inductanceWithTolerance.set_minimum(roundFloat(targetInductance * 0.9, 10));

        if (maximumInductance > 0) {
            // Ensure maximum is not smaller than nominal (can happen in edge cases)
            if (maximumInductance >= targetInductance) {
                inductanceWithTolerance.set_maximum(roundFloat(maximumInductance, 10));
            } else {
                // If max < nominal, use nominal + 20% as max
                inductanceWithTolerance.set_maximum(roundFloat(targetInductance * 1.2, 10));
            }
        } else {
            // No DCM constraint, allow 20% above nominal
            inductanceWithTolerance.set_maximum(roundFloat(targetInductance * 1.2, 10));
        }

        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(MAS::Topology::FLYBACK_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Flyback::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        // Clear per-OP diagnostic vectors so the wizard's Diagnostics table
        // reflects only this run, not any previous one.
        perOpName.clear();
        perOpMode.clear();
        perOpDutyCycle.clear();
        perOpSwitchingFrequency.clear();
        perOpPrimaryAverageCurrent.clear();
        perOpPrimaryPeakToPeak.clear();
        perOpPrimaryPeakCurrent.clear();
        perOpSecondaryPeakCurrent.clear();
        perOpIsCcm.clear();

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto mode = get_mutable_operating_points()[flybackOperatingPointIndex].resolve_mode(get_current_ripple_ratio());

                std::string opName = inputVoltagesNames[inputVoltageIndex];
                if (get_operating_points().size() > 1) {
                    opName += " · OP" + std::to_string(flybackOperatingPointIndex);
                }
                perOpName.push_back(opName);

                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, magnetizingInductance, mode);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(flybackOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> Flyback::process_operating_points(OpenMagnetics::Magnetic magnetic) {
        run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    DesignRequirements AdvancedFlyback::process_design_requirements() {
        // Override the parent Flyback::process_design_requirements()
        // (which requires maximumDutyCycle or maximumDrainSourceVoltage
        // to be set — see Issue M1). The Advanced workflow already
        // pins L and n via desiredInductance / desiredTurnsRatios, so
        // we build the DR directly from those fields. This is the
        // same DR that AdvancedFlyback::process() constructs internally.
        std::vector<double> turnsRatios = get_desired_turns_ratios();
        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        for (auto turnsRatio : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(get_desired_inductance(), 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(MAS::Topology::FLYBACK_CONVERTER);
        return designRequirements;
    }

    Inputs AdvancedFlyback::process() {
        Flyback::run_checks(_assertErrors);

        Inputs inputs;

        double globalNeededInductance = get_desired_inductance();
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
        inductanceWithTolerance.set_nominal(roundFloat(globalNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(MAS::Topology::FLYBACK_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                std::optional<double> customDeadTime = std::nullopt;
                if (get_desired_duty_cycle().size() <= flybackOperatingPointIndex) {
                    throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing duty cycle for flybackOperatingPointIndex: " + std::to_string(flybackOperatingPointIndex));
                }
                double customDutyCycle = get_desired_duty_cycle()[flybackOperatingPointIndex][inputVoltageIndex];

                if (get_desired_dead_time()) {
                    if (get_desired_dead_time()->size() <= flybackOperatingPointIndex) {
                        throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing dead time for flybackOperatingPointIndex: " + std::to_string(flybackOperatingPointIndex));
                    }
                    customDeadTime = get_desired_dead_time().value()[flybackOperatingPointIndex];
                }

                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, globalNeededInductance, std::nullopt, customDutyCycle, customDeadTime);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(flybackOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }
} // namespace OpenMagnetics
