#include "converter_models/FourSwitchBuckBoost.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/NgspiceRunner.h"
#include "support/Utils.h"
#include <cfloat>
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <sstream>
#include "support/Exceptions.h"
#include "converter_models/Topology.h"

namespace OpenMagnetics {

    // ============================================================
    //   Static analytical helpers
    // ============================================================

    FourSwitchBuckBoost::Region FourSwitchBuckBoost::classify_region(
        double Vin, double Vo, double hysteresisRatio,
        double minOnTime, double minOffTime, double Fs)
    {
        if (Vin <= 0.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::classify_region: Vin must be > 0");
        }
        if (Vo <= 0.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::classify_region: Vo must be > 0");
        }
        if (Fs <= 0.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::classify_region: Fs must be > 0");
        }
        const double M = Vo / Vin;
        const double minOnDuty  = minOnTime  * Fs;
        const double minOffDuty = minOffTime * Fs;
        // Buck if Vo/Vin sufficiently below 1 AND D_buck respects min on-time.
        if (M < 1.0 - hysteresisRatio) {
            const double D_buck = M;
            if (D_buck > minOnDuty && (1.0 - D_buck) > minOffDuty) {
                return Region::BUCK;
            }
        }
        // Boost if Vo/Vin sufficiently above 1 AND D_boost respects min on-time.
        if (M > 1.0 + hysteresisRatio) {
            const double D_boost = 1.0 - 1.0 / M;
            if (D_boost > minOnDuty && (1.0 - D_boost) > minOffDuty) {
                return Region::BOOST;
            }
        }
        return Region::BUCK_BOOST;
    }

    double FourSwitchBuckBoost::compute_buck_duty(double Vin, double Vo, double efficiency, double maximumDutyCycle) {
        if (Vin <= 0.0 || efficiency <= 0.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_buck_duty: Vin and efficiency must be > 0");
        }
        const double D = Vo / (Vin * efficiency);
        if (D >= 1.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_buck_duty: D >= 1 — converter cannot regulate (Vo too high for buck)");
        }
        const double dutyTolerance = 0.01;
        if (D >= maximumDutyCycle - dutyTolerance) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_buck_duty: D=" + std::to_string(D) +
                " exceeds maximumDutyCycle " + std::to_string(maximumDutyCycle) +
                " — reduce Vo, raise Vin, or raise maximumDutyCycle");
        }
        return D;
    }

    double FourSwitchBuckBoost::compute_boost_duty(double Vin, double Vo, double efficiency, double maximumDutyCycle) {
        if (Vo <= 0.0 || efficiency <= 0.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_boost_duty: Vo and efficiency must be > 0");
        }
        const double D = 1.0 - (Vin * efficiency) / Vo;
        if (D >= 1.0 || D <= 0.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_boost_duty: D out of range — converter cannot regulate (Vin too high for boost)");
        }
        const double dutyTolerance = 0.01;
        if (D >= maximumDutyCycle - dutyTolerance) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_boost_duty: D=" + std::to_string(D) +
                " exceeds maximumDutyCycle " + std::to_string(maximumDutyCycle) +
                " — reduce Vo, raise Vin, or raise maximumDutyCycle");
        }
        return D;
    }

    void FourSwitchBuckBoost::compute_buck_boost_duties(
        double Vin, double Vo, double efficiency, MAS::TransitionMode mode,
        double& D_buck_out, double& D_boost_out,
        double maximumDutyCycle)
    {
        if (Vin <= 0.0 || Vo <= 0.0 || efficiency <= 0.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_buck_boost_duties: Vin, Vo, η must be > 0");
        }
        if (mode == MAS::TransitionMode::SIMULTANEOUS) {
            // Mode 1: M = D/(1−D), single duty applied to Q1+Q3.
            const double Veff = Vin * efficiency;
            const double D = Vo / (Veff + Vo);
            D_buck_out  = D;
            D_boost_out = D;
        } else {
            // Mode 2 (SPLIT_PWM, LM5176/LT8390 default): two duties per LM5176
            // §8.3.4. Solve M = D_buck / (1−D_boost) with the controller
            // constraint D_buck + D_boost = 1 + δ (δ = 0.05 typical).
            constexpr double delta = 0.05;
            const double Veff = Vin * efficiency;
            const double M = Vo / Veff;
            // From D_buck = M·(1−D_boost) and D_buck + D_boost = 1 + δ:
            //   M·(1−D_boost) + D_boost = 1 + δ
            //   M − M·D_boost + D_boost = 1 + δ
            //   D_boost·(1 − M) = 1 + δ − M
            //   D_boost = (1 + δ − M) / (1 − M)         (M ≠ 1)
            // Special case M = 1: D_buck = D_boost = (1+δ)/2.
            if (std::abs(1.0 - M) < 1e-9) {
                D_buck_out  = (1.0 + delta) / 2.0;
                D_boost_out = (1.0 + delta) / 2.0;
            } else {
                D_boost_out = (1.0 + delta - M) / (1.0 - M);
                D_buck_out  = M * (1.0 - D_boost_out);
                // Throw on out-of-range duties (no silent clamp — would
                // hide a real sizing failure where Vin/Vo demand an
                // unrealisable duty pair). Boost branch is the one that
                // approaches 1; buck branch is bounded above by
                // (1+delta)/2 ≈ 0.525 in normal SPLIT_PWM operation.
                const double dutyTolerance = 0.01;
                if (D_boost_out >= maximumDutyCycle - dutyTolerance) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "FourSwitchBuckBoost::compute_buck_boost_duties: D_boost=" +
                        std::to_string(D_boost_out) +
                        " exceeds maximumDutyCycle " + std::to_string(maximumDutyCycle) +
                        " (SPLIT_PWM, M=" + std::to_string(M) +
                        ") — reduce Vo, raise Vin, or raise maximumDutyCycle");
                }
                if (D_boost_out <= 0.05) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "FourSwitchBuckBoost::compute_buck_boost_duties: D_boost=" +
                        std::to_string(D_boost_out) +
                        " <= 0.05 (SPLIT_PWM, M=" + std::to_string(M) +
                        ") — converter is in pure buck region, classify_region should not have selected BUCK_BOOST");
                }
                if (D_buck_out <= 0.05 || D_buck_out >= maximumDutyCycle - dutyTolerance) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "FourSwitchBuckBoost::compute_buck_boost_duties: D_buck=" +
                        std::to_string(D_buck_out) +
                        " out of [0.05, " + std::to_string(maximumDutyCycle - dutyTolerance) +
                        "] (SPLIT_PWM, M=" + std::to_string(M) + ")");
                }
            }
        }
    }

    double FourSwitchBuckBoost::compute_inductor_buck_region(
        double Vin, double Vo, double Iout, double Fs, double rippleRatio)
    {
        if (Vin <= 0 || Vo <= 0 || Iout <= 0 || Fs <= 0 || rippleRatio <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_inductor_buck_region: all inputs must be > 0");
        }
        if (Vo >= Vin) {
            // Region not valid; return 0 (caller picks max with boost branch).
            return 0.0;
        }
        // L = Vo · (Vin − Vo) / (K · Fs · Vin · Iout),  K = rippleRatio
        return Vo * (Vin - Vo) / (rippleRatio * Iout * Fs * Vin);
    }

    double FourSwitchBuckBoost::compute_inductor_boost_region(
        double Vin, double Vo, double Iout, double Fs, double rippleRatio)
    {
        if (Vin <= 0 || Vo <= 0 || Iout <= 0 || Fs <= 0 || rippleRatio <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_inductor_boost_region: all inputs must be > 0");
        }
        if (Vo <= Vin) {
            return 0.0;
        }
        // L = Vin² · (Vo − Vin) / (K · Fs · Iout · Vo²)
        return (Vin * Vin) * (Vo - Vin) / (rippleRatio * Iout * Fs * Vo * Vo);
    }

    double FourSwitchBuckBoost::compute_worst_case_inductance(
        double VinMin, double VinMax, double Vo, double Iout, double Fs, double rippleRatio)
    {
        const double L_buck  = compute_inductor_buck_region(VinMax, Vo, Iout, Fs, rippleRatio);
        const double L_boost = compute_inductor_boost_region(VinMin, Vo, Iout, Fs, rippleRatio);
        return std::max(L_buck, L_boost);
    }

    double FourSwitchBuckBoost::compute_inductor_dc_bias_boost(double Vin, double Vo, double Iout) {
        if (Vin <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost::compute_inductor_dc_bias_boost: Vin must be > 0");
        }
        return Iout * Vo / Vin;
    }

    // ============================================================
    //   Constructors
    // ============================================================

    FourSwitchBuckBoost::FourSwitchBuckBoost(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    AdvancedFourSwitchBuckBoost::AdvancedFourSwitchBuckBoost(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    // ============================================================
    //   Per-OP analytical worker
    // ============================================================

    OperatingPoint FourSwitchBuckBoost::process_operating_point_for_input_voltage(
        double inputVoltage, const TopologyExcitation& excitation, double inductance)
    {
        OperatingPoint operatingPoint;
        const double switchingFrequency = excitation.get_switching_frequency();
        const double outputVoltage = excitation.get_output_voltages()[0];
        // Bidirectional support: |Iout| sets magnitude; sign is preserved
        // for inductor-current waveform polarity (reverse power flow). The
        // analytical equations below operate on |Iout| then re-sign at the
        // end if get_bidirectional() is true and Iout < 0.
        const double Iout_signed = excitation.get_output_currents()[0];
        const double Iout = std::abs(Iout_signed);
        const bool bidir = get_bidirectional().value_or(false);
        const int currentSign = (bidir && Iout_signed < 0) ? -1 : 1;

        double efficiency = get_efficiency().value_or(1.0);
        const double hysteresis = get_transition_hysteresis_ratio().value_or(transitionHysteresisDefault);
        const MAS::TransitionMode tmode = get_transition_mode().value_or(MAS::TransitionMode::SPLIT_PWM);
        const int phaseCount = static_cast<int>(get_phase_count().value_or(int64_t(1)));
        if (phaseCount < 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "FourSwitchBuckBoost: phaseCount must be >= 1");
        }
        // Per-phase inductor sees Iout/N.
        const double Iout_per_phase = Iout / static_cast<double>(phaseCount);

        const Region region = classify_region(inputVoltage, outputVoltage, hysteresis,
                                              minOnTime, minOffTime, switchingFrequency);
        lastRegion = region;

        double D_buck = 0.0, D_boost = 0.0;
        double iL_avg = 0.0;       // per-phase inductor avg current
        double dIL    = 0.0;       // peak-to-peak ripple
        double iL_peak = 0.0;
        double iL_rms  = 0.0;
        double M = 0.0;
        double iIn_avg = 0.0;
        double q1Rms = 0.0, q2Rms = 0.0, q3Rms = 0.0, q4Rms = 0.0;

        // Helper: triangular RMS = sqrt(I_avg² + (ΔI/2)²/3) for a continuous
        // triangular ramp; for a "switch-on duty D" pulse of triangular
        // current with avg I_avg over the on-window, RMS over full period
        // is sqrt(D · (I_avg² + (ΔI/2)²/3)).
        auto switchRms = [](double duty, double i_avg_on, double dI) {
            if (duty <= 0.0) return 0.0;
            return std::sqrt(duty * (i_avg_on * i_avg_on + (dI * dI) / 12.0));
        };

        const double period = 1.0 / switchingFrequency;
        double tOn = 0.0;            // primary PWM on-time

        // --- Region-dispatched analytical solver ---
        Waveform currentWaveform;
        Waveform voltageWaveform;
        WaveformLabel currentLabel = WaveformLabel::TRIANGULAR;
        WaveformLabel voltageLabel = WaveformLabel::RECTANGULAR;
        double primaryVoltagePtp = 0.0;
        double dutyForWaveform = 0.5;

        if (region == Region::BUCK) {
            D_buck = compute_buck_duty(inputVoltage, outputVoltage, efficiency, maximumDutyCycle.value_or(0.95));
            D_boost = 0.0;
            tOn = D_buck * period;
            M = D_buck;
            iL_avg = Iout_per_phase;
            dIL = (inputVoltage - outputVoltage) * tOn / inductance;
            iL_peak = iL_avg + dIL / 2.0;
            iL_rms  = std::sqrt(iL_avg * iL_avg + (dIL * dIL) / 12.0);
            iIn_avg = iL_avg * D_buck;     // input current = duty-weighted L current
            // Switches: Q1 PWM(D_buck), Q2 PWM(1−D_buck), Q3 ON, Q4 OFF.
            q1Rms = switchRms(D_buck,        iL_avg, dIL);
            q2Rms = switchRms(1.0 - D_buck,  iL_avg, dIL);
            q3Rms = std::sqrt(iL_avg * iL_avg + (dIL * dIL) / 12.0);
            q4Rms = 0.0;

            primaryVoltagePtp = inputVoltage;          // L sees +Vin−Vo / −Vo, pp = Vin
            dutyForWaveform   = D_buck;
            currentLabel = WaveformLabel::TRIANGULAR;
            voltageLabel = WaveformLabel::RECTANGULAR;
        }
        else if (region == Region::BOOST) {
            D_boost = compute_boost_duty(inputVoltage, outputVoltage, efficiency, maximumDutyCycle.value_or(0.95));
            D_buck = 1.0;
            tOn = D_boost * period;
            M = 1.0 / (1.0 - D_boost);
            iL_avg = Iout_per_phase / (1.0 - D_boost);     // = Iout · Vo/Vin
            dIL = inputVoltage * tOn / inductance;
            iL_peak = iL_avg + dIL / 2.0;
            iL_rms  = std::sqrt(iL_avg * iL_avg + (dIL * dIL) / 12.0);
            iIn_avg = iL_avg;                              // L current = source current always
            // Switches: Q1 ON, Q2 OFF, Q3 PWM(1−D_boost), Q4 PWM(D_boost).
            q1Rms = std::sqrt(iL_avg * iL_avg + (dIL * dIL) / 12.0);
            q2Rms = 0.0;
            q3Rms = switchRms(1.0 - D_boost, iL_avg, dIL);
            q4Rms = switchRms(D_boost,       iL_avg, dIL);

            primaryVoltagePtp = outputVoltage;       // L sees +Vin / -(Vo−Vin), pp ≈ Vo
            dutyForWaveform   = D_boost;
            currentLabel = WaveformLabel::TRIANGULAR;
            voltageLabel = WaveformLabel::RECTANGULAR;
        }
        else { // BUCK_BOOST
            compute_buck_boost_duties(inputVoltage, outputVoltage, efficiency, tmode,
                                      D_buck, D_boost, maximumDutyCycle.value_or(0.95));
            if (tmode == MAS::TransitionMode::SIMULTANEOUS) {
                // Mode 1: D = Vo/(Vin+Vo); on-time both legs ON.
                const double D = D_buck;
                tOn = D * period;
                M = D / (1.0 - D);
                iL_avg = Iout_per_phase / (1.0 - D);
                dIL = inputVoltage * tOn / inductance;
            } else {
                // Mode 2 split-PWM. Approximate ΔI as the boost-region ripple
                // over the (1−D_boost) interval where Q3 is ON. The actual
                // four-sub-interval waveform is complex; tier-3 v1 model
                // captures the dominant first-order behaviour.
                M = D_buck / std::max(1.0 - D_boost, 1e-6);
                iL_avg = Iout_per_phase * std::max(1.0, outputVoltage / inputVoltage);
                tOn = D_buck * period;
                dIL = (inputVoltage * D_buck * period) / inductance;
            }
            iL_peak = iL_avg + dIL / 2.0;
            iL_rms  = std::sqrt(iL_avg * iL_avg + (dIL * dIL) / 12.0);
            iIn_avg = iL_avg * D_buck;
            q1Rms = switchRms(D_buck,        iL_avg, dIL);
            q2Rms = switchRms(1.0 - D_buck,  iL_avg, dIL);
            q3Rms = switchRms(1.0 - D_boost, iL_avg, dIL);
            q4Rms = switchRms(D_boost,       iL_avg, dIL);

            primaryVoltagePtp = inputVoltage + outputVoltage;  // L sees +Vin / −Vo
            dutyForWaveform   = D_buck;
            currentLabel = WaveformLabel::TRIANGULAR;
            voltageLabel = WaveformLabel::RECTANGULAR;
        }

        // --- DCM check ---
        const double R_load = std::max(outputVoltage / std::max(Iout, 1e-6), 1e-3);
        const double K = 2.0 * inductance * switchingFrequency / R_load;
        double K_crit = 0.0;
        switch (region) {
            case Region::BUCK:       K_crit = 1.0 - D_buck; break;
            case Region::BOOST:      K_crit = D_boost * (1.0 - D_boost) * (1.0 - D_boost); break;
            case Region::BUCK_BOOST: { const double D = D_buck; K_crit = (1.0 - D) * (1.0 - D); break; }
        }
        const bool isCcm = (K > K_crit);

        // --- Output cap sizing (output ripple from cap charge in BOOST: Iout·D_boost/(Co·Fs)) ---
        const double rippleRatioVo = get_output_voltage_ripple_ratio().value_or(outputRippleRatioDefault);
        const double dVo_target = std::max(rippleRatioVo * outputVoltage, 1e-3);
        // Use the larger of buck (ΔI_L/(8·Co·Fs)) and boost (Iout·D_boost/(Co·Fs)) constraints.
        const double Co_buck  = (region == Region::BUCK)
                              ? dIL / (8.0 * dVo_target * switchingFrequency)
                              : 0.0;
        const double Co_boost = (D_boost > 0.0)
                              ? Iout * D_boost / (dVo_target * switchingFrequency)
                              : 0.0;
        const double Co_sized = std::max({Co_buck, Co_boost, 1e-6});

        // --- Input cap sizing (input ripple from buck pulsed input: Iout·D·(1−D)/(Cin·Fs·rippleRatio)) ---
        const double rippleRatioVin = 0.05;     // 5% target Vin ripple
        const double dVin_target = std::max(rippleRatioVin * inputVoltage, 1e-3);
        const double Cin_buck = (region == Region::BUCK)
                              ? Iout * D_buck * (1.0 - D_buck) / (dVin_target * switchingFrequency)
                              : 0.0;
        const double Cin_sized = std::max({Cin_buck, 1e-6});

        // --- Diagnostics ---
        lastDutyCycleBuck          = D_buck;
        lastDutyCycleBoost         = D_boost;
        lastConversionRatio        = M;
        lastInductorAverageCurrent = iL_avg;
        lastInductorPeakToPeak     = dIL;
        lastInductorPeakCurrent    = iL_peak;
        lastInductorRmsCurrent     = iL_rms;
        lastInputCurrentAverage    = iIn_avg;
        lastInputCurrentRipple     = dIL;     // approximate
        lastQ1RmsCurrent           = q1Rms;
        lastQ2RmsCurrent           = q2Rms;
        lastQ3RmsCurrent           = q3Rms;
        lastQ4RmsCurrent           = q4Rms;
        lastDcmK                   = K;
        lastDcmKcrit               = K_crit;
        lastIsCcm                  = isCcm;
        lastSizedInductance        = inductance;
        lastSizedInputCapacitance  = Cin_sized;
        lastSizedOutputCapacitance = Co_sized;
        lastOutputVoltageRipple    = (Co_sized > 0)
                                     ? Iout * std::max(D_boost, 0.0) / (Co_sized * switchingFrequency)
                                     : 0.0;

        // --- Inductor winding waveforms ---
        currentWaveform = Inputs::create_waveform(
            currentLabel, dIL, switchingFrequency, dutyForWaveform,
            currentSign * iL_avg, 0);
        voltageWaveform = Inputs::create_waveform(
            voltageLabel, primaryVoltagePtp, switchingFrequency, dutyForWaveform, 0, 0);

        auto excitationOut = complete_excitation(currentWaveform, voltageWaveform,
                                                  switchingFrequency, "Inductor");
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitationOut);

        // Capture extra-component waveforms for downstream get_extra_components_inputs().
        extraInductorVoltageWaveforms.push_back(voltageWaveform);
        extraInductorCurrentWaveforms.push_back(currentWaveform);

        // Input cap sees the difference between source DC and switch-pulled current:
        // approximate as a rectangular ripple.
        Waveform vCinWf = Inputs::create_waveform(
            WaveformLabel::RECTANGULAR, dVin_target, switchingFrequency, dutyForWaveform,
            inputVoltage, 0);
        Waveform iCinWf = Inputs::create_waveform(
            WaveformLabel::RECTANGULAR, iL_avg, switchingFrequency, dutyForWaveform, 0, 0);
        extraInputCapVoltageWaveforms.push_back(vCinWf);
        extraInputCapCurrentWaveforms.push_back(iCinWf);

        // Output cap sees the rectifier-side pulsed current (largest in BOOST).
        Waveform vCoWf = Inputs::create_waveform(
            WaveformLabel::RECTANGULAR, lastOutputVoltageRipple, switchingFrequency,
            dutyForWaveform, outputVoltage, 0);
        Waveform iCoWf = Inputs::create_waveform(
            WaveformLabel::RECTANGULAR, iL_avg, switchingFrequency, dutyForWaveform, -Iout, 0);
        extraOutputCapVoltageWaveforms.push_back(vCoWf);
        extraOutputCapCurrentWaveforms.push_back(iCoWf);

        OperatingConditions conditions;
        conditions.set_ambient_temperature(excitation.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    // ============================================================
    //   run_checks / process_design_requirements
    // ============================================================

    bool FourSwitchBuckBoost::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) return false;
            throw InvalidInputException(ErrorCode::MISSING_DATA,
                "FourSwitchBuckBoost: at least one operating point is needed");
        }
        if (!get_input_voltage().get_nominal() && !get_input_voltage().get_maximum()
            && !get_input_voltage().get_minimum())
        {
            if (!assert) return false;
            throw InvalidInputException(ErrorCode::MISSING_DATA,
                "FourSwitchBuckBoost: no input voltage introduced");
        }
        return true;
    }

    DesignRequirements FourSwitchBuckBoost::process_design_requirements() {
        const double VinMin = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        const double VinMax = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument(
                "FourSwitchBuckBoost::process_design_requirements: missing both "
                "currentRippleRatio and maximumSwitchCurrent");
        }
        const double rippleRatio = get_current_ripple_ratio().value_or(inductorRippleRatioDefault);
        const int phaseCount = static_cast<int>(get_phase_count().value_or(int64_t(1)));

        double maxNeededL = 0.0;
        for (const auto& op : get_operating_points()) {
            const double Vo = op.get_output_voltages()[0];
            const double Iout = std::abs(op.get_output_currents()[0]) / static_cast<double>(phaseCount);
            const double Fs = op.get_switching_frequency();
            const double L = compute_worst_case_inductance(VinMin, VinMax, Vo, Iout, Fs, rippleRatio);
            maxNeededL = std::max(maxNeededL, L);
        }
        if (maxNeededL <= 0.0) {
            throw std::invalid_argument(
                "FourSwitchBuckBoost::process_design_requirements: derived inductance is non-positive");
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        DimensionWithTolerance Lwt;
        Lwt.set_minimum(roundFloat(maxNeededL, 10));
        designRequirements.set_magnetizing_inductance(Lwt);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::FOUR_SWITCH_BUCK_BOOST_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> FourSwitchBuckBoost::process_operating_points(
        const std::vector<double>& /*turnsRatios*/, double magnetizingInductance)
    {
        extraInductorVoltageWaveforms.clear();
        extraInductorCurrentWaveforms.clear();
        extraInputCapVoltageWaveforms.clear();
        extraInputCapCurrentWaveforms.clear();
        extraOutputCapVoltageWaveforms.clear();
        extraOutputCapCurrentWaveforms.clear();

        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        for (size_t i = 0; i < inputVoltages.size(); ++i) {
            for (size_t opIdx = 0; opIdx < get_operating_points().size(); ++opIdx) {
                auto op = process_operating_point_for_input_voltage(
                    inputVoltages[i], get_operating_points()[opIdx], magnetizingInductance);
                std::string name = inputVoltagesNames[i] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(opIdx);
                }
                op.set_name(name);
                operatingPoints.push_back(op);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> FourSwitchBuckBoost::process_operating_points(Magnetic magnetic) {
        run_checks(_assertErrors);
        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel
            .calculate_inductance_from_number_turns_and_gapping(
                magnetic.get_mutable_core(), magnetic.get_mutable_coil())
            .get_magnetizing_inductance().get_nominal().value();
        return process_operating_points(magnetic.get_turns_ratios(), magnetizingInductance);
    }

    Inputs AdvancedFourSwitchBuckBoost::process() {
        FourSwitchBuckBoost::run_checks(_assertErrors);
        Inputs inputs;
        const double L = get_desired_inductance();

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        inputs.set_design_requirements(process_design_requirements());

        for (size_t i = 0; i < inputVoltages.size(); ++i) {
            for (size_t opIdx = 0; opIdx < get_operating_points().size(); ++opIdx) {
                auto op = process_operating_point_for_input_voltage(
                    inputVoltages[i], get_operating_points()[opIdx], L);
                std::string name = inputVoltagesNames[i] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(opIdx);
                }
                op.set_name(name);
                inputs.get_mutable_operating_points().push_back(op);
            }
        }
        return inputs;
    }

    DesignRequirements AdvancedFourSwitchBuckBoost::process_design_requirements() {
        // Issue M1: build DR directly from desiredInductance. Do NOT chain
        // to FourSwitchBuckBoost::process_design_requirements() — the parent
        // sizes inductance from current ripple, defeating the "Advanced"
        // override semantics.
        const double L = get_desired_inductance();
        if (L <= 0.0) {
            throw std::invalid_argument(
                "AdvancedFourSwitchBuckBoost::process_design_requirements: desiredInductance must be > 0");
        }

        DesignRequirements dr;
        dr.get_mutable_turns_ratios().clear();
        DimensionWithTolerance Lwt;
        Lwt.set_nominal(roundFloat(L, 10));
        dr.set_magnetizing_inductance(Lwt);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        dr.set_isolation_sides(isolationSides);
        dr.set_topology(Topologies::FOUR_SWITCH_BUCK_BOOST_CONVERTER);
        return dr;
    }

    // ============================================================
    //   get_extra_components_inputs — L (single inductor), Cin, Co
    // ============================================================

    std::vector<std::variant<Inputs, CAS::Inputs>> FourSwitchBuckBoost::get_extra_components_inputs(
        ExtraComponentsMode mode, std::optional<Magnetic> magnetic)
    {
        if (mode == ExtraComponentsMode::REAL && !magnetic.has_value()) {
            throw std::invalid_argument(
                "FourSwitchBuckBoost::get_extra_components_inputs: mode REAL requires a designed magnetic");
        }
        if (lastSizedInductance <= 0.0 || extraInductorVoltageWaveforms.empty()) {
            throw std::runtime_error(
                "FourSwitchBuckBoost::get_extra_components_inputs: call process_operating_points() first");
        }

        const double fsw = get_operating_points()[0].get_switching_frequency();
        std::vector<std::variant<Inputs, CAS::Inputs>> result;

        // L (single inductor, 1 winding).
        {
            Inputs masInputs;
            DesignRequirements dr;
            DimensionWithTolerance L;
            L.set_nominal(lastSizedInductance);
            dr.set_magnetizing_inductance(L);
            dr.set_name("inductor");
            dr.set_topology(Topologies::FOUR_SWITCH_BUCK_BOOST_CONVERTER);
            dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
            masInputs.set_design_requirements(dr);

            std::vector<OperatingPoint> masOps;
            for (size_t i = 0; i < extraInductorVoltageWaveforms.size(); ++i) {
                OperatingPoint op;
                auto exc = complete_excitation(extraInductorCurrentWaveforms[i],
                                               extraInductorVoltageWaveforms[i], fsw, "Inductor");
                op.get_mutable_excitations_per_winding().push_back(exc);
                masOps.push_back(op);
            }
            masInputs.set_operating_points(masOps);
            result.emplace_back(std::move(masInputs));
        }

        // Cin (input capacitor).
        {
            CAS::Inputs casInputs;
            CAS::DesignRequirements dr;
            CAS::DimensionWithTolerance C;
            C.set_nominal(lastSizedInputCapacitance);
            dr.set_capacitance(C);
            double peakV = 0.0;
            for (const auto& wf : extraInputCapVoltageWaveforms)
                for (double v : wf.get_data()) peakV = std::max(peakV, std::abs(v));
            dr.set_rated_voltage(peakV * 1.2);
            dr.set_role(CAS::Application::INPUT_FILTER);
            dr.set_name("inputCapacitor");
            casInputs.set_design_requirements(dr);

            std::vector<CAS::TwoTerminalOperatingPoint> casOps;
            for (size_t i = 0; i < extraInputCapVoltageWaveforms.size(); ++i) {
                CAS::TwoTerminalOperatingPoint op;
                CAS::OperatingPointExcitation exc;
                exc.set_frequency(fsw);
                CAS::SignalDescriptor vSig;
                CAS::Waveform vWfm;
                vWfm.set_data(extraInputCapVoltageWaveforms[i].get_data());
                if (extraInputCapVoltageWaveforms[i].get_time())
                    vWfm.set_time(extraInputCapVoltageWaveforms[i].get_time().value());
                vSig.set_waveform(vWfm);
                exc.set_voltage(vSig);
                CAS::SignalDescriptor iSig;
                CAS::Waveform iWfm;
                iWfm.set_data(extraInputCapCurrentWaveforms[i].get_data());
                if (extraInputCapCurrentWaveforms[i].get_time())
                    iWfm.set_time(extraInputCapCurrentWaveforms[i].get_time().value());
                iSig.set_waveform(iWfm);
                exc.set_current(iSig);
                op.set_excitation(exc);
                casOps.push_back(op);
            }
            casInputs.set_operating_points(casOps);
            result.emplace_back(std::move(casInputs));
        }

        // Co (output capacitor).
        {
            CAS::Inputs casInputs;
            CAS::DesignRequirements dr;
            CAS::DimensionWithTolerance C;
            C.set_nominal(lastSizedOutputCapacitance);
            dr.set_capacitance(C);
            double peakV = 0.0;
            for (const auto& wf : extraOutputCapVoltageWaveforms)
                for (double v : wf.get_data()) peakV = std::max(peakV, std::abs(v));
            dr.set_rated_voltage(peakV * 1.2);
            dr.set_role(CAS::Application::OUTPUT_FILTER);
            dr.set_name("outputCapacitor");
            casInputs.set_design_requirements(dr);

            std::vector<CAS::TwoTerminalOperatingPoint> casOps;
            for (size_t i = 0; i < extraOutputCapVoltageWaveforms.size(); ++i) {
                CAS::TwoTerminalOperatingPoint op;
                CAS::OperatingPointExcitation exc;
                exc.set_frequency(fsw);
                CAS::SignalDescriptor vSig;
                CAS::Waveform vWfm;
                vWfm.set_data(extraOutputCapVoltageWaveforms[i].get_data());
                if (extraOutputCapVoltageWaveforms[i].get_time())
                    vWfm.set_time(extraOutputCapVoltageWaveforms[i].get_time().value());
                vSig.set_waveform(vWfm);
                exc.set_voltage(vSig);
                CAS::SignalDescriptor iSig;
                CAS::Waveform iWfm;
                iWfm.set_data(extraOutputCapCurrentWaveforms[i].get_data());
                if (extraOutputCapCurrentWaveforms[i].get_time())
                    iWfm.set_time(extraOutputCapCurrentWaveforms[i].get_time().value());
                iSig.set_waveform(iWfm);
                exc.set_current(iSig);
                op.set_excitation(exc);
                casOps.push_back(op);
            }
            casInputs.set_operating_points(casOps);
            result.emplace_back(std::move(casInputs));
        }

        return result;
    }

    // ============================================================
    //   SPICE — region-aware netlist
    // ============================================================
    //
    // Node convention:
    //   vin_dc → Vin_sense → vin_p
    //   vin_p → Cin → 0
    //   vin_p → S_Q1 (HS) → sw1
    //   sw1   → S_Q2 (LS) → 0
    //   sw1   → Vl_sense → l_in → L → l_out → Vsw2_sense → sw2
    //   sw2   → S_Q3 (HS) → vout
    //   sw2   → S_Q4 (LS) → 0
    //   vout  → Co → 0;  vout → Vout_sense → out_load → Rload → 0
    //
    // Region-aware gate emission:
    //   BUCK:        Vpwm_q1 PULSE(D_buck), Vpwm_q2 PULSE(1−D_buck)
    //                Vpwm_q3 = +pwmHigh DC, Vpwm_q4 = 0 DC
    //   BOOST:       Vpwm_q1 = +pwmHigh DC, Vpwm_q2 = 0 DC
    //                Vpwm_q3 PULSE(1−D_boost), Vpwm_q4 PULSE(D_boost)
    //   BUCK_BOOST:  all four PWM with appropriate phase alignment
    //
    // RC snubbers on both SW nodes; IC pre-charge on Co.
    // ============================================================

    std::string FourSwitchBuckBoost::generate_ngspice_circuit(
        double inductance, size_t inputVoltageIndex, size_t operatingPointIndex)
    {
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames_;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames_);
        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("FourSwitchBuckBoost::generate_ngspice_circuit: inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("FourSwitchBuckBoost::generate_ngspice_circuit: operatingPointIndex out of range");
        }

        const double inputVoltage = inputVoltages[inputVoltageIndex];
        const auto opPoint = get_operating_points()[operatingPointIndex];
        const double outputVoltage = opPoint.get_output_voltages()[0];
        const double outputCurrent = std::abs(opPoint.get_output_currents()[0]);
        if (outputCurrent <= 0.0) {
            throw std::invalid_argument(
                "FourSwitchBuckBoost::generate_ngspice_circuit: outputCurrent must be > 0 "
                "(received " + std::to_string(outputCurrent) +
                "); cannot derive load resistance.");
        }
        const double switchingFrequency = opPoint.get_switching_frequency();
        const double efficiency = get_efficiency().value_or(1.0);
        const double hysteresis = get_transition_hysteresis_ratio().value_or(transitionHysteresisDefault);
        const MAS::TransitionMode tmode = get_transition_mode().value_or(MAS::TransitionMode::SPLIT_PWM);

        const Region region = classify_region(inputVoltage, outputVoltage, hysteresis,
                                              minOnTime, minOffTime, switchingFrequency);
        double D_buck = 0.0, D_boost = 0.0;
        if (region == Region::BUCK) {
            D_buck = compute_buck_duty(inputVoltage, outputVoltage, efficiency, maximumDutyCycle.value_or(0.95));
        } else if (region == Region::BOOST) {
            D_boost = compute_boost_duty(inputVoltage, outputVoltage, efficiency, maximumDutyCycle.value_or(0.95));
        } else {
            compute_buck_boost_duties(inputVoltage, outputVoltage, efficiency, tmode, D_buck, D_boost, maximumDutyCycle.value_or(0.95));
        }

        const SpiceSimulationConfig cfg = spice_config();
        const double loadResistance = std::max(outputVoltage / outputCurrent, 1e-3);
        const double period = 1.0 / switchingFrequency;
        const int periodsToExtract = get_num_periods_to_extract();
        const int settlingPeriods  = get_num_steady_state_periods();
        const int numPeriodsTotal  = settlingPeriods + periodsToExtract;
        const double simTime = numPeriodsTotal * period;
        const double startTime = settlingPeriods * period;
        const double stepTime = period / cfg.samplesPerPeriod;
        const double tOnBuck  = D_buck  * period;
        const double tOnBoost = D_boost * period;

        std::ostringstream c;
        c << "* 4-Switch Buck-Boost Converter - Generated by OpenMagnetics\n";
        c << "* Vin=" << inputVoltage << "V, Vo=" << outputVoltage << "V, Iout=" << outputCurrent << "A\n";
        c << "* Fs=" << (switchingFrequency / 1e3) << "kHz, L=" << (inductance * 1e6) << "uH, "
          << "region=" << (region == Region::BUCK ? "BUCK"
                            : (region == Region::BOOST ? "BOOST" : "BUCK_BOOST"))
          << ", D_buck=" << D_buck << ", D_boost=" << D_boost << "\n";
        c << "* transitionMode=" << (tmode == MAS::TransitionMode::SIMULTANEOUS ? "simultaneous" : "splitPwm")
          << ", controlMode=" << static_cast<int>(get_control_mode().value_or(MAS::ControlMode::PEAK_CURRENT)) << "\n\n";

        // DC input
        c << "* DC Input + sense\n";
        c << "Vin vin_dc 0 " << inputVoltage << "\n";
        c << "Vin_sense vin_dc vin_p 0\n\n";

        // Models
        c << "* Switch model\n";
        c << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
          << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n\n";

        // PWM gate drives — region-aware emission per FOUR_SWITCH_BUCK_BOOST_PLAN.md §A.6.
        c << "* Region-aware gate drives (region=" << (region == Region::BUCK ? "BUCK"
                            : (region == Region::BOOST ? "BOOST" : "BUCK_BOOST")) << ")\n";
        auto emitPulse = [&](const std::string& name, const std::string& node, double tOn) {
            c << name << " " << node << " 0 PULSE(0 " << cfg.pwmHigh << " 0 "
              << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
              << tOn << " " << period << std::fixed << ")\n";
        };
        auto emitDcHigh = [&](const std::string& name, const std::string& node) {
            c << name << " " << node << " 0 " << cfg.pwmHigh << "\n";
        };
        auto emitDcLow = [&](const std::string& name, const std::string& node) {
            c << name << " " << node << " 0 0\n";
        };
        if (region == Region::BUCK) {
            emitPulse("Vpwm_q1", "pwm_q1", tOnBuck);
            // Q2 anti-phase: avoid shoot-through with deadtime; for ideal SW
            // we simply invert (1-D_buck) on-time starting at tOnBuck.
            c << "Vpwm_q2 pwm_q2 0 PULSE(" << cfg.pwmHigh << " 0 0 "
              << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
              << tOnBuck << " " << period << std::fixed << ")\n";
            emitDcHigh("Vpwm_q3", "pwm_q3");
            emitDcLow("Vpwm_q4", "pwm_q4");
        } else if (region == Region::BOOST) {
            emitDcHigh("Vpwm_q1", "pwm_q1");
            emitDcLow("Vpwm_q2", "pwm_q2");
            // Q3 anti-phase to Q4
            c << "Vpwm_q3 pwm_q3 0 PULSE(" << cfg.pwmHigh << " 0 0 "
              << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
              << tOnBoost << " " << period << std::fixed << ")\n";
            emitPulse("Vpwm_q4", "pwm_q4", tOnBoost);
        } else { // BUCK_BOOST
            emitPulse("Vpwm_q1", "pwm_q1", tOnBuck);
            c << "Vpwm_q2 pwm_q2 0 PULSE(" << cfg.pwmHigh << " 0 0 "
              << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
              << tOnBuck << " " << period << std::fixed << ")\n";
            c << "Vpwm_q3 pwm_q3 0 PULSE(" << cfg.pwmHigh << " 0 0 "
              << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
              << tOnBoost << " " << period << std::fixed << ")\n";
            emitPulse("Vpwm_q4", "pwm_q4", tOnBoost);
        }
        c << "\n";

        // Input capacitor
        const double Cin_val = std::max(lastSizedInputCapacitance, 1e-6);
        c << "* Input capacitor\n";
        c << "Cin vin_p 0 " << std::scientific << Cin_val << std::fixed
          << " IC=" << inputVoltage << "\n\n";

        // Buck-leg (Q1 HS, Q2 LS) — node sw1
        c << "* Buck-leg switches (Q1 HS, Q2 LS) + RC snubber on SW1\n";
        c << "S_Q1 vin_p sw1 pwm_q1 0 SW1\n";
        c << "S_Q2 sw1 0     pwm_q2 0 SW1\n";
        c << "Rsnub_sw1 sw1 sn_sw1 " << cfg.snubR << "\n";
        c << "Csnub_sw1 sn_sw1 0 " << std::scientific << cfg.snubC << std::fixed << "\n\n";

        // Inductor with current sense
        c << "* Inductor with current sense (Vl_sense in series at the SW1 side)\n";
        c << "Vl_sense sw1 l_in 0\n";
        c << "L1 l_in l_out " << std::scientific << inductance << std::fixed
          << " IC=" << lastInductorAverageCurrent << "\n";
        c << "Bvpri_diff vpri_diff 0 V=V(l_in)-V(l_out)\n\n";

        // Boost-leg (Q3 HS, Q4 LS) — node sw2
        c << "* Vsw2_sense sense + boost-leg switches (Q3 HS to vout, Q4 LS to 0) + snubber\n";
        c << "Vsw2_sense l_out sw2 0\n";
        c << "S_Q3 sw2 vout pwm_q3 0 SW1\n";
        c << "S_Q4 sw2 0    pwm_q4 0 SW1\n";
        c << "Rsnub_sw2 sw2 sn_sw2 " << cfg.snubR << "\n";
        c << "Csnub_sw2 sn_sw2 0 " << std::scientific << cfg.snubC << std::fixed << "\n\n";

        // Output cap and load
        c << "* Output capacitor + resistive load with Vout_sense\n";
        c << "Cout vout 0 " << std::scientific << cfg.outputCapacitance << std::fixed
          << " IC=" << (outputVoltage * cfg.outputCapInitialChargeFraction) << "\n";
        c << "Vout_sense vout out_load 0\n";
        c << "Rload out_load 0 " << loadResistance << "\n\n";

        // Transient
        c << "* Transient analysis (UIC enables .ic)\n";
        c << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime
          << std::fixed << " UIC\n\n";

        c << "* Output signals\n";
        c << ".save v(vpri_diff) v(vin_dc) v(sw1) v(sw2) v(vout) "
             "i(Vin_sense) i(Vl_sense) i(Vsw2_sense) i(Vout_sense)\n\n";

        c << ".options RELTOL=" << cfg.relTol
          << " ABSTOL=" << std::scientific << cfg.absTol
          << " VNTOL=" << cfg.vnTol << std::fixed
          << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        c << ".options METHOD=" << cfg.method << " TRTOL=" << cfg.trTol << "\n";
        c << ".ic v(vout)=" << outputVoltage << "\n\n";
        c << ".end\n";
        return c.str();
    }

    std::vector<OperatingPoint> FourSwitchBuckBoost::simulate_and_extract_operating_points(double inductance) {
        std::vector<OperatingPoint> operatingPoints;
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("FourSwitchBuckBoost::simulate_and_extract_operating_points: ngspice not available");
        }

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        for (size_t i = 0; i < inputVoltages.size(); ++i) {
            for (size_t opIdx = 0; opIdx < get_operating_points().size(); ++opIdx) {
                auto opPoint = get_operating_points()[opIdx];
                std::string netlist = generate_ngspice_circuit(inductance, i, opIdx);
                double switchingFrequency = opPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);
                if (!simResult.success) {
                    throw std::runtime_error("FSBB SPICE simulation failed: " + simResult.errorMessage);
                }

                NgspiceRunner::WaveformNameMapping mapping = {
                    {{"voltage", "vpri_diff"}, {"current", "vl_sense#branch"}}
                };
                std::vector<std::string> windingNames = {"Inductor"};
                std::vector<bool> flipCurrentSign = {false};

                OperatingPoint op = NgspiceRunner::extract_operating_point(
                    simResult, mapping, switchingFrequency, windingNames,
                    opPoint.get_ambient_temperature(), flipCurrentSign);
                std::string name = inputVoltagesNames[i] + " input volt. (simulated)";
                if (get_operating_points().size() > 1) {
                    name += " op. point " + std::to_string(opIdx);
                }
                op.set_name(name);
                operatingPoints.push_back(op);
            }
        }
        return operatingPoints;
    }

    std::vector<ConverterWaveforms> FourSwitchBuckBoost::simulate_and_extract_topology_waveforms(
        double inductance, size_t numberOfPeriods)
    {
        std::vector<ConverterWaveforms> results;
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("FourSwitchBuckBoost::simulate_and_extract_topology_waveforms: ngspice not available");
        }
        const int originalNumPeriodsToExtract = get_num_periods_to_extract();
        set_num_periods_to_extract(static_cast<int>(numberOfPeriods));

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        for (size_t i = 0; i < inputVoltages.size(); ++i) {
            for (size_t opIdx = 0; opIdx < get_operating_points().size(); ++opIdx) {
                auto opPoint = get_operating_points()[opIdx];
                std::string netlist = generate_ngspice_circuit(inductance, i, opIdx);
                double switchingFrequency = opPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = numberOfPeriods;
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);
                if (!simResult.success) {
                    throw std::runtime_error("FSBB SPICE simulation failed: " + simResult.errorMessage);
                }

                std::map<std::string, size_t> nameToIndex;
                for (size_t k = 0; k < simResult.waveformNames.size(); ++k) {
                    std::string lower = simResult.waveformNames[k];
                    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                    nameToIndex[lower] = k;
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
                std::string name = inputVoltagesNames[i] + " input";
                if (get_operating_points().size() > 1) {
                    name += " op. point " + std::to_string(opIdx);
                }
                wf.set_operating_point_name(name);
                wf.set_input_voltage(getWaveform("vin_dc"));
                wf.set_input_current(getWaveform("vin_sense#branch"));
                wf.get_mutable_output_voltages().push_back(getWaveform("vout"));
                wf.get_mutable_output_currents().push_back(getWaveform("vout_sense#branch"));
                results.push_back(wf);
            }
        }

        set_num_periods_to_extract(originalNumPeriodsToExtract);
        return results;
    }

} // namespace OpenMagnetics
