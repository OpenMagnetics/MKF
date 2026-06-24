#include "converter_models/Cuk.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/NgspiceRunner.h"
#include "support/Utils.h"
#include <cfloat>
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include "support/Exceptions.h"
#include "converter_models/Topology.h"

namespace OpenMagnetics {

    // ============================================================
    //   Static analytical helpers (CUK_PLAN.md §2.13, §4.1–§4.4)
    // ============================================================

    double Cuk::calculate_duty_cycle(double inputVoltage, double outputVoltageMagnitude, double diodeVoltageDrop, double efficiency, double turnsRatio, double maximumDutyCycle) {
        // Cuk CCM, ideal:
        //   V1/V2 (n=1):     M(D) = -D/(1-D)
        //   V3 isolated:     M(D) = -D / ((1-D) · n),  n = Np/Ns (Flyback convention)
        //
        // Solving |Vo| · (1-D) · n = (Vin·η - Vd·(1-D)) · D for D and treating
        // the small Vd · (1-D) term as a forward-drop correction:
        //   D = (|Vo|+Vd) · n / (Vin·η + (|Vo|+Vd) · n)
        // For n=1 reduces to V1's expression (verified algebraically).
        if (inputVoltage <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_duty_cycle: input voltage must be > 0");
        }
        if (outputVoltageMagnitude <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_duty_cycle: |Vo| must be > 0");
        }
        if (turnsRatio <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_duty_cycle: turnsRatio must be > 0");
        }
        double effVin = inputVoltage * efficiency;
        double reflectedVo = (outputVoltageMagnitude + diodeVoltageDrop) * turnsRatio;
        double dutyCycle = reflectedVo / (reflectedVo + effVin);
        // Explicit configurable maxD gate (no silent clamp). 1 % tolerance
        // absorbs design-requirements rounding.
        constexpr double dutyTolerance = 0.01;
        if (dutyCycle > maximumDutyCycle * (1.0 + dutyTolerance)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Cuk::calculate_duty_cycle: duty cycle " + std::to_string(dutyCycle) +
                " exceeds maximumDutyCycle " + std::to_string(maximumDutyCycle) +
                " — converter would lose regulation; reduce |Vo|, raise Vin, or relax maximumDutyCycle.");
        }
        return dutyCycle;
    }

    double Cuk::calculate_conversion_ratio(double dutyCycle) {
        return -dutyCycle / (1.0 - dutyCycle);
    }

    double Cuk::calculate_coupling_cap_voltage(double inputVoltage, double dutyCycle) {
        return inputVoltage / (1.0 - dutyCycle);
    }

    double Cuk::calculate_l1_min(double inputVoltage, double dutyCycle, double deltaIL1, double switchingFrequency) {
        if (deltaIL1 <= 0 || switchingFrequency <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_l1_min: deltaIL1 and fsw must be > 0");
        }
        return inputVoltage * dutyCycle / (deltaIL1 * switchingFrequency);
    }

    double Cuk::calculate_l2_min(double outputVoltageMagnitude, double dutyCycle, double deltaIL2, double switchingFrequency) {
        if (deltaIL2 <= 0 || switchingFrequency <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_l2_min: deltaIL2 and fsw must be > 0");
        }
        return outputVoltageMagnitude * (1.0 - dutyCycle) / (deltaIL2 * switchingFrequency);
    }

    double Cuk::calculate_c1_min(double outputCurrent, double dutyCycle, double deltaVC1, double switchingFrequency) {
        if (deltaVC1 <= 0 || switchingFrequency <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_c1_min: deltaVC1 and fsw must be > 0");
        }
        return outputCurrent * dutyCycle / (deltaVC1 * switchingFrequency);
    }

    double Cuk::calculate_dcm_K(double L1, double L2, double switchingFrequency, double loadResistance) {
        if (L1 + L2 <= 0 || loadResistance <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_dcm_K: bad L1/L2 or R");
        }
        double Le = (L1 * L2) / (L1 + L2);
        return 2.0 * Le * switchingFrequency / loadResistance;
    }

    double Cuk::calculate_rhp_zero_frequency(double loadResistance, double dutyCycle, double L2) {
        if (L2 <= 0 || dutyCycle <= 0 || dutyCycle >= 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_rhp_zero_frequency: bad D or L2");
        }
        // ωRHP ≈ R·(1-D)² / (D·L2);  fRHP = ωRHP / (2π)
        double omegaRhp = loadResistance * std::pow(1.0 - dutyCycle, 2) / (dutyCycle * L2);
        return omegaRhp / (2.0 * std::numbers::pi);
    }

    double Cuk::compute_l1_for_isolated() const {
        // Mirrors the L1 sizing logic in process_design_requirements (V1 path):
        // worst-case L1 across all OPs at maximum Vin.  We require an explicit
        // ripple budget — no silent default — because V3 isolated Cuk has no
        // primary-magnetic L1 design requirement to fall back on.
        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument(
                "Cuk::compute_l1_for_isolated: V3 isolated requires either "
                "currentRippleRatio or maximumSwitchCurrent to size the input "
                "choke L1 (no fallback per CLAUDE.md)");
        }
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double efficiency = require_input(get_efficiency(), "Cuk", "efficiency");
        double turnsRatio = get_turns_ratio().value_or(1.0);

        double maximumDeltaIL1 = 0.0;
        if (get_current_ripple_ratio()) {
            double rippleRatio = get_current_ripple_ratio().value();
            for (const auto& op : get_operating_points()) {
                double Iout = op.get_output_currents()[0];
                double Vo   = std::abs(op.get_output_voltages()[0]);
                double D    = calculate_duty_cycle(maximumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, turnsRatio, maximumDutyCycle.value_or(0.95));
                double IL1avg = Iout * D / ((1.0 - D) * turnsRatio * efficiency);
                maximumDeltaIL1 = std::max(maximumDeltaIL1, rippleRatio * IL1avg);
            }
        }
        if (get_maximum_switch_current()) {
            double IsMax = get_maximum_switch_current().value();
            for (const auto& op : get_operating_points()) {
                double Iout = op.get_output_currents()[0];
                double Vo   = std::abs(op.get_output_voltages()[0]);
                double D    = calculate_duty_cycle(minimumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, turnsRatio, maximumDutyCycle.value_or(0.95));
                double IL1avg = Iout * D / ((1.0 - D) * turnsRatio * efficiency);
                double IL2avg = Iout;
                // Switch-current peak: IS_pk = IL1avg + IL2avg/n + ripples/2.
                double residual = IsMax - (IL1avg + IL2avg / turnsRatio) - 0.5 * 0.30 * IL2avg / turnsRatio;
                double deltaIL1 = std::max(2.0 * residual, 0.0);
                maximumDeltaIL1 = std::max(maximumDeltaIL1, deltaIL1);
            }
        }
        if (maximumDeltaIL1 <= 0) {
            throw std::invalid_argument(
                "Cuk::compute_l1_for_isolated: derived ΔIL1 budget is non-positive");
        }
        double maximumNeededInductance = 0.0;
        for (const auto& op : get_operating_points()) {
            double switchingFrequency = op.get_switching_frequency();
            double Vo = std::abs(op.get_output_voltages()[0]);
            double D  = calculate_duty_cycle(maximumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, turnsRatio, maximumDutyCycle.value_or(0.95));
            double L1 = calculate_l1_min(maximumInputVoltage, D, maximumDeltaIL1, switchingFrequency);
            maximumNeededInductance = std::max(maximumNeededInductance, L1);
        }
        return maximumNeededInductance;
    }

    // ============================================================
    //   Constructors
    // ============================================================

    Cuk::Cuk(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        MAS::from_json(migrated, *static_cast<MAS::Cuk*>(this));
    }

    AdvancedCuk::AdvancedCuk(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    // ============================================================
    //   Per-operating-point analytical worker
    // ============================================================

    OperatingPoint Cuk::process_operating_points_for_input_voltage(double inputVoltage, const MAS::CukOperatingPoint& outputOperatingPoint, double inductanceL1) {
        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double outputVoltageMag   = std::abs(outputOperatingPoint.get_output_voltages()[0]);
        double outputCurrent      = outputOperatingPoint.get_output_currents()[0];
        double diodeVoltageDrop   = get_diode_voltage_drop();
        double efficiency = require_input(get_efficiency(), "Cuk", "efficiency");

        const bool isIsolated = get_isolated().value_or(false);
        const double turnsRatio = isIsolated ? get_turns_ratio().value_or(1.0) : 1.0;
        if (isIsolated && turnsRatio <= 0) {
            throw std::invalid_argument(
                "Cuk::process_operating_points_for_input_voltage: turnsRatio must "
                "be > 0 for isolated Cuk; received " + std::to_string(turnsRatio));
        }
        lastTurnsRatio = turnsRatio;

        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltageMag, diodeVoltageDrop, efficiency, turnsRatio, maximumDutyCycle.value_or(0.95));

        // Common derived quantities (V1/V2/V3).
        // Power balance: Vin·η·IL1avg = |Vo|·Iout  ⇒  IL1avg = |Vo|·Iout/(Vin·η).
        // For V1 (n=1, η=1) this simplifies to Iout·D/(1-D); the V3 form uses
        // |Vo|/Vin = D/((1-D)·n) so IL1avg = Iout·D/((1-D)·n·η).
        double IL2avg     = outputCurrent;
        double IL1avg     = outputCurrent * dutyCycle / ((1.0 - dutyCycle) * turnsRatio * efficiency);
        double VCa        = inputVoltage / (1.0 - dutyCycle);            // primary-side coupling-cap DC (V3) / VC1 (V1)
        double VCb        = isIsolated ? outputVoltageMag / dutyCycle    // secondary-side coupling-cap DC (V3 only)
                                       : 0.0;
        // For V3 the inductanceL1 argument is repurposed as Lm (transformer
        // magnetizing inductance) — process_design_requirements set it that
        // way. The actual input-choke L1 is sized separately.
        double L1Choke    = isIsolated ? compute_l1_for_isolated() : inductanceL1;
        double deltaIL1   = L1Choke > 0 ? (inputVoltage * dutyCycle) / (L1Choke * switchingFrequency)
                                        : 0.0;
        double deltaIL2_target = std::max(l2RipplePct * IL2avg, 1e-6);
        double inductanceL2    = calculate_l2_min(outputVoltageMag, dutyCycle, deltaIL2_target, switchingFrequency);
        double deltaIL2        = (outputVoltageMag * (1.0 - dutyCycle)) / (inductanceL2 * switchingFrequency);

        // ---- Common diagnostics ----
        lastDutyCycle = dutyCycle;
        lastConversionRatio = -dutyCycle / ((1.0 - dutyCycle) * turnsRatio);
        lastInputInductorAverage = IL1avg;
        lastOutputInductorAverage = IL2avg;
        lastInputInductorRipple = deltaIL1;
        lastOutputInductorRipple = deltaIL2;
        // Switch (S1) blocks VCa during OFF (Cuk classical result; isolation does
        // not change switch stress on primary). Diode peak reverse:
        //   V1/V2: VS,peak = VD,peak = Vin + |Vo| = VCa.
        //   V3:    VS,peak = VCa = Vin/(1-D); VD,peak = VCb = |Vo|/D.
        lastSwitchPeakVoltage = isIsolated ? VCa : (inputVoltage + outputVoltageMag);
        lastDiodePeakReverseVoltage = isIsolated ? VCb : lastSwitchPeakVoltage;
        double IS_peak = IL1avg + IL2avg / std::max(turnsRatio, 1e-12) + (deltaIL1 + deltaIL2 / std::max(turnsRatio, 1e-12)) / 2.0;
        lastSwitchPeakCurrent = IS_peak;
        lastDiodePeakCurrent  = IL2avg + deltaIL2 / 2.0;
        double loadResistance = (outputCurrent > 0) ? outputVoltageMag / outputCurrent : 0.0;
        lastDcmK = calculate_dcm_K(L1Choke, inductanceL2, switchingFrequency, loadResistance);
        lastDcmKcrit = std::pow(1.0 - dutyCycle, 2);
        lastIsCcm = (lastDcmK > lastDcmKcrit);
        lastRhpZeroFrequency = calculate_rhp_zero_frequency(loadResistance, dutyCycle, inductanceL2);
        lastRecommendedLoopBandwidth = lastRhpZeroFrequency / 5.0;
        lastSizedL2 = inductanceL2;

        // Output cap Co sized from ΔIL2 / (8·fsw·ΔVo) (LC filter rule).
        double deltaVo_target = std::max(coRipplePct * outputVoltageMag, 1e-3);
        double outputCapacitance = deltaIL2 / (8.0 * switchingFrequency * deltaVo_target);
        outputCapacitance = std::max(outputCapacitance, 1e-6);
        lastSizedCo = outputCapacitance;
        double deltaVo_actual = deltaIL2 / (8.0 * switchingFrequency * outputCapacitance);

        if (isIsolated) {
            // ============================================================
            //   V3 isolated: primary excitation = transformer Lp winding
            // ============================================================
            //
            // Lm (= inductanceL1 argument repurposed in V3) is the magnetizing
            // inductance returned by process_design_requirements above for the
            // isolated path. Magnetizing current ripple (Erickson §6):
            //     ΔIm = VCa · D · T / Lm   with Im,avg = 0 (symmetric).
            // V·s balance on Lp: VLp_ON · D + VLp_OFF · (1-D) = 0
            //                  ⇒ VLp_OFF = -VCa · D/(1-D).
            // Peak-to-peak winding voltage = VCa - VLp_OFF = VCa / (1-D)
            //                              = Vin / (1-D)².
            const double Lm = inductanceL1;   // sized by process_design_requirements
            const double VCa_DC = VCa;
            double deltaIm = (Lm > 0) ? (VCa_DC * dutyCycle) / (Lm * switchingFrequency) : 0.0;
            double primaryVppLp = VCa_DC / (1.0 - dutyCycle);

            lastSizedLm   = Lm;
            lastCouplingCapVoltage = VCa_DC;     // VCa, primary side
            // For RMS estimate use the same √(D(1-D))(IL1avg+IL2avg/n) form,
            // valid for the primary-side coupling cap which carries the same AC.
            lastCouplingCapRmsCurrent = std::sqrt(dutyCycle * (1.0 - dutyCycle)) *
                                        (IL1avg + IL2avg / turnsRatio);

            // Size Ca and Cb so their ripple stays within c1RipplePct of the DC level.
            double deltaVCa_target = std::max(c1RipplePct * VCa_DC, 1e-3);
            double couplingCa = (IL2avg / turnsRatio) * (1.0 - dutyCycle) / (deltaVCa_target * switchingFrequency);
            double deltaVCb_target = std::max(c1RipplePct * VCb, 1e-3);
            double couplingCb = IL2avg * dutyCycle / (deltaVCb_target * switchingFrequency);
            lastSizedCa = couplingCa;
            lastSizedCb = couplingCb;
            // V1-shape diagnostics that have no V3 counterpart are zeroed for clarity.
            lastSizedC1 = 0.0;

            // ---- Primary excitation = Lp (magnetizing) ----
            {
                Waveform vLp = Inputs::create_waveform(
                    WaveformLabel::RECTANGULAR, primaryVppLp, switchingFrequency, dutyCycle, 0.0, 0);
                Waveform iLp = Inputs::create_waveform(
                    WaveformLabel::TRIANGULAR, deltaIm, switchingFrequency, dutyCycle, 0.0, 0);
                auto excitation = complete_excitation(iLp, vLp, switchingFrequency, "Primary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }
            // ---- Secondary excitation = Ls (turns_ratio scaled) ----
            {
                double primaryVppLs = primaryVppLp / turnsRatio;
                double deltaIs      = deltaIm * turnsRatio;       // ampere-turn balance
                Waveform vLs = Inputs::create_waveform(
                    WaveformLabel::RECTANGULAR, primaryVppLs, switchingFrequency, dutyCycle, 0.0, 0);
                Waveform iLs = Inputs::create_waveform(
                    WaveformLabel::TRIANGULAR, deltaIs, switchingFrequency, dutyCycle, 0.0, 0);
                auto excitation = complete_excitation(iLs, vLs, switchingFrequency, "Secondary");
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }

            // ---- Stash extras: L1 (input choke), L2, Ca, Cb, Co ----
            {
                // L1 voltage: +Vin during ON, Vin-VCa during OFF (= V1's L1).
                Waveform vL1 = Inputs::create_waveform(
                    WaveformLabel::RECTANGULAR, VCa_DC, switchingFrequency, dutyCycle, 0.0, 0);
                Waveform iL1 = Inputs::create_waveform(
                    WaveformLabel::TRIANGULAR, deltaIL1, switchingFrequency, dutyCycle, IL1avg, 0);
                extraL1VoltageWaveforms.push_back(vL1);
                extraL1CurrentWaveforms.push_back(iL1);
            }
            {
                Waveform vL2 = Inputs::create_waveform(
                    WaveformLabel::RECTANGULAR, VCb, switchingFrequency, 1.0 - dutyCycle, 0.0, 0);
                Waveform iL2 = Inputs::create_waveform(
                    WaveformLabel::TRIANGULAR, deltaIL2, switchingFrequency, 1.0 - dutyCycle, IL2avg, 0);
                extraL2VoltageWaveforms.push_back(vL2);
                extraL2CurrentWaveforms.push_back(iL2);
            }
            {
                double deltaVCa_actual = (couplingCa > 0)
                    ? (IL2avg / turnsRatio) * (1.0 - dutyCycle) / (couplingCa * switchingFrequency)
                    : 0.0;
                Waveform vCa = Inputs::create_waveform(
                    WaveformLabel::TRIANGULAR, deltaVCa_actual, switchingFrequency, dutyCycle, VCa_DC, 0);
                // Ca AC current = transformer primary AC current on the Ca side
                // (rectangular ±IL1avg over ON / -(IL2avg/n) over OFF). Use the
                // same pp shape as V1's C1, scaled appropriately.
                double iCapp = IL1avg + IL2avg / turnsRatio;
                Waveform iCa = Inputs::create_waveform(
                    WaveformLabel::RECTANGULAR, iCapp, switchingFrequency, 1.0 - dutyCycle, 0.0, 0);
                extraCaVoltageWaveforms.push_back(vCa);
                extraCaCurrentWaveforms.push_back(iCa);
            }
            {
                double deltaVCb_actual = (couplingCb > 0)
                    ? IL2avg * dutyCycle / (couplingCb * switchingFrequency)
                    : 0.0;
                Waveform vCb = Inputs::create_waveform(
                    WaveformLabel::TRIANGULAR, deltaVCb_actual, switchingFrequency, dutyCycle, VCb, 0);
                double iCbpp = IL1avg * turnsRatio + IL2avg;     // reflected to secondary
                Waveform iCb = Inputs::create_waveform(
                    WaveformLabel::RECTANGULAR, iCbpp, switchingFrequency, 1.0 - dutyCycle, 0.0, 0);
                extraCbVoltageWaveforms.push_back(vCb);
                extraCbCurrentWaveforms.push_back(iCb);
            }
            {
                Waveform vCo = Inputs::create_waveform(
                    WaveformLabel::TRIANGULAR, deltaVo_actual, switchingFrequency, 1.0 - dutyCycle, -outputVoltageMag, 0);
                Waveform iCo = Inputs::create_waveform(
                    WaveformLabel::TRIANGULAR, deltaIL2, switchingFrequency, 1.0 - dutyCycle, 0.0, 0);
                extraCoVoltageWaveforms.push_back(vCo);
                extraCoCurrentWaveforms.push_back(iCo);
            }

            OperatingConditions conditions;
            conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
            conditions.set_cooling(std::nullopt);
            operatingPoint.set_conditions(conditions);
            // Per-OP diagnostic snapshot.
            perOpDutyCycle.push_back(lastDutyCycle);
            perOpConversionRatio.push_back(lastConversionRatio);
            perOpCouplingCapVoltage.push_back(lastCouplingCapVoltage);
            perOpInputInductorAverage.push_back(lastInputInductorAverage);
            perOpOutputInductorAverage.push_back(lastOutputInductorAverage);
            perOpInputInductorRipple.push_back(lastInputInductorRipple);
            perOpOutputInductorRipple.push_back(lastOutputInductorRipple);
            perOpSwitchPeakVoltage.push_back(lastSwitchPeakVoltage);
            perOpSwitchPeakCurrent.push_back(lastSwitchPeakCurrent);
            perOpDiodePeakReverseVoltage.push_back(lastDiodePeakReverseVoltage);
            perOpDiodePeakCurrent.push_back(lastDiodePeakCurrent);
            perOpCouplingCapRmsCurrent.push_back(lastCouplingCapRmsCurrent);
            perOpIsCcm.push_back(lastIsCcm);
            perOpSizedCa.push_back(lastSizedCa);
            perOpSizedCb.push_back(lastSizedCb);
            perOpSizedCo.push_back(lastSizedCo);
            perOpRhpZeroFrequency.push_back(lastRhpZeroFrequency);

            return operatingPoint;
        }

        // ============================================================
        //   V1 / V2 (non-isolated) path — preserved verbatim
        // ============================================================
        const double VC1 = VCa;   // alias: in V1/V2 the single coupling cap is C1 = Ca
        double deltaVC1_target = std::max(c1RipplePct * VC1, 1e-3);
        double couplingCap    = calculate_c1_min(outputCurrent, dutyCycle, deltaVC1_target, switchingFrequency);

        // Primary inductor (L1) waveforms — TRIANGULAR around IL1avg.
        // Voltage swings between +Vin (ON) and Vin - VC1 (OFF, negative).
        double primaryVoltageMaximum = inputVoltage;             // ON
        double primaryVoltageMinimum = inputVoltage - VC1;       // OFF (negative)
        double primaryVoltagePeakToPeak = primaryVoltageMaximum - primaryVoltageMinimum;  // = VC1
        double primaryCurrentMinimum = IL1avg - deltaIL1 / 2.0;

        lastCouplingCapVoltage = VC1;
        lastCouplingCapRmsCurrent = std::sqrt(dutyCycle * (1.0 - dutyCycle)) * (IL1avg + IL2avg);
        lastSizedC1 = couplingCap;

        // Primary winding excitation (L1)
        {
            Waveform currentWaveform = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaIL1, switchingFrequency, dutyCycle, IL1avg, 0);
            Waveform voltageWaveform = Inputs::create_waveform(
                WaveformLabel::RECTANGULAR, primaryVoltagePeakToPeak, switchingFrequency, dutyCycle, 0, 0);

            (void) primaryCurrentMinimum;  // emitted via diagnostics; not asserted here

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // ---- Stash extra-component waveforms for get_extra_components_inputs ----
        // L2 voltage: ON value = -VC1+|Vo| (negative), OFF value = +|Vo|.
        // create_waveform(RECTANGULAR, pp, f, dc, off) emits a waveform that
        // is "high" during dc·T and "low" during (1-dc)·T; in Cuk L2 is high
        // (positive) during the OFF sub-interval, so pass dutyCycle=(1-D).
        // peakToPeak = max-min = |Vo| - (-VC1+|Vo|) = VC1.  Mean = 0 by V·s
        // balance, so offset=0.
        {
            Waveform vL2 = Inputs::create_waveform(
                WaveformLabel::RECTANGULAR, VC1, switchingFrequency, 1.0 - dutyCycle, 0.0, 0);
            Waveform iL2 = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaIL2, switchingFrequency, 1.0 - dutyCycle, IL2avg, 0);
            extraL2VoltageWaveforms.push_back(vL2);
            extraL2CurrentWaveforms.push_back(iL2);
        }

        // C1 voltage: triangular around VC1 with peak-to-peak ΔVC1.
        // C1 current: rectangular between -IL2avg (during ON, D·T) and
        // +IL1avg (during OFF, (1-D)·T).  pp = IL1avg + IL2avg = Iout/(1-D).
        {
            double deltaVC1_actual = (couplingCap > 0)
                ? (outputCurrent * dutyCycle) / (couplingCap * switchingFrequency)
                : 0.0;
            Waveform vC1 = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaVC1_actual, switchingFrequency, dutyCycle, VC1, 0);
            double iC1pp = IL1avg + IL2avg;
            Waveform iC1 = Inputs::create_waveform(
                WaveformLabel::RECTANGULAR, iC1pp, switchingFrequency, 1.0 - dutyCycle, 0.0, 0);
            extraC1VoltageWaveforms.push_back(vC1);
            extraC1CurrentWaveforms.push_back(iC1);
        }

        // Co voltage: tiny triangular ripple ΔVo around the signed mean -|Vo|.
        // Co current: AC component of iL2 (its triangular ripple) around 0 mean.
        {
            Waveform vCo = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaVo_actual, switchingFrequency, 1.0 - dutyCycle, -outputVoltageMag, 0);
            Waveform iCo = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaIL2, switchingFrequency, 1.0 - dutyCycle, 0.0, 0);
            extraCoVoltageWaveforms.push_back(vCo);
            extraCoCurrentWaveforms.push_back(iCo);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        // Per-OP diagnostic snapshot.
        perOpDutyCycle.push_back(lastDutyCycle);
        perOpConversionRatio.push_back(lastConversionRatio);
        perOpCouplingCapVoltage.push_back(lastCouplingCapVoltage);
        perOpInputInductorAverage.push_back(lastInputInductorAverage);
        perOpOutputInductorAverage.push_back(lastOutputInductorAverage);
        perOpInputInductorRipple.push_back(lastInputInductorRipple);
        perOpOutputInductorRipple.push_back(lastOutputInductorRipple);
        perOpSwitchPeakVoltage.push_back(lastSwitchPeakVoltage);
        perOpSwitchPeakCurrent.push_back(lastSwitchPeakCurrent);
        perOpDiodePeakReverseVoltage.push_back(lastDiodePeakReverseVoltage);
        perOpDiodePeakCurrent.push_back(lastDiodePeakCurrent);
        perOpCouplingCapRmsCurrent.push_back(lastCouplingCapRmsCurrent);
        perOpIsCcm.push_back(lastIsCcm);
        perOpSizedCa.push_back(lastSizedCa);
        perOpSizedCb.push_back(lastSizedCb);
        perOpSizedCo.push_back(lastSizedCo);
        perOpRhpZeroFrequency.push_back(lastRhpZeroFrequency);

        return operatingPoint;
    }

    bool Cuk::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) return false;
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        if (!get_input_voltage().get_nominal() && !get_input_voltage().get_maximum() && !get_input_voltage().get_minimum()) {
            if (!assert) return false;
            throw InvalidInputException(ErrorCode::MISSING_DATA, "No input voltage introduced");
        }
        return true;
    }

    DesignRequirements Cuk::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double efficiency = require_input(get_efficiency(), "Cuk", "efficiency");

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Cuk::process_design_requirements: missing both currentRippleRatio and maximumSwitchCurrent");
        }

        // Worst-case L1: largest L1 needed across all OPs at maximum Vin
        // (largest Vin → largest D·Vin/fsw product).
        double maximumDeltaIL1 = 0.0;
        if (get_current_ripple_ratio()) {
            double rippleRatio = get_current_ripple_ratio().value();
            for (const auto& op : get_operating_points()) {
                double Iout = op.get_output_currents()[0];
                double Vo   = std::abs(op.get_output_voltages()[0]);
                double D    = calculate_duty_cycle(maximumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, 1.0, maximumDutyCycle.value_or(0.95));
                double IL1avg = Iout * D / (1.0 - D);
                maximumDeltaIL1 = std::max(maximumDeltaIL1, rippleRatio * IL1avg);
            }
        }
        if (get_maximum_switch_current()) {
            double IsMax = get_maximum_switch_current().value();
            for (const auto& op : get_operating_points()) {
                double Iout = op.get_output_currents()[0];
                double Vo   = std::abs(op.get_output_voltages()[0]);
                double D    = calculate_duty_cycle(minimumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, 1.0, maximumDutyCycle.value_or(0.95));
                double IL1avg = Iout * D / (1.0 - D);
                // IS_peak ≈ IL1avg + IL2avg + (ΔIL1 + ΔIL2)/2; treat ΔIL2 as
                // ~30 % of IL2avg (the L2 sizing default), so the L1 ripple
                // budget that just hits IS_peak = IsMax is:
                double IL2avg = Iout;
                double residual = IsMax - (IL1avg + IL2avg) - 0.5 * 0.30 * IL2avg;
                double deltaIL1 = std::max(2.0 * residual, 0.0);
                maximumDeltaIL1 = std::max(maximumDeltaIL1, deltaIL1);
            }
        }
        if (maximumDeltaIL1 <= 0) {
            throw std::invalid_argument("Cuk::process_design_requirements: derived ΔIL1 budget is non-positive");
        }

        double maximumNeededInductance = 0.0;
        for (const auto& op : get_operating_points()) {
            double switchingFrequency = op.get_switching_frequency();
            double Vo = std::abs(op.get_output_voltages()[0]);
            double D  = calculate_duty_cycle(maximumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, 1.0, maximumDutyCycle.value_or(0.95));
            double L1 = calculate_l1_min(maximumInputVoltage, D, maximumDeltaIL1, switchingFrequency);
            maximumNeededInductance = std::max(maximumNeededInductance, L1);
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        const bool isCoupled  = get_coupled_inductor().value_or(false);
        const bool isIsolated = get_isolated().value_or(false);
        if (isCoupled && isIsolated) {
            throw std::invalid_argument(
                "Cuk::process_design_requirements: coupledInductor and isolated "
                "are mutually exclusive (V2 uses a single coupled inductor on "
                "one core; V3 uses a separate isolation transformer).");
        }

        std::vector<IsolationSide> isolationSides;
        if (isIsolated) {
            // V3 isolated: the "primary" magnetic returned here is the
            // 2-winding isolation transformer (Lp:Ls). The DC-blocking input
            // inductor L1 and output inductor L2 are still present but are
            // emitted as ancillaries via get_extra_components_inputs.
            //
            // The magnetizing inductance Lm referenced from process_design_requirements
            // is sized so that the magnetizing-current ripple ΔIm ≤ 0.20·IL1avg
            // across all OPs at minimum Vin (worst case for ΔIm/IL1avg ratio).
            // ΔIm = VCa · D · T / Lm, with VCa = Vin / (1-D).
            const double turnsRatio = get_turns_ratio().value_or(1.0);
            if (turnsRatio <= 0) {
                throw std::invalid_argument(
                    "Cuk::process_design_requirements: turnsRatio must be > 0 "
                    "for isolated Cuk; received " + std::to_string(turnsRatio));
            }
            double maximumLm = 0.0;
            for (const auto& op : get_operating_points()) {
                double switchingFrequency = op.get_switching_frequency();
                double Vo   = std::abs(op.get_output_voltages()[0]);
                double Iout = op.get_output_currents()[0];
                double D    = calculate_duty_cycle(minimumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, turnsRatio, maximumDutyCycle.value_or(0.95));
                double IL1avg = Iout * D / ((1.0 - D) * turnsRatio * efficiency);
                double VCa    = minimumInputVoltage / (1.0 - D);
                double deltaIm_target = std::max(0.20 * IL1avg, 1e-6);
                double Lm     = VCa * D / (deltaIm_target * switchingFrequency);
                maximumLm = std::max(maximumLm, Lm);
            }
            // Override the L1-based magnetizing-inductance requirement: in V3
            // the "primary magnetic" is the transformer, sized by Lm.
            DimensionWithTolerance lmWithTolerance;
            lmWithTolerance.set_minimum(roundFloat(maximumLm, 10));
            designRequirements.set_magnetizing_inductance(lmWithTolerance);

            // 2-winding transformer: turns_ratios = {Np/Ns}, isolation_sides = {PRIMARY, SECONDARY}.
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 4));
            designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{turnsRatioWithTolerance});
            isolationSides.push_back(get_isolation_side_from_index(0));
            isolationSides.push_back(get_isolation_side_from_index(1));
        } else if (isCoupled) {
            // V2 coupled-inductor: two windings on the same core, no galvanic
            // isolation. Erickson zero-ripple condition: k * sqrt(L2/L1) = 1.
            // We size both windings identically (N1:N2 = 1:1, L1 = L2) so the
            // canonical k -> 1 case yields exact ripple cancellation. The user
            // can perturb the coupling coefficient via couplingCoefficient.
            isolationSides.push_back(get_isolation_side_from_index(0));
            isolationSides.push_back(get_isolation_side_from_index(0));
            DimensionWithTolerance unityRatio;
            unityRatio.set_nominal(1.0);
            designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{unityRatio});
        } else {
            isolationSides.push_back(get_isolation_side_from_index(0));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(MAS::Topology::CUK_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Cuk::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        (void) turnsRatios;  // Cuk V1 has no transformer; turns ratios are unused.

        // Reset extra-component waveform stash so a re-run does not accumulate.
        extraL1VoltageWaveforms.clear();
        extraL1CurrentWaveforms.clear();
        extraL2VoltageWaveforms.clear();
        extraL2CurrentWaveforms.clear();
        extraC1VoltageWaveforms.clear();
        extraC1CurrentWaveforms.clear();
        extraCaVoltageWaveforms.clear();
        extraCaCurrentWaveforms.clear();
        extraCbVoltageWaveforms.clear();
        extraCbCurrentWaveforms.clear();
        extraCoVoltageWaveforms.clear();
        extraCoCurrentWaveforms.clear();

        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        // Clear per-OP diagnostic vectors so the wizard table reflects this run only.
        perOpName.clear();
        perOpDutyCycle.clear();
        perOpConversionRatio.clear();
        perOpCouplingCapVoltage.clear();
        perOpInputInductorAverage.clear();
        perOpOutputInductorAverage.clear();
        perOpInputInductorRipple.clear();
        perOpOutputInductorRipple.clear();
        perOpSwitchPeakVoltage.clear();
        perOpSwitchPeakCurrent.clear();
        perOpDiodePeakReverseVoltage.clear();
        perOpDiodePeakCurrent.clear();
        perOpCouplingCapRmsCurrent.clear();
        perOpIsCcm.clear();
        perOpSizedCa.clear();
        perOpSizedCb.clear();
        perOpSizedCo.clear();
        perOpRhpZeroFrequency.clear();


        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                std::string opName = inputVoltagesNames[inputVoltageIndex];
                if (get_operating_points().size() > 1) {
                    opName += " · OP" + std::to_string(opIndex);
                }
                perOpName.push_back(opName);
                auto operatingPoint = process_operating_points_for_input_voltage(
                    inputVoltage, get_operating_points()[opIndex], magnetizingInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(opIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> Cuk::process_operating_points(OpenMagnetics::Magnetic magnetic) {
        run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(
            magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();

        std::vector<double> turnsRatios = magnetic.get_turns_ratios();

        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedCuk::process() {
        Cuk::run_checks(_assertErrors);

        // Reset extra-component waveform stash so a re-run does not accumulate.
        extraL1VoltageWaveforms.clear();
        extraL1CurrentWaveforms.clear();
        extraL2VoltageWaveforms.clear();
        extraL2CurrentWaveforms.clear();
        extraC1VoltageWaveforms.clear();
        extraC1CurrentWaveforms.clear();
        extraCaVoltageWaveforms.clear();
        extraCaCurrentWaveforms.clear();
        extraCbVoltageWaveforms.clear();
        extraCbCurrentWaveforms.clear();
        extraCoVoltageWaveforms.clear();
        extraCoCurrentWaveforms.clear();

        Inputs inputs;
        double maximumNeededInductance = get_desired_inductance();

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        DesignRequirements designRequirements;
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(MAS::Topology::CUK_CONVERTER);
        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(
                    inputVoltage, get_operating_points()[opIndex], maximumNeededInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(opIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }
        return inputs;
    }

    // ============================================================
    //   SPICE netlist (V1 non-isolated)
    // ============================================================
    //
    // CRITICAL per CUK_PLAN.md §6:
    //   * IC pre-charge for L1, L2, C1, Co (4 reactive elements) —
    //     without these the 4th-order LC tank diverges at cold-start.
    //   * METHOD=GEAR — trapezoidal mishandles the hard-switching edges
    //     on a 4th-order tank.
    //   * Small parasitic R baked into L1, L2 (DCR ≈ 50 mΩ) and ESR
    //     into C1, Co (≈ 5 mΩ). Pure reactive networks make ngspice
    //     unhappy.
    //   * Snubbers across S and D.
    //   * Bridge-voltage probe (sw_probe = V(node_A)).
    // ============================================================

    std::string Cuk::generate_ngspice_circuit(double inductanceL1, size_t inputVoltageIndex, size_t operatingPointIndex) {
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames_;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames_);

        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("Cuk::generate_ngspice_circuit: inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("Cuk::generate_ngspice_circuit: operatingPointIndex out of range");
        }

        double inputVoltage = inputVoltages[inputVoltageIndex];
        auto opPoint = get_operating_points()[operatingPointIndex];

        double outputVoltageMag   = std::abs(opPoint.get_output_voltages()[0]);
        double outputCurrent      = opPoint.get_output_currents()[0];
        double switchingFrequency = opPoint.get_switching_frequency();
        double diodeVoltageDrop   = get_diode_voltage_drop();
        double efficiency = require_input(get_efficiency(), "Cuk", "efficiency");

        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltageMag, diodeVoltageDrop, efficiency, 1.0, maximumDutyCycle.value_or(0.95));

        // Internally-sized L2, C1, Co (mirror process_operating_points_for_input_voltage)
        double IL1avg = outputCurrent * dutyCycle / (1.0 - dutyCycle);
        double IL2avg = outputCurrent;
        double VC1    = calculate_coupling_cap_voltage(inputVoltage, dutyCycle);
        double deltaIL2_target = std::max(l2RipplePct * IL2avg, 1e-6);
        double inductanceL2    = calculate_l2_min(outputVoltageMag, dutyCycle, deltaIL2_target, switchingFrequency);
        double deltaVC1_target = std::max(c1RipplePct * VC1, 1e-3);
        double couplingCap     = calculate_c1_min(outputCurrent, dutyCycle, deltaVC1_target, switchingFrequency);
        // Output cap: ΔVo / |Vo| ≤ coRipplePct.  Use the standard
        // continuous-current LC ripple formula  Co = ΔIL2 / (8·fsw·ΔVo).
        double deltaIL2 = (outputVoltageMag * (1.0 - dutyCycle)) / (inductanceL2 * switchingFrequency);
        double deltaVo  = std::max(coRipplePct * outputVoltageMag, 1e-3);
        double outputCapacitance = deltaIL2 / (8.0 * switchingFrequency * deltaVo);
        outputCapacitance = std::max(outputCapacitance, 1e-6);  // sanity floor for SPICE conditioning

        if (outputCurrent <= 0.0) {
            throw std::invalid_argument(
                "Cuk::generate_ngspice_circuit: outputCurrent must be > 0 "
                "(received " + std::to_string(outputCurrent) +
                "); cannot derive load resistance R_load = |Vo| / Iout.");
        }
        double loadResistance = outputVoltageMag / outputCurrent;

        const SpiceSimulationConfig cfg = spice_config();

        // V4/V5 variant flags. bidirectional implies synchronous regardless of the
        // explicit value (per schema: "Implies synchronous=true").
        const bool isBidirectional = get_bidirectional().value_or(false);
        const bool isSynchronous   = isBidirectional || get_synchronous().value_or(false);
        const bool isCoupled       = get_coupled_inductor().value_or(false);
        const bool isIsolated      = get_isolated().value_or(false);
        if (isCoupled && isIsolated) {
            throw std::invalid_argument(
                "Cuk::generate_ngspice_circuit: coupledInductor and isolated "
                "are mutually exclusive (V2 vs V3)");
        }
        const double turnsRatio    = isIsolated ? get_turns_ratio().value_or(1.0) : 1.0;
        if (isIsolated && turnsRatio <= 0) {
            throw std::invalid_argument(
                "Cuk::generate_ngspice_circuit: turnsRatio must be > 0 for "
                "isolated Cuk; received " + std::to_string(turnsRatio));
        }
        // For V3 the inductanceL1 argument is the transformer Lm; the input
        // choke L1 is sized internally.  For V1/V2 inductanceL1 is L1 directly.
        const double L1Choke   = isIsolated ? compute_l1_for_isolated() : inductanceL1;
        const double Lm        = isIsolated ? inductanceL1 : 0.0;
        // V3 secondary-side cap DC voltage and per-cap sizing (mirror analytical).
        const double VCb_DC    = isIsolated ? outputVoltageMag / dutyCycle : 0.0;
        const double IL1avgV3  = isIsolated
            ? outputCurrent * dutyCycle / ((1.0 - dutyCycle) * turnsRatio * efficiency)
            : IL1avg;
        const double couplingCa = isIsolated
            ? (outputCurrent / turnsRatio) * (1.0 - dutyCycle) /
              (std::max(c1RipplePct * VC1, 1e-3) * switchingFrequency)
            : 0.0;
        const double couplingCb = isIsolated
            ? outputCurrent * dutyCycle /
              (std::max(c1RipplePct * VCb_DC, 1e-3) * switchingFrequency)
            : 0.0;

        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;

        int periodsToExtract = get_num_periods_to_extract();
        int settlingPeriods  = get_num_steady_state_periods();
        const int numPeriodsTotal = settlingPeriods + periodsToExtract;
        double simTime  = numPeriodsTotal * period;
        double startTime = settlingPeriods * period;
        // Tighter step at high Fs (LM2611 1.4 MHz fixture): plan §6.2.
        double samplesPerPeriod = (switchingFrequency >= 1e6) ? 500.0 : cfg.samplesPerPeriod;
        double stepTime = period / samplesPerPeriod;

        // Parasitic series R for SPICE conditioning (CUK_PLAN.md §6.2).
        // Ideal reference: zero winding DCR + cap ESR (consistent with the other ideal decks).
        const double dcrL1 = 0.0;
        const double dcrL2 = 0.0;
        const double esrC1 = 0.0;
        const double esrCo = 0.0;

        // Cuk node convention (CUK_PLAN.md §6.1):
        //   vin_dc → Vin_sense → l1_in →[L1]→ node_A
        //   node_A →[Vc1_sense]→ node_C →[C1]→ node_B
        //   node_B →[L2]→ out_sense →[Vout_sense]→ vout_load_node
        //   S1: node_A → 0 (PWM)
        //   D1: node_B → 0 (anode at node_B, cathode at GND).
        //       During S1 ON, V(node_B) ≈ -VC1 (large negative) so D1 is
        //       reverse-biased (anode << cathode). During S1 OFF, D1
        //       conducts to clamp V(node_B) ≈ -Vd, freewheeling L2.
        //   Vo polarity: vout_load_node is NEGATIVE; load tied between
        //               vout_load_node and 0 with positive current
        //               flowing OUT of vout_load_node into the load
        //               (into 0). i(Vout_sense) reports this current.

        circuit << "* Cuk Converter ("
                << (isIsolated ? "isolated V3" : "non-isolated V1")
                << (isSynchronous ? "/V4 synchronous" : "")
                << (isBidirectional ? "/V5 bidirectional" : "")
                << (isCoupled ? "/V2 coupled-inductor" : "")
                << ") - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, |Vo|=" << outputVoltageMag << "V (sign: negative), "
                << "f=" << (switchingFrequency / 1e3) << "kHz, D=" << (dutyCycle * 100) << " pct";
        if (isIsolated) {
            circuit << ", n=Np/Ns=" << turnsRatio;
        }
        circuit << "\n";
        circuit << "* L1=" << (L1Choke * 1e6) << "uH, L2=" << (inductanceL2 * 1e6) << "uH, ";
        if (isIsolated) {
            circuit << "Lm=" << (Lm * 1e6) << "uH, Ca=" << (couplingCa * 1e6)
                    << "uF, Cb=" << (couplingCb * 1e6) << "uF, ";
        } else {
            circuit << "C1=" << (couplingCap * 1e6) << "uF, ";
        }
        circuit << "Co=" << (outputCapacitance * 1e6) << "uF, "
                << "Iout=" << outputCurrent << "A\n\n";

        // DC Input + sense
        circuit << "* DC Input + sense\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n";
        circuit << "Vin_sense vin_dc l1_in 0\n\n";

        // L1 with DCR (split as L + R for bake-in; use named series RDCR_l1)
        circuit << "* L1 input inductor (with DCR)\n";
        circuit << "L1 l1_in l1_dcr_mid " << std::scientific << L1Choke << std::fixed
                << " ic=" << IL1avgV3 << "\n";
        circuit << "Rdcr_l1 l1_dcr_mid node_A " << dcrL1 << "\n";
        circuit << "Vl1_sense node_A node_A_int 0\n\n";

        // S1 PWM switch (from node_A to GND)
        circuit << "* PWM low-side switch on node_A\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                << tOn << " " << period << std::fixed << ")\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
        circuit << "S1 node_A_int 0 pwm_ctrl 0 SW1\n";
        circuit << "Rsnub_s1 node_A_int 0 " << cfg.snubR << "\n"
                << "Csnub_s1 node_A_int snub_s1_int " << std::scientific << cfg.snubC << std::fixed << "\n"
                << "Rsnub_s1b snub_s1_int 0 " << require_spice_field(cfg.snubDampR, "snubDampR") << "\n\n";

        // Coupling section. V1/V2: single C1 between node_A_int and node_B.
        // V3 isolated: split into Ca → Lp ⟂ Ls → Cb with the transformer
        // providing galvanic isolation. The dot orientation on Lp is on the
        // tx_pri_pos side (so positive winding voltage corresponds to current
        // flowing into tx_pri_pos); on Ls the dot is on tx_sec_pos. With this
        // convention Vsec = -n·Vpri / 1 (negative for our sign convention with
        // n=Np/Ns, since the secondary is wound oppositely on the same core
        // and we want |Vo| < 0 on the load).
        if (isIsolated) {
            circuit << "* V3 isolated: Ca + transformer (Lp/Ls, K=0.999) + Cb\n";
            circuit << "Vc1_sense node_A_int node_C 0\n";
            circuit << "Ca node_C node_C_esr " << std::scientific << couplingCa << std::fixed
                    << " ic=" << VC1 << "\n";
            circuit << "Rca_esr node_C_esr tx_pri_pos " << esrC1 << "\n";
            // Transformer windings. Both reference 0; coupling K1 between Lp and Ls.
            // Lp inductance = Lm; Ls inductance = Lm / n² (turns-square scaling).
            const double Ls = Lm / (turnsRatio * turnsRatio);
            circuit << "Lp tx_pri_pos 0 " << std::scientific << Lm << std::fixed << " ic=0\n";
            circuit << "Ls 0 tx_sec_pos " << std::scientific << Ls << std::fixed << " ic=0\n";
            circuit << "K1 Lp Ls 0.999\n";
            circuit << "Cb tx_sec_pos cb_esr " << std::scientific << couplingCb << std::fixed
                    << " ic=" << VCb_DC << "\n";
            circuit << "Rcb_esr cb_esr node_B " << esrC1 << "\n\n";
        } else {
            // C1 coupling cap (from node_A_int through ESR to node_B). In the
            // ON sub-interval VC1 polarity is +(node_A) - (node_B) ≈ Vin/(1-D).
            circuit << "* C1 coupling cap (with ESR)\n";
            circuit << "Vc1_sense node_A_int node_C 0\n";
            circuit << "C1 node_C node_C_esr " << std::scientific << couplingCap << std::fixed
                    << " ic=" << VC1 << "\n";
            circuit << "Rc1_esr node_C_esr node_B " << esrC1 << "\n\n";
        }

        // D1 freewheel diode (V1/V2/V3 default) or S2 synchronous switch (V4/V5).
        // Conduction direction is identical: from node_B (positive during S1 OFF
        // because L2 forces continuity) into d_cath (= GND-via-sense). We always
        // route through the d_cath / Vd_sense node so the post-sim probes that
        // monitor i(Vd_sense) report the rectifier current uniformly across
        // variants.
        circuit << "Vd_sense d_cath 0 0\n";
        if (isSynchronous) {
            circuit << "* S2 synchronous rectifier (complementary PWM, no dead-time)\n";
            // pwm_ctrl_inv is high during the (1-D)·T sub-interval. We start it
            // high at t=0, go low at tOn (S1 turns on), back high at period.
            // Because PULSE() expresses ON-time, swap the levels and align edges.
            circuit << "Vpwm_inv pwm_ctrl_inv 0 PULSE(" << cfg.pwmHigh << " 0 0 "
                    << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                    << tOn << " " << period << std::fixed << ")\n";
            circuit << ".model SW2 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                    << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
            circuit << "S2 node_B d_cath pwm_ctrl_inv 0 SW2\n";
            circuit << "Rsnub_d1 node_B 0 " << cfg.snubR << "\n"
                    << "Csnub_d1 node_B snub_d1_int " << std::scientific << cfg.snubC << std::fixed << "\n"
                    << "Rsnub_d1b snub_d1_int 0 " << require_spice_field(cfg.snubDampR, "snubDampR") << "\n\n";
        } else {
            circuit << "* D1 freewheel diode\n";
            circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                    << " RS=" << cfg.diodeRS << std::fixed << ")\n";
            circuit << "D1 node_B d_cath DIDEAL\n";
            circuit << "Rsnub_d1 node_B 0 " << cfg.snubR << "\n"
                    << "Csnub_d1 node_B snub_d1_int " << std::scientific << cfg.snubC << std::fixed << "\n"
                    << "Rsnub_d1b snub_d1_int 0 " << require_spice_field(cfg.snubDampR, "snubDampR") << "\n\n";
        }

        // L2 output inductor (from node_B to vout_load_node, with DCR)
        circuit << "* L2 output inductor (with DCR)\n";
        circuit << "L2 node_B l2_dcr_mid " << std::scientific << inductanceL2 << std::fixed
                << " ic=" << IL2avg << "\n";
        circuit << "Rdcr_l2 l2_dcr_mid out_sense " << dcrL2 << "\n";
        circuit << "Vl2_sense out_sense vout_load_node 0\n\n";

        // V2 coupled-inductor: add mutual-inductance element K1 between L1 and L2.
        // ngspice uses the first listed node as the dotted terminal; with our
        // node ordering (L1 l1_in -> ..., L2 node_B -> ...) the L1 dot is on the
        // source side and the L2 dot is on the switching-node side, which is
        // the orientation Erickson §6 uses for input-current ripple cancellation.
        if (isCoupled) {
            // Default coupling coefficient: 0.999 (extremely tight, close to the
            // ideal zero-ripple limit). The user may override via the schema.
            double k = get_coupling_coefficient().value_or(0.999);
            if (k <= 0.0 || k >= 1.0) {
                throw std::invalid_argument(
                    "Cuk::generate_ngspice_circuit: couplingCoefficient must be "
                    "in (0, 1); received " + std::to_string(k));
            }
            circuit << "* V2 coupled-inductor: K1 couples L1 and L2 on the same core\n";
            circuit << "K1 L1 L2 " << k << "\n\n";
        }

        // Output cap and load. Vo is negative — initialise the cap and
        // load with sign matching that. The load resistor is connected
        // between vout_load_node (negative) and 0 (GND), so Iload
        // = (vout_load_node - 0)/Rload < 0 (current sinks from GND
        // into the load and back to vout_load_node, i.e. positive
        // current flows OUT of vout_load_node from the user's
        // perspective, since iLoad < 0 means it flows from 0 to
        // vout_load_node).
        circuit << "* Output cap (with ESR) and resistive load\n";
        circuit << "Cout vout_load_node co_esr " << std::scientific << outputCapacitance << std::fixed
                << " ic=" << (-outputVoltageMag) << "\n";
        circuit << "Rco_esr co_esr 0 " << esrCo << "\n";
        circuit << "Vout_sense vout_load_node vout_load 0\n";
        circuit << "Rload vout_load 0 " << loadResistance << "\n\n";

        // Bridge-voltage probe (CUK_PLAN.md §6.2): switch-node voltage
        circuit << "* Bridge / switch-node probe\n";
        circuit << "Bsw sw_probe 0 V=V(node_A)\n";
        // Differential L1 winding voltage probe (Boost-style §5.0):
        // V(l1_in) - V(node_A) gives the winding voltage with sign
        // matching i(Vin_sense) (current flows vin_dc → l1_in → ... → node_A).
        circuit << "Bvpri_diff vpri_diff 0 V=V(l1_in)-V(node_A)\n\n";

        // Transient analysis with UIC (initial conditions)
        circuit << "* Transient Analysis (UIC enables .ic above)\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime
                << std::fixed << " UIC\n\n";

        // Save signals (split per Boost.cpp §5.0):
        //   Winding-port stream: vpri_diff + Vin_sense#branch (= IL1)
        //   Converter-port stream: vin_dc + sw_probe + vout_load_node
        //                          + i(Vin_sense) + i(Vout_sense)
        circuit << "* Output signals\n";
        circuit << ".save v(vpri_diff) v(vin_dc) v(node_A) v(node_B) v(sw_probe) v(vout_load_node) "
                << "i(Vin_sense) i(Vl1_sense) i(Vc1_sense) i(Vl2_sense) i(Vd_sense) i(Vout_sense)\n\n";

        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::fixed
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        circuit << ".options METHOD=" << cfg.method << " TRTOL=" << cfg.trTol << "\n";

        circuit << ".end\n";
        return circuit.str();
    }

    std::vector<OperatingPoint> Cuk::simulate_and_extract_operating_points(double inductanceL1) {
        std::vector<OperatingPoint> operatingPoints;

        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("Cuk::simulate_and_extract_operating_points: ngspice is not available");
        }

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto opPoint = get_operating_points()[opIndex];
                std::string netlist = generate_ngspice_circuit(inductanceL1, inputVoltageIndex, opIndex);
                double switchingFrequency = opPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);
                if (!simResult.success) {
                    throw std::runtime_error("Cuk SPICE simulation failed: " + simResult.errorMessage);
                }

                NgspiceRunner::WaveformNameMapping waveformMapping = {
                    {{"voltage", "vpri_diff"}, {"current", "vin_sense#branch"}}
                };
                std::vector<std::string> windingNames = {"Primary"};
                std::vector<bool> flipCurrentSign = {false};

                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult, waveformMapping, switchingFrequency, windingNames,
                    opPoint.get_ambient_temperature(), flipCurrentSign);

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

    std::vector<ConverterWaveforms> Cuk::simulate_and_extract_topology_waveforms(double inductanceL1, size_t numberOfPeriods) {
        std::vector<ConverterWaveforms> results;

        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("Cuk::simulate_and_extract_topology_waveforms: ngspice is not available");
        }

        int originalNumPeriodsToExtract = get_num_periods_to_extract();
        set_num_periods_to_extract(static_cast<int>(numberOfPeriods));

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            (void) inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto opPoint = get_operating_points()[opIndex];

                std::string netlist = generate_ngspice_circuit(inductanceL1, inputVoltageIndex, opIndex);
                double switchingFrequency = opPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = numberOfPeriods;
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);
                if (!simResult.success) {
                    throw std::runtime_error("Cuk SPICE simulation failed: " + simResult.errorMessage);
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

                wf.set_input_voltage(getWaveform("vin_dc"));
                wf.set_input_current(getWaveform("vin_sense#branch"));
                wf.get_mutable_output_voltages().push_back(getWaveform("vout_load_node"));
                // Iout via Vout_sense branch — current flowing from the
                // negative output node into the load.
                wf.get_mutable_output_currents().push_back(getWaveform("vout_sense#branch"));

                results.push_back(wf);
            }
        }

        set_num_periods_to_extract(originalNumPeriodsToExtract);
        return results;
    }

    // ============================================================
    //   get_extra_components_inputs (CUK_PLAN.md §9 / §13)
    // ============================================================
    //
    // V1 ancillaries: L2 inductor, C1 coupling cap, Co output cap.
    //
    // Mirrors PhaseShiftedHalfBridge / Llc patterns:
    //   * Inductor → MAS::Inputs (DesignRequirements.magnetizing_inductance
    //     + per-OP per-winding excitation built from the stashed waveforms).
    //   * Capacitor → CAS::Inputs (DesignRequirements.capacitance,
    //     rated_voltage, role; one TwoTerminalOperatingPoint per stored
    //     waveform pair).
    //
    // Both modes (IDEAL and REAL) currently emit the same nominal sizing —
    // V1 does not absorb leakage / ESR into the main magnetic, so REAL is
    // identical to IDEAL.  When V2 (coupled-inductor) lands, REAL will
    // adjust the ancillary L2 down by the leakage already provided by the
    // primary 2-winding magnetic.
    std::vector<std::variant<Inputs, CAS::Inputs>> Cuk::get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic)
    {
        if (mode == ExtraComponentsMode::REAL && !magnetic.has_value()) {
            throw std::invalid_argument("Cuk::get_extra_components_inputs: mode REAL requires a designed magnetic");
        }
        const bool isIsolated = get_isolated().value_or(false);

        // V1/V2 ancillaries: L2, C1, Co  (3 entries).
        // V3 ancillaries: L1, L2, Ca, Cb, Co  (5 entries).
        if (lastSizedL2 <= 0.0 || lastSizedCo <= 0.0
            || extraL2VoltageWaveforms.empty()
            || extraCoVoltageWaveforms.empty())
        {
            throw std::runtime_error("Cuk::get_extra_components_inputs: call process_operating_points() first");
        }
        if (isIsolated) {
            if (lastSizedCa <= 0.0 || lastSizedCb <= 0.0
                || extraL1VoltageWaveforms.empty()
                || extraCaVoltageWaveforms.empty()
                || extraCbVoltageWaveforms.empty())
            {
                throw std::runtime_error("Cuk::get_extra_components_inputs: V3 stash incomplete; "
                                         "call process_operating_points() first");
            }
        } else if (lastSizedC1 <= 0.0 || extraC1VoltageWaveforms.empty()) {
            throw std::runtime_error("Cuk::get_extra_components_inputs: V1 stash incomplete; "
                                     "call process_operating_points() first");
        }

        const size_t nOps = extraL2VoltageWaveforms.size();
        auto checkSync = [nOps](const std::vector<Waveform>& wfms, const char* label) {
            if (wfms.size() != nOps) {
                throw std::runtime_error(std::string("Cuk::get_extra_components_inputs: ") +
                    label + " waveform count out of sync");
            }
        };
        checkSync(extraL2CurrentWaveforms, "L2 current");
        checkSync(extraCoVoltageWaveforms, "Co voltage");
        checkSync(extraCoCurrentWaveforms, "Co current");
        if (isIsolated) {
            checkSync(extraL1VoltageWaveforms, "L1 voltage");
            checkSync(extraL1CurrentWaveforms, "L1 current");
            checkSync(extraCaVoltageWaveforms, "Ca voltage");
            checkSync(extraCaCurrentWaveforms, "Ca current");
            checkSync(extraCbVoltageWaveforms, "Cb voltage");
            checkSync(extraCbCurrentWaveforms, "Cb current");
        } else {
            checkSync(extraC1VoltageWaveforms, "C1 voltage");
            checkSync(extraC1CurrentWaveforms, "C1 current");
        }

        // Use the first OP's switching frequency as the representative
        // frequency for ancillary characterisation (same convention as
        // PhaseShiftedHalfBridge::get_extra_components_inputs).
        const double fsw = get_operating_points()[0].get_switching_frequency();

        std::vector<std::variant<Inputs, CAS::Inputs>> result;

        // ---- Inductor builder ------------------------------------
        auto buildInductor = [&](const std::string& name, double sizedL,
                                 const std::vector<Waveform>& vWfms,
                                 const std::vector<Waveform>& iWfms)
        {
            Inputs masInputs;
            DesignRequirements dr;
            DimensionWithTolerance L;
            L.set_nominal(sizedL);
            dr.set_magnetizing_inductance(L);
            dr.set_name(name);
            dr.set_topology(MAS::Topology::CUK_CONVERTER);
            dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
            dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
            masInputs.set_design_requirements(dr);

            std::vector<OperatingPoint> masOps;
            for (size_t i = 0; i < vWfms.size(); ++i) {
                OperatingPoint op;
                auto excitation = complete_excitation(iWfms[i], vWfms[i], fsw, "Primary");
                op.get_mutable_excitations_per_winding().push_back(excitation);
                masOps.push_back(op);
            }
            masInputs.set_operating_points(masOps);
            result.emplace_back(std::move(masInputs));
        };

        // ---- Capacitor builder ------------------------------------
        auto buildCapacitor = [&](const std::string& name, double sizedC,
                                  const std::vector<Waveform>& vWfms,
                                  const std::vector<Waveform>& iWfms,
                                  CAS::Application role)
        {
            CAS::Inputs casInputs;
            CAS::DesignRequirements dr;
            CAS::DimensionWithTolerance capacitance;
            capacitance.set_nominal(sizedC);
            dr.set_capacitance(capacitance);
            double peakV = 0.0;
            for (const auto& wfm : vWfms) {
                for (double v : wfm.get_data()) peakV = std::max(peakV, std::abs(v));
            }
            dr.set_rated_voltage(peakV * 1.2);  // 20 % derating margin
            dr.set_role(role);
            casInputs.set_design_requirements(dr);

            std::vector<CAS::TwoTerminalOperatingPoint> casOps;
            for (size_t i = 0; i < vWfms.size(); ++i) {
                CAS::TwoTerminalOperatingPoint op;
                CAS::OperatingPointExcitation excitation;
                excitation.set_frequency(fsw);

                CAS::SignalDescriptor vSig;
                CAS::Waveform vWfm;
                vWfm.set_data(vWfms[i].get_data());
                if (vWfms[i].get_time()) vWfm.set_time(vWfms[i].get_time().value());
                vSig.set_waveform(vWfm);
                excitation.set_voltage(vSig);

                CAS::SignalDescriptor iSig;
                CAS::Waveform iWfm;
                iWfm.set_data(iWfms[i].get_data());
                if (iWfms[i].get_time()) iWfm.set_time(iWfms[i].get_time().value());
                iSig.set_waveform(iWfm);
                excitation.set_current(iSig);

                op.set_excitation(excitation);
                casOps.push_back(op);
            }
            casInputs.set_operating_points(casOps);
            result.emplace_back(std::move(casInputs));
        };

        if (isIsolated) {
            // V3: L1 (input choke), L2, Ca, Cb, Co — transformer is the "primary"
            // magnetic returned via process_design_requirements, not an extra.
            const double L1Choke = compute_l1_for_isolated();
            buildInductor("inputInductor",  L1Choke,    extraL1VoltageWaveforms, extraL1CurrentWaveforms);
            buildInductor("outputInductor", lastSizedL2, extraL2VoltageWaveforms, extraL2CurrentWaveforms);
            buildCapacitor("primaryCouplingCapacitor",   lastSizedCa,
                           extraCaVoltageWaveforms, extraCaCurrentWaveforms, CAS::Application::DC_LINK);
            buildCapacitor("secondaryCouplingCapacitor", lastSizedCb,
                           extraCbVoltageWaveforms, extraCbCurrentWaveforms, CAS::Application::DC_LINK);
            buildCapacitor("outputCapacitor",            lastSizedCo,
                           extraCoVoltageWaveforms, extraCoCurrentWaveforms, CAS::Application::OUTPUT_FILTER);
        } else {
            // V1/V2: L2, C1, Co — L1 IS the primary magnetic (DR.magnetizingInductance).
            buildInductor("outputInductor", lastSizedL2, extraL2VoltageWaveforms, extraL2CurrentWaveforms);
            buildCapacitor("couplingCapacitor", lastSizedC1,
                           extraC1VoltageWaveforms, extraC1CurrentWaveforms, CAS::Application::DC_LINK);
            buildCapacitor("outputCapacitor",   lastSizedCo,
                           extraCoVoltageWaveforms, extraCoCurrentWaveforms, CAS::Application::OUTPUT_FILTER);
        }

        return result;
    }

} // namespace OpenMagnetics
