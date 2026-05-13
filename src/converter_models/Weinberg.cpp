#include "converter_models/Weinberg.h"
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
    //   Static analytical helpers
    // ============================================================

    int Weinberg::detect_operating_regime(double dutyCycle) {
        constexpr double kBoundaryEps = 1e-6;
        if (dutyCycle < 0.5 - kBoundaryEps) return 0;       // buck-like
        if (dutyCycle > 0.5 + kBoundaryEps) return 2;       // boost-like
        return 1;                                            // boundary
    }

    double Weinberg::calculate_overlap_fraction(double dutyCycle) {
        return std::max(0.0, 2.0 * dutyCycle - 1.0);
    }

    double Weinberg::calculate_conversion_ratio_boost(double dutyCycle, double turnsRatio) {
        // Boost regime: M = 1 / (2 · n · (1 − D))            [Weinberg PESC 1985]
        if (turnsRatio <= 0.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_conversion_ratio_boost: turnsRatio must be > 0");
        }
        double oneMinusD = 1.0 - dutyCycle;
        if (oneMinusD <= 1e-6) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_conversion_ratio_boost: D too close to 1 — singular gain");
        }
        return 1.0 / (2.0 * turnsRatio * oneMinusD);
    }

    double Weinberg::calculate_conversion_ratio_buck(double dutyCycle, double turnsRatio) {
        if (turnsRatio <= 0.0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_conversion_ratio_buck: turnsRatio must be > 0");
        }
        return 2.0 * dutyCycle / turnsRatio;
    }

    double Weinberg::calculate_switch_peak_voltage_classic(double inputVoltage, double dutyCycle) {
        // V1 push-pull primary: V_Q,pk = 2 · Vin / (1 − D)   [clamped by D3]
        double oneMinusD = 1.0 - dutyCycle;
        if (oneMinusD <= 1e-6) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_switch_peak_voltage_classic: D too close to 1");
        }
        return 2.0 * inputVoltage / oneMinusD;
    }

    double Weinberg::calculate_switch_peak_voltage_bridge(double inputVoltage) {
        return inputVoltage;
    }

    double Weinberg::calculate_duty_cycle(double inputVoltage, double outputVoltage,
                                          double turnsRatio, double diodeVoltageDrop,
                                          double efficiency, double maximumDutyCycle) {
        // Solve M(D) = (Vout + Vd) / (Vin · η) for D.
        // Try the boost regime first (Weinberg's canonical operating regime,
        // D > 0.5):  M_boost = 1 / (2·n·(1−D))
        //   → D = 1 − 1 / (2·n·M_boost)
        // If the resulting D > 0.5 it is the valid boost solution; else the
        // converter is operating in the buck regime (D < 0.5):
        //   M_buck = 2·D / n  → D = n · M_buck / 2
        if (inputVoltage <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_duty_cycle: input voltage must be > 0");
        }
        if (outputVoltage <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_duty_cycle: outputVoltage must be > 0");
        }
        if (turnsRatio <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_duty_cycle: turnsRatio must be > 0");
        }
        double vinEff = inputVoltage * efficiency;
        double M = (outputVoltage + diodeVoltageDrop) / vinEff;
        double dBoost = 1.0 - 1.0 / (2.0 * turnsRatio * M);
        if (dBoost > 0.5) {
            const double dutyTolerance = 0.01;
            if (dBoost >= maximumDutyCycle - dutyTolerance) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Weinberg::calculate_duty_cycle: D=" + std::to_string(dBoost) +
                    " exceeds maximumDutyCycle " + std::to_string(maximumDutyCycle) +
                    " — reduce Vo, raise Vin, or raise maximumDutyCycle");
            }
            return dBoost;
        }
        double dBuck = 0.5 * turnsRatio * M;
        if (dBuck >= 0.5) {
            // Boundary case: exactly D = 0.5 satisfies both branches at M = 1/n.
            return 0.5;
        }
        if (dBuck <= 1e-3) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_duty_cycle: D=" + std::to_string(dBuck) +
                " <= 0.001 — converter would lose regulation; raise Vo or lower Vin");
        }
        return dBuck;
    }

    double Weinberg::calculate_l1_min(double inputVoltage, double dutyCycle, double deltaIL1, double switchingFrequency) {
        // Boost regime (D > 0.5): L1 ≥ Vin · (2D − 1) · Tsw / (4 · ΔI_L1)
        // (per WEINBERG_PLAN.md S1; the 4 in the denominator reflects two
        // overlap intervals per period and the per-winding ripple budget.)
        // For D ≤ 0.5 the input inductor sees only the dead/non-overlap
        // intervals — fall back to the classical boost  L1 = Vin·D·Tsw / ΔI.
        if (deltaIL1 <= 0 || switchingFrequency <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_l1_min: deltaIL1 and fsw must be > 0");
        }
        double overlap = calculate_overlap_fraction(dutyCycle);
        double dEffective = (overlap > 0.0) ? overlap : dutyCycle;
        return inputVoltage * dEffective / (4.0 * deltaIL1 * switchingFrequency);
    }

    double Weinberg::calculate_dcm_K_boost(double L1, double switchingFrequency, double dutyCycle, double loadResistance) {
        // K_crit (boost) = 2 · L1 · fsw · (1 − D)² / R                   [WEINBERG_PLAN B1]
        if (L1 <= 0 || switchingFrequency <= 0 || loadResistance <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Weinberg::calculate_dcm_K_boost: L1, fsw, R must all be > 0");
        }
        double oneMinusD = 1.0 - dutyCycle;
        return 2.0 * L1 * switchingFrequency * oneMinusD * oneMinusD / loadResistance;
    }

    double Weinberg::calculate_rhp_zero_freq(double loadResistance, double dutyCycle, double turnsRatio, double L1) {
        // ω_RHPZ ≈ R · (1 − D)² / (n² · L1)                              [Erickson §8]
        if (loadResistance <= 0 || turnsRatio <= 0 || L1 <= 0) {
            return 0.0;
        }
        double oneMinusD = 1.0 - dutyCycle;
        double omega = loadResistance * oneMinusD * oneMinusD / (turnsRatio * turnsRatio * L1);
        return omega / (2.0 * std::numbers::pi);
    }

    // ============================================================
    //   Constructors
    // ============================================================

    Weinberg::Weinberg(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    AdvancedWeinberg::AdvancedWeinberg(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    // ============================================================
    //   Per-operating-point analytical worker
    // ============================================================

    OperatingPoint Weinberg::process_operating_points_for_input_voltage(double inputVoltage,
        const TopologyExcitation& outputOperatingPoint,
        double turnsRatio,
        double magnetizingInductance)
    {
        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double outputVoltage      = outputOperatingPoint.get_output_voltages()[0];
        double outputCurrent      = outputOperatingPoint.get_output_currents()[0];
        double diodeVoltageDrop   = get_diode_voltage_drop();
        double efficiency = 1.0;
        if (get_efficiency()) efficiency = get_efficiency().value();

        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, turnsRatio,
                                                diodeVoltageDrop, efficiency, maximumDutyCycle.value_or(0.95));
        int regime = detect_operating_regime(dutyCycle);
        double overlap = calculate_overlap_fraction(dutyCycle);

        double M = (regime == 2) ? calculate_conversion_ratio_boost(dutyCycle, turnsRatio)
                                 : calculate_conversion_ratio_buck(dutyCycle, turnsRatio);

        // Average input current from output power and efficiency.
        double inputCurrent = (M > 0) ? outputCurrent / (efficiency * M) : 0.0;
        double IL1perWinding = inputCurrent / 2.0;          // 1:1 input coupled-L splits Iin

        double inductanceL1 = magnetizingInductance;
        double deltaIL1 = inductanceL1 > 0 ? (inputVoltage * std::max(overlap, dutyCycle))
                                              / (inductanceL1 * switchingFrequency)
                                            : 0.0;

        // Variant decides the switch-V_DS expression; default to V1 classic.
        const MAS::Variant variant = get_variant().value_or(MAS::Variant::CLASSIC);
        double switchPeakVoltage = (variant == MAS::Variant::BRIDGE)
            ? calculate_switch_peak_voltage_bridge(inputVoltage)
            : calculate_switch_peak_voltage_classic(inputVoltage, dutyCycle);

        // Switch peak current — non-overlap interval dominates.
        double switchPeakCurrent = inputCurrent + deltaIL1 / 2.0;

        // CT-FW rectifier diode peak reverse voltage = 2·Vout (per CT-FW
        // standard); diode peak current ≈ output current (continuous via
        // Co + reflected L1).
        double diodePeakReverseVoltage = 2.0 * outputVoltage;
        double diodePeakCurrent = outputCurrent + deltaIL1 / (2.0 * std::max(turnsRatio, 1e-9));

        // Energy-recovery diode D3 average current (V1 only; rough estimate).
        double energyRecoveryAvgCurrent = (variant == MAS::Variant::CLASSIC)
            ? inputCurrent * std::max(2.0 * dutyCycle - 1.0, 0.0) * leakageInductanceFraction
            : 0.0;

        // Output cap sized from continuous-current ripple model.
        double loadResistance = (outputCurrent > 0) ? outputVoltage / outputCurrent : 0.0;
        double deltaVo_target = std::max(coRipplePct * outputVoltage, 1e-3);
        double outputCapacitance = (loadResistance > 0)
            ? outputCurrent * dutyCycle / (deltaVo_target * switchingFrequency)
            : 1e-6;
        outputCapacitance = std::max(outputCapacitance, 1e-6);
        double deltaVo_actual = outputCurrent * dutyCycle / (outputCapacitance * switchingFrequency);

        // ---- Diagnostics ----
        lastDutyCycle = dutyCycle;
        lastConversionRatio = M;
        lastOperatingRegime = regime;
        lastOverlapFraction = overlap;
        lastSwitchPeakVoltage = switchPeakVoltage;
        lastSwitchPeakCurrent = switchPeakCurrent;
        lastDiodePeakReverseVoltage = diodePeakReverseVoltage;
        lastDiodePeakCurrent = diodePeakCurrent;
        lastEnergyRecoveryAvgCurrent = energyRecoveryAvgCurrent;
        lastInputInductorAverage = IL1perWinding;
        lastInputInductorRipple = deltaIL1;
        lastMagnetizingRipple = deltaIL1 / std::max(turnsRatio, 1e-9);
        lastFluxImbalanceMargin = 0.0;            // symmetric drive assumed in v1
        lastRhpZeroFrequency = (regime == 2)
            ? calculate_rhp_zero_freq(loadResistance, dutyCycle, turnsRatio, inductanceL1)
            : 0.0;
        if (loadResistance > 0 && inductanceL1 > 0) {
            lastDcmK = calculate_dcm_K_boost(inductanceL1, switchingFrequency, dutyCycle, loadResistance);
            lastDcmKcrit = 1.0;                    // K_crit normalised so K > 1 → CCM in boost
            lastIsCcm = (lastDcmK > 1.0);
        } else {
            lastDcmK = 0.0;
            lastDcmKcrit = 0.0;
            lastIsCcm = true;
        }
        lastSizedL1 = inductanceL1;
        lastSizedCo = outputCapacitance;
        lastOutputVoltageRipple = deltaVo_actual;

        // ---- Primary-winding excitation (combined primary, 2-winding xfmr).
        // The SPICE harness probes one half-winding (Vpri_sense_a in Lpri_a),
        // so the analytical waveform must match that shape. Three regimes:
        //  • Buck (D < 0.5):  Q1-only pulse over [0, D·T], 0 elsewhere.
        //  • Boundary (D=0.5): Q1-only pulse over [0, T/2], 0 over [T/2, T].
        //  • Boost (D > 0.5): overlap segments where both Q1, Q2 conduct
        //                     and half-A carries ~iL1/2; Q1-only segment
        //                     in the middle carries full iL1; half-A is OFF
        //                     during Q2-only segment in the second half.
        const double period = 1.0 / switchingFrequency;
        const double iL1_low  = std::max(inputCurrent - 0.5 * deltaIL1, 0.0);
        const double iL1_high = inputCurrent + 0.5 * deltaIL1;
        const double iL1_half_low  = 0.5 * iL1_low;
        const double iL1_half_high = 0.5 * iL1_high;
        {
            Waveform iPri;
            std::vector<double> iData;
            std::vector<double> iTime;
            if (regime == 2) {
                // Boost: 4 segments — overlap1 / Q1-only / overlap2 / Q2-only.
                const double tOverlap = (2.0 * dutyCycle - 1.0) * period / 2.0;
                const double t1 = tOverlap;
                const double t2 = period / 2.0;
                const double t3 = period / 2.0 + tOverlap;
                iData = { iL1_half_low,  iL1_half_high,           // overlap-1 ramp
                          iL1_low,       iL1_high,                // Q1-only ramp
                          iL1_half_high, iL1_half_low,            // overlap-2 ramp (down)
                          0.0,           0.0 };                   // Q2-only (off)
                iTime = { 0.0, t1,
                          t1, t2,
                          t2, t3,
                          t3, period };
            } else {
                // Buck or boundary: pulse during Q1-on, 0 elsewhere.
                const double tOn = dutyCycle * period;
                iData = { iL1_low, iL1_high, 0.0, 0.0 };
                iTime = { 0.0, tOn, tOn, period };
            }
            iPri.set_ancillary_label(WaveformLabel::CUSTOM);
            iPri.set_data(iData);
            iPri.set_time(iTime);

            Waveform vPri = Inputs::create_waveform(
                WaveformLabel::BIPOLAR_RECTANGULAR, 2.0 * inputVoltage,
                switchingFrequency, 0.5, 0.0, 0);
            auto excitation = complete_excitation(iPri, vPri, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // ---- L1 input coupled-L (per winding) waveforms ----
        // Voltage on each L1 winding swings between +Vin/2 (during non-
        // overlap, when only the opposite primary half is conducting and
        // L1 is being charged through the coupling) and a discharge level;
        // we model it as a triangular voltage with peak-to-peak Vin and a
        // mean equal to (Vin/M)/2 reflected via the input current path.
        // Current: triangular around IL1perWinding, ripple ΔI_L1.
        {
            Waveform vL1 = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, inputVoltage, switchingFrequency, 0.5, 0.0, 0);
            Waveform iL1 = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaIL1, switchingFrequency, 0.5,
                IL1perWinding, 0);
            extraL1VoltageWaveforms.push_back(vL1);
            extraL1CurrentWaveforms.push_back(iL1);
        }

        // ---- Co waveforms ----
        {
            Waveform vCo = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaVo_actual, switchingFrequency, 0.5,
                outputVoltage, 0);
            // Co AC current = rectified secondary − Iout; small triangular
            // approximation around 0 mean with peak-to-peak ≈ ΔIL1/n.
            double iCoPP = deltaIL1 / std::max(turnsRatio, 1e-9);
            Waveform iCo = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, iCoPP, switchingFrequency, 0.5, 0.0, 0);
            extraCoVoltageWaveforms.push_back(vCo);
            extraCoCurrentWaveforms.push_back(iCo);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    bool Weinberg::run_checks(bool assert) {
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

    DesignRequirements Weinberg::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double diodeVoltageDrop    = get_diode_voltage_drop();
        double efficiency = 1.0;
        if (get_efficiency()) efficiency = get_efficiency().value();

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument(
                "Weinberg::process_design_requirements: missing both currentRippleRatio "
                "and maximumSwitchCurrent — at least one is needed to size L1");
        }

        // ----- Turns ratio (single combined Pri:Sec ratio) -----
        // Size n at MAXIMUM Vin (worst case for boost-regime D-min): require
        // D ≥ 0.55 at maxVin so the converter stays in the boost regime
        // across the input-voltage range.  D = 1 − 1/(2·n·M)  →
        //    n = 1 / (2 · M · (1 − D_target))      where M = (Vo + Vd)/(Vin·η)
        constexpr double dTarget = 0.55;
        double maximumTurnsRatio = 0.0;
        for (const auto& op : get_operating_points()) {
            double Vo = op.get_output_voltages()[0];
            double M = (Vo + diodeVoltageDrop) / (maximumInputVoltage * efficiency);
            double n = 1.0 / (2.0 * M * (1.0 - dTarget));
            maximumTurnsRatio = std::max(maximumTurnsRatio, n);
        }
        if (maximumTurnsRatio <= 0.0) {
            throw std::invalid_argument(
                "Weinberg::process_design_requirements: derived turnsRatio is non-positive");
        }

        // ----- L1 sizing — worst case: minimum Vin, maximum Iout -----
        double maximumDeltaIL1 = 0.0;
        if (get_current_ripple_ratio()) {
            double rippleRatio = get_current_ripple_ratio().value();
            for (const auto& op : get_operating_points()) {
                double Iout = op.get_output_currents()[0];
                double Vo   = op.get_output_voltages()[0];
                double D    = calculate_duty_cycle(minimumInputVoltage, Vo, maximumTurnsRatio,
                                                   diodeVoltageDrop, efficiency, maximumDutyCycle.value_or(0.95));
                double M    = calculate_conversion_ratio_boost(D, maximumTurnsRatio);
                double Iin  = Iout / (efficiency * M);
                double IL1perWinding = Iin / 2.0;
                maximumDeltaIL1 = std::max(maximumDeltaIL1, rippleRatio * IL1perWinding);
            }
        }
        if (get_maximum_switch_current()) {
            double IsMax = get_maximum_switch_current().value();
            for (const auto& op : get_operating_points()) {
                double Iout = op.get_output_currents()[0];
                double Vo   = op.get_output_voltages()[0];
                double D    = calculate_duty_cycle(minimumInputVoltage, Vo, maximumTurnsRatio,
                                                   diodeVoltageDrop, efficiency, maximumDutyCycle.value_or(0.95));
                double M    = calculate_conversion_ratio_boost(D, maximumTurnsRatio);
                double Iin  = Iout / (efficiency * M);
                double residual = IsMax - Iin;
                double deltaIL1 = std::max(2.0 * residual, 0.0);
                maximumDeltaIL1 = std::max(maximumDeltaIL1, deltaIL1);
            }
        }
        if (maximumDeltaIL1 <= 0) {
            throw std::invalid_argument(
                "Weinberg::process_design_requirements: derived ΔI_L1 budget is non-positive");
        }

        double maximumNeededInductance = 0.0;
        for (const auto& op : get_operating_points()) {
            double switchingFrequency = op.get_switching_frequency();
            double Vo = op.get_output_voltages()[0];
            double D  = calculate_duty_cycle(minimumInputVoltage, Vo, maximumTurnsRatio,
                                             diodeVoltageDrop, efficiency, maximumDutyCycle.value_or(0.95));
            double L1 = calculate_l1_min(minimumInputVoltage, D, maximumDeltaIL1, switchingFrequency);
            maximumNeededInductance = std::max(maximumNeededInductance, L1);
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        DimensionWithTolerance turnsRatioWithTolerance;
        turnsRatioWithTolerance.set_nominal(roundFloat(maximumTurnsRatio, 3));
        designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        // 2-winding xfmr (combined Pri, combined Sec) — 2 isolation sides.
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        isolationSides.push_back(get_isolation_side_from_index(1));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::WEINBERG_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Weinberg::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        if (turnsRatios.empty()) {
            throw std::invalid_argument(
                "Weinberg::process_operating_points: turnsRatios is empty (need n = Np/Ns)");
        }
        double turnsRatio = turnsRatios[0];

        extraL1VoltageWaveforms.clear();
        extraL1CurrentWaveforms.clear();
        extraCoVoltageWaveforms.clear();
        extraCoCurrentWaveforms.clear();

        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            double inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(
                    inputVoltage, get_operating_points()[opIndex], turnsRatio, magnetizingInductance);

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

    std::vector<OperatingPoint> Weinberg::process_operating_points(OpenMagnetics::Magnetic magnetic) {
        run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(
            magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();

        std::vector<double> turnsRatios = magnetic.get_turns_ratios();

        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedWeinberg::process() {
        Weinberg::run_checks(_assertErrors);

        extraL1VoltageWaveforms.clear();
        extraL1CurrentWaveforms.clear();
        extraCoVoltageWaveforms.clear();
        extraCoCurrentWaveforms.clear();

        Inputs inputs;
        double maximumNeededInductance = get_desired_inductance();
        double turnsRatio = get_desired_turns_ratio();
        if (turnsRatio <= 0.0) {
            throw std::invalid_argument("AdvancedWeinberg::process: desiredTurnsRatio must be > 0");
        }

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        DesignRequirements designRequirements;
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        DimensionWithTolerance turnsRatioWithTolerance;
        turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 3));
        designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        isolationSides.push_back(get_isolation_side_from_index(1));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::WEINBERG_CONVERTER);
        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            double inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(
                    inputVoltage, get_operating_points()[opIndex], turnsRatio, maximumNeededInductance);

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
    //   SPICE netlist
    // ============================================================
    //
    // V1 (classic) node convention:
    //   vin_dc → Vin_sense → vin_p
    //   vin_p →[L1a]→ priCT_a   (input coupled-L upper winding, ic=Iin/2)
    //   vin_p →[L1b]→ priCT_b   (input coupled-L lower winding, ic=Iin/2)
    //   K_in  L1a L1b 0.999     (1:1 input coupled inductor)
    //
    //   priCT_a →[Lpri_a]→[Llk_pa]→ drainQ1
    //   priCT_b →[Lpri_b]→[Llk_pb]→ drainQ2
    //   secCT  →[Lsec_a]→[Llk_sa]→ diodePos
    //   secCT  →[Lsec_b]→[Llk_sb]→ diodeNeg
    //   K_xfmr Lpri_a Lpri_b Lsec_a Lsec_b 0.99
    //
    //   drainQ1 →[D3a]→ vin_p   (energy-recovery clamp, V1 only)
    //   drainQ2 →[D3b]→ vin_p
    //   drainQ1 →[Q1: SW1 to GND, gated by pwm1]
    //   drainQ2 →[Q2: SW1 to GND, gated by pwm2]
    //
    //   diodePos →[D_pos]→ out_node
    //   diodeNeg →[D_neg]→ out_node
    //   secCT → 0    (CT secondary mid-point grounded)
    //
    //   out_node → Co → 0;  out_node →[Vout_sense]→ Rload → 0
    //
    // V2 (bridge): replace the push-pull primary with an H-bridge.  The
    // primary has a SINGLE winding Lpri; the input coupled-L still has 2
    // windings 1:1.  Switches QA1/QA2 (Vin → priLeft → 0) and QB1/QB2
    // (Vin → priRight → 0) form the bridge.  D3 is omitted (the body
    // diodes of the bridge handle the leakage path).
    // ============================================================

    std::string Weinberg::generate_ngspice_circuit(double turnsRatio, double magnetizingInductance,
        size_t inputVoltageIndex, size_t operatingPointIndex)
    {
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames_;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames_);

        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("Weinberg::generate_ngspice_circuit: inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("Weinberg::generate_ngspice_circuit: operatingPointIndex out of range");
        }

        double inputVoltage = inputVoltages[inputVoltageIndex];
        auto opPoint = get_operating_points()[operatingPointIndex];

        double outputVoltage      = opPoint.get_output_voltages()[0];
        double outputCurrent      = opPoint.get_output_currents()[0];
        double switchingFrequency = opPoint.get_switching_frequency();
        double diodeVoltageDrop   = get_diode_voltage_drop();
        double efficiency = 1.0;
        if (get_efficiency()) efficiency = get_efficiency().value();

        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, turnsRatio,
                                                diodeVoltageDrop, efficiency, maximumDutyCycle.value_or(0.95));
        int regime = detect_operating_regime(dutyCycle);
        double M = (regime == 2) ? calculate_conversion_ratio_boost(dutyCycle, turnsRatio)
                                 : calculate_conversion_ratio_buck(dutyCycle, turnsRatio);
        double inputCurrent = (M > 0) ? outputCurrent / (efficiency * M) : 0.0;

        if (outputCurrent <= 0.0) {
            throw std::invalid_argument(
                "Weinberg::generate_ngspice_circuit: outputCurrent must be > 0 "
                "(received " + std::to_string(outputCurrent) +
                "); cannot derive load resistance R_load = Vo / Iout.");
        }
        double loadResistance = outputVoltage / outputCurrent;

        double inductanceL1 = magnetizingInductance;
        double Lpri_half = inductanceL1;                  // mirror Lpri = L1 magnitude
        double Lsec_half = Lpri_half / (turnsRatio * turnsRatio);
        double Llk_pri   = Lpri_half * leakageInductanceFraction;
        double Llk_sec   = Lsec_half * leakageInductanceFraction;
        Lpri_half -= Llk_pri;
        Lsec_half -= Llk_sec;
        double couplingMain  = get_coupling_coefficient_main().value_or(0.9999);
        double couplingInput = get_coupling_coefficient_input().value_or(0.999);
        if (couplingMain <= 0.0 || couplingMain >= 1.0) {
            throw std::invalid_argument(
                "Weinberg::generate_ngspice_circuit: couplingCoefficientMain must be in (0,1)");
        }
        if (couplingInput <= 0.0 || couplingInput >= 1.0) {
            throw std::invalid_argument(
                "Weinberg::generate_ngspice_circuit: couplingCoefficientInput must be in (0,1)");
        }

        double deltaVo_target = std::max(coRipplePct * outputVoltage, 1e-3);
        double outputCapacitance = outputCurrent * dutyCycle / (deltaVo_target * switchingFrequency);
        outputCapacitance = std::max(outputCapacitance, 1e-6);

        const SpiceSimulationConfig cfg = spice_config();
        const bool isSynchronous = get_synchronous_rectifier().value_or(false);
        const MAS::Variant variant = get_variant().value_or(MAS::Variant::CLASSIC);
        const bool isBridge = (variant == MAS::Variant::BRIDGE);

        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        double halfPeriod = period / 2.0;

        int periodsToExtract = get_num_periods_to_extract();
        int settlingPeriods  = get_num_steady_state_periods();
        const int numPeriodsTotal = settlingPeriods + periodsToExtract;
        double simTime  = numPeriodsTotal * period;
        double startTime = settlingPeriods * period;
        double samplesPerPeriod = (switchingFrequency >= 1e6) ? 500.0 : cfg.samplesPerPeriod;
        double stepTime = period / samplesPerPeriod;

        const double dcrL1 = 50e-3;
        const double esrCo = 5e-3;

        circuit << "* Weinberg Converter ("
                << (isBridge ? "V2 bridge primary" : "V1 classic push-pull")
                << (isSynchronous ? "/synchronous SR" : "/CT-FW diode")
                << ") - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vo=" << outputVoltage << "V, "
                << "f=" << (switchingFrequency / 1e3) << "kHz, D=" << (dutyCycle * 100) << " pct, "
                << "n=" << turnsRatio << ", regime="
                << (regime == 0 ? "buck" : (regime == 1 ? "boundary" : "boost")) << "\n";
        circuit << "* L1=" << (inductanceL1 * 1e6) << "uH, Lpri_half=" << (Lpri_half * 1e6)
                << "uH, Lsec_half=" << (Lsec_half * 1e6) << "uH, Co=" << (outputCapacitance * 1e6)
                << "uF, Iout=" << outputCurrent << "A, Iin=" << inputCurrent << "A\n\n";

        // DC Input + sense
        circuit << "* DC Input + sense\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n";
        circuit << "Vin_sense vin_dc vin_p 0\n\n";

        // Models
        circuit << "* Switch and diode models\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                << " RON=0.01 ROFF=1Meg\n";
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::fixed << ")\n\n";

        // PWM sources — Q1 phase 0, Q2 phase Tsw/2.  Overlap (when D > 0.5)
        // arises naturally from the pulse widths.
        circuit << "* PWM sources (Q1 phase=0, Q2 phase=Tsw/2)\n";
        circuit << "Vpwm1 pwm1 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                << tOn << " " << period << std::fixed << ")\n";
        circuit << "Vpwm2 pwm2 0 PULSE(0 " << cfg.pwmHigh << " "
                << std::scientific << halfPeriod << " " << cfg.pwmRise << " " << cfg.pwmFall << " "
                << tOn << " " << period << std::fixed << ")\n\n";

        // Input coupled inductor L1 (1:1, two windings on same core).
        double IL1 = inputCurrent / 2.0;
        circuit << "* Input coupled inductor L1 (1:1, two windings, k=0.999)\n";
        circuit << "L1a vin_p l1a_dcr_mid " << std::scientific << inductanceL1 << std::fixed
                << " ic=" << IL1 << "\n";
        circuit << "Rdcr_l1a l1a_dcr_mid priCT_a " << dcrL1 << "\n";
        circuit << "L1b vin_p l1b_dcr_mid " << std::scientific << inductanceL1 << std::fixed
                << " ic=" << IL1 << "\n";
        circuit << "Rdcr_l1b l1b_dcr_mid priCT_b " << dcrL1 << "\n";
        circuit << "K_in L1a L1b " << couplingInput << "\n\n";

        if (isBridge) {
            // V2 H-bridge primary.  We collapse priCT_a and priCT_b to a
            // single primary node by tying them together (the H-bridge
            // sees the *summed* L1a+L1b current through one Lpri winding).
            circuit << "* V2 bridge primary: tie priCT_a and priCT_b at bridge top\n";
            circuit << "Vbridge_top priCT_a priCT_b 0\n\n";

            circuit << "* H-bridge (4 switches: SA1/SA2 on left leg, SB1/SB2 on right leg)\n";
            // SA1: priCT_a → priLeft (closed when pwm1 HIGH)
            // SA2: priLeft → 0       (closed when pwm2 HIGH)
            // SB1: priCT_a → priRight (closed when pwm2 HIGH)
            // SB2: priRight → 0       (closed when pwm1 HIGH)
            circuit << "SA1 priCT_a priLeft  pwm1 0 SW1\n";
            circuit << "SA2 priLeft 0       pwm2 0 SW1\n";
            circuit << "SB1 priCT_a priRight pwm2 0 SW1\n";
            circuit << "SB2 priRight 0       pwm1 0 SW1\n";
            // Snubbers
            circuit << "Rsnub_sa1 priCT_a priLeft "  << cfg.snubR << "\n";
            circuit << "Csnub_sa1 priLeft sn_sa1 "   << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "Rsnub_sa1b sn_sa1 priCT_a 0.001\n";
            circuit << "Rsnub_sb2 priRight 0 "       << cfg.snubR << "\n";
            circuit << "Csnub_sb2 priRight sn_sb2 " << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "Rsnub_sb2b sn_sb2 0 0.001\n\n";

            // Single primary winding Lpri across priLeft → priRight.
            circuit << "* Main transformer primary (single winding for V2)\n";
            circuit << "Vpri_sense priLeft pri_sense_mid 0\n";
            circuit << "Llk_pri pri_sense_mid pri_top " << std::scientific << Llk_pri << std::fixed << "\n";
            circuit << "Lpri pri_top priRight " << std::scientific << (2.0 * Lpri_half) << std::fixed << "\n\n";

            // For SPICE coupling we still use 4 windings on the secondary
            // side (Lsec_a, Lsec_b for CT-FW), but couple them to Lpri only.
            circuit << "* Main transformer secondary (CT-FW, mixed dot convention: a@diode, b@CT)\n";
            circuit << "Lsec_a secCT_a_mid secCT " << std::scientific << Lsec_half << std::fixed << "\n";
            circuit << "Llk_sa secCT_a_mid diodePos " << std::scientific << Llk_sec << std::fixed << "\n";
            circuit << "Lsec_b secCT secCT_b_mid " << std::scientific << Lsec_half << std::fixed << "\n";
            circuit << "Llk_sb secCT_b_mid diodeNeg " << std::scientific << Llk_sec << std::fixed << "\n";
            // V2 main xfmr coupling: 3 windings → 3 pairwise K statements.
            circuit << "K_pri_sa Lpri Lsec_a " << couplingMain << "\n";
            circuit << "K_pri_sb Lpri Lsec_b " << couplingMain << "\n";
            circuit << "K_sa_sb Lsec_a Lsec_b " << couplingMain << "\n\n";

            circuit << "* Bridge / Vab probe (V_DS of one switch = Vin)\n";
            circuit << "Bvab vab 0 V=V(priLeft)-V(priRight)\n\n";
        } else {
            // V1 push-pull primary.
            // Dot convention: real CT-PP windings are wound in opposite
            // physical directions, so dots fall on OPPOSITE ends:
            //   Lpri_a: dot at drainQ1 (switch) — "Lpri_a drainQ1 pri_a_top"
            //   Lpri_b: dot at pri_b_top (CT)   — "Lpri_b pri_b_top drainQ2"
            //   Lsec_a: dot at secCT (CT)       — "Lsec_a secCT secCT_a_mid"
            //   Lsec_b: dot at secCT_b_mid (D)  — "Lsec_b secCT_b_mid secCT"
            // With ALL pairwise K positive this gives:
            //   S1 ON  → drainQ2 = +2·Vin, secCT_a_mid = +n·Vin (D_pos ON,
            //            D_neg reverse), secCT_b_mid = −n·Vin
            //   S2 ON  → drainQ1 = +2·Vin, secCT_b_mid = +n·Vin (D_neg ON,
            //            D_pos reverse), secCT_a_mid = −n·Vin
            // Anything else either reverse-biases both diodes or short-
            // circuits the off-side switch through the magnetizing path.
            circuit << "* Main transformer primary (CT push-pull, opposite-dot windings)\n";
            circuit << "Vpri_sense_a priCT_a pri_a_sense_mid 0\n";
            circuit << "Llk_pa pri_a_sense_mid pri_a_top " << std::scientific << Llk_pri << std::fixed << "\n";
            circuit << "Lpri_a drainQ1 pri_a_top " << std::scientific << Lpri_half << std::fixed << "\n";
            circuit << "Vpri_sense_b priCT_b pri_b_sense_mid 0\n";
            circuit << "Llk_pb pri_b_sense_mid pri_b_top " << std::scientific << Llk_pri << std::fixed << "\n";
            circuit << "Lpri_b pri_b_top drainQ2 " << std::scientific << Lpri_half << std::fixed << "\n\n";

            circuit << "* Main transformer secondary (CT-FW, opposite-dot windings)\n";
            circuit << "Lsec_a secCT secCT_a_mid " << std::scientific << Lsec_half << std::fixed << "\n";
            circuit << "Llk_sa secCT_a_mid diodePos " << std::scientific << Llk_sec << std::fixed << "\n";
            circuit << "Lsec_b secCT_b_mid secCT " << std::scientific << Lsec_half << std::fixed << "\n";
            circuit << "Llk_sb secCT_b_mid diodeNeg " << std::scientific << Llk_sec << std::fixed << "\n";
            // V1 main xfmr coupling: 4 windings → 6 pairwise K statements.
            circuit << "K_pa_pb Lpri_a Lpri_b " << couplingMain << "\n";
            circuit << "K_pa_sa Lpri_a Lsec_a " << couplingMain << "\n";
            circuit << "K_pa_sb Lpri_a Lsec_b " << couplingMain << "\n";
            circuit << "K_pb_sa Lpri_b Lsec_a " << couplingMain << "\n";
            circuit << "K_pb_sb Lpri_b Lsec_b " << couplingMain << "\n";
            circuit << "K_sa_sb Lsec_a Lsec_b " << couplingMain << "\n\n";

            circuit << "* Push-pull switches S1, S2 (2 SW1) + snubbers\n";
            circuit << "S1 drainQ1 0 pwm1 0 SW1\n";
            circuit << "S2 drainQ2 0 pwm2 0 SW1\n";
            circuit << "Rsnub_s1 drainQ1 0 " << cfg.snubR << "\n";
            circuit << "Csnub_s1 drainQ1 sn_s1 " << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "Rsnub_s1b sn_s1 0 0.001\n";
            circuit << "Rsnub_s2 drainQ2 0 " << cfg.snubR << "\n";
            circuit << "Csnub_s2 drainQ2 sn_s2 " << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "Rsnub_s2b sn_s2 0 0.001\n\n";

            // NOTE: Schreuders' D3 energy-recovery diodes (drain→Vin) would
            // clamp the legitimate 2·Vin/(1−D) push-pull drain swing to Vin+0.7
            // and destroy power transfer. The literature topology relies on an
            // auxiliary clamp winding (or a separate clamp capacitor) that we
            // do not model here. The RC snubbers above suppress the leakage
            // spike adequately for SPICE convergence; D3 is therefore omitted
            // from V1 SPICE. Future enhancement: add aux-winding clamp.
            circuit << "* (Energy-recovery D3 omitted — see comment in Weinberg.cpp)\n\n";

            circuit << "* Bridge / Vab probe (V_DS of S1 = 2·Vin/(1−D))\n";
            circuit << "Bvab vab 0 V=V(drainQ1)-V(drainQ2)\n\n";
        }

        // Secondary CT mid-point grounded.
        circuit << "* Secondary center-tap returns to GND\n";
        circuit << "Vsec_ct_sense secCT 0 0\n\n";

        // Secondary rectifier — diodes (default) or SR MOSFETs (sync).
        if (isSynchronous) {
            circuit << "* Synchronous rectifiers S_pos, S_neg (gated by primary PWM)\n";
            circuit << ".model SW2 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                    << " RON=0.01 ROFF=1Meg\n";
            // S_pos conducts when Q1 conducts (positive half).
            circuit << "Vsec_pos_sense diodePos sec_pos_in 0\n";
            circuit << "S_pos sec_pos_in out_node pwm1 0 SW2\n";
            circuit << "Vsec_neg_sense diodeNeg sec_neg_in 0\n";
            circuit << "S_neg sec_neg_in out_node pwm2 0 SW2\n";
        } else {
            circuit << "* CT-FW rectifier diodes D_pos, D_neg\n";
            circuit << "Vsec_pos_sense diodePos sec_pos_in 0\n";
            circuit << "D_pos sec_pos_in out_node DIDEAL\n";
            circuit << "Vsec_neg_sense diodeNeg sec_neg_in 0\n";
            circuit << "D_neg sec_neg_in out_node DIDEAL\n";
        }
        circuit << "Rsnub_dp diodePos out_node " << cfg.snubR << "\n";
        circuit << "Csnub_dp diodePos sn_dp " << std::scientific << cfg.snubC << std::fixed << "\n";
        circuit << "Rsnub_dpb sn_dp out_node 0.001\n";
        circuit << "Rsnub_dn diodeNeg out_node " << cfg.snubR << "\n";
        circuit << "Csnub_dn diodeNeg sn_dn " << std::scientific << cfg.snubC << std::fixed << "\n";
        circuit << "Rsnub_dnb sn_dn out_node 0.001\n\n";

        // Output cap and load
        circuit << "* Output cap (with ESR) and resistive load\n";
        circuit << "Cout out_node co_esr " << std::scientific << outputCapacitance << std::fixed
                << " ic=" << outputVoltage << "\n";
        circuit << "Rco_esr co_esr 0 " << esrCo << "\n";
        circuit << "Vout_sense out_node out_load 0\n";
        circuit << "Rload out_load 0 " << loadResistance << "\n\n";

        // Probes
        circuit << "* Probes\n";
        circuit << "Bvpri_diff vpri_diff 0 V=V(vab)\n\n";

        circuit << "* Transient analysis (UIC enables .ic)\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime
                << std::fixed << " UIC\n\n";

        circuit << "* Output signals\n";
        circuit << ".save v(vpri_diff) v(vin_dc) v(vab) v(out_node) "
                << "i(Vin_sense) i(Vout_sense)";
        if (isBridge) {
            circuit << " i(Vpri_sense)";
        } else {
            circuit << " i(Vpri_sense_a) i(Vpri_sense_b)";
        }
        circuit << "\n\n";

        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::fixed
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        circuit << ".options METHOD=" << cfg.method << " TRTOL=" << cfg.trTol << "\n";

        circuit << ".end\n";
        return circuit.str();
    }

    std::vector<OperatingPoint> Weinberg::simulate_and_extract_operating_points(double turnsRatio, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;

        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("Weinberg::simulate_and_extract_operating_points: ngspice is not available");
        }

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        const MAS::Variant variant = get_variant().value_or(MAS::Variant::CLASSIC);
        const bool isBridge = (variant == MAS::Variant::BRIDGE);
        const std::string priCurrentBranch = isBridge ? "vpri_sense#branch" : "vpri_sense_a#branch";

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto opPoint = get_operating_points()[opIndex];

                std::string netlist = generate_ngspice_circuit(turnsRatio, magnetizingInductance,
                                                               inputVoltageIndex, opIndex);
                double switchingFrequency = opPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);
                if (!simResult.success) {
                    throw std::runtime_error("Weinberg SPICE simulation failed: " + simResult.errorMessage);
                }

                NgspiceRunner::WaveformNameMapping waveformMapping = {
                    {{"voltage", "vpri_diff"}, {"current", priCurrentBranch}}
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

    std::vector<ConverterWaveforms> Weinberg::simulate_and_extract_topology_waveforms(double turnsRatio, double magnetizingInductance, size_t numberOfPeriods) {
        std::vector<ConverterWaveforms> results;

        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("Weinberg::simulate_and_extract_topology_waveforms: ngspice is not available");
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

                std::string netlist = generate_ngspice_circuit(turnsRatio, magnetizingInductance,
                                                               inputVoltageIndex, opIndex);
                double switchingFrequency = opPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = numberOfPeriods;
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);
                if (!simResult.success) {
                    throw std::runtime_error("Weinberg SPICE simulation failed: " + simResult.errorMessage);
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
                wf.get_mutable_output_voltages().push_back(getWaveform("out_node"));
                wf.get_mutable_output_currents().push_back(getWaveform("vout_sense#branch"));

                results.push_back(wf);
            }
        }

        set_num_periods_to_extract(originalNumPeriodsToExtract);
        return results;
    }

    // ============================================================
    //   get_extra_components_inputs — L1 (input coupled-L), Co
    // ============================================================

    std::vector<std::variant<Inputs, CAS::Inputs>> Weinberg::get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic)
    {
        if (mode == ExtraComponentsMode::REAL && !magnetic.has_value()) {
            throw std::invalid_argument("Weinberg::get_extra_components_inputs: mode REAL requires a designed magnetic");
        }

        if (lastSizedL1 <= 0.0 || lastSizedCo <= 0.0
            || extraL1VoltageWaveforms.empty()
            || extraCoVoltageWaveforms.empty())
        {
            throw std::runtime_error("Weinberg::get_extra_components_inputs: call process_operating_points() first");
        }

        const size_t nOps = extraL1VoltageWaveforms.size();
        auto checkSync = [nOps](const std::vector<Waveform>& wfms, const char* label) {
            if (wfms.size() != nOps) {
                throw std::runtime_error(std::string("Weinberg::get_extra_components_inputs: ") +
                    label + " waveform count out of sync");
            }
        };
        checkSync(extraL1CurrentWaveforms, "L1 current");
        checkSync(extraCoVoltageWaveforms, "Co voltage");
        checkSync(extraCoCurrentWaveforms, "Co current");

        const double fsw = get_operating_points()[0].get_switching_frequency();

        std::vector<std::variant<Inputs, CAS::Inputs>> result;

        // L1 input coupled inductor — emit as a 2-winding 1:1 magnetic
        // (constrained turnsRatio = 1).  Both windings carry the same DC
        // bias and the same triangular ripple (symmetric drive).
        {
            Inputs masInputs;
            DesignRequirements dr;
            DimensionWithTolerance L;
            L.set_nominal(lastSizedL1);
            dr.set_magnetizing_inductance(L);
            dr.set_name("inputCoupledInductor");
            dr.set_topology(Topologies::WEINBERG_CONVERTER);
            DimensionWithTolerance unityRatio;
            unityRatio.set_nominal(1.0);
            dr.set_turns_ratios(std::vector<DimensionWithTolerance>{unityRatio});
            dr.set_isolation_sides(std::vector<IsolationSide>{
                IsolationSide::PRIMARY, IsolationSide::PRIMARY});
            masInputs.set_design_requirements(dr);

            std::vector<OperatingPoint> masOps;
            for (size_t i = 0; i < extraL1VoltageWaveforms.size(); ++i) {
                OperatingPoint op;
                auto exA = complete_excitation(extraL1CurrentWaveforms[i],
                                               extraL1VoltageWaveforms[i], fsw, "L1a");
                auto exB = complete_excitation(extraL1CurrentWaveforms[i],
                                               extraL1VoltageWaveforms[i], fsw, "L1b");
                op.get_mutable_excitations_per_winding().push_back(exA);
                op.get_mutable_excitations_per_winding().push_back(exB);
                masOps.push_back(op);
            }
            masInputs.set_operating_points(masOps);
            result.emplace_back(std::move(masInputs));
        }

        // Co output filter capacitor.
        {
            CAS::Inputs casInputs;
            CAS::DesignRequirements dr;
            CAS::DimensionWithTolerance capacitance;
            capacitance.set_nominal(lastSizedCo);
            dr.set_capacitance(capacitance);
            double peakV = 0.0;
            for (const auto& wfm : extraCoVoltageWaveforms) {
                for (double v : wfm.get_data()) peakV = std::max(peakV, std::abs(v));
            }
            dr.set_rated_voltage(peakV * 1.2);
            dr.set_role(CAS::Application::OUTPUT_FILTER);
            dr.set_name("outputCapacitor");
            casInputs.set_design_requirements(dr);

            std::vector<CAS::TwoTerminalOperatingPoint> casOps;
            for (size_t i = 0; i < extraCoVoltageWaveforms.size(); ++i) {
                CAS::TwoTerminalOperatingPoint op;
                CAS::OperatingPointExcitation excitation;
                excitation.set_frequency(fsw);

                CAS::SignalDescriptor vSig;
                CAS::Waveform vWfm;
                vWfm.set_data(extraCoVoltageWaveforms[i].get_data());
                if (extraCoVoltageWaveforms[i].get_time())
                    vWfm.set_time(extraCoVoltageWaveforms[i].get_time().value());
                vSig.set_waveform(vWfm);
                excitation.set_voltage(vSig);

                CAS::SignalDescriptor iSig;
                CAS::Waveform iWfm;
                iWfm.set_data(extraCoCurrentWaveforms[i].get_data());
                if (extraCoCurrentWaveforms[i].get_time())
                    iWfm.set_time(extraCoCurrentWaveforms[i].get_time().value());
                iSig.set_waveform(iWfm);
                excitation.set_current(iSig);

                op.set_excitation(excitation);
                casOps.push_back(op);
            }
            casInputs.set_operating_points(casOps);
            result.emplace_back(std::move(casInputs));
        }

        return result;
    }

} // namespace OpenMagnetics
