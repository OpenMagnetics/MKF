#include "converter_models/Sepic.h"
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

    double Sepic::calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency, double maximumDutyCycle) {
        // SEPIC CCM, ideal:    M(D) = +D / (1 - D)
        // With η and Vd:       D = (Vo + Vd) / (Vin·η + Vo + Vd)        [TI SNVA168E]
        if (inputVoltage <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Sepic::calculate_duty_cycle: input voltage must be > 0");
        }
        if (outputVoltage <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Sepic::calculate_duty_cycle: Vo must be > 0");
        }
        double effVin = inputVoltage * efficiency;
        double dutyCycle = (outputVoltage + diodeVoltageDrop) / (effVin + outputVoltage + diodeVoltageDrop);
        // Explicit configurable maxD gate (no silent clamp). 1 % tolerance
        // for design-requirements rounding.
        constexpr double dutyTolerance = 0.01;
        if (dutyCycle > maximumDutyCycle * (1.0 + dutyTolerance)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Sepic::calculate_duty_cycle: duty cycle " + std::to_string(dutyCycle) +
                " exceeds maximumDutyCycle " + std::to_string(maximumDutyCycle) +
                " — converter would lose regulation; reduce Vo, raise Vin, or relax maximumDutyCycle.");
        }
        return dutyCycle;
    }

    double Sepic::calculate_conversion_ratio(double dutyCycle) {
        return dutyCycle / (1.0 - dutyCycle);
    }

    double Sepic::calculate_l1_min(double inputVoltage, double dutyCycle, double deltaIL1, double switchingFrequency) {
        if (deltaIL1 <= 0 || switchingFrequency <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Sepic::calculate_l1_min: deltaIL1 and fsw must be > 0");
        }
        return inputVoltage * dutyCycle / (deltaIL1 * switchingFrequency);
    }

    double Sepic::calculate_l2_min(double inputVoltage, double dutyCycle, double deltaIL2, double switchingFrequency) {
        // ΔIL2 = Vin·D / (L2·fsw)  (note: Vin, not Vo — both inductors see ±Vin during ON in SEPIC).
        if (deltaIL2 <= 0 || switchingFrequency <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Sepic::calculate_l2_min: deltaIL2 and fsw must be > 0");
        }
        return inputVoltage * dutyCycle / (deltaIL2 * switchingFrequency);
    }

    double Sepic::calculate_cs_min(double outputCurrent, double dutyCycle, double deltaVCs, double switchingFrequency) {
        // Cs ≥ Iout · D / (ΔVCs · fsw)        [TI SNVA168E]
        if (deltaVCs <= 0 || switchingFrequency <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Sepic::calculate_cs_min: deltaVCs and fsw must be > 0");
        }
        return outputCurrent * dutyCycle / (deltaVCs * switchingFrequency);
    }

    double Sepic::calculate_dcm_K(double L1, double L2, double switchingFrequency, double loadResistance) {
        if (L1 + L2 <= 0 || loadResistance <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Sepic::calculate_dcm_K: bad L1/L2 or R");
        }
        double Le = (L1 * L2) / (L1 + L2);
        return 2.0 * Le * switchingFrequency / loadResistance;
    }

    // ============================================================
    //   Constructors
    // ============================================================

    Sepic::Sepic(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    AdvancedSepic::AdvancedSepic(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    // ============================================================
    //   Per-operating-point analytical worker
    // ============================================================

    OperatingPoint Sepic::process_operating_points_for_input_voltage(double inputVoltage, const TopologyExcitation& outputOperatingPoint, double inductanceL1) {
        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double outputVoltage      = outputOperatingPoint.get_output_voltages()[0];
        double outputCurrent      = outputOperatingPoint.get_output_currents()[0];
        double diodeVoltageDrop   = get_diode_voltage_drop();
        double efficiency = require_input(get_efficiency(), "Sepic", "efficiency");

        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, diodeVoltageDrop, efficiency, maximumDutyCycle.value_or(0.95));

        // Common derived quantities (V1 / V2). Power balance:
        //   Vin·η·IL1avg = Vo·Iout    ⇒   IL1avg = Vo·Iout / (Vin·η).
        //   Lossless:  IL1avg = Iout·D/(1-D).
        double IL2avg     = outputCurrent;
        double IL1avg     = outputCurrent * dutyCycle / ((1.0 - dutyCycle) * efficiency);
        double VCs        = inputVoltage;                              // steady state
        double deltaIL1   = inductanceL1 > 0 ? (inputVoltage * dutyCycle) / (inductanceL1 * switchingFrequency)
                                             : 0.0;
        double deltaIL2_target = std::max(l2RipplePct * IL2avg, 1e-6);
        double inductanceL2    = calculate_l2_min(inputVoltage, dutyCycle, deltaIL2_target, switchingFrequency);
        double deltaIL2        = (inputVoltage * dutyCycle) / (inductanceL2 * switchingFrequency);

        double deltaVCs_target = std::max(csRipplePct * VCs, 1e-3);
        double couplingCap     = calculate_cs_min(outputCurrent, dutyCycle, deltaVCs_target, switchingFrequency);

        // ---- Common diagnostics ----
        lastDutyCycle = dutyCycle;
        lastConversionRatio = dutyCycle / (1.0 - dutyCycle);
        lastInputInductorAverage = IL1avg;
        lastOutputInductorAverage = IL2avg;
        lastInputInductorRipple = deltaIL1;
        lastOutputInductorRipple = deltaIL2;
        lastSwitchPeakVoltage = inputVoltage + outputVoltage;          // VS,peak = Vin + Vo
        double IS_peak = IL1avg + IL2avg + (deltaIL1 + deltaIL2) / 2.0;
        lastSwitchPeakCurrent = IS_peak;
        lastDiodePeakReverseVoltage = lastSwitchPeakVoltage;
        lastDiodePeakCurrent  = IS_peak;
        lastCouplingCapVoltage = VCs;
        lastCouplingCapRmsCurrent = (dutyCycle < 1.0 - 1e-9)
            ? outputCurrent * std::sqrt(dutyCycle / (1.0 - dutyCycle))
            : 0.0;
        double loadResistance = (outputCurrent > 0) ? outputVoltage / outputCurrent : 0.0;
        lastDcmK = calculate_dcm_K(inductanceL1, inductanceL2, switchingFrequency, loadResistance);
        lastDcmKcrit = std::pow(1.0 - dutyCycle, 2);
        lastIsCcm = (lastDcmK > lastDcmKcrit);
        lastSizedL2 = inductanceL2;
        lastSizedCs = couplingCap;

        // Output cap Co sized from ΔIL_load / (8·fsw·ΔVo). In SEPIC the load
        // sees only AC ripple of iD (rectified iL1+iL2) minus the steady DC
        // through Rload — to first order, ΔIo ≈ ΔIL2.
        double deltaVo_target = std::max(coRipplePct * outputVoltage, 1e-3);
        double outputCapacitance = deltaIL2 / (8.0 * switchingFrequency * deltaVo_target);
        outputCapacitance = std::max(outputCapacitance, 1e-6);
        lastSizedCo = outputCapacitance;
        double deltaVo_actual = deltaIL2 / (8.0 * switchingFrequency * outputCapacitance);

        // Primary inductor (L1) waveforms — TRIANGULAR around IL1avg.
        // Voltage swings between +Vin (ON) and -Vo (OFF).
        // peakToPeak = Vin - (-Vo) = Vin + Vo.  V·s mean = 0 by balance, offset=0.
        double primaryVoltagePeakToPeak = inputVoltage + outputVoltage;
        {
            Waveform currentWaveform = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaIL1, switchingFrequency, dutyCycle, IL1avg, 0);
            Waveform voltageWaveform = Inputs::create_waveform(
                WaveformLabel::RECTANGULAR, primaryVoltagePeakToPeak, switchingFrequency, dutyCycle, 0, 0);
            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // ---- Stash extra-component waveforms for get_extra_components_inputs ----
        // L2 voltage: ON value = -VCs = -Vin, OFF value = +Vo.  pp = Vin + Vo,
        // mean = 0 by V·s balance.  L2 is "high" (positive Vo) during OFF
        // sub-interval, so pass dutyCycle = (1-D).
        {
            Waveform vL2 = Inputs::create_waveform(
                WaveformLabel::RECTANGULAR, primaryVoltagePeakToPeak, switchingFrequency, 1.0 - dutyCycle, 0.0, 0);
            Waveform iL2 = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaIL2, switchingFrequency, 1.0 - dutyCycle, IL2avg, 0);
            extraL2VoltageWaveforms.push_back(vL2);
            extraL2CurrentWaveforms.push_back(iL2);
        }

        // Cs voltage: triangular around VCs = Vin with peak-to-peak ΔVCs.
        // Cs current: rectangular between -IL1avg (during ON, D·T — Cs charges
        // from L1 path; L2 holds steady through GND) and +IL2avg (during OFF —
        // Cs discharges as iL2 + iL1 both flow through it). Net pp = IL1+IL2.
        {
            double deltaVCs_actual = (couplingCap > 0)
                ? (outputCurrent * dutyCycle) / (couplingCap * switchingFrequency)
                : 0.0;
            Waveform vCs = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaVCs_actual, switchingFrequency, dutyCycle, VCs, 0);
            double iCspp = IL1avg + IL2avg;
            Waveform iCs = Inputs::create_waveform(
                WaveformLabel::RECTANGULAR, iCspp, switchingFrequency, 1.0 - dutyCycle, 0.0, 0);
            extraCsVoltageWaveforms.push_back(vCs);
            extraCsCurrentWaveforms.push_back(iCs);
        }

        // Co voltage: tiny triangular ripple ΔVo around the positive mean +Vo.
        // Co current: AC component of iD (= iL1+iL2 during OFF, 0 during ON).
        // To first order treat as iL2 ripple around 0 mean.
        {
            Waveform vCo = Inputs::create_waveform(
                WaveformLabel::TRIANGULAR, deltaVo_actual, switchingFrequency, 1.0 - dutyCycle, outputVoltage, 0);
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
        perOpSizedCs.push_back(lastSizedCs);
        perOpSizedCo.push_back(lastSizedCo);

        return operatingPoint;
    }

    bool Sepic::run_checks(bool assert) {
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

    DesignRequirements Sepic::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double efficiency = require_input(get_efficiency(), "Sepic", "efficiency");

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Sepic::process_design_requirements: missing both currentRippleRatio and maximumSwitchCurrent");
        }

        // Worst-case L1: largest L1 needed across all OPs at maximum Vin.
        double maximumDeltaIL1 = 0.0;
        if (get_current_ripple_ratio()) {
            double rippleRatio = get_current_ripple_ratio().value();
            for (const auto& op : get_operating_points()) {
                double Iout = op.get_output_currents()[0];
                double Vo   = op.get_output_voltages()[0];
                double D    = calculate_duty_cycle(maximumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, maximumDutyCycle.value_or(0.95));
                double IL1avg = Iout * D / (1.0 - D);
                maximumDeltaIL1 = std::max(maximumDeltaIL1, rippleRatio * IL1avg);
            }
        }
        if (get_maximum_switch_current()) {
            double IsMax = get_maximum_switch_current().value();
            for (const auto& op : get_operating_points()) {
                double Iout = op.get_output_currents()[0];
                double Vo   = op.get_output_voltages()[0];
                double D    = calculate_duty_cycle(minimumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, maximumDutyCycle.value_or(0.95));
                double IL1avg = Iout * D / (1.0 - D);
                double IL2avg = Iout;
                // IS_peak ≈ IL1+IL2 + (ΔIL1+ΔIL2)/2; treat ΔIL2 as ~30 % of IL2avg
                // (the L2 sizing default), so the L1 ripple budget is:
                double residual = IsMax - (IL1avg + IL2avg) - 0.5 * 0.30 * IL2avg;
                double deltaIL1 = std::max(2.0 * residual, 0.0);
                maximumDeltaIL1 = std::max(maximumDeltaIL1, deltaIL1);
            }
        }
        if (maximumDeltaIL1 <= 0) {
            throw std::invalid_argument("Sepic::process_design_requirements: derived ΔIL1 budget is non-positive");
        }

        double maximumNeededInductance = 0.0;
        for (const auto& op : get_operating_points()) {
            double switchingFrequency = op.get_switching_frequency();
            double Vo = op.get_output_voltages()[0];
            double D  = calculate_duty_cycle(maximumInputVoltage, Vo, get_diode_voltage_drop(), efficiency, maximumDutyCycle.value_or(0.95));
            double L1 = calculate_l1_min(maximumInputVoltage, D, maximumDeltaIL1, switchingFrequency);
            maximumNeededInductance = std::max(maximumNeededInductance, L1);
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        const bool isCoupled = get_coupled_inductor().value_or(false);
        std::vector<IsolationSide> isolationSides;
        if (isCoupled) {
            // V2 coupled-inductor: two windings on the same core, no galvanic
            // isolation. Both share isolationSide=PRIMARY (no isolation).
            // 1:1 turns ratio is the canonical zero-ripple choice (k·√(L2/L1)=1
            // → for L1=L2, k=1).
            isolationSides.push_back(get_isolation_side_from_index(0));
            isolationSides.push_back(get_isolation_side_from_index(0));
            DimensionWithTolerance unityRatio;
            unityRatio.set_nominal(1.0);
            designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{unityRatio});
        } else {
            isolationSides.push_back(get_isolation_side_from_index(0));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::SEPIC_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Sepic::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        (void) turnsRatios;  // SEPIC has no transformer; turns ratios are unused.

        // Reset extra-component waveform stash so a re-run does not accumulate.
        extraL2VoltageWaveforms.clear();
        extraL2CurrentWaveforms.clear();
        extraCsVoltageWaveforms.clear();
        extraCsCurrentWaveforms.clear();
        extraCoVoltageWaveforms.clear();
        extraCoCurrentWaveforms.clear();

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
        perOpSizedCs.clear();
        perOpSizedCo.clear();

        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

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

    std::vector<OperatingPoint> Sepic::process_operating_points(OpenMagnetics::Magnetic magnetic) {
        run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(
            magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();

        std::vector<double> turnsRatios = magnetic.get_turns_ratios();

        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedSepic::process() {
        Sepic::run_checks(_assertErrors);

        extraL2VoltageWaveforms.clear();
        extraL2CurrentWaveforms.clear();
        extraCsVoltageWaveforms.clear();
        extraCsCurrentWaveforms.clear();
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
        designRequirements.set_topology(Topologies::SEPIC_CONVERTER);
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
    //   SPICE netlist
    // ============================================================
    //
    // SEPIC node convention:
    //   vin_dc → Vin_sense → l1_in →[L1]→ node_A
    //   node_A →[Vc_sense]→ node_C →[Cs]→ node_B
    //   node_A → S1 → 0  (PWM low-side switch)
    //   node_B →[L2]→ 0   (L2 to GND, NOT to vout — SEPIC distinguishing feature)
    //   node_B → D1 → vout  (anode at node_B, cathode at vout; vout > 0)
    //   vout → Cout → 0;  vout → Rload → 0
    // ============================================================

    std::string Sepic::generate_ngspice_circuit(double inductanceL1, size_t inputVoltageIndex, size_t operatingPointIndex) {
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames_;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames_);

        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("Sepic::generate_ngspice_circuit: inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("Sepic::generate_ngspice_circuit: operatingPointIndex out of range");
        }

        double inputVoltage = inputVoltages[inputVoltageIndex];
        auto opPoint = get_operating_points()[operatingPointIndex];

        double outputVoltage      = opPoint.get_output_voltages()[0];
        double outputCurrent      = opPoint.get_output_currents()[0];
        double switchingFrequency = opPoint.get_switching_frequency();
        double diodeVoltageDrop   = get_diode_voltage_drop();
        double efficiency = require_input(get_efficiency(), "Sepic", "efficiency");

        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, diodeVoltageDrop, efficiency, maximumDutyCycle.value_or(0.95));

        double IL1avg = outputCurrent * dutyCycle / (1.0 - dutyCycle);
        double IL2avg = outputCurrent;
        double VCs    = inputVoltage;
        double deltaIL2_target = std::max(l2RipplePct * IL2avg, 1e-6);
        double inductanceL2    = calculate_l2_min(inputVoltage, dutyCycle, deltaIL2_target, switchingFrequency);
        double deltaVCs_target = std::max(csRipplePct * VCs, 1e-3);
        double couplingCap     = calculate_cs_min(outputCurrent, dutyCycle, deltaVCs_target, switchingFrequency);
        double deltaIL2 = (inputVoltage * dutyCycle) / (inductanceL2 * switchingFrequency);
        double deltaVo  = std::max(coRipplePct * outputVoltage, 1e-3);
        double outputCapacitance = deltaIL2 / (8.0 * switchingFrequency * deltaVo);
        outputCapacitance = std::max(outputCapacitance, 1e-6);

        if (outputCurrent <= 0.0) {
            throw std::invalid_argument(
                "Sepic::generate_ngspice_circuit: outputCurrent must be > 0 "
                "(received " + std::to_string(outputCurrent) +
                "); cannot derive load resistance R_load = Vo / Iout.");
        }
        double loadResistance = outputVoltage / outputCurrent;

        const SpiceSimulationConfig cfg = spice_config();

        const bool isSynchronous = get_synchronous_rectifier().value_or(false);
        const bool isCoupled     = get_coupled_inductor().value_or(false);

        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;

        int periodsToExtract = get_num_periods_to_extract();
        int settlingPeriods  = get_num_steady_state_periods();
        const int numPeriodsTotal = settlingPeriods + periodsToExtract;
        double simTime  = numPeriodsTotal * period;
        double startTime = settlingPeriods * period;
        double samplesPerPeriod = (switchingFrequency >= 1e6) ? 500.0 : cfg.samplesPerPeriod;
        double stepTime = period / samplesPerPeriod;

        // Ideal reference: zero winding DCR + cap ESR (consistent with the other converters' ideal
        // decks, which use pure L/C). Real parasitics belong to a designed magnetic / real cap, not
        // to the "ideal" topology reference. (Were 50 mOhm / 5 mOhm hardcoded.)
        const double dcrL1 = 0.0;
        const double dcrL2 = 0.0;
        const double esrCs = 0.0;
        const double esrCo = 0.0;

        circuit << "* SEPIC Converter ("
                << (isCoupled ? "V2 coupled-inductor" : "V1 uncoupled")
                << (isSynchronous ? "/synchronous" : "")
                << ") - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vo=" << outputVoltage << "V, "
                << "f=" << (switchingFrequency / 1e3) << "kHz, D=" << (dutyCycle * 100) << " pct\n";
        circuit << "* L1=" << (inductanceL1 * 1e6) << "uH, L2=" << (inductanceL2 * 1e6) << "uH, "
                << "Cs=" << (couplingCap * 1e6) << "uF, Co=" << (outputCapacitance * 1e6) << "uF, "
                << "Iout=" << outputCurrent << "A\n\n";

        // DC Input + sense
        circuit << "* DC Input + sense\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n";
        circuit << "Vin_sense vin_dc l1_in 0\n\n";

        // L1 with DCR
        circuit << "* L1 input inductor (with DCR)\n";
        circuit << "L1 l1_in l1_dcr_mid " << std::scientific << inductanceL1 << std::fixed
                << " ic=" << IL1avg << "\n";
        circuit << "Rdcr_l1 l1_dcr_mid node_A " << dcrL1 << "\n";
        circuit << "Vl1_sense node_A node_A_int 0\n\n";

        // S1 PWM switch (from node_A_int to GND)
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

        // Cs coupling cap (node_A_int → node_B with sense + ESR)
        circuit << "* Cs coupling cap (with ESR; +plate at node_A side, IC=Vin)\n";
        circuit << "Vc_sense node_A_int node_C 0\n";
        circuit << "Cs node_C node_C_esr " << std::scientific << couplingCap << std::fixed
                << " ic=" << VCs << "\n";
        circuit << "Rcs_esr node_C_esr node_B " << esrCs << "\n\n";

        // L2 — bottom terminal at GND, top at node_B (SEPIC distinguishing).
        circuit << "* L2 inductor (bottom-side at GND; with DCR)\n";
        circuit << "L2 0 l2_top " << std::scientific << inductanceL2 << std::fixed
                << " ic=" << IL2avg << "\n";
        circuit << "Rdcr_l2 l2_top node_B " << dcrL2 << "\n";
        circuit << "Vl2_sense node_B node_B_int 0\n\n";

        // Rectifier between node_B_int and vout. D1 (default) or S2 (synchronous).
        circuit << "Vd_sense d_anode 0 0\n";  // unused for SEPIC routing — kept for parity
        if (isSynchronous) {
            circuit << "* S2 synchronous rectifier (complementary PWM, no dead-time)\n";
            circuit << "Vpwm_inv pwm_ctrl_inv 0 PULSE(" << cfg.pwmHigh << " 0 0 "
                    << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                    << tOn << " " << period << std::fixed << ")\n";
            circuit << ".model SW2 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                    << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
            circuit << "Vrect_sense node_B_int rect_in 0\n";
            circuit << "S2 rect_in vout pwm_ctrl_inv 0 SW2\n";
            circuit << "Rsnub_d1 rect_in 0 " << cfg.snubR << "\n"
                    << "Csnub_d1 rect_in snub_d1_int " << std::scientific << cfg.snubC << std::fixed << "\n"
                    << "Rsnub_d1b snub_d1_int 0 " << require_spice_field(cfg.snubDampR, "snubDampR") << "\n\n";
        } else {
            circuit << "* D1 rectifier diode\n";
            circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                    << " RS=" << cfg.diodeRS << std::fixed << ")\n";
            circuit << "Vrect_sense node_B_int rect_in 0\n";
            circuit << "D1 rect_in vout DIDEAL\n";
            circuit << "Rsnub_d1 rect_in 0 " << cfg.snubR << "\n"
                    << "Csnub_d1 rect_in snub_d1_int " << std::scientific << cfg.snubC << std::fixed << "\n"
                    << "Rsnub_d1b snub_d1_int 0 " << require_spice_field(cfg.snubDampR, "snubDampR") << "\n\n";
        }

        // V2 coupled-inductor: K_L1L2 between L1 and L2.
        if (isCoupled) {
            double k = get_coupling_coefficient().value_or(0.999);
            if (k <= 0.0 || k >= 1.0) {
                throw std::invalid_argument(
                    "Sepic::generate_ngspice_circuit: couplingCoefficient must be "
                    "in (0, 1); received " + std::to_string(k));
            }
            circuit << "* V2 coupled-inductor: K_L1L2 couples L1 and L2 on the same core\n";
            circuit << "K_L1L2 L1 L2 " << k << "\n\n";
        }

        // Output cap (positive Vo) and resistive load.
        circuit << "* Output cap (with ESR) and resistive load\n";
        circuit << "Cout vout co_esr " << std::scientific << outputCapacitance << std::fixed
                << " ic=" << outputVoltage << "\n";
        circuit << "Rco_esr co_esr 0 " << esrCo << "\n";
        circuit << "Vout_sense vout vout_load 0\n";
        circuit << "Rload vout_load 0 " << loadResistance << "\n\n";

        // Switch-node + L1 winding voltage probes
        circuit << "* Bridge / switch-node probe\n";
        circuit << "Bsw sw_probe 0 V=V(node_A)\n";
        circuit << "Bvpri_diff vpri_diff 0 V=V(l1_in)-V(node_A)\n\n";

        // Transient analysis with UIC (initial conditions)
        circuit << "* Transient Analysis (UIC enables .ic above)\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime
                << std::fixed << " UIC\n\n";

        circuit << "* Output signals\n";
        circuit << ".save v(vpri_diff) v(vin_dc) v(node_A) v(node_B) v(sw_probe) v(vout) "
                << "i(Vin_sense) i(Vl1_sense) i(Vc_sense) i(Vl2_sense) i(Vrect_sense) i(Vout_sense)\n\n";

        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::fixed
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        circuit << ".options METHOD=" << cfg.method << " TRTOL=" << cfg.trTol << "\n";

        circuit << ".end\n";
        return circuit.str();
    }

    std::vector<OperatingPoint> Sepic::simulate_and_extract_operating_points(double inductanceL1) {
        std::vector<OperatingPoint> operatingPoints;

        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("Sepic::simulate_and_extract_operating_points: ngspice is not available");
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
                    throw std::runtime_error("Sepic SPICE simulation failed: " + simResult.errorMessage);
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

    std::vector<ConverterWaveforms> Sepic::simulate_and_extract_topology_waveforms(double inductanceL1, size_t numberOfPeriods) {
        std::vector<ConverterWaveforms> results;

        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("Sepic::simulate_and_extract_topology_waveforms: ngspice is not available");
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
                    throw std::runtime_error("Sepic SPICE simulation failed: " + simResult.errorMessage);
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
                wf.get_mutable_output_voltages().push_back(getWaveform("vout"));
                wf.get_mutable_output_currents().push_back(getWaveform("vout_sense#branch"));

                results.push_back(wf);
            }
        }

        set_num_periods_to_extract(originalNumPeriodsToExtract);
        return results;
    }

    // ============================================================
    //   get_extra_components_inputs — L2 / Cs / Co
    // ============================================================

    std::vector<std::variant<Inputs, CAS::Inputs>> Sepic::get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic)
    {
        if (mode == ExtraComponentsMode::REAL && !magnetic.has_value()) {
            throw std::invalid_argument("Sepic::get_extra_components_inputs: mode REAL requires a designed magnetic");
        }

        if (lastSizedL2 <= 0.0 || lastSizedCs <= 0.0 || lastSizedCo <= 0.0
            || extraL2VoltageWaveforms.empty()
            || extraCsVoltageWaveforms.empty()
            || extraCoVoltageWaveforms.empty())
        {
            throw std::runtime_error("Sepic::get_extra_components_inputs: call process_operating_points() first");
        }

        const size_t nOps = extraL2VoltageWaveforms.size();
        auto checkSync = [nOps](const std::vector<Waveform>& wfms, const char* label) {
            if (wfms.size() != nOps) {
                throw std::runtime_error(std::string("Sepic::get_extra_components_inputs: ") +
                    label + " waveform count out of sync");
            }
        };
        checkSync(extraL2CurrentWaveforms, "L2 current");
        checkSync(extraCsVoltageWaveforms, "Cs voltage");
        checkSync(extraCsCurrentWaveforms, "Cs current");
        checkSync(extraCoVoltageWaveforms, "Co voltage");
        checkSync(extraCoCurrentWaveforms, "Co current");

        const double fsw = get_operating_points()[0].get_switching_frequency();

        std::vector<std::variant<Inputs, CAS::Inputs>> result;

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
            dr.set_topology(Topologies::SEPIC_CONVERTER);
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
            dr.set_rated_voltage(peakV * 1.2);
            dr.set_role(role);
            dr.set_name(name);
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

        // V1/V2: L2, Cs, Co — L1 IS the primary magnetic.
        buildInductor("outputInductor", lastSizedL2, extraL2VoltageWaveforms, extraL2CurrentWaveforms);
        buildCapacitor("couplingCapacitor", lastSizedCs,
                       extraCsVoltageWaveforms, extraCsCurrentWaveforms, CAS::Application::DC_LINK);
        buildCapacitor("outputCapacitor",   lastSizedCo,
                       extraCoVoltageWaveforms, extraCoCurrentWaveforms, CAS::Application::OUTPUT_FILTER);

        return result;
    }

} // namespace OpenMagnetics
