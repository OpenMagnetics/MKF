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

    double Cuk::calculate_duty_cycle(double inputVoltage, double outputVoltageMagnitude, double diodeVoltageDrop, double efficiency) {
        // Cuk CCM, ideal:                M(D)  = -D/(1-D)
        // With diode drop and efficiency: |Vo| = (Vin·η - Vd·(1-D))·D/(1-D)
        // Ignoring the small Vd term (it's a 2nd-order correction in the
        // already-approximated η model), solve |Vo| = Vin·η·D/(1-D) for D:
        //                                  D = |Vo| / (|Vo| + Vin·η)
        // Then add a small Vd-driven correction: D += Vd / (Vin·η + Vd) · (1-D)
        // For the test fixtures (Vd ≤ 0.6 V, Vin ≥ 3 V) the correction is
        // < 5 %; ignore for analytical anchor and let SPICE see the real
        // diode drop.
        if (inputVoltage <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_duty_cycle: input voltage must be > 0");
        }
        if (outputVoltageMagnitude <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cuk::calculate_duty_cycle: |Vo| must be > 0");
        }
        double effVin = inputVoltage * efficiency;
        double dutyCycle = (outputVoltageMagnitude + diodeVoltageDrop) / (outputVoltageMagnitude + diodeVoltageDrop + effVin);
        if (dutyCycle >= 0.95) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Cuk::calculate_duty_cycle: duty cycle " + std::to_string(dutyCycle) +
                " >= 0.95 — converter would lose regulation; reduce |Vo| or raise Vin");
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

    // ============================================================
    //   Constructors
    // ============================================================

    Cuk::Cuk(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    AdvancedCuk::AdvancedCuk(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    // ============================================================
    //   Per-operating-point analytical worker
    // ============================================================

    OperatingPoint Cuk::process_operating_points_for_input_voltage(double inputVoltage, const BaseOperatingPoint& outputOperatingPoint, double inductanceL1) {
        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double outputVoltageMag   = std::abs(outputOperatingPoint.get_output_voltages()[0]);
        double outputCurrent      = outputOperatingPoint.get_output_currents()[0];
        double diodeVoltageDrop   = get_diode_voltage_drop();
        double efficiency = 1.0;
        if (get_efficiency()) efficiency = get_efficiency().value();

        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltageMag, diodeVoltageDrop, efficiency);

        // Internally-sized L2, C1 (V1 only — see CUK_PLAN.md §13).
        double IL2avg     = outputCurrent;
        double IL1avg     = outputCurrent * dutyCycle / (1.0 - dutyCycle);
        double deltaIL1   = inductanceL1 > 0 ? (inputVoltage * dutyCycle) / (inductanceL1 * switchingFrequency)
                                             : 0.0;
        double deltaIL2_target = std::max(l2RipplePct * IL2avg, 1e-6);
        double inductanceL2    = calculate_l2_min(outputVoltageMag, dutyCycle, deltaIL2_target, switchingFrequency);
        double deltaIL2        = (outputVoltageMag * (1.0 - dutyCycle)) / (inductanceL2 * switchingFrequency);

        double VC1            = calculate_coupling_cap_voltage(inputVoltage, dutyCycle);
        double deltaVC1_target = std::max(c1RipplePct * VC1, 1e-3);
        double couplingCap    = calculate_c1_min(outputCurrent, dutyCycle, deltaVC1_target, switchingFrequency);

        // Primary inductor (L1) waveforms — TRIANGULAR around IL1avg.
        // Voltage swings between +Vin (ON) and Vin - VC1 (OFF, negative).
        double primaryVoltageMaximum = inputVoltage;             // ON
        double primaryVoltageMinimum = inputVoltage - VC1;       // OFF (negative)
        double primaryVoltagePeakToPeak = primaryVoltageMaximum - primaryVoltageMinimum;  // = VC1
        double primaryCurrentMinimum = IL1avg - deltaIL1 / 2.0;

        // ---- Diagnostics ----
        lastDutyCycle = dutyCycle;
        lastConversionRatio = calculate_conversion_ratio(dutyCycle);
        lastCouplingCapVoltage = VC1;
        lastInputInductorAverage = IL1avg;
        lastOutputInductorAverage = IL2avg;
        lastInputInductorRipple = deltaIL1;
        lastOutputInductorRipple = deltaIL2;
        lastSwitchPeakVoltage = inputVoltage + outputVoltageMag;
        lastDiodePeakReverseVoltage = lastSwitchPeakVoltage;
        double IS_peak = IL1avg + IL2avg + (deltaIL1 + deltaIL2) / 2.0;
        lastSwitchPeakCurrent = IS_peak;
        lastDiodePeakCurrent = IS_peak;
        lastCouplingCapRmsCurrent = std::sqrt(dutyCycle * (1.0 - dutyCycle)) * (IL1avg + IL2avg);
        double loadResistance = (outputCurrent > 0) ? outputVoltageMag / outputCurrent : 0.0;
        lastDcmK = calculate_dcm_K(inductanceL1, inductanceL2, switchingFrequency, loadResistance);
        lastDcmKcrit = std::pow(1.0 - dutyCycle, 2);
        lastIsCcm = (lastDcmK > lastDcmKcrit);
        lastRhpZeroFrequency = calculate_rhp_zero_frequency(loadResistance, dutyCycle, inductanceL2);
        lastRecommendedLoopBandwidth = lastRhpZeroFrequency / 5.0;
        lastSizedL2 = inductanceL2;
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

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

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
        double efficiency = 1.0;
        if (get_efficiency()) efficiency = get_efficiency().value();

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
                double D    = calculate_duty_cycle(maximumInputVoltage, Vo, get_diode_voltage_drop(), efficiency);
                double IL1avg = Iout * D / (1.0 - D);
                maximumDeltaIL1 = std::max(maximumDeltaIL1, rippleRatio * IL1avg);
            }
        }
        if (get_maximum_switch_current()) {
            double IsMax = get_maximum_switch_current().value();
            for (const auto& op : get_operating_points()) {
                double Iout = op.get_output_currents()[0];
                double Vo   = std::abs(op.get_output_voltages()[0]);
                double D    = calculate_duty_cycle(minimumInputVoltage, Vo, get_diode_voltage_drop(), efficiency);
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
            double D  = calculate_duty_cycle(maximumInputVoltage, Vo, get_diode_voltage_drop(), efficiency);
            double L1 = calculate_l1_min(maximumInputVoltage, D, maximumDeltaIL1, switchingFrequency);
            maximumNeededInductance = std::max(maximumNeededInductance, L1);
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::CUK_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Cuk::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        (void) turnsRatios;  // Cuk V1 has no transformer; turns ratios are unused.
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
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
        designRequirements.set_topology(Topologies::CUK_CONVERTER);
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
        double efficiency = 1.0;
        if (get_efficiency()) efficiency = get_efficiency().value();

        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltageMag, diodeVoltageDrop, efficiency);

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
        const double dcrL1 = 50e-3;          // 50 mΩ DCR
        const double dcrL2 = 50e-3;
        const double esrC1 = 5e-3;           // 5 mΩ ESR
        const double esrCo = 5e-3;

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

        circuit << "* Cuk Converter (non-isolated, V1) - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, |Vo|=" << outputVoltageMag << "V (sign: negative), "
                << "f=" << (switchingFrequency / 1e3) << "kHz, D=" << (dutyCycle * 100) << " pct\n";
        circuit << "* L1=" << (inductanceL1 * 1e6) << "uH, L2=" << (inductanceL2 * 1e6) << "uH, "
                << "C1=" << (couplingCap * 1e6) << "uF, Co=" << (outputCapacitance * 1e6) << "uF, "
                << "Iout=" << outputCurrent << "A\n\n";

        // DC Input + sense
        circuit << "* DC Input + sense\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n";
        circuit << "Vin_sense vin_dc l1_in 0\n\n";

        // L1 with DCR (split as L + R for bake-in; use named series RDCR_l1)
        circuit << "* L1 input inductor (with DCR)\n";
        circuit << "L1 l1_in l1_dcr_mid " << std::scientific << inductanceL1 << std::fixed
                << " ic=" << IL1avg << "\n";
        circuit << "Rdcr_l1 l1_dcr_mid node_A " << dcrL1 << "\n";
        circuit << "Vl1_sense node_A node_A_int 0\n\n";

        // S1 PWM switch (from node_A to GND)
        circuit << "* PWM low-side switch on node_A\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                << tOn << " " << period << std::fixed << ")\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH << "\n";
        circuit << "S1 node_A_int 0 pwm_ctrl 0 SW1\n";
        circuit << "Rsnub_s1 node_A_int 0 " << cfg.snubR << "\n"
                << "Csnub_s1 node_A_int snub_s1_int " << std::scientific << cfg.snubC << std::fixed << "\n"
                << "Rsnub_s1b snub_s1_int 0 0.001\n\n";

        // C1 coupling cap (from node_A_int through ESR to node_B). In the
        // ON sub-interval VC1 polarity is +(node_A) - (node_B) ≈ Vin/(1-D).
        circuit << "* C1 coupling cap (with ESR)\n";
        circuit << "Vc1_sense node_A_int node_C 0\n";
        circuit << "C1 node_C node_C_esr " << std::scientific << couplingCap << std::fixed
                << " ic=" << VC1 << "\n";
        circuit << "Rc1_esr node_C_esr node_B " << esrC1 << "\n\n";

        // D1 freewheel diode (anode at node_B, cathode at GND via sense)
        circuit << "* D1 freewheel diode\n";
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::fixed << ")\n";
        circuit << "Vd_sense d_cath 0 0\n";
        circuit << "D1 node_B d_cath DIDEAL\n";
        circuit << "Rsnub_d1 node_B 0 " << cfg.snubR << "\n"
                << "Csnub_d1 node_B snub_d1_int " << std::scientific << cfg.snubC << std::fixed << "\n"
                << "Rsnub_d1b snub_d1_int 0 0.001\n\n";

        // L2 output inductor (from node_B to vout_load_node, with DCR)
        circuit << "* L2 output inductor (with DCR)\n";
        circuit << "L2 node_B l2_dcr_mid " << std::scientific << inductanceL2 << std::fixed
                << " ic=" << IL2avg << "\n";
        circuit << "Rdcr_l2 l2_dcr_mid out_sense " << dcrL2 << "\n";
        circuit << "Vl2_sense out_sense vout_load_node 0\n\n";

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

} // namespace OpenMagnetics
