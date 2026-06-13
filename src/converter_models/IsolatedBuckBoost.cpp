#include "converter_models/IsolatedBuckBoost.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/CircuitSimulatorInterface.h"
#include "support/Utils.h"
#include <cfloat>
#include <limits>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double IsolatedBuckBoost::calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency) {
        auto dutyCycle = outputVoltage / (inputVoltage + outputVoltage) * efficiency;
        if (dutyCycle >= 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Duty cycle must be smaller than 1");
        }
        // Explicit maxD gate (no silent clamp). 1 % tolerance for
        // design-requirements rounding. Same pattern as Flyback,
        // forward family, Buck, Boost, IsolatedBuck.
        if (maximumDutyCycle.has_value()) {
            const double dMax = maximumDutyCycle.value();
            constexpr double dutyTolerance = 0.01;
            if (dutyCycle > dMax * (1.0 + dutyTolerance)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "IsolatedBuckBoost: required dutyCycle " + std::to_string(dutyCycle) +
                    " exceeds maximumDutyCycle " + std::to_string(dMax) +
                    " at Vin=" + std::to_string(inputVoltage) +
                    ", Vout=" + std::to_string(outputVoltage) +
                    ". Raise Vin_min, lower Vout, or relax maximumDutyCycle.");
            }
        }
        return dutyCycle;
    }

    IsolatedBuckBoost::IsolatedBuckBoost(const json& j) {
        from_json(j, *this);
    }

    AdvancedIsolatedBuckBoost::AdvancedIsolatedBuckBoost(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint IsolatedBuckBoost::processOperatingPointsForInputVoltage(double inputVoltage, IsolatedBuckBoostOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double primaryOutputVoltage = outputOperatingPoint.get_output_voltages()[0];
        double primaryOutputCurrent = outputOperatingPoint.get_output_currents()[0];
        double diodeVoltageDrop = get_diode_voltage_drop();
        double totalReflectedSecondaryCurrent = 0;
        for (size_t secondaryIndex = 0; secondaryIndex + 1 < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
            totalReflectedSecondaryCurrent += outputOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex];
        }

        double efficiency = 1;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }


        auto dutyCycle = calculate_duty_cycle(inputVoltage, primaryOutputVoltage, efficiency);

        auto tOn = dutyCycle / switchingFrequency;

        auto primaryCurrentPeakToPeak = (inputVoltage * primaryOutputVoltage) / (inputVoltage + primaryOutputVoltage) / (switchingFrequency * inductance);

        // Fly-Buck-Boost: continuous magnetizing current (K=1, Bimag probe sums
        // i(Lpri) + Σ i(Lsec)/n). In CCM, each cap's charge balance gives
        // ⟨i_winding⟩_off · (1−D) = I_out_winding, so ⟨i_pri⟩_off = I_out_pri/(1−D)
        // and ⟨i_sec⟩_off = I_out_sec/(1−D). Summing referred to primary:
        //   ⟨i_mag⟩ = (I_out_pri + Σ I_out_sec/n) / (1−D)
        // (Same buck-boost current scaling as Erickson §5.2 applied to all windings
        // sharing the magnetizing flux.) Earlier code dropped the /(1−D) factor,
        // under-estimating I_mag by (1−D)× — D=0.294 design 1: 0.567 A vs true 0.802 A.
        auto primaryCurrentAverage = (primaryOutputCurrent + totalReflectedSecondaryCurrent) / (1.0 - dutyCycle);

        // Fly-Buck-Boost primary winding voltage swing:
        //   ON  (S1 closed):    V(pri_in) = +Vin
        //   OFF (Dpri conducts): V(pri_in) clamps to V(vpri_out) − Vd ≈ −|Vout_pri|
        // → V_max = +Vin, V_min = −(Vout_pri + Vd), ptp = Vin + Vout_pri + Vd.
        auto primaryVoltaveMaximum = inputVoltage;
        auto primaryVoltaveMinimum = -(primaryOutputVoltage + diodeVoltageDrop);
        auto primaryVoltavePeaktoPeak = primaryVoltaveMaximum - primaryVoltaveMinimum;

        // Per-OP analytical diagnostics (mirrors IsolatedBuck.cpp).
        // R_eff_pri uses the asymptotic Im_avg = Iout + Σ Iout_sec/n
        // (the same quantity the model uses to centre the triangular
        // primary current waveform), giving R_eff_pri = Vpri / Im_avg.
        // CCM criterion: K = 2·L·Fs / R_eff_pri ≥ K_crit = (1−D)²
        // (Erickson §5.2 buck-boost). This is the standard buck-boost
        // boundary applied to the magnetizing inductance referred to
        // the primary side.
        double rEffPri = (primaryCurrentAverage > 0)
            ? (primaryVoltaveMaximum / primaryCurrentAverage)
            : std::numeric_limits<double>::infinity();
        double kFactor = 2.0 * inductance * switchingFrequency / rEffPri;
        double kCrit   = (1.0 - dutyCycle) * (1.0 - dutyCycle);
        lastDutyCycle                 = dutyCycle;
        lastMagnetizingCurrentRipple  = primaryCurrentPeakToPeak;
        lastPrimaryAverageCurrent     = primaryCurrentAverage;
        lastPrimaryPeakCurrent        = primaryCurrentAverage + primaryCurrentPeakToPeak / 2;
        lastIsCcm                     = (kFactor >= kCrit);
        lastSecondaryPeakCurrent      = 0.0;  // updated in the secondaries loop below

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, primaryCurrentAverage, 0);
            voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, 0);

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // Secondaries
        for (size_t secondaryIndex = 0; secondaryIndex + 1 < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
            Waveform currentWaveform;
            Waveform voltageWaveform;
            double secondaryOutputCurrent = outputOperatingPoint.get_output_currents()[secondaryIndex + 1];

            auto secondaryCurrentMaximum = (1 + dutyCycle) / (1 - dutyCycle) * secondaryOutputCurrent - secondaryOutputCurrent;
            auto secondaryCurrentMinimum = 0;
            auto secondaryCurrentPeakToPeak = secondaryCurrentMaximum - secondaryCurrentMinimum;

            // Track worst-case secondary peak across all secondaries
            // for diagnostic readout.
            double secPeak = secondaryOutputCurrent + secondaryCurrentMaximum;
            if (secPeak > lastSecondaryPeakCurrent) {
                lastSecondaryPeakCurrent = secPeak;
            }

            // Secondary winding levels follow the primary winding (+Vin during ON,
            // clamped to -(Vout_pri + Vd) during OFF, see the comment block above)
            // scaled by 1/n with inverted dot. The diode only conducts during OFF,
            // so its drop applies to the maximum (OFF) level, not to the blocking
            // (ON) level. Volt-second balance: D*(-Vin/n) + (1-D)*((Vout_pri+Vd)/n) = 0.
            auto secondaryVoltaveMaximum = (primaryOutputVoltage + diodeVoltageDrop) / turnsRatios[secondaryIndex] - diodeVoltageDrop;
            auto secondaryVoltaveMinimum = -inputVoltage / turnsRatios[secondaryIndex];
            auto secondaryVoltavePeaktoPeak = secondaryVoltaveMaximum - secondaryVoltaveMinimum;

            currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, secondaryCurrentPeakToPeak, switchingFrequency, 1.0 - dutyCycle, secondaryOutputCurrent, 0, tOn);
            voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, secondaryVoltavePeaktoPeak, switchingFrequency, 1.0 - dutyCycle, 0, 0, tOn);

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Secondary " + std::to_string(secondaryIndex));
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        // Per-OP diagnostic snapshot.
        perOpDutyCycle.push_back(lastDutyCycle);
        perOpMagnetizingCurrentRipple.push_back(lastMagnetizingCurrentRipple);
        perOpPrimaryAverageCurrent.push_back(lastPrimaryAverageCurrent);
        perOpPrimaryPeakCurrent.push_back(lastPrimaryPeakCurrent);
        perOpSecondaryPeakCurrent.push_back(lastSecondaryPeakCurrent);
        perOpIsCcm.push_back(lastIsCcm);

        return operatingPoint;
    }

    bool IsolatedBuckBoost::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        for (size_t isolatedbuckBoostOperatingPointIndex = 1; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
            if (get_operating_points()[isolatedbuckBoostOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[isolatedbuckBoostOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
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

    DesignRequirements IsolatedBuckBoost::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double efficiency = 1;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Missing both current ripple ratio and maximum swtich current");
        }

        // Turns ratio calculation
        size_t numOutputVoltages = get_operating_points()[0].get_output_voltages().size();
        if (numOutputVoltages < 2) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "IsolatedBuckBoost requires at least 2 output voltages (primary + secondary)");
        }
        std::vector<double> turnsRatios(numOutputVoltages - 1, 0);
        for (auto isolatedbuckBoostOperatingPoint : get_operating_points()) {
            double primaryVoltage = isolatedbuckBoostOperatingPoint.get_output_voltages()[0];
            for (size_t secondaryIndex = 0; secondaryIndex + 1 < isolatedbuckBoostOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                auto turnsRatio = primaryVoltage / (isolatedbuckBoostOperatingPoint.get_output_voltages()[secondaryIndex + 1] + get_diode_voltage_drop());
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
                for (size_t secondaryIndex = 0; secondaryIndex + 1 < isolatedbuckOperatingPoint.get_output_currents().size(); ++secondaryIndex) {
                    totalReflectedSecondaryCurrent += isolatedbuckOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex];
                }
                totalCurrentInPoint += totalReflectedSecondaryCurrent;
                maximumOutputCurrent = std::max(maximumOutputCurrent, totalCurrentInPoint);
            }

            maximumCurrentRiple = currentRippleRatio * maximumOutputCurrent;
        }

        if (get_maximum_switch_current()) {
            // Worst case across operating points is the SMALLEST ripple budget,
            // which produces the largest needed inductance below.
            double worstCaseRippleBudget = std::numeric_limits<double>::max();
            for (auto isolatedbuckBoostOperatingPoint : get_operating_points()) {
                auto primaryCurrent = isolatedbuckBoostOperatingPoint.get_output_currents()[0];
                double totalReflectedSecondaryCurrent = 0;
                for (size_t secondaryIndex = 0; secondaryIndex + 1 < isolatedbuckBoostOperatingPoint.get_output_currents().size(); ++secondaryIndex) {
                    totalReflectedSecondaryCurrent += isolatedbuckBoostOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex];
                }
                double primaryOutputVoltage = isolatedbuckBoostOperatingPoint.get_output_voltages()[0];
                auto dutyCycle = calculate_duty_cycle(minimumInputVoltage, primaryOutputVoltage, efficiency);
                worstCaseRippleBudget = std::min(worstCaseRippleBudget, get_maximum_switch_current().value() - (primaryCurrent + totalReflectedSecondaryCurrent) / (1 - dutyCycle));
            }
            if (!get_operating_points().empty()) {
                maximumCurrentRiple = worstCaseRippleBudget;
            }
        }

        double maximumNeededInductance = 0;
        for (size_t isolatedbuckBoostOperatingPointIndex = 0; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
            auto isolatedbuckBoostOperatingPoint = get_operating_points()[isolatedbuckBoostOperatingPointIndex];
            double primaryVoltage = isolatedbuckBoostOperatingPoint.get_output_voltages()[0];
            auto switchingFrequency = isolatedbuckBoostOperatingPoint.get_switching_frequency();
            auto desiredInductance = primaryVoltage * maximumInputVoltage / (primaryVoltage + maximumInputVoltage) / (2 * maximumCurrentRiple* switchingFrequency);

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
        designRequirements.set_topology(Topologies::ISOLATED_BUCK_BOOST_CONVERTER);
        return designRequirements;
    }

    std::vector<MAS::OperatingPoint> IsolatedBuckBoost::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        extraLoVoltageWaveforms.clear();
        extraLoCurrentWaveforms.clear();
        extraLoInductances.clear();

        // IsolatedBuckBoost: the primary winding IS the converter inductor.
        // Secondary windings are flyback-type outputs and do not require a separate output inductor Lo.
        // No Lo waveforms are captured here.

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



        // Clear per-OP diagnostic vectors so the wizard table reflects this run only.
        perOpName.clear();
        perOpDutyCycle.clear();
        perOpMagnetizingCurrentRipple.clear();
        perOpPrimaryAverageCurrent.clear();
        perOpPrimaryPeakCurrent.clear();
        perOpSecondaryPeakCurrent.clear();
        perOpIsCcm.clear();
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t isolatedbuckBoostOperatingPointIndex = 0; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
                std::string opName = inputVoltagesNames[inputVoltageIndex];
                if (get_operating_points().size() > 1) {
                    opName += " · OP" + std::to_string(isolatedbuckBoostOperatingPointIndex);
                }
                perOpName.push_back(opName);
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[isolatedbuckBoostOperatingPointIndex], turnsRatios, magnetizingInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(isolatedbuckBoostOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<std::variant<Inputs, CAS::Inputs>> IsolatedBuckBoost::get_extra_components_inputs(
        ExtraComponentsMode /*mode*/,
        std::optional<Magnetic> /*magnetic*/)
    {
        // IsolatedBuckBoost: the transformer primary winding acts as the converter inductor.
        // Secondary windings are flyback-type outputs and do not require a separate output inductor Lo.
        return {};
    }

    std::vector<MAS::OperatingPoint> IsolatedBuckBoost::process_operating_points(Magnetic magnetic) {
        IsolatedBuckBoost::run_checks(_assertErrors);

        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedIsolatedBuckBoost::process() {
        IsolatedBuckBoost::run_checks(_assertErrors);

        Inputs inputs;

        double maximumNeededInductance = get_desired_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;


        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
            inputVoltagesNames.push_back("Nom.");
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
            inputVoltagesNames.push_back("Max.");
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
            inputVoltagesNames.push_back("Min.");
        }

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
        designRequirements.set_topology(Topologies::ISOLATED_BUCK_BOOST_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t isolatedbuckBoostOperatingPointIndex = 0; isolatedbuckBoostOperatingPointIndex < get_operating_points().size(); ++isolatedbuckBoostOperatingPointIndex) {
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[isolatedbuckBoostOperatingPointIndex], turnsRatios, maximumNeededInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(isolatedbuckBoostOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

    std::string IsolatedBuckBoost::generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {

        // Minimal spice_config() consolidation — PULSE rise/fall only.
        const auto cfg = spice_config();

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
        const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;  // Steady state + extraction
        double simTime = numPeriodsTotal * period;
        double startTime = numSteadyStatePeriods * period;  // Start extracting after steady state
        double stepTime = period / 200;
        
        circuit << "* Isolated Buck-Boost Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vout=" << primaryOutputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, " << numSecondaries << " secondaries\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";

        // §8a.5 input-current sense: a zero-V source upstream of S1's
        // drain isolates the true switch-channel current. IsolatedBuck-
        // Boost is flyback-like: only S1 connects vin_dc to the
        // primary; the primary diode Dpri returns to a NEGATIVE rail
        // (vpri_out < 0), not to vin_dc, so no return path pollutes the
        // input-current measurement. Sign: + terminal is vin_dc, - is
        // q1_drain, so positive i(Vq1_sense) means current flowing OUT
        // of vin_dc INTO S1 (non-negative in CCM/DCM).
        circuit << "* §8a.5 input-current sense\n";
        circuit << "Vq1_sense vin_dc q1_drain 0\n\n";

        // PWM Switch
        circuit << "* PWM Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << cfg.pwmRise << " " << cfg.pwmFall
                << " " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                << " RON=" << cfg.swModelRON
                << " ROFF=" << std::scientific << cfg.swModelROFF << std::defaultfloat << "\n";
        circuit << "S1 q1_drain pri_p pwm_ctrl 0 SW1\n\n";
        
        // Primary current sense - Vpri_sense only needed for switch node connectivity
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense pri_p pri_in 0\n\n";

        circuit << "* Coupled Inductor (Primary = flyback-type primary winding)\n";
        circuit << "Lpri pri_in 0 " << std::scientific << magnetizingInductance << std::fixed << "\n";
        
        // Secondary windings: dot at gnd (0) side, OPPOSITE to primary dot at pri_in.
        // SPICE coupled-inductor convention: the "dot" terminal is the FIRST node of
        // each L<x> element. To get flyback polarity (secondary blocked during ON,
        // conducting during OFF), we must list 0 as the first node of Lsec so the
        // dot lands at the grounded end. Earlier wiring `Lsec<i> sec<i>_in 0` placed
        // both dots on the same side (top), producing forward-converter behaviour and
        // breaking the magnetizing-current probe used in the analytical-vs-SPICE
        // PtP comparison (NRMSE ~0.9 — discontinuous Bimag at switch transitions).
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double secondaryInductance = magnetizingInductance / (turnsRatios[secIdx] * turnsRatios[secIdx]);
            circuit << "Lsec" << secIdx << " 0 sec" << secIdx << "_in " << std::scientific << secondaryInductance << std::fixed << "\n";
        }
        
        // Couple primary to each secondary
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "Kpri_sec" << secIdx << " Lpri Lsec" << secIdx << " 1\n";
        }
        // Couple secondaries to each other (for proper cross-regulation)
        for (size_t i = 0; i < numSecondaries; ++i) {
            for (size_t j = i + 1; j < numSecondaries; ++j) {
                circuit << "Ksec" << i << "_" << j << " Lsec" << i << " Lsec" << j << " 1\n";
            }
        }
        circuit << "\n";
        
        // Diode model - relaxed for convergence
        circuit << "* Diode model\n";
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::defaultfloat;
        if (!cfg.diodeExtra.empty()) circuit << " " << cfg.diodeExtra;
        circuit << ")\n\n";
        
        // Primary output: inverting (Fly-Buck-Boost) topology per TI SNVAA84
        // ("negative output for fly-buck-boost"). Canonical inverting buck-boost:
        //   Vin → S1 → pri_in → Lpri → 0   (during ON, i_L ramps from low to high)
        //   vpri_out (negative cap) → Dpri → pri_in   (during OFF, Dpri conducts;
        //     V(pri_in) clamps to V(vpri_out)≈−|Vout|, inductor discharges into cap,
        //     making vpri_out MORE negative)
        // Earlier wiring (Dpri pri_in→vpri_rect with IC=+Vout) is non-physical: during
        // OFF V(pri_in) goes negative, so a forward-biased Dpri would need vpri_rect
        // even more negative; with +5V IC it was permanently reverse-biased, the
        // primary cap never charged, all power went through Dsec only, and the
        // converter ran in DCM with NRMSE 17-36% on the magnetizing-current shape.
        circuit << "* Primary Output Stage (inverting Buck-Boost; Vout < 0)\n";
        circuit << "Dpri vpri_rect pri_in DIDEAL\n";
        circuit << "Vpri_out_sense vpri_rect vpri_out 0\n";
        double primaryLoadResistance = primaryOutputVoltage / primaryOutputCurrent;
        // 50 mΩ output-cap ESR (physical electrolytic spec) damps the
        // primary-inductor / output-cap LC ring that otherwise took
        // ~1000+ switching periods to settle. With ESR in place the
        // ring damps in <500 periods. The ESR is small enough that the
        // steady-state cap voltage is unaffected to <0.1 %.
        circuit << "Cpri vpri_out vpri_out_esr 100u IC=" << (-primaryOutputVoltage) << "\n";
        circuit << "Rpri_esr vpri_out_esr 0 0.1\n";
        circuit << "Rload_pri vpri_out 0 " << primaryLoadResistance << "\n\n";

        // Secondary output stages
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "* Secondary " << secIdx << " output stage\n";
            circuit << "Rsec" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_node 0.01\n";
            circuit << "Dsec" << secIdx << " sec" << secIdx << "_node sec" << secIdx << "_rect DIDEAL\n";
            circuit << "Vsec_sense" << secIdx << " sec" << secIdx << "_rect vout" << secIdx << " 0\n";

            double outputVoltage = opPoint.get_output_voltages()[secIdx + 1];
            double outputCurrent = opPoint.get_output_currents()[secIdx + 1];
            double loadResistance = outputVoltage / outputCurrent;
            // Same ESR + cap damping pattern as the primary output stage.
            circuit << "Cout" << secIdx << " vout" << secIdx << " vout" << secIdx << "_esr 100u IC=" << outputVoltage << "\n";
            circuit << "Rout" << secIdx << "_esr vout" << secIdx << "_esr 0 0.1\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }
        
        // Magnetizing current probe:
        // I_mag = I(Lpri) during ON time (secondary blocked)
        //       + sum_k( I(Lsec_k) / turnsRatio_k ) during OFF time (primary switch open)
        // Combined gives a continuous triangular waveform matching the analytical TRIANGULAR model.
        circuit << "* Magnetizing current probe (continuous triangular waveform)\n";
        circuit << "Bimag 0 imag_p I = I(Lpri)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " + I(Lsec" << secIdx << ")/" << std::fixed << std::setprecision(6) << turnsRatios[secIdx];
        }
        circuit << "\n";
        circuit << "Vimag_sense imag_p 0 0\n\n";

        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << "\n\n";

        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(vin_dc) i(Vin) i(Vq1_sense) v(pri_in) i(Vpri_sense) i(Vimag_sense) i(Lpri)"
                << " v(vpri_out) i(Vpri_out_sense)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) i(Vsec_sense" << secIdx << ") v(vout" << secIdx << ")";
        }
        circuit << "\n\n";
        
        // Options with convergence aids for flyback topology
        // RSHUNT: prevents floating nodes (1 TΩ to ground from each node)
        // RSERIES: small series resistance for inductors (0.1 mΩ)
        // ITL4=1000: more iterations for transient convergence
        // Solver tolerances from cfg; RSHUNT/RSERIES kept hardcoded
        // (topology-specific flyback-convergence aids, not generic
        // SpiceSimulationConfig fields).
        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::defaultfloat
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4
                << " RSHUNT=1e12 RSERIES=1e-4\n";
        
        // Initial conditions - use .nodeset for better DC convergence
        circuit << ".nodeset v(pri_in)=0 v(vpri_out)=" << (-primaryOutputVoltage) << "\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << ".nodeset v(vout" << secIdx << ")=" << opPoint.get_output_voltages()[secIdx + 1] << "\n";
        }
        circuit << "\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> IsolatedBuckBoost::simulate_and_extract_operating_points(
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
                auto ibbOpPoint = get_operating_points()[opIndex];
                
                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = ibbOpPoint.get_switching_frequency();
                
                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);

                if (!simResult.success) {
                    throw std::runtime_error("Simulation failed: " + simResult.errorMessage);
                }
                
                // Define waveform name mapping
                NgspiceRunner::WaveformNameMapping waveformMapping;

                // Primary winding: magnetizing current probe (Bimag = I(Lpri) + ΣI(Lsec_k)/n_k).
                // This continuous triangular waveform matches the analytical TRIANGULAR
                // model used in process_operating_points (primary winding current is
                // modeled as the magnetizing current, not the switch / inductor current).
                waveformMapping.push_back({{"voltage", "pri_in"}, {"current", "vimag_sense#branch"}});

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
                // With the secondary dot now at gnd (see Lsec orientation comment
                // above), the magnetizing-current probe matches the analytical sign
                // convention without flipping. Secondaries do not need flipping either —
                // the FLYBACK_PRIMARY analytical waveform sits on the same sign axis as
                // i(Vsec_sense<k>) (current flowing from sec_rect into vout during OFF).
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    ibbOpPoint.get_ambient_temperature(),
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
    
    std::vector<ConverterWaveforms> IsolatedBuckBoost::simulate_and_extract_topology_waveforms(const std::vector<double>& turnsRatios,
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
            double switchingFrequency = opPoint.get_switching_frequency();
            
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
            
            // §5.1 / §8a.5 converter-port stream:
            //   input_voltage = v(vin_dc) (real DC rail, constant Vin).
            //     NOT v(pri_in) — that's the post-switch / winding-top
            //     node which pulses to ~Vin during ON and to ≈−|Vout|
            //     during OFF (Dpri clamps to vpri_out). v(pri_in) is a
            //     winding-port quantity. Bug A doesn't apply — we were
            //     already on vin_dc.
            //   input_current = i(Vq1_sense) — zero-V sense source
            //     upstream of S1's drain. Only S1 connects vin_dc to
            //     the primary winding; Dpri returns to the negative
            //     vpri_out rail (not vin_dc), so i(Vq1_sense) captures
            //     the full converter draw with no contamination from
            //     the primary diode loop. i(Vin) (the prior probe,
            //     sign-flipped) is numerically equivalent in the ideal
            //     case, but the §8a.5 template inserts a dedicated
            //     ammeter so future snubbers/bypass branches don't
            //     silently pollute the input-current reading. Sign:
            //     + terminal of Vq1_sense is vin_dc, so positive
            //     i(Vq1_sense) is current flowing OUT of vin_dc INTO
            //     S1 (non-negative in CCM/DCM).
            //   Bug C: Lpri is wired `Lpri pri_in 0` — far terminal IS
            //     GND, so v(pri_in) is a true differential winding-only
            //     voltage (used by simulate_and_extract_operating_points
            //     for the winding-port stream). Each Lsec<i> is wired
            //     `Lsec<i> 0 sec<i>_in` — dot at GND, so v(sec<i>_in)
            //     is also winding-only. No E-source needed on either
            //     side. (The post-rectifier vout<i> nodes are NEVER
            //     used as a winding-voltage probe.)
            wf.set_input_voltage(getWaveform("vin_dc"));

            Waveform iIn = getWaveform("vq1_sense#branch");
            if (iIn.get_data().empty()) {
                throw std::runtime_error(
                    "IsolatedBuckBoost simulate_and_extract_topology_waveforms: "
                    "missing i(Vq1_sense) — netlist/save mismatch");
            }

            // ngspice SW1 model produces narrow di/dt spikes (~10^5 A)
            // around switching edges that are a numerical artefact of
            // the discontinuous RON/ROFF transition, not a physical
            // bound. Clamp samples to ±2·max(|i(Vpri_sense)|) — the
            // primary current sets the actual conduction-current
            // scale; this preserves all real samples while discarding
            // the artefact spikes that would poison the DC mean used
            // downstream by §5.1.
            //
            // numerical guard against ngspice SW-model di/dt spikes
            // (~10^5 A), not a physical bound.
            {
                auto& d = iIn.get_mutable_data();
                Waveform iPri = getWaveform("vpri_sense#branch");
                double iMax = 0.0;
                for (double v : iPri.get_data()) iMax = std::max(iMax, std::abs(v));
                if (iMax > 0.0) {
                    const double clamp = 2.0 * iMax;
                    for (auto& v : d) {
                        if (v >  clamp) v =  clamp;
                        if (v < -clamp) v = -clamp;
                    }
                }
            }
            wf.set_input_current(iIn);
            
            // Output port: index 0 = primary buck-boost output (vpri_out,
            // ~Vin·D/(1-D)), indices 1..N = transformer secondaries
            // (vout0..vout(N-1)). Currents reconstructed from Vout/Rload
            // — the ammeters sit in the cap-return path so their branch
            // current averages to zero in steady state. Mirrors §5.1
            // template across all isolated topologies.
            {
                const double Vo_pri = opPoint.get_output_voltages()[0];
                const double Io_pri = opPoint.get_output_currents()[0];
                const double Rload_pri = std::max(Vo_pri / Io_pri, 1e-3);
                wf.get_mutable_output_voltages().push_back(getWaveform("vpri_out"));
                Waveform ioPriWf = getWaveform("vpri_out");
                auto& ioPriData = ioPriWf.get_mutable_data();
                for (auto& v : ioPriData) v = v / Rload_pri;
                wf.get_mutable_output_currents().push_back(ioPriWf);
            }
            for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
                const double Vo_sec = opPoint.get_output_voltages()[secIdx + 1];
                const double Io_sec = opPoint.get_output_currents()[secIdx + 1];
                const double Rload_sec = std::max(Vo_sec / Io_sec, 1e-3);
                std::string si = std::to_string(secIdx);
                wf.get_mutable_output_voltages().push_back(getWaveform("vout" + si));
                Waveform ioSecWf = getWaveform("vout" + si);
                auto& ioSecData = ioSecWf.get_mutable_data();
                for (auto& v : ioSecData) v = v / Rload_sec;
                wf.get_mutable_output_currents().push_back(ioSecWf);
            }
            
            results.push_back(wf);
        }
    }
    
    // Restore original value
    set_num_periods_to_extract(originalNumPeriodsToExtract);
    
    return results;
}
} // namespace OpenMagnetics
