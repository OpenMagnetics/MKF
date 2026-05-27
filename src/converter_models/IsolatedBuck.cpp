#include "converter_models/IsolatedBuck.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double IsolatedBuck::calculate_duty_cycle(double inputVoltage, double outputVoltage, double efficiency) {
        auto dutyCycle = outputVoltage / inputVoltage * efficiency;
        if (dutyCycle >= 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Duty cycle must be smaller than 1");
        }
        // Explicit maxD gate (no silent clamp). 1 % tolerance absorbs
        // design-requirements rounding. Same pattern as Flyback, forward
        // family, Buck, Boost.
        if (maximumDutyCycle.has_value()) {
            const double dMax = maximumDutyCycle.value();
            constexpr double dutyTolerance = 0.01;
            if (dutyCycle > dMax * (1.0 + dutyTolerance)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "IsolatedBuck: required dutyCycle " + std::to_string(dutyCycle) +
                    " exceeds maximumDutyCycle " + std::to_string(dMax) +
                    " at Vin=" + std::to_string(inputVoltage) +
                    ", V_pri=" + std::to_string(outputVoltage) +
                    ". Raise Vin_min, lower V_pri, or relax maximumDutyCycle.");
            }
        }
        return dutyCycle;
    }

    IsolatedBuck::IsolatedBuck(const json& j) {
        from_json(j, *this);
    }

    AdvancedIsolatedBuck::AdvancedIsolatedBuck(const json& j) {
        from_json(j, *this);
    }

    OperatingPoint IsolatedBuck::process_operating_point_for_input_voltage(double inputVoltage, const IsolatedBuckOperatingPoint& outputOperatingPoint, const std::vector<double>& turnsRatios, double inductance) {

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
        auto period = 1.0 / switchingFrequency;

        auto tOn = dutyCycle * period;
        auto magnetizingCurrentRipple = (inputVoltage - primaryOutputVoltage) * tOn / inductance;
        // Magnetizing-current bounds (continuous symmetric triangle anchored
        // on Iout_pri + Σ Iout_sec/n, i.e. i_mag avg). i_mag = i_pri during
        // ON-time (secondary blocked), so primary peak = i_mag_max.
        auto magnetizingCurrentAverage = primaryOutputCurrent + totalReflectedSecondaryCurrent;
        auto magnetizingCurrentMax     = magnetizingCurrentAverage + magnetizingCurrentRipple / 2;
        auto magnetizingCurrentMin     = magnetizingCurrentAverage - magnetizingCurrentRipple / 2;
        // OFF-time secondary inductor current (avg over OFF-window) by
        // Cout_sec charge balance: ∫i_sec dt = Iout_sec · T → i_sec_off_avg
        // = Iout_sec / (1-D). Reflected into the primary winding this is
        // subtracted from i_mag during the OFF interval (i_pri = i_mag −
        // i_sec/n, dot convention). At switching edges, i_sec snaps
        // 0 ↔ I_sec_off_avg/n which is a STEP on the primary winding
        // current (the magnetizing flux is continuous, individual winding
        // currents are not — this is intrinsic to K=1 coupled inductors).
        auto reflectedSecondaryOffsetOff = 0.0;
        for (size_t secondaryIndex = 0; secondaryIndex + 1 < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
            const double Isec = outputOperatingPoint.get_output_currents()[secondaryIndex + 1];
            reflectedSecondaryOffsetOff += (Isec / (1.0 - dutyCycle)) / turnsRatios[secondaryIndex];
        }
        // Bookkeeping aliases kept for downstream code (peak / pp).
        auto primaryCurrentMaximum     = magnetizingCurrentMax;
        auto primaryCurrentMinimum     = magnetizingCurrentMin - reflectedSecondaryOffsetOff;
        (void)primaryCurrentMinimum;   // exposed via lastPrimary* diagnostics only

        auto primaryVoltageMaximum = inputVoltage - primaryOutputVoltage;
        auto primaryVoltageMinimum = -primaryOutputVoltage;
        auto primaryVoltagePeaktoPeak = primaryVoltageMaximum - primaryVoltageMinimum;

        // Populate per-OP analytical diagnostics. Exposed via
        // get_last_*() accessors so tests / Web frontend can verify the
        // analytical model without re-deriving the equations. CCM check
        // uses the standard buck criterion on the primary (magnetizing)
        // inductor: K = 2·L·Fs/R_load_pri ≥ 1−D.
        double rLoadPri = (primaryOutputCurrent > 0)
            ? (primaryOutputVoltage / primaryOutputCurrent)
            : std::numeric_limits<double>::infinity();
        double kFactor = 2.0 * inductance * switchingFrequency / rLoadPri;
        double kCrit   = 1.0 - dutyCycle;
        lastDutyCycle                 = dutyCycle;
        lastMagnetizingCurrentRipple  = magnetizingCurrentRipple;
        lastPrimaryAverageCurrent     = primaryOutputCurrent + totalReflectedSecondaryCurrent;
        lastPrimaryPeakCurrent        = primaryCurrentMaximum;
        lastIsCcm                     = (kFactor >= kCrit);
        lastSecondaryPeakCurrent      = 0.0;  // updated in the secondaries loop below

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            // PRIMARY WINDING current is the Lpri inductor current and has
            // a PIECEWISE shape (NOT a symmetric triangle) because at K=1
            // each winding current can step at switching events while the
            // magnetizing flux (and hence i_mag) stays continuous.
            //
            // ON-time  (S1 closed, secondary diode reverse-biased):
            //   i_sec = 0 → i_pri = i_mag, ramps from i_mag_min → i_mag_max
            //   at slope (Vin − Vout_pri)/Lpri.
            //
            // At t = D·T (ON→OFF): secondary diode commutates ON.
            //   i_sec jumps 0 → I_sec_off_avg ≈ Iout_sec/(1−D) (Cout_sec
            //   charge balance). Primary winding loses i_sec/n of its
            //   current: STEP DOWN by Σ Iout_sec/((1−D)·n).
            //
            // OFF-time (S1 open, S2/synchronous freewheel):
            //   i_pri = i_mag(t) − Σ i_sec/n, ramps down at the SAME slope
            //   −Vout_pri/Lpri as i_mag (since the reflected offset is
            //   approximately constant over the OFF window for small
            //   Cout-driven i_sec ripple).
            //
            // At t = T (OFF→ON): secondary diode reverse-biases.
            //   i_sec snaps back to 0, i_pri STEPS UP by the same offset
            //   and lands on i_mag_min, ready for the next on-ramp.
            //
            // Period-average check (KCL at vpri_out):
            //   ⟨i_pri⟩ = D·(i_mag_min + i_mag_max)/2
            //             + (1−D)·[(i_mag_min + i_mag_max)/2 − I_off]
            //           = ⟨i_mag⟩ − (1−D)·I_off
            //           = (Iout_pri + Σ Iout_sec/n) − Σ Iout_sec/n
            //           = Iout_pri ✓
            //
            // Replaces the older symmetric-TRIANGULAR(magnetizingRipple,
            // centred on Iout_pri) shape, which ignored the off-time
            // reflected-secondary step and produced ≥30 % NRMSE versus
            // SPICE on flybuck reference designs (LM5160 / LM5017 / LM5160-Q1
            // — see TestIsolatedBuckReferenceDesignsPtp.cpp).
            const double iPriOnStart   = magnetizingCurrentMin;
            const double iPriOnEnd     = magnetizingCurrentMax;
            const double iPriOffStart  = magnetizingCurrentMax - reflectedSecondaryOffsetOff;
            const double iPriOffEnd    = magnetizingCurrentMin - reflectedSecondaryOffsetOff;
            {
                std::vector<double> data = {iPriOnStart, iPriOnEnd, iPriOffStart, iPriOffEnd};
                std::vector<double> time = {0.0,         tOn,       tOn,         period};
                currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                currentWaveform.set_data(data);
                currentWaveform.set_time(time);
            }
            voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, 0);

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

            // Track the worst-case secondary peak across all secondaries
            // for diagnostic readout.
            double secPeak = secondaryOutputCurrent + secondaryCurrentMaximum;
            if (secPeak > lastSecondaryPeakCurrent) {
                lastSecondaryPeakCurrent = secPeak;
            }

            // Flybuck-mode secondary winding voltage (NOT forward-mode).
            //
            // Topology refresher: in IsolatedBuck (flybuck), the primary IS
            // the buck inductor. The secondary winding is COUPLED into that
            // same magnetic, with a flyback-style rectifier on the output.
            // The secondary diode therefore conducts during the FREEWHEEL
            // (OFF) interval — the opposite of a forward converter.
            //
            // Per CONVERTER_MODELS_GOLDEN_GUIDE.md §5.0, this signal lands in
            // excitations_per_winding[secondary].voltage and is checked by
            // Test_VoltSecondBalance_IsolatedBuck. Pre-fix this code emitted
            // the forward-mode formula (max during ON = +(Vin−Vp)/n−Vd, min
            // during OFF = −Vp/n+Vd), which violates volt-second balance by
            // ~89 % of peak — observable as the failing analytical check.
            //
            //   ON  (S1 closed): V_pri,winding = +(Vin − Vp_pri)
            //                    V_sec,winding = −(Vin − Vp_pri)/n
            //                    (dot of Lsec is at GND in the netlist, so
            //                     it inverts; secondary diode reverse-biased)
            //   OFF (freewheel): V_pri,winding = −Vp_pri
            //                    V_sec,winding = +Vp_pri/n
            //                    Diode then conducts and the winding terminal
            //                    clamps to Vout_sec + Vd, i.e. the winding-
            //                    side voltage seen by the inductor port is
            //                    +Vp_pri/n minus the resistive Vd loss.
            //
            // Sanity (ideal, Vd→0, D = Vp_pri/Vin):
            //   D·V_on + (1−D)·V_off
            //     = D·[−(Vin−Vp)/n] + (1−D)·[Vp/n]
            //     = (1/n)·[−D·Vin + D·Vp + Vp − D·Vp]
            //     = (1/n)·[Vp − D·Vin]  =  0.    ✓
            auto secondaryVoltageMinimum = -(inputVoltage - primaryOutputVoltage) / turnsRatios[secondaryIndex];
            auto secondaryVoltageMaximum =  primaryOutputVoltage / turnsRatios[secondaryIndex] - diodeVoltageDrop;

            // Current
            {
                std::vector<double> data = {
                    0,
                    0,
                    secondaryOutputCurrent + secondaryCurrentMinimum,
                    secondaryOutputCurrent + secondaryCurrentMaximum
                };
                std::vector<double> time = {
                    0,
                    tOn,
                    tOn,
                    period
                };
                currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                currentWaveform.set_data(data);
                currentWaveform.set_time(time);
            }
            // Voltage
            {
                std::vector<double> data = {
                    secondaryVoltageMinimum,
                    secondaryVoltageMinimum,
                    secondaryVoltageMaximum,
                    secondaryVoltageMaximum
                };
                std::vector<double> time = {
                    0,
                    tOn,
                    tOn,
                    period
                };
                voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
                voltageWaveform.set_data(data);
                voltageWaveform.set_time(time);
            }

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

    bool IsolatedBuck::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        for (size_t isolatedbuckOperatingPointIndex = 1; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
            if (get_operating_points()[isolatedbuckOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[isolatedbuckOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
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

    DesignRequirements IsolatedBuck::process_design_requirements() {
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Missing both current ripple ratio and maximum switch current");
        }

        // Turns ratio calculation
        size_t numOutputVoltages = get_operating_points()[0].get_output_voltages().size();
        if (numOutputVoltages < 2) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "IsolatedBuck requires at least 2 output voltages (primary + secondary)");
        }
        std::vector<double> turnsRatios(numOutputVoltages - 1, 0);
        for (auto isolatedbuckOperatingPoint : get_operating_points()) {
            double primaryVoltage = isolatedbuckOperatingPoint.get_output_voltages()[0];
            for (size_t secondaryIndex = 0; secondaryIndex + 1 < isolatedbuckOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                auto turnsRatio = primaryVoltage / (isolatedbuckOperatingPoint.get_output_voltages()[secondaryIndex + 1] + get_diode_voltage_drop());
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
            for (auto isolatedbuckOperatingPoint : get_operating_points()) {
                auto primaryCurrent = isolatedbuckOperatingPoint.get_output_currents()[0];
                double totalReflectedSecondaryCurrent = 0;
                for (size_t secondaryIndex = 0; secondaryIndex + 1 < isolatedbuckOperatingPoint.get_output_currents().size(); ++secondaryIndex) {
                    totalReflectedSecondaryCurrent += isolatedbuckOperatingPoint.get_output_currents()[secondaryIndex + 1] / turnsRatios[secondaryIndex];
                }

                maximumCurrentRiple = (get_maximum_switch_current().value() - primaryCurrent - totalReflectedSecondaryCurrent) * 2;
            }
        }

        double maximumNeededInductance = 0;
        for (size_t isolatedbuckOperatingPointIndex = 0; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
            auto isolatedbuckOperatingPoint = get_operating_points()[isolatedbuckOperatingPointIndex];
            double primaryVoltage = isolatedbuckOperatingPoint.get_output_voltages()[0];
            auto switchingFrequency = isolatedbuckOperatingPoint.get_switching_frequency();
            auto desiredInductance =  (maximumInputVoltage - primaryVoltage) * primaryVoltage / (maximumInputVoltage * switchingFrequency * maximumCurrentRiple);
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
        designRequirements.set_topology(Topologies::ISOLATED_BUCK_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> IsolatedBuck::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        // Clear per-OP diagnostic vectors so the wizard table reflects this run only.
        perOpName.clear();
        perOpDutyCycle.clear();
        perOpMagnetizingCurrentRipple.clear();
        perOpPrimaryAverageCurrent.clear();
        perOpPrimaryPeakCurrent.clear();
        perOpSecondaryPeakCurrent.clear();
        perOpIsCcm.clear();


        extraLoVoltageWaveforms.clear();
        extraLoCurrentWaveforms.clear();
        extraLoInductances.clear();

        // IsolatedBuck (Fly-Buck): the primary winding IS the inductor.
        // Secondary windings are flyback-coupled outputs; they don't require a separate Lo.
        // No Lo waveforms are captured here.

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t isolatedbuckOperatingPointIndex = 0; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
                std::string opName = inputVoltagesNames[inputVoltageIndex];
                if (get_operating_points().size() > 1) {
                    opName += " · OP" + std::to_string(isolatedbuckOperatingPointIndex);
                }
                perOpName.push_back(opName);
                auto operatingPoint = process_operating_point_for_input_voltage(inputVoltage, get_operating_points()[isolatedbuckOperatingPointIndex], turnsRatios, magnetizingInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(isolatedbuckOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<std::variant<Inputs, CAS::Inputs>> IsolatedBuck::get_extra_components_inputs(
        ExtraComponentsMode /*mode*/,
        std::optional<Magnetic> /*magnetic*/)
    {
        // IsolatedBuck (Fly-Buck): the transformer primary winding acts as the converter inductor.
        // Secondary windings are flyback-type outputs and do not require a separate output inductor Lo.
        return {};
    }

    std::vector<OperatingPoint> IsolatedBuck::process_operating_points(Magnetic magnetic) {
        IsolatedBuck::run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedIsolatedBuck::process() {
        IsolatedBuck::run_checks(_assertErrors);

        Inputs inputs;

        double maximumNeededInductance = get_desired_inductance();
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
        inductanceWithTolerance.set_nominal(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::ISOLATED_BUCK_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t isolatedbuckOperatingPointIndex = 0; isolatedbuckOperatingPointIndex < get_operating_points().size(); ++isolatedbuckOperatingPointIndex) {
                auto operatingPoint = process_operating_point_for_input_voltage(inputVoltage, get_operating_points()[isolatedbuckOperatingPointIndex], turnsRatios, maximumNeededInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(isolatedbuckOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

    std::vector<ConverterWaveforms> IsolatedBuck::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
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
            //     NOT v(sw_node) (the switch node swings between vin_dc
            //     during S1 on-time and 0 during S2 on-time / synchronous
            //     freewheel; using it would be a square wave). Bug A
            //     doesn't apply — we were already on vin_dc.
            //   input_current = i(Vq1_sense) — zero-V sense source
            //     upstream of S1's drain. S2 (sync rectifier) sits
            //     between sw_node and GND, never sourcing from vin_dc,
            //     so Vq1_sense captures the full converter draw.
            //     i(Vin) (prior probe) was equivalent in the ideal case,
            //     but the §8a.5 template inserts a dedicated sense
            //     source so future snubber / bypass additions don't
            //     silently pollute the input-current reading. Sign:
            //     + terminal of Vq1_sense is vin_dc, so positive
            //     i(Vq1_sense) is current flowing OUT of vin_dc INTO
            //     S1 (non-negative in CCM/DCM).
            //   Bug C: primary winding voltage probe is Bvpri_diff =
            //     v(pri_in)-v(vpri_out) — a true differential winding
            //     voltage. Secondary windings: each Lsec<i> is wired
            //     `0 sec<i>_in`, so v(sec<i>_in) is already referenced
            //     to GND = the other Lsec terminal, hence already
            //     winding-only (no E-source needed). The post-rectifier
            //     v(sec<i>_rect)/v(vout<i>) nodes are NOT used for
            //     winding-voltage extraction (per
            //     simulate_and_extract_operating_points mapping above).
            wf.set_input_voltage(getWaveform("vin_dc"));

            Waveform iIn = getWaveform("vq1_sense#branch");
            if (iIn.get_data().empty()) {
                throw std::runtime_error(
                    "IsolatedBuck simulate_and_extract_topology_waveforms: "
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

            // Output[0] = primary buck output (non-isolated, across Cpri).
            // Reconstruct Iout = Vout/Rload (DC by design).
            const auto& voutsNom = opPoint.get_output_voltages();
            const auto& ioutsNom = opPoint.get_output_currents();
            const size_t numSecondaries = turnsRatios.size();
            if (voutsNom.size() != numSecondaries + 1 || ioutsNom.size() != numSecondaries + 1) {
                throw std::runtime_error(
                    "IsolatedBuck::simulate_and_extract_topology_waveforms: "
                    "operating point output_voltages/currents size (" +
                    std::to_string(voutsNom.size()) + "/" + std::to_string(ioutsNom.size()) +
                    ") inconsistent with numSecondaries+1 (" + std::to_string(numSecondaries + 1) + ")");
            }

            // Primary buck output rail.
            wf.get_mutable_output_voltages().push_back(getWaveform("vpri_out"));
            {
                Waveform priIout = getWaveform("vpri_out");
                const double rLoad = voutsNom[0] / ioutsNom[0];
                for (auto& v : priIout.get_mutable_data()) v = v / rLoad;
                wf.get_mutable_output_currents().push_back(priIout);
            }

            // Each isolated secondary output cap.
            for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
                const std::string voutNode = "vout" + std::to_string(secIdx);
                wf.get_mutable_output_voltages().push_back(getWaveform(voutNode));
                Waveform secIout = getWaveform(voutNode);
                const double rLoad = voutsNom[secIdx + 1] / ioutsNom[secIdx + 1];
                for (auto& v : secIout.get_mutable_data()) v = v / rLoad;
                wf.get_mutable_output_currents().push_back(secIout);
            }

            results.push_back(wf);
        }
    }

    // Restore original value
    set_num_periods_to_extract(originalNumPeriodsToExtract);

    return results;
}

    std::string IsolatedBuck::generate_ngspice_circuit(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {

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
        const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;
        double simTime = numPeriodsTotal * period;
        double startTime = numSteadyStatePeriods * period;
        double stepTime = period / 200;

        circuit << "* Isolated Buck (Flybuck) Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vout_pri=" << primaryOutputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* Lmag=" << (magnetizingInductance*1e6) << "uH, " << numSecondaries << " secondaries\n\n";

        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";

        // §8a.5 input-current sense: a zero-V source upstream of S1's
        // drain isolates the true switch-channel current. The
        // synchronous-rectifier S2 sits sw_node -> 0 and does NOT draw
        // from vin_dc, so the only path from vin_dc is through S1; a
        // single ammeter there captures the converter draw. Sign:
        // + terminal is vin_dc, - is q1_drain, so positive i(Vq1_sense)
        // means current flowing OUT of vin_dc INTO S1 (non-negative in
        // CCM/DCM).
        circuit << "* §8a.5 input-current sense\n";
        circuit << "Vq1_sense vin_dc q1_drain 0\n\n";

        // High-side switch (synchronous buck)
        circuit << "* High-side Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << std::scientific << cfg.pwmRise << " " << cfg.pwmFall
                << std::defaultfloat << " " << tOn << " " << period << ")\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
        circuit << "S1 q1_drain sw_node pwm_ctrl 0 SW1\n\n";

        // Low-side switch (synchronous rectification)
        circuit << "* Low-side Switch (Synchronous Rectifier)\n";
        circuit << "Vpwm_inv pwm_inv 0 PULSE(" << cfg.pwmHigh << " 0 0 "
                << std::scientific << cfg.pwmRise << " " << cfg.pwmFall
                << std::defaultfloat << " " << tOn << " " << period << ")\n";
        circuit << ".model SW2 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
        circuit << "S2 sw_node 0 pwm_inv 0 SW2\n\n";

        // Coupled inductor (Primary = buck inductor)
        // Using proper buck topology: series inductor with parallel DC path resistor
        circuit << "* Coupled Inductor (Primary = Buck Inductor)\n";
        
        // Primary current sense - measure current entering the inductor
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense sw_node pri_in 0\n";
        
        // Primary inductor: from switch node to output (series inductor like buck converter)
        // This is the proper topology - inductor in series with output, not to ground
        circuit << "Lpri pri_in vpri_out " << std::scientific << magnetizingInductance << std::fixed << "\n\n";

        // Across-winding differential voltage probe for the primary.
        // The raw v(pri_in) is the SWITCH-node voltage (≈Vin when S1 ON,
        // ≈0 when S2 ON / synchronous freewheel) — NOT the actual
        // voltage across Lpri, which is V(pri_in)−V(vpri_out). Same
        // class of bug fixed in PushPull commit 4ffd3c28: without this
        // correction MagneticFilterAreaProduct's powerMean estimate
        // is far too small (it uses |V·I| over the saved waveform), the
        // AP filter degrades to a no-op, and CoreAdviser advances cores
        // that physically cannot fit any wire — the failure mode behind
        // TestMagneticAdviser.cpp:5008.
        //
        // Sign convention: Lpri's first SPICE node = dot, so dot is at
        // pri_in. The current sense Vpri_sense (sw_node → pri_in)
        // reports positive i when current flows INTO the dot. Defining
        // V_winding = V(pri_in) − V(vpri_out) (dot to non-dot) keeps
        // passive sign convention against this current — i.e.
        // avg(V·I) > 0 represents power flowing INTO the inductor.
        circuit << "* Across-winding differential voltage probe (passive sign)\n";
        circuit << "Bvpri_diff vpri_diff 0 V=V(pri_in)-V(vpri_out)\n\n";

        // Secondary windings
        // For flyback action, secondary windings referenced to ground
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double secondaryInductance = magnetizingInductance / (turnsRatios[secIdx] * turnsRatios[secIdx]);
            // Secondary winding: ground to sec_in
            circuit << "Lsec" << secIdx << " 0 sec" << secIdx << "_in " << std::scientific << secondaryInductance << std::fixed << "\n";
        }

        // Couple primary to each secondary with high coupling
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            // K=1 (perfect coupling) matches the analytical model's
            // ideal-transformer assumption. With K<1 the induced
            // secondary voltage (K · V_pri / n) is below the nominal
            // secondary rail by a factor (1−K)/n, so for low-voltage
            // secondaries (e.g. design 2: 5/9 = 0.55 V nominal) the
            // secondary diode NEVER reaches forward bias against its
            // Cout_sec IC=Vout_sec_nom initial condition; the Cout
            // then discharges through Rload_sec to ~0 V and breaks
            // the entire flybuck operation in SPICE. K=1 also matches
            // Flyback's coupling (already validated DAB-quality with
            // K=1). The trade-off is that with K=1 the individual
            // winding currents can step at switching events (only
            // total ampere-turns is constrained continuous, not each
            // winding) — see piecewise analytical model below for
            // how this is handled.
            circuit << "Kpri_sec" << secIdx << " Lpri Lsec" << secIdx << " 1\n";
        }
        // Couple secondaries to each other. Must also be K=1 for matrix
        // consistency: with Kpri_secA = Kpri_secB = 1, the K-matrix forces
        // Ksec_secA_secB = 1 too (any lesser K violates the Cauchy-Schwarz
        // bound on the inductance matrix and makes the coupled-inductor
        // ngspice element non-positive-definite → simulation aborts with
        // "No output file generated by ngspice").
        for (size_t i = 0; i < numSecondaries; ++i) {
            for (size_t j = i + 1; j < numSecondaries; ++j) {
                circuit << "Ksec" << i << "_" << j << " Lsec" << i << " Lsec" << j << " 1\n";
            }
        }
        circuit << "\n";

        // Diode model
        circuit << "* Diode model\n";
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::defaultfloat;
        if (!cfg.diodeExtra.empty()) circuit << " " << cfg.diodeExtra;
        circuit << ")\n\n";

        // Primary output (buck output) - directly from inductor output
        // In this topology, vpri_out is the output node (after the inductor)
        circuit << "* Primary Output Stage (Buck - non-isolated)\n";
        double primaryLoadResistance = primaryOutputVoltage / primaryOutputCurrent;
        circuit << "Cpri vpri_out 0 100u IC=" << primaryOutputVoltage << "\n";
        circuit << "Rload_pri vpri_out 0 " << primaryLoadResistance << "\n\n";

        // Secondary output stages (isolated outputs)
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << "* Secondary " << secIdx << " output stage (isolated)\n";
            // Add series resistance between winding and diode (models winding resistance)
            circuit << "Rsec" << secIdx << " sec" << secIdx << "_in sec" << secIdx << "_node 0.01\n";
            // Secondary rectifier: diode from secondary winding node to rectified node
            // Anode at secX_node (winding side), cathode at secX_rect (output side)
            circuit << "Dsec" << secIdx << " sec" << secIdx << "_node sec" << secIdx << "_rect DIDEAL\n";
            circuit << "Vsec_sense" << secIdx << " sec" << secIdx << "_rect vout" << secIdx << " 0\n";

            double outputVoltage = opPoint.get_output_voltages()[secIdx + 1];
            double outputCurrent = opPoint.get_output_currents()[secIdx + 1];
            double loadResistance = outputVoltage / outputCurrent;

            // Cout startup IC. The secondary diode forward-biases only when
            // V(sec_in)_peak ≈ Vpri_out/n exceeds V(vout) + Vd. When the
            // analytical peak is strictly below the spec'd outputVoltage,
            // even an ideal diode (Vd→0) cannot forward-bias against an
            // IC=outputVoltage pre-charge, so the secondary loop carries only
            // numerical noise — the failure mode documented in
            // Heaviside/docs/MKF-HANDOFF.md §2 (FFT amplitudes ~1000× too
            // small → CoilMesher empty-harmonic). In that case seed Cout at
            // 0 V so the diode conducts during transient. For marginal specs
            // where peak ≈ outputVoltage (e.g. LM5160 24→12V, n=1) keep
            // IC=outputVoltage: the cap stays pre-charged through the short
            // PtP sim window and downstream gates are unaffected.
            double secondaryVoltagePeak = primaryOutputVoltage / turnsRatios[secIdx];
            double coutIC = outputVoltage;
            if (secondaryVoltagePeak < outputVoltage) {
                coutIC = 0.0;
            }
            circuit << "Cout" << secIdx << " vout" << secIdx << " 0 100u IC=" << coutIC << "\n";
            circuit << "Rload" << secIdx << " vout" << secIdx << " 0 " << loadResistance << "\n\n";
        }

        // Transient Analysis with UIC (Use Initial Conditions)
        // UIC skips the DC operating point calculation and uses the specified initial conditions
        // This is necessary for circuits with inductors that don't have DC paths to ground
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << " UIC\n\n";

        // Save signals
        // For Flybuck topology:
        // - pri_in: inductor input (after current sense)
        // - vpri_out: primary output (inductor output / non-isolated output)
        // - voutX: secondary outputs (isolated outputs)
        // - secX_in: secondary winding inputs
        // - secX_rect: secondary rectified node (before output cap)
        circuit << "* Output signals\n";
        circuit << ".save v(sw_node) v(pri_in) v(vpri_out) v(vpri_diff) i(Vpri_sense) v(vin_dc) i(Vin) i(Vq1_sense)";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            circuit << " v(sec" << secIdx << "_in) v(sec" << secIdx << "_rect) i(Vsec_sense" << secIdx << ") v(vout" << secIdx << ")";
        }
        circuit << "\n\n";

        // Options (matching other converters for convergence)
        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::defaultfloat
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        circuit << ".ic v(vpri_out)=" << primaryOutputVoltage << "\n";
        for (size_t secIdx = 0; secIdx < numSecondaries; ++secIdx) {
            double outputVoltage = opPoint.get_output_voltages()[secIdx + 1];
            double secondaryVoltagePeak = primaryOutputVoltage / turnsRatios[secIdx];
            double coutIC = (secondaryVoltagePeak < outputVoltage) ? 0.0 : outputVoltage;
            circuit << ".ic v(vout" << secIdx << ")=" << coutIC << "\n";
        }
        circuit << "\n";

        circuit << ".end\n";

        return circuit.str();
    }

    std::vector<OperatingPoint> IsolatedBuck::simulate_and_extract_operating_points(
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
                auto ibOpPoint = get_operating_points()[opIndex];

                std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, inputVoltageIndex, opIndex);

                double switchingFrequency = ibOpPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);

                if (!simResult.success) {
                    throw std::runtime_error("ngspice simulation failed: " + simResult.errorMessage);
                }

                // Define waveform name mapping
                NgspiceRunner::WaveformNameMapping waveformMapping;
                
                // Primary winding — across-winding differential voltage
                // V(pri_in)−V(vpri_out) (passive sign w.r.t. Vpri_sense).
                // Was previously v(pri_in) which is the SWITCH-node
                // voltage and led to ~3-orders-of-magnitude AP filter
                // underestimate (3606-class bug; see SPICE comment).
                waveformMapping.push_back({{"voltage", "vpri_diff"}, {"current", "vpri_sense#branch"}});
                
                // Secondary windings: probe the across-Lsec voltage. Each
                // Lsec<N> in the netlist is `Lsec<N> 0 sec<N>_in`, so its
                // far terminal is GND and v(sec<N>_in) IS the differential
                // winding voltage (per §5.0). Previously this used
                // sec<N>_rect — the post-diode/post-Cout rectified node —
                // which is downstream of the rectifier and a §5.0 violation
                // (~98 % normAvg violation pre-fix).
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
                    ibOpPoint.get_ambient_temperature(),
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

} // namespace OpenMagnetics
