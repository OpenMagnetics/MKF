#include "converter_models/Buck.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/NgspiceRunner.h"
#include "support/Utils.h"
#include <cfloat>
#include <algorithm>
#include <map>
#include <sstream>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    double Buck::calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency) {
        auto dutyCycle = (outputVoltage + diodeVoltageDrop) / ((inputVoltage + diodeVoltageDrop) * efficiency);
        if (dutyCycle >= 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Duty cycle must be smaller than 1");
        }
        // Explicit maxD gate (no silent clamp). 1 % tolerance absorbs
        // design-requirements rounding. Same pattern as Flyback and the
        // forward family.
        if (maximumDutyCycle.has_value()) {
            const double dMax = maximumDutyCycle.value();
            constexpr double dutyTolerance = 0.01;
            if (dutyCycle > dMax * (1.0 + dutyTolerance)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Buck: required dutyCycle " + std::to_string(dutyCycle) +
                    " exceeds maximumDutyCycle " + std::to_string(dMax) +
                    " at Vin=" + std::to_string(inputVoltage) +
                    ", Vout=" + std::to_string(outputVoltage) +
                    ". Raise Vin_min, lower Vout, or relax maximumDutyCycle.");
            }
        }
        return dutyCycle;
    }

    Buck::Buck(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    AdvancedBuck::AdvancedBuck(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    OperatingPoint Buck::process_operating_points_for_input_voltage(double inputVoltage, const BaseOperatingPoint& outputOperatingPoint, double inductance) {
        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double outputVoltage = outputOperatingPoint.get_output_voltages()[0];
        double outputCurrent = outputOperatingPoint.get_output_currents()[0];
        double diodeVoltageDrop = get_diode_voltage_drop();
        double efficiency = require_input(get_efficiency(), "Buck", "efficiency");

        auto dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, diodeVoltageDrop, efficiency);

        auto tOn = dutyCycle / switchingFrequency;
        auto primaryCurrentPeakToPeak = (inputVoltage - outputVoltage) * tOn / inductance;
        auto minimumCurrent = outputCurrent - primaryCurrentPeakToPeak / 2;
        auto primaryVoltageMinimum = -outputVoltage -diodeVoltageDrop;
        auto primaryVoltageMaximum = inputVoltage - outputVoltage;
        auto primaryVoltagePeaktoPeak = primaryVoltageMaximum - primaryVoltageMinimum;

        // ---- Diagnostic capture (analytical, last-call) ----
        // Populated for whichever Vin this method is called with; the
        // outer process_operating_points() iterates Nom/Min/Max so the
        // values reflect the LAST call. Tests that need a specific Vin
        // should call this method directly.
        const double period = 1.0 / switchingFrequency;
        lastDutyCycle = dutyCycle;
        lastInductorAverageCurrent = outputCurrent;
        lastInductorPeakToPeak = primaryCurrentPeakToPeak;
        lastIsCcm = (minimumCurrent >= 0);
        if (lastIsCcm) {
            lastPeakInductorCurrent = outputCurrent + primaryCurrentPeakToPeak / 2.0;
            lastConductionRatio = 1.0;
        }

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            if (minimumCurrent < 0) {
                tOn = sqrt(2 * outputCurrent * inductance * (outputVoltage + diodeVoltageDrop) / (switchingFrequency * (inputVoltage - outputVoltage) * (inputVoltage + diodeVoltageDrop)));
                auto tOff = tOn * ((inputVoltage + diodeVoltageDrop) / (outputVoltage + diodeVoltageDrop) - 1);
                auto deadTime = 1.0 / switchingFrequency - tOn - tOff;
                primaryCurrentPeakToPeak = (inputVoltage - outputVoltage) * tOn / inductance;
                outputCurrent = primaryCurrentPeakToPeak / 2;

                // DCM: the waveform timing must use the recomputed tOn, not the
                // CCM duty cycle (which is larger and at light load even pushed
                // the triangle peak past period - deadTime, producing a
                // non-monotonic time vector).
                double dcmDutyCycle = tOn * switchingFrequency;
                currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR_WITH_DEADTIME, primaryCurrentPeakToPeak, switchingFrequency, dcmDutyCycle, outputCurrent, deadTime);
                // Voltage: explicit levels and times. A single duty parameter
                // cannot express both the tOn transition and the physical level
                // split in DCM (RECTANGULAR_WITH_DEADTIME would violate
                // volt-second balance over the conduction interval).
                voltageWaveform.set_data(std::vector<double>{primaryVoltageMaximum, primaryVoltageMaximum, primaryVoltageMinimum, primaryVoltageMinimum, 0, 0});
                voltageWaveform.set_time(std::vector<double>{0, tOn, tOn, tOn + tOff, tOn + tOff, 1.0 / switchingFrequency});
                voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);

                // DCM diagnostic capture: peak current is the full ΔIL (waveform
                // returns to zero each cycle). Conduction ratio < 1.
                lastInductorPeakToPeak = primaryCurrentPeakToPeak;
                lastPeakInductorCurrent = primaryCurrentPeakToPeak;
                lastInductorAverageCurrent = outputCurrent;  // = ΔIL/2 per area-balance
                lastConductionRatio = (tOn + tOff) / period;
            }
            else {
                currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, outputCurrent, 0);
                voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, 0);
            }

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        // Per-OP diagnostic snapshot for the wizard's multi-column table.
        perOpDutyCycle.push_back(lastDutyCycle);
        perOpInductorAverageCurrent.push_back(lastInductorAverageCurrent);
        perOpInductorPeakToPeak.push_back(lastInductorPeakToPeak);
        perOpPeakInductorCurrent.push_back(lastPeakInductorCurrent);
        perOpIsCcm.push_back(lastIsCcm);
        perOpConductionRatio.push_back(lastConductionRatio);

        return operatingPoint;
    }

    bool Buck::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
        }
        if (!get_input_voltage().get_nominal() && !get_input_voltage().get_maximum() && !get_input_voltage().get_minimum()) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA, "No input voltage introduced");
        }

        return true;
    }

    DesignRequirements Buck::process_design_requirements() {
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Missing both current ripple ratio and maximum switch current");
        }

        double maximumCurrentRiple = 0;
        if (get_current_ripple_ratio()) {
            double currentRippleRatio = get_current_ripple_ratio().value();
            double maximumOutputCurrent = 0;

            for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
                auto buckOperatingPoint = get_operating_points()[buckOperatingPointIndex];
                maximumOutputCurrent = std::max(maximumOutputCurrent, buckOperatingPoint.get_output_currents()[0]);
            }

            maximumCurrentRiple = currentRippleRatio * maximumOutputCurrent;
        }

        if (get_maximum_switch_current()) {
            double maximumSwitchCurrent = get_maximum_switch_current().value();
            double maximumOutputCurrent = 0;

            for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
                auto buckOperatingPoint = get_operating_points()[buckOperatingPointIndex];
                maximumOutputCurrent = std::max(maximumOutputCurrent, buckOperatingPoint.get_output_currents()[0]);
            }

            maximumCurrentRiple = (maximumSwitchCurrent - maximumOutputCurrent) * 2;
        }

        double maximumNeededInductance = 0;
        for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
            auto buckOperatingPoint = get_operating_points()[buckOperatingPointIndex];
            auto switchingFrequency = buckOperatingPoint.get_switching_frequency();
            auto outputVoltage = buckOperatingPoint.get_output_voltages()[0];
            auto desiredInductance = outputVoltage * (maximumInputVoltage - outputVoltage) / (maximumCurrentRiple * switchingFrequency * maximumInputVoltage);
            maximumNeededInductance = std::max(maximumNeededInductance, desiredInductance);
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::BUCK_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Buck::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        // Clear per-OP diagnostic vectors so the wizard table reflects this run only.
        perOpName.clear();
        perOpDutyCycle.clear();
        perOpInductorAverageCurrent.clear();
        perOpInductorPeakToPeak.clear();
        perOpPeakInductorCurrent.clear();
        perOpIsCcm.clear();
        perOpConductionRatio.clear();

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
                std::string opName = inputVoltagesNames[inputVoltageIndex];
                if (get_operating_points().size() > 1) {
                    opName += " · OP" + std::to_string(buckOperatingPointIndex);
                }
                perOpName.push_back(opName);
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[buckOperatingPointIndex], magnetizingInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(buckOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> Buck::process_operating_points(OpenMagnetics::Magnetic magnetic) {
        run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    DesignRequirements AdvancedBuck::process_design_requirements() {
        // Advanced workflow: the user already knows the L they want;
        // set it as nominal in DesignRequirements so the CoreAdviser
        // targets it instead of the parent Buck's ripple-floor minimum.
        DesignRequirements designRequirements;
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(get_desired_inductance(), 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::BUCK_CONVERTER);
        return designRequirements;
    }

    Inputs AdvancedBuck::process() {
        Buck::run_checks(_assertErrors);

        Inputs inputs;

        double maximumNeededInductance = get_desired_inductance();

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;


        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        DesignRequirements designRequirements;

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::BUCK_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t buckOperatingPointIndex = 0; buckOperatingPointIndex < get_operating_points().size(); ++buckOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[buckOperatingPointIndex], maximumNeededInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(buckOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

    std::string Buck::generate_ngspice_circuit(
        double inductance,
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
        auto opPoint = get_operating_points()[operatingPointIndex];
        
        double outputVoltage = opPoint.get_output_voltages()[0];
        double outputCurrent = opPoint.get_output_currents()[0];
        double switchingFrequency = opPoint.get_switching_frequency();
        double diodeVoltageDrop = get_diode_voltage_drop();
        double efficiency = require_input(get_efficiency(), "Buck", "efficiency");
        
        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, diodeVoltageDrop, efficiency);

        // SPICE configuration knobs (registered defaults, see
        // Topology.cpp::spice_simulation_defaults()). Per-instance
        // overrides apply via Topology::set_spice_config().
        const SpiceSimulationConfig cfg = spice_config();

        // Build netlist
        std::ostringstream circuit;
        double period = 1.0 / switchingFrequency;
        double tOn = period * dutyCycle;
        
        // Simulation: run 10x the extraction periods for settling, then extract the last N periods
        int periodsToExtract = get_num_periods_to_extract();
    int numSteadyStatePeriods = get_num_steady_state_periods();
    const int numPeriodsTotal = numSteadyStatePeriods + periodsToExtract;
    double simTime = numPeriodsTotal * period;
    double startTime = numSteadyStatePeriods * period;
        double stepTime = period / cfg.samplesPerPeriod;
        
        circuit << "* Buck Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vout=" << outputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* L=" << (inductance*1e6) << "uH, Iout=" << outputCurrent << "A\n\n";
        
        // DC Input.  The input source itself does NOT measure source-side
        // current correctly via i(Vin) once we want it under our control,
        // so we insert a 0V sense source between the source and the rest
        // of the circuit.  i(Vin_sense) is positive when current flows
        // from the source into the converter (load convention), which is
        // exactly what `wf.set_input_current` expects below.
        // Buck-specific: the inductor is on the OUTPUT side (after the
        // switch), so the inductor sense (vl_sense#branch) measures the
        // OUTPUT inductor current, NOT the source-side input current.
        // Source-side current only flows during the switch ON window.
        circuit << "* DC Input\n";
        circuit << "Vin vin_src 0 " << inputVoltage << "\n";
        circuit << "Vin_sense vin_src vin_dc 0\n\n";

        // §8a.5 input-current sense: a zero-V source upstream of S1's
        // drain isolates the true switch-channel current. Buck has no
        // return path from elsewhere into vin_dc (D1 freewheels to GND,
        // not back to vin_dc), so i(Vin_sense) and the switch current
        // are equal in the ideal case — but the snubber is wired
        // vin_dc <-> sw and DOES carry small charge/discharge current
        // through Vin_sense. Per the §8a.5 template, the snubber and
        // any future bypass branches stay on vin_dc and are excluded
        // from the converter input-draw measurement by routing S1 from
        // a separate q1_drain node and inserting Vq1_sense between
        // vin_dc and q1_drain. Sign: + terminal is vin_dc, - is
        // q1_drain, so positive i(Vq1_sense) means current flowing OUT
        // of vin_dc INTO S1 (non-negative in CCM/DCM).
        circuit << "* §8a.5 input-current sense\n";
        circuit << "Vq1_sense vin_dc q1_drain 0\n\n";

        // PWM High-side Switch (ideal). PULSE timing emitted in scientific
        // notation so sub-µs values aren't silently rounded by std::fixed
        // (precision 6) to the wrong duty cycle / Fs (cf. Boost.cpp PWM
        // precision bug fix; same pattern previously latent here).
        circuit << "* PWM High-side Switch\n";
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                << tOn << " " << period << std::fixed
                << ")\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
        circuit << "S1 q1_drain sw pwm_ctrl 0 SW1\n";
        // Snubber RC across the switch — absorbs di/dt at hard-switching
        // events, improves convergence in CCM/DCM transitions. Stays on
        // vin_dc (NOT q1_drain) so its current is excluded from
        // i(Vq1_sense) per §8a.5.
        circuit << "Rsnub_s1 vin_dc sw " << cfg.snubR << "\n"
                << "Csnub_s1 vin_dc sw " << std::scientific << cfg.snubC << std::fixed << "\n\n";
        
        // Freewheeling Diode (low-side)
        circuit << "* Freewheeling Diode\n";
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::fixed << ")\n";
        circuit << "D1 0 sw DIDEAL\n\n";
        
        // Inductor with current sense
        circuit << "* Inductor with current sense\n";
        circuit << "Vl_sense sw l_in 0\n";
        circuit << "L1 l_in vout " << std::scientific << inductance << std::fixed << "\n";
        // §5.0: across-winding differential probe. The inductor's "other"
        // terminal (vout) is NOT GND. v(sw) avg ≈ D·Vin = Vout, which
        // collapses the bipolar swing in the winding-port stream.
        // Differential V(l_in)-V(vout) gives the correct ±Vin/Vout swing
        // with avg = 0. Vl_sense = (sw → l_in) so i(Vl_sense) > 0 flows
        // INTO L1's l_in terminal (= dot). V_winding = V(l_in)-V(vout)
        // is +ve at the dot ⇒ avg(V·I) > 0 (passive convention).
        circuit << "Bvpri_diff vpri_diff 0 V=V(l_in)-V(vout)\n\n";
        
        // Output capacitor and load
        if (outputCurrent <= 0.0) {
            throw std::invalid_argument(
                "Buck::generate_ngspice_circuit: outputCurrent must be > 0 "
                "(received " + std::to_string(outputCurrent) +
                "); cannot derive load resistance R_load = Vout / Iout.");
        }
        double loadResistance = outputVoltage / outputCurrent;
        circuit << "* Output Filter and Load\n";
        circuit << "Cout vout 0 " << std::scientific << cfg.outputCapacitance << std::fixed
                << " IC=" << (outputVoltage * cfg.outputCapInitialChargeFraction) << "\n";
        circuit << "Rload vout 0 " << loadResistance << "\n\n";
        
        // Transient Analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime << " " << startTime << std::fixed << "\n\n";
        
        // Save signals — split per §5.0:
        //   - Winding-port: v(vpri_diff) i(Vl_sense)
        //   - Converter-port: v(vin_dc) v(sw) v(vout) i(Vl_sense)
        circuit << "* Output signals\n";
        circuit << ".save v(vpri_diff) v(vin_dc) v(sw) v(vout) i(Vl_sense) i(Vin_sense) i(Vq1_sense)\n\n";
        
        // Solver options for convergence in switching circuits.
        // METHOD=GEAR + larger TRTOL handles the hard-switching event at
        // each PWM edge robustly; higher ITL4 lets the solver retry
        // CCM/DCM boundary transients.
        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::fixed
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        circuit << ".options METHOD=" << cfg.method << " TRTOL=" << cfg.trTol << "\n";
        circuit << ".ic v(vout)=" << outputVoltage << "\n\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> Buck::simulate_and_extract_operating_points(double inductance) {
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto buckOpPoint = get_operating_points()[opIndex];
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(inductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = buckOpPoint.get_switching_frequency();
                
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
                
                // §5.0: winding-port stream — differential voltage across
                // L1 terminals (l_in → vout), passive-convention current.
                NgspiceRunner::WaveformNameMapping waveformMapping = {
                    {{"voltage", "vpri_diff"}, {"current", "vl_sense#branch"}}
                };
                
                std::vector<std::string> windingNames = {"Primary"};
                std::vector<bool> flipCurrentSign = {false};
                
                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult,
                    waveformMapping,
                    switchingFrequency,
                    windingNames,
                    buckOpPoint.get_ambient_temperature(),
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
        
        return operatingPoints;
    }

    std::vector<ConverterWaveforms> Buck::simulate_and_extract_topology_waveforms(double inductance, size_t numberOfPeriods) {
    
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
            
            std::string netlist = generate_ngspice_circuit(inductance, inputVoltageIndex, opIndex);
            double switchingFrequency = opPoint.get_switching_frequency();
            const double outputVoltage = opPoint.get_output_voltages()[0];
            const double outputCurrent = opPoint.get_output_currents()[0];
            
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
            //     NOT v(sw) — the switch-node is a square wave between
            //     0 and Vin within a switching period and would be a
            //     wrong "Input V" reading; bug A doesn't apply here
            //     because we were already on vin_dc, but the §8a.5
            //     audit confirms the choice.
            //   input_current = i(Vq1_sense) — zero-V sense source
            //     upstream of S1's drain. Captures only the on-time
            //     switch-channel current. i(Vin_sense) (the prior probe)
            //     ALSO includes the snubber Rsnub/Csnub charge-discharge
            //     current (the snubber stays on vin_dc by design), so
            //     it overstates the converter input draw by the snubber
            //     contribution at each switching edge. Sign: + terminal
            //     of Vq1_sense is vin_dc, so positive i(Vq1_sense) is
            //     current flowing OUT of vin_dc INTO S1 (non-negative
            //     in CCM/DCM).
            //   Bug C: n/a — Buck is non-isolated, only one winding
            //     (Lpri = L1 wired sw->l_in via Vl_sense -> vout). The
            //     across-winding probe Bvpri_diff = v(l_in)-v(vout)
            //     already provides the true winding-only differential
            //     voltage for the winding-port stream.
            wf.set_input_voltage(getWaveform("vin_dc"));

            Waveform iIn = getWaveform("vq1_sense#branch");
            if (iIn.get_data().empty()) {
                throw std::runtime_error(
                    "Buck simulate_and_extract_topology_waveforms: "
                    "missing i(Vq1_sense) — netlist/save mismatch");
            }

            // ngspice SW1 model produces narrow di/dt spikes (~10^5 A)
            // around switching edges that are a numerical artefact of
            // the discontinuous RON/ROFF transition, not a physical
            // bound. Clamp samples to ±2·max(|i(Vl_sense)|) — the
            // inductor current sets the actual conduction-current
            // scale (during ON it IS the switch current); this
            // preserves all real samples while discarding the artefact
            // spikes that would poison the DC mean used downstream by
            // §5.1.
            //
            // numerical guard against ngspice SW-model di/dt spikes
            // (~10^5 A), not a physical bound.
            {
                auto& d = iIn.get_mutable_data();
                Waveform iL = getWaveform("vl_sense#branch");
                double iMax = 0.0;
                for (double v : iL.get_data()) iMax = std::max(iMax, std::abs(v));
                if (iMax > 0.0) {
                    const double clamp = 2.0 * iMax;
                    for (auto& v : d) {
                        if (v >  clamp) v =  clamp;
                        if (v < -clamp) v = -clamp;
                    }
                }
            }
            wf.set_input_current(iIn);

            wf.get_mutable_output_voltages().push_back(getWaveform("vout"));
            // Reconstruct Iout(t) = Vout(t) / Rload (DC by design).
            Waveform ioutWf = getWaveform("vout");
            auto& ioutData = ioutWf.get_mutable_data();
            const double rLoad = outputVoltage / outputCurrent;
            for (auto& v : ioutData) v = v / rLoad;
            wf.get_mutable_output_currents().push_back(ioutWf);
            
            results.push_back(wf);
        }
    }
    
    // Restore original value
    set_num_periods_to_extract(originalNumPeriodsToExtract);
    
    return results;
}

} // namespace OpenMagnetics
