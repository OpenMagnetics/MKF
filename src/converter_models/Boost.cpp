#include "converter_models/Boost.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/NgspiceRunner.h"
#include "support/Utils.h"
#include <cfloat>
#include <algorithm>
#include <map>
#include <sstream>
#include "support/Exceptions.h"
#include "converter_models/Topology.h"

namespace OpenMagnetics {

    double Boost::calculate_duty_cycle(double inputVoltage, double outputVoltage, double diodeVoltageDrop, double efficiency) {
        auto dutyCycle = 1 - inputVoltage * efficiency / (outputVoltage + diodeVoltageDrop);
        if (dutyCycle >= 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Duty cycle must be smaller than 1");
        }
        // Explicit maxD gate (no silent clamp). 1 % tolerance absorbs
        // design-requirements rounding. Same pattern as Flyback, forward
        // family, Buck.
        if (maximumDutyCycle.has_value()) {
            const double dMax = maximumDutyCycle.value();
            constexpr double dutyTolerance = 0.01;
            if (dutyCycle > dMax * (1.0 + dutyTolerance)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Boost: required dutyCycle " + std::to_string(dutyCycle) +
                    " exceeds maximumDutyCycle " + std::to_string(dMax) +
                    " at Vin=" + std::to_string(inputVoltage) +
                    ", Vout=" + std::to_string(outputVoltage) +
                    ". Raise Vin_min, lower Vout, or relax maximumDutyCycle.");
            }
        }
        return dutyCycle;
    }

    Boost::Boost(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    AdvancedBoost::AdvancedBoost(const json& j) {
        json migrated = j;
        migrate_operating_point_json(migrated);
        from_json(migrated, *this);
    }

    OperatingPoint Boost::process_operating_points_for_input_voltage(double inputVoltage, const BaseOperatingPoint& outputOperatingPoint, double inductance) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.get_switching_frequency();
        double outputVoltage = outputOperatingPoint.get_output_voltages()[0];
        double outputCurrent = outputOperatingPoint.get_output_currents()[0];
        double diodeVoltageDrop = get_diode_voltage_drop();
        double efficiency = 1;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }


        auto dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, diodeVoltageDrop, efficiency);

        auto tOn = dutyCycle / switchingFrequency;
        auto primaryCurrentPeakToPeak = inputVoltage * tOn / inductance;
        auto primaryCurrentAverage = outputCurrent * (outputVoltage + diodeVoltageDrop) / inputVoltage;
        auto primaryCurrentMinimum = primaryCurrentAverage - primaryCurrentPeakToPeak / 2;

        auto primaryVoltageMinimum = inputVoltage - outputVoltage - diodeVoltageDrop;
        auto primaryVoltageMaximum = inputVoltage;
        auto primaryVoltagePeaktoPeak = primaryVoltageMaximum - primaryVoltageMinimum;

        // ---- Diagnostics (CCM defaults; DCM overwrites below) ----
        lastDutyCycle = dutyCycle;
        lastInductorAverageCurrent = primaryCurrentAverage;
        lastInductorPeakToPeak = primaryCurrentPeakToPeak;
        lastPeakInductorCurrent = primaryCurrentAverage + primaryCurrentPeakToPeak / 2;
        lastIsCcm = (primaryCurrentMinimum >= 0);
        lastConductionRatio = 1.0;

        // Primary
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            if (primaryCurrentMinimum < 0) {
                tOn = sqrt(2 * outputCurrent * inductance * (outputVoltage + diodeVoltageDrop - inputVoltage) / (switchingFrequency * pow(inputVoltage, 2)));
                auto tOff = tOn * ((outputVoltage + diodeVoltageDrop) / (outputVoltage + diodeVoltageDrop - inputVoltage) - 1);
                auto deadTime = 1.0 / switchingFrequency - tOn - tOff;
                outputCurrent = primaryCurrentPeakToPeak / 2;

                currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR_WITH_DEADTIME, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, outputCurrent, deadTime);
                voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);

                // ---- DCM diagnostics ----
                // Peak current in DCM = ΔIL_pp (current ramps from 0 to peak
                // during t_on, then back to 0 during t_off, then 0 for t_dt).
                lastIsCcm = false;
                lastInductorPeakToPeak = primaryCurrentPeakToPeak;
                lastPeakInductorCurrent = primaryCurrentPeakToPeak;
                lastInductorAverageCurrent = outputCurrent;  // local var was reassigned to ΔIL/2
                lastConductionRatio = (tOn + tOff) * switchingFrequency;
            }
            else {
                currentWaveform = Inputs::create_waveform(WaveformLabel::TRIANGULAR, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, primaryCurrentAverage, 0);
                voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, 0);
            }

            auto excitation = complete_excitation(currentWaveform, voltageWaveform, switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    bool Boost::run_checks(bool assert) {
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

    DesignRequirements Boost::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);
        double efficiency = 1;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }

        if (!get_current_ripple_ratio() && !get_maximum_switch_current()) {
            throw std::invalid_argument("Missing both current ripple ratio and maximum switch current");
        }

        double maximumCurrentRiple = 0;
        if (get_current_ripple_ratio()) {
            double currentRippleRatio = get_current_ripple_ratio().value();
            double maximumOutputCurrent = 0;

            for (size_t boostOperatingPointIndex = 0; boostOperatingPointIndex < get_operating_points().size(); ++boostOperatingPointIndex) {
                auto boostOperatingPoint = get_operating_points()[boostOperatingPointIndex];
                maximumOutputCurrent = std::max(maximumOutputCurrent, boostOperatingPoint.get_output_currents()[0]);
            }

            maximumCurrentRiple = currentRippleRatio * maximumOutputCurrent;
        }

        if (get_maximum_switch_current()) {
            double maximumSwitchCurrent = get_maximum_switch_current().value();

            for (size_t boostOperatingPointIndex = 0; boostOperatingPointIndex < get_operating_points().size(); ++boostOperatingPointIndex) {
                auto boostOperatingPoint = get_operating_points()[boostOperatingPointIndex];
                auto dutyCycle = calculate_duty_cycle(minimumInputVoltage, boostOperatingPoint.get_output_voltages()[0], get_diode_voltage_drop(), efficiency);
                auto currentRiple = (maximumSwitchCurrent - boostOperatingPoint.get_output_currents()[0] / (1.0 - dutyCycle)) * 2;
                maximumCurrentRiple = std::max(maximumCurrentRiple, currentRiple);
            }

        }

        double maximumNeededInductance = 0;
        for (size_t boostOperatingPointIndex = 0; boostOperatingPointIndex < get_operating_points().size(); ++boostOperatingPointIndex) {
            auto boostOperatingPoint = get_operating_points()[boostOperatingPointIndex];
            auto switchingFrequency = boostOperatingPoint.get_switching_frequency();
            auto outputVoltage = boostOperatingPoint.get_output_voltages()[0];
            auto desiredInductance = maximumInputVoltage * (outputVoltage - maximumInputVoltage) / (maximumCurrentRiple * switchingFrequency * outputVoltage);
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
        designRequirements.set_topology(Topologies::BOOST_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Boost::process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t boostOperatingPointIndex = 0; boostOperatingPointIndex < get_operating_points().size(); ++boostOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[boostOperatingPointIndex], magnetizingInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(boostOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> Boost::process_operating_points(OpenMagnetics::Magnetic magnetic) {
        run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();
        
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    Inputs AdvancedBoost::process() {
        Boost::run_checks(_assertErrors);

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
        designRequirements.set_topology(Topologies::BOOST_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t boostOperatingPointIndex = 0; boostOperatingPointIndex < get_operating_points().size(); ++boostOperatingPointIndex) {
                auto operatingPoint = process_operating_points_for_input_voltage(inputVoltage, get_operating_points()[boostOperatingPointIndex], maximumNeededInductance);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(boostOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

    std::string Boost::generate_ngspice_circuit(
        double inductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames_;
    Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames_);
        
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
        double efficiency = 1.0;
        if (get_efficiency()) {
            efficiency = get_efficiency().value();
        }
        
        double dutyCycle = calculate_duty_cycle(inputVoltage, outputVoltage, diodeVoltageDrop, efficiency);

        // Resolve SPICE configuration (per-instance override, else registered default).
        // See SpiceSimulationConfig in Topology.h for field semantics.
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
        
        circuit << "* Boost Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vout=" << outputVoltage << "V, f=" << (switchingFrequency/1e3) << "kHz, D=" << (dutyCycle*100) << " pct\n";
        circuit << "* L=" << (inductance*1e6) << "uH, Iout=" << outputCurrent << "A\n\n";
        
        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";
        
        // Inductor with current sense (between input and switch node)
        circuit << "* Inductor with current sense\n";
        circuit << "Vl_sense vin_dc l_in 0\n";
        circuit << "L1 l_in sw " << std::scientific << inductance << std::fixed << "\n";
        // §5.0: across-winding differential probe. The inductor's "other"
        // terminal (sw) is NOT GND — it bounces between 0 (S1 closed) and
        // Vout (S1 open via D1). v(l_in)=Vin is constant DC, so saving it
        // as the winding voltage collapses the bipolar swing entirely
        // (1.00 → 0.00 normalised volt-second). The differential probe
        // V(l_in)-V(sw) gives the correct ±Vin-Vout swing with avg=0.
        // Sign matches passive convention against i(Vl_sense): current
        // flows vin_dc → l_in (INTO the dot), V_winding = V(l_in)-V(sw)
        // positive at the dot ⇒ avg(V·I) > 0.
        circuit << "Bvpri_diff vpri_diff 0 V=V(l_in)-V(sw)\n\n";
        
        // PWM Low-side Switch (ideal)
        circuit << "* PWM Low-side Switch\n";
        // tOn and period emitted in scientific notation. With the default
        // std::fixed / precision 6 active in the surrounding stream, sub-µs
        // values get rounded to 6 decimal places (e.g. 2.5e-6 → "0.000003"),
        // which silently shifts both Fs and duty cycle. Scientific format
        // preserves full mantissa so the simulator runs at the requested OP.
        circuit << "Vpwm pwm_ctrl 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                << tOn << " " << period << std::fixed
                << ")\n";
        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH << "\n";
        circuit << "S1 sw 0 pwm_ctrl 0 SW1\n";
        // Snubber RC across the switch — absorbs di/dt at hard-switching
        // events and improves ngspice convergence in CCM/DCM transitions
        // and at high duty (Vout >> Vin). Defaults R=100 Ω, C=100 pF
        // (τ = 10 ns); see SpiceSimulationConfig for override knobs.
        circuit << "Rsnub_s1 sw 0 " << cfg.snubR << "\n"
                << "Csnub_s1 sw 0 " << std::scientific << cfg.snubC << std::fixed << "\n\n";

        // Output Diode
        circuit << "* Output Diode\n";
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::fixed << ")\n";
        circuit << "D1 sw vout DIDEAL\n\n";

        // Output capacitor and load
        if (outputCurrent <= 0.0) {
            throw std::invalid_argument(
                "Boost::generate_ngspice_circuit: outputCurrent must be > 0 "
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

        // Save signals — split into two streams per §5.0:
        //   - Winding-port stream: v(vpri_diff) i(Vl_sense)
        //     (consumed by simulate_and_extract_operating_points →
        //      excitations_per_winding)
        //   - Converter-port stream: v(vin_dc) v(sw) v(vout) i(Vl_sense)
        //     (consumed by simulate_and_extract_topology_waveforms →
        //      ConverterWaveforms; output rails are intentionally DC)
        circuit << "* Output signals\n";
        circuit << ".save v(vpri_diff) v(vin_dc) v(sw) v(vout) i(Vl_sense)\n\n";

        // Solver options for convergence in switching circuits.
        // METHOD=GEAR + larger TRTOL handles the hard-switching event at
        // each PWM edge robustly; higher ITL4 lets the solver retry
        // CCM/DCM boundary transients. Pattern mirrors Dab.cpp.
        circuit << ".options RELTOL=" << cfg.relTol
                << " ABSTOL=" << std::scientific << cfg.absTol
                << " VNTOL=" << cfg.vnTol << std::fixed
                << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        circuit << ".options METHOD=" << cfg.method << " TRTOL=" << cfg.trTol << "\n";
        circuit << ".ic v(vout)=" << outputVoltage << "\n\n";
        
        circuit << ".end\n";
        
        return circuit.str();
    }
    
    std::vector<OperatingPoint> Boost::simulate_and_extract_operating_points(double inductance) {
        
        std::vector<OperatingPoint> operatingPoints;
        
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }
        
        // Get input voltages
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        Topology::collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        
        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto boostOpPoint = get_operating_points()[opIndex];
                
                // Generate circuit
                std::string netlist = generate_ngspice_circuit(inductance, inputVoltageIndex, opIndex);
                
                double switchingFrequency = boostOpPoint.get_switching_frequency();
                
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
                // the inductor terminals (vpri_diff = V(l_in)-V(sw)),
                // current via the in-line zero-volt sense source.
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
                    boostOpPoint.get_ambient_temperature(),
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
std::vector<ConverterWaveforms> Boost::simulate_and_extract_topology_waveforms(double inductance, size_t numberOfPeriods) {
    
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
        (void) inputVoltages[inputVoltageIndex];
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
                if (it != nameToIndex.end()) {
                    return simResult.waveforms[it->second];
                }
                return Waveform();
            };
            
            ConverterWaveforms wf;
            wf.set_switching_frequency(switchingFrequency);
            std::string name = inputVoltagesNames[inputVoltageIndex] + " input";
            if (get_operating_points().size() > 1) {
                name += " op. point " + std::to_string(opIndex);
            }
            wf.set_operating_point_name(name);
            
            // Topology / converter-port signals (§5.0): the converter's
            // electrical ports — Vin (DC source), Iin (inductor current
            // = source current for Boost), Vout (rectified output),
            // Iout (load current). Output current is reconstructed as
            // v(vout)/Rload; ngspice's diode/cap structure does not
            // expose Iout directly without an extra Vsense source on
            // the load.
            wf.set_input_voltage(getWaveform("vin_dc"));
            wf.set_input_current(getWaveform("vl_sense#branch"));
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
