#include "converter_models/PowerFactorCorrection.h"
#include "converter_models/PfcControllerDesign.h"
#include "converter_models/PfcControllerSubcircuits.h"
#include "processors/NgspiceRunner.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <numbers>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    PowerFactorCorrection::PowerFactorCorrection(const json& j) : Topology(j) {
        // Required fields per the MAS schema. Guard up front so users see
        // PFC-specific messages instead of a less specific j.at() throw.
        if (!j.contains("inputVoltage")) {
            throw std::runtime_error("PowerFactorCorrection: 'inputVoltage' is required");
        }
        if (!j.contains("outputVoltage")) {
            throw std::runtime_error("PowerFactorCorrection: 'outputVoltage' is required");
        }
        if (!j.contains("outputPower")) {
            throw std::runtime_error("PowerFactorCorrection: 'outputPower' is required");
        }
        if (!j.contains("switchingFrequency")) {
            throw std::runtime_error("PowerFactorCorrection: 'switchingFrequency' is required");
        }
        if (!j.contains("ambientTemperature")) {
            throw std::runtime_error("PowerFactorCorrection: 'ambientTemperature' is required");
        }

        // Delegate to the MAS-generated parser. It maps the optional
        // 'mode' field through the PfcModes enum (continuousConductionMode,
        // discontinuousConductionMode, criticalConductionMode, transitionMode).
        MAS::from_json(j, static_cast<MAS::PowerFactorCorrection&>(*this));

        // MKF-only knob (not in the MAS schema).
        if (j.contains("numberOfPeriods")) {
            _numberOfPeriods = j["numberOfPeriods"].get<int>();
        }
    }

    void PowerFactorCorrection::validate_topology_variant() const {
        const auto variant = get_topology_variant_or_default();
        switch (variant) {
            case PfcTopologyVariants::BOOST:
                // Single-phase boost — no companion-field requirements.
                return;

            case PfcTopologyVariants::TOTEM_POLE: {
                // Bridgeless totem-pole. CCM mode requires GaN/SiC switches
                // (no Si-MOSFET body-diode reverse-recovery tolerable).
                // DCM/CrCM/Transition modes are fine with any switch tech.
                auto m = MAS::PowerFactorCorrection::get_mode();
                if (m.has_value() && m.value() == PfcModes::CONTINUOUS_CONDUCTION_MODE
                    && !get_wide_bandgap_switch_or_default()) {
                    throw std::runtime_error(
                        "PowerFactorCorrection: totemPole topologyVariant in "
                        "CONTINUOUS_CONDUCTION_MODE requires wideBandgapSwitch=true "
                        "(GaN or SiC). Si MOSFETs are not viable due to body-diode "
                        "reverse-recovery losses (Erickson §17, ON Semi AND8016).");
                }
                return;
            }

            case PfcTopologyVariants::INTERLEAVED_BOOST: {
                const int64_t n = get_number_of_phases_or_default();
                if (n < 2) {
                    throw std::runtime_error(
                        "PowerFactorCorrection: interleavedBoost topologyVariant "
                        "requires numberOfPhases >= 2 (set via "
                        "set_number_of_phases() or the JSON field 'numberOfPhases'). "
                        "Got: " + std::to_string(n));
                }
                if (n > 3) {
                    throw std::runtime_error(
                        "PowerFactorCorrection: interleavedBoost numberOfPhases "
                        "must be 2 or 3 (got " + std::to_string(n) + "). "
                        "N>3 phase-shifted PFCs are not industry-standard.");
                }
                return;
            }

            case PfcTopologyVariants::VIENNA:
                throw std::runtime_error(
                    "PowerFactorCorrection: vienna is not implemented as a PFC "
                    "topologyVariant code branch. The Vienna rectifier is a "
                    "3-phase, 3-level topology with its own dedicated model — "
                    "see VIENNA_PLAN.md.");

            case PfcTopologyVariants::BRIDGELESS:
            case PfcTopologyVariants::SEMI_BRIDGELESS:
            case PfcTopologyVariants::BUCK:
            case PfcTopologyVariants::BUCK_BOOST:
            case PfcTopologyVariants::SEPIC:
            case PfcTopologyVariants::CUK:
                throw std::runtime_error(
                    "PowerFactorCorrection: topologyVariant is not yet "
                    "implemented in the MKF engineering layer. Currently "
                    "supported: BOOST, TOTEM_POLE, INTERLEAVED_BOOST.");
        }
    }

    bool PowerFactorCorrection::run_checks(bool assert) {
        bool valid = true;

        // Topology-variant gate FIRST so unsupported variants surface here
        // rather than from a downstream calc_inductance call.
        try {
            validate_topology_variant();
        } catch (const std::exception& e) {
            if (assert) throw;
            valid = false;
        }

        auto inputVoltage = get_input_voltage();
        double outputVoltage = get_output_voltage();
        double efficiency   = get_efficiency_required();

        // Check that output voltage is greater than peak input voltage
        double vinPeakMax = resolve_dimensional_values(inputVoltage, DimensionalValues::MAXIMUM) * std::sqrt(2);
        if (outputVoltage <= vinPeakMax) {
            if (assert) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "PFC output voltage must be greater than peak input voltage. Vout: " +
                    std::to_string(outputVoltage) + " <= Vin_peak_max: " + std::to_string(vinPeakMax));
            }
            valid = false;
        }

        // Check efficiency is reasonable
        if (efficiency <= 0 || efficiency > 1) {
            if (assert) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Efficiency must be between 0 and 1. Got: " + std::to_string(efficiency));
            }
            valid = false;
        }

        return valid;
    }

    double PowerFactorCorrection::calculate_duty_cycle(double vinPeak, double /*vout*/) {
        // For boost PFC: D = 1 - Vin/Vout
        // Accounting for diode drop: D = 1 - Vin/(Vout + Vd)
        // Totem-Pole replaces the boost diode with a synchronous switch, so
        // the conduction-path drop is essentially zero (RDS(on)·I, not Vf).
        const double vd = (get_topology_variant_or_default() == PfcTopologyVariants::TOTEM_POLE)
                              ? 0.0
                              : get_diode_voltage_drop_required();
        return 1.0 - vinPeak / (get_output_voltage() + vd);
    }

    // Per-phase power (Pout/N for INTERLEAVED_BOOST, Pout otherwise). Used by
    // the inductance-sizing formulas: each per-phase inductor carries 1/N of
    // the total line current, so its peak current — and therefore its sizing
    // current — scales by 1/N. Total stored magnetic energy across all N
    // inductors equals the single-phase boost: ½·L_phase·I_phase² · N
    // = ½·(N·L₁)·(I₁/N)²·N = ½·L₁·I₁².
    static double per_phase_power(const PowerFactorCorrection& pfc) {
        if (pfc.get_topology_variant_or_default() == PfcTopologyVariants::INTERLEAVED_BOOST) {
            return pfc.get_output_power() /
                   static_cast<double>(pfc.get_number_of_phases_or_default());
        }
        return pfc.get_output_power();
    }

    double PowerFactorCorrection::calculate_inductance_ccm() {
        validate_topology_variant();
        // For CCM PFC, worst case (maximum inductance requirement) is at minimum input voltage
        // at the peak of the AC line (where current is maximum)

        auto inputVoltage = get_input_voltage();
        double vinRmsMin  = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);
        double outputVoltage = get_output_voltage();

        // Calculate duty cycle at minimum input voltage peak
        double D = calculate_duty_cycle(vinPeakMin, outputVoltage);

        // Calculate input power (accounting for efficiency); per-phase for INTERLEAVED_BOOST
        double pinAvg = per_phase_power(*this) / get_efficiency_required();

        // Calculate average input current at minimum voltage
        if (vinRmsMin <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "PFC: minimum RMS input voltage must be positive");
        }
        double iinAvg = pinAvg / vinRmsMin;

        // Peak inductor current at line peak
        double iLpeak = iinAvg * std::sqrt(2);

        // Calculate ripple current based on ratio
        double deltaI = iLpeak * get_current_ripple_ratio_required();

        // CCM inductance formula
        double L = vinPeakMin * D / (deltaI * get_switching_frequency());

        return L;
    }

    double PowerFactorCorrection::calculate_inductance_dcm() {
        validate_topology_variant();
        auto inputVoltage = get_input_voltage();
        double vinRmsMin  = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);

        double D      = calculate_duty_cycle(vinPeakMin, get_output_voltage());
        double pinAvg = per_phase_power(*this) / get_efficiency_required();

        double L = pow(vinPeakMin, 2) * pow(D, 2) / (2 * pinAvg * get_switching_frequency());
        return L;
    }

    double PowerFactorCorrection::calculate_inductance_crcm() {
        validate_topology_variant();
        auto inputVoltage = get_input_voltage();
        double vinRmsMin  = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);

        double D      = calculate_duty_cycle(vinPeakMin, get_output_voltage());
        double pinAvg = per_phase_power(*this) / get_efficiency_required();
        double iinAvg = pinAvg / vinRmsMin;
        double iLpeak = iinAvg * std::sqrt(2) * 2;  // Peak current is 2x average in CrCM

        double L = vinPeakMin * D / (iLpeak * get_switching_frequency());
        return L;
    }

    double PowerFactorCorrection::calculate_peak_current(double vinPeak, double inductance) {
        double D = calculate_duty_cycle(vinPeak, get_output_voltage());

        double pinAvg = get_output_power() / get_efficiency_required();
        double iinRms = pinAvg / (vinPeak / std::sqrt(2));
        double iAvg = iinRms * std::sqrt(2);

        double deltaI = vinPeak * D / (inductance * get_switching_frequency());
        return iAvg + deltaI / 2;
    }

    std::string PowerFactorCorrection::determine_actual_mode(double inductance) {
        double L_critical = calculate_inductance_crcm();
        const double tolerance = 0.05;

        if (inductance > L_critical * (1.0 + tolerance)) {
            return "Continuous Conduction Mode";
        } else if (inductance < L_critical * (1.0 - tolerance)) {
            return "Discontinuous Conduction Mode";
        } else {
            return "Critical Conduction Mode";
        }
    }

    DesignRequirements PowerFactorCorrection::process_design_requirements() {
        validate_topology_variant();
        DesignRequirements designRequirements;

        // PFC inductor is a single winding - no turns ratio
        designRequirements.get_mutable_turns_ratios().clear();

        // Single winding
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(IsolationSide::PRIMARY);
        designRequirements.set_isolation_sides(isolationSides);

        // Set application - PFC uses POWER application for low-loss materials
        designRequirements.set_application(Application::POWER);
        // Tag the topology so downstream advisers (CoreAdviser, MagneticFilter,
        // CoilAdviser) can distinguish PFC from generic boost/buck inductors and
        // apply the correct routing — same as CMC/DMC.
        designRequirements.set_topology(Topologies::POWER_FACTOR_CORRECTION);

        // Calculate required inductance based on mode
        std::string mode = get_mode_string();
        double inductance;
        if (mode == "Continuous Conduction Mode") {
            inductance = calculate_inductance_ccm();
        } else if (mode == "Discontinuous Conduction Mode") {
            inductance = calculate_inductance_dcm();
        } else if (mode == "Critical Conduction Mode" || mode == "Transition Mode") {
            inductance = calculate_inductance_crcm();
        } else {
            inductance = calculate_inductance_ccm();
        }

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(inductance);
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        return designRequirements;
    }

    std::vector<OperatingPoint> PowerFactorCorrection::process_operating_points(const std::vector<double>& /*turnsRatios*/, double magnetizingInductance) {
        validate_topology_variant();
        std::vector<OperatingPoint> operatingPoints;

        std::string mode = get_mode_string();
        const auto variant = get_topology_variant_or_default();
        // Totem-pole has no bridge — inductor sees AC voltage / bipolar
        // current. Boost / Interleaved-Boost both have an input bridge → unipolar.
        const bool bipolar = (variant == PfcTopologyVariants::TOTEM_POLE);

        // Calculate required inductance if not provided
        double L = magnetizingInductance;
        if (L <= 0) {
            if (mode == "Continuous Conduction Mode") {
                L = calculate_inductance_ccm();
            } else if (mode == "Discontinuous Conduction Mode") {
                L = calculate_inductance_dcm();
            } else {
                L = calculate_inductance_crcm();
            }
        }

        auto inputVoltage = get_input_voltage();
        double vinRmsMin  = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);

        double outputVoltage     = get_output_voltage();
        double outputPower       = per_phase_power(*this);
        double switchingFrequency= get_switching_frequency();
        double lineFrequency     = get_line_frequency_required();
        // Totem-pole: synchronous low-side switch replaces the boost diode;
        // conduction-path Vf is essentially zero.
        double diodeVoltageDrop  = bipolar ? 0.0 : get_diode_voltage_drop_required();
        double efficiency        = get_efficiency_required();
        double ambientTemperature= get_ambient_temperature();

        double pinAvg     = outputPower / efficiency;
        double iinRmsAvg  = pinAvg / vinRmsMin;
        double iLinePeak  = iinRmsAvg * std::sqrt(2);

        double mainsPeriod     = 1.0 / lineFrequency;
        double switchingPeriod = 1.0 / switchingFrequency;

        int actualPeriods = (_numberOfPeriods > 0) ? _numberOfPeriods : 2;
        double totalTime  = mainsPeriod * actualPeriods;
        double timeStep   = switchingPeriod / 4.0;
        size_t numPoints  = static_cast<size_t>(totalTime / timeStep) + 1;

        std::vector<double> currentData;
        std::vector<double> voltageData;
        std::vector<double> timeData;
        currentData.reserve(numPoints);
        voltageData.reserve(numPoints);
        timeData.reserve(numPoints);

        for (size_t i = 0; i < numPoints; ++i) {
            double t = i * timeStep;
            timeData.push_back(t);

            double theta = 2.0 * std::numbers::pi * t / mainsPeriod;

            // Bridged (boost / interleaved-boost) → rectified |sin|;
            // bridgeless totem-pole → true sine (bipolar inductor current/voltage).
            const double vinShape = bipolar ? std::sin(theta) : std::abs(std::sin(theta));
            const double vinAbs   = std::abs(vinShape);
            double vinInst        = vinPeakMin * vinShape;
            double vinAbsInst     = vinPeakMin * vinAbs;
            // Floor on |Vin| to keep duty bounded near line zero-crossings.
            if (vinAbsInst < vinPeakMin * 0.05) {
                vinAbsInst = vinPeakMin * 0.05;
                if (bipolar) vinInst = std::copysign(vinAbsInst, vinShape);
                else         vinInst = vinAbsInst;
            }

            double D = 1.0 - vinAbsInst / (outputVoltage + diodeVoltageDrop);
            if (D < 0.01) D = 0.01;
            if (D > 0.95) D = 0.95;

            double iAvgInst = iLinePeak * vinShape;     // signed for totem-pole
            double deltaI   = vinAbsInst * D / (L * switchingFrequency);

            double tInSwitchCycle = std::fmod(t, switchingPeriod);
            double switchPhase    = tInSwitchCycle / switchingPeriod;

            double ripple;
            if (switchPhase < D) {
                ripple = -deltaI/2 + deltaI * (switchPhase / D);
            } else {
                ripple = deltaI/2 - deltaI * ((switchPhase - D) / (1 - D));
            }

            currentData.push_back(iAvgInst + ripple);

            if (switchPhase < D) {
                voltageData.push_back(vinInst);
            } else {
                // Off-time: inductor sees Vin − Vout (signed); for totem-pole
                // negative half-cycle the polarity of Vout is also reversed.
                const double voutSigned = bipolar ? std::copysign(outputVoltage, vinShape)
                                                  : outputVoltage;
                voltageData.push_back(vinInst - voutSigned - diodeVoltageDrop);
            }
        }

        Waveform currentWaveform;
        currentWaveform.set_data(currentData);
        currentWaveform.set_time(timeData);

        SignalDescriptor current;
        current.set_waveform(currentWaveform);
        auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, lineFrequency);
        current.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform, lineFrequency));
        auto currentProcessed = Inputs::calculate_processed_data(currentWaveform, lineFrequency, true);
        currentProcessed.set_label(WaveformLabel::CUSTOM);
        current.set_processed(currentProcessed);

        Waveform voltageWaveform;
        voltageWaveform.set_data(voltageData);
        voltageWaveform.set_time(timeData);

        SignalDescriptor voltage;
        voltage.set_waveform(voltageWaveform);
        auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, lineFrequency);
        voltage.set_harmonics(Inputs::calculate_harmonics_data(sampledVoltageWaveform, lineFrequency));
        auto voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, lineFrequency, true);
        voltageProcessed.set_label(WaveformLabel::CUSTOM);
        voltage.set_processed(voltageProcessed);

        OperatingPointExcitation excitation;
        excitation.set_current(current);
        excitation.set_frequency(lineFrequency);
        excitation.set_voltage(voltage);

        OperatingPoint operatingPoint;
        operatingPoint.set_excitations_per_winding({excitation});
        operatingPoint.get_mutable_conditions().set_ambient_temperature(ambientTemperature);
        operatingPoint.set_name("PFC_" +
            std::string(bipolar ? "TotemPole_" :
                        (variant == PfcTopologyVariants::INTERLEAVED_BOOST
                            ? "InterleavedBoost_" : "Boost_")) +
            "Line_Period_" + std::to_string(actualPeriods) + "_Periods");

        operatingPoints.push_back(operatingPoint);

        return operatingPoints;
    }

    std::string PowerFactorCorrection::generate_ngspice_circuit(double inductance,
                                                                 double /*dcResistance*/,
                                                                 double simulationTime,
                                                                 double timeStep) {
        validate_topology_variant();
        std::ostringstream netlist;

        const auto variant = get_topology_variant_or_default();
        const char* variantName = "Boost";
        switch (variant) {
            case PfcTopologyVariants::TOTEM_POLE:        variantName = "Totem-Pole";        break;
            case PfcTopologyVariants::INTERLEAVED_BOOST: variantName = "Interleaved Boost"; break;
            default: break;
        }
        netlist << "* PFC " << variantName << " Converter - Ideal Behavioral Model\n";
        netlist << "* Generated by OpenMagnetics\n";
        netlist << "* Models ideal PFC with sinusoidal current envelope + switching ripple\n\n";

        auto inputVoltage = get_input_voltage();
        double vinRms;
        if (inputVoltage.get_nominal().has_value()) {
            vinRms = inputVoltage.get_nominal().value();
        } else {
            vinRms = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        }
        double vinPeak = vinRms * std::sqrt(2);
        double vout    = get_output_voltage();
        double pout    = get_output_power();
        double switchingFrequency = get_switching_frequency();
        double lineFrequency      = get_line_frequency_required();

        double iPeak = std::sqrt(2) * pout / vinRms;
        double tSw   = 1.0 / switchingFrequency;

        netlist << ".param vin_peak=" << vinPeak << "\n";
        netlist << ".param vout=" << vout << "\n";
        netlist << ".param fline=" << lineFrequency << "\n";
        netlist << ".param fsw=" << switchingFrequency << "\n";
        netlist << ".param L=" << inductance << "\n";
        netlist << ".param i_peak=" << iPeak << "\n\n";

        netlist << "* Rectified AC Input\n";
        netlist << "B_vin vin_rect 0 V=vin_peak*abs(sin(2*3.14159265*fline*time))\n\n";

        netlist << "* Ideal current envelope (sinusoidal, in phase with voltage)\n";
        netlist << "B_ienv i_env 0 V=i_peak*abs(sin(2*3.14159265*fline*time))\n\n";

        netlist << "* Instantaneous duty cycle\n";
        netlist << "B_duty duty 0 V=1-V(vin_rect)/vout\n\n";

        netlist << "* Current ripple amplitude\n";
        netlist << "B_rip ripple 0 V=V(vin_rect)*V(duty)/(L*fsw)/2\n\n";

        netlist << "* Sawtooth for triangular switching ripple\n";
        netlist << "V_saw saw 0 PULSE(-1 1 0 " << tSw/2 << " " << tSw/2 << " 1n " << tSw << ")\n\n";

        netlist << "* Total inductor current (envelope + ripple)\n";
        netlist << "B_iL i_L 0 V=V(i_env)+V(ripple)*V(saw)\n\n";

        netlist << "* Analysis\n";
        netlist << ".tran " << timeStep << " " << simulationTime << " 0 " << timeStep << "\n";
        netlist << ".save v(vin_rect) v(i_env) v(i_L) v(ripple)\n";
        netlist << ".end\n";

        return netlist.str();
    }

    PfcSimulationWaveforms PowerFactorCorrection::simulate_and_extract_waveforms(double inductance,
                                                                                   double /*dcResistance*/,
                                                                                   int numberOfCycles) {
        validate_topology_variant();
        PfcSimulationWaveforms waveforms;
        const auto variant = get_topology_variant_or_default();
        const bool bipolar = (variant == PfcTopologyVariants::TOTEM_POLE);
        double switchingFrequency = get_switching_frequency();
        double lineFrequency      = get_line_frequency_required();
        waveforms.switchingFrequency = switchingFrequency;
        waveforms.lineFrequency      = lineFrequency;

        auto inputVoltage = get_input_voltage();
        double vinRms;
        if (inputVoltage.get_nominal().has_value()) {
            vinRms = inputVoltage.get_nominal().value();
        } else {
            vinRms = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        }
        double vinPeak = vinRms * std::sqrt(2);
        double vout    = get_output_voltage();
        // Per-phase power for INTERLEAVED_BOOST.
        double pout    = per_phase_power(*this);

        double iPeak = std::sqrt(2) * pout / vinRms;

        double linePeriod      = 1.0 / lineFrequency;
        double switchingPeriod = 1.0 / switchingFrequency;
        double simulationTime  = numberOfCycles * linePeriod;

        double timeStep   = switchingPeriod / 4.0;
        size_t numPoints  = static_cast<size_t>(simulationTime / timeStep) + 1;

        waveforms.time.resize(numPoints);
        waveforms.inputVoltage.resize(numPoints);
        waveforms.inductorCurrent.resize(numPoints);
        waveforms.inductorVoltage.resize(numPoints);
        waveforms.currentEnvelope.resize(numPoints);
        waveforms.currentRipple.resize(numPoints);

        double omega_line = 2.0 * M_PI * lineFrequency;

        for (size_t i = 0; i < numPoints; i++) {
            double t = i * timeStep;
            waveforms.time[i] = t;

            // Bridged → unipolar |sin|; bridgeless totem-pole → bipolar sin.
            const double shape = bipolar ? std::sin(omega_line * t)
                                          : std::abs(std::sin(omega_line * t));
            const double absShape = std::abs(shape);
            double vin = vinPeak * shape;
            waveforms.inputVoltage[i] = vin;

            double iEnv = iPeak * shape;
            waveforms.currentEnvelope[i] = iEnv;

            double duty = 1.0 - (vinPeak * absShape) / vout;
            if (duty < 0) duty = 0;
            if (duty > 1) duty = 1;

            double rippleAmplitude = (vinPeak * absShape) * duty / (2.0 * inductance * switchingFrequency);
            waveforms.currentRipple[i] = rippleAmplitude;

            double tInSwitchingCycle = std::fmod(t, switchingPeriod);
            double sawtoothPhase = tInSwitchingCycle / switchingPeriod;
            double triangular;
            if (sawtoothPhase < 0.5) {
                triangular = -1.0 + 4.0 * sawtoothPhase;
            } else {
                triangular = 3.0 - 4.0 * sawtoothPhase;
            }

            waveforms.inductorCurrent[i] = iEnv + rippleAmplitude * triangular;

            // Off-time inductor voltage = Vin − Vout (signed for totem-pole).
            const double voutSigned = bipolar ? std::copysign(vout, shape) : vout;
            if (sawtoothPhase < 0.5) {
                waveforms.inductorVoltage[i] = vin;
            } else {
                waveforms.inductorVoltage[i] = vin - voutSigned;
            }
        }

        waveforms.inputCurrent = waveforms.inductorCurrent;

        double realPower = 0, vRms = 0, iRms = 0;
        for (size_t i = 0; i < numPoints; i++) {
            realPower += waveforms.inputVoltage[i] * waveforms.inputCurrent[i];
            vRms      += waveforms.inputVoltage[i] * waveforms.inputVoltage[i];
            iRms      += waveforms.inputCurrent[i] * waveforms.inputCurrent[i];
        }
        realPower /= numPoints;
        vRms = std::sqrt(vRms / numPoints);
        iRms = std::sqrt(iRms / numPoints);

        double apparentPower = vRms * iRms;
        waveforms.powerFactor = (apparentPower > 0) ? realPower / apparentPower : 1.0;
        waveforms.efficiency  = 1.0;
        waveforms.operatingPointName = "PFC_analytical";

        return waveforms;
    }

    std::vector<OperatingPoint> PowerFactorCorrection::simulate_and_extract_operating_points(double inductance,
                                                                                               double dcResistance) {
        std::vector<OperatingPoint> operatingPoints;

        int numPeriods = (_numberOfPeriods > 0) ? _numberOfPeriods : 2;
        PfcSimulationWaveforms waveforms = simulate_and_extract_waveforms(inductance, dcResistance, numPeriods);

        if (waveforms.inductorCurrent.empty() || waveforms.time.empty()) {
            return process_operating_points({}, inductance);
        }

        double switchingFrequency = get_switching_frequency();
        double lineFrequency      = get_line_frequency_required();
        double ambientTemperature = get_ambient_temperature();

        Waveform currentWaveform;
        currentWaveform.set_time(waveforms.time);
        currentWaveform.set_data(waveforms.inductorCurrent);

        SignalDescriptor current;
        current.set_waveform(currentWaveform);

        try {
            auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
            current.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency));
            auto currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true);
            currentProcessed.set_label(WaveformLabel::CUSTOM);
            current.set_processed(currentProcessed);
        } catch (const std::exception&) {
            Processed processed;
            if (!waveforms.inductorCurrent.empty()) {
                double iMax = *std::max_element(waveforms.inductorCurrent.begin(), waveforms.inductorCurrent.end());
                double iMin = *std::min_element(waveforms.inductorCurrent.begin(), waveforms.inductorCurrent.end());
                double iAvg = std::accumulate(waveforms.inductorCurrent.begin(), waveforms.inductorCurrent.end(), 0.0) /
                              waveforms.inductorCurrent.size();
                processed.set_peak(iMax);
                processed.set_peak_to_peak(iMax - iMin);
                processed.set_average(iAvg);

                double sumSq = 0;
                for (double val : waveforms.inductorCurrent) {
                    sumSq += val * val;
                }
                processed.set_rms(std::sqrt(sumSq / waveforms.inductorCurrent.size()));
            }
            current.set_processed(processed);
        }

        SignalDescriptor voltage;
        if (!waveforms.inductorVoltage.empty() && waveforms.inductorVoltage.size() == waveforms.time.size()) {
            Waveform voltageWaveform;
            voltageWaveform.set_time(waveforms.time);
            voltageWaveform.set_data(waveforms.inductorVoltage);
            voltage.set_waveform(voltageWaveform);

            try {
                auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
                voltage.set_harmonics(Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency));
                auto voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true);
                voltageProcessed.set_label(WaveformLabel::CUSTOM);
                voltage.set_processed(voltageProcessed);
            } catch (const std::exception&) {
                // Just use raw waveform without harmonics
            }
        }

        OperatingPointExcitation excitation;
        excitation.set_current(current);
        excitation.set_frequency(lineFrequency);
        excitation.set_voltage(voltage);

        OperatingPoint operatingPoint;
        operatingPoint.set_excitations_per_winding({excitation});
        operatingPoint.get_mutable_conditions().set_ambient_temperature(ambientTemperature);
        operatingPoint.set_name("PFC_simulated_" + std::to_string(static_cast<int>(lineFrequency)) + "Hz");

        operatingPoints.push_back(operatingPoint);

        return operatingPoints;
    }

    std::vector<ConverterWaveforms>
    PowerFactorCorrection::simulate_and_extract_topology_waveforms(double inductance,
                                                                    size_t numberOfCycles) {
        validate_topology_variant();
        std::vector<ConverterWaveforms> results;

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);
        if (inputVoltages.empty()) {
            throw std::runtime_error(
                "PowerFactorCorrection::simulate_and_extract_topology_waveforms: "
                "no input voltage (nominal/min/max) populated");
        }

        const double switchingFrequency = get_switching_frequency();
        const double lineFrequency      = get_line_frequency_required();
        const double vout               = get_output_voltage();
        const double pout               = get_output_power();
        const double efficiency         = get_efficiency_required();

        // Bulk-capacitance is required for the §5.1 converter-port output
        // model (without it the DC bus would be modelled as flat, which
        // hides the twice-line-frequency ripple that drives bulk-cap
        // sizing). Per the repo "no fallbacks — throw" rule, throw if unset.
        auto bulkCapOpt = MAS::PowerFactorCorrection::get_bulk_capacitance();
        if (!bulkCapOpt.has_value()) {
            throw std::runtime_error(
                "PowerFactorCorrection::simulate_and_extract_topology_waveforms: "
                "bulkCapacitance is required (set via set_bulk_capacitance() or "
                "the JSON field 'bulkCapacitance'). Typical sizing: 1–4 µF/W for "
                "400 V bus / 20 ms hold-up.");
        }
        const double bulkCapacitance = bulkCapOpt.value();

        const double linePeriod      = 1.0 / lineFrequency;
        const double switchingPeriod = 1.0 / switchingFrequency;
        const int    actualCycles    = std::max<size_t>(numberOfCycles, 1);
        const double simulationTime  = actualCycles * linePeriod;
        const double timeStep        = switchingPeriod / 4.0;
        const size_t numPoints       = static_cast<size_t>(simulationTime / timeStep) + 1;
        const double omegaLine       = 2.0 * M_PI * lineFrequency;

        for (size_t vinIdx = 0; vinIdx < inputVoltages.size(); ++vinIdx) {
            const double vinRms  = inputVoltages[vinIdx];
            const double vinPeak = vinRms * std::sqrt(2.0);
            // Average input current (RMS, ideal PF=1): I_in_rms = P_in / V_in_rms,
            // P_in = P_out / η. Peak of the |sin| envelope is √2·I_in_rms.
            const double pinAvg    = pout / efficiency;
            const double iLinePeak = std::sqrt(2.0) * pinAvg / vinRms;
            // Load current (mean DC; the instantaneous value tracks Vbus(t)/Rload below).
            const double iLoad = pout / vout;
            const double rLoad = vout / iLoad;

            // ──────────────────────────────────────────────────────────────
            // DC bus ripple — linearised state-variable model.
            //
            // Bulk-cap ODE:  Cbus · dVbus/dt = i_diode(t) - Vbus(t)/Rload
            // where the cycle-averaged diode current
            //   i_diode_avg(t) = i_line(t)·(1−D(t)) = (Pin_avg / Vbus)·(1 − cos 2ωt)
            // Linearising around Vbus ≈ Vbus_nom (small ripple, ω >> 1/(Rload·Cbus)):
            //   ΔV(t) ≈ −[Pin_avg / (2·ω_line·Cbus·Vbus_nom)] · sin(2·ω_line·t)
            // so the peak-to-peak bus ripple is
            //   ΔVbus_pp ≈ Pout / (π · f_line · Cbus · Vbus)
            // This is the textbook "low-frequency bus ripple" formula used to
            // size Cbus for both ripple-current rating and the 2·f_line hum-loop
            // shaping in the voltage loop (Erickson §17, TI SLUA754).
            //
            // Mean of Vbus(t) remains Vbus_nom (sinusoid is zero-mean), so
            // ConverterPortChecks::check_pfc_ports still passes.
            // ──────────────────────────────────────────────────────────────
            const double rippleAmpPk = pinAvg / (2.0 * omegaLine * bulkCapacitance * vout);

            std::vector<double> time(numPoints);
            std::vector<double> vinData(numPoints);
            std::vector<double> iinData(numPoints);
            std::vector<double> voutData(numPoints);
            std::vector<double> ioutData(numPoints);

            for (size_t i = 0; i < numPoints; ++i) {
                const double t       = i * timeStep;
                time[i]              = t;
                const double absSin  = std::abs(std::sin(omegaLine * t));
                const double vinInst = vinPeak * absSin;
                vinData[i]           = vinInst;

                // Triangular switching ripple superimposed on the rectified-sine envelope.
                const double iEnv  = iLinePeak * absSin;
                double duty = 1.0 - vinInst / vout;
                if (duty < 0.0) duty = 0.0;
                if (duty > 1.0) duty = 1.0;
                const double rippleAmp = vinInst * duty / (2.0 * inductance * switchingFrequency);
                const double phase     = std::fmod(t, switchingPeriod) / switchingPeriod;
                const double tri       = (phase < 0.5) ? (-1.0 + 4.0 * phase)
                                                       : ( 3.0 - 4.0 * phase);
                iinData[i] = iEnv + rippleAmp * tri;

                const double vbus = vout - rippleAmpPk * std::sin(2.0 * omegaLine * t);
                voutData[i] = vbus;
                ioutData[i] = vbus / rLoad;
            }

            ConverterWaveforms wf;
            wf.set_switching_frequency(switchingFrequency);
            wf.set_operating_point_name(inputVoltagesNames[vinIdx] + " input (PFC)");

            Waveform vinWf;  vinWf.set_data(vinData);   vinWf.set_time(time);
            Waveform iinWf;  iinWf.set_data(iinData);   iinWf.set_time(time);
            Waveform voutWf; voutWf.set_data(voutData); voutWf.set_time(time);
            Waveform ioutWf; ioutWf.set_data(ioutData); ioutWf.set_time(time);

            wf.set_input_voltage(vinWf);
            wf.set_input_current(iinWf);
            wf.get_mutable_output_voltages().push_back(voutWf);
            wf.get_mutable_output_currents().push_back(ioutWf);

            results.push_back(wf);
        }

        return results;
    }

    // ─────────────────────────────────────────────────────────────────────
    // Switching boost-PFC netlist + ngspice runner.
    //
    // Unlike `generate_ngspice_circuit`, which emits a behavioural-source
    // synthesis of the analytical model, the methods below build a true
    // switching boost converter (B_vin rectified-sine source → L1 → ideal
    // SW S1 → ideal D1 → Cbus → Rload) driven by the native
    // average-current-mode controller documented in
    // `PfcControllerDesign.h` (six analog blocks: voltage error amp,
    // RMS² feed-forward, multiplier, current error amp, sawtooth, PWM
    // comparator). All controller component values come from
    // `derive_pfc_controller_tuning()` in closed form from the converter
    // spec; subcircuits (OPAMP_IDEAL, COMPARATOR_IDEAL) come from
    // `PfcControllerSubcircuits.h`. Initial conditions warm-start every
    // state variable to its analytical steady-state value to skip the
    // multi-line-cycle Cbus charging transient.
    //
    // The intent is to drive an INDEPENDENT (switching-circuit) check of
    // the PFC analytical model, used by the §3.D Phase 6 industry
    // reference-design PtP gates.
    // ─────────────────────────────────────────────────────────────────────

    namespace {
        // sci(n): force scientific notation regardless of stream state. Used
        // to keep .param strings stable across compilers / locales (some
        // glibc locales render e.g. 0.0001 as "0,0001" by default).
        std::string sci(double v) {
            std::ostringstream o;
            o.imbue(std::locale::classic());
            o << std::scientific << std::setprecision(8) << v;
            return o.str();
        }
    } // namespace

    std::string PowerFactorCorrection::generate_ngspice_switching_circuit(
            double inductance,
            int numberOfLineCycles) {
        validate_topology_variant();
        const auto variant = get_topology_variant_or_default();
        if (variant != PfcTopologyVariants::BOOST) {
            throw std::runtime_error(
                "PowerFactorCorrection::generate_ngspice_switching_circuit: "
                "only the BOOST variant is supported by the native switching "
                "controller; TOTEM_POLE / INTERLEAVED_BOOST require their own "
                "switching netlists (out of scope here).");
        }
        if (!std::isfinite(inductance) || inductance <= 0.0) {
            throw std::invalid_argument(
                "PowerFactorCorrection::generate_ngspice_switching_circuit: "
                "inductance must be positive and finite (got " +
                std::to_string(inductance) + ")");
        }
        if (numberOfLineCycles < 1) numberOfLineCycles = 1;

        // ── Spec resolution ─────────────────────────────────────────────
        auto inputVoltage = get_input_voltage();
        double vrmsNom;
        if (inputVoltage.get_nominal().has_value()) {
            vrmsNom = inputVoltage.get_nominal().value();
        } else {
            vrmsNom = resolve_dimensional_values(inputVoltage, DimensionalValues::NOMINAL);
        }
        const double vbusNom = get_output_voltage();
        const double pout    = get_output_power();
        const double fsw     = get_switching_frequency();
        const double fline   = get_line_frequency_required();
        if (pout <= 0.0) {
            throw std::invalid_argument(
                "PowerFactorCorrection::generate_ngspice_switching_circuit: "
                "output power must be positive (got " + std::to_string(pout) + ")");
        }
        auto bulkCapOpt = MAS::PowerFactorCorrection::get_bulk_capacitance();
        if (!bulkCapOpt.has_value()) {
            throw std::runtime_error(
                "PowerFactorCorrection::generate_ngspice_switching_circuit: "
                "bulkCapacitance is required for the switching netlist "
                "(set via set_bulk_capacitance() or JSON 'bulkCapacitance'). "
                "Typical sizing for a boost-PFC: 1–2 µF/W.");
        }
        const double cbus  = bulkCapOpt.value();
        const double vpk   = std::sqrt(2.0) * vrmsNom;
        const double rload = vbusNom * vbusNom / pout;
        const double tsw   = 1.0 / fsw;
        const double tline = 1.0 / fline;

        // ── Controller tuning (closed-form from spec) ──────────────────
        const PfcControllerTuning t = derive_pfc_controller_tuning(
            vrmsNom, vbusNom, pout, fsw, fline, inductance);

        std::ostringstream c;
        c.imbue(std::locale::classic());

        // ── Header ──────────────────────────────────────────────────────
        c << "* PFC Boost Switching Converter (native average-current-mode controller)\n";
        c << "* Generated by OpenMagnetics — see PfcControllerDesign.h\n";
        c << "* Vrms=" << vrmsNom << " V, Vbus=" << vbusNom << " V, Pout=" << pout
          << " W, Fsw=" << (fsw/1e3) << " kHz, fline=" << fline
          << " Hz, L=" << (inductance*1e6) << " uH, Cbus=" << (cbus*1e6)
          << " uF, Rload=" << rload << " Ohm\n";
        c << "* Cycles to simulate: " << numberOfLineCycles
          << " (line period = " << (tline*1e3) << " ms)\n\n";

        // ── Parameters ──────────────────────────────────────────────────
        c << ".param vpk="              << sci(vpk)              << "\n";
        c << ".param vbus_nom="         << sci(vbusNom)          << "\n";
        c << ".param fline="            << sci(fline)            << "\n";
        c << ".param fsw="              << sci(fsw)              << "\n";
        c << ".param tsw="              << sci(tsw)              << "\n";
        c << ".param l_ind="            << sci(inductance)       << "\n";
        c << ".param cbus="             << sci(cbus)             << "\n";
        c << ".param rload="            << sci(rload)            << "\n\n";

        c << ".param vref="             << sci(t.vref)           << "\n";
        c << ".param k_div="            << sci(t.k_div)          << "\n";
        c << ".param vea_min="          << sci(t.vea_min)        << "\n";
        c << ".param vea_max="          << sci(t.vea_max)        << "\n";
        c << ".param gv_rz="            << sci(t.gv_rz)          << "\n";
        c << ".param gv_rin="           << sci(t.gv_rz)          << "  ; mid-band gain = 1\n";
        c << ".param gv_cz="            << sci(t.gv_cz)          << "\n";
        c << ".param gv_cp="            << sci(t.gv_cp)          << "\n\n";

        c << ".param ff_r="             << sci(t.ff_r)           << "\n";
        c << ".param ff_c="             << sci(t.ff_c)           << "\n";
        c << ".param vrms_ff_floor="    << sci(t.vrms_ff_floor)  << "\n";
        c << ".param ic_vff="           << sci(t.ic_vff)         << "\n\n";

        c << ".param g_mul="            << sci(t.g_mul)          << "\n";
        c << ".param rs_sense="         << sci(t.rs_sense)       << "\n\n";

        c << ".param gi_rz="            << sci(t.gi_rz)          << "\n";
        c << ".param gi_rin="           << sci(t.gi_rin)         << "\n";
        c << ".param gi_cz="            << sci(t.gi_cz)          << "\n";
        c << ".param gi_cp="            << sci(t.gi_cp)          << "\n\n";

        c << ".param v_pk_saw="         << sci(t.pwm_v_pk_saw)   << "\n";
        c << ".param v_high="           << sci(t.pwm_v_high)     << "\n";
        c << ".param pwm_t_rise="       << sci(t.pwm_t_rise)     << "\n";
        c << ".param ic_vbus="          << sci(t.ic_vbus)        << "\n";
        c << ".param ic_vea="           << sci(t.ic_vea)         << "\n\n";

        // ── Subckt prelude (OPAMP_IDEAL + COMPARATOR_IDEAL) ─────────────
        c << "* ─── Reusable subcircuits ──────────────────────────────────\n";
        c << spice_subckt_prelude_pfc_controller();
        c << "\n";

        // ── Power stage ─────────────────────────────────────────────────
        c << "* ─── Power stage (Boost) ───────────────────────────────────\n";
        c << "B_vin     vin_rect 0     V=vpk*abs(sin(2*3.141592653589793*fline*time))\n";
        c << "Vl_sense  vin_rect l_in  0          ; in-line current sense\n";
        c << "L1        l_in     sw    {l_ind}\n";
        c << ".model SW1 SW (VT=0.5 VH=0.1 RON=10m ROFF=1MEG)\n";
        c << "S1        sw  0    gate 0 SW1\n";
        c << "Rsnub_s1  sw  0    100\n";
        c << "Csnub_s1  sw  0    100p\n";
        c << ".model DIDEAL D (IS=1e-12 RS=1m N=1)\n";
        c << "D1        sw  vbus DIDEAL\n";
        c << "Cout      vbus 0   {cbus}  IC={ic_vbus}\n";
        c << "Rload     vbus 0   {rload}\n\n";

        // ── Voltage error amp (Block 1) ─────────────────────────────────
        c << "* ─── Voltage error amp (Gv) — type-II ─────────────────────\n";
        c << "B_div     vbus_div 0  V=V(vbus)*k_div\n";
        c << "V_ref_v   vref_v   0  {vref}\n";
        c << "R_in_v    vbus_div vea_n  {gv_rin}\n";
        c << "R_fb_z    vea_n    vea_z  {gv_rz}\n";
        c << "C_fb_z    vea_z    vea    {gv_cz} IC=0\n";
        c << "C_fb_p    vea_n    vea    {gv_cp} IC=0\n";
        c << "X_op_v    vref_v   vea_n  vea OPAMP_IDEAL "
             "params: A0=1e5 GBW=1e6 VSSPOS={vea_max} VSSNEG={vea_min}\n\n";

        // ── Feed-forward (Block 2) ──────────────────────────────────────
        c << "* ─── RMS^2 feed-forward — two LP stages + squarer ─────────\n";
        c << "R_ff1     vin_rect vff1   {ff_r}\n";
        c << "C_ff1     vff1     0      {ff_c} IC={ic_vff}\n";
        c << "R_ff2     vff1     vff2   {ff_r}\n";
        c << "C_ff2     vff2     0      {ff_c} IC={ic_vff}\n";
        c << "B_squarer vrms_ff  0      V=V(vff2)*V(vff2)\n\n";

        // ── Multiplier (Block 3) ────────────────────────────────────────
        c << "* ─── Multiplier — i_ref = G_mul · vea · vin_rect / vrms_ff ─\n";
        c << "B_iref    i_ref    0  "
             "V=g_mul*V(vea)*V(vin_rect)/max(V(vrms_ff),vrms_ff_floor)*rs_sense\n\n";

        // ── Current sense + error amp (Block 4) ─────────────────────────
        c << "* ─── Current sense + Gi (type-II) ─────────────────────────\n";
        c << "B_isense  i_sense  0  V=I(Vl_sense)*rs_sense\n";
        c << "R_in_i    i_sense  ic_n   {gi_rin}\n";
        c << "R_fb_zi   ic_n     ic_z   {gi_rz}\n";
        c << "C_fb_zi   ic_z     vc_i   {gi_cz} IC=0\n";
        c << "C_fb_pi   ic_n     vc_i   {gi_cp} IC=0\n";
        c << "X_op_i    i_ref    ic_n   vc_i OPAMP_IDEAL "
             "params: A0=1e5 GBW=10e6 VSSPOS={v_pk_saw} VSSNEG=0\n\n";

        // ── Oscillator + PWM comparator (Blocks 5–6) ────────────────────
        c << "* ─── Sawtooth + PWM comparator ────────────────────────────\n";
        c << "V_saw     saw 0   PULSE(0 {v_pk_saw} 0 {tsw-pwm_t_rise} "
             "{pwm_t_rise} {pwm_t_rise} {tsw})\n";
        c << "B_gate    gate 0  V=(V(vc_i) > V(saw)) ? v_high : 0\n\n";

        // ── Initial conditions ──────────────────────────────────────────
        c << "* ─── Initial conditions (warm start to analytical SS) ────\n";
        c << ".ic v(vbus)={ic_vbus} v(vea)={ic_vea} "
             "v(vff1)={ic_vff} v(vff2)={ic_vff}\n\n";

        // ── Transient analysis + solver options ─────────────────────────
        const double simTime = numberOfLineCycles * tline;
        const double maxStep = tsw / 40.0;
        c << "* ─── Analysis ──────────────────────────────────────────────\n";
        c << ".tran " << sci(maxStep) << " " << sci(simTime) << " 0 "
          << sci(maxStep) << " uic\n";
        c << ".save i(vl_sense) v(vin_rect) v(vbus) v(vea) v(vrms_ff) "
             "v(i_ref) v(i_sense) v(vc_i) v(saw) v(gate)\n";
        c << ".options METHOD=GEAR TRTOL=7 RELTOL=1e-3 ABSTOL=1e-9 VNTOL=1e-6\n";
        c << ".options ITL1=500 ITL4=200\n";
        c << ".end\n";

        return c.str();
    }

    OperatingPoint PowerFactorCorrection::simulate_with_ngspice_switching(
            double inductance,
            int numberOfLineCycles) {
        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error(
                "PowerFactorCorrection::simulate_with_ngspice_switching: "
                "ngspice is not available on this system");
        }

        const std::string netlist =
            generate_ngspice_switching_circuit(inductance, numberOfLineCycles);

        SimulationConfig config;
        // Switching netlist already specifies the full transient window
        // (numberOfLineCycles · Tline). We want the entire raw waveform
        // back; the caller windows it to the cycle of interest.
        config.frequency        = get_line_frequency_required();
        config.extractOnePeriod = false;
        config.numberOfPeriods  = static_cast<size_t>(std::max(1, numberOfLineCycles));
        config.keepTempFiles    = false;
        // 3 line cycles × thousands of switching events: raise timeout.
        config.timeout          = 300.0;

        auto sim = runner.run_simulation(netlist, config);
        if (!sim.success) {
            throw std::runtime_error(
                "PowerFactorCorrection::simulate_with_ngspice_switching: "
                "ngspice simulation failed: " + sim.errorMessage);
        }

        // Locate the inductor-current waveform — ngspice typically returns
        // i(vl_sense) (or vl_sense#branch in some builds). Case-insensitive
        // suffix match on "vl_sense" hits both forms.
        auto matchesILSense = [](const std::string& n) {
            std::string s = n;
            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            return s.find("vl_sense") != std::string::npos;
        };
        const Waveform* iLwf = nullptr;
        for (size_t k = 0; k < sim.waveformNames.size(); ++k) {
            if (matchesILSense(sim.waveformNames[k])) {
                iLwf = &sim.waveforms[k];
                break;
            }
        }
        if (!iLwf) {
            std::ostringstream e;
            e << "PowerFactorCorrection::simulate_with_ngspice_switching: "
                 "inductor-current waveform 'i(vl_sense)' not found in "
                 "ngspice output. Available names:";
            for (const auto& n : sim.waveformNames) e << " " << n;
            throw std::runtime_error(e.str());
        }

        // Wrap into an OperatingPoint with one excitation carrying i_L.
        OperatingPoint op;
        OperatingPointExcitation exc;
        SignalDescriptor curSig;
        curSig.set_waveform(*iLwf);
        exc.set_current(curSig);
        op.get_mutable_excitations_per_winding().push_back(exc);
        return op;
    }

} // namespace OpenMagnetics
