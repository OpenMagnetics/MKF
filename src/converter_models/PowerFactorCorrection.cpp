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

            case PfcTopologyVariants::BRIDGELESS:
            case PfcTopologyVariants::SEMI_BRIDGELESS:
                // Bridgeless / semi-bridgeless boost. These remove (or halve)
                // the input rectifier bridge to cut conduction loss, but the
                // BOOST INDUCTOR itself sees the identical rectified-sine
                // current envelope as a classic boost PFC — the only
                // difference is which return-path devices conduct. From the
                // magnetic-component standpoint they are boost: same
                // inductance sizing (calculate_inductance_*), same
                // |sin|-enveloped unipolar waveform synthesis, same switching
                // netlist. Unlike totem-pole (which is also bridgeless but
                // puts the inductor on the AC side → bipolar current), these
                // keep the inductor on the rectified side, so no
                // wide-bandgap-switch requirement and no companion fields.
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

            case PfcTopologyVariants::SEPIC:
            case PfcTopologyVariants::CUK:
                // Buck-boost-class with a series input inductor (continuous
                // input current). Analytical sizing/waveforms are supported
                // via the buck-boost duty relation + the shared L1 = Vin·D/
                // (ΔI·fsw) form (see is_buck_boost_class). Their 4th-order
                // power stage has no native switching controller yet, so
                // generate_ngspice_switching_circuit still rejects them
                // (no switching-PtP). DCM sizing is not yet validated for
                // these and is rejected in calculate_inductance_dcm.
                return;

            case PfcTopologyVariants::BUCK:
            case PfcTopologyVariants::BUCK_BOOST:
                // Input current is discontinuous (inductor not in series with
                // the line), so these are not boost-/SEPIC-class and have no
                // shared sizing path. Unsupported.
                throw std::runtime_error(
                    "PowerFactorCorrection: buck / buckBoost PFC is not "
                    "implemented (discontinuous input current — no shared "
                    "series-inductor sizing path). Supported: BOOST, "
                    "BRIDGELESS, SEMI_BRIDGELESS, TOTEM_POLE, "
                    "INTERLEAVED_BOOST, SEPIC, CUK.");
        }
    }

    double PowerFactorCorrection::suggested_bulk_capacitance(double vBus,
                                                              double vBusMin,
                                                              double pOut,
                                                              double holdupTime) {
        auto check = [](const char* name, double v) {
            if (!std::isfinite(v) || v <= 0.0) {
                throw std::invalid_argument(
                    std::string("PowerFactorCorrection::suggested_bulk_capacitance: ")
                    + name + " must be positive and finite (got "
                    + std::to_string(v) + ")");
            }
        };
        check("vBus",       vBus);
        check("vBusMin",    vBusMin);
        check("pOut",       pOut);
        check("holdupTime", holdupTime);
        if (vBusMin >= vBus) {
            throw std::invalid_argument(
                "PowerFactorCorrection::suggested_bulk_capacitance: "
                "vBusMin must be strictly less than vBus (got vBus=" +
                std::to_string(vBus) + ", vBusMin=" + std::to_string(vBusMin) + ")");
        }
        return 2.0 * pOut * holdupTime / (vBus * vBus - vBusMin * vBusMin);
    }

    double PowerFactorCorrection::get_vrms_for_simulation() const {
        auto inputVoltage = const_cast<PowerFactorCorrection*>(this)->get_input_voltage();
        if (inputVoltage.get_nominal().has_value()) {
            return inputVoltage.get_nominal().value();
        }
        const double v = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        if (!std::isfinite(v) || v <= 0.0) {
            throw std::runtime_error(
                "PowerFactorCorrection::get_vrms_for_simulation: inputVoltage "
                "has neither a nominal value nor a usable minimum");
        }
        return v;
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

        // Boost-family PFCs can only step UP, so Vout must exceed the peak
        // input. SEPIC/Ćuk are buck-boost-class and may regulate a bus below
        // the peak line, so this constraint does not apply to them.
        if (!is_buck_boost_class()) {
            double vinPeakMax = resolve_dimensional_values(inputVoltage, DimensionalValues::MAXIMUM) * std::sqrt(2);
            if (outputVoltage <= vinPeakMax) {
                if (assert) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "PFC output voltage must be greater than peak input voltage. Vout: " +
                        std::to_string(outputVoltage) + " <= Vin_peak_max: " + std::to_string(vinPeakMax));
                }
                valid = false;
            }
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
        const double vout = get_output_voltage();
        const double vd   = effective_diode_voltage_drop();
        if (is_buck_boost_class()) {
            // SEPIC / Ćuk (buck-boost conversion ratio):
            //   Vout/Vin = D/(1−D)  ⇒  D = (Vout+Vd)/(Vin+Vout+Vd).
            // Consistent with Sepic::calculate_conversion_ratio().
            return (vout + vd) / (vinPeak + vout + vd);
        }
        // Boost family: D = 1 − Vin/(Vout+Vd).
        // effective_diode_voltage_drop() centralises the TOTEM_POLE → 0 policy.
        return 1.0 - vinPeak / (vout + vd);
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
        // The DCM closed form below is the boost discontinuous-mode result
        // (L = Vin²·D² / (2·P·fsw)). For SEPIC/Ćuk, DCM involves the effective
        // inductance Le = L1‖L2 and a different energy relation, which is not
        // validated here — reject rather than return a boost number.
        if (is_buck_boost_class()) {
            throw std::runtime_error(
                "PowerFactorCorrection: DCM inductance sizing is not yet "
                "validated for SEPIC/Ćuk (requires the Le = L1‖L2 model). "
                "Use CONTINUOUS_CONDUCTION_MODE or CRITICAL_CONDUCTION_MODE.");
        }
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
        double D    = calculate_duty_cycle(vinPeak, get_output_voltage());
        double fsw  = get_switching_frequency();
        std::string mode = determine_actual_mode(inductance);

        if (mode == "Discontinuous Conduction Mode") {
            // DCM peak current is set by the energy delivered per switching
            // cycle, not by the average envelope: ipk_sw = Vin·D·Tsw / L.
            // The CCM "envelope + ΔI/2" formula does NOT apply — it would
            // under-report the peak by a large factor.
            return vinPeak * D / (inductance * fsw);
        }

        // CCM (and CrCM, which is the CCM/DCM boundary): peak = envelope avg
        // at line peak (√2·iinRms) + ΔI/2.
        double pinAvg = per_phase_power(*this) / get_efficiency_required();
        double iinRms = pinAvg / (vinPeak / std::sqrt(2));
        double iAvg   = iinRms * std::sqrt(2);
        double deltaI = vinPeak * D / (inductance * fsw);
        return iAvg + deltaI / 2;
    }

    std::string PowerFactorCorrection::determine_actual_mode(double inductance) {
        double L_critical = calculate_inductance_crcm();
        // ±5 % band around L_crit labels the result as Critical (boundary)
        // mode. The exact figure is a labelling convenience, not a physics
        // boundary — at L = L_crit the converter operates at the CCM/DCM
        // boundary by definition. The 5 % band absorbs floating-point /
        // tolerance noise from the sizing formulas (`calculate_inductance_*`
        // multiply/divide √2, π, η — a few % relative drift is normal) so
        // that a value sized as "CrCM" doesn't get re-labelled CCM/DCM
        // purely by rounding. Tightening below ~1 % causes spurious flips;
        // widening past ~10 % overlaps real CCM/DCM operating regions.
        // This band is consistent with the band used in TI SLUA479 fig. 6
        // (industry-conventional CrCM tolerance window).
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
            throw std::runtime_error(
                "PowerFactorCorrection::process_design_requirements: "
                "unsupported mode '" + mode + "'");
        }

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(inductance);
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        // Populate the computed-inductance / actual-mode diagnostics now so
        // the wizard can show them even on an analytical-only call (before
        // process_operating_points runs).
        computedInductance  = inductance;
        computedActualMode  = determine_actual_mode(inductance);

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
        // SEPIC/Ćuk: unipolar input current but a buck-boost OFF-time inductor
        // voltage (−(Vout+Vd)), handled in the waveform loop below.
        const bool buckBoostClass = is_buck_boost_class();

        // Calculate required inductance if not provided
        double L = magnetizingInductance;
        if (L <= 0) {
            if (mode == "Continuous Conduction Mode") {
                L = calculate_inductance_ccm();
            } else if (mode == "Discontinuous Conduction Mode") {
                L = calculate_inductance_dcm();
            } else if (mode == "Critical Conduction Mode" || mode == "Transition Mode") {
                L = calculate_inductance_crcm();
            } else {
                throw std::runtime_error(
                    "PowerFactorCorrection::process_operating_points: "
                    "unsupported mode '" + mode + "'");
            }
        }

        auto inputVoltage = get_input_voltage();
        double vinRmsMin  = resolve_dimensional_values(inputVoltage, DimensionalValues::MINIMUM);
        double vinPeakMin = vinRmsMin * std::sqrt(2);

        double outputVoltage     = get_output_voltage();
        double outputPower       = per_phase_power(*this);
        double switchingFrequency= get_switching_frequency();
        double lineFrequency     = get_line_frequency_required();
        double diodeVoltageDrop  = effective_diode_voltage_drop();
        double efficiency        = get_efficiency_required();
        double ambientTemperature= get_ambient_temperature();

        double pinAvg     = outputPower / efficiency;
        double iinRmsAvg  = pinAvg / vinRmsMin;
        double iLinePeak  = iinRmsAvg * std::sqrt(2);

        // Populate per-OP diagnostics. Worst-case is at minimum input voltage
        // peak-of-line (where current is max → ΔI_L is max and D is max).
        // For DCM the peak-current formula above already assumes continuous
        // envelope; the lastInductorRipple reflects the peak-of-line value.
        {
            double Dpeak    = calculate_duty_cycle(vinPeakMin, outputVoltage);
            double deltaIpk = vinPeakMin * Dpeak / (L * switchingFrequency);
            computedInductance      = L;
            computedActualMode      = mode;
            lastDutyCyclePeak       = Dpeak;
            lastPeakInductorCurrent = iLinePeak + deltaIpk / 2.0;
            lastInductorRipple      = deltaIpk;
            lastLineRmsCurrent      = iinRmsAvg;
            lastInputPower          = pinAvg;

            // Single-shot per-OP surface (one entry, worst case at min-Vin peak-of-line).
            perOpName.clear();
            perOpDutyCyclePeak.clear();
            perOpPeakInductorCurrent.clear();
            perOpInductorRipple.clear();
            perOpLineRmsCurrent.clear();
            perOpInputPower.clear();
            perOpName.push_back("Worst-case");
            perOpDutyCyclePeak.push_back(lastDutyCyclePeak);
            perOpPeakInductorCurrent.push_back(lastPeakInductorCurrent);
            perOpInductorRipple.push_back(lastInductorRipple);
            perOpLineRmsCurrent.push_back(lastLineRmsCurrent);
            perOpInputPower.push_back(lastInputPower);
        }

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
            const double vinInst    = vinPeakMin * vinShape;
            const double vinAbsInst = vinPeakMin * vinAbs;

            // Switching duty — TOPOLOGY-AWARE (must match
            // simulate_and_extract_waveforms and the OFF-time voltage selected
            // below). Boost family: D = 1 − Vin/(Vout+Vd). SEPIC/Ćuk: the
            // buck-boost ratio D = (Vout+Vd)/(Vin+Vout+Vd). Using the boost duty
            // for SEPIC/Ćuk would violate L1 volt-second balance
            // (Vin·D ≠ (Vout+Vd)·(1−D)) — see the OFF-time branch comment.
            //
            // Analytical-model boundary at line zero-crossings: as Vin→0 the
            // boost duty → 1 and ΔI = Vin·D·Tsw/L → 0; the model produces no
            // usable sample at exact Vin = 0 (a boost can't sustain output
            // current with zero input — it relies on Cbus dump-out, which this
            // analytical surface does not model). The physical duty bound
            // D ∈ [0, 1] holds; we clip there (NOT [0.01, 0.95]) because [0, 1]
            // is the true physical range. A real PFC controller simply disables
            // gating during the few µs of near-zero Vin per cycle.
            double D = buckBoostClass
                ? (outputVoltage + diodeVoltageDrop) / (vinAbsInst + outputVoltage + diodeVoltageDrop)
                : 1.0 - vinAbsInst / (outputVoltage + diodeVoltageDrop);
            if (D < 0.0) D = 0.0;
            if (D > 1.0) D = 1.0;

            double iAvgInst = iLinePeak * vinShape;     // signed for totem-pole
            double deltaI   = vinAbsInst * D / (L * switchingFrequency);

            // Integer-arithmetic switching-cycle phase (samples per switching
            // period = 4 since timeStep = Tsw/4). Using std::fmod(t, Tsw)
            // over ~80k samples accumulates float-rounding drift large
            // enough to shift the triangle's vertex by a fraction of a
            // sample late in the simulation.
            constexpr size_t samplesPerSwCycle = 4;
            double switchPhase = static_cast<double>(i % samplesPerSwCycle)
                                 / static_cast<double>(samplesPerSwCycle);

            double ripple;
            if (switchPhase < D) {
                ripple = -deltaI/2 + deltaI * (switchPhase / D);
            } else {
                ripple = deltaI/2 - deltaI * ((switchPhase - D) / (1 - D));
            }

            currentData.push_back(iAvgInst + ripple);

            if (switchPhase < D) {
                // ON-time: every supported variant's series inductor sees +Vin.
                voltageData.push_back(vinInst);
            } else if (buckBoostClass) {
                // SEPIC/Ćuk OFF-time: the input inductor L1 sees −(Vout+Vd),
                // independent of Vin (volt-second balance Vin·D = (Vout+Vd)·
                // (1−D) gives the buck-boost ratio). Unipolar input current,
                // so no sign flip.
                voltageData.push_back(-(outputVoltage + diodeVoltageDrop));
            } else {
                // Boost-family OFF-time: inductor sees Vin − Vout (signed); for
                // totem-pole's negative half-cycle the polarity of Vout flips.
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
        // PFC is a single-winding inductor: the inductor current IS the
        // magnetizing current. Set it explicitly so downstream code does
        // not try to re-derive it by integrating the switched voltage
        // waveform — that integration is sampled at only 4 points per
        // switching cycle over ~80k samples (line × switching scales),
        // which accumulates enough numerical drift to produce a
        // ~1000× overestimate of peak magnetizing current and breaks
        // every gap / energy calculation downstream.
        excitation.set_magnetizing_current(current);

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

        const double vinRms = get_vrms_for_simulation();
        const double vinPeak = vinRms * std::sqrt(2);
        double vout    = get_output_voltage();
        double pout    = get_output_power();
        double switchingFrequency = get_switching_frequency();
        double lineFrequency      = get_line_frequency_required();

        // Pin = Pout / η — the sizing path uses Pin (Erickson §17.3); use the
        // same here so the analytical and behavioural surfaces agree.
        const double pin = pout / get_efficiency_required();
        double iPeak = std::sqrt(2) * pin / vinRms;
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
        const bool buckBoostClass = is_buck_boost_class();
        const double vd = effective_diode_voltage_drop();
        double switchingFrequency = get_switching_frequency();
        double lineFrequency      = get_line_frequency_required();
        waveforms.switchingFrequency = switchingFrequency;
        waveforms.lineFrequency      = lineFrequency;

        const double vinRms = get_vrms_for_simulation();
        const double vinPeak = vinRms * std::sqrt(2);
        double vout    = get_output_voltage();
        // Per-phase power for INTERLEAVED_BOOST.
        double pout    = per_phase_power(*this);

        // Pin = Pout / η — sizing uses Pin, so the envelope amplitude here
        // must too. Using Pout would under-represent the line current by 1/η.
        const double pin = pout / get_efficiency_required();
        double iPeak = std::sqrt(2) * pin / vinRms;

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

            // Duty is bounded to its physical range [0, 1]. Near the line
            // zero-crossing the analytical model degrades (D → 1, ripple → 0);
            // a real controller stops gating there. See note in
            // process_operating_points for the full rationale.
            // Topology-aware: SEPIC/Ćuk use the buck-boost ratio, boost family
            // the boost ratio. (The series inductor sees +Vin during ON in
            // both, so the ripple amplitude below is the same form.)
            const double vinInst = vinPeak * absShape;
            double duty = buckBoostClass
                ? (vout + vd) / (vinInst + vout + vd)
                : 1.0 - vinInst / (vout + vd);
            if (duty < 0) duty = 0;
            if (duty > 1) duty = 1;

            double rippleAmplitude = vinInst * duty / (2.0 * inductance * switchingFrequency);
            waveforms.currentRipple[i] = rippleAmplitude;

            // Integer-arithmetic phase (timeStep = Tsw/4 → 4 samples/cycle).
            constexpr size_t samplesPerSwCycle = 4;
            double sawtoothPhase = static_cast<double>(i % samplesPerSwCycle)
                                   / static_cast<double>(samplesPerSwCycle);
            double triangular;
            if (sawtoothPhase < 0.5) {
                triangular = -1.0 + 4.0 * sawtoothPhase;
            } else {
                triangular = 3.0 - 4.0 * sawtoothPhase;
            }

            waveforms.inductorCurrent[i] = iEnv + rippleAmplitude * triangular;

            if (sawtoothPhase < 0.5) {
                // ON-time: series inductor sees +Vin.
                waveforms.inductorVoltage[i] = vin;
            } else if (buckBoostClass) {
                // SEPIC/Ćuk OFF-time: L1 sees −(Vout+Vd), independent of Vin.
                waveforms.inductorVoltage[i] = -(vout + vd);
            } else {
                // Boost-family OFF-time = Vin − Vout (signed for totem-pole).
                const double voutSigned = bipolar ? std::copysign(vout, shape) : vout;
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
        // Reflect the spec-level efficiency used everywhere else in the
        // sizing path (Pin = Pout/η). The analytical synthesis does not
        // model loss explicitly, but reporting 1.0 here would silently
        // contradict the rest of the model — use the spec value.
        waveforms.efficiency  = get_efficiency_required();
        waveforms.operatingPointName = "PFC_analytical";

        return waveforms;
    }

    std::vector<OperatingPoint> PowerFactorCorrection::simulate_and_extract_operating_points(double inductance,
                                                                                               double dcResistance) {
        std::vector<OperatingPoint> operatingPoints;

        int numPeriods = (_numberOfPeriods > 0) ? _numberOfPeriods : 2;
        PfcSimulationWaveforms waveforms = simulate_and_extract_waveforms(inductance, dcResistance, numPeriods);

        // Per CLAUDE.md: do not silently fall back to the analytical path
        // when the synthesis returned no samples. That hides a bug.
        if (waveforms.inductorCurrent.empty() || waveforms.time.empty()) {
            throw std::runtime_error(
                "PowerFactorCorrection::simulate_and_extract_operating_points: "
                "simulate_and_extract_waveforms returned empty data (no samples)");
        }

        double switchingFrequency = get_switching_frequency();
        double lineFrequency      = get_line_frequency_required();
        double ambientTemperature = get_ambient_temperature();

        Waveform currentWaveform;
        currentWaveform.set_time(waveforms.time);
        currentWaveform.set_data(waveforms.inductorCurrent);

        SignalDescriptor current;
        current.set_waveform(currentWaveform);

        // No try/catch swallow: if harmonics extraction fails, propagate the
        // exception so callers see the real problem rather than a silently-
        // analytical-only Processed.
        auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
        current.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency));
        auto currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true);
        currentProcessed.set_label(WaveformLabel::CUSTOM);
        current.set_processed(currentProcessed);

        SignalDescriptor voltage;
        if (!waveforms.inductorVoltage.empty() && waveforms.inductorVoltage.size() == waveforms.time.size()) {
            Waveform voltageWaveform;
            voltageWaveform.set_time(waveforms.time);
            voltageWaveform.set_data(waveforms.inductorVoltage);
            voltage.set_waveform(voltageWaveform);

            auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
            voltage.set_harmonics(Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency));
            auto voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true);
            voltageProcessed.set_label(WaveformLabel::CUSTOM);
            voltage.set_processed(voltageProcessed);
        }

        OperatingPointExcitation excitation;
        excitation.set_current(current);
        excitation.set_frequency(lineFrequency);
        excitation.set_voltage(voltage);
        // See note in process_operating_points (analytical path): PFC is a
        // single-winding inductor, so the inductor current IS the
        // magnetizing current. Setting it explicitly avoids the integrate-
        // voltage fallback in pre_process_inputs producing a wildly wrong
        // peak from the sparsely-sampled switched waveform.
        excitation.set_magnetizing_current(current);

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
                // Physical bound only — D ∈ [0, 1]. Near line zero-crossing
                // the analytical model degenerates; see the rationale note
                // in process_operating_points. Topology-aware duty: SEPIC/Ćuk
                // use the buck-boost ratio (Vd folded into vout here, which is
                // 0 for the §5.1 ideal-port view).
                double duty = is_buck_boost_class()
                    ? vout / (vinInst + vout)
                    : 1.0 - vinInst / vout;
                if (duty < 0.0) duty = 0.0;
                if (duty > 1.0) duty = 1.0;
                const double rippleAmp = vinInst * duty / (2.0 * inductance * switchingFrequency);
                // Integer-arithmetic phase (timeStep = Tsw/4 → 4 samples/cycle).
                constexpr size_t samplesPerSwCycle = 4;
                const double phase = static_cast<double>(i % samplesPerSwCycle)
                                     / static_cast<double>(samplesPerSwCycle);
                const double tri   = (phase < 0.5) ? (-1.0 + 4.0 * phase)
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
        // SPICE-side knobs from spice_config(). PFC registry entry
        // (Topology.cpp) matches this file's historical hardcoded
        // netlist byte-for-byte — behaviour-preserving.
        const auto cfg = spice_config();

        validate_topology_variant();
        const auto variant = get_topology_variant_or_default();
        // The native average-current-mode controller drives a unipolar boost
        // power stage. Boost, bridgeless and semi-bridgeless all present the
        // identical rectified-sine inductor current, so they share this
        // netlist verbatim. Interleaved boost is N independent boost cells;
        // we emit the netlist for ONE phase (per-phase power + per-phase
        // inductance, both already returned by process_design_requirements),
        // which is what the per-phase magnetic must be verified against.
        // TOTEM_POLE takes the genuine bipolar branch below: its inductor sits
        // on the AC side and carries BIPOLAR current, so it gets a floating
        // line source + a 4-switch totem-pole power stage (HF leg with the
        // active boost switch + rectifying diode role-swapped by line polarity,
        // and a line-frequency LF polarity leg). The average-current controller
        // runs in the rectified domain (senses |i_L|, multiplier on |vline|) and
        // the PWM is routed to the active HF switch by polarity — exactly how
        // industry totem-pole controllers operate. Only one active switch + one
        // diode conduct at any instant, so it converges like the boost.
        switch (variant) {
            case PfcTopologyVariants::BOOST:
            case PfcTopologyVariants::BRIDGELESS:
            case PfcTopologyVariants::SEMI_BRIDGELESS:
            case PfcTopologyVariants::INTERLEAVED_BOOST:
            case PfcTopologyVariants::TOTEM_POLE:
                break;  // supported (boost family unipolar; totem-pole bipolar)
            default:
                throw std::runtime_error(
                    "PowerFactorCorrection::generate_ngspice_switching_circuit: "
                    "only the boost family (BOOST, BRIDGELESS, SEMI_BRIDGELESS, "
                    "INTERLEAVED_BOOST) and TOTEM_POLE are supported by the "
                    "native switching controller. SEPIC/CUK/BUCK/BUCK_BOOST/"
                    "VIENNA need their own switching netlists (out of scope here).");
        }
        // Bipolar totem-pole vs unipolar boost family. Selects the power stage,
        // the current-sense rectification, the PWM routing, and the saved
        // bipolar line node below.
        const bool bipolar = (variant == PfcTopologyVariants::TOTEM_POLE);
        if (!std::isfinite(inductance) || inductance <= 0.0) {
            throw std::invalid_argument(
                "PowerFactorCorrection::generate_ngspice_switching_circuit: "
                "inductance must be positive and finite (got " +
                std::to_string(inductance) + ")");
        }
        if (numberOfLineCycles < 1) numberOfLineCycles = 1;

        // ── Spec resolution ─────────────────────────────────────────────
        const double vrmsNom = get_vrms_for_simulation();
        const double vbusNom = get_output_voltage();
        // Per-phase power: Pout for boost/bridgeless/semi-bridgeless, Pout/N
        // for interleaved (one cell carries 1/N of the load). The inductance
        // argument is likewise the per-phase value, so the emitted single-cell
        // boost netlist is self-consistent.
        const double pout    = per_phase_power(*this);
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
        const char* variantLabel = "Boost";
        switch (variant) {
            case PfcTopologyVariants::BRIDGELESS:        variantLabel = "Bridgeless (boost-equivalent inductor)"; break;
            case PfcTopologyVariants::SEMI_BRIDGELESS:   variantLabel = "Semi-bridgeless (boost-equivalent inductor)"; break;
            case PfcTopologyVariants::INTERLEAVED_BOOST: variantLabel = "Interleaved boost (single per-phase cell)"; break;
            case PfcTopologyVariants::TOTEM_POLE:        variantLabel = "Totem-pole (bipolar 4-switch stage)"; break;
            default: break;
        }
        c << "* PFC " << variantLabel << " Switching Converter (native average-current-mode controller)\n";
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
        c << ".param ic_vea="           << sci(t.ic_vea)         << "\n";
        c << ".param ic_vc_i="          << sci(t.ic_vc_i)        << "\n";
        c << ".param t_ss_release="     << sci(t.t_ss_release)   << "\n\n";

        // ── Subckt prelude (OPAMP_IDEAL + COMPARATOR_IDEAL) ─────────────
        c << "* ─── Reusable subcircuits ──────────────────────────────────\n";
        c << spice_subckt_prelude_pfc_controller();
        c << "\n";

        // ── Power stage ─────────────────────────────────────────────────
        // SW and diode models are shared by both branches.
        const auto emitSwModel = [&]() {
            c << ".model SW1 SW (VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
              << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << ")\n";
        };
        const auto emitDiodeModel = [&]() {
            c << ".model DIDEAL D (IS=" << std::scientific << cfg.diodeIS
              << " RS=" << cfg.diodeRS << std::defaultfloat;
            if (!cfg.diodeExtra.empty()) c << " " << cfg.diodeExtra;
            c << ")\n";
        };
        if (!bipolar) {
            c << "* ─── Power stage (Boost) ───────────────────────────────────\n";
            c << "B_vin     vin_rect 0     V=vpk*abs(sin(2*3.141592653589793*fline*time))\n";
            c << "Vl_sense  vin_rect l_in  0          ; in-line current sense\n";
            c << "L1        l_in     sw    {l_ind}\n";
            emitSwModel();
            c << "S1        sw  0    gate 0 SW1\n";
            // Series-RC snubber across the switch (R in series with C, NOT both
            // shunted to 0).  When S1 is OFF and D1 is reverse-biased (vin_rect <
            // vbus), the inductor current must reach 0 before the cycle ends; if
            // the snubber R were directly to ground (the "parallel-RC" pattern
            // used elsewhere in this codebase) it would carry vpk/Rsnub of DC
            // current at every line-cycle peak — for Rsnub=100 Ω that's >3 A,
            // dwarfing the 0.6 A i_ref and locking the current EA into hard
            // anti-windup at vc_i = 0.  In a real boost the inductor *can* fully
            // demagnetize each cycle (DCM near the line zero-crossings); in CCM
            // the diode keeps conducting into vbus.  The series-RC pattern blocks
            // DC and only damps switching-edge ringing, which is the textbook
            // role of an RC snubber (Erickson §A.2; Basso 2008 §4.5).
            c << "Rsnub_s1  sw   snub_n  " << cfg.snubR << "\n";
            c << "Csnub_s1  snub_n 0     " << std::scientific << cfg.snubC
              << std::defaultfloat << "\n";
            emitDiodeModel();
            c << "D1        sw  vbus DIDEAL\n";
            c << "Cout      vbus 0   {cbus}  IC={ic_vbus}\n";
            c << "Rload     vbus 0   {rload}\n\n";
        } else {
            // Genuine bipolar totem-pole. The boost inductor is in series with
            // the AC line (no input bridge), so it carries a SIGNED sinusoidal
            // current. Two half-bridge legs share the bus:
            //   • HF leg (mid_hf): switched at fsw. In each line half-cycle ONE
            //     of its switches is the active boost switch and the OPPOSITE
            //     leg's diode rectifies into the bus — so only one switch + one
            //     diode conduct at a time, just like the boost. Positive half:
            //     S_hf_lo is the boost switch, D_hf_hi rectifies (mid_hf>vbus).
            //     Negative half: S_hf_hi is the boost switch, D_hf_lo rectifies
            //     (mid_hf<0).
            //   • LF leg (mid_lf = the AC return): switched at the LINE rate,
            //     tying the return to gnd (positive half) or vbus (negative
            //     half). Exactly one LF switch is ON at all times (no overlap,
            //     no gap), so mid_lf is always hard-driven.
            // vin_rect / vline are ground-referenced copies used only by the
            // controller (feed-forward, multiplier, polarity); the power path
            // is the FLOATING source B_vac between vac and mid_lf.
            c << "* ─── Power stage (Totem-pole, bipolar) ─────────────────────\n";
            c << "B_vline   vline 0      V=vpk*sin(2*3.141592653589793*fline*time)\n";
            c << "B_vin_rect vin_rect 0  V=abs(V(vline))   ; rectified line for control\n";
            c << "B_vac     vac mid_lf   V=vpk*sin(2*3.141592653589793*fline*time)\n";
            c << "Vl_sense  vac  l_in    0          ; signed in-line current sense\n";
            c << "L1        l_in mid_hf  {l_ind}\n";
            emitSwModel();
            emitDiodeModel();
            c << "S_hf_hi   mid_hf vbus  g_hi 0 SW1\n";
            c << "S_hf_lo   mid_hf 0     g_lo 0 SW1\n";
            c << "D_hf_hi   mid_hf vbus  DIDEAL    ; rectifies positive-half boost\n";
            c << "D_hf_lo   0 mid_hf     DIDEAL    ; rectifies negative-half boost\n";
            c << "S_lf_hi   mid_lf vbus  g_lf_hi 0 SW1\n";
            c << "S_lf_lo   mid_lf 0     g_lf_lo 0 SW1\n";
            // Series-RC snubber on the HF midpoint (same role/rationale as the
            // boost snubber above — blocks DC, damps switching-edge ringing).
            c << "Rsnub_s1  mid_hf snub_n  " << cfg.snubR << "\n";
            c << "Csnub_s1  snub_n 0     " << std::scientific << cfg.snubC
              << std::defaultfloat << "\n";
            c << "Cout      vbus 0   {cbus}  IC={ic_vbus}\n";
            c << "Rload     vbus 0   {rload}\n\n";
        }

        // ── Voltage error amp (Block 1) ─────────────────────────────────
        c << "* ─── Voltage error amp (Gv) — type-II ─────────────────────\n";
        c << "B_div     vbus_div 0  V=V(vbus)*k_div\n";
        // Soft-start: ramp the EA reference linearly from
        //   vref(0)  = Vpk · k_div   (= ic_vbus · k_div, so error is 0 at t=0+)
        //   vref(t)  = vref          (= 2.5 V) for t ≥ t_ss_release
        // This grows the bus-voltage target gradually from the bridge
        // equilibrium (Vpk) up to Vbus_nom, which keeps the inrush
        // current that flows into Cbus during the boost-up bounded.
        // Without the ramp, the EA must close from a 75-V target gap,
        // pumping >3 A through L1 each line peak — the current EA reads
        // i_sense >> i_ref and winds C_fb_zi negative within the first
        // line cycle, dead-locking vc_i at 0 forever (gate OFF).
        // Per Basso 2012 Ch. 6 + TI SLUA479 (UCC28070 SS architecture).
        c << "B_vref_v  vref_v   0  V=(vpk*k_div)+(vref-vpk*k_div)*min(1, time/t_ss_release)\n";
        c << "R_in_v    vbus_div vea_n  {gv_rin}\n";
        c << "R_fb_z    vea_n    vea_z  {gv_rz}\n";
        c << "C_fb_z    vea_z    vea    {gv_cz} IC=0\n";
        c << "C_fb_p    vea_n    vea    {gv_cp} IC=0\n";
        c << "X_op_v    vref_v   vea_n  vea OPAMP_IDEAL "
             "params: A0=1e3 GBW=10e6 VSSPOS={vea_max} VSSNEG={vea_min} IC_FILT={ic_vea}\n\n";

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
        // Totem-pole inductor current is bipolar; the average-current loop
        // regulates its MAGNITUDE (the LF leg + HF role-swap handle polarity),
        // so sense |i_L|. Boost current is already unipolar.
        if (bipolar)
            c << "B_isense  i_sense  0  V=abs(I(Vl_sense))*rs_sense\n";
        else
            c << "B_isense  i_sense  0  V=I(Vl_sense)*rs_sense\n";
        c << "R_in_i    i_sense  ic_n   {gi_rin}\n";
        c << "R_fb_zi   ic_n     ic_z   {gi_rz}\n";
        // Pre-bias the integrator capacitor C_fb_zi to ic_vc_i. ic_z sits at
        // virtual ground (≈ V(i_ref) ≈ 0 at startup), and V_cap = V(ic_z) -
        // V(vc_i). For vc_i to start at +ic_vc_i, V_cap must start at
        // -ic_vc_i. Without this pre-bias, the integrator state is washed
        // out within ~1 µs (the opamp's internal R_pole·C_pole) and the
        // current EA cannot remember its operating point through the
        // first inrush event of the line cycle (which would otherwise
        // wind C_fb_zi irreversibly negative — gate dead-locked OFF).
        c << "C_fb_zi   ic_z     vc_i   {gi_cz} IC=-{ic_vc_i}\n";
        c << "C_fb_pi   ic_n     vc_i   {gi_cp} IC=0\n";
        c << "X_op_i    i_ref    ic_n   vc_i OPAMP_IDEAL "
             // Current EA: fc_i = fsw/10 — internal opamp pole must sit
             // ≥ 10·fc_i. A0=1e3, GBW=100 MHz → fp = 100 kHz. IC_FILT
             // intentionally NOT set on the current EA — it is washed
             // out in 1.6 µs by the natural feedback dynamics; the real
             // integrator state lives on C_fb_zi above (see IC=-{ic_vc_i}).
             "params: A0=1e3 GBW=100e6 VSSPOS={v_pk_saw} VSSNEG=0\n\n";

        // ── Oscillator + PWM comparator (Blocks 5–6) ────────────────────
        c << "* ─── Sawtooth + PWM comparator ────────────────────────────\n";
        c << "V_saw     saw 0   PULSE(0 {v_pk_saw} 0 {tsw-pwm_t_rise} "
             "{pwm_t_rise} {pwm_t_rise} {tsw})\n";
        if (!bipolar) {
            c << "B_gate    gate    0  V=(V(vc_i) > V(saw)) ? v_high : 0\n\n";
        } else {
            // Bipolar totem-pole modulator. The current loop produces ONE PWM
            // boost command (vc_i vs saw); a line-polarity selector routes it to
            // whichever HF switch is the active boost switch this half-cycle, and
            // drives the LF leg to tie the AC return to the correct rail. pol_pos
            // = 1 on the positive line half-cycle, 0 on the negative.
            //   positive half: S_hf_lo = PWM (boost), D_hf_hi rectifies, LF→gnd
            //   negative half: S_hf_hi = PWM (boost), D_hf_lo rectifies, LF→vbus
            c << "B_polpos  pol_pos 0  V=(V(vline) >= 0) ? 1 : 0\n";
            c << "B_pwm     pwm     0  V=(V(vc_i) > V(saw)) ? 1 : 0\n";
            c << "B_g_lo    g_lo    0  V=V(pol_pos)*V(pwm)*v_high\n";
            c << "B_g_hi    g_hi    0  V=(1-V(pol_pos))*V(pwm)*v_high\n";
            c << "B_g_lf_lo g_lf_lo 0  V=V(pol_pos)*v_high\n";
            c << "B_g_lf_hi g_lf_hi 0  V=(1-V(pol_pos))*v_high\n\n";
        }

        // ── Initial conditions ──────────────────────────────────────────
        // Note: vea / vc_i / vrms_ff are B-source outputs and CANNOT be
        // ICed directly — `.ic v(vea)=...` is silently dropped by ngspice.
        // The opamp warm-start is applied via the IC_FILT subckt parameter
        // on the X_op_v / X_op_i instances above (sets the internal `filt`
        // capacitor — that IS a state node). vrms_ff = vff2² follows
        // automatically from vff1 / vff2 capacitor ICs.
        c << "* ─── Initial conditions (warm start to analytical SS) ────\n";
        c << ".ic v(vbus)={ic_vbus} v(vff1)={ic_vff} v(vff2)={ic_vff}\n\n";

        // ── Transient analysis + solver options ─────────────────────────
        const double simTime = numberOfLineCycles * tline;
        const double maxStep = tsw / 40.0;
        c << "* ─── Analysis ──────────────────────────────────────────────\n";
        c << ".tran " << sci(maxStep) << " " << sci(simTime) << " 0 "
          << sci(maxStep) << " uic\n";
        if (!bipolar) {
            c << ".save i(vl_sense) v(vin_rect) v(vbus) v(vea) v(vrms_ff) "
                 "v(i_ref) v(i_sense) v(vc_i) v(vref_v) v(saw) v(gate)\n";
        } else {
            // Bipolar: also save v(vline) (the signed AC line — the genuine
            // input-port voltage and the signed power-balance partner of the
            // bipolar i(vl_sense)) and the routed gate signals.
            c << ".save i(vl_sense) v(vline) v(vin_rect) v(vbus) v(vea) v(vrms_ff) "
                 "v(i_ref) v(i_sense) v(vc_i) v(vref_v) v(saw) v(pwm) "
                 "v(g_lo) v(g_hi)\n";
        }
        // std::defaultfloat after std::scientific block — see the
        // IsolatedBuck commit (6f795fef) for why std::fixed would break.
        c << ".options METHOD=" << cfg.method << " TRTOL=" << cfg.trTol
          << " RELTOL=" << cfg.relTol
          << " ABSTOL=" << std::scientific << cfg.absTol
          << " VNTOL=" << cfg.vnTol << std::defaultfloat << "\n";
        c << ".options ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
        c << ".end\n";

        return c.str();
    }

    OperatingPoint PowerFactorCorrection::simulate_with_ngspice_switching(
            double inductance,
            int numberOfLineCycles,
            bool trimToLastLineCycle) {
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

        // ── Optional trim to last full line cycle ───────────────────────
        // The simulator always returns the full transient (soft-start +
        // any number of line cycles).  Downstream consumers — Painter's
        // composite SVG renderer, THD/PF analysers, PtP comparators —
        // almost always want only the *steady-state* portion: one full
        // line period, starting at a zero-crossing so the half-sine
        // envelope is complete.  When trimToLastLineCycle is true (the
        // default) we slice every saved waveform to t ∈ [t_end − Tline,
        // t_end] before packing it into the OperatingPoint.  Set the
        // parameter to false to retain the full simulation history (the
        // diagnostic harness in TestPfcPhase4Explore.cpp does this so it
        // can also report startup metrics).
        const double tline = 1.0 / get_line_frequency_required();
        std::vector<Waveform> trimmedWaveforms;
        std::vector<std::string> trimmedNames;
        if (trimToLastLineCycle) {
            // Pick the time vector from the inductor-current waveform —
            // every saved signal shares the same simulator time grid.
            auto tOpt = iLwf->get_time();
            if (!tOpt.has_value() || tOpt->empty()) {
                throw std::runtime_error(
                    "PowerFactorCorrection::simulate_with_ngspice_switching:"
                    " inductor-current waveform has no time vector — "
                    "cannot trim to last line cycle");
            }
            const auto& tFull = tOpt.value();
            const double tEnd = tFull.back();
            const double tBeg = std::max(tFull.front(), tEnd - tline);
            // Find first sample with t >= tBeg.
            auto itBeg = std::lower_bound(tFull.begin(), tFull.end(), tBeg);
            const size_t iBeg = std::distance(tFull.begin(), itBeg);

            trimmedWaveforms.reserve(sim.waveforms.size());
            trimmedNames.reserve(sim.waveforms.size());
            for (size_t k = 0; k < sim.waveforms.size(); ++k) {
                const auto& src = sim.waveforms[k];
                auto srcT = src.get_time();
                const auto& srcD = src.get_data();
                // ngspice's NgspiceRunner only attaches a time vector to the
                // dedicated `time` waveform — every other signal shares the
                // same grid implicitly.  Broadcast `tFull` to every data
                // signal so downstream consumers (Painter, ConverterPort
                // Checks::check_pfc_switching_ports) can pull `.get_time()`
                // off any waveform.
                Waveform w;
                if (srcT.has_value() && srcT->size() == srcD.size()) {
                    std::vector<double> tNew(srcT->begin() + iBeg, srcT->end());
                    std::vector<double> dNew(srcD.begin()  + iBeg, srcD.end());
                    w.set_time(tNew);
                    w.set_data(dNew);
                } else if (srcD.size() == tFull.size()) {
                    std::vector<double> tNew(tFull.begin() + iBeg, tFull.end());
                    std::vector<double> dNew(srcD.begin()  + iBeg, srcD.end());
                    w.set_time(tNew);
                    w.set_data(dNew);
                } else {
                    // Truly mismatched — keep as-is rather than fabricating.
                    w = src;
                }
                trimmedWaveforms.push_back(std::move(w));
                trimmedNames.push_back(sim.waveformNames[k]);
            }
        } else {
            // Untrimmed path: still need the per-signal time vector for
            // downstream consumers — broadcast the time vector from any
            // signal that carries one (typically the dedicated `time`
            // waveform).
            std::vector<double> tShared;
            for (const auto& wf : sim.waveforms) {
                auto t = wf.get_time();
                if (t.has_value() && !t->empty()) {
                    tShared = t.value();
                    break;
                }
            }
            trimmedWaveforms.reserve(sim.waveforms.size());
            trimmedNames.reserve(sim.waveforms.size());
            for (size_t k = 0; k < sim.waveforms.size(); ++k) {
                Waveform w = sim.waveforms[k];
                auto t = w.get_time();
                if ((!t.has_value() || t->empty()) && !tShared.empty()
                    && w.get_data().size() == tShared.size()) {
                    w.set_time(tShared);
                }
                trimmedWaveforms.push_back(std::move(w));
                trimmedNames.push_back(sim.waveformNames[k]);
            }
        }

        // Re-resolve iL after trimming so the OperatingPoint and downstream
        // signal lookups all share the same (trimmed or full) time axis.
        const Waveform* iLtrim = nullptr;
        for (size_t k = 0; k < trimmedNames.size(); ++k) {
            if (matchesILSense(trimmedNames[k])) {
                iLtrim = &trimmedWaveforms[k];
                break;
            }
        }

        // Wrap into an OperatingPoint with four diagnostic "windings" so
        // that downstream consumers (e.g. Painter::paint_operating_point_
        // waveforms, ConverterPortChecks::check_pfc_switching_ports) can
        // render every loop signal in a single composite SVG and apply the
        // standard converter-port DC-stream gates.  ngspice's saved nodes
        // are mapped onto the OperatingPointExcitation V/I pair as follows:
        //
        //   Winding 0  InputPort       V = vin_rect   I = i(vl_sense)
        //   Winding 1  PowerStage      V = vbus       I = i(vl_sense)
        //   Winding 2  VoltageLoop     V = vea        I = i_ref
        //   Winding 3  CurrentLoop     V = vc_i       I = i_sense
        //
        // (The V/I labelling on the controller windings is diagrammatic —
        // vea / vc_i are controller signals, not winding voltages — but it
        // lets the existing OperatingPoint plumbing carry them through
        // unchanged.)  The InputPort winding is what the PFC port-check
        // helper consumes to validate the rectified-line input port (mean
        // ≈ 0.9·Vrms, RMS ≈ Vrms).
        auto findByName = [&](const std::string& key) -> const Waveform* {
            std::string k = key;
            std::transform(k.begin(), k.end(), k.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            for (size_t i = 0; i < trimmedNames.size(); ++i) {
                std::string n = trimmedNames[i];
                std::transform(n.begin(), n.end(), n.begin(),
                               [](unsigned char c){ return std::tolower(c); });
                if (n.find(k) != std::string::npos) return &trimmedWaveforms[i];
            }
            return nullptr;
        };
        auto makeExc = [&](const std::string& name,
                           const Waveform* vWf,
                           const Waveform* iWf) {
            OperatingPointExcitation e;
            e.set_name(name);
            if (vWf) {
                SignalDescriptor s;
                s.set_waveform(*vWf);
                e.set_voltage(s);
            }
            if (iWf) {
                SignalDescriptor s;
                s.set_waveform(*iWf);
                e.set_current(s);
            }
            return e;
        };

        // Totem-pole is bridgeless: its genuine input-port voltage is the
        // SIGNED AC line (v(vline), saved only for that variant), in phase with
        // the bipolar inductor current. The bridged boost family presents the
        // rectified |line| (v(vin_rect)). Pick whichever the netlist saved.
        const bool bipolarInput =
            (get_topology_variant_or_default() == PfcTopologyVariants::TOTEM_POLE);
        const Waveform* inputVoltageWf =
            bipolarInput ? findByName("vline") : findByName("vin_rect");
        OperatingPoint op;
        op.get_mutable_excitations_per_winding().push_back(
            makeExc("InputPort",   inputVoltageWf, iLtrim));
        op.get_mutable_excitations_per_winding().push_back(
            makeExc("PowerStage",  findByName("vbus"), iLtrim));
        op.get_mutable_excitations_per_winding().push_back(
            makeExc("VoltageLoop", findByName("vea"),  findByName("i_ref")));
        op.get_mutable_excitations_per_winding().push_back(
            makeExc("CurrentLoop", findByName("vc_i"), findByName("i_sense")));
        return op;
    }

} // namespace OpenMagnetics
