#include "converter_models/Cllc.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <complex>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    // =========================================================================
    // Construction
    // =========================================================================

    CllcConverter::CllcConverter(const json& j) {
        // Parse base CllcResonant fields
        set_bidirectional(get_stack_optional<bool>(j, "bidirectional"));
        set_efficiency(get_stack_optional<double>(j, "efficiency"));
        set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
        set_max_switching_frequency(j.at("maxSwitchingFrequency").get<double>());
        set_min_switching_frequency(j.at("minSwitchingFrequency").get<double>());
        set_operating_points(j.at("operatingPoints").get<std::vector<CllcOperatingPoint>>());
        set_quality_factor(get_stack_optional<double>(j, "qualityFactor"));
        set_symmetric_design(get_stack_optional<bool>(j, "symmetricDesign"));
    }

    AdvancedCllcConverter::AdvancedCllcConverter(const json& j) {
        from_json(j, *this);
    }

    bool CllcConverter::run_checks(bool assert) {
        return Topology::run_checks_common(this, assert);
    }

    // =========================================================================
    // Resonant Tank Parameter Calculation
    // =========================================================================

    /**
     * Implements the design procedure from Infineon AN-2024-06, Section 2.3.
     *
     * The CLLC resonant tank is defined by:
     *   Primary side:   C1 (series capacitor) + L1 (series inductor, = leakage inductance)
     *   Transformer:    Lm (magnetizing inductance)
     *   Secondary side: L2 (series inductor) + C2 (series capacitor)
     *
     * The resonant frequency is: fr = 1/(2π√(L1·C1)) = 1/(2π√(L2·C2))
     *
     * For a symmetric tank: L2 = L1/n², C2 = n²·C1  (so fr1 = fr2)
     */
    CllcResonantParameters CllcConverter::calculate_resonant_parameters() {
        CllcResonantParameters params;

        double nominalInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::NOMINAL);
        if (nominalInputVoltage == 0) {
            nominalInputVoltage = (resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM) +
                                   resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM)) / 2.0;
        }

        // Find representative output voltage and power from operating points
        // Use the first forward-mode operating point, or the first one available
        double nominalOutputVoltage = 0;
        double nominalOutputPower = 0;
        double representativeFrequency = 0;

        for (const auto& op : get_operating_points()) {
            if (op.get_output_voltages().size() > 0) {
                nominalOutputVoltage = op.get_output_voltages()[0];
                double current = op.get_output_currents().size() > 0 ? op.get_output_currents()[0] : 0;
                nominalOutputPower = nominalOutputVoltage * current;
                representativeFrequency = op.get_switching_frequency();
                break;
            }
        }

        if (nominalOutputVoltage == 0 || nominalOutputPower == 0) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS,
                "CLLC: Cannot determine output voltage/power from operating points");
        }

        // Step 1: Transformer turns ratio (Infineon AN Step 1)
        // n = Vin_nominal / Vout_nominal
        params.turnsRatio = nominalInputVoltage / nominalOutputVoltage;

        // Use a representative resonant frequency:
        // If operating point frequency is within min/max range, use it;
        // otherwise use the geometric mean of min and max switching frequency
        double fr;
        if (representativeFrequency >= get_min_switching_frequency() &&
            representativeFrequency <= get_max_switching_frequency()) {
            fr = representativeFrequency;
        } else {
            fr = sqrt(get_min_switching_frequency() * get_max_switching_frequency());
        }
        params.resonantFrequency = fr;

        // Step 4: Effective AC resistance (FHA) (Infineon AN Step 4)
        // Ro = (8·n²/π²) · (Vout²/Pout) = (8·n²/π²) · R_load
        double n = params.turnsRatio;
        double Rload = nominalOutputVoltage * nominalOutputVoltage / nominalOutputPower;
        double Ro = (8.0 * n * n / (M_PI * M_PI)) * Rload;
        params.equivalentAcResistance = Ro;

        // Step 3: Quality factor Q (use user-specified or default)
        double Q;
        if (get_quality_factor()) {
            Q = get_quality_factor().value();
        } else {
            Q = 0.3;  // Default from Infineon AN recommendation (between 0.2 and 0.4)
        }
        params.qualityFactor = Q;

        // Step 5: Primary resonant capacitor C1 (Infineon AN Step 5)
        // C1 = 1 / (2π·Q·fr·Ro)
        double C1 = 1.0 / (2.0 * M_PI * Q * fr * Ro);
        params.primaryResonantCapacitance = C1;

        // Step 6: Primary resonant inductor L1 (Infineon AN Step 6)
        // L1 = 1 / ((2πfr)²·C1)
        double omega_r = 2.0 * M_PI * fr;
        double L1 = 1.0 / (omega_r * omega_r * C1);
        params.primaryResonantInductance = L1;

        // Step 7: Magnetizing inductance Lm (Infineon AN Step 7)
        // Lm = k · L1
        double k = defaultInductanceRatio;
        double Lm = k * L1;
        params.magnetizingInductance = Lm;
        params.inductanceRatio = k;

        // Steps 8-9: Secondary resonant components
        bool symmetric = get_symmetric_design().value_or(true);
        double a, b;
        if (symmetric) {
            // Symmetric design: a = 1, b = 1
            a = 1.0;
            b = 1.0;
        } else {
            // Asymmetric design: typical values from Infineon AN Step 8
            a = 0.95;
            b = 1.052;
        }
        params.resonantInductorRatio = a;
        params.resonantCapacitorRatio = b;

        // L2 = a·L1/n² (Infineon AN Step 9)
        params.secondaryResonantInductance = a * L1 / (n * n);

        // C2 = n²·b·C1 (Infineon AN Step 9)
        params.secondaryResonantCapacitance = n * n * b * C1;

        return params;
    }

    // =========================================================================
    // FHA Voltage Gain Calculation
    // =========================================================================

    /**
     * Compute the FHA voltage gain |nVout/Vin| at a given switching frequency.
     *
     * From Infineon AN Equation 2 (for ωr1 = ωr2, i.e., symmetric or
     * when primary and secondary resonant frequencies match):
     *
     *   nVout/Vin = 1 / |Denominator|
     *
     * where the denominator in terms of s = jω is given by Equation 2.
     *
     * This implementation uses the impedance-based form:
     *   H(jω) = Zm·Ro / (Z1·Zm + Z1·Z2 + Z1·Ro + Zm·Z2 + Zm·Ro)
     *
     * where:
     *   Z1 = jωL1 + 1/(jωC1)        (primary resonant impedance)
     *   Z2 = n²·(jωL2 + 1/(jωC2))   (secondary resonant impedance referred to primary)
     *   Zm = jωLm                     (magnetizing impedance)
     *   Ro = 8n²/(π²)·R_load         (FHA equivalent AC load resistance)
     */
    double CllcConverter::get_voltage_gain(double switchingFrequency,
                                            const CllcResonantParameters& params) {
        double omega_s = 2.0 * M_PI * switchingFrequency;
        std::complex<double> j(0.0, 1.0);
        std::complex<double> s = j * omega_s;

        double n = params.turnsRatio;
        double L1 = params.primaryResonantInductance;
        double C1 = params.primaryResonantCapacitance;
        double L2 = params.secondaryResonantInductance;
        double C2 = params.secondaryResonantCapacitance;
        double Lm = params.magnetizingInductance;
        double Ro = params.equivalentAcResistance;

        // Primary resonant impedance: Z1 = sL1 + 1/(sC1)
        std::complex<double> Z1 = s * L1 + 1.0 / (s * C1);

        // Secondary resonant impedance referred to primary: Z2 = n²·(sL2 + 1/(sC2))
        std::complex<double> Z2 = n * n * (s * L2 + 1.0 / (s * C2));

        // Magnetizing impedance: Zm = sLm
        std::complex<double> Zm = s * Lm;

        // Transfer function: H = Zm·Ro / (Z1·Zm + Z1·Z2 + Z1·Ro + Zm·Z2 + Zm·Ro)
        std::complex<double> numerator = Zm * Ro;
        std::complex<double> denominator = Z1 * Zm + Z1 * Z2 + Z1 * Ro + Zm * Z2 + Zm * Ro;

        return std::abs(numerator / denominator);
    }

    // =========================================================================
    // Design Requirements
    // =========================================================================

    DesignRequirements CllcConverter::process_design_requirements() {
        auto params = calculate_resonant_parameters();

        DesignRequirements designRequirements;

        // Turns ratio
        designRequirements.get_mutable_turns_ratios().clear();
        DimensionWithTolerance turnsRatioWithTolerance;
        turnsRatioWithTolerance.set_nominal(roundFloat(params.turnsRatio, 2));
        designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);

        // Magnetizing inductance
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(params.magnetizingInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        // Isolation sides: primary and one secondary
        designRequirements.set_isolation_sides(
            Topology::create_isolation_sides(1, false));

        designRequirements.set_topology(Topologies::CLLC_RESONANT_CONVERTER);

        return designRequirements;
    }

    // =========================================================================
    // Analytical Waveform Generation
    // =========================================================================

    /**
     * Generate analytical operating point for a single input voltage and operating condition.
     *
     * The CLLC converter waveforms are modeled using FHA:
     *
     * PRIMARY WINDING:
     *   Voltage: Bipolar rectangular ±Vin with dead time
     *            +Vin for [0, T/2 - td], 0 for dead time, -Vin for [T/2, T - td], 0
     *   Current: Sum of resonant current and magnetizing current
     *            - Magnetizing current: triangular, peak = Vin/(4·Lm·fs)
     *            - Resonant current: sinusoidal at frequency ≈ fr
     *            - Total Ip(t) ≈ I_res·sin(2π·fs·t) + Im(t)
     *
     * SECONDARY WINDING:
     *   Voltage: Bipolar rectangular ±Vout*n (reflected to primary, then /n to secondary)
     *   Current: Resonant current minus magnetizing current, scaled by n
     *            Isec(t) = n · (Ip(t) - Im(t))
     *
     * At resonance (fs = fr):
     *   - Primary current is quasi-sinusoidal
     *   - Secondary current completes exactly one half-sine per half period
     *   - Peak current: Ip_peak ≈ Pout/(Vin·η) · π/2 (fundamental approximation)
     */
    OperatingPoint CllcConverter::process_operating_points_for_input_voltage(
        double inputVoltage,
        const CllcOperatingPoint& cllcOpPoint,
        double turnsRatio,
        double magnetizingInductance,
        const CllcResonantParameters& params) {

        OperatingPoint operatingPoint;
        double switchingFrequency = cllcOpPoint.get_switching_frequency();
        double outputVoltage = cllcOpPoint.get_output_voltages()[0];
        double outputCurrent = cllcOpPoint.get_output_currents()[0];
        double outputPower = outputVoltage * outputCurrent;
        double n = turnsRatio;
        double Lm = magnetizingInductance;

        double period = 1.0 / switchingFrequency;
        double halfPeriod = period / 2.0;
        double td = deadTime;

        // Clamp dead time to a reasonable fraction of the half period
        if (td > halfPeriod * 0.1) {
            td = halfPeriod * 0.1;
        }

        double efficiency = get_efficiency().value_or(0.95);

        // Peak magnetizing current (triangular, linear ramp over half period minus dead time)
        // Im ramps from -Im_peak to +Im_peak over (T/2 - td), with Vin across Lm
        // Vin = Lm · dI/dt → Im_peak = Vin · (T/2 - td) / (2·Lm)
        double Im_peak = inputVoltage * (halfPeriod - td) / (2.0 * Lm);

        // RMS input current from power balance: Pin = Pout/η
        // For sinusoidal current: Irms = Ip_peak / √2
        // Pin = Vin · Irms · (2√2/π)  [for fundamental of square wave voltage]
        // So Ip_peak ≈ (π/2) · Pout / (Vin · η)  [approximate for FHA]
        double inputPower = outputPower / efficiency;
        double Ip_rms = inputPower / inputVoltage;  // approximate DC equivalent
        double Ip_resonant_peak = Ip_rms * M_PI / 2.0;  // FHA fundamental peak

        // Ensure resonant peak is larger than magnetizing peak (power transfer happens)
        if (Ip_resonant_peak < Im_peak * 1.2) {
            Ip_resonant_peak = Im_peak * 1.2;
        }

        // Secondary resonant current peak (referred to secondary side)
        // Isec = n · (Ip - Im), so secondary peak ≈ n · (Ip_res_peak - Im_peak)
        // But at the moment of peak resonant current, Im ≈ 0 (at midpoint)
        // So secondary peak ≈ n · Ip_resonant_peak (approximately)
        // More precisely, Isec_peak = n · Ip_resonant_peak when Im is small
        double Isec_peak = n * Ip_resonant_peak;

        // Number of sample points for waveform construction
        const int numPoints = 200;
        double dt = period / numPoints;

        // =====================================================================
        // PRIMARY WINDING
        // =====================================================================
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            std::vector<double> currentData(numPoints + 1);
            std::vector<double> voltageData(numPoints + 1);
            std::vector<double> timeData(numPoints + 1);

            for (int i = 0; i <= numPoints; ++i) {
                double t = i * dt;
                timeData[i] = t;

                // --- Primary voltage: bipolar rectangular with dead time ---
                // [0, T/2-td]: +Vin
                // [T/2-td, T/2]: 0 (dead time)
                // [T/2, T-td]: -Vin
                // [T-td, T]: 0 (dead time)
                double tMod = fmod(t, period);
                if (tMod < halfPeriod - td) {
                    voltageData[i] = inputVoltage;
                } else if (tMod < halfPeriod) {
                    voltageData[i] = 0.0;
                } else if (tMod < period - td) {
                    voltageData[i] = -inputVoltage;
                } else {
                    voltageData[i] = 0.0;
                }

                // --- Primary current: sinusoidal resonant + triangular magnetizing ---
                // Resonant component (approximately sinusoidal at switching frequency)
                double I_resonant = Ip_resonant_peak * sin(2.0 * M_PI * switchingFrequency * t);

                // Magnetizing component (triangular)
                // Ramps linearly from -Im_peak to +Im_peak over first half period,
                // then from +Im_peak to -Im_peak over second half period
                double I_mag;
                if (tMod < halfPeriod) {
                    I_mag = -Im_peak + 2.0 * Im_peak * tMod / halfPeriod;
                } else {
                    I_mag = Im_peak - 2.0 * Im_peak * (tMod - halfPeriod) / halfPeriod;
                }

                currentData[i] = I_resonant + I_mag;
            }

            currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
            currentWaveform.set_data(currentData);
            currentWaveform.set_time(timeData);

            voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
            voltageWaveform.set_data(voltageData);
            voltageWaveform.set_time(timeData);

            auto excitation = complete_excitation(currentWaveform, voltageWaveform,
                                                   switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // =====================================================================
        // SECONDARY WINDING
        // =====================================================================
        {
            Waveform currentWaveform;
            Waveform voltageWaveform;

            std::vector<double> currentData(numPoints + 1);
            std::vector<double> voltageData(numPoints + 1);
            std::vector<double> timeData(numPoints + 1);

            // Secondary voltage referred to secondary side: ±Vout (output voltage)
            // which is ±Vin/n from the primary perspective
            double secVoltage = outputVoltage;

            for (int i = 0; i <= numPoints; ++i) {
                double t = i * dt;
                timeData[i] = t;

                // --- Secondary voltage: bipolar rectangular ±Vout with dead time ---
                double tMod = fmod(t, period);
                if (tMod < halfPeriod - td) {
                    voltageData[i] = secVoltage;
                } else if (tMod < halfPeriod) {
                    voltageData[i] = 0.0;
                } else if (tMod < period - td) {
                    voltageData[i] = -secVoltage;
                } else {
                    voltageData[i] = 0.0;
                }

                // --- Secondary current ---
                // Isec(t) = n · (Ip(t) - Im(t)) = n · Iresonant(t)
                // The resonant current on the secondary is the resonant component only
                // (magnetizing current stays in the primary, doesn't transfer)
                double I_resonant = Ip_resonant_peak * sin(2.0 * M_PI * switchingFrequency * t);
                currentData[i] = n * I_resonant;
            }

            currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
            currentWaveform.set_data(currentData);
            currentWaveform.set_time(timeData);

            voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
            voltageWaveform.set_data(voltageData);
            voltageWaveform.set_time(timeData);

            auto excitation = complete_excitation(currentWaveform, voltageWaveform,
                                                   switchingFrequency, "Secondary 0");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(cllcOpPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    // =========================================================================
    // Process Operating Points (all input voltages × all operating points)
    // =========================================================================

    std::vector<OperatingPoint> CllcConverter::process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {

        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        double n = turnsRatios[0];

        // Calculate resonant parameters for ngspice and waveform generation
        CllcResonantParameters params = calculate_resonant_parameters();
        // Override with actual values provided
        params.turnsRatio = n;
        params.magnetizingInductance = magnetizingInductance;
        // Recalculate L1 from Lm and k
        params.primaryResonantInductance = magnetizingInductance / params.inductanceRatio;
        double L1 = params.primaryResonantInductance;
        double omega_r = 2.0 * M_PI * params.resonantFrequency;
        params.primaryResonantCapacitance = 1.0 / (omega_r * omega_r * L1);
        bool symmetric = get_symmetric_design().value_or(true);
        double a = symmetric ? 1.0 : 0.95;
        double b = symmetric ? 1.0 : 1.052;
        params.secondaryResonantInductance = a * L1 / (n * n);
        params.secondaryResonantCapacitance = n * n * b * params.primaryResonantCapacitance;

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto cllcOpPoint = get_operating_points()[opIndex];

                auto operatingPoint = process_operating_points_for_input_voltage(
                    inputVoltage, cllcOpPoint, n, magnetizingInductance, params);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(opIndex);
                }
                std::string flowStr = (cllcOpPoint.get_power_flow() == CllcPowerFlow::FORWARD) ?
                    " (Forward)" : " (Reverse)";
                name += flowStr;
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> CllcConverter::process_operating_points(Magnetic magnetic) {
        CllcConverter::run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel
            .calculate_inductance_from_number_turns_and_gapping(
                magnetic.get_mutable_core(), magnetic.get_mutable_coil())
            .get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();

        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    // =========================================================================
    // Ngspice Circuit Generation
    // =========================================================================

    /**
     * Generate a complete ngspice netlist for the CLLC converter.
     *
     * Circuit topology (forward mode):
     *
     *  Vin ─┬─ S1 ─┬─ S3 ─┬─ GND          GND ─┬─ Ds1 ─┬─ Ds3 ─┬─ Vout+
     *       │      A      │                     │       C       │
     *       │      │      │                     │       │       │
     *       └─ S2 ─┘─ S4 ─┘                    └─ Ds2 ─┘─ Ds4 ─┘
     *              B                                    D
     *
     *  A ── C1 ── L1 ── Lpri ──── B
     *                    ║ Lm
     *  C ── C2 ── L2 ── Lsec ──── D
     *
     * S1,S4 are controlled together (diagonal pair 1)
     * S2,S3 are controlled together (diagonal pair 2)
     * 50% duty cycle with dead time between transitions
     */
    std::string CllcConverter::generate_ngspice_circuit(
        double turnsRatio,
        const CllcResonantParameters& params,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {

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
        double switchingFrequency = opPoint.get_switching_frequency();
        double outputVoltage = opPoint.get_output_voltages()[0];
        double outputCurrent = opPoint.get_output_currents()[0];
        double n = turnsRatio;

        double L1 = params.primaryResonantInductance;
        double C1 = params.primaryResonantCapacitance;
        double L2 = params.secondaryResonantInductance;
        double C2 = params.secondaryResonantCapacitance;
        double Lm = params.magnetizingInductance;

        double period = 1.0 / switchingFrequency;
        double halfPeriod = period / 2.0;
        double td = deadTime;

        // Ensure dead time doesn't exceed reasonable portion of half period
        if (td > halfPeriod * 0.05) {
            td = halfPeriod * 0.05;
        }

        double tOn = halfPeriod - td;

        // Simulation timing
        int periodsToExtract = get_num_periods_to_extract();
        int steadyPeriods = get_num_steady_state_periods();
        int totalPeriods = steadyPeriods + periodsToExtract;
        double simTime = totalPeriods * period;
        double startTime = steadyPeriods * period;
        double stepTime = period / 200.0;

        // Load resistance
        double Rload = outputVoltage / outputCurrent;

        // Secondary inductance for coupled inductor model
        double Lsec = Lm / (n * n);  // Lm referred to secondary for coupling

        std::ostringstream circuit;

        circuit << "* CLLC Bidirectional Resonant Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vout=" << outputVoltage
                << "V, f=" << (switchingFrequency/1e3) << "kHz\n";
        circuit << "* n=" << n << ", L1=" << (L1*1e6) << "uH, C1=" << (C1*1e9) << "nF\n";
        circuit << "* Lm=" << (Lm*1e6) << "uH, L2=" << (L2*1e6) << "uH, C2=" << (C2*1e9) << "nF\n\n";

        // DC Input
        circuit << "* DC Input\n";
        circuit << "Vin vin_p 0 " << inputVoltage << "\n\n";

        // Switch models
        circuit << "* Switch and diode models\n";
        circuit << ".model SW1 SW VT=2.5 VH=0.5\n";
        circuit << ".model DIDEAL D(IS=1e-14 RS=0.01 CJO=1e-12)\n\n";

        // PWM control signals for full bridge
        // Pair 1: S1 (high-side leg A), S4 (low-side leg B) - first half cycle
        // Pair 2: S2 (low-side leg A), S3 (high-side leg B) - second half cycle
        circuit << "* PWM control signals (complementary pairs with dead time)\n";
        circuit << "Vpwm1 pwm1 0 PULSE(0 5 0 10n 10n " << std::scientific << tOn
                << " " << period << std::fixed << ")\n";
        circuit << "Vpwm2 pwm2 0 PULSE(0 5 " << std::scientific << halfPeriod
                << " 10n 10n " << tOn << " " << period << std::fixed << ")\n\n";

        // Primary full bridge
        circuit << "* Primary Full Bridge\n";
        circuit << "* Leg A: S1 (high-side), S2 (low-side)\n";
        circuit << "S1 vin_p node_a pwm1 0 SW1\n";
        circuit << "S2 node_a 0 pwm2 0 SW1\n";
        circuit << "* Leg B: S3 (high-side), S4 (low-side)\n";
        circuit << "S3 vin_p node_b pwm2 0 SW1\n";
        circuit << "S4 node_b 0 pwm1 0 SW1\n\n";

        // Primary current sense
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense node_a pri_c1_in 0\n\n";

        // Primary resonant tank: C1 in series with L1
        circuit << "* Primary Resonant Tank (C1 series with L1)\n";
        circuit << "C_res1 pri_c1_in pri_l1_in " << std::scientific << C1 << std::fixed << "\n";
        circuit << "L_res1 pri_l1_in pri_trafo_in " << std::scientific << L1 << std::fixed << "\n\n";

        // Transformer: coupled inductors (Lpri = Lm, Lsec = Lm/n²)
        // The magnetizing inductance IS the primary inductor of the transformer
        circuit << "* Transformer (coupled inductors)\n";
        circuit << "Lpri pri_trafo_in node_b " << std::scientific << Lm << std::fixed << "\n";
        circuit << "Lsec sec_trafo_p sec_trafo_n " << std::scientific << Lsec << std::fixed << "\n";
        circuit << "Kpri_sec Lpri Lsec 0.9999\n\n";

        // Secondary resonant tank: L2 in series with C2
        circuit << "* Secondary Resonant Tank (L2 series with C2)\n";
        circuit << "Vsec_sense sec_trafo_p sec_l2_in 0\n";
        circuit << "L_res2 sec_l2_in sec_c2_in " << std::scientific << L2 << std::fixed << "\n";
        circuit << "C_res2 sec_c2_in node_c " << std::scientific << C2 << std::fixed << "\n\n";

        // Secondary node_d = sec_trafo_n
        circuit << "* Secondary bridge reference\n";
        circuit << "Vd_ref sec_trafo_n node_d 0\n\n";

        // Secondary full bridge rectifier (diodes for simplicity)
        // In forward mode: Ds1,Ds4 conduct first half; Ds2,Ds3 conduct second half
        circuit << "* Secondary Full Bridge Rectifier (diodes)\n";
        circuit << "Ds1 node_c vout_p DIDEAL\n";
        circuit << "Ds2 vout_n node_c DIDEAL\n";
        circuit << "Ds3 node_d vout_p DIDEAL\n";
        circuit << "Ds4 vout_n node_d DIDEAL\n\n";

        // Snubber resistors for convergence
        circuit << "* Snubber resistors for convergence\n";
        circuit << "Rsnub1 node_c vout_p 1MEG\n";
        circuit << "Rsnub2 vout_n node_c 1MEG\n";
        circuit << "Rsnub3 node_d vout_p 1MEG\n";
        circuit << "Rsnub4 vout_n node_d 1MEG\n\n";

        // Ground reference for secondary
        circuit << "* Secondary ground reference\n";
        circuit << "Vgnd_sec vout_n 0 0\n\n";

        // Output filter and load
        circuit << "* Output filter and load\n";
        circuit << "Cout vout_p vout_n 100u IC=" << outputVoltage << "\n";
        circuit << "Rload vout_p vout_n " << Rload << "\n\n";

        // Transient analysis
        circuit << "* Transient Analysis\n";
        circuit << ".tran " << std::scientific << stepTime << " " << simTime
                << " " << startTime << std::fixed << " UIC\n\n";

        // Save signals
        circuit << "* Output signals\n";
        circuit << ".save v(pri_trafo_in) v(node_b) i(Vpri_sense)"
                << " v(sec_trafo_p) v(sec_trafo_n) i(Vsec_sense)"
                << " v(vout_p) v(vout_n)\n\n";

        // Options
        circuit << ".options RELTOL=0.003 ABSTOL=1e-8 VNTOL=1e-5 TRTOL=10 ITL1=500 ITL4=100\n";
        circuit << ".ic v(vout_p)=" << outputVoltage << "\n\n";

        circuit << ".end\n";

        return circuit.str();
    }

    // =========================================================================
    // Simulation-based Operating Point Extraction
    // =========================================================================

    std::vector<OperatingPoint> CllcConverter::simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {

        std::vector<OperatingPoint> operatingPoints;

        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        double n = turnsRatios[0];
        CllcResonantParameters params = calculate_resonant_parameters();
        params.turnsRatio = n;
        params.magnetizingInductance = magnetizingInductance;
        params.primaryResonantInductance = magnetizingInductance / params.inductanceRatio;
        double L1 = params.primaryResonantInductance;
        double omega_r = 2.0 * M_PI * params.resonantFrequency;
        params.primaryResonantCapacitance = 1.0 / (omega_r * omega_r * L1);
        bool symmetric = get_symmetric_design().value_or(true);
        double a = symmetric ? 1.0 : 0.95;
        double b = symmetric ? 1.0 : 1.052;
        params.secondaryResonantInductance = a * L1 / (n * n);
        params.secondaryResonantCapacitance = n * n * b * params.primaryResonantCapacitance;

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto cllcOpPoint = get_operating_points()[opIndex];

                std::string netlist = generate_ngspice_circuit(n, params, inputVoltageIndex, opIndex);
                double switchingFrequency = cllcOpPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);

                if (!simResult.success) {
                    throw std::runtime_error("CLLC Simulation failed: " + simResult.errorMessage);
                }

                // Map waveform names to winding excitations
                NgspiceRunner::WaveformNameMapping waveformMapping;

                // Primary winding: voltage across transformer primary = v(pri_trafo_in) - v(node_b)
                // Current through primary = i(Vpri_sense)
                waveformMapping.push_back({
                    {"voltage", "pri_trafo_in"},
                    {"current", "vpri_sense#branch"}
                });

                // Secondary winding: voltage = v(sec_trafo_p) - v(sec_trafo_n)
                // Current = i(Vsec_sense)
                waveformMapping.push_back({
                    {"voltage", "sec_trafo_p"},
                    {"current", "vsec_sense#branch"}
                });

                std::vector<std::string> windingNames = {"Primary", "Secondary"};
                std::vector<bool> flipCurrentSign = {false, false};

                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult, waveformMapping, switchingFrequency,
                    windingNames, cllcOpPoint.get_ambient_temperature(),
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

    // =========================================================================
    // Simulation-based Topology Waveform Extraction
    // =========================================================================

    std::vector<ConverterWaveforms> CllcConverter::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {

        std::vector<ConverterWaveforms> results;

        NgspiceRunner runner;
        if (!runner.is_available()) {
            throw std::runtime_error("ngspice is not available for simulation");
        }

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        double n = turnsRatios[0];
        CllcResonantParameters params = calculate_resonant_parameters();
        params.turnsRatio = n;
        params.magnetizingInductance = magnetizingInductance;
        params.primaryResonantInductance = magnetizingInductance / params.inductanceRatio;
        double L1 = params.primaryResonantInductance;
        double omega_r = 2.0 * M_PI * params.resonantFrequency;
        params.primaryResonantCapacitance = 1.0 / (omega_r * omega_r * L1);
        bool symmetric = get_symmetric_design().value_or(true);
        double a = symmetric ? 1.0 : 0.95;
        double b = symmetric ? 1.0 : 1.052;
        params.secondaryResonantInductance = a * L1 / (n * n);
        params.secondaryResonantCapacitance = n * n * b * params.primaryResonantCapacitance;

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto opPoint = get_operating_points()[opIndex];

                std::string netlist = generate_ngspice_circuit(n, params, inputVoltageIndex, opIndex);
                double switchingFrequency = opPoint.get_switching_frequency();

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);
                if (!simResult.success) {
                    throw std::runtime_error("CLLC Simulation failed: " + simResult.errorMessage);
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

                wf.set_input_voltage(getWaveform("pri_trafo_in"));
                wf.set_input_current(getWaveform("vpri_sense#branch"));
                wf.get_mutable_output_voltages().push_back(getWaveform("vout_p"));
                wf.get_mutable_output_currents().push_back(getWaveform("vsec_sense#branch"));

                results.push_back(wf);
            }
        }

        return results;
    }

    // =========================================================================
    // Advanced CLLC Converter
    // =========================================================================

    Inputs AdvancedCllcConverter::process() {
        CllcConverter::run_checks(_assertErrors);

        Inputs inputs;

        double magnetizingInductance = get_desired_magnetizing_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();

        if (turnsRatios.empty()) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS,
                "CLLC: desiredTurnsRatios must not be empty");
        }

        // Build design requirements
        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();

        for (auto tr : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(tr, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(magnetizingInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        designRequirements.set_isolation_sides(
            Topology::create_isolation_sides(1, false));
        designRequirements.set_topology(Topologies::CLLC_RESONANT_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        // Generate operating points
        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        CllcResonantParameters params = calculate_resonant_parameters();
        params.turnsRatio = turnsRatios[0];
        params.magnetizingInductance = magnetizingInductance;

        // Override resonant components if user specified them
        if (get_desired_resonant_inductance_primary()) {
            params.primaryResonantInductance = get_desired_resonant_inductance_primary().value();
        } else {
            params.primaryResonantInductance = magnetizingInductance / params.inductanceRatio;
        }

        double L1 = params.primaryResonantInductance;

        if (get_desired_resonant_capacitance_primary()) {
            params.primaryResonantCapacitance = get_desired_resonant_capacitance_primary().value();
        } else {
            double omega_r = 2.0 * M_PI * params.resonantFrequency;
            params.primaryResonantCapacitance = 1.0 / (omega_r * omega_r * L1);
        }

        double n = turnsRatios[0];
        bool symmetric = get_symmetric_design().value_or(true);
        double a = symmetric ? 1.0 : 0.95;
        double b = symmetric ? 1.0 : 1.052;

        if (get_desired_resonant_inductance_secondary()) {
            params.secondaryResonantInductance = get_desired_resonant_inductance_secondary().value();
        } else {
            params.secondaryResonantInductance = a * L1 / (n * n);
        }

        if (get_desired_resonant_capacitance_secondary()) {
            params.secondaryResonantCapacitance = get_desired_resonant_capacitance_secondary().value();
        } else {
            params.secondaryResonantCapacitance = n * n * b * params.primaryResonantCapacitance;
        }

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto cllcOpPoint = get_operating_points()[opIndex];

                auto operatingPoint = process_operating_points_for_input_voltage(
                    inputVoltage, cllcOpPoint, n, magnetizingInductance, params);

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

} // namespace OpenMagnetics
