#include "converter_models/Dab.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "support/Exceptions.h"

namespace OpenMagnetics {

// =========================================================================
// Construction
// =========================================================================
Dab::Dab(const json& j) {
    from_json(j, *static_cast<DualActiveBridge*>(this));
}

AdvancedDab::AdvancedDab(const json& j) {
    from_json(j, *this);
}

// =========================================================================
// Static helper: Power transfer (SPS modulation)
// =========================================================================
// P = N * V1 * V2 * phi * (pi - |phi|) / (2 * pi^2 * Fs * L)
// Reference: [1] TI TIDA-010054 Eq.6, [2] Demetriades Ch.6
// =========================================================================
double Dab::compute_power(double V1, double V2, double N,
                          double phi, double Fs, double L) {
    return N * V1 * V2 * phi * (M_PI - std::abs(phi))
           / (2.0 * M_PI * M_PI * Fs * L);
}

// =========================================================================
// Static helper: Series inductance for desired power at given phase shift
// =========================================================================
// L = N * V1 * V2 * phi * (pi - phi) / (2 * pi^2 * Fs * P)
// =========================================================================
double Dab::compute_series_inductance(double V1, double V2, double N,
                                      double phi, double Fs, double P) {
    if (P <= 0) return 1e-3; // fallback
    return N * V1 * V2 * phi * (M_PI - std::abs(phi))
           / (2.0 * M_PI * M_PI * Fs * P);
}

// =========================================================================
// Static helper: Phase shift for desired power with given inductance
// =========================================================================
// phi = (pi/2) * (1 - sqrt(1 - 8*Fs*L*P / (N*V1*V2)))
// Reference: [1] TI TIDA-010054 Eq.16
// =========================================================================
double Dab::compute_phase_shift(double V1, double V2, double N,
                                double Fs, double L, double P) {
    double discriminant = 1.0 - 8.0 * Fs * L * P / (N * V1 * V2);
    if (discriminant < 0) {
        // Power exceeds maximum transferable power
        return M_PI / 2.0; // Maximum phase shift
    }
    return (M_PI / 2.0) * (1.0 - std::sqrt(discriminant));
}

// =========================================================================
// Static helper: Voltage conversion ratio
// =========================================================================
double Dab::compute_voltage_ratio(double V1, double V2, double N) {
    return N * V2 / V1;
}

// =========================================================================
// Static helper: Inductor current at switching instants
// =========================================================================
// Reference: [1] TI TIDA-010054 Eq.7-10
//   d = N * V2 / V1
//   Ibase = V1 / (2*pi*Fs*L)
//   i1 = 0.5 * (2*phi - (1-d)*pi) * Ibase
//   i2 = 0.5 * (2*d*phi + (1-d)*pi) * Ibase
//
// i1 = current at t=0 (start of positive half-cycle, before phase shift)
// i2 = current at t=phi/(2*pi*Fs) (at the phase shift instant)
// =========================================================================
void Dab::compute_switching_currents(double V1, double V2, double N,
                                     double phi, double Fs, double L,
                                     double& i1, double& i2) {
    double d = N * V2 / V1;
    double Ibase = V1 / (2.0 * M_PI * Fs * L);
    i1 = 0.5 * (2.0 * phi - (1.0 - d) * M_PI) * Ibase;
    i2 = 0.5 * (2.0 * d * phi + (1.0 - d) * M_PI) * Ibase;
}

// =========================================================================
// Static helper: Primary RMS current
// =========================================================================
// Reference: [1] TI TIDA-010054 Eq.14
//   Ip_rms = sqrt(1/3 * (i1^2 + i2^2 + (1 - 2*phi/pi)*i1*i2))
// =========================================================================
double Dab::compute_primary_rms_current(double i1, double i2, double phi) {
    double factor = 1.0 - 2.0 * phi / M_PI;
    return std::sqrt((i1 * i1 + i2 * i2 + factor * i1 * i2) / 3.0);
}

// =========================================================================
// Static helper: ZVS check
// =========================================================================
// Reference: [1] TI TIDA-010054 Eq.11-12
//   Primary ZVS:   phi > (1 - 1/d) * pi/2
//   Secondary ZVS: phi > (1 - d) * pi/2
// =========================================================================
bool Dab::check_zvs_primary(double phi, double d) {
    if (d <= 0) return false;
    double phi_min = (1.0 - 1.0 / d) * M_PI / 2.0;
    return phi > phi_min;
}

bool Dab::check_zvs_secondary(double phi, double d) {
    double phi_min = (1.0 - d) * M_PI / 2.0;
    return phi > phi_min;
}


// =========================================================================
// Input validation
// =========================================================================
bool Dab::run_checks(bool assertErrors) {
    _assertErrors = assertErrors;
    bool ok = true;

    auto& ops = get_operating_points();
    if (ops.empty()) {
        if (assertErrors) throw std::runtime_error("DAB: no operating points");
        return false;
    }

    for (size_t i = 0; i < ops.size(); i++) {
        auto& op = ops[i];
        if (op.get_output_voltages().empty() || op.get_output_currents().empty()) {
            if (assertErrors) throw std::runtime_error("DAB: OP missing voltages/currents");
            ok = false;
        }
        if (op.get_output_voltages().size() != op.get_output_currents().size()) {
            if (assertErrors) throw std::runtime_error("DAB: voltage/current count mismatch");
            ok = false;
        }
        double fsw = op.get_switching_frequency();
        if (fsw <= 0) {
            if (assertErrors) throw std::runtime_error("DAB: invalid switching frequency");
            ok = false;
        }
        // Phase shift should be in range (-90, 90) degrees for SPS
        double phi_deg = op.get_phase_shift();
        if (std::abs(phi_deg) > 90.0) {
            if (assertErrors) throw std::runtime_error("DAB: phase shift out of range (|phi| > 90 deg)");
            ok = false;
        }
    }
    return ok;
}


// =========================================================================
// Design Requirements
// =========================================================================
//
// Design flow:
// 1. N = V1_nom / V2_nom (turns ratio)
// 2. Compute output power from output voltage and current
// 3. If series_inductance is given, use it; otherwise compute from power & phase shift
// 4. If phase shift is given, use it; otherwise compute from power & inductance
// 5. Compute magnetizing inductance Lm (large, to keep mag. current small)
// 6. Build DesignRequirements with turns ratio, Lm, and optionally leakage = L
// =========================================================================
DesignRequirements Dab::process_design_requirements() {
    auto& inputVoltage = get_input_voltage();
    double Vin_nom = inputVoltage.get_nominal().value_or(
        (inputVoltage.get_minimum().value_or(0) + inputVoltage.get_maximum().value_or(0)) / 2.0);

    auto& ops = get_operating_points();
    double mainOutputVoltage = ops[0].get_output_voltages()[0];
    double mainOutputCurrent = ops[0].get_output_currents()[0];
    double mainOutputPower = mainOutputVoltage * mainOutputCurrent;
    double Fs = ops[0].get_switching_frequency();

    // 1. Turns ratio: N = V1_nom / V2_nom
    double N = Vin_nom / mainOutputVoltage;

    // Turns ratios for all secondaries
    std::vector<double> turnsRatios;
    turnsRatios.push_back(N);
    for (size_t i = 1; i < ops[0].get_output_voltages().size(); i++) {
        turnsRatios.push_back(Vin_nom / ops[0].get_output_voltages()[i]);
    }

    // 2. Phase shift (from operating point, converted from degrees to radians)
    double phi_deg = ops[0].get_phase_shift();
    double phi_rad = phi_deg * M_PI / 180.0;

    // 3. Series inductance
    double L;
    if (get_series_inductance().has_value() && get_series_inductance().value() > 0) {
        L = get_series_inductance().value();
        // If phase shift is zero or very small, compute it from power
        if (std::abs(phi_rad) < 1e-6 && mainOutputPower > 0) {
            phi_rad = compute_phase_shift(Vin_nom, mainOutputVoltage, N, Fs, L, mainOutputPower);
        }
    } else {
        // Compute L from power and phase shift
        if (std::abs(phi_rad) < 1e-6) {
            // Default phase shift: target ~20-30 degrees for good controllability
            phi_rad = 25.0 * M_PI / 180.0;
        }
        L = compute_series_inductance(Vin_nom, mainOutputVoltage, N, phi_rad, Fs, mainOutputPower);
    }

    computedSeriesInductance = L;
    computedPhaseShift = phi_rad;

    // 4. Magnetizing inductance
    //    Lm should be large enough that magnetizing current is small (~5-10% of load current)
    //    Im_peak = V1 / (4 * Fs * Lm)
    //    Target: Im_peak < 0.1 * I_load_primary
    //
    //    I_load_primary (approximate) = P / V1
    //    So Lm > V1 / (4 * Fs * 0.1 * P / V1) = V1^2 / (0.4 * Fs * P)
    //
    //    We use a factor of 20x the series inductance as a minimum,
    //    or compute from the 10% magnetizing current criterion
    double I_load_pri = mainOutputPower / Vin_nom;
    double Im_target = 0.1 * I_load_pri; // 10% of load current
    double Lm_from_current = (Im_target > 0) ? Vin_nom / (4.0 * Fs * Im_target) : 20.0 * L;
    double Lm_from_ratio = 20.0 * L; // At least 20x series inductance
    double Lm = std::max(Lm_from_current, Lm_from_ratio);

    computedMagnetizingInductance = Lm;

    // 5. Build DesignRequirements
    DesignRequirements designRequirements;
    designRequirements.get_mutable_turns_ratios().clear();
    for (auto n : turnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(roundFloat(n, 2));
        designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }

    DimensionWithTolerance inductanceWithTolerance;
    inductanceWithTolerance.set_nominal(roundFloat(Lm, 10));
    designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

    // If using leakage inductance as series inductor, request leakage = L
    if (get_use_leakage_inductance().value_or(false)) {
        std::vector<DimensionWithTolerance> leakageReqs;
        DimensionWithTolerance lrTol;
        lrTol.set_nominal(roundFloat(L, 10));
        leakageReqs.push_back(lrTol);
        designRequirements.set_leakage_inductance(leakageReqs);
    }

    return designRequirements;
}


// =========================================================================
// Process operating points for all input voltages
// =========================================================================
std::vector<OperatingPoint> Dab::process_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    std::vector<OperatingPoint> result;
    auto& inputVoltage = get_input_voltage();
    auto& ops = get_operating_points();

    // Collect input voltages: nominal, minimum, maximum
    std::vector<double> inputVoltages;
    if (inputVoltage.get_nominal().has_value())
        inputVoltages.push_back(inputVoltage.get_nominal().value());
    if (inputVoltage.get_minimum().has_value())
        inputVoltages.push_back(inputVoltage.get_minimum().value());
    if (inputVoltage.get_maximum().has_value())
        inputVoltages.push_back(inputVoltage.get_maximum().value());

    // Remove duplicates
    std::sort(inputVoltages.begin(), inputVoltages.end());
    inputVoltages.erase(std::unique(inputVoltages.begin(), inputVoltages.end()), inputVoltages.end());

    for (double Vin : inputVoltages) {
        auto op = process_operating_point_for_input_voltage(
            Vin, ops[0], turnsRatios, magnetizingInductance);
        result.push_back(op);
    }
    return result;
}

std::vector<OperatingPoint> Dab::process_operating_points(Magnetic magnetic) {
    auto req = process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    return process_operating_points(turnsRatios, Lm);
}


// =========================================================================
// CORE WAVEFORM GENERATION - Analytical piecewise linear model
// =========================================================================
//
// Reference: [1] TI TIDA-010054 Section 2.3.2, Eq. 2-5, 7-8
//            [2] Demetriades Chapter 6, Figure 6.3
//
// The inductor current in a DAB is piecewise linear with 4 segments per period.
// The primary and secondary bridge voltages are square waves phase-shifted by phi.
//
// Positive half-cycle (0 <= t <= Ts/2):
//
// Primary bridge voltage:  Vab = +V1  for 0 <= t <= Ts/2
// Secondary bridge voltage: Vcd = -V2 for 0 <= t < t_phi  (interval 1)
//                           Vcd = +V2 for t_phi <= t < Ts/2 (interval 2)
//
// Inductor voltage (referred to primary, Vab - N*Vcd):
//   Interval 1: vL = V1 + N*V2  --> current ramps UP steeply
//   Interval 2: vL = V1 - N*V2  --> current ramps (up if d<1, down if d>1)
//
// Negative half-cycle by antisymmetry: iL(t + Ts/2) = -iL(t)
//
// Transformer primary excitation for MAS:
//   Voltage across transformer primary = Vab (the H-bridge 1 output square wave)
//   This is the voltage that drives the flux in the core.
//   Current through primary = iL(t)
//
// Transformer secondary excitation for MAS:
//   Voltage across secondary = Vcd (the H-bridge 2 output square wave)
//   Current through secondary = N * iL(t) (reflected to secondary)
//
// Magnetizing current:
//   The transformer magnetizing inductance Lm sees the primary voltage Vab.
//   Im(t) is a triangular wave: slope = +V1/Lm for positive half, -V1/Lm for negative.
//   Im_peak = V1 / (4 * Fs * Lm)
//
// =========================================================================
OperatingPoint Dab::process_operating_point_for_input_voltage(
    double inputVoltage,
    const DabOperatingPoint& dabOpPoint,
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    OperatingPoint operatingPoint;

    double Fs = dabOpPoint.get_switching_frequency();
    double V1 = inputVoltage;
    double V2 = dabOpPoint.get_output_voltages()[0];
    double N = turnsRatios[0]; // Primary-to-secondary turns ratio
    double Lm = magnetizingInductance;
    double L = computedSeriesInductance;

    // Phase shift: from operating point (degrees) or computed
    double phi_deg = dabOpPoint.get_phase_shift();
    double phi_rad;
    if (std::abs(phi_deg) > 1e-6) {
        phi_rad = phi_deg * M_PI / 180.0;
    } else {
        phi_rad = computedPhaseShift;
    }

    double period = 1.0 / Fs;
    double Thalf = period / 2.0;
    double t_phi = std::abs(phi_rad) / (2.0 * M_PI * Fs); // Time for phase shift

    // Sign of phi determines power flow direction
    double phi_sign = (phi_rad >= 0) ? 1.0 : -1.0;
    double phi_abs = std::abs(phi_rad);

    // Compute switching currents
    double i1, i2;
    compute_switching_currents(V1, V2, N, phi_abs, Fs, L, i1, i2);
    // For negative phase shift (reverse power flow), negate currents
    if (phi_sign < 0) {
        i1 = -i1;
        i2 = -i2;
    }

    // Sampling
    const int N_samples = 256; // Samples per half-period
    double dt = Thalf / N_samples;

    // Build full-period waveforms
    int totalSamples = 2 * N_samples + 1;
    std::vector<double> time_full(totalSamples);
    std::vector<double> iL_full(totalSamples);      // Inductor current
    std::vector<double> Vab_full(totalSamples);      // Primary bridge voltage
    std::vector<double> Vcd_full(totalSamples);      // Secondary bridge voltage
    std::vector<double> Im_full(totalSamples);       // Magnetizing current

    // Voltage slopes
    double slope1 = (V1 + N * V2) / L;  // Interval 1: Vab=+V1, Vcd=-V2
    double slope2 = (V1 - N * V2) / L;  // Interval 2: Vab=+V1, Vcd=+V2
    double Im_slope = V1 / Lm;          // Magnetizing current slope
    double Im_peak = V1 / (4.0 * Fs * Lm);

    // Positive half-cycle (0 <= t <= Thalf)
    for (int k = 0; k <= N_samples; ++k) {
        double t = k * dt;
        time_full[k] = t;

        // Primary bridge: Vab = +V1
        Vab_full[k] = V1;

        // Secondary bridge: Vcd depends on phase shift
        if (t < t_phi) {
            // Interval 1: secondary hasn't switched yet
            // For positive phi (forward): Vcd = -V2
            Vcd_full[k] = -phi_sign * V2;
            iL_full[k] = i1 + phi_sign * slope1 * t;
        } else {
            // Interval 2: secondary has switched
            // For positive phi (forward): Vcd = +V2
            Vcd_full[k] = phi_sign * V2;
            iL_full[k] = i2 + phi_sign * slope2 * (t - t_phi);
        }

        // Magnetizing current: triangular, starts at -Im_peak at t=0
        Im_full[k] = -Im_peak + Im_slope * t;
    }

    // Negative half-cycle by antisymmetry
    for (int k = 1; k <= N_samples; ++k) {
        time_full[N_samples + k] = Thalf + k * dt;
        iL_full[N_samples + k] = -iL_full[k];
        Vab_full[N_samples + k] = -Vab_full[k];
        Vcd_full[N_samples + k] = -Vcd_full[k];
        Im_full[N_samples + k] = -Im_full[k];
    }

    // ---- Primary winding excitation ----
    // Current: iL(t) (the inductor current flows through the transformer primary)
    // Voltage: Vab(t) (the full-bridge output voltage is across the transformer primary
    //          + series inductor. The voltage across just the transformer primary is
    //          Vab - L*diL/dt. However, for MAS operating point, we use the voltage
    //          that appears across the transformer magnetizing inductance, which drives
    //          the flux. In practice, if use_leakage_inductance is true, the series L
    //          is the leakage, and the voltage across Lm is Vab - L*diL/dt.
    //          If use_leakage_inductance is false, the transformer primary voltage
    //          is approximately Vab (the external inductor drops the difference).
    //          For the MAS design, we provide the bipolar rectangular waveform
    //          that the core sees, which is approximately Vab for the transformer
    //          primary voltage.
    {
        // For the transformer core flux calculation, the relevant voltage is
        // across the magnetizing inductance. We compute it properly:
        // V_Lm = Lm * dIm/dt = Lm * (V1/Lm) = V1 during positive half
        // This is simply the bipolar rectangular wave +/- V1
        // (same as Vab since the series L drop doesn't affect the core voltage
        //  in the ideal transformer model)
        std::vector<double> V_primary(totalSamples);
        for (int k = 0; k < totalSamples; ++k) {
            V_primary[k] = Vab_full[k];
        }

        Waveform currentWaveform;
        currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        currentWaveform.set_data(iL_full);
        currentWaveform.set_time(time_full);

        Waveform voltageWaveform;
        voltageWaveform.set_ancillary_label(WaveformLabel::BIPOLAR_RECTANGULAR);
        voltageWaveform.set_data(V_primary);
        voltageWaveform.set_time(time_full);

        auto excitation = complete_excitation(currentWaveform, voltageWaveform,
                                              Fs, "Primary");
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    // ---- Secondary winding excitation(s) ----
    for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
        double n = turnsRatios[secIdx];
        double Vout = dabOpPoint.get_output_voltages()[secIdx];

        std::vector<double> iSecData(totalSamples);
        std::vector<double> vSecData(totalSamples);

        for (int k = 0; k < totalSamples; ++k) {
            // Secondary current = N * iL (reflected from primary)
            iSecData[k] = n * iL_full[k];
            // Secondary voltage = Vcd (H-bridge 2 output)
            vSecData[k] = Vcd_full[k];
        }

        Waveform secCurrentWfm;
        secCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        secCurrentWfm.set_data(iSecData);
        secCurrentWfm.set_time(time_full);

        Waveform secVoltageWfm;
        secVoltageWfm.set_ancillary_label(WaveformLabel::BIPOLAR_RECTANGULAR);
        secVoltageWfm.set_data(vSecData);
        secVoltageWfm.set_time(time_full);

        auto excitation = complete_excitation(secCurrentWfm, secVoltageWfm,
                                              Fs,
                                              "Secondary " + std::to_string(secIdx));
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    // Operating conditions
    OperatingConditions conditions;
    conditions.set_ambient_temperature(dabOpPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    return operatingPoint;
}


// =========================================================================
// SPICE Circuit Generation
// =========================================================================
//
// Generates an NgSpice netlist for the DAB converter.
//
// Circuit topology:
//   V1 (DC) -> Full Bridge 1 (Q1-Q4) -> L_series -> Transformer (Np:Ns)
//            -> Full Bridge 2 (Q5-Q8, phase-shifted) -> V2 / R_load
//
// The transformer is modeled as a coupled inductor pair:
//   L_primary (magnetizing) coupled to L_secondary with k ≈ 1
//   The leakage is modeled separately as L_series
//
// PWM signals:
//   Primary bridge: 50% duty, period = 1/Fs
//   Secondary bridge: 50% duty, phase-shifted by phi from primary
//
// =========================================================================
std::string Dab::generate_ngspice_circuit(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t inputVoltageIndex,
    size_t operatingPointIndex)
{
    auto& inputVoltageSpec = get_input_voltage();
    auto& ops = get_operating_points();

    // Select input voltage
    std::vector<double> inputVoltages;
    if (inputVoltageSpec.get_nominal().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_nominal().value());
    if (inputVoltageSpec.get_minimum().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_minimum().value());
    if (inputVoltageSpec.get_maximum().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_maximum().value());

    double V1 = inputVoltages[std::min(inputVoltageIndex, inputVoltages.size() - 1)];
    auto& dabOp = ops[std::min(operatingPointIndex, ops.size() - 1)];

    double Fs = dabOp.get_switching_frequency();
    double period = 1.0 / Fs;
    double halfPeriod = period / 2.0;
    double deadTime = computedDeadTime;
    double tOn = halfPeriod - deadTime;

    double V2 = dabOp.get_output_voltages()[0];
    double Iout = dabOp.get_output_currents()[0];
    double N = turnsRatios[0];

    double L = computedSeriesInductance;
    double Lm = magnetizingInductance;

    // Phase shift time delay
    double phi_deg = dabOp.get_phase_shift();
    double phi_rad = (std::abs(phi_deg) > 1e-6) ? phi_deg * M_PI / 180.0 : computedPhaseShift;
    double phaseDelay = std::abs(phi_rad) / (2.0 * M_PI * Fs);

    // Simulation timing
    int periodsToExtract = get_num_periods_to_extract();
    int steadyStatePeriods = get_num_steady_state_periods();
    int numPeriodsTotal = steadyStatePeriods + periodsToExtract;
    double simTime = numPeriodsTotal * period;
    double startTime = steadyStatePeriods * period;
    double stepTime = period / 500;

    std::ostringstream circuit;

    circuit << "* Dual Active Bridge (DAB) Converter - Generated by OpenMagnetics\n";
    circuit << "* V1=" << V1 << "V, V2=" << V2 << "V, Fs=" << (Fs/1e3)
            << "kHz, phi=" << phi_deg << "deg\n";
    circuit << "* N=" << N << ", L=" << (L*1e6) << "uH, Lm=" << (Lm*1e6) << "uH\n\n";

    // Switch and diode models
    circuit << ".model SW1 SW(Ron=10m Roff=10Meg Vt=2.5)\n";
    circuit << ".model DIDEAL D(Is=1e-14 N=0.001)\n\n";

    // DC input voltage
    circuit << "Vdc1 vin_dc1 0 " << V1 << "\n\n";

    // ==== PRIMARY FULL BRIDGE (Q1-Q4) ====
    circuit << "* Primary Full Bridge\n";
    circuit << "Vpwm_p1 pwm_p1 0 PULSE(0 5 0 10n 10n "
            << std::scientific << tOn << " " << period << std::fixed << ")\n";
    circuit << "Vpwm_p2 pwm_p2 0 PULSE(0 5 "
            << std::scientific << halfPeriod << std::fixed
            << " 10n 10n "
            << std::scientific << tOn << " " << period << std::fixed << ")\n\n";

    // Q1, Q4: driven by pwm_p1 (positive half)
    circuit << "S1 vin_dc1 bridge_a1 pwm_p1 0 SW1\n";
    circuit << "D1 0 bridge_a1 DIDEAL\n";
    circuit << "S4 bridge_b1 0 pwm_p1 0 SW1\n";
    circuit << "D4 bridge_b1 vin_dc1 DIDEAL\n\n";

    // Q2, Q3: driven by pwm_p2 (negative half)
    circuit << "S2 bridge_a1 0 pwm_p2 0 SW1\n";
    circuit << "D2 bridge_a1 vin_dc1 DIDEAL\n";
    circuit << "S3 vin_dc1 bridge_b1 pwm_p2 0 SW1\n";
    circuit << "D3 0 bridge_b1 DIDEAL\n\n";

    // Sense primary current
    circuit << "Vpri_sense bridge_a1 pri_out 0\n\n";

    // ==== SERIES INDUCTANCE ====
    circuit << "* Series inductance (leakage + external)\n";
    circuit << "L_series pri_out trafo_pri " << std::scientific << L << "\n\n";

    // ==== TRANSFORMER ====
    // Model as coupled inductors with coupling coefficient k
    // Lp = Lm (magnetizing), Ls = Lm/N^2
    // k = 1 (ideal coupling, leakage modeled separately)
    circuit << "* Transformer Np:Ns = " << N << ":1\n";
    double Ls_sec = Lm / (N * N);
    circuit << "L_pri trafo_pri bridge_b1 " << std::scientific << Lm << "\n";
    circuit << "L_sec trafo_sec_a bridge_sec_b " << std::scientific << Ls_sec << "\n";
    circuit << "K_trafo L_pri L_sec 0.9999\n\n";

    // ==== SECONDARY FULL BRIDGE (Q5-Q8) ====
    circuit << "* Secondary Full Bridge (phase-shifted by phi)\n";
    double secDelay = phaseDelay;

    circuit << "Vpwm_s1 pwm_s1 0 PULSE(0 5 "
            << std::scientific << secDelay << std::fixed
            << " 10n 10n "
            << std::scientific << tOn << " " << period << std::fixed << ")\n";
    circuit << "Vpwm_s2 pwm_s2 0 PULSE(0 5 "
            << std::scientific << (halfPeriod + secDelay) << std::fixed
            << " 10n 10n "
            << std::scientific << tOn << " " << period << std::fixed << ")\n\n";

    // DC output (modeled as voltage source for validation, or as R_load)
    circuit << "Vdc2 vin_dc2 0 " << V2 << "\n\n";

    // Q5, Q8: driven by pwm_s1
    circuit << "S5 vin_dc2 trafo_sec_a pwm_s1 0 SW1\n";
    circuit << "D5 0 trafo_sec_a DIDEAL\n";
    circuit << "S8 bridge_sec_b 0 pwm_s1 0 SW1\n";
    circuit << "D8 bridge_sec_b vin_dc2 DIDEAL\n\n";

    // Q6, Q7: driven by pwm_s2
    circuit << "S6 trafo_sec_a 0 pwm_s2 0 SW1\n";
    circuit << "D6 trafo_sec_a vin_dc2 DIDEAL\n";
    circuit << "S7 vin_dc2 bridge_sec_b pwm_s2 0 SW1\n";
    circuit << "D7 0 bridge_sec_b DIDEAL\n\n";

    // Simulation commands
    circuit << ".tran " << std::scientific << stepTime
            << " " << simTime << " " << startTime << "\n";
    circuit << ".end\n";

    return circuit.str();
}


// =========================================================================
// SPICE simulation wrappers
// =========================================================================
std::vector<OperatingPoint> Dab::simulate_and_extract_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    // Placeholder: run NgSpice and extract operating points
    // For now, fall back to analytical
    return process_operating_points(turnsRatios, magnetizingInductance);
}

std::vector<ConverterWaveforms> Dab::simulate_and_extract_topology_waveforms(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t numberOfPeriods)
{
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
                throw std::runtime_error("DAB simulation failed: " + simResult.errorMessage);
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
            
            wf.set_input_voltage(getWaveform("v(pri_out)"));
            wf.set_input_current(getWaveform("i(vpri_sense)"));
            
            if (!turnsRatios.empty()) {
                wf.get_mutable_output_voltages().push_back(getWaveform("v(trafo_sec_a)"));
                wf.get_mutable_output_currents().push_back(getWaveform("i(vdc2)"));
            }
            
            results.push_back(wf);
        }
    }
    
    // Restore original value
    set_num_periods_to_extract(originalNumPeriodsToExtract);
    
    return results;
}


// =========================================================================
// AdvancedDab
// =========================================================================
Inputs AdvancedDab::process() {
    auto designRequirements = process_design_requirements();

    // Override turns ratios
    designRequirements.get_mutable_turns_ratios().clear();
    for (auto n : desiredTurnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(n);
        designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }

    // Override magnetizing inductance
    DimensionWithTolerance LmTol;
    LmTol.set_nominal(desiredMagnetizingInductance);
    designRequirements.set_magnetizing_inductance(LmTol);

    // Override series inductance if specified
    if (desiredSeriesInductance.has_value()) {
        set_computed_series_inductance(desiredSeriesInductance.value());
    }

    auto ops = process_operating_points(desiredTurnsRatios, desiredMagnetizingInductance);

    Inputs inputs;
    inputs.set_design_requirements(designRequirements);
    inputs.set_operating_points(ops);

    return inputs;
}

} // namespace OpenMagnetics
