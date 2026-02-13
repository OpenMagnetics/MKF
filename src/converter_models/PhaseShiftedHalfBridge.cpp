#include "converter_models/PhaseShiftedHalfBridge.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "support/Exceptions.h"

namespace OpenMagnetics {

Pshb::Pshb(const json& j) {
    from_json(j, *static_cast<PhaseShiftFullBridge*>(this));
}

AdvancedPshb::AdvancedPshb(const json& j) {
    from_json(j, *this);
}

// =========================================================================
// Static: effective duty cycle from phase shift
// D_eff = phaseShift(deg) / 180
// =========================================================================
double Pshb::compute_effective_duty_cycle(double phaseShiftDeg) {
    return std::abs(phaseShiftDeg) / 180.0;
}

// =========================================================================
// Static: output voltage  (NOTE: Vin here is the FULL input voltage;
//         the Vin/2 factor is applied internally)
// =========================================================================
double Pshb::compute_output_voltage(double Vin, double Deff, double n,
                                    double Vd, PsfbRectifierType rectType) {
    double Vhb = Vin * BRIDGE_VOLTAGE_FACTOR;  // Vin/2
    switch (rectType) {
        case PsfbRectifierType::CENTER_TAPPED:
            return Vhb * Deff / n - Vd;
        case PsfbRectifierType::CURRENT_DOUBLER:
            return Vhb * Deff / (2.0 * n) - Vd;
        case PsfbRectifierType::FULL_BRIDGE:
        default:
            return Vhb * Deff / n - 2.0 * Vd;
    }
}

// =========================================================================
// Static: turns ratio for target output voltage
// =========================================================================
double Pshb::compute_turns_ratio(double Vin, double Vo, double Deff,
                                 double Vd, PsfbRectifierType rectType) {
    double Vhb = Vin * BRIDGE_VOLTAGE_FACTOR;
    switch (rectType) {
        case PsfbRectifierType::CENTER_TAPPED:
            return Vhb * Deff / (Vo + Vd);
        case PsfbRectifierType::CURRENT_DOUBLER:
            return Vhb * Deff / (2.0 * (Vo + Vd));
        case PsfbRectifierType::FULL_BRIDGE:
        default:
            return Vhb * Deff / (Vo + 2.0 * Vd);
    }
}

// =========================================================================
// Static: output inductor
// Lo = Vo * (1 - Deff) / (Fs * dIo)
// =========================================================================
double Pshb::compute_output_inductance(double Vo, double Deff, double Fs,
                                       double Io, double rippleRatio) {
    double dIo = rippleRatio * Io;
    if (dIo <= 0) return 1e-3;
    return Vo * (1.0 - Deff) / (Fs * dIo);
}

// =========================================================================
// Static: primary RMS current (simplified)
// =========================================================================
double Pshb::compute_primary_rms_current(double Io, double n, double Deff) {
    return (Io / n) * std::sqrt(Deff);
}


// =========================================================================
// Input validation
// =========================================================================
bool Pshb::run_checks(bool assertErrors) {
    _assertErrors = assertErrors;
    bool ok = true;
    auto& ops = get_operating_points();
    if (ops.empty()) {
        if (assertErrors) throw std::runtime_error("PSHB: no operating points");
        return false;
    }
    for (size_t i = 0; i < ops.size(); i++) {
        auto& op = ops[i];
        if (op.get_output_voltages().empty() || op.get_output_currents().empty()) {
            if (assertErrors) throw std::runtime_error("PSHB: OP missing voltages/currents");
            ok = false;
        }
        double phi = op.get_phase_shift();
        if (phi < 0 || phi > 180.0) {
            if (assertErrors) throw std::runtime_error("PSHB: phase shift out of range [0,180]");
            ok = false;
        }
        if (op.get_switching_frequency() <= 0) {
            if (assertErrors) throw std::runtime_error("PSHB: invalid switching frequency");
            ok = false;
        }
    }
    return ok;
}


// =========================================================================
// Design Requirements
// =========================================================================
DesignRequirements Pshb::process_design_requirements() {
    auto& inputVoltage = get_input_voltage();
    double Vin_nom = inputVoltage.get_nominal().value_or(
        (inputVoltage.get_minimum().value_or(0) + inputVoltage.get_maximum().value_or(0)) / 2.0);

    auto& ops = get_operating_points();
    double Vo = ops[0].get_output_voltages()[0];
    double Io = ops[0].get_output_currents()[0];
    double Fs = ops[0].get_switching_frequency();
    double phi_deg = ops[0].get_phase_shift();

    PsfbRectifierType rectType = get_rectifier_type().value_or(PsfbRectifierType::CENTER_TAPPED);
    double Vd = 0.6;
    computedDiodeVoltageDrop = Vd;

    // Effective duty cycle
    double Deff;
    if (phi_deg > 1e-6) {
        Deff = compute_effective_duty_cycle(phi_deg);
    } else {
        Deff = 0.75;  // Higher default Deff for HB (compensates for halved voltage)
    }
    computedEffectiveDutyCycle = Deff;

    // Turns ratio  (compute_turns_ratio already applies Vin/2 factor)
    double n = compute_turns_ratio(Vin_nom, Vo, Deff, Vd, rectType);

    std::vector<double> turnsRatios;
    turnsRatios.push_back(n);
    for (size_t i = 1; i < ops[0].get_output_voltages().size(); i++) {
        double Voi = ops[0].get_output_voltages()[i];
        turnsRatios.push_back(compute_turns_ratio(Vin_nom, Voi, Deff, Vd, rectType));
    }

    // Output inductance
    double Lo;
    if (get_output_inductance().has_value() && get_output_inductance().value() > 0) {
        Lo = get_output_inductance().value();
    } else {
        double rippleRatio = 0.3;
        Lo = compute_output_inductance(Vo, Deff, Fs, Io, rippleRatio);
    }
    computedOutputInductance = Lo;

    // Series inductance
    double Vhb = Vin_nom * BRIDGE_VOLTAGE_FACTOR;
    double Lr;
    if (get_series_inductance().has_value() && get_series_inductance().value() > 0) {
        Lr = get_series_inductance().value();
    } else {
        double t_loss = 0.02 / Fs;
        double Ip_min = Io / (4.0 * n);
        Lr = (Ip_min > 0) ? Vhb * t_loss / (2.0 * Ip_min) : 2e-6;
        Lr = std::max(Lr, 1e-7);
    }
    computedSeriesInductance = Lr;

    // Magnetizing inductance
    // Im_peak = Vhb * Deff / (4 * Fs * Lm)
    double Io_pri = Io / n;
    double Im_target = 0.1 * Io_pri;
    double Lm;
    if (Im_target > 0) {
        Lm = Vhb * Deff / (4.0 * Fs * Im_target);
    } else {
        Lm = 20.0 * Lr;
    }
    Lm = std::max(Lm, 20.0 * Lr);
    computedMagnetizingInductance = Lm;

    // Build DesignRequirements
    DesignRequirements designRequirements;
    designRequirements.get_mutable_turns_ratios().clear();
    for (auto tr : turnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(roundFloat(tr, 2));
        designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }

    DimensionWithTolerance inductanceWithTolerance;
    inductanceWithTolerance.set_nominal(roundFloat(Lm, 10));
    designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

    if (get_use_leakage_inductance().value_or(false)) {
        std::vector<DimensionWithTolerance> leakageReqs;
        DimensionWithTolerance lrTol;
        lrTol.set_nominal(roundFloat(Lr, 10));
        leakageReqs.push_back(lrTol);
        designRequirements.set_leakage_inductance(leakageReqs);
    }

    return designRequirements;
}


// =========================================================================
// Process operating points
// =========================================================================
std::vector<OperatingPoint> Pshb::process_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    std::vector<OperatingPoint> result;
    auto& inputVoltage = get_input_voltage();

    std::vector<double> inputVoltages;
    if (inputVoltage.get_nominal().has_value())
        inputVoltages.push_back(inputVoltage.get_nominal().value());
    if (inputVoltage.get_minimum().has_value())
        inputVoltages.push_back(inputVoltage.get_minimum().value());
    if (inputVoltage.get_maximum().has_value())
        inputVoltages.push_back(inputVoltage.get_maximum().value());

    std::sort(inputVoltages.begin(), inputVoltages.end());
    inputVoltages.erase(std::unique(inputVoltages.begin(), inputVoltages.end()), inputVoltages.end());

    for (double Vin : inputVoltages) {
        auto op = process_operating_point_for_input_voltage(
            Vin, get_operating_points()[0], turnsRatios, magnetizingInductance);
        result.push_back(op);
    }
    return result;
}

std::vector<OperatingPoint> Pshb::process_operating_points(Magnetic magnetic) {
    auto req = process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    return process_operating_points(turnsRatios, Lm);
}


// =========================================================================
// CORE WAVEFORM GENERATION - Phase-Shifted Half Bridge
// =========================================================================
//
// The PSHB transformer voltage is a 3-level waveform, identical in SHAPE
// to the PSFB but with HALF the amplitude:
//
//   +(Vin/2) during power transfer
//   0        during freewheeling
//   -(Vin/2) during opposite power transfer
//   0        during opposite freewheeling
//
// Primary current shape is the same as PSFB:
//   i_pri = Io_ref ± Im(t)   (trapezoidal / ramp)
//
// The half-bridge nature means:
//   - For the same output power, primary current is ~2x that of full-bridge
//     (because Vpri is halved, so n must be halved, increasing Io/n)
//   - Split capacitors must carry the full AC component of primary current
//   - Only two switches to drive (simpler gate drive)
//
// =========================================================================
OperatingPoint Pshb::process_operating_point_for_input_voltage(
    double inputVoltage,
    const PsfbOperatingPoint& pshbOpPoint,
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    OperatingPoint operatingPoint;

    double Fs = pshbOpPoint.get_switching_frequency();
    double Vin = inputVoltage;
    double Vhb = Vin * BRIDGE_VOLTAGE_FACTOR;  // Vin/2
    double Vo = pshbOpPoint.get_output_voltages()[0];
    double Io = pshbOpPoint.get_output_currents()[0];
    double n = turnsRatios[0];
    double Lm = magnetizingInductance;

    double phi_deg = pshbOpPoint.get_phase_shift();
    double Deff = (phi_deg > 1e-6) ? compute_effective_duty_cycle(phi_deg) : computedEffectiveDutyCycle;

    double period = 1.0 / Fs;
    double Thalf = period / 2.0;
    double t_power = Deff * Thalf;
    double t_free = Thalf - t_power;

    // Primary reflected load current
    double Io_ref = Io / n;

    // Magnetizing current peak (uses Vhb, not Vin)
    double Im_peak = Vhb * Deff / (4.0 * Fs * Lm);

    const int N_samples = 256;
    double dt = Thalf / N_samples;
    int totalSamples = 2 * N_samples + 1;

    std::vector<double> time_full(totalSamples);
    std::vector<double> Vpri_full(totalSamples);
    std::vector<double> Ipri_full(totalSamples);

    int n_power = (int)(t_power / dt);
    if (n_power > N_samples) n_power = N_samples;

    double i_start = Io_ref - Im_peak;
    double i_end_power = Io_ref + Im_peak;

    // Positive half-cycle
    for (int k = 0; k <= N_samples; ++k) {
        double t = k * dt;
        time_full[k] = t;

        if (k <= n_power) {
            Vpri_full[k] = Vhb;  // +Vin/2
            double frac = (n_power > 0) ? (double)k / n_power : 0;
            Ipri_full[k] = i_start + (i_end_power - i_start) * frac;
        } else {
            Vpri_full[k] = 0.0;
            Ipri_full[k] = i_end_power;
        }
    }

    // Negative half-cycle by antisymmetry
    for (int k = 1; k <= N_samples; ++k) {
        time_full[N_samples + k] = Thalf + k * dt;
        Ipri_full[N_samples + k] = -Ipri_full[k];
        Vpri_full[N_samples + k] = -Vpri_full[k];
    }

    // ---- Primary winding excitation ----
    {
        Waveform currentWaveform;
        currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        currentWaveform.set_data(Ipri_full);
        currentWaveform.set_time(time_full);

        Waveform voltageWaveform;
        voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        voltageWaveform.set_data(Vpri_full);
        voltageWaveform.set_time(time_full);

        auto excitation = complete_excitation(currentWaveform, voltageWaveform,
                                              Fs, "Primary");
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    // ---- Secondary winding excitation(s) ----
    for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
        double ni = turnsRatios[secIdx];

        std::vector<double> iSecData(totalSamples);
        std::vector<double> vSecData(totalSamples);

        for (int k = 0; k < totalSamples; ++k) {
            vSecData[k] = Vpri_full[k] / ni;
            iSecData[k] = ni * Ipri_full[k];
        }

        Waveform secCurrentWfm;
        secCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        secCurrentWfm.set_data(iSecData);
        secCurrentWfm.set_time(time_full);

        Waveform secVoltageWfm;
        secVoltageWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        secVoltageWfm.set_data(vSecData);
        secVoltageWfm.set_time(time_full);

        auto excitation = complete_excitation(secCurrentWfm, secVoltageWfm,
                                              Fs,
                                              "Secondary " + std::to_string(secIdx));
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    OperatingConditions conditions;
    conditions.set_ambient_temperature(pshbOpPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    return operatingPoint;
}


// =========================================================================
// SPICE Circuit Generation - Half Bridge
// =========================================================================
std::string Pshb::generate_ngspice_circuit(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t inputVoltageIndex,
    size_t operatingPointIndex)
{
    auto& inputVoltageSpec = get_input_voltage();
    auto& ops = get_operating_points();

    std::vector<double> inputVoltages;
    if (inputVoltageSpec.get_nominal().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_nominal().value());
    if (inputVoltageSpec.get_minimum().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_minimum().value());
    if (inputVoltageSpec.get_maximum().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_maximum().value());

    double Vin = inputVoltages[std::min(inputVoltageIndex, inputVoltages.size() - 1)];
    auto& pshbOp = ops[std::min(operatingPointIndex, ops.size() - 1)];

    double Fs = pshbOp.get_switching_frequency();
    double period = 1.0 / Fs;
    double halfPeriod = period / 2.0;
    double deadTime = computedDeadTime;
    double tOn = halfPeriod - deadTime;
    double Vo = pshbOp.get_output_voltages()[0];
    double Io = pshbOp.get_output_currents()[0];
    double n = turnsRatios[0];
    double Lm = magnetizingInductance;
    double Lr = computedSeriesInductance;
    double Lo = computedOutputInductance;

    double phi_deg = pshbOp.get_phase_shift();
    double Deff = (phi_deg > 1e-6) ? compute_effective_duty_cycle(phi_deg) : computedEffectiveDutyCycle;

    int periodsToExtract = get_num_periods_to_extract();
    int steadyStatePeriods = get_num_steady_state_periods();
    int numPeriodsTotal = steadyStatePeriods + periodsToExtract;
    double simTime = numPeriodsTotal * period;
    double startTime = steadyStatePeriods * period;
    double stepTime = period / 500;

    std::ostringstream circuit;

    circuit << "* Phase-Shifted Half Bridge (PSHB) Converter\n";
    circuit << "* Vin=" << Vin << "V, Vo=" << Vo << "V, Fs=" << (Fs/1e3)
            << "kHz, Deff=" << Deff << "\n";
    circuit << "* n=" << n << ", Lr=" << (Lr*1e6) << "uH, Lm=" << (Lm*1e6)
            << "uH, Lo=" << (Lo*1e6) << "uH\n";
    circuit << "* Bridge voltage factor = 0.5 (half-bridge)\n\n";

    circuit << ".model SW1 SW(Ron=10m Roff=10Meg Vt=2.5)\n";
    circuit << ".model DIDEAL D(Is=1e-14 N=0.001)\n\n";

    circuit << "Vdc vin_dc 0 " << Vin << "\n\n";

    // Split capacitor divider
    circuit << "* Split capacitor bus (Vin/2 each)\n";
    circuit << "C_split_hi vin_dc mid_cap 100u ic=" << (Vin/2.0) << "\n";
    circuit << "C_split_lo mid_cap 0 100u ic=" << (Vin/2.0) << "\n\n";

    // Half-bridge: QA (high-side), QB (low-side)
    circuit << "* Half-bridge leg: QA (hi), QB (lo)\n";
    circuit << "Vpwm_A pwm_A 0 PULSE(0 5 0 10n 10n "
            << std::scientific << tOn << " " << period << std::fixed << ")\n";
    circuit << "Vpwm_B pwm_B 0 PULSE(0 5 "
            << std::scientific << halfPeriod << std::fixed
            << " 10n 10n "
            << std::scientific << tOn << " " << period << std::fixed << ")\n\n";

    circuit << "SA vin_dc mid_sw pwm_A 0 SW1\n";
    circuit << "DA 0 mid_sw DIDEAL\n";
    circuit << "SB mid_sw 0 pwm_B 0 SW1\n";
    circuit << "DB mid_sw vin_dc DIDEAL\n\n";

    // Primary current sense
    circuit << "Vpri_sense mid_sw pri_lr 0\n\n";

    // Series inductance (leakage + external)
    circuit << "L_series pri_lr trafo_pri " << std::scientific << Lr << "\n\n";

    // Transformer: primary between switch mid-point and cap mid-point
    double Ls_sec = Lm / (n * n);
    circuit << "* Transformer primary between mid_sw and mid_cap\n";
    circuit << "L_pri trafo_pri mid_cap " << std::scientific << Lm << "\n";
    circuit << "L_sec sec_a sec_b " << std::scientific << Ls_sec << "\n";
    circuit << "K_trafo L_pri L_sec 0.9999\n\n";

    // Output rectifier (full-bridge diode)
    circuit << "* Output full-bridge rectifier\n";
    circuit << "D_r1 sec_a out_rect DIDEAL\n";
    circuit << "D_r2 sec_b out_rect DIDEAL\n";
    circuit << "D_r3 out_gnd sec_a DIDEAL\n";
    circuit << "D_r4 out_gnd sec_b DIDEAL\n\n";

    // Output filter
    circuit << "L_out out_rect out_node " << std::scientific << Lo << "\n";
    circuit << "R_load out_node out_gnd " << std::fixed << (Vo / Io) << "\n";
    circuit << "C_out out_node out_gnd 100u\n\n";

    circuit << ".tran " << std::scientific << stepTime
            << " " << simTime << " " << startTime << "\n";
    circuit << ".end\n";

    return circuit.str();
}


std::vector<OperatingPoint> Pshb::simulate_and_extract_operating_points(
    const std::vector<double>& turnsRatios, double magnetizingInductance)
{
    return process_operating_points(turnsRatios, magnetizingInductance);
}

std::vector<ConverterWaveforms> Pshb::simulate_and_extract_topology_waveforms(
    const std::vector<double>& turnsRatios, double magnetizingInductance)
{
    return {};
}


// =========================================================================
// AdvancedPshb
// =========================================================================
Inputs AdvancedPshb::process() {
    auto designRequirements = process_design_requirements();

    designRequirements.get_mutable_turns_ratios().clear();
    for (auto n : desiredTurnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(n);
        designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }

    DimensionWithTolerance LmTol;
    LmTol.set_nominal(desiredMagnetizingInductance);
    designRequirements.set_magnetizing_inductance(LmTol);

    if (desiredSeriesInductance.has_value())
        set_computed_series_inductance(desiredSeriesInductance.value());
    if (desiredOutputInductance.has_value())
        set_computed_output_inductance(desiredOutputInductance.value());

    auto ops = process_operating_points(desiredTurnsRatios, desiredMagnetizingInductance);

    Inputs inputs;
    inputs.set_design_requirements(designRequirements);
    inputs.set_operating_points(ops);
    return inputs;
}

} // namespace OpenMagnetics
