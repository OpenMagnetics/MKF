#include "converter_models/PhaseShiftedHalfBridge.h"
#include "converter_models/PwmBridgeSolver.h"
#include "physical_models/LeakageInductance.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include "support/Exceptions.h"

namespace OpenMagnetics {

Pshb::Pshb(const json& j) {
    from_json(j, *static_cast<PhaseShiftedHalfBridge*>(this));
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
                                    double Vd, BRectifierType rectType) {
    double Vhb = Vin * BRIDGE_VOLTAGE_FACTOR;  // Vin/2
    switch (rectType) {
        case BRectifierType::CENTER_TAPPED:
            return Vhb * Deff / n - Vd;
        case BRectifierType::CURRENT_DOUBLER:
            return Vhb * Deff / (2.0 * n) - Vd;
        case BRectifierType::FULL_BRIDGE:
        default:
            return Vhb * Deff / n - 2.0 * Vd;
    }
}

// =========================================================================
// Static: turns ratio for target output voltage
// =========================================================================
double Pshb::compute_turns_ratio(double Vin, double Vo, double Deff,
                                 double Vd, BRectifierType rectType) {
    double Vhb = Vin * BRIDGE_VOLTAGE_FACTOR;
    switch (rectType) {
        case BRectifierType::CENTER_TAPPED:
            return Vhb * Deff / (Vo + Vd);
        case BRectifierType::CURRENT_DOUBLER:
            return Vhb * Deff / (2.0 * (Vo + Vd));
        case BRectifierType::FULL_BRIDGE:
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
// Static: diode forward drop at a given conduction current
// Matches the SPICE diode model used in generate_ngspice_circuit:
//     .model DIDEAL D(IS=1e-12 RS=0.005 BV=1000 IBV=1e-12)
// Vd(I) = Vt · ln(I / IS) + I · RS,    Vt = kT/q ≈ 0.025 V at 300 K
// RS=0.005 matches a representative Schottky (Vf ≈ 0.45 V at 40 A);
// RS=0.05 gave a physically wrong 3.25 V drop at 50 A.
// =========================================================================
double Pshb::compute_diode_drop_at_current(double Io) {
    constexpr double DIODE_IS = 1e-12;
    constexpr double DIODE_RS = 0.005;
    constexpr double VT       = 0.025;
    if (Io <= 0)
        throw std::runtime_error("PSHB: compute_diode_drop_at_current requires Io > 0");
    return VT * std::log(Io / DIODE_IS) + Io * DIODE_RS;
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

    BRectifierType rectType = get_rectifier_type().value_or(BRectifierType::CENTER_TAPPED);
    // Diode drop computed from the same model used in the SPICE netlist;
    // the prior fixed 0.6 V default caused ~25 % Vo droop at 50 A.
    double Vd = compute_diode_drop_at_current(Io);
    computedDiodeVoltageDrop = Vd;

    double D_cmd = (phi_deg > 1e-6) ? compute_effective_duty_cycle(phi_deg) : 0.75;
    double Vhb = Vin_nom * BRIDGE_VOLTAGE_FACTOR;

    // Lr seed (see PSFB for rationale).
    double Lr;
    if (get_series_inductance().has_value() && get_series_inductance().value() > 0) {
        Lr = get_series_inductance().value();
    } else {
        double n_seed = compute_turns_ratio(Vin_nom, Vo, D_cmd, Vd, rectType);
        double Lr_dcl_cap = (Io > 0)
            ? 0.02 * std::max(n_seed, 0.1) * Vhb / (4.0 * Io * Fs)
            : 2e-6;
        Lr = std::min(2e-6, Lr_dcl_cap);
        Lr = std::max(Lr, 1e-7);
    }
    computedSeriesInductance = Lr;

    // Iterate n / D_eff with duty-cycle-loss correction (Vbus = Vin/2 here).
    double n = compute_turns_ratio(Vin_nom, Vo, D_cmd, Vd, rectType);
    double Deff = D_cmd;
    for (int iter = 0; iter < 8; ++iter) {
        double dcl = PwmBridgeSolver::compute_duty_cycle_loss(Vhb, Lr, Io, n, Fs);
        Deff = PwmBridgeSolver::compute_effective_duty_cycle(D_cmd, dcl);
        double n_new = (Deff > 1e-3)
                       ? compute_turns_ratio(Vin_nom, Vo, Deff, Vd, rectType)
                       : n;
        if (std::abs(n_new - n) < 1e-3 * std::max(n, 1.0)) { n = n_new; break; }
        n = n_new;
    }
    computedEffectiveDutyCycle = Deff;

    std::vector<double> turnsRatios;
    turnsRatios.push_back(n);
    for (size_t i = 1; i < ops[0].get_output_voltages().size(); i++) {
        double Voi = ops[0].get_output_voltages()[i];
        turnsRatios.push_back(compute_turns_ratio(Vin_nom, Voi, Deff, Vd, rectType));
    }

    double Lo;
    if (get_output_inductance().has_value() && get_output_inductance().value() > 0) {
        Lo = get_output_inductance().value();
    } else {
        double rippleRatio = 0.3;
        Lo = compute_output_inductance(Vo, Deff, Fs, Io, rippleRatio);
    }
    computedOutputInductance = Lo;

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

    designRequirements.set_topology(MAS::Topology::PHASE_SHIFTED_HALF_BRIDGE_CONVERTER);
    designRequirements.set_isolation_sides(
        Topology::create_isolation_sides(ops[0].get_output_currents().size(), false));

    return designRequirements;
}


// =========================================================================
// Process operating points
// =========================================================================
std::vector<OperatingPoint> Pshb::process_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    extraLoVoltageWaveforms.clear();
    extraLoCurrentWaveforms.clear();
    extraLo2VoltageWaveforms.clear();
    extraLo2CurrentWaveforms.clear();
    extraLrVoltageWaveforms.clear();
    extraLrCurrentWaveforms.clear();

    perOpName.clear();
    perOpDutyCycleLoss.clear();
    perOpEffectiveDutyCycle.clear();
    perOpZvsMarginLagging.clear();
    perOpZvsLoadThreshold.clear();
    perOpResonantTransitionTime.clear();
    perOpPrimaryPeakCurrent.clear();

    std::vector<OperatingPoint> result;
    auto& inputVoltage = get_input_voltage();

    std::vector<std::pair<double, std::string>> labelled;
    if (inputVoltage.get_nominal().has_value())
        labelled.push_back({inputVoltage.get_nominal().value(), "Nominal"});
    if (inputVoltage.get_minimum().has_value())
        labelled.push_back({inputVoltage.get_minimum().value(), "Min"});
    if (inputVoltage.get_maximum().has_value())
        labelled.push_back({inputVoltage.get_maximum().value(), "Max"});
    std::sort(labelled.begin(), labelled.end(),
        [](const auto& a, const auto& b) { return a.first < b.first; });
    auto last = std::unique(labelled.begin(), labelled.end(),
        [](const auto& a, const auto& b) { return a.first == b.first; });
    labelled.erase(last, labelled.end());

    for (const auto& [Vin, name] : labelled) {
        auto op = process_operating_point_for_input_voltage(
            Vin, get_operating_points()[0], turnsRatios, magnetizingInductance);
        result.push_back(op);

        perOpName.push_back(name);
        perOpDutyCycleLoss.push_back(lastDutyCycleLoss);
        perOpEffectiveDutyCycle.push_back(lastEffectiveDutyCycle);
        perOpZvsMarginLagging.push_back(lastZvsMarginLagging);
        perOpZvsLoadThreshold.push_back(lastZvsLoadThreshold);
        perOpResonantTransitionTime.push_back(lastResonantTransitionTime);
        perOpPrimaryPeakCurrent.push_back(lastPrimaryPeakCurrent);
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
// CORE WAVEFORM GENERATION — PSHB Analytical Model (sub-interval solver)
// =========================================================================
//
// Same three-segment model as PSFB but with Vbus = Vin/2 (the split-cap
// midpoint sets the primary swing to ±Vin/2). All physics — duty-cycle
// loss, ZVS energy balance, output-inductor ripple — carry over with
// the bus voltage substitution.
//
// Refs: Pinheiro & Barbi, IEEE TPE 8(4) 1993; Sabate APEC 1990 (duty
// loss derivation reused with Vbus → Vin/2).
// =========================================================================
OperatingPoint Pshb::process_operating_point_for_input_voltage(
    double inputVoltage,
    const MAS::PshbOperatingPoint&pshbOpPoint,
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    OperatingPoint operatingPoint;

    double Fs = pshbOpPoint.get_switching_frequency();
    // Guard unconditionally before any 1/Fs division (direct callers bypass
    // run_checks — same rationale as PSFB).
    if (Fs <= 0)
        throw std::runtime_error("PSHB: switching frequency must be positive");
    double Vin = inputVoltage;
    double Vhb = Vin * BRIDGE_VOLTAGE_FACTOR;
    double Vo  = pshbOpPoint.get_output_voltages()[0];
    double Io  = pshbOpPoint.get_output_currents()[0];
    double n   = turnsRatios[0];
    double Lm  = magnetizingInductance;
    if (computedSeriesInductance <= 0)
        throw std::runtime_error("PSHB: series (resonant) inductance must be positive; "
                                 "it was not computed/provided");
    double Lr  = computedSeriesInductance;

    BRectifierType rectType = get_rectifier_type().value_or(BRectifierType::CENTER_TAPPED);

    double phi_deg = pshbOpPoint.get_phase_shift();
    double D_cmd = (phi_deg > 1e-6) ? compute_effective_duty_cycle(phi_deg)
                                    : computedEffectiveDutyCycle;
    if (D_cmd <= 0)
        throw std::runtime_error("PSHB: effective duty cycle is non-positive; provide a phase "
                                 "shift or valid design requirements (no default duty substituted)");

    double dcl_duty = PwmBridgeSolver::compute_duty_cycle_loss(Vhb, Lr, Io, n, Fs);
    double Deff = PwmBridgeSolver::compute_effective_duty_cycle(D_cmd, dcl_duty);
    lastDutyCycleLoss = dcl_duty;
    lastEffectiveDutyCycle = Deff;

    double period = 1.0 / Fs;
    double Thalf = period / 2.0;

    // Im integrates +Vhb during BOTH commutation and active intervals
    // (total +Vhb time per half-cycle = D_cmd·Thalf), so use D_cmd.
    double Im_peak = Vhb * D_cmd / (4.0 * Fs * Lm);

    double Lo = computedOutputInductance;
    double dILo = (Lo > 0) ? Vo * (1.0 - Deff) / (Fs * Lo) : 0.0;
    double Io_in_inductor = (rectType == BRectifierType::CURRENT_DOUBLER) ? Io / 2.0 : Io;
    double ILo_min = Io_in_inductor - dILo / 2.0;
    double ILo_max = Io_in_inductor + dILo / 2.0;

    auto segs = PwmBridgeSolver::build_first_half_cycle(period, Vhb, D_cmd, dcl_duty);

    double t_dcl    = std::min(dcl_duty * Thalf, Thalf);
    double t_active = std::min(D_cmd * Thalf,    Thalf);
    if (t_active < t_dcl) t_active = t_dcl;
    double Vsec_pk = Vhb / n;

    const int N_samples = 256;
    double dt = Thalf / N_samples;
    int totalSamples = 2 * N_samples + 1;

    std::vector<double> time_full(totalSamples);
    std::vector<double> Vpri_full(totalSamples);
    std::vector<double> Ipri_full(totalSamples);
    std::vector<double> Im_full(totalSamples);
    std::vector<double> ILo_full(totalSamples);
    std::vector<double> VLo_full(totalSamples);

    double Im0 = -Im_peak;

    std::vector<double> Im_at_seg_start(segs.size());
    double Im_running = Im0;
    for (size_t si = 0; si < segs.size(); ++si) {
        Im_at_seg_start[si] = Im_running;
        double dur = segs[si].t_end - segs[si].t_start;
        Im_running += segs[si].Vp * dur / Lm;
    }

    auto seg_index_at = [&](double t) -> size_t {
        for (size_t si = 0; si < segs.size(); ++si) {
            if (t < segs[si].t_end - 1e-15) return si;
        }
        return segs.size() - 1;
    };

    for (int k = 0; k <= N_samples; ++k) {
        double t = k * dt;
        time_full[k] = t;
        size_t si = seg_index_at(t);
        const auto& s = segs[si];
        double Vp = s.Vp;
        Vpri_full[k] = Vp;

        Im_full[k] = Im_at_seg_start[si] + Vp * (t - s.t_start) / Lm;
        if (Im_full[k] >  Im_peak) Im_full[k] =  Im_peak;
        if (Im_full[k] < -Im_peak) Im_full[k] = -Im_peak;

        // Output-inductor — exact piecewise linear: trough at t_dcl, peak
        // at t_active (see PSFB analogous code).
        double slope_neg = (Lo > 0) ? -Vo / Lo : 0.0;
        double slope_pos = (Lo > 0) ?  (Vsec_pk - Vo) / Lo : 0.0;
        double ILo_at_zero = ILo_min - slope_neg * t_dcl;
        if (t < t_dcl) {
            ILo_full[k] = ILo_at_zero + slope_neg * t;
            VLo_full[k] = -Vo;
        } else if (t < t_active) {
            ILo_full[k] = ILo_min + slope_pos * (t - t_dcl);
            VLo_full[k] = Vsec_pk - Vo;
        } else {
            ILo_full[k] = ILo_max + slope_neg * (t - t_active);
            VLo_full[k] = -Vo;
        }

        if (s.is_commutation && (s.t_end - s.t_start) > 1e-15) {
            double cfrac = (t - s.t_start) / (s.t_end - s.t_start);
            double ipri_start = -ILo_max / n;
            double ipri_end   =  ILo_min / n;
            Ipri_full[k] = ipri_start + (ipri_end - ipri_start) * cfrac + Im_full[k];
        } else if (s.is_active) {
            Ipri_full[k] = ILo_full[k] / n + Im_full[k];
        } else {
            // Freewheel: PSHB has an NPC-style bridge that shorts the
            // primary through one inner switch + the mid-cap stack during
            // freewheel.  Effective resistance ≈ 2·RON_mosfet (one inner
            // switch + one outer switch path), so:
            //     τ = Lr / (2 · RON)
            // For SPICE RON=0.01 Ω and typical Lr = 1-10 µH this gives
            // τ = 50-500 µs ≫ dur_fw (~1.5 µs at 100 kHz), so the primary
            // current is essentially flat through freewheel with a small
            // exponential droop.  The old "tau = 1.0·dur_fw" linear decay
            // ran from ILo_max/n + Im_peak to Im_peak across the *entire*
            // freewheel interval — too aggressive; was tuned against a
            // SPICE model with ~25 % Vo droop, where the coincidental
            // amplitude mismatch was masking the shape mismatch.
            constexpr double RON_PER_SWITCH = 0.01;   // matches SPICE SW1 model
            double tau_fw = (Lr > 0) ? Lr / (2.0 * RON_PER_SWITCH) : 1.0;
            double t_fw = t - t_active;
            double ipri_start_fw = ILo_max / n + Im_peak;
            double ipri_end_fw   = Im_peak;
            double f = 1.0 - std::exp(-t_fw / tau_fw);
            Ipri_full[k] = ipri_start_fw + (ipri_end_fw - ipri_start_fw) * f;
        }
    }

    for (int k = 1; k <= N_samples; ++k) {
        time_full[N_samples + k] = Thalf + k * dt;
        Ipri_full[N_samples + k] = -Ipri_full[k];
        Vpri_full[N_samples + k] = -Vpri_full[k];
        Im_full[N_samples + k]   = -Im_full[k];
        ILo_full[N_samples + k]  = ILo_full[k];
        VLo_full[N_samples + k]  = VLo_full[k];
    }

    double Ip_pk = std::abs(ILo_max / n + Im_peak);
    lastPrimaryPeakCurrent = Ip_pk;
    lastZvsMarginLagging = PwmBridgeSolver::compute_zvs_margin_lagging(
        Lr, Ip_pk, mosfetCoss, Vhb);
    double Ip_zvs_min = PwmBridgeSolver::compute_zvs_primary_current_min(
        Lr, mosfetCoss, Vhb);
    double Io_thr = (Ip_zvs_min - Im_peak) * n;
    lastZvsLoadThreshold = std::max(Io_thr, 0.0);
    lastResonantTransitionTime = PwmBridgeSolver::compute_resonant_transition_time(
        Lr, mosfetCoss);

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

    for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
        double ni = turnsRatios[secIdx];
        if (ni <= 0) continue;

        if (rectType == BRectifierType::CENTER_TAPPED) {
            double Vsec_pk = Vhb / ni;
            std::vector<double> v1(totalSamples), i1(totalSamples);
            std::vector<double> v2(totalSamples), i2(totalSamples);
            for (int k = 0; k < totalSamples; ++k) {
                bool positive_half = (k <= N_samples);
                double vpri = Vpri_full[k];
                double iLo_k = ILo_full[k];
                if (positive_half) {
                    v1[k] = (vpri > 0) ? Vsec_pk : (vpri < 0 ? -Vsec_pk : 0.0);
                    i1[k] = (vpri > 0) ? iLo_k : 0.0;
                    v2[k] = -v1[k]; i2[k] = 0.0;
                } else {
                    v2[k] = (vpri < 0) ? Vsec_pk : (vpri > 0 ? -Vsec_pk : 0.0);
                    i2[k] = (vpri < 0) ? iLo_k : 0.0;
                    v1[k] = -v2[k]; i1[k] = 0.0;
                }
            }
            Waveform i1Wfm; i1Wfm.set_ancillary_label(WaveformLabel::CUSTOM);
            i1Wfm.set_data(i1); i1Wfm.set_time(time_full);
            Waveform v1Wfm; v1Wfm.set_ancillary_label(WaveformLabel::CUSTOM);
            v1Wfm.set_data(v1); v1Wfm.set_time(time_full);
            operatingPoint.get_mutable_excitations_per_winding().push_back(
                complete_excitation(i1Wfm, v1Wfm, Fs,
                                    "Secondary " + std::to_string(secIdx) + "a"));
            Waveform i2Wfm; i2Wfm.set_ancillary_label(WaveformLabel::CUSTOM);
            i2Wfm.set_data(i2); i2Wfm.set_time(time_full);
            Waveform v2Wfm; v2Wfm.set_ancillary_label(WaveformLabel::CUSTOM);
            v2Wfm.set_data(v2); v2Wfm.set_time(time_full);
            operatingPoint.get_mutable_excitations_per_winding().push_back(
                complete_excitation(i2Wfm, v2Wfm, Fs,
                                    "Secondary " + std::to_string(secIdx) + "b"));
        } else {
            std::vector<double> vSec(totalSamples), iSec(totalSamples);
            for (int k = 0; k < totalSamples; ++k) {
                vSec[k] = Vpri_full[k] / ni;
                double sign = (Vpri_full[k] > 1e-6) ? 1.0
                            : (Vpri_full[k] < -1e-6) ? -1.0 : 0.0;
                iSec[k] = sign * ILo_full[k];
            }
            Waveform iWfm; iWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            iWfm.set_data(iSec); iWfm.set_time(time_full);
            Waveform vWfm; vWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            vWfm.set_data(vSec); vWfm.set_time(time_full);
            operatingPoint.get_mutable_excitations_per_winding().push_back(
                complete_excitation(iWfm, vWfm, Fs,
                                    "Secondary " + std::to_string(secIdx)));
        }
    }

    OperatingConditions conditions;
    conditions.set_ambient_temperature(pshbOpPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    {
        Waveform loVWfm;
        loVWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        loVWfm.set_data(VLo_full);
        loVWfm.set_time(time_full);
        Waveform loIWfm;
        loIWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        loIWfm.set_data(ILo_full);
        loIWfm.set_time(time_full);
        extraLoVoltageWaveforms.push_back(loVWfm);
        extraLoCurrentWaveforms.push_back(loIWfm);

        if (rectType == BRectifierType::CURRENT_DOUBLER) {
            std::vector<double> ILo2(totalSamples), VLo2(totalSamples);
            for (int k = 0; k < totalSamples; ++k) {
                int kshift = (k + N_samples) % (totalSamples - 1);
                ILo2[k] = ILo_full[kshift];
                VLo2[k] = VLo_full[kshift];
            }
            Waveform lo2VWfm; lo2VWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            lo2VWfm.set_data(VLo2); lo2VWfm.set_time(time_full);
            Waveform lo2IWfm; lo2IWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            lo2IWfm.set_data(ILo2); lo2IWfm.set_time(time_full);
            extraLo2VoltageWaveforms.push_back(lo2VWfm);
            extraLo2CurrentWaveforms.push_back(lo2IWfm);
        }
    }
    {
        Waveform lrVWfm;
        lrVWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        lrVWfm.set_data(Vpri_full);
        lrVWfm.set_time(time_full);
        Waveform lrIWfm;
        lrIWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        lrIWfm.set_data(Ipri_full);
        lrIWfm.set_time(time_full);
        extraLrVoltageWaveforms.push_back(lrVWfm);
        extraLrCurrentWaveforms.push_back(lrIWfm);
    }

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
    // All SPICE-side knobs are pulled from spice_config(). The PSHB
    // registry entry (Topology.cpp) was split out from the previously
    // shared `psb` entry and tuned to match this file's historical
    // hardcoded netlist byte-for-byte; this refactor is behaviour-
    // preserving.
    const auto cfg = spice_config();

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
    double stepTime = period / 200;

    BRectifierType rectType = get_rectifier_type().value_or(BRectifierType::CENTER_TAPPED);

    size_t numOutputs = std::min(turnsRatios.size(), pshbOp.get_output_voltages().size());
    if (numOutputs == 0) numOutputs = 1;

    std::ostringstream circuit;

    circuit << "* Phase-Shifted Half Bridge (PSHB) — three-level Pinheiro-Barbi NPC\n";
    circuit << "* Vin=" << Vin << "V (primary swings +Vin/2 / 0 / -Vin/2 via 3-level NPC), Vo=" << Vo
            << "V, Fs=" << (Fs/1e3) << "kHz, Deff=" << Deff
            << ", outputs=" << numOutputs << "\n";
    circuit << "* n=" << n << ", Lr=" << (Lr*1e6) << "uH, Lm=" << (Lm*1e6)
            << "uH, Lo=" << (Lo*1e6) << "uH, rect=" << static_cast<int>(rectType) << "\n";
    circuit << "* Bridge: single 3-level NPC leg (4 stacked S-switches + 2 NPC clamp diodes\n";
    circuit << "* + per-switch RC snubbers). Same modeling approach as DAB/PSFB (real S-switches\n";
    circuit << "* with PWM gates and DIDEAL antiparallel body diodes). DC1/DC2 anchor the inner\n";
    circuit << "* nodes to mid_cap (Vin/2) during the freewheel intervals when only one inner\n";
    circuit << "* switch conducts. Transformer primary connects between point 'a' (leg midpoint)\n";
    circuit << "* and 'b' (split-cap midpoint), so vab swings +Vin/2 / 0 / -Vin/2 / 0.\n";
    circuit << "* Reference: Pinheiro & Barbi, IEEE TPE 8(4) 1993; Barbi & Pottker, Springer 2018.\n\n";

    const auto bridgeMode = get_bridge_simulation_mode();

    // PWM SCHEME — Pinheiro-Barbi 1993 / Barbi & Pottker 2018, ch. 9.
    // Single 3-level NPC leg between vin_dc and 0:
    //
    //   vin_dc --[S1]-- nH --[S2]-- a --[S3]-- nL --[S4]-- 0
    //                    ^                       ^
    //               DC1: mid_cap → nH      DC2: nL → mid_cap
    //
    // Transformer primary (with Lr in series) connects between point 'a'
    // (the leg midpoint) and 'b' = mid_cap (the split-cap midpoint).
    //
    // Gating (Deff = active duty per half-period, t_act = Deff·Thalf):
    //   S1 (outer top): ON  during [0, t_act]                    → a = Vin
    //   S2 (inner top): ON  during [0, halfPeriod]               → first half
    //   S3 (inner bot): ON  during [halfPeriod, period]          → second half
    //   S4 (outer bot): ON  during [halfPeriod, halfPeriod+t_act] → a = 0
    //
    // Resulting waveform at point 'a' (relative to ground):
    //   [0,         t_act]              : S1+S2 on  → a = Vin       (active +)
    //   [t_act,     Thalf]              : only S2 on → a = Vin/2     (clamp via DC2)
    //   [Thalf,     Thalf+t_act]        : S3+S4 on  → a = 0          (active −)
    //   [Thalf+t_act, Ts]               : only S3 on → a = Vin/2     (clamp via DC1)
    //
    // vab = a − mid_cap = +Vin/2, 0, −Vin/2, 0 — matches the analytical
    // model (Vhb = Vin/2, Deff active fraction per half-cycle).
    const double t_act = Deff * halfPeriod;

    circuit << ".model DIDEAL D(IS=" << cfg.diodeIS << " RS=" << cfg.diodeRS;
    if (!cfg.diodeExtra.empty()) {
        circuit << " " << cfg.diodeExtra;
    }
    circuit << ")\n";
    // SW1 is referenced by BOTH bridge modes — VOLTAGE_CONTROLLED_SWITCH uses
    // it for the four NPC switches, and BEHAVIORAL_PULSE also instantiates
    // S1..S4 SW1 elements (lines ~794-800) as the §8a.5-correct probe
    // insertion path. Without this model definition the netlist parses as
    // "can't find model 'sw1'" and the simulation aborts before .tran even
    // starts, surfacing on the frontend as an unrecoverable "Simulated" hang.
    circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
            << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
    circuit << "\n";

    circuit << "Vdc vin_dc 0 " << Vin << "\n\n";

    if (bridgeMode == BridgeSimulationMode::BEHAVIORAL_PULSE) {
        // -----------------------------------------------------------------
        // Behavioral PWL bridge — fast path (MagneticAdviser default).
        //
        // The 3-level NPC waveform at point 'a' is:
        //   [0,             t_act]              : Vin       (S1+S2 active +)
        //   [t_act,         halfPeriod]         : Vin/2     (clamp via DC2)
        //   [halfPeriod,    halfPeriod+t_act]   : 0         (S3+S4 active −)
        //   [halfPeriod+t_act, period]          : Vin/2     (clamp via DC1)
        // We emit this as a piecewise-linear independent V-source on
        // bridge_a with deadTime-length slew between levels (approximating
        // the ZVS resonant transitions that SW1 mode would resolve at
        // sub-step granularity), and we hard-pin mid_cap to Vin/2 (the
        // analytical model already treats the split-cap midpoint as a
        // stiff reference, so no information is lost). No SW1 / NPC
        // diodes / snubbers / split-cap dynamics — saves ~10× ngspice
        // runtime in WASM (root cause of MKF M2).
        const double slew = computedDeadTime;
        if (slew <= 0.0) {
            throw std::runtime_error(
                "Pshb BEHAVIORAL_PULSE bridge: computedDeadTime must be > 0 "
                "to define the leg-voltage slew interval. Got "
                + std::to_string(slew));
        }
        if (2.0 * slew >= t_act) {
            throw std::runtime_error(
                "Pshb BEHAVIORAL_PULSE bridge: 2*deadTime ("
                + std::to_string(2.0 * slew)
                + " s) must be < t_act ("
                + std::to_string(t_act)
                + " s). Reduce dead time or increase Deff.");
        }
        if (2.0 * slew >= (halfPeriod - t_act)) {
            throw std::runtime_error(
                "Pshb BEHAVIORAL_PULSE bridge: 2*deadTime ("
                + std::to_string(2.0 * slew)
                + " s) must be < (halfPeriod - t_act) ("
                + std::to_string(halfPeriod - t_act)
                + " s). Reduce dead time or decrease Deff.");
        }

        const double halfVin = 0.5 * Vin;

        // SPICE_PROBE_HANDOFF.md Bug 2 fix:
        //
        // The old form was `Vbridge_a leg_a_src 0 PWL(...)` (an
        // independent 3-level source) + `Vq1_sense leg_a_src bridge_a
        // 0`. Vdc was decoupled from the bridge by construction, so
        // `i(Vdc) = 0` and `i(Vq1_sense)` measured the leg-internal
        // current, not the actual DC-bus draw the §8a.5 probe rule (B)
        // requires. The half-bridge split-cap mid-point made this
        // approximation worse than PSFB's.
        //
        // Replace the independent PWL leg with a 4-switch NPC leg
        // (mirrors the SW1 branch below): split capacitors create the
        // mid_cap = Vin/2 reference, S1/S4 are the outer switches
        // routed through `Vq1_sense` on `vin_dc`, S2/S3 are the inner
        // switches that clamp `bridge_a` to mid_cap during the
        // non-active half of each cycle. Body diodes and snubbers are
        // intentionally omitted to keep the BEHAVIORAL path lean — only
        // S1's drain runs through `Vq1_sense`, so the measured
        // i(Vq1_sense) is the outer-top switch channel current, which
        // is the §8a.5-compliant input-current probe.
        circuit << "C_split_hi vin_dc mid_cap 470u ic=" << halfVin << "\n";
        circuit << "C_split_lo mid_cap 0 470u ic=" << halfVin << "\n\n";

        // PWM signals: S1 active for t_act in first half, S4 active for
        // t_act in second half (outer switches); S2 active throughout
        // the first half, S3 throughout the second half (inner clamp).
        auto pulse = [&](const std::string& name, double delay, double width) {
            circuit << "V" << name << " " << name << " 0 PULSE(0 5 "
                    << std::scientific << delay << " " << slew << " " << slew
                    << " " << width << " " << period << std::fixed << ")\n";
        };
        pulse("pwm_S1", 0.0,        t_act);
        pulse("pwm_S2", 0.0,        halfPeriod);
        pulse("pwm_S3", halfPeriod, halfPeriod);
        pulse("pwm_S4", halfPeriod, t_act);

        circuit << "Vq1_sense vin_dc qa_drain 0\n";
        circuit << "S1 qa_drain nH pwm_S1 0 SW1\n";
        circuit << "S2 nH bridge_a pwm_S2 0 SW1\n";
        circuit << "S3 bridge_a nL pwm_S3 0 SW1\n";
        circuit << "S4 nL 0 pwm_S4 0 SW1\n";
        // NPC clamp: nH/nL pulled to mid_cap during the off-half so
        // bridge_a sees Vin/2 instead of floating.
        circuit << "DC1 mid_cap nH DIDEAL\n";
        circuit << "DC2 nL mid_cap DIDEAL\n\n";
    }
    else {
        // -----------------------------------------------------------------
        // Voltage-controlled-switch bridge — high-fidelity path.
        // -----------------------------------------------------------------
        // Split-capacitor midpoint — provides Vin/2 reference for the NPC
        // clamp diodes and the second terminal of the transformer primary
        // AC interface (point 'b').
        circuit << "C_split_hi vin_dc mid_cap 470u ic=" << (Vin/2.0) << "\n";
        circuit << "C_split_lo mid_cap 0 470u ic=" << (Vin/2.0) << "\n\n";

        auto pulse = [&](const std::string& name, double delay, double width) {
            circuit << "V" << name << " " << name << " 0 PULSE(0 5 "
                    << std::scientific << delay
                    << " " << cfg.pwmRise << " " << cfg.pwmFall << " "
                    << width << " " << period << std::fixed << ")\n";
        };

        pulse("pwm_S1", 0.0,        t_act);          // outer top (S1)
        pulse("pwm_S2", 0.0,        halfPeriod);     // inner top (S2): first half
        pulse("pwm_S3", halfPeriod, halfPeriod);     // inner bot (S3): second half
        pulse("pwm_S4", halfPeriod, t_act);          // outer bot (S4)
        circuit << "\n";

        // §8a.5 input-current probe — Vq1_sense is a 0-V ammeter inserted
        // between vin_dc and the outer-top switch S1's drain. i(Vq1_sense)
        // reports the clean high-side switch current, free of the snubber
        // RC spikes that contaminate i(Vdc). D1 (body diode) and the S1
        // snubber stay on the original vin_dc node so they don't feed the
        // ammeter — only the S1 channel current is measured. Used in
        // extractors as the converter's input current.
        circuit << "Vq1_sense vin_dc qa_drain 0\n";
        circuit << "S1 qa_drain nH pwm_S1 0 SW1\n";
        circuit << "D1 nH vin_dc DIDEAL\n";
        circuit << "S2 nH bridge_a pwm_S2 0 SW1\n";
        circuit << "D2 bridge_a nH DIDEAL\n";
        circuit << "S3 bridge_a nL pwm_S3 0 SW1\n";
        circuit << "D3 nL bridge_a DIDEAL\n";
        circuit << "S4 nL 0 pwm_S4 0 SW1\n";
        circuit << "D4 0 nL DIDEAL\n";
        circuit << "DC1 mid_cap nH DIDEAL\n";
        circuit << "DC2 nL mid_cap DIDEAL\n";
        circuit << "Rsnub_S1 vin_dc nH " << cfg.snubR
                << "\nCsnub_S1 vin_dc nH " << cfg.snubC << "\n";
        circuit << "Rsnub_S2 nH bridge_a " << cfg.snubR
                << "\nCsnub_S2 nH bridge_a " << cfg.snubC << "\n";
        circuit << "Rsnub_S3 bridge_a nL " << cfg.snubR
                << "\nCsnub_S3 bridge_a nL " << cfg.snubC << "\n";
        circuit << "Rsnub_S4 nL 0 " << cfg.snubR
                << "\nCsnub_S4 nL 0 " << cfg.snubC << "\n\n";
    }

    // Primary current sense (in series with the transformer leg) and
    // differential bridge-output probe vab = v(bridge_a) − v(mid_cap).
    circuit << "Vpri_sense bridge_a pri_lr 0\n";
    circuit << "Evab vab 0 bridge_a mid_cap 1\n\n";

    circuit << "L_series pri_lr trafo_pri " << std::scientific << Lr
            << " IC=" << (std::max(Io / n, 1.0) / 4.0) << "\n\n";

    // Transformer with multi-secondary K-matrix (primary between trafo_pri and mid_cap).
    constexpr double K_COUPLING = 0.9999;
    circuit << "L_pri trafo_pri mid_cap " << std::scientific << Lm << "\n";
    for (size_t i = 0; i < numOutputs; ++i) {
        double n_i = turnsRatios[i];
        if (n_i <= 0) n_i = 1.0;
        double Ls_i = Lm / (n_i * n_i);
        std::string si = std::to_string(i + 1);
        circuit << "L_sec_o" << si << " sec_a_o" << si << " sec_b_o" << si
                << " " << std::scientific << Ls_i << "\n";
    }
    int kCount = 0;
    for (size_t i = 0; i < numOutputs; ++i) {
        ++kCount;
        circuit << "K" << kCount << " L_pri L_sec_o" << (i + 1)
                << " " << K_COUPLING << "\n";
    }
    for (size_t i = 0; i < numOutputs; ++i) {
        for (size_t j = i + 1; j < numOutputs; ++j) {
            ++kCount;
            circuit << "K" << kCount
                    << " L_sec_o" << (i + 1)
                    << " L_sec_o" << (j + 1)
                    << " " << K_COUPLING << "\n";
        }
    }
    circuit << "\n";

    // ---- Differential winding-voltage probes (§8a.5 fix) ----
    //
    // The magnetic-view operating-point excitations require the voltage
    // ACROSS each transformer winding (the EMF that drives core flux),
    // not the converter-side terminal voltage. `vab` lumps together
    // L_series + L_pri (the primary winding alone), and `out_node_o*` is
    // the post-rectifier DC bus that shows ≈Vo instead of the bipolar
    // ±Vhb/n square that actually crosses the secondary winding.
    //
    // E-source probes here expose the actual L_pri / L_sec_o* terminal
    // voltages with the correct dot-convention sign. Polarity: the
    // primary winding is `L_pri trafo_pri mid_cap` (dot at trafo_pri),
    // so v_pri_w = v(trafo_pri) - v(mid_cap). Each secondary winding is
    // `L_sec_o<i> sec_a_o<i> sec_b_o<i>` with the dot at sec_a_o<i>.
    circuit << "Evpri_w vpri_w 0 trafo_pri mid_cap 1\n";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        circuit << "Evsec_w_o" << si << " vsec_w_o" << si
                << " 0 sec_a_o" << si << " sec_b_o" << si << " 1\n";
    }
    circuit << "\n";

    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        double Vo_i = pshbOp.get_output_voltages()[i];
        double Io_i = (i < pshbOp.get_output_currents().size())
                      ? pshbOp.get_output_currents()[i] : Io;
        if (Io_i <= 0) Io_i = 1.0;
        double Rload_i = std::max(Vo_i / Io_i, 1e-3);
        double Cout_i = 100e-6;

        circuit << "Vsec1_sense_o" << si << " sec_a_o" << si << " rec_a_o" << si << " 0\n";
        circuit << "Vsec2_sense_o" << si << " sec_b_o" << si << " rec_b_o" << si << " 0\n";

        if (rectType == BRectifierType::FULL_BRIDGE) {
            circuit << "D_r1_o" << si << " rec_a_o" << si << " out_rect_o" << si << " DIDEAL\n";
            circuit << "D_r2_o" << si << " rec_b_o" << si << " out_rect_o" << si << " DIDEAL\n";
            circuit << "D_r3_o" << si << " out_gnd_o" << si << " rec_a_o" << si << " DIDEAL\n";
            circuit << "D_r4_o" << si << " out_gnd_o" << si << " rec_b_o" << si << " DIDEAL\n";
            circuit << "L_out_o" << si << " out_rect_o" << si << " out_node_o" << si
                    << " " << std::scientific << Lo << " IC=" << Io_i << "\n";
        } else if (rectType == BRectifierType::CURRENT_DOUBLER) {
            double Io_half = Io_i / 2.0;
            circuit << "L_out1_o" << si << " rec_a_o" << si << " out_node_o" << si
                    << " " << std::scientific << Lo << " IC=" << Io_half << "\n";
            circuit << "L_out2_o" << si << " rec_b_o" << si << " out_node_o" << si
                    << " " << std::scientific << Lo << " IC=" << Io_half << "\n";
            circuit << "D_r1_o" << si << " out_gnd_o" << si << " rec_a_o" << si << " DIDEAL\n";
            circuit << "D_r2_o" << si << " out_gnd_o" << si << " rec_b_o" << si << " DIDEAL\n";
        } else {
            circuit << "* Center-tapped rectifier (output " << si << ")\n";
            circuit << "D_r1_o" << si << " rec_a_o" << si << " out_rect_o" << si << " DIDEAL\n";
            circuit << "D_r2_o" << si << " rec_b_o" << si << " out_rect_o" << si << " DIDEAL\n";
            circuit << "Vct_o" << si << " out_gnd_o" << si << " sec_ct_o" << si << " 0\n";
            circuit << "Rct_o" << si << " sec_ct_o" << si << " sec_b_o" << si << " 1u\n";
            circuit << "L_out_o" << si << " out_rect_o" << si << " out_node_o" << si
                    << " " << std::scientific << Lo << " IC=" << Io_i << "\n";
        }
        circuit << "C_out_o" << si << " out_node_o" << si << " out_gnd_o" << si
                << " " << std::scientific << Cout_i << " IC=" << Vo_i << "\n";
        circuit << "R_load_o" << si << " out_node_o" << si << " out_gnd_o" << si
                << " " << std::fixed << Rload_i << "\n";
        circuit << "Vout_sense_o" << si << " out_gnd_o" << si << " 0 0\n\n";
    }

    circuit << ".tran " << std::scientific << stepTime
            << " " << simTime << " " << startTime << "\n\n";

    circuit << ".save v(vin_dc) v(vab) v(bridge_a) v(mid_cap) v(nH) v(nL) v(vpri_w)"
            << " i(Vpri_sense) i(Vdc) i(Vq1_sense)";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        circuit << " v(out_node_o" << si << ")"
                << " v(vsec_w_o" << si << ")"
                << " i(Vsec1_sense_o" << si << ")"
                << " i(Vsec2_sense_o" << si << ")"
                << " i(Vout_sense_o" << si << ")";
    }
    circuit << "\n\n";

    circuit << ".options RELTOL=" << cfg.relTol
            << " ABSTOL=" << cfg.absTol
            << " VNTOL=" << cfg.vnTol
            << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
    circuit << ".options METHOD=" << cfg.method
            << " TRTOL=" << cfg.trTol << "\n";
    circuit << ".ic v(mid_cap)=" << (Vin/2.0);
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        double Vo_i = pshbOp.get_output_voltages()[i];
        circuit << " v(out_node_o" << si << ")=" << Vo_i;
    }
    circuit << "\n";

    // Nodeset hints for rectifier nodes — helps ngspice find correct DC
    // operating point so secondary diodes become forward biased.
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        double Vo_i = pshbOp.get_output_voltages()[i];
        double Vd_fwd = Vo_i + 0.7;  // ~Vo + diode drop
        circuit << ".nodeset v(rec_a_o" << si << ")=" << Vd_fwd
                << " v(rec_b_o" << si << ")=" << Vd_fwd
                << " v(out_rect_o" << si << ")=" << Vo_i
                << " v(sec_a_o" << si << ")=" << Vd_fwd
                << " v(sec_b_o" << si << ")=" << Vd_fwd << "\n";
    }

    circuit << "\n";

    circuit << ".end\n";

    return circuit.str();
}


std::vector<OperatingPoint> Pshb::simulate_and_extract_operating_points(
    const std::vector<double>& turnsRatios, double magnetizingInductance)
{
    std::vector<OperatingPoint> operatingPoints;
    NgspiceRunner runner;
    if (!runner.is_available()) {
        std::cerr << "[PSHB] ngspice not available — falling back to analytical "
                     "operating-point generation" << std::endl;
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    if (computedSeriesInductance <= 0) {
        process_design_requirements();
    }

    std::vector<double> inputVoltages;
    std::vector<std::string> inputVoltageNames;
    collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltageNames);

    auto& ops = get_operating_points();

    for (size_t vinIdx = 0; vinIdx < inputVoltages.size(); ++vinIdx) {
        for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
            double switchingFrequency = ops[opIdx].get_switching_frequency();
            if (switchingFrequency <= 0) {
                throw std::runtime_error("PSHB: operating point has invalid switching frequency");
            }

            std::string netlist = generate_ngspice_circuit(
                turnsRatios, magnetizingInductance, vinIdx, opIdx);

            SimulationConfig config;
            config.frequency = switchingFrequency;
            config.extractOnePeriod = true;
            config.numberOfPeriods = static_cast<size_t>(get_num_periods_to_extract());
            config.keepTempFiles = false;

            auto simResult = runner.run_simulation(netlist, config);
            if (!simResult.success) {
                std::cerr << "[PSHB] ngspice run failed: " << simResult.errorMessage
                          << " — falling back to analytical for this op." << std::endl;
                auto analytical = process_operating_point_for_input_voltage(
                    inputVoltages[vinIdx], ops[opIdx], turnsRatios, magnetizingInductance);
                analytical.set_name(inputVoltageNames[vinIdx] + " input ("
                    + std::to_string(static_cast<int>(inputVoltages[vinIdx])) + "V) [analytical]");
                operatingPoints.push_back(analytical);
                continue;
            }

            size_t numOutputs = std::min(turnsRatios.size(),
                                         ops[opIdx].get_output_voltages().size());
            if (numOutputs == 0) numOutputs = 1;

            NgspiceRunner::WaveformNameMapping waveformMapping;
            std::vector<std::string> windingNames;
            std::vector<bool> flipCurrentSign;

            // §8a.5 — use the differential winding-voltage probes
            // (vpri_w, vsec_w_o*) added in generate_ngspice_circuit.
            // These expose the actual EMF across each transformer
            // winding — NOT vab (which lumps L_series + L_pri) and
            // NOT out_node_o* (which is the post-rectifier DC rail).
            // Using those wrong probes here produced a primary trace
            // with a series-inductor ramp and secondary traces stuck
            // at ≈Vo, neither of which is a physical winding voltage.
            waveformMapping.push_back({{"voltage", "vpri_w"},
                                       {"current", "vpri_sense#branch"}});
            windingNames.push_back("Primary");
            flipCurrentSign.push_back(false);

            for (size_t outIdx = 0; outIdx < numOutputs; ++outIdx) {
                std::string si = std::to_string(outIdx + 1);
                waveformMapping.push_back({{"voltage", "vsec_w_o" + si},
                                           {"current", "vsec1_sense_o" + si + "#branch"}});
                windingNames.push_back("Secondary " + std::to_string(outIdx));
                flipCurrentSign.push_back(false);
            }

            OperatingPoint op = NgspiceRunner::extract_operating_point(
                simResult, waveformMapping, switchingFrequency, windingNames,
                ops[opIdx].get_ambient_temperature(), flipCurrentSign);

            op.set_name(inputVoltageNames[vinIdx] + " input ("
                + std::to_string(static_cast<int>(inputVoltages[vinIdx])) + "V)");
            operatingPoints.push_back(op);
        }
    }

    return operatingPoints;
}

std::vector<ConverterWaveforms> Pshb::simulate_and_extract_topology_waveforms(
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
                throw std::runtime_error("PSHB simulation failed: " + simResult.errorMessage);
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
            
            // §8a.5 converter-port stream — input voltage is the DC
            // source rail v(vin_dc) (the converter's external supply,
            // constant in steady state), NOT vab (bridge-midpoint AC
            // node, lumps the bridge + L_series + L_pri tank).
            //
            // Input current is the leg-top high-side switch current
            // i(Vq1_sense), measured by the in-series 0-V ammeter
            // inserted upstream of S1 (SW1 mode) or upstream of the
            // bridge midpoint in BEHAVIORAL_PULSE mode. It is positive
            // when current flows from vin_dc into the primary tank,
            // i.e. the converter is drawing from the bus.
            //
            // Why not i(Vdc)? The convergence-aid RC snubbers between
            // vin_dc and bridge nodes inject huge dV/dt-driven spikes
            // (~10^5 A peaks) at every switch transition that swamp
            // the real bus current — and in BEHAVIORAL_PULSE mode
            // i(Vdc) is identically zero (Vdc is decoupled from the
            // independent Vbridge_a / Vmid_cap stiff sources).
            // Why not i(Vpri_sense)? Vpri_sense sits inside the
            // primary tank (between bridge_a and L_series), so it
            // measures the bipolar tank current that averages to zero
            // — not the converter's actual input-port draw.
            //
            // Numerical-artifact clamp: ngspice's SW model transitions
            // instantaneously between Roff and Ron, so when S1 closes
            // the di/dt is unbounded and the resampled trace can hold
            // transient spikes ~10^5 A that are pure simulation noise
            // (NOT measured current). Physically, |i(S1)| <= |i(L_pri)|
            // because S1 only conducts the primary tank current while
            // ON. So we clamp the raw S1 trace to ±2·max|i(Vpri_sense)|
            // (2× headroom for the switching-instant overshoot we want
            // to keep visible). This is a numerical guard against the
            // ngspice idealised-switch di/dt artifact, not a physical
            // bound on the current.
            wf.set_input_voltage(getWaveform("vin_dc"));
            Waveform iQ1 = getWaveform("vq1_sense#branch");
            Waveform iPri = getWaveform("vpri_sense#branch");
            std::vector<double>& iQ1Data = iQ1.get_mutable_data();
            const std::vector<double>& iPriData = iPri.get_data();
            double iPriMax = 0.0;
            for (double v : iPriData) iPriMax = std::max(iPriMax, std::abs(v));
            const double clampLimit = 2.0 * iPriMax;
            if (clampLimit > 0.0) {
                for (auto& v : iQ1Data) {
                    if (v >  clampLimit) v =  clampLimit;
                    if (v < -clampLimit) v = -clampLimit;
                }
            }
            wf.set_input_current(iQ1);

            size_t numOutputs = std::min(turnsRatios.size(), opPoint.get_output_voltages().size());
            if (numOutputs == 0) numOutputs = 1;
            for (size_t outIdx = 0; outIdx < numOutputs; ++outIdx) {
                std::string si = std::to_string(outIdx + 1);
                wf.get_mutable_output_voltages().push_back(getWaveform("out_node_o" + si));
                // Reconstruct Iout(t) = Vout(t)/Rload (DC by design).
                // Vout_sense in this netlist sits in the cap-return path
                // so its branch current averages to zero in steady state.
                // Mirrors Buck.cpp / IsolatedBuck.cpp / AHB / PSFB §5.1.
                const double Vo_i = opPoint.get_output_voltages()[outIdx];
                const double Io_i = (outIdx < opPoint.get_output_currents().size())
                                    ? opPoint.get_output_currents()[outIdx]
                                    : opPoint.get_output_currents()[0];
                const double Rload_i = std::max(Vo_i / Io_i, 1e-3);
                Waveform ioutWf = getWaveform("out_node_o" + si);
                auto& ioutData = ioutWf.get_mutable_data();
                for (auto& v : ioutData) v = v / Rload_i;
                wf.get_mutable_output_currents().push_back(ioutWf);
            }

            results.push_back(wf);
        }
    }
    
    // Restore original value
    set_num_periods_to_extract(originalNumPeriodsToExtract);
    
    return results;
}


// =========================================================================
// AdvancedPshb
// =========================================================================
Inputs AdvancedPshb::process() {
    auto designRequirements = process_design_requirements();

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

DesignRequirements AdvancedPshb::process_design_requirements() {
    // Issue M1: chain to parent PDR then override turns ratios and
    // magnetizing inductance with user-provided desired* values.
    auto designRequirements = Pshb::process_design_requirements();

    designRequirements.get_mutable_turns_ratios().clear();
    for (auto n : desiredTurnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(n);
        designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }

    DimensionWithTolerance LmTol;
    LmTol.set_minimum(desiredMagnetizingInductance);
    designRequirements.set_magnetizing_inductance(LmTol);

    return designRequirements;
}


// =========================================================================
// Extra components inputs (output inductor Lo and series/ZVS inductor Lr)
// =========================================================================
std::vector<std::variant<Inputs, CAS::Inputs>> Pshb::get_extra_components_inputs(
    ExtraComponentsMode mode,
    std::optional<Magnetic> magnetic)
{
    if (mode == ExtraComponentsMode::REAL && !magnetic.has_value())
        throw std::invalid_argument("Pshb::get_extra_components_inputs: mode REAL requires a designed magnetic");

    if (computedOutputInductance <= 0 || extraLoVoltageWaveforms.empty())
        throw std::runtime_error("Pshb::get_extra_components_inputs: call process() first");

    std::vector<std::variant<Inputs, CAS::Inputs>> result;
    size_t nOps = extraLoVoltageWaveforms.size();
    double fsw = get_operating_points()[0].get_switching_frequency();

    // ---- Output inductor Lo ----
    {
        Inputs masInputs;
        DesignRequirements dr;

        DimensionWithTolerance loInductance;
        loInductance.set_nominal(computedOutputInductance);
        dr.set_magnetizing_inductance(loInductance);
        dr.set_name("outputInductor");
        dr.set_topology(MAS::Topology::PHASE_SHIFTED_HALF_BRIDGE_CONVERTER);
        dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
        dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
        masInputs.set_design_requirements(dr);

        std::vector<OperatingPoint> masOps;
        for (size_t i = 0; i < nOps; ++i) {
            OperatingPoint op;
            auto excitation = complete_excitation(
                extraLoCurrentWaveforms[i], extraLoVoltageWaveforms[i], fsw, "Primary");
            op.get_mutable_excitations_per_winding().push_back(excitation);
            masOps.push_back(op);
        }
        masInputs.set_operating_points(masOps);
        result.emplace_back(std::move(masInputs));
    }

    // ---- Second output inductor (current-doubler only) ----
    if (!extraLo2VoltageWaveforms.empty()) {
        Inputs masInputs;
        DesignRequirements dr;

        DimensionWithTolerance lo2Inductance;
        lo2Inductance.set_nominal(computedOutputInductance);
        dr.set_magnetizing_inductance(lo2Inductance);
        dr.set_name("outputInductor2");
        dr.set_topology(MAS::Topology::PHASE_SHIFTED_HALF_BRIDGE_CONVERTER);
        dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
        dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
        masInputs.set_design_requirements(dr);

        std::vector<OperatingPoint> masOps;
        for (size_t i = 0; i < extraLo2VoltageWaveforms.size(); ++i) {
            OperatingPoint op;
            op.get_mutable_excitations_per_winding().push_back(
                complete_excitation(extraLo2CurrentWaveforms[i],
                                    extraLo2VoltageWaveforms[i], fsw, "Primary"));
            masOps.push_back(op);
        }
        masInputs.set_operating_points(masOps);
        result.emplace_back(std::move(masInputs));
    }

    // ---- Series (ZVS) inductor Lr ----
    double Lr = computedSeriesInductance;
    double Llk = 0.0;

    if (mode == ExtraComponentsMode::REAL) {
        auto leakageOutput = LeakageInductance().calculate_leakage_inductance_all_windings(
            magnetic.value(), fsw);
        auto perWinding = leakageOutput.get_leakage_inductance_per_winding();
        if (!perWinding.empty())
            Llk = resolve_dimensional_values(perWinding[0]);
    }

    double Lr_external = (mode == ExtraComponentsMode::IDEAL)
        ? Lr
        : ((Llk < Lr) ? (Lr - Llk) : 0.0);

    if (Lr_external > 0.0 && !extraLrVoltageWaveforms.empty()) {
        Inputs masInputs;
        DesignRequirements dr;

        DimensionWithTolerance lrInductance;
        lrInductance.set_nominal(Lr_external);
        dr.set_magnetizing_inductance(lrInductance);
        dr.set_name("seriesInductor");
        dr.set_topology(MAS::Topology::PHASE_SHIFTED_HALF_BRIDGE_CONVERTER);
        dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
        dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
        masInputs.set_design_requirements(dr);

        std::vector<OperatingPoint> masOps;
        for (size_t i = 0; i < nOps; ++i) {
            OperatingPoint op;
            auto excitation = complete_excitation(
                extraLrCurrentWaveforms[i], extraLrVoltageWaveforms[i], fsw, "Primary");
            op.get_mutable_excitations_per_winding().push_back(excitation);
            masOps.push_back(op);
        }
        masInputs.set_operating_points(masOps);
        result.emplace_back(std::move(masInputs));
    }

    return result;
}

} // namespace OpenMagnetics
