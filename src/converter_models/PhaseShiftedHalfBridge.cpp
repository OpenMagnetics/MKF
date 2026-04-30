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
    double Vd = 0.6;
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

    designRequirements.set_topology(Topologies::PHASE_SHIFTED_HALF_BRIDGE_CONVERTER);
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
    const PshbOperatingPoint& pshbOpPoint,
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    OperatingPoint operatingPoint;

    double Fs = pshbOpPoint.get_switching_frequency();
    double Vin = inputVoltage;
    double Vhb = Vin * BRIDGE_VOLTAGE_FACTOR;
    double Vo  = pshbOpPoint.get_output_voltages()[0];
    double Io  = pshbOpPoint.get_output_currents()[0];
    double n   = turnsRatios[0];
    double Lm  = magnetizingInductance;
    double Lr  = computedSeriesInductance > 0 ? computedSeriesInductance : 1e-6;

    BRectifierType rectType = get_rectifier_type().value_or(BRectifierType::CENTER_TAPPED);

    double phi_deg = pshbOpPoint.get_phase_shift();
    double D_cmd = (phi_deg > 1e-6) ? compute_effective_duty_cycle(phi_deg)
                                    : computedEffectiveDutyCycle;
    if (D_cmd <= 0) D_cmd = 0.75;

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
            // PSHB single-leg HB has NO physical freewheel state — the two
            // primary switches are complementary at 50 % duty, so Vab is
            // always ±Vbus, never zero (apart from sub-µs dead time). The
            // 3-level shape we generate analytically is an approximation
            // that the SPICE single-leg netlist does not produce; this is
            // documented in the topology-disambiguation note in the
            // header docstring. Hold Ipri at end-of-active value here as
            // the least-bad approximation of the "freewheel" interval.
            Ipri_full[k] = ILo_max / n + Im_peak;
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
    double stepTime = period / 200;

    BRectifierType rectType = get_rectifier_type().value_or(BRectifierType::CENTER_TAPPED);

    size_t numOutputs = std::min(turnsRatios.size(), pshbOp.get_output_voltages().size());
    if (numOutputs == 0) numOutputs = 1;

    std::ostringstream circuit;

    circuit << "* Phase-Shifted Half Bridge (PSHB) — single-leg + split-cap\n";
    circuit << "* Vin=" << Vin << "V (primary swings ±Vin/2), Vo=" << Vo << "V, Fs="
            << (Fs/1e3) << "kHz, Deff=" << Deff << ", outputs=" << numOutputs << "\n";
    circuit << "* n=" << n << ", Lr=" << (Lr*1e6) << "uH, Lm=" << (Lm*1e6)
            << "uH, Lo=" << (Lo*1e6) << "uH, rect=" << static_cast<int>(rectType) << "\n\n";

    circuit << ".model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg\n";
    circuit << ".model DIDEAL D(IS=1e-12 RS=0.05)\n\n";

    circuit << "Vdc vin_dc 0 " << Vin << "\n\n";

    // Split capacitor divider — clamps mid_cap to Vin/2 so the primary
    // (between mid_sw and mid_cap) sees ±Vin/2.
    circuit << "C_split_hi vin_dc mid_cap 470u ic=" << (Vin/2.0) << "\n";
    circuit << "C_split_lo mid_cap 0 470u ic=" << (Vin/2.0) << "\n";
    circuit << "Rbal_hi vin_dc mid_cap 100k\n";
    circuit << "Rbal_lo mid_cap 0 100k\n\n";

    // Single half-bridge leg (QA-QB): 50% complementary duty.
    circuit << "Vpwm_A pwm_A 0 PULSE(0 5 0 10n 10n "
            << std::scientific << tOn << " " << period << std::fixed << ")\n";
    circuit << "Vpwm_B pwm_B 0 PULSE(0 5 "
            << std::scientific << halfPeriod << std::fixed
            << " 10n 10n "
            << std::scientific << tOn << " " << period << std::fixed << ")\n";
    circuit << "SA vin_dc mid_sw pwm_A 0 SW1\n";
    circuit << "DA 0 mid_sw DIDEAL\n";
    circuit << "SB mid_sw 0 pwm_B 0 SW1\n";
    circuit << "DB mid_sw vin_dc DIDEAL\n";
    circuit << "Rsnub_QA vin_dc mid_sw 1k\nCsnub_QA vin_dc mid_sw 1n\n";
    circuit << "Rsnub_QB mid_sw 0 1k\nCsnub_QB mid_sw 0 1n\n\n";

    // Primary current sense + primary winding voltage probe.
    // The *transformer* primary sits between mid_sw and mid_cap via Lr,
    // and swings ±Vin/2.
    circuit << "Vpri_sense mid_sw pri_lr 0\n";
    circuit << "Evab vab 0 mid_sw mid_cap 1\n\n";

    circuit << "L_series pri_lr trafo_pri " << std::scientific << Lr << "\n\n";

    // Transformer with multi-secondary K-matrix.
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
                    << " " << std::scientific << Lo << "\n";
        } else if (rectType == BRectifierType::CURRENT_DOUBLER) {
            circuit << "L_out1_o" << si << " rec_a_o" << si << " out_node_o" << si
                    << " " << std::scientific << Lo << "\n";
            circuit << "L_out2_o" << si << " rec_b_o" << si << " out_node_o" << si
                    << " " << std::scientific << Lo << "\n";
            circuit << "D_r1_o" << si << " out_gnd_o" << si << " rec_a_o" << si << " DIDEAL\n";
            circuit << "D_r2_o" << si << " out_gnd_o" << si << " rec_b_o" << si << " DIDEAL\n";
        } else {
            circuit << "* Center-tapped rectifier (output " << si << ")\n";
            circuit << "D_r1_o" << si << " rec_a_o" << si << " out_rect_o" << si << " DIDEAL\n";
            circuit << "D_r2_o" << si << " rec_b_o" << si << " out_rect_o" << si << " DIDEAL\n";
            circuit << "Vct_o" << si << " out_gnd_o" << si << " sec_ct_o" << si << " 0\n";
            circuit << "Rct_o" << si << " sec_ct_o" << si << " sec_b_o" << si << " 1u\n";
            circuit << "L_out_o" << si << " out_rect_o" << si << " out_node_o" << si
                    << " " << std::scientific << Lo << "\n";
        }
        circuit << "C_out_o" << si << " out_node_o" << si << " out_gnd_o" << si
                << " " << std::scientific << Cout_i << "\n";
        circuit << "R_load_o" << si << " out_node_o" << si << " out_gnd_o" << si
                << " " << std::fixed << Rload_i << "\n";
        circuit << "Vout_sense_o" << si << " out_gnd_o" << si << " 0 0\n\n";
    }

    circuit << ".tran " << std::scientific << stepTime
            << " " << simTime << " " << startTime << "\n\n";

    circuit << ".save v(vab) v(mid_sw) v(mid_cap) i(Vpri_sense) i(Vdc)";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        circuit << " v(out_node_o" << si << ")"
                << " i(Vsec1_sense_o" << si << ")"
                << " i(Vsec2_sense_o" << si << ")"
                << " i(Vout_sense_o" << si << ")";
    }
    circuit << "\n\n";

    circuit << ".options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500\n";
    circuit << ".options METHOD=GEAR TRTOL=7\n";
    circuit << ".ic v(mid_cap)=" << (Vin/2.0);
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        double Vo_i = pshbOp.get_output_voltages()[i];
        circuit << " v(out_node_o" << si << ")=" << Vo_i;
    }
    circuit << "\n\n";

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

            waveformMapping.push_back({{"voltage", "vab"},
                                       {"current", "vpri_sense#branch"}});
            windingNames.push_back("Primary");
            flipCurrentSign.push_back(false);

            for (size_t outIdx = 0; outIdx < numOutputs; ++outIdx) {
                std::string si = std::to_string(outIdx + 1);
                waveformMapping.push_back({{"voltage", "out_node_o" + si},
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
            
            wf.set_input_voltage(getWaveform("vab"));
            wf.set_input_current(getWaveform("vpri_sense#branch"));

            size_t numOutputs = std::min(turnsRatios.size(), opPoint.get_output_voltages().size());
            if (numOutputs == 0) numOutputs = 1;
            for (size_t outIdx = 0; outIdx < numOutputs; ++outIdx) {
                std::string si = std::to_string(outIdx + 1);
                wf.get_mutable_output_voltages().push_back(getWaveform("out_node_o" + si));
                wf.get_mutable_output_currents().push_back(getWaveform("vout_sense_o" + si + "#branch"));
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
        dr.set_topology(Topologies::PHASE_SHIFTED_HALF_BRIDGE_CONVERTER);
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
        dr.set_topology(Topologies::PHASE_SHIFTED_HALF_BRIDGE_CONVERTER);
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
        dr.set_topology(Topologies::PHASE_SHIFTED_HALF_BRIDGE_CONVERTER);
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
