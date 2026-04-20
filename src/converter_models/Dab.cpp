#include "converter_models/Dab.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdexcept>
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
// =========================================================================
// Static helpers: General piecewise-linear voltage waveforms
// =========================================================================
//
// All modulation types (SPS/EPS/DPS/TPS) share the same bridge structure.
// Inner phase shifts D1, D2 (radians) control the bridge duty cycles.
// Outer phase shift phi controls power transfer.
//
// Primary bridge output Vab(theta), theta in [0, 2pi):
//   [0,   D1)     : 0
//   [D1,  pi)     : +V1
//   [pi,  pi+D1)  : 0
//   [pi+D1, 2pi)  : -V1
//
// Secondary bridge output Vcd(theta), shifted by phi, inner shift D2:
//   theta_sec = (theta - phi) mod 2pi
//   [0,   D2)     : 0
//   [D2,  pi)     : +V2
//   [pi,  pi+D2)  : 0
//   [pi+D2, 2pi)  : -V2
//
// Inductor voltage (referred to primary): vL = Vab - N*Vcd
// Inductor current slope: diL/dtheta = vL / (L * 2*pi*Fs)
//
// Initial condition from half-wave antisymmetry (iL(theta+pi) = -iL(theta)):
//   delta = integral_0^pi  vL/L d(theta) / (2*pi*Fs)
//   iL(0) = -delta/2
//
// Power: P = (1/(2*pi)) * integral_0^{2pi} Vab(theta) * iL(theta) d(theta)
// =========================================================================

static double Vab_at(double theta, double V1, double D1) {
    theta = std::fmod(theta, 2.0 * M_PI);
    if (theta < 0) theta += 2.0 * M_PI;
    if (theta < D1)            return 0.0;
    if (theta < M_PI)          return  V1;
    if (theta < M_PI + D1)     return 0.0;
    return -V1;
}

static double Vcd_at(double theta, double V2, double D2, double phi) {
    double ts = std::fmod(theta - phi, 2.0 * M_PI);
    if (ts < 0) ts += 2.0 * M_PI;
    if (ts < D2)            return 0.0;
    if (ts < M_PI)          return  V2;
    if (ts < M_PI + D2)     return 0.0;
    return -V2;
}

// =========================================================================
// Closed-form sub-interval propagator (replaces the previous Euler scheme)
// =========================================================================
//
// All four DAB modulation modes (SPS / EPS / DPS / TPS) share the same
// underlying structure: within one switching period, both bridge voltages
// Vab(θ) and Vcd(θ) are piecewise-constant. The breakpoints come only from
// the known constants D1, D2 and φ. Inside any sub-interval [θ_k, θ_{k+1}]:
//
//   Vab(θ) = const ⇒ dIm/dθ = Vab/(Lm·ω) = const ⇒ Im(θ) is exactly LINEAR
//   vL = Vab − N·Vcd = const ⇒ diL/dθ = vL/(L·ω) = const ⇒ iL(θ) is exactly LINEAR
//
// where ω = 2π·Fs (so dt = dθ/ω).
//
// Strategy:
//   1. Build the sorted list of unique sub-interval boundary angles in [0, 2π).
//   2. The number of sub-intervals is at most 8 (TPS): {0, D1, π, π+D1,
//      φ, φ+D2, φ+π, φ+π+D2}, all reduced mod 2π and de-duplicated.
//   3. iL(0) and Im(0) are fixed by half-wave antisymmetry x(π) = −x(0):
//        iL(0) = − ½ · ∫₀^π (vL/L) dt   = − ½ · Σ_{k: θ_k<π} (vL_k/L)·Δt_k
//        Im(0) = − ½ · ∫₀^π (Vab/Lm) dt = − ½ · Σ_{k: θ_k<π} (Vab_k/Lm)·Δt_k
//   4. Forward-propagate iL and Im exactly through each sub-interval. The
//      result at each θ_k is exact, and within a sub-interval the values
//      are linear interpolations of the endpoints (also exact).
//
// The same routine is used by compute_power_general (which integrates Vab·iL
// over the period exactly) and by process_operating_point_for_input_voltage
// (which samples the resulting waveform onto the 256-point output grid).
// =========================================================================

struct DabSubInterval {
    double theta_start;   // radians, [0, 2π]
    double theta_end;
    double Vab;           // primary bridge voltage in this sub-interval
    double Vcd;           // secondary bridge voltage in this sub-interval (referred to its own side)
};

// Build the sorted, deduplicated list of sub-interval boundary angles in
// [0, 2π], including both 0 and 2π. Returns at most 9 entries (8 segments).
static std::vector<double> dab_boundary_angles(double D1, double D2, double phi) {
    auto wrap = [](double a) {
        a = std::fmod(a, 2.0 * M_PI);
        if (a < 0) a += 2.0 * M_PI;
        return a;
    };
    std::vector<double> raw = {
        0.0,
        wrap(D1),
        wrap(M_PI),
        wrap(M_PI + D1),
        wrap(phi),
        wrap(phi + D2),
        wrap(phi + M_PI),
        wrap(phi + M_PI + D2),
    };
    std::sort(raw.begin(), raw.end());
    std::vector<double> out;
    out.reserve(raw.size() + 1);
    constexpr double EPS = 1e-9;
    for (double a : raw) {
        if (out.empty() || a - out.back() > EPS) out.push_back(a);
    }
    if (out.empty() || out.back() < 2.0 * M_PI - EPS) out.push_back(2.0 * M_PI);
    return out;
}

// Build the sub-interval list for the full period. Each entry has the
// constant Vab and Vcd values that apply over [theta_start, theta_end].
static std::vector<DabSubInterval> dab_build_subintervals(
    double V1, double V2, double D1, double D2, double phi)
{
    auto angles = dab_boundary_angles(D1, D2, phi);
    std::vector<DabSubInterval> segs;
    segs.reserve(angles.size());
    for (size_t k = 0; k + 1 < angles.size(); ++k) {
        double a0 = angles[k];
        double a1 = angles[k + 1];
        double mid = 0.5 * (a0 + a1);
        DabSubInterval seg;
        seg.theta_start = a0;
        seg.theta_end   = a1;
        seg.Vab = Vab_at(mid, V1, D1);
        seg.Vcd = Vcd_at(mid, V2, D2, phi);
        segs.push_back(seg);
    }
    return segs;
}

// Closed-form initial conditions for iL and Im at θ=0, enforced by the
// half-wave antisymmetry constraint x(π) = −x(0). Uses the analytical
// piecewise-linear integral.
static void dab_initial_conditions(const std::vector<DabSubInterval>& segs,
                                   double N, double L, double Lm,
                                   double Fs, double& iL0, double& Im0)
{
    // ∫₀^π (vL/L) dt and ∫₀^π (Vab/Lm) dt
    double inv_iL = 1.0 / (L * 2.0 * M_PI * Fs);
    double inv_Im = 1.0 / (Lm * 2.0 * M_PI * Fs);
    double sum_iL = 0.0;
    double sum_Im = 0.0;
    for (const auto& s : segs) {
        if (s.theta_start >= M_PI) break;
        double end = std::min(s.theta_end, M_PI);
        double dtheta = end - s.theta_start;
        if (dtheta <= 0) continue;
        double vL = s.Vab - N * s.Vcd;
        sum_iL += vL * inv_iL * dtheta;
        sum_Im += s.Vab * inv_Im * dtheta;
    }
    iL0 = -0.5 * sum_iL;
    Im0 = -0.5 * sum_Im;
}

// Forward-propagate iL and Im exactly through every sub-interval. Returns
// arrays of (theta, iL, Im) at every sub-interval boundary, length =
// segs.size() + 1. By construction the last entry equals -first entry to
// within floating-point round-off (antisymmetry holds exactly when N=integer
// segments cover [0, 2π], which they always do here).
static void dab_propagate_period(
    const std::vector<DabSubInterval>& segs,
    double N, double L, double Lm, double Fs,
    double iL0, double Im0,
    std::vector<double>& boundary_theta,
    std::vector<double>& boundary_iL,
    std::vector<double>& boundary_Im)
{
    double inv_iL = 1.0 / (L * 2.0 * M_PI * Fs);
    double inv_Im = 1.0 / (Lm * 2.0 * M_PI * Fs);
    boundary_theta.clear();
    boundary_iL.clear();
    boundary_Im.clear();
    boundary_theta.reserve(segs.size() + 1);
    boundary_iL.reserve(segs.size() + 1);
    boundary_Im.reserve(segs.size() + 1);

    double iL = iL0;
    double Im = Im0;
    boundary_theta.push_back(0.0);
    boundary_iL.push_back(iL);
    boundary_Im.push_back(Im);
    for (const auto& s : segs) {
        double dtheta = s.theta_end - s.theta_start;
        double vL = s.Vab - N * s.Vcd;
        iL += vL * inv_iL * dtheta;
        Im += s.Vab * inv_Im * dtheta;
        boundary_theta.push_back(s.theta_end);
        boundary_iL.push_back(iL);
        boundary_Im.push_back(Im);
    }
}

// Sample iL or Im at angle θ ∈ [0, 2π] given the boundary arrays. Each
// sub-interval is exactly linear in θ (because Vab/Vcd are constant), so
// linear interpolation between boundary samples is mathematically exact.
static double dab_sample(double theta,
                         const std::vector<double>& bnd_theta,
                         const std::vector<double>& bnd_y)
{
    if (bnd_theta.empty()) return 0.0;
    if (theta <= bnd_theta.front()) return bnd_y.front();
    if (theta >= bnd_theta.back())  return bnd_y.back();
    // Binary search for the right segment
    auto it = std::upper_bound(bnd_theta.begin(), bnd_theta.end(), theta);
    size_t hi = static_cast<size_t>(it - bnd_theta.begin());
    size_t lo = hi - 1;
    double t = (theta - bnd_theta[lo]) / (bnd_theta[hi] - bnd_theta[lo]);
    return bnd_y[lo] + t * (bnd_y[hi] - bnd_y[lo]);
}

double Dab::compute_power_general(double V1, double V2, double N,
                                   double phi, double D1, double D2,
                                   double Fs, double L) {
    // Closed-form average power over one period using piecewise-constant Vab
    // and piecewise-linear iL. Within each sub-interval [θ_k, θ_{k+1}]:
    //   contribution = Vab_k · ∫_{θ_k}^{θ_{k+1}} iL(θ) dθ / (2π)
    // and ∫ iL(θ) dθ = ½ (iL_start + iL_end) · Δθ  (because iL is linear).
    auto segs = dab_build_subintervals(V1, V2, D1, D2, phi);
    // We need Lm too in dab_initial_conditions (for Im), but the power
    // computation only depends on iL — pass a dummy Lm > 0.
    constexpr double Lm_dummy = 1.0;
    double iL0, Im0;
    dab_initial_conditions(segs, N, L, Lm_dummy, Fs, iL0, Im0);
    std::vector<double> bnd_theta, bnd_iL, bnd_Im;
    dab_propagate_period(segs, N, L, Lm_dummy, Fs, iL0, Im0,
                         bnd_theta, bnd_iL, bnd_Im);

    double power = 0.0;
    for (size_t k = 0; k < segs.size(); ++k) {
        double dtheta = segs[k].theta_end - segs[k].theta_start;
        double iL_avg = 0.5 * (bnd_iL[k] + bnd_iL[k + 1]);
        power += segs[k].Vab * iL_avg * dtheta;
    }
    return power / (2.0 * M_PI);
}

// Binary-search for phi given desired power (for EPS/DPS/TPS).
// Returns the achieved phi; if the requested power exceeds what the topology
// can deliver at this L/Fs/voltage operating point, saturates at π/2 and
// emits a warning so the caller can react.
double Dab::compute_phase_shift_general(double V1, double V2, double N,
                                         double D1, double D2,
                                         double Fs, double L, double P) {
    double P_max = compute_power_general(V1, V2, N, M_PI / 2.0, D1, D2, Fs, L);
    if (P >= P_max) {
        std::cerr << "[DAB] requested power " << P << "W exceeds achievable max "
                  << P_max << "W at this operating point — saturating phi at π/2"
                  << std::endl;
        return M_PI / 2.0;
    }
    double phi_lo = 0.0;
    double phi_hi = M_PI / 2.0;
    for (int iter = 0; iter < 60; ++iter) {
        double phi_mid = (phi_lo + phi_hi) / 2.0;
        double P_mid = compute_power_general(V1, V2, N, phi_mid, D1, D2, Fs, L);
        if (P_mid < P) phi_lo = phi_mid;
        else           phi_hi = phi_mid;
        if ((phi_hi - phi_lo) < 1e-8) break;
    }
    return (phi_lo + phi_hi) / 2.0;
}


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
        // Phase shift bound: |phi| <= 90 deg (SPS theoretical maximum).
        double phi_deg = op.get_inner_phase_shift3().value_or(0.0);
        if (std::abs(phi_deg) > 90.0) {
            if (assertErrors) throw std::runtime_error("DAB: phase shift out of range (|phi| > 90 deg)");
            ok = false;
        }

        // Inner phase shifts: 0 <= D < 90 deg (90 collapses the half-cycle to 0 width).
        double D1_deg = op.get_inner_phase_shift1().value_or(0.0);
        double D2_deg = op.get_inner_phase_shift2().value_or(0.0);
        if (D1_deg < 0.0 || D1_deg >= 90.0) {
            if (assertErrors) throw std::runtime_error("DAB: innerPhaseShift1 must be in [0, 90) degrees");
            ok = false;
        }
        if (D2_deg < 0.0 || D2_deg >= 90.0) {
            if (assertErrors) throw std::runtime_error("DAB: innerPhaseShift2 must be in [0, 90) degrees");
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
    double mainOutputPower = 0.0;
    for (size_t i = 0; i < ops[0].get_output_voltages().size() && i < ops[0].get_output_currents().size(); ++i) {
        mainOutputPower += ops[0].get_output_voltages()[i] * ops[0].get_output_currents()[i];
    }
    if (mainOutputPower <= 0) {
        mainOutputPower = mainOutputVoltage * mainOutputCurrent;
    }
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
    double phi_deg = ops[0].get_inner_phase_shift3().value_or(0.0);
    double phi_rad = phi_deg * M_PI / 180.0;

    // For EPS/DPS/TPS: get inner phase shifts from operating point
    double D1_rad = Dab::get_D1_rad(ops[0]);
    double D2_rad = Dab::get_D2_rad(ops[0]);

    // 3. Series inductance
    double L;
    if (get_series_inductance().has_value() && get_series_inductance().value() > 0) {
        L = get_series_inductance().value();
        // If phase shift is zero or very small, compute it from power
        if (std::abs(phi_rad) < 1e-6 && mainOutputPower > 0) {
            auto modTypeOpt = ops[0].get_modulation_type();
            bool isSPS = !modTypeOpt.has_value() || modTypeOpt.value() == ModulationType::SPS;
            if (isSPS) {
                phi_rad = compute_phase_shift(Vin_nom, mainOutputVoltage, N, Fs, L, mainOutputPower);
            } else {
                phi_rad = compute_phase_shift_general(Vin_nom, mainOutputVoltage, N,
                                                       D1_rad, D2_rad, Fs, L, mainOutputPower);
            }
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
    //    Lm should be large enough that magnetizing current ripple is manageable.
    //    Im_peak = V1 / (4 * Fs * Lm)
    //    Target: Im_peak < 30% of primary load current (conservative enough for ZVS,
    //    practical enough for standard ferrite cores at typical power densities).
    //
    //    I_load_primary (approximate) = P / V1
    //    So Lm > V1 / (4 * Fs * 0.3 * P / V1) = V1^2 / (1.2 * Fs * P)
    //
    //    We use a factor of 10x the series inductance as a minimum.
    double I_load_pri = mainOutputPower / Vin_nom;
    double Im_target = 0.3 * I_load_pri; // 30% of load current
    double Lm_from_current = (Im_target > 0) ? Vin_nom / (4.0 * Fs * Im_target) : 10.0 * L;
    double Lm_from_ratio = 10.0 * L; // At least 10x series inductance
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
    inductanceWithTolerance.set_minimum(roundFloat(Lm, 10));
    designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

    // If using leakage inductance as series inductor, request leakage = L
    if (get_use_leakage_inductance().value_or(false)) {
        std::vector<DimensionWithTolerance> leakageReqs;
        DimensionWithTolerance lrTol;
        lrTol.set_nominal(roundFloat(L, 10));
        leakageReqs.push_back(lrTol);
        designRequirements.set_leakage_inductance(leakageReqs);
    }

    designRequirements.set_topology(Topologies::DUAL_ACTIVE_BRIDGE_CONVERTER);
    designRequirements.set_isolation_sides(
        Topology::create_isolation_sides(ops[0].get_output_currents().size(), false));

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

    // Sort and label so callers can disambiguate by Vin level.
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
            Vin, ops[0], turnsRatios, magnetizingInductance);
        op.set_name(name + " input (" + std::to_string(static_cast<int>(Vin)) + "V)");
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
    if (Fs <= 0)
        throw std::runtime_error("DAB: switching frequency must be positive");
    if (dabOpPoint.get_output_voltages().empty())
        throw std::runtime_error("DAB: operating point has no output voltages");
    if (turnsRatios.empty())
        throw std::runtime_error("DAB: turnsRatios is empty");

    double V1 = inputVoltage;
    double N_main = turnsRatios[0];  // Primary-to-secondary turns ratio for output 0
    double V2_main = dabOpPoint.get_output_voltages()[0];
    double Lm = magnetizingInductance;
    double L = computedSeriesInductance;
    if (Lm <= 0) throw std::runtime_error("DAB: magnetizing inductance must be positive");
    if (L <= 0) throw std::runtime_error("DAB: series inductance must be positive (run process_design_requirements first)");

    // Phase shift: signed (degrees) → radians. Negative phi = reverse power flow.
    double phi_deg = dabOpPoint.get_inner_phase_shift3().value_or(0.0);
    double phi_rad = (std::abs(phi_deg) > 1e-6)
                   ? phi_deg * M_PI / 180.0
                   : computedPhaseShift;

    double period = 1.0 / Fs;
    double Thalf = period / 2.0;

    double D1_rad = Dab::get_D1_rad(dabOpPoint);
    double D2_rad = Dab::get_D2_rad(dabOpPoint);

    // ---- Detect modulation type for diagnostics ----
    // SPS: D1 = D2 = 0; EPS: exactly one of {D1, D2} > 0; DPS: D1 = D2 > 0;
    // TPS: D1 ≠ D2 and both > 0.
    constexpr double D_TOL = 1e-6;
    int mod_type = 0;  // SPS
    if (D1_rad > D_TOL || D2_rad > D_TOL) {
        if (D1_rad > D_TOL && D2_rad < D_TOL)      mod_type = 1;  // EPS (primary)
        else if (D1_rad < D_TOL && D2_rad > D_TOL) mod_type = 1;  // EPS (secondary)
        else if (std::abs(D1_rad - D2_rad) < D_TOL) mod_type = 2; // DPS
        else                                        mod_type = 3; // TPS
    }
    lastModulationType = mod_type;
    lastPhaseShiftRad = phi_rad;

    // ---- ZVS margin diagnostics ----
    double d_ratio = (V1 > 0) ? (N_main * V2_main / V1) : 1.0;
    lastVoltageConversionRatio = d_ratio;
    double phi_min_pri = (d_ratio > 0) ? (1.0 - 1.0 / d_ratio) * M_PI / 2.0 : 0.0;
    double phi_min_sec = (1.0 - d_ratio) * M_PI / 2.0;
    lastZvsMarginPrimary = std::abs(phi_rad) - phi_min_pri;
    lastZvsMarginSecondary = std::abs(phi_rad) - phi_min_sec;

    // ---- Closed-form sub-interval propagator (replaces Euler) ----
    // Build the sub-interval list for one period using the analytical
    // breakpoints {0, D1, π, π+D1, φ, φ+D2, φ+π, φ+π+D2}. iL and Im are
    // EXACTLY linear within each sub-interval, so the propagator is exact.
    auto segs = dab_build_subintervals(V1, V2_main, D1_rad, D2_rad, phi_rad);
    double iL_start = 0.0, Im_start = 0.0;
    dab_initial_conditions(segs, N_main, L, Lm, Fs, iL_start, Im_start);

    std::vector<double> bnd_theta, bnd_iL, bnd_Im;
    dab_propagate_period(segs, N_main, L, Lm, Fs, iL_start, Im_start,
                         bnd_theta, bnd_iL, bnd_Im);

    // Expose sub-interval boundary angles for diagnostics.
    lastSubIntervalTimes = bnd_theta;

    // Sample onto the existing 256-points-per-half-cycle output grid.
    const int N_samples = 256;
    int totalSamples = 2 * N_samples + 1;
    std::vector<double> time_full(totalSamples);
    std::vector<double> iL_full(totalSamples);
    std::vector<double> Vab_full(totalSamples);
    std::vector<double> Im_full(totalSamples);

    double dt = Thalf / N_samples;
    double angle_per_time = 2.0 * M_PI * Fs;

    for (int k = 0; k < totalSamples; ++k) {
        double t = k * dt;
        double theta = std::fmod(t * angle_per_time, 2.0 * M_PI);
        if (theta < 0) theta += 2.0 * M_PI;
        time_full[k] = t;
        Vab_full[k]  = Vab_at(theta, V1, D1_rad);
        iL_full[k]   = dab_sample(theta, bnd_theta, bnd_iL);
        Im_full[k]   = dab_sample(theta, bnd_theta, bnd_Im);
    }

    // ---- Primary winding excitation ----
    // Total primary current = iL (ideal-transformer) + Im (magnetizing).
    {
        std::vector<double> V_primary(totalSamples);
        std::vector<double> I_primary(totalSamples);
        for (int k = 0; k < totalSamples; ++k) {
            V_primary[k] = Vab_full[k];
            I_primary[k] = iL_full[k] + Im_full[k];
        }

        Waveform currentWaveform;
        currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        currentWaveform.set_data(I_primary);
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
    // Each output i has its own secondary bridge driving Vcd_i (amplitude V2_i)
    // through turns ratio n_i. The primary tank current iL_full was computed
    // against output 0 (the regulated/main output). For multi-output DAB this
    // is a single-loop approximation: each secondary sees a current
    // proportional to its load share, scaled by its turns ratio.
    //
    // **Multi-output limitation**: this load-share projection is only exact
    // when all outputs have the same V2 (and therefore the same n_i) AND the
    // AC current splits equally across secondaries — neither of which is
    // generally true for an under-specified multi-output topology where
    // per-secondary leakage inductance isn't given. For magnetic-component
    // design with multi-output DAB, treat the per-output current waveforms
    // as approximate; trust the primary winding waveform (which is the one
    // that determines the transformer's primary copper and core loss).
    if (turnsRatios.size() > 1) {
        static thread_local bool warned = false;
        if (!warned) {
            std::cerr << "[DAB] multi-output configuration detected — per-output "
                         "secondary current waveforms are an approximation "
                         "(load-share projection of the primary tank current). "
                         "The primary waveform and the total power are accurate; "
                         "individual secondary current magnitudes may be off by "
                         "20-100% depending on the V2 spread and the actual "
                         "(unspecified) per-secondary leakage." << std::endl;
            warned = true;
        }
    }
    double total_g = 0.0;
    for (size_t i = 0; i < turnsRatios.size() &&
                       i < dabOpPoint.get_output_voltages().size() &&
                       i < dabOpPoint.get_output_currents().size(); ++i) {
        double Vo_i = dabOpPoint.get_output_voltages()[i];
        double Io_i = dabOpPoint.get_output_currents()[i];
        if (Vo_i > 0 && Io_i > 0) total_g += Io_i / Vo_i;
    }
    if (total_g <= 0) total_g = 1.0;

    for (size_t secIdx = 0; secIdx < turnsRatios.size(); ++secIdx) {
        double n_i = turnsRatios[secIdx];
        if (n_i <= 0) continue;
        double V2_i = (secIdx < dabOpPoint.get_output_voltages().size())
                      ? dabOpPoint.get_output_voltages()[secIdx]
                      : V2_main;
        double Io_i = (secIdx < dabOpPoint.get_output_currents().size())
                      ? dabOpPoint.get_output_currents()[secIdx]
                      : 0.0;
        double share = (V2_i > 0 && Io_i > 0)
                       ? (Io_i / V2_i) / total_g
                       : (secIdx == 0 ? 1.0 : 0.0);

        std::vector<double> iSecData(totalSamples);
        std::vector<double> vSecData(totalSamples);

        for (int k = 0; k < totalSamples; ++k) {
            // Secondary current: scale primary iL by share, then by ampere-turn
            // balance the physical secondary current = n_i * (primary share).
            iSecData[k] = n_i * iL_full[k] * share;
            // Secondary bridge voltage: Vcd_i has amplitude V2_i and the same
            // phase-shifted square waveform.
            double t = time_full[k];
            double theta = t * angle_per_time;
            vSecData[k] = Vcd_at(theta, V2_i, D2_rad, phi_rad);
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

    // Inner phase shifts (degrees → radians) and outer phase shift.
    double D1_deg = dabOp.get_inner_phase_shift1().value_or(0.0);
    double D2_deg;
    if (dabOp.get_inner_phase_shift2().has_value()) {
        D2_deg = dabOp.get_inner_phase_shift2().value();
    } else {
        auto modTypeOpt = dabOp.get_modulation_type();
        D2_deg = (modTypeOpt.has_value() && modTypeOpt.value() == ModulationType::DPS) ? D1_deg : 0.0;
    }
    double D1_rad = D1_deg * M_PI / 180.0;
    double D2_rad = D2_deg * M_PI / 180.0;

    double phi_deg = dabOp.get_inner_phase_shift3().value_or(0.0);
    double phi_rad = (std::abs(phi_deg) > 1e-6) ? phi_deg * M_PI / 180.0 : computedPhaseShift;
    double phaseDelay = phi_rad / (2.0 * M_PI * Fs);  // signed

    // Time-domain inner phase shift offsets
    double t_D1 = D1_rad / (2.0 * M_PI * Fs);
    double t_D2 = D2_rad / (2.0 * M_PI * Fs);

    // Number of secondary outputs (multi-output: each output gets its own
    // secondary winding + active bridge driven by the same secondary PWMs).
    size_t numOutputs = std::min(turnsRatios.size(), dabOp.get_output_voltages().size());
    if (numOutputs == 0) numOutputs = 1;

    double V2_main = dabOp.get_output_voltages()[0];

    double L = computedSeriesInductance;
    double Lm = magnetizingInductance;

    // Simulation timing
    int periodsToExtract = get_num_periods_to_extract();
    int steadyStatePeriods = get_num_steady_state_periods();
    int numPeriodsTotal = steadyStatePeriods + periodsToExtract;
    double simTime = numPeriodsTotal * period;
    double startTime = steadyStatePeriods * period;
    double stepTime = period / 200;

    std::ostringstream circuit;

    circuit << "* Dual Active Bridge (DAB) Converter - Generated by OpenMagnetics\n";
    circuit << "* V1=" << V1 << "V, Vout=" << V2_main << "V, Fs=" << (Fs/1e3)
            << "kHz, phi=" << phi_deg << "deg, D1=" << D1_deg
            << "deg, D2=" << D2_deg << "deg, outputs=" << numOutputs << "\n";
    circuit << "* L_series=" << (L*1e6) << "uH, Lm=" << (Lm*1e6) << "uH\n\n";

    // Switch and diode models. Slightly softer body diode (RS=0.05) and
    // higher SW hysteresis (VH=0.8) reduce switching-event chatter.
    circuit << ".model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg\n";
    circuit << ".model DIDEAL D(IS=1e-12 RS=0.05)\n\n";

    // DC input voltage
    circuit << "Vdc1 vin_dc1 0 " << V1 << "\n\n";

    // ==== PRIMARY FULL BRIDGE (independent left/right legs) ====
    // For SPS (D1=0): the two legs are diagonally synchronised → standard
    //   square-wave Vab = ±V1 with 50% duty.
    // For EPS/DPS/TPS (D1>0): the LEFT leg is delayed by t_D1 relative to
    //   the RIGHT leg, creating zero-voltage plateaus during the transitions.
    //
    //   Right leg (Q3 top, Q4 bot):
    //     Q4 ON during [0, halfPeriod)             → bridge_b1 = 0
    //     Q3 ON during [halfPeriod, period)        → bridge_b1 = V1
    //   Left leg (Q1 top, Q2 bot):
    //     Q1 ON during [t_D1, t_D1 + halfPeriod)   → bridge_a1 = V1
    //     Q2 ON during [t_D1 + halfPeriod, +)      → bridge_a1 = 0
    //
    // Negative outer phi: delay the entire primary bridge by |phi|/(2π·Fs)
    // so the relative phase (secondary − primary) becomes negative.
    double priBaseDelay = (phaseDelay < 0) ? -phaseDelay : 0.0;
    double secBaseDelay = (phaseDelay > 0) ?  phaseDelay : 0.0;

    auto pulse = [&](const std::string& name, double delay) {
        circuit << "V" << name << " " << name << " 0 PULSE(0 5 "
                << std::scientific << delay << std::fixed
                << " 10n 10n "
                << std::scientific << tOn << " " << period << std::fixed << ")\n";
    };

    circuit << "* Primary Full Bridge — independent leg PWMs (supports SPS/EPS/DPS/TPS)\n";
    // Right leg (drives bridge_b1)
    pulse("pwm_p_r1", priBaseDelay);                          // Q4: bridge_b1 → 0
    pulse("pwm_p_r2", priBaseDelay + halfPeriod);             // Q3: bridge_b1 → V1
    // Left leg (drives bridge_a1, delayed by t_D1)
    pulse("pwm_p_l1", priBaseDelay + t_D1);                   // Q1: bridge_a1 → V1
    pulse("pwm_p_l2", priBaseDelay + t_D1 + halfPeriod);      // Q2: bridge_a1 → 0
    circuit << "\n";

    // Q1: top of left leg
    circuit << "S1 vin_dc1 bridge_a1 pwm_p_l1 0 SW1\n";
    circuit << "D1 0 bridge_a1 DIDEAL\n";
    // Q2: bottom of left leg
    circuit << "S2 bridge_a1 0 pwm_p_l2 0 SW1\n";
    circuit << "D2 bridge_a1 vin_dc1 DIDEAL\n";
    // Q3: top of right leg
    circuit << "S3 vin_dc1 bridge_b1 pwm_p_r2 0 SW1\n";
    circuit << "D3 0 bridge_b1 DIDEAL\n";
    // Q4: bottom of right leg
    circuit << "S4 bridge_b1 0 pwm_p_r1 0 SW1\n";
    circuit << "D4 bridge_b1 vin_dc1 DIDEAL\n";
    // Snubber RC across each switch — absorbs di/dt at switching events,
    // critical for ngspice convergence with 4-position switching (EPS/DPS/TPS).
    circuit << "Rsnub_q1 vin_dc1 bridge_a1 1k\nCsnub_q1 vin_dc1 bridge_a1 1n\n";
    circuit << "Rsnub_q2 bridge_a1 0 1k\nCsnub_q2 bridge_a1 0 1n\n";
    circuit << "Rsnub_q3 vin_dc1 bridge_b1 1k\nCsnub_q3 vin_dc1 bridge_b1 1n\n";
    circuit << "Rsnub_q4 bridge_b1 0 1k\nCsnub_q4 bridge_b1 0 1n\n\n";

    // Sense primary current (0V source in series between bridge_a1 and pri_out)
    circuit << "Vpri_sense bridge_a1 pri_out 0\n";
    // Differential bridge output Vab = V(pri_out) - V(bridge_b1) via VCVS
    circuit << "Evab vab 0 pri_out bridge_b1 1\n\n";

    // ==== SERIES INDUCTANCE ====
    circuit << "* Series inductance (leakage + external)\n";
    circuit << "L_series pri_out trafo_pri " << std::scientific << L << "\n\n";

    // ==== TRANSFORMER (multi-secondary) ====
    // One primary winding (Lm), one secondary per output. All windings share
    // the same magnetic core, so the coupling matrix must include EVERY
    // pair: each (L_pri, L_sec_oN) and each (L_sec_oN, L_sec_oM) for N≠M.
    // Without the inter-secondary K statements, ngspice treats the model as
    // independent transformers sharing only L_pri, which gives wildly wrong
    // currents (factor of ~50× off in 3-output stress tests).
    circuit << "* Transformer with " << numOutputs << " secondary winding(s)\n";
    circuit << "L_pri trafo_pri bridge_b1 " << std::scientific << Lm << "\n";
    constexpr double K_COUPLING = 0.9999;
    for (size_t i = 0; i < numOutputs; ++i) {
        double n_i = turnsRatios[i];
        if (n_i <= 0) n_i = 1.0;
        double Ls_i = Lm / (n_i * n_i);
        std::string si = std::to_string(i + 1);
        circuit << "L_sec_o" << si << " trafo_sec_a_o" << si
                << " bridge_sec_b_o" << si << " " << std::scientific << Ls_i << "\n";
    }
    // K statements: primary ↔ each secondary
    int kCount = 0;
    for (size_t i = 0; i < numOutputs; ++i) {
        ++kCount;
        circuit << "K" << kCount << " L_pri L_sec_o" << (i + 1)
                << " " << K_COUPLING << "\n";
    }
    // K statements: every pair of secondaries (N choose 2)
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

    // ==== SECONDARY FULL BRIDGES (per output) ====
    // Same independent-leg-PWM scheme as the primary, with the secondary
    // base delay = t_phi (positive phi → secondary lags primary) plus the
    // inner phase shift t_D2.
    circuit << "* Secondary bridge PWMs (D2=" << D2_deg << "deg)\n";
    pulse("pwm_s_r1", secBaseDelay);
    pulse("pwm_s_r2", secBaseDelay + halfPeriod);
    pulse("pwm_s_l1", secBaseDelay + t_D2);
    pulse("pwm_s_l2", secBaseDelay + t_D2 + halfPeriod);
    circuit << "\n";

    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        double V2_i = dabOp.get_output_voltages()[i];
        double Iout_i = dabOp.get_output_currents()[i];
        if (Iout_i <= 0) Iout_i = 1.0;
        double Rload_i = V2_i / Iout_i;
        double Cout_i = 47e-6;

        // Per-secondary current senses
        circuit << "Vsec1_sense_o" << si << " trafo_sec_a_o" << si << " sec_a_o" << si << " 0\n";
        circuit << "Vsec2_sense_o" << si << " bridge_sec_b_o" << si << " sec_b_o" << si << " 0\n";

        // Per-output secondary H-bridge
        std::string vdc2 = "vin_dc2_o" + si;
        // Q5: top of left leg
        circuit << "S5_o" << si << " " << vdc2 << " sec_a_o" << si << " pwm_s_l1 0 SW1\n";
        circuit << "D5_o" << si << " 0 sec_a_o" << si << " DIDEAL\n";
        // Q6: bottom of left leg
        circuit << "S6_o" << si << " sec_a_o" << si << " 0 pwm_s_l2 0 SW1\n";
        circuit << "D6_o" << si << " sec_a_o" << si << " " << vdc2 << " DIDEAL\n";
        // Q7: top of right leg
        circuit << "S7_o" << si << " " << vdc2 << " sec_b_o" << si << " pwm_s_r2 0 SW1\n";
        circuit << "D7_o" << si << " 0 sec_b_o" << si << " DIDEAL\n";
        // Q8: bottom of right leg
        circuit << "S8_o" << si << " sec_b_o" << si << " 0 pwm_s_r1 0 SW1\n";
        circuit << "D8_o" << si << " sec_b_o" << si << " " << vdc2 << " DIDEAL\n";
        // Snubber RC across each switch
        circuit << "Rsnub_q5_o" << si << " " << vdc2 << " sec_a_o" << si << " 1k\n";
        circuit << "Csnub_q5_o" << si << " " << vdc2 << " sec_a_o" << si << " 1n\n";
        circuit << "Rsnub_q6_o" << si << " sec_a_o" << si << " 0 1k\n";
        circuit << "Csnub_q6_o" << si << " sec_a_o" << si << " 0 1n\n";
        circuit << "Rsnub_q7_o" << si << " " << vdc2 << " sec_b_o" << si << " 1k\n";
        circuit << "Csnub_q7_o" << si << " " << vdc2 << " sec_b_o" << si << " 1n\n";
        circuit << "Rsnub_q8_o" << si << " sec_b_o" << si << " 0 1k\n";
        circuit << "Csnub_q8_o" << si << " sec_b_o" << si << " 0 1n\n";

        // Realistic Rload + Cout sink
        circuit << "Resr_o" << si << " " << vdc2 << " vout_cap_o" << si << " 0.05\n";
        circuit << "Cout_o" << si << " vout_cap_o" << si << " vout_neg_o" << si
                << " " << std::scientific << Cout_i << "\n";
        circuit << "Rload_o" << si << " vout_cap_o" << si << " vout_neg_o" << si
                << " " << Rload_i << "\n";
        circuit << "Vsec_sense_o" << si << " vout_neg_o" << si << " 0 0\n\n";
    }

    // Transient Analysis
    circuit << ".tran " << std::scientific << stepTime
            << " " << simTime << " " << startTime << "\n\n";

    // Saved signals
    circuit << ".save v(vab) v(trafo_pri) v(bridge_b1) v(vin_dc1)"
               " i(Vpri_sense) i(Vdc1)";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        circuit << " v(sec_a_o" << si << ") v(sec_b_o" << si << ")"
                << " v(vin_dc2_o" << si << ") v(vout_cap_o" << si << ")"
                << " i(Vsec1_sense_o" << si << ") i(Vsec2_sense_o" << si << ")"
                << " i(Vsec_sense_o" << si << ")";
    }
    circuit << "\n\n";

    // Solver options for convergence in switching circuits.
    // METHOD=GEAR + larger TRTOL handle the 4-position switching robustly;
    // higher ITL4 lets the solver retry hard switching transients.
    circuit << ".options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500\n";
    circuit << ".options METHOD=GEAR TRTOL=7\n";
    // Initial conditions: pre-charge each output capacitor to its target voltage.
    circuit << ".ic";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        double V2_i = dabOp.get_output_voltages()[i];
        circuit << " v(vout_cap_o" << si << ")=" << V2_i
                << " v(vin_dc2_o" << si << ")=" << V2_i;
    }
    circuit << "\n\n";

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
    std::vector<OperatingPoint> operatingPoints;
    NgspiceRunner runner;
    if (!runner.is_available()) {
        std::cerr << "[DAB] ngspice not available — falling back to analytical "
                     "operating-point generation" << std::endl;
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    // Make sure design quantities (computedSeriesInductance, computedPhaseShift)
    // are populated. process_design_requirements is idempotent.
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
                throw std::runtime_error("DAB: operating point has invalid switching frequency");
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
                std::cerr << "[DAB] ngspice run failed: " << simResult.errorMessage
                          << " — falling back to analytical for this op." << std::endl;
                auto analytical = process_operating_point_for_input_voltage(
                    inputVoltages[vinIdx], ops[opIdx], turnsRatios, magnetizingInductance);
                analytical.set_name(inputVoltageNames[vinIdx] + " input ("
                    + std::to_string(static_cast<int>(inputVoltages[vinIdx])) + "V) [analytical]");
                operatingPoints.push_back(analytical);
                continue;
            }

            // Build the per-winding waveform mapping for the multi-secondary
            // netlist. The primary uses the differential v(vab) sense source.
            // Each output gets its own per-half winding sense (Vsec1_sense_oN
            // and Vsec2_sense_oN), but for a single Secondary winding per
            // output we map to the first half-sense and the sec_a node.
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
                waveformMapping.push_back({{"voltage", "sec_a_o" + si},
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
    
    // Use nominal input voltage only (index 0) for topology waveforms — same as LLC.
    // Min/max variants are covered by simulate_and_extract_operating_points.
    for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
        const size_t inputVoltageIndex = 0;
        {
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
            
            // Names are normalized by NgspiceRunner::parse_raw_file:
            //   v(node)   -> node
            //   i(source) -> source#branch
            // input_voltage: differential bridge output v(vab) (bipolar ±V1).
            // input_current: AC tank current via Vpri_sense (bipolar). The
            //   actual DC supply current is also available as i(Vdc1) for
            //   callers that need it.
            wf.set_input_voltage(getWaveform("vab"));
            wf.set_input_current(getWaveform("vpri_sense#branch"));

            // Multi-secondary output: each output has its own vout_cap_oN and
            // vsec_sense_oN.
            size_t numOutputs = std::min(turnsRatios.size(), opPoint.get_output_voltages().size());
            if (numOutputs == 0) numOutputs = 1;
            for (size_t outIdx = 0; outIdx < numOutputs; ++outIdx) {
                std::string si = std::to_string(outIdx + 1);
                wf.get_mutable_output_voltages().push_back(getWaveform("vout_cap_o" + si));
                wf.get_mutable_output_currents().push_back(getWaveform("vsec_sense_o" + si + "#branch"));
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
    LmTol.set_minimum(desiredMagnetizingInductance);
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
