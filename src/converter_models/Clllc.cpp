// Version: 2.0 — Phase A: 4-state linear-ODE analytical solver (forward direction).
//
// STATUS: Phase A of CLLLC v2 (per CLLLC_PLAN.md §A.5).
//
// What works in this version
// --------------------------
//  * process_design_requirements (already in v1)
//  * Static analytical helpers (already in v1)
//  * process_operating_points + process_operating_point_for_input_voltage
//    via a one-shot affine-propagator steady-state solver on a 4-state
//    linear ODE [i_Lr1, i_Lr2, v_Cr1, v_Cr2] (i_Lm derived as i_Lr1 - i_Lr2/n).
//  * Per-OP diagnostics: peak / RMS currents, magnetizing peak, Cr1/Cr2 peak
//    voltages, ZVS margins on BOTH bridges, current-sharing ratio, operating
//    region (above/at/below resonance).
//  * Forward + reverse direction (reverse uses the symmetric-tank reflection
//    trick: swap which side is "active" and re-use the same solver).
//
// Pending for later Phase(s)
// --------------------------
//  * SPICE netlist generation (`generate_ngspice_circuit`) — STILL throws.
//  * `simulate_and_extract_*` SPICE wrappers — STILL throw.
//  * `get_extra_components_inputs` — STILL throws (needs Cr1/Cr2 RMS + peak
//    waveforms; the analytical solver here populates extraCr*VoltageWaveforms,
//    but the CAS-Inputs builder is left for the next phase).
//  * AdvancedClllc::process — STILL throws.
//  * Below-resonance operation: v1 assumes both bridges flip synchronously at
//    the period boundary, which is exact at and above the primary resonant
//    frequency. Below resonance, i_Lr2 zero-crosses inside the half period and
//    the LV bridge would commutate early. The solver reports |i_Lr2| at the
//    period boundary as `lastZvsMarginSecondaryLagging`; if the user's
//    operating point is well below resonance this margin loses physical
//    meaning. Documented limit: fsw / fr ≥ 0.95.
//
// Solver design rationale (per OpenMagnetics' no-shortcuts rule)
// --------------------------------------------------------------
// The ODE between bridge transitions is purely linear with constant
// coefficients (LCLCL tank driven by step bridge voltages). Its half-period
// propagator P : x0 → x(Thalf) is therefore affine: P(x0) = M·x0 + g for some
// 4×4 matrix M and 4-vector g. We measure M and g by running the propagator
// from five carefully chosen initial conditions (zero plus the four basis
// vectors), then solve the half-wave-antisymmetric steady-state condition
// P(x0) = -x0, i.e. (M+I)·x0 = -g, by 4×4 Gaussian elimination. This produces
// an exact (modulo RK4 truncation) steady state in 6 RK4 passes — no Newton,
// no Picard, no fallbacks. A small parasitic ESR is added to Lr1 and Lr2
// (1 mΩ each, far smaller than any physical wire) so M has spectral radius
// strictly < 1; without this guard, M+I would be singular for an idealised
// lossless tank at exact resonance.

#include "converter_models/Clllc.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace OpenMagnetics {

// =====================================================================
// Construction
// =====================================================================
Clllc::Clllc(const json& j) {
    from_json(j, *static_cast<ClllcResonant*>(this));
}

AdvancedClllc::AdvancedClllc(const json& j) {
    from_json(j, *this);
}

// =====================================================================
// Bridge / frequency / ratio accessors
// =====================================================================
double Clllc::get_primary_bridge_voltage_factor() const {
    auto bt = get_bridge_type_primary();
    if (bt.has_value() && bt.value() == LlcBridgeType::FULL_BRIDGE) return 1.0;
    return 0.5;
}

double Clllc::get_secondary_bridge_voltage_factor() const {
    auto bt = get_bridge_type_secondary();
    if (bt.has_value() && bt.value() == LlcBridgeType::FULL_BRIDGE) return 1.0;
    return 0.5;
}

double Clllc::get_effective_resonant_frequency() const {
    if (get_primary_resonant_frequency().has_value())
        return get_primary_resonant_frequency().value();
    return std::sqrt(get_min_switching_frequency() * get_max_switching_frequency());
}

double Clllc::get_effective_inductance_ratio_k() const {
    auto kBase = MAS::ClllcResonant::get_inductance_ratio_k();
    if (kBase.has_value() && kBase.value() > 0) return kBase.value();
    return computedInductanceRatioK > 0 ? computedInductanceRatioK : 5.0;
}

double Clllc::get_effective_tank_symmetry_ratio() const {
    auto s = get_tank_symmetry_ratio();
    if (s.has_value()) return s.value();
    return 1.0;  // symmetric by default — the headline CLLLC property.
}

// =====================================================================
// Static analytical helpers (these CAN be unit-tested in isolation
// before the solver is finished — they're pure math.)
// =====================================================================
double Clllc::compute_primary_resonant_frequency(double Lr1, double Cr1) {
    if (Lr1 <= 0 || Cr1 <= 0)
        throw std::invalid_argument(
            "Clllc::compute_primary_resonant_frequency: Lr1 and Cr1 must be > 0");
    return 1.0 / (2.0 * M_PI * std::sqrt(Lr1 * Cr1));
}

double Clllc::compute_magnetizing_resonant_frequency_fwd(
        double Lr1, double Lm, double Cr1) {
    if (Lr1 <= 0 || Lm <= 0 || Cr1 <= 0)
        throw std::invalid_argument(
            "Clllc::compute_magnetizing_resonant_frequency_fwd: positive args required");
    return 1.0 / (2.0 * M_PI * std::sqrt((Lr1 + Lm) * Cr1));
}

double Clllc::compute_magnetizing_resonant_frequency_rev(
        double Lr2, double Lm, double n, double Cr2) {
    if (Lr2 <= 0 || Lm <= 0 || n <= 0 || Cr2 <= 0)
        throw std::invalid_argument(
            "Clllc::compute_magnetizing_resonant_frequency_rev: positive args required");
    // From the LV-side reference frame: Lm reflected as Lm/n^2.
    return 1.0 / (2.0 * M_PI * std::sqrt((Lr2 + Lm / (n * n)) * Cr2));
}

double Clllc::compute_fha_gain(double f_n, double K, double Q,
                                double a, double b) {
    if (f_n <= 0 || K <= 0)
        throw std::invalid_argument("Clllc::compute_fha_gain: f_n and K must be > 0");
    // Symmetric (a = b = 1) reduces to LLC FHA shape; asymmetric per
    // IEEE 10362332 / MDPI Electronics 12(7):1605 §3 closed form.
    // v1 uses the symmetric reduction and applies a small (a+b)/2 weighting
    // for the asymmetric branch — refine when the asymmetric solver lands.
    double ab = 0.5 * (a + b);
    double term1 = 1.0 + ab * K * (1.0 - 1.0 / (f_n * f_n));
    double term2 = Q * (f_n - 1.0 / f_n);
    double denom = std::sqrt(term1 * term1 + term2 * term2);
    if (denom <= 0)
        throw std::runtime_error("Clllc::compute_fha_gain: denominator collapsed");
    return 1.0 / denom;
}

double Clllc::compute_inductance_ratio_K(double Lm, double Lr1) {
    if (Lr1 <= 0 || Lm <= 0)
        throw std::invalid_argument("Clllc::compute_inductance_ratio_K: Lm and Lr1 must be > 0");
    return Lm / Lr1;
}

double Clllc::compute_quality_factor(double Lr1, double Cr1, double R_load_reflected) {
    if (Lr1 <= 0 || Cr1 <= 0 || R_load_reflected <= 0)
        throw std::invalid_argument("Clllc::compute_quality_factor: positive args required");
    return std::sqrt(Lr1 / Cr1) / R_load_reflected;
}

double Clllc::compute_zvs_lm_max(double T_dead, double C_oss, double Fs) {
    if (T_dead <= 0 || C_oss <= 0 || Fs <= 0)
        throw std::invalid_argument("Clllc::compute_zvs_lm_max: positive args required");
    // Huang SLUP263 §6 (extended for CLLLC mirror): Lm <= T_dead / (16*Coss*Fs)
    return T_dead / (16.0 * C_oss * Fs);
}

double Clllc::compute_turns_ratio(double V_HV, double V_LV) {
    if (V_HV <= 0 || V_LV <= 0)
        throw std::invalid_argument("Clllc::compute_turns_ratio: positive args required");
    return V_HV / V_LV;
}

double Clllc::compute_symmetric_lr2(double Lr1, double n) {
    if (Lr1 <= 0 || n <= 0)
        throw std::invalid_argument("Clllc::compute_symmetric_lr2: positive args required");
    return Lr1 / (n * n);
}

double Clllc::compute_symmetric_cr2(double Cr1, double n) {
    if (Cr1 <= 0 || n <= 0)
        throw std::invalid_argument("Clllc::compute_symmetric_cr2: positive args required");
    return Cr1 * (n * n);
}

// =====================================================================
// run_checks — validates the input deck (analogue of Llc::run_checks)
// =====================================================================
bool Clllc::run_checks(bool assertErrors) {
    _assertErrors = assertErrors;
    bool ok = true;
    auto& ops = get_operating_points();
    if (ops.empty()) {
        if (assertErrors) throw std::runtime_error("CLLLC: no operating points");
        return false;
    }
    for (size_t i = 0; i < ops.size(); ++i) {
        auto& op = ops[i];
        if (op.get_output_voltages().empty() || op.get_output_currents().empty()) {
            if (assertErrors) throw std::runtime_error("CLLLC: OP missing voltages/currents");
            ok = false;
        }
        if (op.get_output_voltages().size() != op.get_output_currents().size()) {
            if (assertErrors) throw std::runtime_error("CLLLC: voltage/current count mismatch");
            ok = false;
        }
        double fsw = op.get_switching_frequency();
        if (fsw < get_min_switching_frequency() * 0.99 ||
            fsw > get_max_switching_frequency() * 1.01) {
            if (assertErrors) throw std::runtime_error("CLLLC: fsw out of range");
            ok = false;
        }
    }
    auto sym = get_effective_tank_symmetry_ratio();
    if (sym <= 0) {
        if (assertErrors) throw std::runtime_error("CLLLC: tankSymmetryRatio must be > 0");
        ok = false;
    }
    return ok;
}

// =====================================================================
// Design Requirements — derives Lr1, Cr1, Lm, Lr2, Cr2 and turns ratio.
// Symmetric: Lr2 = Lr1/n^2, Cr2 = Cr1*n^2.
// =====================================================================
DesignRequirements Clllc::process_design_requirements() {
    auto& ops = get_operating_points();
    if (ops.empty())
        throw std::runtime_error("CLLLC: at least one operating point required");

    auto& hv = get_high_voltage_bus_voltage();
    auto& lv = get_low_voltage_bus_voltage();
    if (!hv.get_nominal().has_value() &&
        !(hv.get_minimum().has_value() && hv.get_maximum().has_value()))
        throw std::runtime_error("CLLLC: highVoltageBusVoltage requires nominal or min+max");
    if (!lv.get_nominal().has_value() &&
        !(lv.get_minimum().has_value() && lv.get_maximum().has_value()))
        throw std::runtime_error("CLLLC: lowVoltageBusVoltage requires nominal or min+max");

    double Vhv_nom = hv.get_nominal().value_or(
        (hv.get_minimum().value_or(0) + hv.get_maximum().value_or(0)) / 2.0);
    double Vlv_nom = lv.get_nominal().value_or(
        (lv.get_minimum().value_or(0) + lv.get_maximum().value_or(0)) / 2.0);

    // Turns ratio sized at nominal so M = 1 at f_r1 (Clllc plan §A.5).
    double k_pri = get_primary_bridge_voltage_factor();
    double k_sec = get_secondary_bridge_voltage_factor();
    double n = (Vhv_nom * k_pri) / (Vlv_nom * k_sec);
    if (n <= 0)
        throw std::runtime_error("CLLLC: derived turns ratio non-positive");
    computedTurnsRatio = n;

    // Pick worst-case OP for tank sizing — max |Pout| across all OPs.
    double Pmax = 0.0;
    double Vlv_design = Vlv_nom;
    for (auto& op : ops) {
        double V = op.get_output_voltages()[0];
        double I = op.get_output_currents()[0];
        double P = std::abs(V * I);
        if (P > Pmax) {
            Pmax = P;
            Vlv_design = V;
        }
    }
    if (Pmax <= 0)
        throw std::runtime_error("CLLLC: no positive-power operating point found");

    double Iload = Pmax / Vlv_design;
    double Rload = Vlv_design / Iload;
    double Rac_pri = (8.0 * n * n) / (M_PI * M_PI) * Rload;

    double Q  = get_quality_factor().value_or(0.4);
    double K  = MAS::ClllcResonant::get_inductance_ratio_k().value_or(5.0);
    double fr = get_effective_resonant_frequency();

    double Zr  = Q * Rac_pri;
    double Lr1 = (userPrimarySeriesInductance > 0)
                  ? userPrimarySeriesInductance
                  : Zr / (2.0 * M_PI * fr);
    double Cr1 = (userPrimaryResonantCapacitance > 0)
                  ? userPrimaryResonantCapacitance
                  : 1.0 / (2.0 * M_PI * fr * Zr);
    double Lm  = K * Lr1;

    // Symmetric tank
    double sym = get_effective_tank_symmetry_ratio();
    double Lr2 = compute_symmetric_lr2(Lr1, n);
    double Cr2 = compute_symmetric_cr2(Cr1, n);
    if (std::abs(sym - 1.0) > 1e-9) {
        // Asymmetric tank: scale Lr2/Cr2 to the requested ratio.
        // sym = (Lr2/(n^2 Lr1)) * (Cr2 n^2 / Cr1). Apply sqrt(sym) to each.
        double s = std::sqrt(sym);
        Lr2 *= s;
        Cr2 *= s;
    }

    computedPrimarySeriesInductance      = Lr1;
    computedPrimaryResonantCapacitance   = Cr1;
    computedSecondarySeriesInductance    = Lr2;
    computedSecondaryResonantCapacitance = Cr2;
    computedMagnetizingInductance        = Lm;
    computedInductanceRatioK             = K;
    computedQualityFactor                = Q;
    computedPrimaryResonantFrequency     = compute_primary_resonant_frequency(Lr1, Cr1);
    computedMagnetizingResonantFreqFwd   =
        compute_magnetizing_resonant_frequency_fwd(Lr1, Lm, Cr1);
    computedMagnetizingResonantFreqRev   =
        compute_magnetizing_resonant_frequency_rev(Lr2, Lm, n, Cr2);

    DesignRequirements designRequirements;
    designRequirements.set_topology(Topologies::CLLLC_RESONANT_CONVERTER);

    std::vector<double> turnsRatios = { n };
    DimensionWithTolerance turnsRatioDim;
    turnsRatioDim.set_nominal(n);
    designRequirements.set_turns_ratios({turnsRatioDim});

    DimensionWithTolerance LmDim;
    LmDim.set_nominal(Lm);
    designRequirements.set_magnetizing_inductance(LmDim);

    return designRequirements;
}

// =====================================================================
// 4-state linear-ODE solver — Phase A core.
//
// State vector x = [i_Lr1, i_Lr2, v_Cr1, v_Cr2]ᵀ. Magnetizing current is
// algebraically derived as i_Lm = i_Lr1 - i_Lr2/n (transformer ampere-turn
// balance). Lm therefore enters the dynamics through the coupled inductor
// matrix below, not as a separate state.
//
// Topology (forward power flow, both bridges full-bridge):
//   u_HV = k_pri · V_HV · sgn(half-cycle)   ← HV bridge active
//   u_LV = k_sec · V_LV · sgn(half-cycle)   ← LV bridge synchronous-rectifying
//
// KVL HV loop: u_HV = Lr1·di_Lr1/dt + r_Lr1·i_Lr1 + v_Cr1 + v_pri
// KVL LV loop: v_pri/n = Lr2·di_Lr2/dt + r_Lr2·i_Lr2 + v_Cr2 + u_LV
// v_pri = Lm·di_Lm/dt = Lm·(di_Lr1/dt - di_Lr2/dt / n)
//
// Eliminating v_pri gives the 2×2 coupled inductance system
//   M_L · [di_Lr1/dt; di_Lr2/dt] = [u_HV - v_Cr1 - r·i_Lr1;
//                                   -v_Cr2 - u_LV - r·i_Lr2]
// where M_L = [[Lr1+Lm, -Lm/n], [-Lm/n, Lr2+Lm/n²]].
// =====================================================================

namespace {

struct ClllcState {
    double iLr1 = 0;
    double iLr2 = 0;
    double vCr1 = 0;
    double vCr2 = 0;
};

struct ClllcParams {
    double Lr1 = 0;
    double Lr2 = 0;
    double Lm  = 0;
    double Cr1 = 0;
    double Cr2 = 0;
    double n   = 1;
    double rLr1 = 1e-3;
    double rLr2 = 1e-3;
    double k_pri = 1.0;
    double k_sec = 1.0;
};

// Compute dx/dt at state x with raw bridge voltages (will be scaled by k_pri / k_sec).
// uHV_raw and uLV_raw are the full-side bridge voltages (V_HV / V_LV with sign).
static ClllcState clllc_dxdt(const ClllcState& x, double uHV_raw, double uLV_raw,
                              const ClllcParams& p) {
    double uHV = p.k_pri * uHV_raw;
    double uLV = p.k_sec * uLV_raw;
    double Lm_n  = p.Lm / p.n;
    double Lm_n2 = p.Lm / (p.n * p.n);
    double det = (p.Lr1 + p.Lm) * (p.Lr2 + Lm_n2) - Lm_n * Lm_n;
    if (det <= 0)
        throw std::runtime_error("Clllc::clllc_dxdt: coupled-L matrix singular");

    double rhs1 = uHV - x.vCr1 - p.rLr1 * x.iLr1;
    double rhs2 = -x.vCr2 - uLV - p.rLr2 * x.iLr2;

    ClllcState d;
    d.iLr1 = ((p.Lr2 + Lm_n2) * rhs1 + Lm_n * rhs2) / det;
    d.iLr2 = (Lm_n * rhs1 + (p.Lr1 + p.Lm) * rhs2) / det;
    d.vCr1 = x.iLr1 / p.Cr1;
    d.vCr2 = x.iLr2 / p.Cr2;
    return d;
}

static inline ClllcState state_axpy(const ClllcState& a, const ClllcState& b, double s) {
    return ClllcState{a.iLr1 + s*b.iLr1, a.iLr2 + s*b.iLr2,
                      a.vCr1 + s*b.vCr1, a.vCr2 + s*b.vCr2};
}

// Single RK4 step of size dt with constant bridge inputs.
static ClllcState rk4_step(const ClllcState& x, double dt, double uHV, double uLV,
                            const ClllcParams& p) {
    ClllcState k1 = clllc_dxdt(x, uHV, uLV, p);
    ClllcState k2 = clllc_dxdt(state_axpy(x, k1, dt * 0.5), uHV, uLV, p);
    ClllcState k3 = clllc_dxdt(state_axpy(x, k2, dt * 0.5), uHV, uLV, p);
    ClllcState k4 = clllc_dxdt(state_axpy(x, k3, dt), uHV, uLV, p);
    return ClllcState{
        x.iLr1 + dt / 6.0 * (k1.iLr1 + 2*k2.iLr1 + 2*k3.iLr1 + k4.iLr1),
        x.iLr2 + dt / 6.0 * (k1.iLr2 + 2*k2.iLr2 + 2*k3.iLr2 + k4.iLr2),
        x.vCr1 + dt / 6.0 * (k1.vCr1 + 2*k2.vCr1 + 2*k3.vCr1 + k4.vCr1),
        x.vCr2 + dt / 6.0 * (k1.vCr2 + 2*k2.vCr2 + 2*k3.vCr2 + k4.vCr2),
    };
}

// Propagate a half-period with both bridges in the +polarity. Caller produces
// the negative half by mirror-symmetry (half-wave antisymmetry of the steady
// state).
static ClllcState propagate_half(ClllcState x0, double Thalf, double VHV, double VLV,
                                  const ClllcParams& p, int N) {
    double dt = Thalf / N;
    ClllcState x = x0;
    for (int k = 0; k < N; ++k) {
        x = rk4_step(x, dt, +VHV, +VLV, p);
    }
    return x;
}

// 4×4 linear solve via Gaussian elimination with partial pivoting. Throws on
// singular matrix.
static void solve4x4(double A[4][4], double b[4], double out[4]) {
    double M[4][5];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) M[i][j] = A[i][j];
        M[i][4] = b[i];
    }
    for (int col = 0; col < 4; ++col) {
        int piv = col;
        double mx = std::abs(M[col][col]);
        for (int r = col + 1; r < 4; ++r) {
            if (std::abs(M[r][col]) > mx) { mx = std::abs(M[r][col]); piv = r; }
        }
        if (mx < 1e-30)
            throw std::runtime_error("Clllc steady-state: M+I singular (lossless tank?)");
        if (piv != col) for (int k = 0; k < 5; ++k) std::swap(M[col][k], M[piv][k]);
        for (int r = col + 1; r < 4; ++r) {
            double f = M[r][col] / M[col][col];
            for (int k = col; k < 5; ++k) M[r][k] -= f * M[col][k];
        }
    }
    for (int i = 3; i >= 0; --i) {
        double s = M[i][4];
        for (int j = i + 1; j < 4; ++j) s -= M[i][j] * out[j];
        out[i] = s / M[i][i];
    }
}

// One-shot steady-state x0 such that P(x0) = -x0 (half-wave antisymmetry).
//
// P is affine: P(x0) = M·x0 + g. We measure g (= P(0)) and the four columns of
// M (= P(e_i) - g) by 5 RK4 passes, then solve (M + I)·x0 = -g.
static ClllcState solve_steady_state(double Thalf, double VHV, double VLV,
                                      const ClllcParams& p, int N) {
    ClllcState zero;
    ClllcState ge = propagate_half(zero, Thalf, VHV, VLV, p, N);
    double g[4] = {ge.iLr1, ge.iLr2, ge.vCr1, ge.vCr2};

    double M[4][4];
    for (int i = 0; i < 4; ++i) {
        ClllcState ei;
        if (i == 0) ei.iLr1 = 1.0;
        else if (i == 1) ei.iLr2 = 1.0;
        else if (i == 2) ei.vCr1 = 1.0;
        else ei.vCr2 = 1.0;
        ClllcState end = propagate_half(ei, Thalf, VHV, VLV, p, N);
        M[0][i] = end.iLr1 - g[0];
        M[1][i] = end.iLr2 - g[1];
        M[2][i] = end.vCr1 - g[2];
        M[3][i] = end.vCr2 - g[3];
    }

    double A[4][4];
    double rhs[4];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) A[i][j] = M[i][j] + (i == j ? 1.0 : 0.0);
        rhs[i] = -g[i];
    }
    double x0[4];
    solve4x4(A, rhs, x0);
    return ClllcState{x0[0], x0[1], x0[2], x0[3]};
}

}  // namespace

// =====================================================================
// process_operating_points — outer loop over input voltages
// =====================================================================
std::vector<OperatingPoint> Clllc::process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
    if (computedPrimarySeriesInductance <= 0 || computedPrimaryResonantCapacitance <= 0)
        process_design_requirements();

    extraLr1CurrentWaveforms.clear();
    extraLr1VoltageWaveforms.clear();
    extraLr2CurrentWaveforms.clear();
    extraLr2VoltageWaveforms.clear();
    extraCr1CurrentWaveforms.clear();
    extraCr1VoltageWaveforms.clear();
    extraCr2CurrentWaveforms.clear();
    extraCr2VoltageWaveforms.clear();
    extraLmCurrentWaveforms.clear();
    extraTimeVectors.clear();

    std::vector<OperatingPoint> result;
    auto& ops = get_operating_points();
    if (ops.empty())
        throw std::runtime_error("Clllc::process_operating_points: no operating points");

    auto& hvSpec = get_high_voltage_bus_voltage();
    auto& lvSpec = get_low_voltage_bus_voltage();

    auto pick_voltages = [](const DimensionWithTolerance& spec) {
        std::vector<std::pair<double, std::string>> out;
        if (spec.get_nominal().has_value())
            out.push_back({spec.get_nominal().value(), "Nominal"});
        if (spec.get_minimum().has_value())
            out.push_back({spec.get_minimum().value(), "Min"});
        if (spec.get_maximum().has_value())
            out.push_back({spec.get_maximum().value(), "Max"});
        std::sort(out.begin(), out.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        auto last = std::unique(out.begin(), out.end(),
                  [](const auto& a, const auto& b) { return a.first == b.first; });
        out.erase(last, out.end());
        return out;
    };

    auto hvVoltages = pick_voltages(hvSpec);
    auto lvVoltages = pick_voltages(lvSpec);
    if (hvVoltages.empty() || lvVoltages.empty())
        throw std::runtime_error(
            "Clllc::process_operating_points: HV/LV bus voltages need nominal or min+max");

    for (auto& op : ops) {
        auto dir = op.get_power_flow_direction();
        bool reverse = dir.has_value() && dir.value() == CllcPowerFlow::REVERSE;
        // Pick the input voltage spec corresponding to the active side.
        const auto& inVoltages = reverse ? lvVoltages : hvVoltages;
        for (const auto& [Vin, name] : inVoltages) {
            auto sample = process_operating_point_for_input_voltage(
                Vin, op, turnsRatios, magnetizingInductance);
            std::ostringstream tag;
            tag << name << (reverse ? " LV" : " HV") << " input ("
                << static_cast<int>(Vin) << " V, "
                << (reverse ? "reverse" : "forward") << ")";
            sample.set_name(tag.str());
            result.push_back(sample);
        }
    }
    return result;
}

// =====================================================================
// process_operating_point_for_input_voltage — solver entry point
// =====================================================================
OperatingPoint Clllc::process_operating_point_for_input_voltage(
        double inputVoltage,
        const ClllcOperatingPoint& op,
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
    OperatingPoint operatingPoint;

    if (op.get_output_voltages().empty() || op.get_output_currents().empty())
        throw std::runtime_error(
            "Clllc::process_operating_point_for_input_voltage: OP missing voltages/currents");

    double fsw = op.get_switching_frequency();
    if (!std::isfinite(fsw) || fsw <= 0)
        throw std::runtime_error(
            "Clllc::process_operating_point_for_input_voltage: switching frequency must be > 0");

    if (computedPrimarySeriesInductance <= 0 || computedPrimaryResonantCapacitance <= 0)
        throw std::runtime_error(
            "Clllc::process_operating_point_for_input_voltage: tank not designed; "
            "call process_design_requirements first");

    double n = (turnsRatios.empty()) ? computedTurnsRatio : turnsRatios[0];
    if (!std::isfinite(n) || n <= 0)
        throw std::runtime_error(
            "Clllc::process_operating_point_for_input_voltage: turns ratio must be > 0");

    auto dirOpt = op.get_power_flow_direction();
    bool reverse = dirOpt.has_value() && dirOpt.value() == CllcPowerFlow::REVERSE;

    // In reverse direction we exploit the symmetric tank: swap "primary" and
    // "secondary" labels for the solver. The schema's outputVoltages/currents
    // describe the LV side; in reverse they become the source side.
    double Vlv = op.get_output_voltages()[0];
    double Vhv = inputVoltage;
    if (reverse) {
        // inputVoltage IS the LV-side bus (active); the HV side is the load.
        Vlv = inputVoltage;
        Vhv = op.get_output_voltages()[0];
    }
    if (Vlv <= 0 || Vhv <= 0)
        throw std::runtime_error("Clllc: HV/LV bus voltages must be positive");

    ClllcParams p;
    if (!reverse) {
        // Forward: HV bridge active, LV bridge SR.
        p.Lr1 = computedPrimarySeriesInductance;
        p.Lr2 = computedSecondarySeriesInductance;
        p.Cr1 = computedPrimaryResonantCapacitance;
        p.Cr2 = computedSecondaryResonantCapacitance;
        p.n   = n;
        p.k_pri = get_primary_bridge_voltage_factor();
        p.k_sec = get_secondary_bridge_voltage_factor();
    } else {
        // Reverse: reflect the system to the LV-side reference frame. With
        // n_rev = 1/n the roles of the primary and secondary tank swap.
        p.Lr1 = computedSecondarySeriesInductance;
        p.Lr2 = computedPrimarySeriesInductance;
        p.Cr1 = computedSecondaryResonantCapacitance;
        p.Cr2 = computedPrimaryResonantCapacitance;
        p.n   = 1.0 / n;
        p.k_pri = get_secondary_bridge_voltage_factor();
        p.k_sec = get_primary_bridge_voltage_factor();
    }
    p.Lm  = (magnetizingInductance > 0) ? magnetizingInductance : computedMagnetizingInductance;
    if (reverse) p.Lm = p.Lm / (n * n);  // reflect Lm to LV-side reference frame
    if (!std::isfinite(p.Lm) || p.Lm <= 0)
        throw std::runtime_error("Clllc: magnetizing inductance must be > 0");
    p.rLr1 = 1e-3;  // 1 mΩ damping — cosmetic; prevents (M+I) singularity at exact LIP
    p.rLr2 = 1e-3;

    double VHV_drv = reverse ? Vlv : Vhv;
    double VLV_drv = reverse ? Vhv : Vlv;

    double Thalf = 0.5 / fsw;
    constexpr int N = 512;
    double dt = Thalf / N;

    ClllcState x0 = solve_steady_state(Thalf, VHV_drv, VLV_drv, p, N);

    // Re-propagate to fill positive half waveforms.
    std::vector<double> tHalf(N + 1), iLr1H(N + 1), iLr2H(N + 1),
                        iLmH(N + 1), vCr1H(N + 1), vCr2H(N + 1);
    auto fill_at = [&](int k, const ClllcState& s) {
        tHalf[k]  = k * dt;
        iLr1H[k]  = s.iLr1;
        iLr2H[k]  = s.iLr2;
        iLmH[k]   = s.iLr1 - s.iLr2 / p.n;
        vCr1H[k]  = s.vCr1;
        vCr2H[k]  = s.vCr2;
    };
    ClllcState x = x0;
    fill_at(0, x);
    for (int k = 1; k <= N; ++k) {
        x = rk4_step(x, dt, +VHV_drv, +VLV_drv, p);
        fill_at(k, x);
    }

    // Half-wave antisymmetry residual (sanity check; should be ~0 to 1e-6 of
    // the peak current by construction).
    double residual = std::sqrt(
        std::pow(x.iLr1 + x0.iLr1, 2) + std::pow(x.iLr2 + x0.iLr2, 2) +
        std::pow(x.vCr1 + x0.vCr1, 2) + std::pow(x.vCr2 + x0.vCr2, 2));
    lastSteadyStateResidual = residual;

    // Build full-period (2N+1 samples) by mirroring.
    int total = 2 * N + 1;
    std::vector<double> tF(total), iLr1F(total), iLr2F(total), iLmF(total),
                        vCr1F(total), vCr2F(total),
                        vPriF(total), vSecF(total);
    for (int k = 0; k <= N; ++k) {
        tF[k] = tHalf[k];
        iLr1F[k] = iLr1H[k]; iLr2F[k] = iLr2H[k]; iLmF[k] = iLmH[k];
        vCr1F[k] = vCr1H[k]; vCr2F[k] = vCr2H[k];
    }
    for (int k = 1; k <= N; ++k) {
        tF[N + k]   = Thalf + k * dt;
        iLr1F[N + k] = -iLr1H[k]; iLr2F[N + k] = -iLr2H[k]; iLmF[N + k] = -iLmH[k];
        vCr1F[N + k] = -vCr1H[k]; vCr2F[N + k] = -vCr2H[k];
    }

    // Recover primary/secondary winding voltages from the ODE: v_pri = Lm·di_Lm/dt.
    for (int k = 0; k < total; ++k) {
        double uHV_full = (k <= N) ? +VHV_drv : -VHV_drv;
        double uLV_full = (k <= N) ? +VLV_drv : -VLV_drv;
        ClllcState xk{iLr1F[k], iLr2F[k], vCr1F[k], vCr2F[k]};
        ClllcState d = clllc_dxdt(xk, uHV_full, uLV_full, p);
        double diLm = d.iLr1 - d.iLr2 / p.n;
        vPriF[k] = p.Lm * diLm;
        vSecF[k] = vPriF[k] / p.n;
    }

    // ----- Diagnostics -----
    double iLr1pk = 0, iLr2pk = 0, iLmpk = 0;
    double iLr1sq = 0, iLr2sq = 0;
    double vCr1pk = 0, vCr2pk = 0;
    for (int k = 0; k < total; ++k) {
        iLr1pk = std::max(iLr1pk, std::abs(iLr1F[k]));
        iLr2pk = std::max(iLr2pk, std::abs(iLr2F[k]));
        iLmpk  = std::max(iLmpk,  std::abs(iLmF[k]));
        vCr1pk = std::max(vCr1pk, std::abs(vCr1F[k]));
        vCr2pk = std::max(vCr2pk, std::abs(vCr2F[k]));
        iLr1sq += iLr1F[k] * iLr1F[k];
        iLr2sq += iLr2F[k] * iLr2F[k];
    }

    // Map solver-frame diagnostics back to the original (primary / secondary)
    // labels that callers expect. In reverse direction the solver's "primary"
    // (active) is actually the LV side, so we swap.
    if (!reverse) {
        lastPrimaryPeakCurrent     = iLr1pk;
        lastSecondaryPeakCurrent   = iLr2pk;
        lastPrimaryRmsCurrent      = std::sqrt(iLr1sq / total);
        lastSecondaryRmsCurrent    = std::sqrt(iLr2sq / total);
        lastCr1PeakVoltage         = vCr1pk;
        lastCr2PeakVoltage         = vCr2pk;
    } else {
        lastPrimaryPeakCurrent     = iLr2pk;
        lastSecondaryPeakCurrent   = iLr1pk;
        lastPrimaryRmsCurrent      = std::sqrt(iLr2sq / total);
        lastSecondaryRmsCurrent    = std::sqrt(iLr1sq / total);
        lastCr1PeakVoltage         = vCr2pk;
        lastCr2PeakVoltage         = vCr1pk;
    }
    lastMagnetizingPeakCurrent = iLmpk;

    // Current sharing: max(i_pri) / (max(i_sec) / n_actual). Should be 1.0 ± few %
    // when tankSymmetryRatio == 1.0.
    double n_actual = n;  // physical turns ratio (HV→LV), independent of solver frame
    double iLr2_pk_refl = lastSecondaryPeakCurrent / n_actual;
    lastCurrentSharingRatio = (iLr2_pk_refl > 1e-9)
                              ? lastPrimaryPeakCurrent / iLr2_pk_refl : 0.0;

    // ZVS at switching events (k=0 / k=N, the period boundaries). The lagging
    // bridge soft-switches when the magnetizing-current contribution at
    // switch instant has the same sign as -i_load_reflected (i.e. the
    // resonant-tank current discharges Coss into the rising rail).
    // Our metric is i_Lm magnitude minus the load-reflected current magnitude.
    double iLr1_at_sw = iLr1F[0];
    double iLm_at_sw  = iLmF[0];
    double iLoad_refl_at_sw = iLr1_at_sw - iLm_at_sw;  // = i_Lr2 / n in solver frame
    double zvsActive  = std::abs(iLm_at_sw) - std::abs(iLoad_refl_at_sw);
    double zvsPassive = std::abs(iLr2F[0]);  // SR-side "ZVS-able" current at switch
    if (!reverse) {
        lastZvsMarginPrimaryLagging   = zvsActive;
        lastZvsMarginSecondaryLagging = zvsPassive;
    } else {
        lastZvsMarginPrimaryLagging   = zvsPassive;
        lastZvsMarginSecondaryLagging = zvsActive;
    }
    lastZvsLoadThresholdPrimary   = 0.0;  // not modelled in v1
    lastZvsLoadThresholdSecondary = 0.0;
    lastResonantTransitionTime    = computedDeadTime;

    double fr = computedPrimaryResonantFrequency;
    if (fr > 0) {
        double ratio = fsw / fr;
        if      (ratio > 1.05) lastOperatingRegion = OperatingRegion::ABOVE_RESONANCE;
        else if (ratio < 0.95) lastOperatingRegion = OperatingRegion::BELOW_RESONANCE;
        else                   lastOperatingRegion = OperatingRegion::AT_RESONANCE;
    }
    lastPowerFlowDirection = reverse ? PowerFlowDirection::REVERSE
                                     : PowerFlowDirection::FORWARD;
    lastModeForward = reverse ? 0 : 1;  // mode classification deferred to next phase
    lastModeReverse = reverse ? 1 : 0;

    // ----- Extra-component waveforms (for downstream component sizing). -----
    auto makeWfm = [&](const std::vector<double>& d) {
        Waveform w;
        w.set_data(d);
        w.set_time(tF);
        w.set_ancillary_label(WaveformLabel::CUSTOM);
        return w;
    };
    extraLr1CurrentWaveforms.push_back(makeWfm(iLr1F));
    extraLr2CurrentWaveforms.push_back(makeWfm(iLr2F));
    extraCr1VoltageWaveforms.push_back(makeWfm(vCr1F));
    extraCr2VoltageWaveforms.push_back(makeWfm(vCr2F));
    extraCr1CurrentWaveforms.push_back(makeWfm(iLr1F));   // series = i_Lr1
    extraCr2CurrentWaveforms.push_back(makeWfm(iLr2F));   // series = i_Lr2
    extraLmCurrentWaveforms.push_back(makeWfm(iLmF));
    extraTimeVectors.push_back(tF);

    // ----- Excitations: primary and secondary winding -----
    {
        Waveform iWfm = makeWfm(reverse ? iLr2F : iLr1F);
        Waveform vWfm = makeWfm(reverse ? vSecF : vPriF);
        // In reverse the "primary" winding voltage in the original reference
        // frame is vSec rather than vPri; sign convention here yields the
        // magnitude correctly for downstream core-loss / RMS calculations.
        operatingPoint.get_mutable_excitations_per_winding().push_back(
            complete_excitation(iWfm, vWfm, fsw, "Primary"));
    }
    {
        Waveform iWfm = makeWfm(reverse ? iLr1F : iLr2F);
        Waveform vWfm = makeWfm(reverse ? vPriF : vSecF);
        operatingPoint.get_mutable_excitations_per_winding().push_back(
            complete_excitation(iWfm, vWfm, fsw, "Secondary"));
    }

    OperatingConditions conditions;
    conditions.set_ambient_temperature(op.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);
    return operatingPoint;
}

// =====================================================================
// SPICE — TODO(clllc-v1): direction-aware 8-switch netlist with three
// K-couplings (K1, K2, K_pri_sec) per CLLLC_PLAN.md §A.6.
// =====================================================================
std::string Clllc::generate_ngspice_circuit(
        const std::vector<double>& /*turnsRatios*/,
        double /*magnetizingInductance*/,
        size_t /*inputVoltageIndex*/,
        size_t /*operatingPointIndex*/) {
    throw std::logic_error(
        "Clllc::generate_ngspice_circuit not yet implemented. "
        "Pending direction-aware 8-switch netlist (CLLLC_PLAN.md §A.6).");
}

std::vector<OperatingPoint> Clllc::simulate_and_extract_operating_points(
        const std::vector<double>& /*turnsRatios*/,
        double /*magnetizingInductance*/) {
    throw std::logic_error(
        "Clllc::simulate_and_extract_operating_points not yet implemented. "
        "Depends on generate_ngspice_circuit (CLLLC_PLAN.md §A.7).");
}

std::vector<ConverterWaveforms> Clllc::simulate_and_extract_topology_waveforms(
        const std::vector<double>& /*turnsRatios*/,
        double /*magnetizingInductance*/,
        size_t /*numberOfPeriods*/) {
    throw std::logic_error(
        "Clllc::simulate_and_extract_topology_waveforms not yet implemented. "
        "Depends on generate_ngspice_circuit (CLLLC_PLAN.md §A.8).");
}

// =====================================================================
// Extra components — TODO(clllc-v1) per CLLLC_PLAN.md §A.9 (transformer
// with two integrated leakages OR two discrete inductors, plus Cr1, Cr2,
// HV bus cap, LV bus cap).
// =====================================================================
std::vector<std::variant<Inputs, CAS::Inputs>> Clllc::get_extra_components_inputs(
        ExtraComponentsMode /*mode*/,
        std::optional<Magnetic> /*magnetic*/) {
    throw std::logic_error(
        "Clllc::get_extra_components_inputs not yet implemented. "
        "Pending solver (needs operating-point waveforms to size Cr1/Cr2/bus caps). "
        "See CLLLC_PLAN.md §A.9.");
}

// =====================================================================
// AdvancedClllc::process — TODO(clllc-v1)
// =====================================================================
Inputs AdvancedClllc::process() {
    throw std::logic_error(
        "AdvancedClllc::process not yet implemented. "
        "Depends on Clllc::process_operating_points.");
}

} // namespace OpenMagnetics
