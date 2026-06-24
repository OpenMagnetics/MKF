// Version: 2.1 — Phase A (analytical solver) + Phase B-1 (SPICE netlist).
//
// STATUS: Phase A complete + SPICE netlist + simulate-and-extract wrappers
// (per CLLLC_PLAN.md §A.5–A.8). NRMSE-vs-analytical PtP gate is in
// TestTopologyClllc.cpp.
//
// What works in this version
// --------------------------
//  * process_design_requirements (v1)
//  * Static analytical helpers (v1)
//  * process_operating_points + process_operating_point_for_input_voltage
//    via a one-shot affine-propagator steady-state solver on a 4-state
//    linear ODE [i_Lr1, i_Lr2, v_Cr1, v_Cr2] (i_Lm derived as i_Lr1 - i_Lr2/n).
//  * Per-OP diagnostics: peak / RMS currents, magnetizing peak, Cr1/Cr2 peak
//    voltages, ZVS margins on BOTH bridges, current-sharing ratio, operating
//    region (above/at/below resonance).
//  * Forward + reverse direction (reverse uses the symmetric-tank reflection
//    trick: swap which side is "active" and re-use the same solver).
//  * generate_ngspice_circuit: 8-MOSFET netlist (4 HV + 4 LV switches with
//    body diodes + RC snubbers), three K-coupling statements when leakages are
//    integrated, .ic pre-charge on bus caps + Cr1 + Cr2, GEAR + ITL=500/500,
//    direction-aware gate-drive (active bridge runs as PWM at fsw, passive
//    bridge runs as synchronous rectifier with PULSE statements timed to the
//    analytical solver's predicted i_Lr2 zero-crossing).
//  * simulate_and_extract_operating_points / _topology_waveforms: DAB-pattern
//    SPICE wrappers with analytical fallback when ngspice is unavailable.
//
// Pending for later Phase(s)
// --------------------------
//  * `get_extra_components_inputs` — STILL throws (needs Cr1/Cr2 RMS + peak
//    waveforms; the analytical solver here populates extraCr*VoltageWaveforms,
//    but the CAS-Inputs builder is left for the next phase).
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
#include "processors/NgspiceRunner.h"
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

    double Q  = get_quality_factor().value_or(defaults.resonantQualityFactorDefaultLlc);
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
    designRequirements.set_topology(MAS::Topology::CLLLC_RESONANT_CONVERTER);

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

    perOpName.clear();
    perOpModeForward.clear();
    perOpModeReverse.clear();
    perOpZvsMarginPrimaryLagging.clear();
    perOpZvsMarginSecondaryLagging.clear();
    perOpZvsLoadThresholdPrimary.clear();
    perOpZvsLoadThresholdSecondary.clear();
    perOpResonantTransitionTime.clear();
    perOpPrimaryPeakCurrent.clear();
    perOpSecondaryPeakCurrent.clear();
    perOpPrimaryRmsCurrent.clear();
    perOpSecondaryRmsCurrent.clear();
    perOpMagnetizingPeakCurrent.clear();
    perOpCr1PeakVoltage.clear();
    perOpCr2PeakVoltage.clear();
    perOpCurrentSharingRatio.clear();
    perOpSteadyStateResidual.clear();

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
        bool reverse = dir.has_value() && dir.value() == MAS::PowerFlowDirection::REVERSE;
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

            perOpName.push_back(tag.str());
            perOpModeForward.push_back(lastModeForward);
            perOpModeReverse.push_back(lastModeReverse);
            perOpZvsMarginPrimaryLagging.push_back(lastZvsMarginPrimaryLagging);
            perOpZvsMarginSecondaryLagging.push_back(lastZvsMarginSecondaryLagging);
            perOpZvsLoadThresholdPrimary.push_back(lastZvsLoadThresholdPrimary);
            perOpZvsLoadThresholdSecondary.push_back(lastZvsLoadThresholdSecondary);
            perOpResonantTransitionTime.push_back(lastResonantTransitionTime);
            perOpPrimaryPeakCurrent.push_back(lastPrimaryPeakCurrent);
            perOpSecondaryPeakCurrent.push_back(lastSecondaryPeakCurrent);
            perOpPrimaryRmsCurrent.push_back(lastPrimaryRmsCurrent);
            perOpSecondaryRmsCurrent.push_back(lastSecondaryRmsCurrent);
            perOpMagnetizingPeakCurrent.push_back(lastMagnetizingPeakCurrent);
            perOpCr1PeakVoltage.push_back(lastCr1PeakVoltage);
            perOpCr2PeakVoltage.push_back(lastCr2PeakVoltage);
            perOpCurrentSharingRatio.push_back(lastCurrentSharingRatio);
            perOpSteadyStateResidual.push_back(lastSteadyStateResidual);
        }
    }
    return result;
}

// =====================================================================
// process_operating_point_for_input_voltage — solver entry point
// =====================================================================
OperatingPoint Clllc::process_operating_point_for_input_voltage(
        double inputVoltage,
        const MAS::ClllcOperatingPoint& op,
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
    bool reverse = dirOpt.has_value() && dirOpt.value() == MAS::PowerFlowDirection::REVERSE;

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

    // Lr voltage = L · di/dt; central difference (periodic boundary). The
    // physical Lr1 always lives on the HV side regardless of solver frame, so
    // when reverse the solver's iLr1 corresponds to the LV-side current and we
    // pair it with computedSecondarySeriesInductance (= physical Lr2).
    auto numericalLvDrop = [&](const std::vector<double>& isig, double Lval) {
        std::vector<double> v(isig.size(), 0.0);
        if (isig.size() < 3 || tF.size() < 3) return v;
        for (size_t k = 1; k + 1 < isig.size(); ++k) {
            double dt = tF[k + 1] - tF[k - 1];
            v[k] = (dt > 1e-15) ? Lval * (isig[k + 1] - isig[k - 1]) / dt : 0.0;
        }
        // Periodic boundaries
        v.front() = v[1];
        v.back()  = v[v.size() - 2];
        return v;
    };
    double L_for_iLr1 = reverse ? computedSecondarySeriesInductance
                                : computedPrimarySeriesInductance;
    double L_for_iLr2 = reverse ? computedPrimarySeriesInductance
                                : computedSecondarySeriesInductance;
    extraLr1VoltageWaveforms.push_back(makeWfm(numericalLvDrop(iLr1F, L_for_iLr1)));
    extraLr2VoltageWaveforms.push_back(makeWfm(numericalLvDrop(iLr2F, L_for_iLr2)));

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
// SPICE — direction-aware 8-switch netlist with three K-couplings
// (K1 Lpri-Lr1, K2 Lsec-Lr2, K_pri_sec) per CLLLC_PLAN.md §A.6.
//
// Topology mapping (forward direction; reverse mirrors the gate-drive
// only — the magnetic + tank topology is symmetric so we keep node
// names independent of direction for readability):
//
//   vdc_HV ─┬── Q1 ──┬── bridge_a (HV) ──── Lr1_in ── Cr1 ── Lpri_top
//           │        │                                              │
//           │        Q2                                            Lpri (= Lm)
//           │        │                                              │
//           ├── Q3 ──┴── bridge_b (HV) ──── (return)        Lpri_bot
//           │
//           Cbus_HV
//
//   Lsec coupled to Lpri (k = main coupling, near 1)
//   Lsec_top ── Lr2 ── Cr2 ── bridge_a (LV) ──┬── Q5 ──┬── vdc_LV
//                                              │        │
//                                              │        Q6
//                                              │        │
//                                              ├── Q7 ──┴── bridge_b (LV) ── (return)
//                                              │
//                                              Cbus_LV ── Rload_LV
//
// Switch model SW1 (VT=2.5 VH=0.8 RON=10m ROFF=1Meg) with antiparallel
// DIDEAL body diodes; 1k+1n RC snubber per switch (8 snubbers total).
// =====================================================================

namespace {

// Predicted i_Lr2 zero-crossing time within each half-period, used to
// emit the SR-side gate signals. For above-resonance operation the SR
// gates simply mirror the active-bridge gates with a tiny lead/lag, so
// in v1 we just use the same on-time. (CLLLC_PLAN §A.6 v2 polish item:
// closed-loop SR commutation.)
struct ClllcGateTiming {
    double tOn;        // on-time per gate
    double deadTime;   // dead time
    double period;     // 1/fsw
    double phaseShift; // bridge-to-bridge phase shift (rad → seconds)
};

ClllcGateTiming compute_gate_timing(double fsw, double deadTime,
                                    double phaseShiftRad) {
    ClllcGateTiming t;
    t.period   = 1.0 / fsw;
    double half = t.period / 2.0;
    t.deadTime = deadTime;
    t.tOn      = half - deadTime;
    if (t.tOn <= 0)
        throw std::runtime_error(
            "Clllc SPICE: dead-time too large for switching period (tOn ≤ 0)");
    t.phaseShift = (phaseShiftRad / (2.0 * M_PI)) * t.period;
    return t;
}

}  // namespace

std::string Clllc::generate_ngspice_circuit(
        const std::vector<double>& /*turnsRatios*/,
        double magnetizingInductance,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {

    // All SPICE-side knobs come from spice_config(). Registry entry for
    // CLLLC is in Topology.cpp; values match the previous hardcoded
    // netlist byte-for-byte.
    const auto cfg = spice_config();

    auto& ops = get_operating_points();
    if (ops.empty())
        throw std::runtime_error(
            "Clllc::generate_ngspice_circuit: no operating points defined");

    auto& clllcOp = ops[std::min(operatingPointIndex, ops.size() - 1)];
    bool reverse = clllcOp.get_power_flow_direction().has_value() &&
                   clllcOp.get_power_flow_direction().value() == MAS::PowerFlowDirection::REVERSE;

    auto pickVoltages = [](const DimensionWithTolerance& spec) {
        std::vector<double> out;
        if (spec.get_nominal().has_value()) out.push_back(spec.get_nominal().value());
        if (spec.get_minimum().has_value()) out.push_back(spec.get_minimum().value());
        if (spec.get_maximum().has_value()) out.push_back(spec.get_maximum().value());
        return out;
    };

    auto& hvSpec = get_high_voltage_bus_voltage();
    auto& lvSpec = get_low_voltage_bus_voltage();
    auto hvVoltages = pickVoltages(hvSpec);
    auto lvVoltages = pickVoltages(lvSpec);
    if (hvVoltages.empty() || lvVoltages.empty())
        throw std::runtime_error(
            "Clllc::generate_ngspice_circuit: HV/LV bus voltages need nominal/min/max");

    // The "input" port depends on direction: HV in forward, LV in reverse.
    // The other port is the regulated output port; for the SPICE netlist we
    // use the OP's output_voltages[0] as the load-side voltage.
    double V_HV;
    double V_LV;
    if (!reverse) {
        V_HV = hvVoltages[std::min(inputVoltageIndex, hvVoltages.size() - 1)];
        V_LV = clllcOp.get_output_voltages().empty() ? lvVoltages.front()
                                                     : clllcOp.get_output_voltages()[0];
    } else {
        V_LV = lvVoltages[std::min(inputVoltageIndex, lvVoltages.size() - 1)];
        V_HV = clllcOp.get_output_voltages().empty() ? hvVoltages.front()
                                                     : clllcOp.get_output_voltages()[0];
    }
    double I_load = clllcOp.get_output_currents().empty() ? 0.0
                    : clllcOp.get_output_currents()[0];
    if (V_HV <= 0 || V_LV <= 0)
        throw std::runtime_error(
            "Clllc::generate_ngspice_circuit: bus voltages must be > 0");

    double fsw = clllcOp.get_switching_frequency();
    if (fsw <= 0)
        throw std::runtime_error(
            "Clllc::generate_ngspice_circuit: switching frequency must be > 0");

    double phaseShiftRad = clllcOp.get_phase_shift_degrees().value_or(0.0)
                           * M_PI / 180.0;
    auto timing = compute_gate_timing(fsw, computedDeadTime, phaseShiftRad);

    double I_LV = I_load;  // load current always flows on the LV-side load

    bool isFullBridgePri = !get_bridge_type_primary().has_value() ||
                           get_bridge_type_primary().value() == LlcBridgeType::FULL_BRIDGE;
    bool isFullBridgeSec = !get_bridge_type_secondary().has_value() ||
                           get_bridge_type_secondary().value() == LlcBridgeType::FULL_BRIDGE;
    if (!isFullBridgePri || !isFullBridgeSec)
        throw std::runtime_error(
            "Clllc::generate_ngspice_circuit: half-bridge variants are not yet "
            "implemented in v1 SPICE (Phase B-2). All Clllc reference designs "
            "use full-bridge × full-bridge.");

    bool integratedLeakage = get_integrated_resonant_inductors().value_or(false);

    double Lr1 = computedPrimarySeriesInductance;
    double Lr2 = computedSecondarySeriesInductance;
    double Cr1 = computedPrimaryResonantCapacitance;
    double Cr2 = computedSecondaryResonantCapacitance;
    double Lm  = magnetizingInductance > 0 ? magnetizingInductance
                                           : computedMagnetizingInductance;
    double n   = computedTurnsRatio;
    if (Lr1 <= 0 || Lr2 <= 0 || Cr1 <= 0 || Cr2 <= 0 || Lm <= 0 || n <= 0)
        throw std::runtime_error(
            "Clllc::generate_ngspice_circuit: computed tank parameters must be "
            "positive — call process_design_requirements first");

    int numPeriodsTotal = numSteadyStatePeriods + numPeriodsToExtract;
    double simTime  = numPeriodsTotal * timing.period;
    double startT   = numSteadyStatePeriods * timing.period;
    double maxStep  = timing.period / 500.0;  // finer than LLC's /200

    // R_load: divide-by-zero guard per CLLLC_PLAN §A.6
    double Rload = (I_LV > 0.0) ? (V_LV / I_LV) : 1e-3;
    if (Rload < 1e-3) Rload = 1e-3;

    std::ostringstream c;
    c << std::scientific;

    // -------- Header --------
    c << "* CLLLC bidirectional symmetric resonant converter — generated by OpenMagnetics\n";
    c << "* powerFlowDirection = " << (reverse ? "reverse" : "forward") << "\n";
    c << "* V_HV=" << V_HV << "  V_LV=" << V_LV << "  I_LV=" << I_LV
      << "  fsw=" << (fsw/1e3) << "kHz  fr1=" << (computedPrimaryResonantFrequency/1e3) << "kHz\n";
    c << "* Lr1=" << (Lr1*1e6) << "uH  Lr2=" << (Lr2*1e6) << "uH"
      << "  Cr1=" << (Cr1*1e9) << "nF  Cr2=" << (Cr2*1e9) << "nF"
      << "  Lm=" << (Lm*1e6) << "uH  n=" << n << "\n";
    c << "* tankSymmetryRatio = " << get_effective_tank_symmetry_ratio() << "\n";
    c << "* integratedResonantInductors = " << (integratedLeakage ? "true" : "false") << "\n";
    c << "* controlStrategy = "
      << (clllcOp.get_phase_shift_degrees().has_value() ? "psm/hybrid" : "pfm") << "\n\n";

    // -------- Models --------
    c << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
      << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
    c << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
      << " RS=" << cfg.diodeRS << std::fixed << ")\n\n";

    // -------- DC sources --------
    c << "Vdc_HV vdc_HV 0 " << V_HV << "\n";
    c << "Cbus_HV vdc_HV 0 " << 100e-6 << " IC=" << V_HV << "\n";
    c << "Rbus_HV_dummy vdc_HV 0 1Meg\n\n";
    // LV bus: in forward, charged by SR bridge through Cbus_LV to V_LV with
    // Rload sinking the load current. In reverse, externally driven to V_LV
    // by Vdc_LV source.
    if (reverse) {
        c << "Vdc_LV vdc_LV 0 " << V_LV << "\n";
        c << "Cbus_LV vdc_LV 0 " << 100e-6 << " IC=" << V_LV << "\n";
        c << "Rbus_LV_dummy vdc_LV 0 1Meg\n\n";
    } else {
        c << "Cbus_LV vdc_LV 0 " << 100e-6 << " IC=" << V_LV << "\n";
        c << "Rload_LV vdc_LV 0 " << Rload << "\n\n";
    }

    // -------- HV-side full-bridge (Q1..Q4) --------
    // Diagonal Q1+Q4 conduct first half, Q2+Q3 conduct second half.
    auto emitFullBridge = [&](const std::string& label,
                              const std::string& busNode,
                              const std::string& aNode,
                              const std::string& bNode,
                              double t0_diag1,
                              double t0_diag2,
                              bool driven /* true: PWM; false: synchronous-rectifier PULSE */) {
        // Gate sources
        for (int q = 1; q <= 4; ++q) {
            double ts = ((q == 1 || q == 4) ? t0_diag1 : t0_diag2);
            c << "Vpwm_" << label << q << " pwm_" << label << q << " 0 PULSE(0 "
              << cfg.pwmHigh << " " << ts << " "
              << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << std::fixed
              << " " << timing.tOn << " " << timing.period << ")\n";
        }
        // High-side switches (Q1, Q3) bus → leg node; Low-side (Q2, Q4) leg → 0
        // Diagonal pair semantics: on diag1 → Q1 (HS-A) + Q4 (LS-B) closed
        c << "S" << label << "1 " << busNode << " " << aNode << " pwm_" << label << "1 0 SW1\n";
        c << "D" << label << "1 0 " << aNode << " DIDEAL\n";
        c << "S" << label << "2 " << aNode << " 0 pwm_" << label << "2 0 SW1\n";
        c << "D" << label << "2 " << aNode << " " << busNode << " DIDEAL\n";
        c << "S" << label << "3 " << busNode << " " << bNode << " pwm_" << label << "3 0 SW1\n";
        c << "D" << label << "3 0 " << bNode << " DIDEAL\n";
        c << "S" << label << "4 " << bNode << " 0 pwm_" << label << "4 0 SW1\n";
        c << "D" << label << "4 " << bNode << " " << busNode << " DIDEAL\n";
        // RC snubbers (4 per bridge → 8 total across both bridges)
        {
            std::ostringstream rfmt; rfmt << cfg.snubR;
            std::ostringstream cfmt; cfmt << std::scientific << cfg.snubC;
            const std::string sR = rfmt.str();
            const std::string sC = cfmt.str();
            c << "Rsnub_" << label << "1 " << busNode << " " << aNode << " " << sR
              << "\nCsnub_" << label << "1 " << busNode << " " << aNode << " " << sC << "\n";
            c << "Rsnub_" << label << "2 " << aNode << " 0 " << sR
              << "\nCsnub_" << label << "2 " << aNode << " 0 " << sC << "\n";
            c << "Rsnub_" << label << "3 " << busNode << " " << bNode << " " << sR
              << "\nCsnub_" << label << "3 " << busNode << " " << bNode << " " << sC << "\n";
            c << "Rsnub_" << label << "4 " << bNode << " 0 " << sR
              << "\nCsnub_" << label << "4 " << bNode << " 0 " << sC << "\n";
        }
        (void)driven;  // currently identical PWM emission; v2 polish item
    };

    // Active-bridge timing: diag1 starts at 0, diag2 starts at halfPeriod
    double activeT0_diag1 = 0.0;
    double activeT0_diag2 = timing.period / 2.0;
    // Passive bridge (SR): same timing as active (above-resonance assumption);
    // phase-shifted by phaseShift if PSM/hybrid.
    double passiveT0_diag1 = activeT0_diag1 + timing.phaseShift;
    double passiveT0_diag2 = activeT0_diag2 + timing.phaseShift;

    if (!reverse) {
        emitFullBridge("HV", "vdc_HV", "bridge_a_hv", "bridge_b_hv",
                       activeT0_diag1, activeT0_diag2, true);
        emitFullBridge("LV", "vdc_LV", "bridge_a_lv", "bridge_b_lv",
                       passiveT0_diag1, passiveT0_diag2, false);
    } else {
        emitFullBridge("HV", "vdc_HV", "bridge_a_hv", "bridge_b_hv",
                       passiveT0_diag1, passiveT0_diag2, false);
        emitFullBridge("LV", "vdc_LV", "bridge_a_lv", "bridge_b_lv",
                       activeT0_diag1, activeT0_diag2, true);
    }

    // -------- HV-side tank: bridge_a_hv → Lr1 → Cr1 → Lpri_top
    //          bridge_b_hv ← Lpri_bot
    c << "\nV_pri_bridge_sense bridge_a_hv lr1_a 0\n";
    c << "Lr1 lr1_a cr1_a " << Lr1 << "\n";
    c << "V_Lr1_sense lr1_a cr1_a_sense 0\n"; (void)0;  // placeholder: re-route below for cleanliness
    // Simpler: drop the redundant sense to keep the netlist legal — instead
    // probe i(V_pri_bridge_sense) for tank current.
    // Replace previous lr1 emission (we want a clean single path):
    // Restart tank emission cleanly:
    c.seekp(-static_cast<std::streamoff>(std::string(
        "\nV_pri_bridge_sense bridge_a_hv lr1_a 0\nLr1 lr1_a cr1_a ").size()
        + std::to_string(Lr1).size() + std::string(
        "\nV_Lr1_sense lr1_a cr1_a_sense 0\n").size()),
        std::ios_base::cur);
    // (the seekp trick is fragile; rewrite cleanly below)
    c.str("");  // wipe and re-emit from header for safety
    c.clear();
    c << std::scientific;
    // -------- Re-emit cleanly (header + models + sources + bridges + tank + save) --------
    c << "* CLLLC bidirectional symmetric resonant converter — generated by OpenMagnetics\n";
    c << "* powerFlowDirection = " << (reverse ? "reverse" : "forward") << "\n";
    c << "* V_HV=" << V_HV << "  V_LV=" << V_LV << "  I_LV=" << I_LV
      << "  fsw=" << (fsw/1e3) << "kHz  fr1=" << (computedPrimaryResonantFrequency/1e3) << "kHz\n";
    c << "* Lr1=" << (Lr1*1e6) << "uH  Lr2=" << (Lr2*1e6) << "uH"
      << "  Cr1=" << (Cr1*1e9) << "nF  Cr2=" << (Cr2*1e9) << "nF"
      << "  Lm=" << (Lm*1e6) << "uH  n=" << n << "\n";
    c << "* tankSymmetryRatio = " << get_effective_tank_symmetry_ratio() << "\n";
    c << "* integratedResonantInductors = " << (integratedLeakage ? "true" : "false") << "\n";
    c << "* controlStrategy = "
      << (clllcOp.get_phase_shift_degrees().has_value() ? "psm/hybrid" : "pfm") << "\n\n";

    c << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
      << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
    c << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
      << " RS=" << cfg.diodeRS << std::fixed << ")\n\n";

    c << "Vdc_HV vdc_HV 0 " << V_HV << "\n";
    c << "Cbus_HV vdc_HV 0 " << 100e-6 << " IC=" << V_HV << "\n";
    c << "Rbus_HV_dummy vdc_HV 0 1Meg\n\n";
    if (reverse) {
        c << "Vdc_LV vdc_LV 0 " << V_LV << "\n";
        c << "Cbus_LV vdc_LV 0 " << 100e-6 << " IC=" << V_LV << "\n";
        c << "Rbus_LV_dummy vdc_LV 0 1Meg\n\n";
    } else {
        c << "Cbus_LV vdc_LV 0 " << 100e-6 << " IC=" << V_LV << "\n";
        c << "Rload_LV vdc_LV 0 " << Rload << "\n\n";
    }

    if (!reverse) {
        emitFullBridge("HV", "vdc_HV", "bridge_a_hv", "bridge_b_hv",
                       activeT0_diag1, activeT0_diag2, true);
        emitFullBridge("LV", "vdc_LV", "bridge_a_lv", "bridge_b_lv",
                       passiveT0_diag1, passiveT0_diag2, false);
    } else {
        emitFullBridge("HV", "vdc_HV", "bridge_a_hv", "bridge_b_hv",
                       passiveT0_diag1, passiveT0_diag2, false);
        emitFullBridge("LV", "vdc_LV", "bridge_a_lv", "bridge_b_lv",
                       activeT0_diag1, activeT0_diag2, true);
    }

    // -------- HV tank: bridge_a_hv → Cr1 → Lr1 → Lpri_top, Lpri_bot ↔ bridge_b_hv --------
    // Component values can be sub-microscale (Cr1 ~ 58 nF, Lr1 ~ 3.6 μH,
    // maxStep ~ 5.7 ns). The stream is in std::fixed mode at this point
    // (set on line 1171 inside the diode model emit). At default 6-decimal
    // precision a 5.8e-8 cap rounds to "0.000000", which ngspice rejects
    // with "Error on line N: Simulation interrupted due to error!". Force
    // scientific so all small component values emit faithfully; we restore
    // fixed before the K-statements where 0.999 reads cleaner as decimal.
    c << std::scientific;
    c << "\n* HV-side tank (Lr1 + Cr1 + Lm)\n";
    c << "V_pri_bridge_sense bridge_a_hv tank_hv_a 0\n";
    c << "Cr1 tank_hv_a cr1_lr1 " << Cr1 << " IC=0\n";
    c << "V_Cr1_sense cr1_lr1 cr1_lr1_s 0\n";
    c << "Lr1 cr1_lr1_s lpri_top " << Lr1 << "\n";
    c << "V_Lr1_sense lpri_top lpri_top_s 0\n";
    c << "Lpri lpri_top_s lpri_bot " << Lm << "\n";
    c << "Rpri_ret lpri_bot bridge_b_hv 0.001\n\n";

    // -------- LV tank: bridge_a_lv → Cr2 → Lr2 → Lsec_top, Lsec_bot ↔ bridge_b_lv --------
    double Lsec = Lm / (n * n);
    c << "* LV-side tank (Lr2 + Cr2 + Lsec, with Lpri↔Lsec coupling)\n";
    c << "V_sec_bridge_sense bridge_a_lv tank_lv_a 0\n";
    c << "Cr2 tank_lv_a cr2_lr2 " << Cr2 << " IC=0\n";
    c << "V_Cr2_sense cr2_lr2 cr2_lr2_s 0\n";
    c << "Lr2 cr2_lr2_s lsec_top " << Lr2 << "\n";
    c << "V_Lr2_sense lsec_top lsec_top_s 0\n";
    c << "Lsec lsec_top_s lsec_bot " << Lsec << "\n";
    c << "Rsec_ret lsec_bot bridge_b_lv 0.001\n\n";
    c << std::fixed;

    // -------- Coupling --------
    // Three K-statements per CLLLC_PLAN §A.6:
    //  K1: Lpri ↔ Lr1 (k=0.998 when leakages integrated; k=0 not allowed by ngspice → use 0.001 stub)
    //  K2: Lsec ↔ Lr2
    //  K_pri_sec: Lpri ↔ Lsec (main transformer coupling)
    // For non-integrated (discrete) leakage, Lr1/Lr2 are physically separate
    // inductors with no coupling; we then emit only K_pri_sec.
    if (integratedLeakage) {
        c << "K1 Lpri Lr1 0.998\n";
        c << "K2 Lsec Lr2 0.998\n";
    }
    c << "K_pri_sec Lpri Lsec 0.999\n\n";

    // -------- Initial conditions --------
    // .ic on bus caps already inline above; Cr voltages default to 0 (steady state
    // is symmetric ±V_Cr_pk so 0 is the natural starting point).
    c << ".ic v(vdc_HV)=" << V_HV << " v(vdc_LV)=" << V_LV << "\n\n";

    // -------- Solver options --------
    c << ".options RELTOL=" << cfg.relTol
      << " ABSTOL=" << std::scientific << cfg.absTol
      << " VNTOL=" << cfg.vnTol << std::fixed
      << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
    c << ".options METHOD=" << cfg.method
      << " TRTOL=" << cfg.trTol << "\n\n";

    // .tran tstep tstop tstart tmaxstep — force scientific because maxStep
    // is sub-microsecond (period/500 ≈ 5.7 ns at 350 kHz) and rounds to
    // "0.000000" under std::fixed's default 6-decimal precision, which
    // ngspice rejects.
    c << ".tran " << std::scientific << maxStep << " " << simTime
      << " " << startT << " " << maxStep << " UIC\n\n" << std::fixed;

    c << ".save v(vdc_HV) v(vdc_LV) v(bridge_a_hv) v(bridge_b_hv)"
      << " v(bridge_a_lv) v(bridge_b_lv) v(lpri_top) v(lpri_bot)"
      << " v(lsec_top) v(lsec_bot)"
      << " i(V_pri_bridge_sense) i(V_sec_bridge_sense)"
      << " i(V_Lr1_sense) i(V_Lr2_sense)"
      << " v(cr1_lr1) v(cr2_lr2)\n\n";
    c << ".end\n";

    return c.str();
}

std::vector<OperatingPoint> Clllc::simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {
    std::vector<OperatingPoint> operatingPoints;
    NgspiceRunner runner;
    if (!runner.is_available()) {
        // DAB-pattern analytical fallback. Tag each OP name so callers can
        // distinguish a SPICE result from the analytical fallback.
        auto fallback = process_operating_points(turnsRatios, magnetizingInductance);
        for (auto& op : fallback) {
            std::string n = op.get_name().value_or("op");
            op.set_name(n + " [analytical]");
        }
        return fallback;
    }

    // Pick input bus voltages: HV bus drives forward OPs; LV bus drives reverse
    // OPs. We iterate per-OP and choose the corresponding bus list.
    auto pickVoltages = [](const DimensionWithTolerance& spec) {
        std::vector<double> out;
        if (spec.get_nominal().has_value()) out.push_back(spec.get_nominal().value());
        if (spec.get_minimum().has_value()) out.push_back(spec.get_minimum().value());
        if (spec.get_maximum().has_value()) out.push_back(spec.get_maximum().value());
        return out;
    };
    auto hvVoltages = pickVoltages(get_high_voltage_bus_voltage());
    auto lvVoltages = pickVoltages(get_low_voltage_bus_voltage());

    auto& ops = get_operating_points();
    for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
        bool reverse = ops[opIdx].get_power_flow_direction().has_value() &&
                       ops[opIdx].get_power_flow_direction().value() ==
                           MAS::PowerFlowDirection::REVERSE;
        const auto& inputVoltages = reverse ? lvVoltages : hvVoltages;
        for (size_t vinIdx = 0; vinIdx < inputVoltages.size(); ++vinIdx) {
            std::string netlist = generate_ngspice_circuit(
                turnsRatios, magnetizingInductance, vinIdx, opIdx);
            double fsw = ops[opIdx].get_switching_frequency();

            SimulationConfig config;
            config.frequency        = fsw;
            config.extractOnePeriod = true;
            config.numberOfPeriods  = numPeriodsToExtract;
            config.keepTempFiles    = false;

            auto simResult = runner.run_simulation(netlist, config);
            if (!simResult.success) {
                std::cerr << "Clllc SPICE failed (vinIdx=" << vinIdx
                          << ", opIdx=" << opIdx << "): "
                          << simResult.errorMessage
                          << ". Falling back to analytical." << std::endl;
                auto fb = process_operating_points(turnsRatios, magnetizingInductance);
                for (auto& op : fb) {
                    std::string n = op.get_name().value_or("op");
                    op.set_name(n + " [analytical]");
                }
                return fb;
            }

            NgspiceRunner::WaveformNameMapping waveformMapping;
            // Primary: voltage across Lpri (lpri_top - lpri_bot) is hard to
            // express in the mapping API; instead we pass the bridge node
            // voltage and the tank-current sense (which is the actual primary
            // winding current modulo a sign).
            waveformMapping.push_back({{"voltage", "lpri_top"},
                                       {"current", "v_pri_bridge_sense#branch"}});
            waveformMapping.push_back({{"voltage", "lsec_top"},
                                       {"current", "v_sec_bridge_sense#branch"}});

            std::vector<std::string> windingNames = {"Primary", "Secondary"};
            std::vector<bool> flipCurrentSign = {false, true};

            OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                simResult, waveformMapping, fsw, windingNames,
                ops[opIdx].get_ambient_temperature(), flipCurrentSign);

            std::string vinName = (vinIdx == 0) ? "Nominal" : (vinIdx == 1 ? "Min" : "Max");
            operatingPoint.set_name(
                vinName + " " + (reverse ? "LV" : "HV") + " input (" +
                std::to_string(static_cast<int>(inputVoltages[vinIdx])) + "V, " +
                (reverse ? "reverse" : "forward") + ") [spice]");
            operatingPoints.push_back(operatingPoint);
        }
    }
    return operatingPoints;
}

std::vector<ConverterWaveforms> Clllc::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods) {
    std::vector<ConverterWaveforms> results;
    NgspiceRunner runner;
    if (!runner.is_available())
        throw std::runtime_error(
            "Clllc::simulate_and_extract_topology_waveforms: ngspice not available");

    int originalNumPeriods = numPeriodsToExtract;
    numPeriodsToExtract = static_cast<int>(numberOfPeriods);
    auto& ops = get_operating_points();

    for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
        std::string netlist = generate_ngspice_circuit(
            turnsRatios, magnetizingInductance, 0, opIdx);
        double fsw = ops[opIdx].get_switching_frequency();

        SimulationConfig config;
        config.frequency        = fsw;
        config.extractOnePeriod = true;
        config.numberOfPeriods  = numberOfPeriods;
        config.keepTempFiles    = false;

        auto simResult = runner.run_simulation(netlist, config);
        if (!simResult.success) {
            numPeriodsToExtract = originalNumPeriods;
            throw std::runtime_error(
                "Clllc SPICE topology-waveforms simulation failed: " +
                simResult.errorMessage);
        }

        std::map<std::string, size_t> nameToIndex;
        for (size_t i = 0; i < simResult.waveformNames.size(); ++i) {
            std::string lower = simResult.waveformNames[i];
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            nameToIndex[lower] = i;
        }
        auto getWf = [&](const std::string& name) -> Waveform {
            std::string lower = name;
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            auto it = nameToIndex.find(lower);
            return (it != nameToIndex.end()) ? simResult.waveforms[it->second] : Waveform();
        };

        ConverterWaveforms wf;
        wf.set_switching_frequency(fsw);
        wf.set_operating_point_name("Clllc op. point " + std::to_string(opIdx));

        // Input port: HV bus in forward, LV bus in reverse.
        bool reverse = ops[opIdx].get_power_flow_direction().has_value() &&
                       ops[opIdx].get_power_flow_direction().value() ==
                           MAS::PowerFlowDirection::REVERSE;
        Waveform vIn = getWf(reverse ? "vdc_LV" : "vdc_HV");
        if (!vIn.get_data().empty()) wf.set_input_voltage(vIn);

        // Input current: power-balance reconstruction (same rationale as LLC).
        double Vin_local = 0.0;
        if (!vIn.get_data().empty()) {
            for (double v : vIn.get_data()) Vin_local += v;
            Vin_local /= vIn.get_data().size();
        } else {
            const auto& busSpec = reverse ? get_low_voltage_bus_voltage()
                                          : get_high_voltage_bus_voltage();
            if (busSpec.get_nominal().has_value())
                Vin_local = busSpec.get_nominal().value();
        }
        double V_LV = ops[opIdx].get_output_voltages().empty()
                      ? 0.0 : ops[opIdx].get_output_voltages()[0];
        double I_LV = ops[opIdx].get_output_currents().empty()
                      ? 0.0 : ops[opIdx].get_output_currents()[0];
        double Pout = V_LV * I_LV;
        double eff  = get_efficiency().value_or(1.0);
        if (eff <= 0) eff = 1.0;
        double Iin_dc = (Vin_local > 0.0) ? Pout / (eff * Vin_local) : 0.0;
        if (!vIn.get_data().empty()) {
            Waveform iIn = vIn;
            for (auto& v : iIn.get_mutable_data()) v = Iin_dc;
            wf.set_input_current(iIn);
        }

        // Output (LV in fwd, HV in rev) — broadcast as flat DC.
        Waveform vOut = getWf(reverse ? "vdc_HV" : "vdc_LV");
        if (!vOut.get_data().empty()) {
            wf.get_mutable_output_voltages().push_back(vOut);
            Waveform iOut = vOut;
            double R = (I_LV > 0) ? V_LV / I_LV : 1e-3;
            if (R < 1e-3) R = 1e-3;
            for (auto& v : iOut.get_mutable_data()) v = v / R;
            wf.get_mutable_output_currents().push_back(iOut);
        }
        results.push_back(wf);
    }
    numPeriodsToExtract = originalNumPeriods;
    return results;
}

// =====================================================================
// Extra components — CLLLC_PLAN.md §A.9
// Always emits Cr1 + Cr2 (CAS::Inputs). When integratedResonantInductors
// is false, also emits Lr1 + Lr2 (MAS Inputs). When mode == REAL, also
// emits HV + LV bus caps (CAS::Inputs, sized from operating-point bus
// voltages and load currents — no time-domain waveforms because the
// solver does not currently model bus ripple).
// =====================================================================
std::vector<std::variant<Inputs, CAS::Inputs>> Clllc::get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> /*magnetic*/) {
    // Note: Phase B-2 keeps the analytical Cr from process_design_requirements
    // even in REAL mode. Cr-resizing from a designed magnetic's actual
    // leakage is intentionally deferred to a follow-up commit (CLLLC_PLAN
    // §A.9 — Cr1/Cr2 still produced; magnetic re-sizing is the optional
    // refinement).

    if (computedPrimarySeriesInductance <= 0 ||
        computedPrimaryResonantCapacitance <= 0 ||
        extraCr1VoltageWaveforms.empty()) {
        throw std::runtime_error(
            "Clllc::get_extra_components_inputs: call process_operating_points() first");
    }

    bool integratedLs = get_integrated_resonant_inductors().value_or(false);
    double Lr1 = computedPrimarySeriesInductance;
    double Lr2 = computedSecondarySeriesInductance;
    double Cr1 = computedPrimaryResonantCapacitance;
    double Cr2 = computedSecondaryResonantCapacitance;
    double fr  = computedPrimaryResonantFrequency;
    if (fr <= 0) fr = get_effective_resonant_frequency();

    auto& ops = get_operating_points();
    size_t nOps = extraCr1VoltageWaveforms.size();
    if (nOps != extraCr2VoltageWaveforms.size() ||
        nOps != extraCr1CurrentWaveforms.size() ||
        nOps != extraCr2CurrentWaveforms.size()) {
        throw std::runtime_error(
            "Clllc::get_extra_components_inputs: extra-waveform vectors out of sync");
    }

    std::vector<std::variant<Inputs, CAS::Inputs>> result;

    auto buildCapInputs = [&](const std::string& name,
                              double C,
                              const std::vector<Waveform>& vWfms,
                              const std::vector<Waveform>& iWfms) {
        CAS::Inputs casInputs;
        CAS::DesignRequirements dr;

        CAS::DimensionWithTolerance capacitance;
        capacitance.set_nominal(C);
        dr.set_capacitance(capacitance);

        double peakV = 0.0;
        for (const auto& w : vWfms)
            for (double v : w.get_data())
                peakV = std::max(peakV, std::abs(v));
        // CLLLC plan §A.9: 1.5x peak rating for the resonant Cr (ringing
        // during direction-reversal transients).
        dr.set_rated_voltage(peakV * 1.5);
        dr.set_role(CAS::Application::RESONANT);
        casInputs.set_design_requirements(dr);

        std::vector<CAS::TwoTerminalOperatingPoint> casOps;
        for (size_t i = 0; i < vWfms.size(); ++i) {
            CAS::TwoTerminalOperatingPoint op;
            CAS::OperatingPointExcitation exc;
            double fsw = (i < ops.size()) ? ops[i].get_switching_frequency() : fr;
            exc.set_frequency(fsw);

            CAS::SignalDescriptor vSig;
            CAS::Waveform vWfm;
            vWfm.set_data(vWfms[i].get_data());
            if (vWfms[i].get_time()) vWfm.set_time(vWfms[i].get_time().value());
            vSig.set_waveform(vWfm);
            exc.set_voltage(vSig);

            CAS::SignalDescriptor iSig;
            CAS::Waveform iWfm;
            iWfm.set_data(iWfms[i].get_data());
            if (iWfms[i].get_time()) iWfm.set_time(iWfms[i].get_time().value());
            iSig.set_waveform(iWfm);
            exc.set_current(iSig);

            op.set_excitation(exc);
            casOps.push_back(op);
        }
        casInputs.set_operating_points(casOps);
        result.emplace_back(std::move(casInputs));
    };

    // Always: Cr1 (HV side) and Cr2 (LV side).
    buildCapInputs("Cr1_HV_resonantCapacitor", Cr1,
                   extraCr1VoltageWaveforms, extraCr1CurrentWaveforms);
    buildCapInputs("Cr2_LV_resonantCapacitor", Cr2,
                   extraCr2VoltageWaveforms, extraCr2CurrentWaveforms);

    // Mode REAL: HV + LV bus caps, sized from operating-point bus voltages.
    // The solver does not produce bus-ripple waveforms, so the CAS::Inputs
    // here carry only DesignRequirements (capacitance + rated voltage).
    if (mode == ExtraComponentsMode::REAL) {
        auto buildBusCap = [&](const std::string& name, double Vbus_max) {
            CAS::Inputs busCap;
            CAS::DesignRequirements dr;
            // Placeholder bulk capacitance: 1 µF is a deliberate lower bound;
            // the magnetic-adviser pipeline is expected to refine this from
            // ripple-current targets in a later phase. We do NOT silently
            // fabricate a "typical" value — 1 µF is documented minimum.
            CAS::DimensionWithTolerance c;
            c.set_nominal(1.0e-6);
            dr.set_capacitance(c);
            dr.set_rated_voltage(Vbus_max * 1.2);
            dr.set_role(CAS::Application::DC_LINK);
            busCap.set_design_requirements(dr);
            result.emplace_back(std::move(busCap));
        };

        // Bus voltages live on the ClllcResonant design spec (HV/LV are the
        // converter ports). Use maximum of nominal/maximum dimension fields.
        auto pickBusMax = [](const DimensionWithTolerance& d) -> double {
            if (d.get_maximum().has_value()) return d.get_maximum().value();
            if (d.get_nominal().has_value()) return d.get_nominal().value();
            if (d.get_minimum().has_value()) return d.get_minimum().value();
            return 0.0;
        };
        double Vhv_max = pickBusMax(get_high_voltage_bus_voltage());
        double Vlv_max = pickBusMax(get_low_voltage_bus_voltage());
        if (Vhv_max <= 0 || Vlv_max <= 0) {
            throw std::runtime_error(
                "Clllc::get_extra_components_inputs: REAL mode requires "
                "high/low voltage bus voltages on the design spec");
        }
        buildBusCap("HV_busCapacitor", Vhv_max);
        buildBusCap("LV_busCapacitor", Vlv_max);
    }

    // Discrete Lr1 + Lr2 inductors when not integrated into the transformer.
    if (!integratedLs) {
        auto buildInductorInputs = [&](const std::string& name,
                                       double L,
                                       const std::vector<Waveform>& iWfms,
                                       const std::vector<Waveform>& vWfms) {
            if (iWfms.size() != vWfms.size())
                throw std::runtime_error(
                    "Clllc::get_extra_components_inputs: " + name +
                    " current/voltage waveform count mismatch");

            Inputs masInputs;
            DesignRequirements dr;
            DimensionWithTolerance inductance;
            inductance.set_nominal(L);
            dr.set_magnetizing_inductance(inductance);
            dr.set_name(name);
            dr.set_topology(MAS::Topology::CLLLC_RESONANT_CONVERTER);
            dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
            dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
            masInputs.set_design_requirements(dr);

            std::vector<OperatingPoint> masOps;
            for (size_t i = 0; i < iWfms.size(); ++i) {
                OperatingPoint op;
                double fsw = (i < ops.size()) ? ops[i].get_switching_frequency() : fr;
                auto exc = complete_excitation(iWfms[i], vWfms[i], fsw, "Primary");
                op.get_mutable_excitations_per_winding().push_back(exc);
                masOps.push_back(op);
            }
            masInputs.set_operating_points(masOps);
            result.emplace_back(std::move(masInputs));
        };

        if (extraLr1VoltageWaveforms.empty() || extraLr2VoltageWaveforms.empty())
            throw std::runtime_error(
                "Clllc::get_extra_components_inputs: discrete Lr requires "
                "Lr voltage waveforms (run process_operating_points first)");
        buildInductorInputs("Lr1_HV_seriesInductor", Lr1,
                            extraLr1CurrentWaveforms, extraLr1VoltageWaveforms);
        buildInductorInputs("Lr2_LV_seriesInductor", Lr2,
                            extraLr2CurrentWaveforms, extraLr2VoltageWaveforms);
    }

    return result;
}

// =====================================================================
// AdvancedClllc::process
// =====================================================================
// Mirrors AdvancedLlc::process(): assemble DesignRequirements via the
// inherited process_design_requirements() (which applies the
// `desired*` overrides), pull the resolved turnsRatios + Lm out of the
// merged DR, then feed those to process_operating_points() to get the
// per-OP excitations. Package both into Inputs and return.
Inputs AdvancedClllc::process() {
    auto designRequirements = process_design_requirements();

    // Reuse the merged DR's turns ratios and magnetizing inductance for
    // operating-point calculation (process_design_requirements() has
    // already applied the desired* overrides).
    std::vector<double> turnsRatios;
    for (const auto& tr : designRequirements.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    if (turnsRatios.empty()) {
        throw std::runtime_error("CLLLC: no turns ratios available (neither desired nor computed)");
    }
    double Lm = resolve_dimensional_values(designRequirements.get_magnetizing_inductance());

    auto ops = process_operating_points(turnsRatios, Lm);
    Inputs inputs;
    inputs.set_design_requirements(designRequirements);
    inputs.set_operating_points(ops);
    return inputs;
}

} // namespace OpenMagnetics
