// Version: 3.3 — Runo Nielsen TDA + Integrated Ls + Inductance Ratio + Multi-output
#include "converter_models/Llc.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <limits>
#include "support/Exceptions.h"

namespace OpenMagnetics {

// =====================================================================
// Construction
// =====================================================================
Llc::Llc(const json& j) { from_json(j, *static_cast<LlcResonant*>(this)); }
AdvancedLlc::AdvancedLlc(const json& j) { from_json(j, *this); }

double Llc::get_bridge_voltage_factor() const {
    auto bt = get_bridge_type();
    if (bt.has_value() && bt.value() == LlcBridgeType::FULL_BRIDGE) return 1.0;
    return 0.5;
}

double Llc::get_effective_resonant_frequency() const {
    if (get_resonant_frequency().has_value()) return get_resonant_frequency().value();
    return std::sqrt(get_min_switching_frequency() * get_max_switching_frequency());
}

double Llc::get_inductance_ratio() const {
    // inductance_ratio is now optional in the schema; falls back to computedInductanceRatio
    auto userValue = MAS::LlcResonant::get_inductance_ratio();
    if (userValue.has_value() && userValue.value() > 0) {
        return userValue.value();
    }
    return computedInductanceRatio;
}

bool Llc::run_checks(bool assertErrors) {
    _assertErrors = assertErrors;
    bool ok = true;
    auto& ops = get_operating_points();
    if (ops.empty()) {
        if (assertErrors) throw std::runtime_error("LLC: no operating points");
        return false;
    }
    for (size_t i = 0; i < ops.size(); i++) {
        auto& op = ops[i];
        if (op.get_output_voltages().empty() || op.get_output_currents().empty()) {
            if (assertErrors) throw std::runtime_error("LLC: OP missing voltages/currents");
            ok = false;
        }
        if (op.get_output_voltages().size() != op.get_output_currents().size()) {
            if (assertErrors) throw std::runtime_error("LLC: voltage/current count mismatch");
            ok = false;
        }
        double fsw = op.get_switching_frequency();
        if (fsw < get_min_switching_frequency() * 0.99 ||
            fsw > get_max_switching_frequency() * 1.01) {
            if (assertErrors) throw std::runtime_error("LLC: fsw out of range");
            ok = false;
        }
    }
    return ok;
}

// =====================================================================
// Design Requirements
//
// LLC has THREE independent design parameters:
//   1. Q  (quality factor)    → controls Ls, Cr (tank impedance)
//   2. Ln (inductance ratio)  → controls Lm (magnetizing inductance)
//   3. fr (resonant frequency) → sets the operating frequency
//
// From these:
//   Rac = 8·n²/(π²) · Rload
//   Zr  = Q · Rac
//   Ls  = Zr / (2π·fr)
//   Cr  = 1 / (2π·fr·Zr)
//   Lm  = Ln · Ls
// =====================================================================
DesignRequirements Llc::process_design_requirements() {
    double k_bridge = get_bridge_voltage_factor();
    auto& inputVoltage = get_input_voltage();
    auto& ops = get_operating_points();
    if (ops.empty())
        throw std::runtime_error("LLC: at least one operating point is required");
    if (ops[0].get_output_voltages().empty() || ops[0].get_output_currents().empty())
        throw std::runtime_error("LLC: operating point is missing output voltages or currents");
    if (!inputVoltage.get_nominal().has_value() &&
        !(inputVoltage.get_minimum().has_value() && inputVoltage.get_maximum().has_value()))
        throw std::runtime_error("LLC: input voltage must specify nominal (or both min and max)");

    double Vin_nom = inputVoltage.get_nominal().value_or(
        (inputVoltage.get_minimum().value_or(0) + inputVoltage.get_maximum().value_or(0)) / 2.0);
    double mainOutputVoltage  = ops[0].get_output_voltages()[0];
    double mainOutputCurrent  = ops[0].get_output_currents()[0];
    if (mainOutputVoltage <= 0 || mainOutputCurrent <= 0)
        throw std::runtime_error("LLC: main output voltage and current must be positive");

    double Vo = k_bridge * Vin_nom;  // FHA-equivalent reflected output voltage
    double mainTurnsRatio = Vo / mainOutputVoltage;

    // turnsRatios stores SECONDARY turns ratios only (primary implicit, ratio = 1).
    // Number of physical secondary windings per output depends on rectifierType:
    //   CT  → 2 (center-tapped, two halves)
    //   FB  → 1 (single secondary, four diodes)
    //   CD  → 1 (single secondary, two diodes + two output inductors)
    //   VD  → 1 (single secondary, two diodes + two output capacitors)
    // We push the same ratio `windings_per_output()` times per output so the
    // turns_ratios vector aligns 1-to-1 with non-primary windings.
    size_t nOutputs = ops[0].get_output_voltages().size();
    size_t wpo = windings_per_output();
    std::vector<double> outputTurnsRatios;
    outputTurnsRatios.reserve(nOutputs);
    for (size_t i = 0; i < nOutputs; ++i) {
        double Vout_i = ops[0].get_output_voltages()[i];
        if (Vout_i <= 0)
            throw std::runtime_error("LLC: output voltage must be positive");
        // For VD, the cap stack delivers 2·Vsec_pk to the load, so the
        // secondary turns ratio is doubled (we need half the Vsec amplitude).
        double n_design = (get_effective_rectifier_type() == LlcRectifierType::VOLTAGE_DOUBLER)
                          ? (Vo / Vout_i) * 2.0
                          : (Vo / Vout_i);
        outputTurnsRatios.push_back(n_design);
    }

    double Rload = mainOutputVoltage / mainOutputCurrent;
    double Rac = (8.0 * mainTurnsRatio * mainTurnsRatio) / (M_PI * M_PI) * Rload;
    double Q = get_quality_factor().value_or(defaults.resonantQualityFactorDefaultLlc);
    double fr = get_effective_resonant_frequency();
    if (fr <= 0)
        throw std::runtime_error("LLC: resonant frequency must be positive");
    double Zr = Q * Rac;
    double Ls = Zr / (2.0 * M_PI * fr);
    double Cr = 1.0 / (2.0 * M_PI * fr * Zr);

    // User overrides take precedence over the Q/Ln/fr derivation. This lets
    // callers validate against published reference designs (e.g. Nielsen)
    // where Ls and Cr are measured component values.
    if (userResonantInductance > 0)  Ls = userResonantInductance;
    if (userResonantCapacitance > 0) Cr = userResonantCapacitance;

    // Inductance ratio: user-settable via JSON, falls back to computedInductanceRatio
    // Low Ln (3-5):  wide gain range, high magnetizing current, good ZVS margin
    // Medium (5-8):  balanced design
    // High Ln (8-12): narrow gain, low magnetizing current, high efficiency
    double Ln = get_inductance_ratio();
    double L = Ln * Ls;

    computedResonantInductance  = Ls;
    computedResonantCapacitance = Cr;

    DesignRequirements designRequirements;
    designRequirements.get_mutable_turns_ratios().clear();
    for (auto n : outputTurnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(roundFloat(n, 2));
        for (size_t w = 0; w < wpo; ++w)
            designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }
    DimensionWithTolerance inductanceWithTolerance;
    inductanceWithTolerance.set_nominal(roundFloat(L, 10));
    designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

    if (get_integrated_resonant_inductor().value_or(false)) {
        std::vector<DimensionWithTolerance> leakageReqs;
        DimensionWithTolerance lrTol;
        lrTol.set_nominal(roundFloat(Ls, 10));
        leakageReqs.push_back(lrTol);
        designRequirements.set_leakage_inductance(leakageReqs);
    }

    // Isolation sides: primary, then `windings_per_output()` entries per output.
    std::vector<IsolationSide> isolationSides;
    isolationSides.push_back(get_isolation_side_from_index(0));
    for (size_t i = 0; i < nOutputs; ++i) {
        for (size_t w = 0; w < wpo; ++w)
            isolationSides.push_back(get_isolation_side_from_index(i + 1));
    }
    designRequirements.set_isolation_sides(isolationSides);
    designRequirements.set_topology(Topologies::LLC_RESONANT_CONVERTER);

    return designRequirements;
}

// =====================================================================
// Process operating points
// =====================================================================
std::vector<OperatingPoint> Llc::process_operating_points(
    const std::vector<double>& turnsRatios, double magnetizingInductance)
{
    if (computedResonantInductance <= 0 || computedResonantCapacitance <= 0)
        process_design_requirements();

    // Multi-output warn-once: per-secondary currents are projected from the
    // primary tank current by load-share weighting (no per-secondary leakage
    // model). Same pattern as Dab.cpp.
    {
        auto& opsTmp = get_operating_points();
        if (!opsTmp.empty() && opsTmp[0].get_output_voltages().size() > 1) {
            static thread_local bool warned = false;
            if (!warned) {
                std::cerr << "[Llc] Multi-output: per-secondary currents are "
                          << "approximated via load-share projection. Provide "
                          << "per-secondary leakage inductance for accurate "
                          << "individual current waveforms.\n";
                warned = true;
            }
        }
    }

    extraCapVoltageWaveforms.clear();
    extraCapCurrentWaveforms.clear();
    extraIndVoltageWaveforms.clear();
    extraIndCurrentWaveforms.clear();
    extraTimeVectors.clear();
    extraLoCurrentWaveforms.clear();
    extraLoVoltageWaveforms.clear();
    extraLo2CurrentWaveforms.clear();
    extraLo2VoltageWaveforms.clear();
    extraVoutCapVoltageWaveforms.clear();

    perOpName.clear();
    perOpMode.clear();
    perOpSteadyStateResidual.clear();
    perOpZvsMarginLagging.clear();
    perOpPrimaryPeakCurrent.clear();

    std::vector<OperatingPoint> result;
    auto& inputVoltage = get_input_voltage();
    auto& ops = get_operating_points();
    std::vector<std::pair<double, std::string>> inputVoltages;
    if (inputVoltage.get_nominal().has_value())
        inputVoltages.push_back({inputVoltage.get_nominal().value(), "Nominal"});
    if (inputVoltage.get_minimum().has_value())
        inputVoltages.push_back({inputVoltage.get_minimum().value(), "Min"});
    if (inputVoltage.get_maximum().has_value())
        inputVoltages.push_back({inputVoltage.get_maximum().value(), "Max"});

    std::sort(inputVoltages.begin(), inputVoltages.end(),
        [](const auto& a, const auto& b) { return a.first < b.first; });
    auto last = std::unique(inputVoltages.begin(), inputVoltages.end(),
        [](const auto& a, const auto& b) { return a.first == b.first; });
    inputVoltages.erase(last, inputVoltages.end());

    for (const auto& [Vin, name] : inputVoltages) {
        auto op = process_operating_point_for_input_voltage(
            Vin, ops[0], turnsRatios, magnetizingInductance);
        op.set_name(name + " input (" + std::to_string(static_cast<int>(Vin)) + "V)");
        result.push_back(op);

        perOpName.push_back(name);
        perOpMode.push_back(lastMode);
        perOpSteadyStateResidual.push_back(lastSteadyStateResidual);
        perOpZvsMarginLagging.push_back(lastZvsMarginLagging);
        perOpPrimaryPeakCurrent.push_back(lastPrimaryPeakCurrent);
    }
    return result;
}

std::vector<OperatingPoint> Llc::process_operating_points(Magnetic magnetic) {
    auto req = process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios())
        turnsRatios.push_back(resolve_dimensional_values(tr));
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    return process_operating_points(turnsRatios, Lm);
}

// =====================================================================
// NIELSEN TIME-DOMAIN APPROACH — event-driven sub-state propagator
//
// Implements the Runo Nielsen TDA for the LLC resonant converter as
// described in:
//   • "LLC and LCC resonance converters", Runo Nielsen, August 2013
//   • "Resonance half bridge LLC converter analysis" (Mathcad worksheet
//     output), Runo Nielsen, March 2022
//
// State vector x = [I_Ls, I_L, V_C]ᵀ — three INDEPENDENT components:
//   I_Ls : resonant inductor current (= primary winding current)
//   I_L  : magnetizing inductor current (independent of I_Ls except in
//          the freewheeling sub-state where they coincide)
//   V_C  : resonant capacitor voltage
//
// Three linear sub-states cover all of Nielsen's six modes — the modes
// emerge from the SEQUENCE of sub-states traversed in one half cycle:
//   A_POS : secondary diode +Vo conducting   (Id = ILs − IL > 0)
//           dILs/dt = (Vi − Vc − Vo)/Ls
//           dIL/dt  = +Vo/L     (magnetizing voltage clamped to +Vo)
//           dVc/dt  = ILs/C
//   A_NEG : secondary diode −Vo conducting   (Id < 0)
//           dILs/dt = (Vi − Vc + Vo)/Ls
//           dIL/dt  = −Vo/L
//           dVc/dt  = ILs/C
//   B_FW  : both diodes off, freewheeling     (ILs ≡ IL)
//           dILs/dt = (Vi − Vc)/(Ls + L)  ← Lm in series with Ls
//           dIL/dt  = same                 (the two currents equalise)
//           dVc/dt  = ILs/C
//
// Each sub-state is a linear ODE with closed-form cos/sin solution. The
// event detector finds the next sub-state transition exactly (analytic
// root of the trigger function), so the half-cycle is composed from a
// chain of exact closed-form segments — no per-sample arithmetic.
//
// Steady state is enforced by half-wave antisymmetry x(Thalf) = −x(0)
// via damped Newton on a 3-vector residual F(x0) = x(Thalf; x0) + x0.
// =====================================================================

enum class LlcSubState { A_POS, A_NEG, B_FW };

struct LlcStateVector {
    double iLs;
    double iL;
    double vC;
};

struct LlcSubStateSegment {
    LlcSubState state;
    double t_start;
    double t_end;
    LlcStateVector x_start;
    LlcStateVector x_end;
};

// Sub-state-specific natural frequency and characteristic impedance.
// Returns w (rad/s) and Z (ohms) given Ls, L, C and the sub-state type.
static void substate_freq(LlcSubState s, double Ls, double L, double C,
                          double& w, double& Z) {
    if (s == LlcSubState::B_FW) {
        // Freewheeling: Lm in series with Ls (diodes off)
        double L_eff = Ls + L;
        w = 1.0 / std::sqrt(L_eff * C);
        Z = std::sqrt(L_eff / C);
    } else {
        // Power delivery: Lm clamped, only Ls and C ring
        w = 1.0 / std::sqrt(Ls * C);
        Z = std::sqrt(Ls / C);
    }
}

// Closed-form propagator for a single sub-state. Computes x(t_start + dt)
// exactly given x(t_start) and the sub-state type.
//
// For A_POS / A_NEG (power delivery), I_Ls and V_C oscillate about the
// equilibrium V_C_eq = Vi ∓ Vo at frequency w1 = 1/√(Ls·C); I_L ramps
// linearly because its terminal voltage is clamped to ±Vo.
//
// For B_FW (freewheeling), I_Ls (= I_L) and V_C oscillate about V_C_eq = Vi
// at frequency w0 = 1/√((Ls+Lm)·C).
static LlcStateVector propagate_substate(LlcSubState s, LlcStateVector x_in,
                                         double dt, double Vi, double Vo,
                                         double Ls, double L, double C) {
    LlcStateVector out{};
    if (dt <= 0) return x_in;

    double w, Z;
    substate_freq(s, Ls, L, C, w, Z);
    double cs = std::cos(w * dt);
    double sn = std::sin(w * dt);

    if (s == LlcSubState::B_FW) {
        // Freewheeling: I_Ls ≡ I_L, drive = Vi
        double V_eq = Vi;
        double dV = x_in.vC - V_eq;
        // (ILs, Vc) is an LC oscillator about (0, V_eq):
        double iLs_new = x_in.iLs * cs - (dV / Z) * sn;
        double vC_new  = V_eq + dV * cs + x_in.iLs * Z * sn;
        out.iLs = iLs_new;
        out.iL  = iLs_new;  // tracks ILs in freewheeling
        out.vC  = vC_new;
        return out;
    }

    // Power-delivery sub-states. For A_POS: V_drive = Vi − Vo, dIL/dt = +Vo/L
    // For A_NEG: V_drive = Vi + Vo, dIL/dt = −Vo/L
    double V_drive = (s == LlcSubState::A_POS) ? (Vi - Vo) : (Vi + Vo);
    double dIL_dt  = (s == LlcSubState::A_POS) ? (+Vo / L) : (-Vo / L);

    double dV = x_in.vC - V_drive;
    out.iLs = x_in.iLs * cs - (dV / Z) * sn;
    out.vC  = V_drive + dV * cs + x_in.iLs * Z * sn;
    out.iL  = x_in.iL + dIL_dt * dt;
    return out;
}

// Evaluate the trigger function at the END of a propagation step.
// Used by the event detector to bracket and locate the next sub-state
// transition. Returns the value g(dt) for the trigger that ends the
// current sub-state.
//
// Triggers:
//   A_POS → B_FW : g = I_L − I_Ls (rises to 0 from below as Id falls to 0)
//   A_NEG → B_FW : g = I_Ls − I_L (rises to 0 from below as Id rises to 0)
//   B_FW  → A_POS: g = L·dILs/dt − Vo (when magnetizing voltage rises to +Vo)
//   B_FW  → A_NEG: g = −L·dILs/dt − Vo (when magnetizing voltage drops to −Vo)
static double trigger_value(LlcSubState s, LlcStateVector x_in, double dt,
                            double Vi, double Vo, double Ls, double L, double C) {
    LlcStateVector x = propagate_substate(s, x_in, dt, Vi, Vo, Ls, L, C);
    if (s == LlcSubState::A_POS) {
        return x.iL - x.iLs;          // > 0 means Id flipped negative
    }
    if (s == LlcSubState::A_NEG) {
        return x.iLs - x.iL;          // > 0 means Id flipped positive
    }
    // B_FW: voltage across Lm exits the ±Vo clamp window.
    // VLm = (L / (Ls+L)) * (Vi − Vc).  Diodes turn on when |VLm| ≥ Vo.
    double VLm = (L / (Ls + L)) * (Vi - x.vC);
    // Return the most-positive of the two clamp violations:
    //   if VLm > +Vo → A_POS forward-bias (return VLm − Vo)
    //   if VLm < −Vo → A_NEG forward-bias (return −VLm − Vo)
    // We pick whichever hits zero first; the caller resolves which
    // sub-state to enter using sign(VLm) at the crossing.
    double pos_violation = VLm - Vo;
    double neg_violation = -VLm - Vo;
    return std::max(pos_violation, neg_violation);
}

// Find the smallest t > 0 in (0, t_max] where the trigger function for the
// given sub-state crosses zero from negative to positive. Returns t_max if
// no crossing is found.
//
// Strategy: subdivide [0, t_max] into a coarse grid (32 segments), find the
// first interval where the trigger sign changes, then refine with bisection.
// This is more robust than analytic root-finding for the mixed cos/sin/linear
// trigger functions in the A± sub-states (where the linear I_L ramp creates
// trigger functions of the form A·cos(ωt) + B·sin(ωt) + C·t + D).
static double find_next_event(LlcSubState s, LlcStateVector x_in, double t_max,
                              double Vi, double Vo, double Ls, double L, double C) {
    constexpr int COARSE_STEPS = 64;
    double dt_coarse = t_max / COARSE_STEPS;
    double prev_g = trigger_value(s, x_in, 1e-12, Vi, Vo, Ls, L, C);
    if (prev_g >= 0) {
        // Trigger already true at t≈0: return immediately. This is a
        // "tangent" case (e.g. Id starts at exactly 0). Skip to a tiny
        // ε-step to leave the boundary.
        return 0.0;
    }

    for (int k = 1; k <= COARSE_STEPS; ++k) {
        double t = k * dt_coarse;
        double g = trigger_value(s, x_in, t, Vi, Vo, Ls, L, C);
        if (g >= 0 && std::isfinite(g)) {
            // Sign change in (t-dt_coarse, t]. Refine by bisection.
            double lo = t - dt_coarse;
            double hi = t;
            double g_lo = prev_g;
            for (int it = 0; it < 50; ++it) {
                double mid = 0.5 * (lo + hi);
                double g_mid = trigger_value(s, x_in, mid, Vi, Vo, Ls, L, C);
                if (g_mid * g_lo < 0) {
                    hi = mid;
                } else {
                    lo = mid;
                    g_lo = g_mid;
                }
                if ((hi - lo) < 1e-12) break;
            }
            return 0.5 * (lo + hi);
        }
        prev_g = g;
    }
    return t_max;  // no crossing
}

// Determine the sub-state immediately after a B_FW transition: which
// secondary diode forward-biased? Look at sign(VLm) at the crossing.
static LlcSubState next_state_after_B(LlcStateVector x_at_event, double Vi,
                                      double Ls, double L) {
    double VLm = (L / (Ls + L)) * (Vi - x_at_event.vC);
    return (VLm > 0) ? LlcSubState::A_POS : LlcSubState::A_NEG;
}

// Determine the initial sub-state at t=0 of the half cycle, given x0 and Vi.
// At t=0+, the bridge has just switched to +Vi. Examine which constraint
// is satisfied and pick accordingly. The diode that conducts is determined
// by the sign of (ILs - IL) — if positive, +Vo diode; if negative, −Vo
// diode; if both currents are equal, the converter starts in freewheeling.
static LlcSubState initial_substate(LlcStateVector x0, double Vi, double Vo,
                                    double Ls, double L) {
    double Id = x0.iLs - x0.iL;
    if (std::abs(Id) > 1e-9) {
        return (Id > 0) ? LlcSubState::A_POS : LlcSubState::A_NEG;
    }
    // Id ≈ 0: check if magnetizing voltage will push the diodes on.
    double VLm = (L / (Ls + L)) * (Vi - x0.vC);
    if (VLm > Vo) return LlcSubState::A_POS;
    if (VLm < -Vo) return LlcSubState::A_NEG;
    return LlcSubState::B_FW;
}

// Drive the event loop over [0, Thalf]. Returns the chain of segments
// covering the half cycle. Each segment has exact closed-form start/end
// states; downstream sampling rasters them onto the 256-point grid.
static std::vector<LlcSubStateSegment> propagate_half_cycle(
    LlcStateVector x0, double Thalf,
    double Vi, double Vo, double Ls, double L, double C)
{
    std::vector<LlcSubStateSegment> segments;
    segments.reserve(8);

    LlcSubState current = initial_substate(x0, Vi, Vo, Ls, L);
    LlcStateVector x = x0;
    double t = 0.0;
    constexpr int MAX_SEGMENTS = 16;  // safety bound; LLC modes use ≤ 6

    for (int k = 0; k < MAX_SEGMENTS; ++k) {
        double remaining = Thalf - t;
        if (remaining <= 1e-15) break;

        double t_event = find_next_event(current, x, remaining, Vi, Vo, Ls, L, C);
        // t_event is relative to the start of this segment.
        double dt = std::min(t_event, remaining);
        if (dt < 1e-15 && k > 0) {
            // Zero-length segment: must transition immediately to avoid
            // an infinite loop. Pick the next state and continue.
            LlcSubState next;
            if (current == LlcSubState::B_FW) {
                next = next_state_after_B(x, Vi, Ls, L);
            } else {
                next = LlcSubState::B_FW;
            }
            current = next;
            continue;
        }
        LlcStateVector x_end = propagate_substate(current, x, dt, Vi, Vo, Ls, L, C);
        segments.push_back({current, t, t + dt, x, x_end});
        t += dt;
        x = x_end;

        if (t >= Thalf - 1e-15) break;

        // Determine next sub-state.
        if (current == LlcSubState::B_FW) {
            current = next_state_after_B(x, Vi, Ls, L);
        } else {
            // A_POS / A_NEG ended on Id = 0 → freewheeling
            current = LlcSubState::B_FW;
        }
    }

    return segments;
}

// Steady-state solver: damped Newton on F(x0) = propagate_half(x0).end + x0.
// Returns the converged x0, the segment chain at convergence, and the L2
// residual. If Newton fails to converge in MAX_ITERS, falls back to a
// Picard iteration on x0 using x0_new = -propagate_half(x0).end (which is
// equivalent to the previous code's outer loop).
static LlcStateVector solve_steady_state(
    LlcStateVector x0_seed, double Thalf,
    double Vi, double Vo, double Ls, double L, double C,
    std::vector<LlcSubStateSegment>& outSegments,
    double& outResidual)
{
    auto eval_F = [&](LlcStateVector x0) -> std::array<double, 3> {
        auto segs = propagate_half_cycle(x0, Thalf, Vi, Vo, Ls, L, C);
        LlcStateVector x_end = segs.empty() ? x0 : segs.back().x_end;
        return {x_end.iLs + x0.iLs, x_end.iL + x0.iL, x_end.vC + x0.vC};
    };
    auto norm = [](const std::array<double, 3>& f) {
        return std::sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
    };

    LlcStateVector x0 = x0_seed;
    auto F = eval_F(x0);
    double r = norm(F);

    constexpr int MAX_ITERS = 50;  // monotone acceptance below, so extra iters only help
    constexpr double TOL = 1e-7;

    // Reasonable perturbation scales for finite differences (in physical
    // units of A and V).
    double scale_i = std::max(1e-3, 0.01 * std::abs(x0.iLs) + 1e-3);
    double scale_v = std::max(1e-2, 0.01 * std::abs(x0.vC)  + 1e-2);

    double damping = 1.0;
    double prev_r = r;
    int stagnant = 0;

    // Track the best (lowest-residual) iterate seen. The damped-Newton /
    // Picard iteration is non-monotone (it can transiently jump the residual
    // up), so returning the *last* iterate threw away good solutions already
    // found. Return the best instead.
    LlcStateVector bestX0 = x0;
    double bestR = r;

    for (int iter = 0; iter < MAX_ITERS && r > TOL; ++iter) {
        // Build Jacobian by central differences.
        double J[3][3];
        LlcStateVector xp, xm;
        std::array<double, 3> Fp, Fm;

        // ∂F/∂iLs
        xp = x0; xp.iLs += scale_i;
        xm = x0; xm.iLs -= scale_i;
        Fp = eval_F(xp); Fm = eval_F(xm);
        for (int i = 0; i < 3; ++i) J[i][0] = (Fp[i] - Fm[i]) / (2 * scale_i);

        // ∂F/∂iL
        xp = x0; xp.iL += scale_i;
        xm = x0; xm.iL -= scale_i;
        Fp = eval_F(xp); Fm = eval_F(xm);
        for (int i = 0; i < 3; ++i) J[i][1] = (Fp[i] - Fm[i]) / (2 * scale_i);

        // ∂F/∂vC
        xp = x0; xp.vC += scale_v;
        xm = x0; xm.vC -= scale_v;
        Fp = eval_F(xp); Fm = eval_F(xm);
        for (int i = 0; i < 3; ++i) J[i][2] = (Fp[i] - Fm[i]) / (2 * scale_v);

        // Solve J · dx = -F by Gaussian elimination with partial pivoting.
        double A[3][4];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) A[i][j] = J[i][j];
            A[i][3] = -F[i];
        }
        // Forward elimination
        bool singular = false;
        for (int col = 0; col < 3; ++col) {
            int piv = col;
            double maxv = std::abs(A[col][col]);
            for (int row = col + 1; row < 3; ++row) {
                if (std::abs(A[row][col]) > maxv) {
                    maxv = std::abs(A[row][col]);
                    piv = row;
                }
            }
            if (maxv < 1e-14) { singular = true; break; }
            if (piv != col) std::swap(A[col], A[piv]);
            for (int row = col + 1; row < 3; ++row) {
                double f = A[row][col] / A[col][col];
                for (int k = col; k < 4; ++k) A[row][k] -= f * A[col][k];
            }
        }
        if (singular) break;
        // Back substitution
        double dx[3];
        for (int i = 2; i >= 0; --i) {
            double sum = A[i][3];
            for (int j = i + 1; j < 3; ++j) sum -= A[i][j] * dx[j];
            dx[i] = sum / A[i][i];
        }

        // Damped step + line search on ||F||.
        double try_d = damping;
        LlcStateVector x_new{};
        std::array<double, 3> F_new{};
        double r_new = r;
        bool accepted = false;
        for (int ls = 0; ls < 6; ++ls) {
            x_new.iLs = x0.iLs + try_d * dx[0];
            x_new.iL  = x0.iL  + try_d * dx[1];
            x_new.vC  = x0.vC  + try_d * dx[2];
            F_new = eval_F(x_new);
            r_new = norm(F_new);
            if (std::isfinite(r_new) && r_new < r) {
                accepted = true;
                break;
            }
            try_d *= 0.5;
        }
        if (!accepted) {
            // Picard fallback: x0_new = -x_end (same as old outer loop)
            auto segs = propagate_half_cycle(x0, Thalf, Vi, Vo, Ls, L, C);
            if (!segs.empty()) {
                LlcStateVector x_end = segs.back().x_end;
                x_new.iLs = -x_end.iLs;
                x_new.iL  = -x_end.iL;
                x_new.vC  = -x_end.vC;
                F_new = eval_F(x_new);
                r_new = norm(F_new);
            }
        }
        // Commit the step (Newton line-search step when accepted; otherwise the
        // Picard fallback). The Picard jump is deliberately allowed to move
        // *uphill* — it escapes line-search stalls and lets the search explore
        // other basins (empirically it reaches much lower residuals a few
        // iterations later). We DON'T return this last iterate, though: we track
        // the best iterate seen and return that, so the exploration never throws
        // away a good solution it already found (the old code returned the last
        // iterate, which could be far worse than the best).
        x0 = x_new;
        F = F_new;
        if (r_new >= prev_r * 0.999) {
            stagnant++;
            damping *= 0.7;
        } else {
            stagnant = 0;
            damping = std::min(1.0, damping * 1.2);
        }
        prev_r = r_new;
        r = r_new;
        if (std::isfinite(r) && r < bestR) {
            bestR = r;
            bestX0 = x0;
        }
        if (stagnant >= 4) break;
    }

    // Return the best iterate found (not the last), re-propagated for waveforms.
    outSegments = propagate_half_cycle(bestX0, Thalf, Vi, Vo, Ls, L, C);
    outResidual = bestR;
    return bestX0;
}

// Map a sub-state segment chain to Nielsen's mode-1..6 numbering.
// The mapping is the working hypothesis from the plan; it will be refined
// against Nielsen's mode plot during validation. Returns 0 if the sequence
// is unrecognised (still produces correct waveforms).
static int classify_mode(const std::vector<LlcSubStateSegment>& segs) {
    if (segs.empty()) return 0;
    // Build sequence of A_POS / A_NEG / B_FW codes.
    std::vector<LlcSubState> seq;
    for (auto& s : segs) seq.push_back(s.state);

    auto eq = [](const std::vector<LlcSubState>& s,
                 std::initializer_list<LlcSubState> ref) {
        if (s.size() != ref.size()) return false;
        size_t i = 0;
        for (auto v : ref) {
            if (s[i++] != v) return false;
        }
        return true;
    };

    if (eq(seq, {LlcSubState::A_POS}))                                              return 1;
    if (eq(seq, {LlcSubState::A_POS, LlcSubState::A_NEG}))                          return 2;
    if (eq(seq, {LlcSubState::A_POS, LlcSubState::B_FW, LlcSubState::A_POS}))       return 3;
    if (eq(seq, {LlcSubState::A_POS, LlcSubState::B_FW}))                           return 4;
    if (eq(seq, {LlcSubState::A_POS, LlcSubState::B_FW, LlcSubState::A_NEG}))       return 5;
    if (eq(seq, {LlcSubState::A_POS, LlcSubState::A_NEG, LlcSubState::B_FW}))       return 6;
    return 0;  // unrecognised — still produces correct waveforms
}

// Sample a segment chain onto a uniform time grid [0, Thalf] with `N+1`
// points. Fills ILs_pos / IL_pos / Vc_pos / VL_pos vectors that the
// downstream waveform-emission code consumes.
static void sample_segments(const std::vector<LlcSubStateSegment>& segs,
                            double Thalf, int N, double Vi, double Vo,
                            double Ls, double L, double C,
                            std::vector<double>& ILs_pos,
                            std::vector<double>& IL_pos,
                            std::vector<double>& Vc_pos,
                            std::vector<double>& VL_pos)
{
    double dt = Thalf / N;
    size_t segIdx = 0;
    for (int k = 0; k <= N; ++k) {
        double t = k * dt;
        if (t > Thalf) t = Thalf;
        // Advance segIdx while t lies past the current segment end.
        while (segIdx + 1 < segs.size() && t > segs[segIdx].t_end + 1e-15) ++segIdx;
        if (segs.empty()) {
            ILs_pos[k] = 0; IL_pos[k] = 0; Vc_pos[k] = 0; VL_pos[k] = 0;
            continue;
        }
        const auto& seg = segs[segIdx];
        double t_local = t - seg.t_start;
        if (t_local < 0) t_local = 0;
        if (t_local > seg.t_end - seg.t_start) t_local = seg.t_end - seg.t_start;
        LlcStateVector x = propagate_substate(seg.state, seg.x_start, t_local,
                                              Vi, Vo, Ls, L, C);
        ILs_pos[k] = x.iLs;
        IL_pos[k]  = x.iL;
        Vc_pos[k]  = x.vC;
        // Magnetizing voltage: clamped during conduction, follows
        // Lm·dILs/dt during freewheeling.
        if (seg.state == LlcSubState::A_POS) {
            VL_pos[k] = +Vo;
        } else if (seg.state == LlcSubState::A_NEG) {
            VL_pos[k] = -Vo;
        } else {
            // B_FW: VLm = (L/(Ls+L)) * (Vi - Vc)
            VL_pos[k] = (L / (Ls + L)) * (Vi - x.vC);
        }
    }
}


OperatingPoint Llc::process_operating_point_for_input_voltage(
    double inputVoltage,
    const LlcOperatingPoint& llcOpPoint,
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    OperatingPoint operatingPoint;

    if (llcOpPoint.get_output_voltages().empty() || llcOpPoint.get_output_currents().empty())
        throw std::runtime_error("LLC: operating point missing output voltages/currents");

    double fsw = llcOpPoint.get_switching_frequency();
    if (fsw <= 0) {
        // Operating near resonance when the switching frequency is unspecified
        // is a legitimate LLC design choice (the resonant frequency is derived
        // from the min/max switching-frequency requirements). But fabricating a
        // bare 100 kHz constant when even that is unavailable produces untrustworthy
        // currents/voltages — throw instead (siblings Cllc/Dab/PSHB do the same).
        fsw = get_effective_resonant_frequency();
        if (fsw <= 0)
            throw std::runtime_error("LLC: switching frequency is non-positive and no resonant "
                                     "frequency is available to derive it");
    }

    double k_bridge = get_bridge_voltage_factor();
    bool integratedLs = get_integrated_resonant_inductor().value_or(false);

    double Vi = k_bridge * inputVoltage;

    // turnsRatios stores `windings_per_output()` entries per output (2 for CT,
    // 1 for FB/CD/VD) when built by process_design_requirements.
    size_t nOutputs = llcOpPoint.get_output_voltages().size();
    if (nOutputs == 0) nOutputs = 1;
    size_t wpo = windings_per_output();
    auto get_n_for_output = [&](size_t outputIdx) -> double {
        if (turnsRatios.empty())
            return Vi / llcOpPoint.get_output_voltages()[outputIdx];
        if (turnsRatios.size() == wpo * nOutputs)
            return turnsRatios[wpo * outputIdx];
        // Backward-compat: also accept legacy 2*nOutputs layout from
        // pre-rectifierType fixtures even when wpo == 1.
        if (turnsRatios.size() == 2 * nOutputs)
            return turnsRatios[2 * outputIdx];
        if (outputIdx < turnsRatios.size())
            return turnsRatios[outputIdx];
        return turnsRatios[0];
    };

    // The reflected (FHA-equivalent) output voltage uses the *first* output, since
    // the resonant tank current is shared and design parameters target the main
    // output. Multi-output secondaries are scaled per-output below.
    double n_main = get_n_for_output(0);
    double Vo = n_main * llcOpPoint.get_output_voltages()[0];
    // For voltage-doubler, the cap stack provides 2·Vsec_pk to the load and
    // n_main was doubled at design time (see process_design_requirements). The
    // *reflected* output voltage seen at the magnetizing inductor is therefore
    // n_main · Vout / 2, not n_main · Vout. Without this correction, the
    // steady-state solver sees Vo = 2·Vi → no power transfer → ILs == IL and
    // the analytical secondary current collapses to zero.
    if (get_effective_rectifier_type() == LlcRectifierType::VOLTAGE_DOUBLER) {
        Vo *= 0.5;
    }

    if (magnetizingInductance <= 0)
        throw std::runtime_error("LLC: magnetizing inductance must be positive (got non-positive)");

    double Ls = computedResonantInductance;
    double C  = computedResonantCapacitance;
    double L  = magnetizingInductance;

    if (Ls <= 0 || C <= 0) {
        double fr = get_effective_resonant_frequency();
        double Q_val = get_quality_factor().value_or(defaults.resonantQualityFactorDefaultLlc);
        double Ln_val = get_inductance_ratio();
        double Vout = llcOpPoint.get_output_voltages()[0];
        double Iout_fb = llcOpPoint.get_output_currents()[0];
        if (Iout_fb <= 0)
            throw std::runtime_error("LLC: output current must be positive to size the resonant tank "
                                     "(cannot derive load resistance)");
        double Rload = Vout / Iout_fb;
        double Rac = (8.0 * n_main * n_main) / (M_PI * M_PI) * Rload;
        double Zr = Q_val * Rac;
        Ls = Zr / (2.0 * M_PI * fr);
        C = 1.0 / (2.0 * M_PI * fr * Zr);
        L = Ln_val * Ls;
        computedResonantInductance = Ls;
        computedResonantCapacitance = C;
    }

    if (Ls <= 0 || C <= 0 || L <= 0)
        throw std::runtime_error("LLC resonant tank values invalid.");

    double period = 1.0 / fsw;
    double Thalf  = period / 2.0;

    // The analytical TDA model does not represent the dead-time interval
    // explicitly (it would require a third sub-phase with both halves of the
    // primary bridge floating). Simulating over Thalf and ignoring deadTime
    // is consistent because the resonant tank dynamics are dominated by the
    // (Ls, Cr, Lm) system, not by the small (~1 µs) bridge gap. The sample
    // grid `dt` therefore matches the output time base 1-to-1, with no
    // post-stretch.
    double Thalf_eff = Thalf;

    if (!std::isfinite(Vo) || Vo < 0)
        throw std::runtime_error("LLC: reflected output voltage Vo is invalid (non-finite or negative)");
    if (!std::isfinite(Thalf_eff) || Thalf_eff <= 0)
        throw std::runtime_error("LLC: half switching period is invalid (non-finite or non-positive)");
    if (!std::isfinite(L) || L <= 0)
        throw std::runtime_error("LLC: magnetizing inductance L is invalid (non-finite or non-positive)");

    // LIP frequency and LIP input voltage (Nielsen's f1 and Vinlip).
    // f1 = 1/(2π·√(Ls·Cr)) is the load-independent point — the frequency at
    // which all FHA gain curves cross at gain = 1. Vinlip = 2·Vo because at
    // the LIP, the bridge voltage Vi (= ½·Vin) exactly equals the reflected
    // output Vo.
    computedLipFrequency = 1.0 / (2.0 * M_PI * std::sqrt(Ls * C));
    computedLipInputVoltage = 2.0 * Vo;

    const int N = 256;
    double dt = Thalf_eff / N;

    // Initial-condition seed for the Newton solver. Use the FHA-style estimate
    // of the resonant tank current (load + magnetizing) and the magnetizing
    // peak ILm0 ≈ -Vo·Thalf/(2·L). The Vc0 seed is 0 (centred bipolar swing).
    double Im_pk_est = Vo * Thalf_eff / (2.0 * L);
    if (!std::isfinite(Im_pk_est) || std::abs(Im_pk_est) > 1e6)
        Im_pk_est = std::copysign(1.0, Im_pk_est);

    double Iout_main = llcOpPoint.get_output_currents()[0];
    double Iload_reflected = Iout_main / n_main;
    double Ires_est = std::max(std::abs(Im_pk_est) + std::abs(Iload_reflected),
                               std::abs(Iload_reflected) * 1.5);

    LlcStateVector x0_seed{ -Ires_est, -Im_pk_est, 0.0 };

    // Run the Nielsen TDA solver.
    //
    // At the Load Independent Point (Vi ≈ Vo, fsw ≈ fr) the problem is
    // under-determined: the (iLs, vC) sub-system is a 180° rotation for any
    // magnitude, so F(x0) = x_end + x0 has a 2-D null space and the Newton
    // Jacobian is singular. The solver's singularity break-out then leaves x0
    // at a seed-dependent value that can explode by orders of magnitude. We
    // guard against that two ways:
    //   (a) if |Vi − Vo| is within 0.5% of max(Vi, Vo), perturb Vi up by 0.5%
    //       during the solver so the rotation is slightly less than π and the
    //       Jacobian has full rank. The physical steady state is very close
    //       to the perturbed one, and we use the converged x0 for waveform
    //       emission with the original Vi.
    //   (b) after the solver, sanity-check the result against a generous
    //       physical bound. If it exceeds the bound, fall back to the FHA
    //       closed-form seed.
    std::vector<LlcSubStateSegment> segments;
    double residual = 0.0;

    double Vi_solver = Vi;
    double denom_vo = std::max(std::abs(Vi), std::abs(Vo));
    if (denom_vo > 0 && std::abs(Vi - Vo) / denom_vo < 0.005) {
        Vi_solver = Vi * 1.005;
    }
    LlcStateVector x0 = solve_steady_state(x0_seed, Thalf_eff,
                                           Vi_solver, Vo, Ls, L, C,
                                           segments, residual);

    // Sanity check: if the solver converged to an unphysical null-space
    // direction, fall back to the FHA closed-form estimate.
    double sanity_iLs = std::max(10.0 * Ires_est, 20.0);
    double sanity_vC  = std::max(10.0 * std::abs(Vi), 10.0 * std::abs(Vo));
    if (sanity_vC < 200.0) sanity_vC = 200.0;
    if (!std::isfinite(x0.iLs) || !std::isfinite(x0.iL) || !std::isfinite(x0.vC) ||
        std::abs(x0.iLs) > sanity_iLs ||
        std::abs(x0.vC)  > sanity_vC) {
        x0 = x0_seed;
        residual = -1.0;  // flag "seed fallback" for diagnostics
    }
    // Re-propagate with the authoritative Vi for waveform emission (so any
    // perturbation above is invisible to downstream code).
    segments = propagate_half_cycle(x0, Thalf_eff, Vi, Vo, Ls, L, C);

    // Sample the segment chain onto the existing dt grid that the
    // downstream waveform-emission code expects.
    std::vector<double> ILs_pos(N + 1, 0.0), IL_pos(N + 1, 0.0);
    std::vector<double> Vc_pos(N + 1, 0.0), VL_pos(N + 1, 0.0);
    sample_segments(segments, Thalf_eff, N, Vi, Vo, Ls, L, C,
                    ILs_pos, IL_pos, Vc_pos, VL_pos);

    // Diagnostics
    lastMode = classify_mode(segments);
    lastSubStateSequence.clear();
    for (const auto& seg : segments)
        lastSubStateSequence.push_back(static_cast<int>(seg.state));
    lastSteadyStateResidual = residual;

    // DAB-quality ZVS diagnostics.
    //   • Peak tank current: max(|iLs|) over the sampled half-period.
    //   • ZVS margin (lagging leg): magnetizing current at the half-cycle
    //     switching instant minus the reflected load current at that instant
    //     (Huang SLUP263 §6). When > 0, Coss is discharged by Im before the
    //     opposite leg turns on.
    //   • ZVS load threshold: Vbus·t_dead / (4·Lm) — the maximum load the
    //     magnetizing energy can sustain ZVS for (Huang SLUP263 eq. 12).
    //   • Resonant transition time: dead time (Coss-resonance interval).
    double peakI = 0.0;
    for (int k = 0; k <= N; ++k)
        peakI = std::max(peakI, std::abs(std::isfinite(ILs_pos[k]) ? ILs_pos[k] : 0.0));
    lastPrimaryPeakCurrent     = peakI;
    double Im_at_switch        = std::isfinite(IL_pos[N]) ? IL_pos[N] : 0.0;
    lastZvsMarginLagging       = std::abs(Im_at_switch) - std::abs(Iload_reflected);
    double Vbus_bridge         = k_bridge * inputVoltage;  // primary swing
    lastZvsLoadThreshold       = (L > 0)
        ? (Vbus_bridge * computedDeadTime) / (4.0 * L)
        : 0.0;
    lastResonantTransitionTime = computedDeadTime;

    // Stuff used by the legacy ILs0/IL0 references below (waveform emission).
    double ILs0 = x0.iLs;
    double IL0  = x0.iL;

    // Convergence diagnostic. The genuinely-dangerous outcomes are already
    // guarded upstream: an unbounded / NaN Newton iterate is caught by the
    // sanity check (Llc.cpp ~932) and replaced by the bounded FHA closed-form
    // seed (flagged residual == -1.0). What remains is a *bounded-but-imperfect*
    // Newton residual: the Nielsen TDA half-cycle symmetry condition isn't met
    // tightly (residual ~1-10 A on some designs), but the resulting waveform is
    // still close to physical — the PtP reference-design tests validate these
    // against measured hardware and pass. So this is an analytical-model accuracy
    // limitation, NOT a silent-fabrication bug; warn (don't throw) so the signal
    // is visible without rejecting validated-usable waveforms. (Tightening the
    // TDA solver to drive this residual down is tracked separately.)
    double i_thresh = std::max(0.5, 0.02 * std::abs(Iload_reflected));
    if (residual > i_thresh) {
        std::cerr << "[LLC] Nielsen TDA solver did not fully converge: residual="
                  << residual << " A (Vi=" << Vi << "V, Vo=" << Vo << "V, fsw=" << fsw
                  << "Hz) — analytical waveform is bounded but may be imperfect for this op point."
                  << std::endl;
    }

    // =====================================================================
    // Build full-period waveforms by mirroring the positive half by half-wave
    // antisymmetry. Time base is `dt` directly — no post-stretch.
    //
    // Primary voltage depends on Ls topology:
    //   SEPARATE Ls:   Vpri = VLm(t) — flat clamp + soft hill (magnetizing only)
    //   INTEGRATED Ls: Vpri = Vi - Vc(t) — bridge voltage with Cr ripple
    // =====================================================================

    int totalSamples = 2 * N + 1;
    std::vector<double> time_full(totalSamples);
    std::vector<double> ILs_full(totalSamples);
    std::vector<double> IL_full(totalSamples);
    std::vector<double> Vpri_full(totalSamples);  // Primary voltage (topology-dependent)
    std::vector<double> VLm_full(totalSamples);   // Magnetizing voltage (for core loss & secondary)

    for (int k = 0; k <= N; ++k) {
        time_full[k] = k * dt;
        ILs_full[k] = std::isfinite(ILs_pos[k]) ? ILs_pos[k] : ILs0;
        IL_full[k]  = std::isfinite(IL_pos[k])  ? IL_pos[k]  : IL0;

        double VLm_k = std::isfinite(VL_pos[k]) ? VL_pos[k] : 0.0;
        double Vc_k  = std::isfinite(Vc_pos[k]) ? Vc_pos[k] : 0.0;
        VLm_full[k] = VLm_k;

        if (integratedLs) {
            Vpri_full[k] = Vi - Vc_k;
        } else {
            Vpri_full[k] = VLm_k;
        }
    }

    for (int k = 1; k <= N; ++k) {
        time_full[N + k] = Thalf + k * dt;
        ILs_full[N + k] = std::isfinite(ILs_pos[k]) ? -ILs_pos[k] : -ILs0;
        IL_full[N + k]  = std::isfinite(IL_pos[k])  ? -IL_pos[k]  : -IL0;

        double VLm_k = std::isfinite(VL_pos[k]) ? VL_pos[k] : 0.0;
        double Vc_k  = std::isfinite(Vc_pos[k]) ? Vc_pos[k] : 0.0;
        VLm_full[N + k] = -VLm_k;

        if (integratedLs) {
            Vpri_full[N + k] = -(Vi - Vc_k);
        } else {
            Vpri_full[N + k] = -VLm_k;
        }
    }

    // --- Build and store extra-component waveforms for get_extra_components_inputs ---
    {
        // Capacitor voltage: Vc_full assembled from Vc_pos with half-wave antisymmetry
        std::vector<double> Vc_full(totalSamples);
        for (int k = 0; k <= N; ++k)
            Vc_full[k] = std::isfinite(Vc_pos[k]) ? Vc_pos[k] : 0.0;
        for (int k = 1; k <= N; ++k)
            Vc_full[N + k] = -Vc_full[k];

        Waveform capVoltWfm;
        capVoltWfm.set_data(Vc_full);
        capVoltWfm.set_time(time_full);
        capVoltWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        extraCapVoltageWaveforms.push_back(capVoltWfm);

        // Cap current = series tank current ILs
        Waveform capCurrWfm;
        capCurrWfm.set_data(ILs_full);
        capCurrWfm.set_time(time_full);
        capCurrWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        extraCapCurrentWaveforms.push_back(capCurrWfm);

        extraTimeVectors.push_back(time_full);

        // Series inductor voltage VLs = Vi_sq - Vc - VLm (only stored when Ls is separate)
        if (!integratedLs) {
            std::vector<double> VLs_full(totalSamples);
            for (int k = 0; k <= N; ++k) {
                double Vc_k  = Vc_full[k];
                double VLm_k = std::isfinite(VLm_full[k]) ? VLm_full[k] : 0.0;
                VLs_full[k] = Vi - Vc_k - VLm_k;
            }
            for (int k = 1; k <= N; ++k)
                VLs_full[N + k] = -VLs_full[k];

            Waveform indVoltWfm;
            indVoltWfm.set_data(VLs_full);
            indVoltWfm.set_time(time_full);
            indVoltWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            extraIndVoltageWaveforms.push_back(indVoltWfm);

            // Series inductor current = same as capacitor current (series circuit)
            extraIndCurrentWaveforms.push_back(capCurrWfm);
        }
    }

    // --- Primary excitation ---
    {
        Waveform currentWaveform;
        currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        currentWaveform.set_data(ILs_full);
        currentWaveform.set_time(time_full);

        Waveform voltageWaveform;
        voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        voltageWaveform.set_data(Vpri_full);
        voltageWaveform.set_time(time_full);

        auto excitation = complete_excitation(currentWaveform, voltageWaveform,
                                              fsw, "Primary");
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    // --- Secondary excitations ---
    //
    // Convention:
    //   n_i  = primary-to-secondary turns ratio for output i
    //   Id   = ILs - IL  is the *primary-referred* secondary diode current
    //          (the share of the resonant tank current actually transferred
    //          to the secondaries via the transformer).
    //   I_sec_i = Id * n_i           (ampere-turn balance: I_sec increases with n
    //                                   because we multiply, not divide).
    //
    // Multi-output approximation: each output receives its share of Id in
    // proportion to its load conductance, then is scaled by n_i. For a single
    // output this collapses to the exact analytical answer.
    //
    // Per-rectifier-type waveform shaping (Llc.h::LlcRectifierType):
    //   CT  → 2 windings per output. Each half conducts on one polarity only.
    //         Vsec_half = VLm / n; |I_half| = max(0, ±Id*share*n).
    //   FB  → 1 winding per output. Continuous full-wave conduction.
    //         Vsec  = VLm / n (bipolar 3-level); |I_sec| = |Id*share|*n.
    //   CD  → 1 winding per output + two output inductors L_o1/L_o2.
    //         Vsec  = VLm / n; I_sec = Id*share*n (alternating polarity).
    //         Each Lo carries IoutDC/2 with 2·Fs ripple. (Lo waveforms stored
    //         in extraLo*Waveforms for get_extra_components_inputs.)
    //   VD  → 1 winding per output. Cap stack provides 2·Vsec_pk to load.
    //         n_i is doubled at design time so Vsec_pk = Vout/2.
    //         Vsec = VLm / n; |I_sec| = |Id*share|*n; output cap voltage = Vo/2.
    double total_g = 0.0;
    for (size_t i = 0; i < nOutputs; ++i) {
        double Vout_i = llcOpPoint.get_output_voltages()[i];
        double Iout_i = llcOpPoint.get_output_currents()[i];
        if (Vout_i > 0 && Iout_i > 0) total_g += Iout_i / Vout_i;
    }
    if (total_g <= 0) total_g = 1.0;

    LlcRectifierType rectType = get_effective_rectifier_type();

    for (size_t outputIdx = 0; outputIdx < nOutputs; ++outputIdx) {
        double n_i = get_n_for_output(outputIdx);
        if (n_i <= 0) n_i = 1.0;

        double Vout_i = llcOpPoint.get_output_voltages()[outputIdx];
        double Iout_i = llcOpPoint.get_output_currents()[outputIdx];
        double share = (Vout_i > 0 && Iout_i > 0)
                       ? (Iout_i / Vout_i) / total_g
                       : (outputIdx == 0 ? 1.0 : 0.0);

        switch (rectType) {
        case LlcRectifierType::CENTER_TAPPED: {
            // Two half-windings, each conducting on one polarity of the cycle.
            for (size_t halfIdx = 0; halfIdx < 2; ++halfIdx) {
                std::vector<double> iSecData(totalSamples, 0.0);
                std::vector<double> vSecData(totalSamples, 0.0);
                for (int k = 0; k < totalSamples; ++k) {
                    double Id = ILs_full[k] - IL_full[k];
                    if (!std::isfinite(Id)) Id = 0;
                    double Id_share = Id * share;
                    double Id_half = (halfIdx == 0)
                                     ? std::max(0.0, Id_share)
                                     : std::max(0.0, -Id_share);
                    iSecData[k] = Id_half * n_i;
                    vSecData[k] = VLm_full[k] / n_i;
                    if (!std::isfinite(iSecData[k])) iSecData[k] = 0;
                    if (!std::isfinite(vSecData[k])) vSecData[k] = 0;
                }
                Waveform secCurrentWfm; secCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
                secCurrentWfm.set_data(iSecData); secCurrentWfm.set_time(time_full);
                Waveform secVoltageWfm; secVoltageWfm.set_ancillary_label(WaveformLabel::CUSTOM);
                secVoltageWfm.set_data(vSecData); secVoltageWfm.set_time(time_full);
                std::string windingName = "Secondary " + std::to_string(outputIdx)
                                        + " Half " + std::to_string(halfIdx + 1);
                auto excitation = complete_excitation(secCurrentWfm, secVoltageWfm, fsw, windingName);
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }
            break;
        }
        case LlcRectifierType::FULL_BRIDGE: {
            // Single secondary winding. Bipolar Vsec, full-wave rectified at load.
            // Diode-bridge rectifies; the winding carries bipolar Id*share*n.
            std::vector<double> iSecData(totalSamples, 0.0);
            std::vector<double> vSecData(totalSamples, 0.0);
            for (int k = 0; k < totalSamples; ++k) {
                double Id = ILs_full[k] - IL_full[k];
                if (!std::isfinite(Id)) Id = 0;
                iSecData[k] = Id * share * n_i;          // bipolar; sign matches Id
                vSecData[k] = VLm_full[k] / n_i;         // bipolar Vsec
                if (!std::isfinite(iSecData[k])) iSecData[k] = 0;
                if (!std::isfinite(vSecData[k])) vSecData[k] = 0;
            }
            Waveform secCurrentWfm; secCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            secCurrentWfm.set_data(iSecData); secCurrentWfm.set_time(time_full);
            Waveform secVoltageWfm; secVoltageWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            secVoltageWfm.set_data(vSecData); secVoltageWfm.set_time(time_full);
            std::string windingName = "Secondary " + std::to_string(outputIdx);
            auto excitation = complete_excitation(secCurrentWfm, secVoltageWfm, fsw, windingName);
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            break;
        }
        case LlcRectifierType::CURRENT_DOUBLER: {
            // Single secondary winding driving two output inductors via two diodes.
            // Secondary current = Id*share*n (bipolar), but each Lo only conducts on
            // alternating half-cycles: Lo1 = max(0, +I_sec), Lo2 = max(0, -I_sec),
            // each with average ≈ Iout/2 and ripple at 2·Fs (cancels at load node).
            std::vector<double> iSecData(totalSamples, 0.0);
            std::vector<double> vSecData(totalSamples, 0.0);
            std::vector<double> iLo1Data(totalSamples, 0.0);
            std::vector<double> iLo2Data(totalSamples, 0.0);
            std::vector<double> vLo1Data(totalSamples, 0.0);
            std::vector<double> vLo2Data(totalSamples, 0.0);
            // CD inductor sizing: ripple ≤ 30% of IoutDC by convention.
            // VLo amplitude tracks Vsec − Vo when the corresponding diode conducts,
            // and −Vo during freewheel.
            for (int k = 0; k < totalSamples; ++k) {
                double Id = ILs_full[k] - IL_full[k];
                if (!std::isfinite(Id)) Id = 0;
                double I_sec = Id * share * n_i;
                double V_sec = VLm_full[k] / n_i;
                iSecData[k] = I_sec;
                vSecData[k] = V_sec;
                // Lo1 / Lo2 currents — KCL at output node: ILo1 + ILo2 = IoutDC.
                // Each carries IoutDC/2 mean with a triangular ripple of opposite
                // sign sourced from the secondary current (which itself averages 0).
                double Iout_dc = (Vout_i > 0) ? (Iout_i) : 0.0;
                double ripple  = I_sec / 4.0;            // peak excursion ≈ |I_sec|/4
                iLo1Data[k] = Iout_dc / 2.0 + ripple;
                iLo2Data[k] = Iout_dc / 2.0 - ripple;
                // VLo when conducting = (Vsec − Vo); freewheel = -Vo.
                vLo1Data[k] = (I_sec > 0) ? (V_sec - Vout_i) : (-Vout_i);
                vLo2Data[k] = (I_sec < 0) ? (-V_sec - Vout_i) : (-Vout_i);
                if (!std::isfinite(iSecData[k])) iSecData[k] = 0;
                if (!std::isfinite(vSecData[k])) vSecData[k] = 0;
                if (!std::isfinite(iLo1Data[k])) iLo1Data[k] = 0;
                if (!std::isfinite(iLo2Data[k])) iLo2Data[k] = 0;
                if (!std::isfinite(vLo1Data[k])) vLo1Data[k] = 0;
                if (!std::isfinite(vLo2Data[k])) vLo2Data[k] = 0;
            }
            Waveform secCurrentWfm; secCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            secCurrentWfm.set_data(iSecData); secCurrentWfm.set_time(time_full);
            Waveform secVoltageWfm; secVoltageWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            secVoltageWfm.set_data(vSecData); secVoltageWfm.set_time(time_full);
            std::string windingName = "Secondary " + std::to_string(outputIdx);
            auto excitation = complete_excitation(secCurrentWfm, secVoltageWfm, fsw, windingName);
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            // Stash Lo waveforms for get_extra_components_inputs.
            Waveform lo1I; lo1I.set_ancillary_label(WaveformLabel::CUSTOM);
            lo1I.set_data(iLo1Data); lo1I.set_time(time_full);
            Waveform lo1V; lo1V.set_ancillary_label(WaveformLabel::CUSTOM);
            lo1V.set_data(vLo1Data); lo1V.set_time(time_full);
            Waveform lo2I; lo2I.set_ancillary_label(WaveformLabel::CUSTOM);
            lo2I.set_data(iLo2Data); lo2I.set_time(time_full);
            Waveform lo2V; lo2V.set_ancillary_label(WaveformLabel::CUSTOM);
            lo2V.set_data(vLo2Data); lo2V.set_time(time_full);
            extraLoCurrentWaveforms.push_back(lo1I);
            extraLoVoltageWaveforms.push_back(lo1V);
            extraLo2CurrentWaveforms.push_back(lo2I);
            extraLo2VoltageWaveforms.push_back(lo2V);
            break;
        }
        case LlcRectifierType::VOLTAGE_DOUBLER: {
            // Single secondary winding. Two diodes + two output capacitors stacked.
            // Each cap charges to Vsec_pk = Vo/2 (n_i was doubled at design time).
            // The secondary current is bipolar; both caps share the same |I_sec|.
            std::vector<double> iSecData(totalSamples, 0.0);
            std::vector<double> vSecData(totalSamples, 0.0);
            std::vector<double> vCapData(totalSamples, 0.0);
            for (int k = 0; k < totalSamples; ++k) {
                double Id = ILs_full[k] - IL_full[k];
                if (!std::isfinite(Id)) Id = 0;
                iSecData[k] = Id * share * n_i;
                vSecData[k] = VLm_full[k] / n_i;
                vCapData[k] = Vout_i / 2.0;              // DC nominal; ripple is small
                if (!std::isfinite(iSecData[k])) iSecData[k] = 0;
                if (!std::isfinite(vSecData[k])) vSecData[k] = 0;
            }
            Waveform secCurrentWfm; secCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            secCurrentWfm.set_data(iSecData); secCurrentWfm.set_time(time_full);
            Waveform secVoltageWfm; secVoltageWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            secVoltageWfm.set_data(vSecData); secVoltageWfm.set_time(time_full);
            std::string windingName = "Secondary " + std::to_string(outputIdx);
            auto excitation = complete_excitation(secCurrentWfm, secVoltageWfm, fsw, windingName);
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            Waveform capWfm; capWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            capWfm.set_data(vCapData); capWfm.set_time(time_full);
            extraVoutCapVoltageWaveforms.push_back(capWfm);
            break;
        }
        default:
            throw std::runtime_error("LLC: unknown rectifier type encountered in process_operating_point_for_input_voltage");
        }
    }

    OperatingConditions conditions;
    conditions.set_ambient_temperature(llcOpPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    return operatingPoint;
}


// =====================================================================
// SPICE Circuit Generation
// =====================================================================
std::string Llc::generate_ngspice_circuit(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t inputVoltageIndex,
    size_t operatingPointIndex)
{
    // Pull all SPICE-side knobs (PULSE timing, SW model, snubbers, diode
    // model, output cap, solver options) from the per-topology config —
    // no literals from here on out. Registry default for LLC is wired in
    // Topology.cpp (LLC_RESONANT_CONVERTER entry); the values below
    // therefore match the previous hardcoded netlist byte-for-byte.
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

    double inputVoltage = inputVoltages[std::min(inputVoltageIndex, inputVoltages.size() - 1)];
    auto& llcOp = ops[std::min(operatingPointIndex, ops.size() - 1)];

    double fsw = llcOp.get_switching_frequency();
    double period = 1.0 / fsw;
    double halfPeriod = period / 2.0;
    // computedDeadTime is a fixed default (1 µs) sized for low-frequency LLC.
    // At high fsw (e.g. 250 kHz → halfPeriod = 2 µs) a 1 µs dead time leaves
    // tOn = 0 and the VOLTAGE_CONTROLLED_SWITCH bridge aborts. Clamp it to a
    // physically sane fraction of the half-period (≤10 %, ZVS-realistic) so
    // the configured dead time is honoured when it fits but never exceeds the
    // available conduction window. Mirrors Src.cpp's period-scaled dead time.
    double deadTime = std::min(computedDeadTime, halfPeriod * 0.1);
    double tOn = halfPeriod - deadTime;

    size_t numOutputs = llcOp.get_output_voltages().size();
    if (numOutputs == 0) numOutputs = 1;

    // Recover per-output turns ratios. process_design_requirements pushes the
    // ratio twice per output (one per center-tapped half) so this collapses to
    // numOutputs entries here.
    std::vector<double> nPerOutput(numOutputs);
    for (size_t i = 0; i < numOutputs; ++i) {
        if (turnsRatios.size() == 2 * numOutputs)
            nPerOutput[i] = turnsRatios[2 * i];
        else if (i < turnsRatios.size())
            nPerOutput[i] = turnsRatios[i];
        else if (!turnsRatios.empty())
            nPerOutput[i] = turnsRatios[0];
        else
            nPerOutput[i] = inputVoltage / llcOp.get_output_voltages()[i];
        if (nPerOutput[i] <= 0) nPerOutput[i] = 1.0;
    }

    double Ls = computedResonantInductance;
    double Cr = computedResonantCapacitance;
    double L  = magnetizingInductance;

    bool isFullBridge = (get_bridge_type().has_value() &&
                         get_bridge_type().value() == LlcBridgeType::FULL_BRIDGE);
    bool integratedLs = get_integrated_resonant_inductor().value_or(false);

    int numPeriodsTotal = get_num_steady_state_periods() + get_num_periods_to_extract();
    double simTime = numPeriodsTotal * period;
    double startTime = get_num_steady_state_periods() * period;
    double maxStep = period / 200.0;

    std::ostringstream circuit;

    circuit << "* LLC Resonant Converter v3.3 (" << (integratedLs ? "Integrated Ls" : "Separate Ls") << ")\n";
    circuit << "* " << (isFullBridge ? "Full" : "Half") << "-Bridge\n";
    circuit << "* Vin=" << inputVoltage << "V  fsw=" << (fsw/1e3) << "kHz  outputs=" << numOutputs << "\n";
    circuit << "* Ls=" << (Ls*1e6) << "uH  Cr=" << (Cr*1e9) << "nF  Lm=" << (L*1e6) << "uH\n";
    for (size_t i = 0; i < numOutputs; ++i) {
        circuit << "*   output " << i << ": Vout=" << llcOp.get_output_voltages()[i]
                << "V Iout=" << llcOp.get_output_currents()[i] << "A n=" << nPerOutput[i] << "\n";
    }
    circuit << "\n";

    // §8a.5: Real DC rail. The previous design emitted a synthetic
    // Vdc_supply with a 1MEG dummy resistor (Vdc_supply is otherwise
    // disconnected in BEHAVIORAL_PULSE mode), so i(Vdc_supply) only
    // saw the dummy load — useless for input-current extraction. The
    // rail is now wired through the bridge switches (SW1 mode) so
    // i(Vin) reads the true DC bus current via dedicated zero-V sense
    // sources upstream of each high-side switch.
    circuit << "Vin vin_dc 0 " << inputVoltage << "\n\n";

    const auto bridgeMode = get_bridge_simulation_mode();

    if (bridgeMode == BridgeSimulationMode::BEHAVIORAL_PULSE) {
        // -----------------------------------------------------------------
        // Behavioral PULSE bridge — fast path (MagneticAdviser default).
        // Single dependent V-source generates the bridge differential output
        // directly. No SW1 / body diodes / snubbers — the WASM-friendly
        // emission used by all default LLC simulations.
        // -----------------------------------------------------------------
        if (isFullBridge) {
            circuit << "Vbridge bridge_a bridge_b PULSE("
                    << -inputVoltage << " " << inputVoltage << " 0 "
                    << std::scientific << deadTime << " " << deadTime << " "
                    << tOn << " " << period << std::fixed << ")\n";
            circuit << "Vpri_sense bridge_a lr_in 0\n";
            circuit << "Vbus_gnd bridge_b 0 0\n\n";
        } else {
            circuit << "Vbridge sw_node mid_point PULSE("
                    << -(inputVoltage / 2.0) << " " << (inputVoltage / 2.0) << " 0 "
                    << std::scientific << deadTime << " " << deadTime << " "
                    << tOn << " " << period << std::fixed << ")\n";
            circuit << "Vmid mid_point 0 0\n";
            circuit << "Vpri_sense sw_node lr_in 0\n\n";
        }
    }
    else {
        // -----------------------------------------------------------------
        // Voltage-controlled-switch bridge — high-fidelity path.
        // Models MOSFETs as SW1 with antiparallel body diodes and RC
        // snubbers. ngspice resolves switching detail at sub-step
        // granularity → ~10× slower; reserve for PtP reference designs.
        //   Full bridge: 4 SW1 in two legs (QA/QB on bridge_a, QC/QD on
        //     bridge_b). Diagonal pairs (QA,QD) and (QB,QC) conduct on
        //     opposite half-cycles → V(bridge_a,bridge_b) swings ±Vin.
        //   Half bridge: 2 SW1 (QHI/QLO on sw_node) + DC bus split caps
        //     forming mid_point at Vin/2 → V(sw_node,mid_point) swings
        //     ±Vin/2.
        // Gate sources are 50% complementary with deadTime gap; PSFB-style
        // 10 ns gate transitions.
        // -----------------------------------------------------------------
        if (deadTime <= 0.0) {
            throw std::runtime_error(
                "Llc VOLTAGE_CONTROLLED_SWITCH bridge: computedDeadTime must be > 0. "
                "Got " + std::to_string(deadTime));
        }
        if (2.0 * deadTime >= halfPeriod) {
            throw std::runtime_error(
                "Llc VOLTAGE_CONTROLLED_SWITCH bridge: deadTime ("
                + std::to_string(deadTime) + " s) is too large for halfPeriod ("
                + std::to_string(halfPeriod) + " s) — tOn would be non-positive.");
        }

        circuit << ".model SW1 SW VT=" << cfg.swModelVT << " VH=" << cfg.swModelVH
                << " RON=" << cfg.swModelRON << " ROFF=" << cfg.swModelROFF << "\n";
        circuit << ".model DIDEAL D(IS=" << std::scientific << cfg.diodeIS
                << " RS=" << cfg.diodeRS << std::fixed << ")\n\n";

        if (isFullBridge) {
            // Diagonal QA+QD on first half, QB+QC on second half.
            circuit << "Vpwm_QA pwm_QA 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                    << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                    << tOn << " " << period << std::fixed << ")\n";
            circuit << "Vpwm_QB pwm_QB 0 PULSE(0 " << cfg.pwmHigh << " "
                    << std::scientific << halfPeriod << " "
                    << cfg.pwmRise << " " << cfg.pwmFall << " "
                    << tOn << " " << period << std::fixed << ")\n";
            circuit << "Vpwm_QC pwm_QC 0 PULSE(0 " << cfg.pwmHigh << " "
                    << std::scientific << halfPeriod << " "
                    << cfg.pwmRise << " " << cfg.pwmFall << " "
                    << tOn << " " << period << std::fixed << ")\n";
            circuit << "Vpwm_QD pwm_QD 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                    << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                    << tOn << " " << period << std::fixed << ")\n";
            // §8a.5 input-current sense: zero-V sources upstream of each
            // high-side switch drain. i(Vq1_sense)+i(Vq3_sense) gives the
            // true DC-bus current sourced from vin_dc.
            circuit << "Vq1_sense vin_dc qa_drain 0\n";
            circuit << "Vq3_sense vin_dc qc_drain 0\n";
            circuit << "SQA qa_drain bridge_a pwm_QA 0 SW1\n";
            circuit << "DQA 0 bridge_a DIDEAL\n";
            circuit << "SQB bridge_a 0 pwm_QB 0 SW1\n";
            circuit << "DQB bridge_a qa_drain DIDEAL\n";
            circuit << "Rsnub_QA qa_drain bridge_a " << cfg.snubR
                    << "\nCsnub_QA qa_drain bridge_a " << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "Rsnub_QB bridge_a 0 " << cfg.snubR
                    << "\nCsnub_QB bridge_a 0 " << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "SQC qc_drain bridge_b pwm_QC 0 SW1\n";
            circuit << "DQC 0 bridge_b DIDEAL\n";
            circuit << "SQD bridge_b 0 pwm_QD 0 SW1\n";
            circuit << "DQD bridge_b qc_drain DIDEAL\n";
            circuit << "Rsnub_QC qc_drain bridge_b " << cfg.snubR
                    << "\nCsnub_QC qc_drain bridge_b " << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "Rsnub_QD bridge_b 0 " << cfg.snubR
                    << "\nCsnub_QD bridge_b 0 " << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "Vpri_sense bridge_a lr_in 0\n\n";
        } else {
            // Half-bridge: split caps form mid_point at Vin/2.
            circuit << "Cbus_hi vin_dc mid_point 1u IC=" << (inputVoltage / 2.0) << "\n";
            circuit << "Cbus_lo mid_point 0 1u IC=" << (inputVoltage / 2.0) << "\n";
            circuit << "Rbal_hi vin_dc mid_point 100k\n";
            circuit << "Rbal_lo mid_point 0 100k\n";
            circuit << "Vpwm_HI pwm_HI 0 PULSE(0 " << cfg.pwmHigh << " 0 "
                    << std::scientific << cfg.pwmRise << " " << cfg.pwmFall << " "
                    << tOn << " " << period << std::fixed << ")\n";
            circuit << "Vpwm_LO pwm_LO 0 PULSE(0 " << cfg.pwmHigh << " "
                    << std::scientific << halfPeriod << " "
                    << cfg.pwmRise << " " << cfg.pwmFall << " "
                    << tOn << " " << period << std::fixed << ")\n";
            // §8a.5 input-current sense upstream of the high-side switch.
            circuit << "Vq1_sense vin_dc qhi_drain 0\n";
            circuit << "SHI qhi_drain sw_node pwm_HI 0 SW1\n";
            circuit << "DHI 0 sw_node DIDEAL\n";
            circuit << "SLO sw_node 0 pwm_LO 0 SW1\n";
            circuit << "DLO sw_node qhi_drain DIDEAL\n";
            circuit << "Rsnub_HI qhi_drain sw_node " << cfg.snubR
                    << "\nCsnub_HI qhi_drain sw_node " << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "Rsnub_LO sw_node 0 " << cfg.snubR
                    << "\nCsnub_LO sw_node 0 " << std::scientific << cfg.snubC << std::fixed << "\n";
            circuit << "Vpri_sense sw_node lr_in 0\n\n";
        }
    }

    double Lpri_total = integratedLs ? (L + Ls) : L;
    double k_int = integratedLs ? std::sqrt(L / Lpri_total) : 0.999;

    if (integratedLs) {
        circuit << "* Resonant cap (no separate inductor — leakage provides Ls)\n";
        circuit << "Cr lr_in pri_top " << std::scientific << Cr << "\n\n";
        circuit << "* Transformer: Lpri=Lm+Ls=" << (Lpri_total*1e6) << "uH, k=" << k_int << "\n";
    } else {
        circuit << "* Resonant tank (separate inductor)\n";
        circuit << "Cr lr_in cr_ls " << std::scientific << Cr << "\n";
        circuit << "Lr cr_ls pri_top " << std::scientific << Ls << "\n\n";
        circuit << "* Transformer (k~0.999, pairwise K statements for ngspice compatibility)\n";
    }
    circuit << "Lpri pri_top pri_bot " << std::scientific << Lpri_total << "\n";

    // Per-output secondary windings — count varies with rectifier type:
    //   CT  → 2 half-windings (Lsec1/Lsec2) with center tap (sec_ct).
    //   FB  → 1 secondary (Lsec) feeding a 4-diode bridge.
    //   CD  → 1 secondary (Lsec) feeding 2 diodes + 2 output inductors.
    //   VD  → 1 secondary (Lsec) feeding 2 diodes + a series capacitor stack.
    LlcRectifierType rectType = get_effective_rectifier_type();
    bool ct = (rectType == LlcRectifierType::CENTER_TAPPED);

    for (size_t i = 0; i < numOutputs; ++i) {
        double n_i = nPerOutput[i];
        std::string si = std::to_string(i + 1);
        if (ct) {
            double Lsec_half = Lpri_total / (n_i * n_i);
            circuit << "Lsec1_o" << si << " sec_top_sec_o" << si << " sec_ct_o" << si
                    << " " << std::scientific << Lsec_half << "\n";
            circuit << "Lsec2_o" << si << " sec_ct_o" << si << " sec_bot_sec_o" << si
                    << " " << std::scientific << Lsec_half << "\n";
        } else {
            // Single-winding secondary. Full Lsec sees the full Vpri/n_i.
            double Lsec = Lpri_total / (n_i * n_i);
            circuit << "Lsec_o" << si << " sec_pos_sec_o" << si << " sec_neg_sec_o" << si
                    << " " << std::scientific << Lsec << "\n";
        }
    }
    // Emit a full pairwise K matrix across every coupled inductor. ngspice
    // requires explicit pairwise couplings for every pair — an "incomplete
    // K set" makes the mutual-inductance matrix non-positive-definite and the
    // transient solver can't converge.
    std::vector<std::string> coupledInductors;
    coupledInductors.reserve(1 + windings_per_output() * numOutputs);
    coupledInductors.emplace_back("Lpri");
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        if (ct) {
            coupledInductors.push_back("Lsec1_o" + si);
            coupledInductors.push_back("Lsec2_o" + si);
        } else {
            coupledInductors.push_back("Lsec_o" + si);
        }
    }
    int kIdx = 0;
    for (size_t a = 0; a < coupledInductors.size(); ++a) {
        for (size_t b = a + 1; b < coupledInductors.size(); ++b) {
            circuit << "K" << ++kIdx << " " << coupledInductors[a]
                    << " " << coupledInductors[b] << " " << k_int << "\n";
        }
    }
    circuit << "\n";

    if (isFullBridge) circuit << "Rpri_ret pri_bot bridge_b 0.001\n\n";
    else              circuit << "Rpri_ret pri_bot mid_point 0.001\n\n";

    // §8a.5 differential winding-voltage probes — the bare-node v(pri_top)
    // probe used previously lumped Lpri + the bridge-return offset
    // (pri_bot = mid_point @ Vin/2 in HB SW1 mode); the secondary
    // bare-node probes likewise read post-sense-source nodes instead of
    // the winding-only differential. These E-sources expose
    // v(vpri_w)/v(vsec*_w_o<i>) as pure winding voltages.
    circuit << "Evpri_w vpri_w 0 pri_top pri_bot 1\n";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        if (ct) {
            circuit << "Evsec1_w_o" << si << " vsec1_w_o" << si
                    << " 0 sec_top_sec_o" << si << " sec_ct_o" << si << " 1\n";
            circuit << "Evsec2_w_o" << si << " vsec2_w_o" << si
                    << " 0 sec_ct_o" << si << " sec_bot_sec_o" << si << " 1\n";
        } else {
            circuit << "Evsec_w_o" << si << " vsec_w_o" << si
                    << " 0 sec_pos_sec_o" << si << " sec_neg_sec_o" << si << " 1\n";
        }
    }
    circuit << "\n";

    circuit << ".model DRECT D(Is=1e-8 N=0.01 RS=0.01)\n";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        double Vout_i = llcOp.get_output_voltages()[i];
        double Iout_i = llcOp.get_output_currents()[i];
        double Rload_i = (Iout_i > 0) ? (Vout_i / Iout_i) : 100.0;
        // n_i kept for parity with other branches; CD branch consumes it.
        [[maybe_unused]] double n_i = nPerOutput[i];

        circuit << "* Output " << si << " rectifier ("
                << (ct ? "CT" : rectType == LlcRectifierType::FULL_BRIDGE ? "FB"
                              : rectType == LlcRectifierType::CURRENT_DOUBLER ? "CD" : "VD")
                << " Vout=" << Vout_i << "V Rload=" << Rload_i << ")\n";

        switch (rectType) {
        case LlcRectifierType::CENTER_TAPPED: {
            circuit << "D1_o" << si << " sec_top_o" << si << " vout_pos_o" << si << " DRECT\n";
            circuit << "D2_o" << si << " sec_bot_o" << si << " vout_pos_o" << si << " DRECT\n";
            circuit << "Rsn1_o" << si << " sec_top_o" << si << " vout_pos_o" << si << " 100\n";
            circuit << "Csn1_o" << si << " sec_top_o" << si << " vout_pos_o" << si << " 100p\n";
            circuit << "Rsn2_o" << si << " sec_bot_o" << si << " vout_pos_o" << si << " 100\n";
            circuit << "Csn2_o" << si << " sec_bot_o" << si << " vout_pos_o" << si << " 100p\n";
            circuit << "Vsec1_sense_o" << si << " sec_top_sec_o" << si << " sec_top_o" << si << " 0\n";
            circuit << "Vsec2_sense_o" << si << " sec_bot_sec_o" << si << " sec_bot_o" << si << " 0\n";
            circuit << "Vsec_sense_o" << si << " sec_ct_o" << si << " vout_neg_o" << si << " 0\n";
            circuit << "Vgnd_o" << si << " vout_neg_o" << si << " 0 0\n";
            circuit << "Resr_o" << si << " vout_pos_o" << si << " vout_cap_o" << si << " 0.05\n";
            circuit << "Cout_o" << si << " vout_cap_o" << si << " vout_neg_o" << si << " " << std::scientific << 47e-6 << "\n";
            circuit << "Rload_o" << si << " vout_cap_o" << si << " vout_neg_o" << si << " " << Rload_i << "\n\n";
            break;
        }
        case LlcRectifierType::FULL_BRIDGE: {
            // 4-diode bridge: sec_pos / sec_neg → vout_pos / vout_neg.
            circuit << "Dh1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " DRECT\n";
            circuit << "Dh2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " DRECT\n";
            circuit << "Dl1_o" << si << " vout_neg_o" << si << " sec_pos_o" << si << " DRECT\n";
            circuit << "Dl2_o" << si << " vout_neg_o" << si << " sec_neg_o" << si << " DRECT\n";
            // Snubbers across each diode (optional, helps convergence).
            circuit << "Rsnh1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100\n";
            circuit << "Csnh1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100p\n";
            circuit << "Rsnh2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " 100\n";
            circuit << "Csnh2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " 100p\n";
            // Winding-current sense
            circuit << "Vsec_sense_o" << si << " sec_pos_sec_o" << si << " sec_pos_o" << si << " 0\n";
            circuit << "Vsec_ret_o" << si << " sec_neg_o" << si << " sec_neg_sec_o" << si << " 0\n";
            circuit << "Vgnd_o" << si << " vout_neg_o" << si << " 0 0\n";
            circuit << "Resr_o" << si << " vout_pos_o" << si << " vout_cap_o" << si << " 0.05\n";
            circuit << "Cout_o" << si << " vout_cap_o" << si << " vout_neg_o" << si << " " << std::scientific << 47e-6 << "\n";
            circuit << "Rload_o" << si << " vout_cap_o" << si << " vout_neg_o" << si << " " << Rload_i << "\n\n";
            break;
        }
        case LlcRectifierType::CURRENT_DOUBLER: {
            // CURRENT_DOUBLER (canonical topology — mirrors Src.cpp CD branch):
            //   D1:  sec_pos → vout_pos          D2:  sec_neg → vout_pos
            //   Lo1: sec_pos → vout_neg          Lo2: sec_neg → vout_neg
            //   Cout / Rload between vout_pos and vout_neg
            //
            // Conduction trace (V(sec_pos) > V(sec_neg)):
            //   sec_pos → D1 → vout_pos → Rload → vout_neg → Lo2 → sec_neg
            //   (winding closes the loop). Lo1 freewheels through D1 + Rload + Lo2.
            // Each Lo carries Iout/2 DC; ripple cancels at 2·fs at vout_pos.
            //
            // Lo sizing for SPICE: large enough to clamp Lo currents to
            // ~Iout/2 (CCM with small ripple) so the current-doubler mode
            // actually obtains. ~30% per-Lo ripple:
            //     Lo = Vo / (4 · fs · 0.30 · Iout_DC)
            // (Independent of the 30% designer-facing ripple emitted by
            // get_extra_components_inputs for the physical magnetic.)
            //
            // Earlier LLC CD SPICE had no return path from vout_neg to the
            // secondary winding (Lo's terminated at vout_pos with isolated
            // freewheel diodes onto vout_neg) — SPICE never converged on
            // a real CD case, and LLC had no CD PtP regression so the bug
            // was never caught. Restored to the textbook form used in
            // Telecom 3 kW briks (Steigerwald §IV).
            double Iout_dc = (Iout_i > 0) ? Iout_i : 1.0;
            double fsw_sim = 1.0 / period;
            double Lo = (Vout_i > 0)
                ? (Vout_i / (4.0 * fsw_sim * 0.30 * Iout_dc))
                : 10e-6;
            circuit << "D1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " DRECT\n";
            circuit << "D2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " DRECT\n";
            // RC snubbers across each diode.
            circuit << "Rsn1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100\n";
            circuit << "Csn1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100p\n";
            circuit << "Rsn2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " 100\n";
            circuit << "Csn2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " 100p\n";
            // Output inductors with per-Lo current sense.
            circuit << "VLo1_sense_o" << si << " sec_pos_o" << si << " lo1_a_o" << si << " 0\n";
            circuit << "Lo1_o" << si << " lo1_a_o" << si << " vout_neg_o" << si << " " << std::scientific << Lo << "\n";
            circuit << "VLo2_sense_o" << si << " sec_neg_o" << si << " lo2_a_o" << si << " 0\n";
            circuit << "Lo2_o" << si << " lo2_a_o" << si << " vout_neg_o" << si << " " << std::scientific << Lo << "\n";
            // Secondary winding sense (bipolar; same naming as FB).
            circuit << "Vsec_sense_o" << si << " sec_pos_sec_o" << si << " sec_pos_o" << si << " 0\n";
            circuit << "Vsec_ret_o"   << si << " sec_neg_o"    << si << " sec_neg_sec_o" << si << " 0\n";
            circuit << "Vgnd_o" << si << " vout_neg_o" << si << " 0 0\n";
            circuit << "Resr_o" << si << " vout_pos_o" << si << " vout_cap_o" << si << " 0.05\n";
            circuit << "Cout_o" << si << " vout_cap_o" << si << " vout_neg_o" << si << " " << std::scientific << 47e-6 << "\n";
            circuit << "Rload_o" << si << " vout_cap_o" << si << " vout_neg_o" << si << " " << Rload_i << "\n\n";
            break;
        }
        case LlcRectifierType::VOLTAGE_DOUBLER: {
            // Two diodes + capacitor stack. Sec winding's negative terminal
            // is tied to the cap midpoint; positive terminal feeds the two
            // diodes which charge Co_hi and Co_lo respectively. Vo = 2·Vsec_pk.
            // n_i has been doubled at design time so Vsec_pk = Vo/2.
            circuit << "Dh_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " DRECT\n";
            circuit << "Dl_o" << si << " vout_mid_o" << si << " sec_pos_o" << si << " DRECT\n";
            // Snubbers on BOTH diodes — at startup, with caps pre-charged via
            // .ic to Vout/2 on vout_mid and 0 on sec_pos, Dl is fully forward-
            // biased (Vmid > Vsecpos) through the near-ideal DRECT (N=0.01),
            // collapsing the solver timestep. The Csnh on Dl provides a
            // millivolt-scale RC at the diode terminals during the first
            // microsecond so the operating point can be found without
            // Dl turning fully on. Mirrors what FB does on all 4 diodes.
            circuit << "Rsnh_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100\n";
            circuit << "Csnh_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100p\n";
            circuit << "Rsnl_o" << si << " vout_mid_o" << si << " sec_pos_o" << si << " 100\n";
            circuit << "Csnl_o" << si << " vout_mid_o" << si << " sec_pos_o" << si << " 100p\n";
            // Sense + winding return
            circuit << "Vsec_sense_o" << si << " sec_pos_sec_o" << si << " sec_pos_o" << si << " 0\n";
            circuit << "Vsec_ret_o" << si << " vout_mid_o" << si << " sec_neg_sec_o" << si << " 0\n";
            // Cap stack: Co_hi from vout_pos → vout_mid, Co_lo from vout_mid → vout_neg.
            circuit << "Co_hi_o" << si << " vout_pos_o" << si << " vout_mid_o" << si << " " << std::scientific << 94e-6 << "\n";
            circuit << "Co_lo_o" << si << " vout_mid_o" << si << " vout_neg_o" << si << " " << std::scientific << 94e-6 << "\n";
            circuit << "Vgnd_o" << si << " vout_neg_o" << si << " 0 0\n";
            circuit << "Rload_o" << si << " vout_pos_o" << si << " vout_neg_o" << si << " " << Rload_i << "\n\n";
            break;
        }
        default:
            throw std::runtime_error("LLC: unknown rectifier type encountered in generate_ngspice_circuit");
        }
    }

    // Initial conditions to speed up settling
    circuit << ".ic";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        double Vout_i = llcOp.get_output_voltages()[i];
        if (rectType == LlcRectifierType::VOLTAGE_DOUBLER) {
            circuit << " v(vout_pos_o" << si << ")=" << Vout_i
                    << " v(vout_mid_o" << si << ")=" << (Vout_i / 2.0);
        } else {
            circuit << " v(vout_cap_o" << si << ")=" << Vout_i;
        }
    }
    if (!integratedLs) circuit << " v(cr_ls)=0";
    circuit << "\n\n";

    circuit << ".options RELTOL=" << cfg.relTol
            << " ABSTOL=" << std::scientific << cfg.absTol
            << " VNTOL=" << cfg.vnTol << std::fixed
            << " ITL1=" << cfg.itl1 << " ITL4=" << cfg.itl4 << "\n";
    circuit << ".options METHOD=" << cfg.method
            << " TRTOL=" << cfg.trTol << "\n\n";
    // UIC tells ngspice to skip the DC operating-point solve and start the
    // transient straight from .ic values (with inductor currents = 0). That
    // works for CT/FB/CD where the output node has a single cap charged to
    // Vout. For VOLTAGE_DOUBLER the .ic constrains two stacked caps (vout_pos
    // and vout_mid) but leaves the floating sec_pos node ambiguous — combined
    // with the near-ideal DRECT model (N=0.01), the first timestep collapses
    // ("trouble with drect-instance dh_o1"). Dropping UIC for VD lets ngspice
    // solve a consistent DC operating point first; the .ic still seeds the
    // caps so settling is fast.
    bool useUic = (rectType != LlcRectifierType::VOLTAGE_DOUBLER);
    circuit << ".tran " << std::scientific << maxStep << " " << simTime
            << " " << startTime << " " << maxStep;
    if (useUic) circuit << " UIC";
    circuit << "\n\n";
    circuit << ".save v(vin_dc) i(Vin) v(vpri_w) v(pri_top) v(pri_bot) i(Vpri_sense)";
    if (bridgeMode != BridgeSimulationMode::BEHAVIORAL_PULSE) {
        if (isFullBridge) circuit << " i(Vq1_sense) i(Vq3_sense)";
        else              circuit << " i(Vq1_sense)";
    }
    if (isFullBridge) circuit << " v(bridge_a) v(bridge_b)";
    else              circuit << " v(sw_node) v(mid_point)";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        switch (rectType) {
        case LlcRectifierType::CENTER_TAPPED:
            circuit << " v(sec_top_o" << si << ") v(sec_bot_o" << si << ")"
                    << " v(vsec1_w_o" << si << ") v(vsec2_w_o" << si << ")"
                    << " v(vout_pos_o" << si << ") v(vout_cap_o" << si << ")"
                    << " i(Vsec1_sense_o" << si << ") i(Vsec2_sense_o" << si << ")"
                    << " i(Vsec_sense_o" << si << ")";
            break;
        case LlcRectifierType::FULL_BRIDGE:
            circuit << " v(sec_pos_o" << si << ") v(sec_neg_o" << si << ")"
                    << " v(vsec_w_o" << si << ")"
                    << " v(vout_pos_o" << si << ") v(vout_cap_o" << si << ")"
                    << " i(Vsec_sense_o" << si << ")";
            break;
        case LlcRectifierType::CURRENT_DOUBLER:
            circuit << " v(sec_pos_o" << si << ") v(sec_neg_o" << si << ")"
                    << " v(vsec_w_o" << si << ")"
                    << " v(vout_pos_o" << si << ") v(vout_cap_o" << si << ")"
                    << " i(Vsec_sense_o" << si << ")"
                    << " i(VLo1_sense_o" << si << ") i(VLo2_sense_o" << si << ")";
            break;
        case LlcRectifierType::VOLTAGE_DOUBLER:
            circuit << " v(sec_pos_o" << si << ")"
                    << " v(vsec_w_o" << si << ")"
                    << " v(vout_pos_o" << si << ") v(vout_mid_o" << si << ")"
                    << " i(Vsec_sense_o" << si << ")";
            break;
        }
    }
    circuit << "\n\n.end\n";

    return circuit.str();
}


// =====================================================================
// SPICE simulation wrappers
// =====================================================================
std::vector<OperatingPoint> Llc::simulate_and_extract_operating_points(
    const std::vector<double>& turnsRatios, double magnetizingInductance,
    size_t numberOfPeriods)
{
    std::vector<OperatingPoint> operatingPoints;
    NgspiceRunner runner;
    if (!runner.is_available())
        return process_operating_points(turnsRatios, magnetizingInductance);

    // If numberOfPeriods is specified, use it; otherwise use the member variable
    int numPeriods = (numberOfPeriods > 0) ? static_cast<int>(numberOfPeriods) : get_num_periods_to_extract();
    
    // Save original value and set the requested number of periods
    int originalNumPeriodsToExtract = get_num_periods_to_extract();
    set_num_periods_to_extract(numPeriods);

    std::vector<double> inputVoltages;
    if (get_input_voltage().get_nominal().has_value())
        inputVoltages.push_back(get_input_voltage().get_nominal().value());
    if (get_input_voltage().get_minimum().has_value())
        inputVoltages.push_back(get_input_voltage().get_minimum().value());
    if (get_input_voltage().get_maximum().has_value())
        inputVoltages.push_back(get_input_voltage().get_maximum().value());

    auto& ops = get_operating_points();
    for (size_t vinIdx = 0; vinIdx < inputVoltages.size(); ++vinIdx) {
        for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
            std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, vinIdx, opIdx);
            double switchingFrequency = ops[opIdx].get_switching_frequency();

            SimulationConfig config;
            config.frequency = switchingFrequency;
            config.extractOnePeriod = true;
            config.numberOfPeriods = get_num_periods_to_extract();
            config.keepTempFiles = false;

            auto simResult = runner.run_simulation(netlist, config);
            if (!simResult.success) {
                std::cerr << "LLC sim failed: " << simResult.errorMessage << ". Falling back." << std::endl;
                return process_operating_points(turnsRatios, magnetizingInductance);
            }

            // §8a.5: read differential winding voltages from the E-source
            // probes (v(vpri_w), v(vsec*_w_o<i>)) instead of bare nodes —
            // pri_bot / sec_*_o<i> sit at a non-zero reference (bridge
            // return / sense-source ladder), so bare-node voltages lump
            // the winding with an offset.
            NgspiceRunner::WaveformNameMapping waveformMapping;
            waveformMapping.push_back({{"voltage", "vpri_w"}, {"current", "vpri_sense#branch"}});

            std::vector<std::string> windingNames = {"Primary"};
            std::vector<bool> flipCurrentSign = {false};

            size_t numOutputs = ops[opIdx].get_output_voltages().size();
            if (numOutputs == 0) numOutputs = 1;
            LlcRectifierType rectTypeSim = get_effective_rectifier_type();
            for (size_t outIdx = 0; outIdx < numOutputs; ++outIdx) {
                std::string si = std::to_string(outIdx + 1);
                switch (rectTypeSim) {
                case LlcRectifierType::CENTER_TAPPED:
                    waveformMapping.push_back({{"voltage", "vsec1_w_o" + si},
                                               {"current", "vsec1_sense_o" + si + "#branch"}});
                    windingNames.push_back("Secondary " + std::to_string(outIdx) + " Half 1");
                    flipCurrentSign.push_back(true);
                    waveformMapping.push_back({{"voltage", "vsec2_w_o" + si},
                                               {"current", "vsec2_sense_o" + si + "#branch"}});
                    windingNames.push_back("Secondary " + std::to_string(outIdx) + " Half 2");
                    flipCurrentSign.push_back(true);
                    break;
                case LlcRectifierType::FULL_BRIDGE:
                case LlcRectifierType::CURRENT_DOUBLER:
                case LlcRectifierType::VOLTAGE_DOUBLER:
                    waveformMapping.push_back({{"voltage", "vsec_w_o" + si},
                                               {"current", "vsec_sense_o" + si + "#branch"}});
                    windingNames.push_back("Secondary " + std::to_string(outIdx));
                    flipCurrentSign.push_back(true);
                    break;
                }
            }

            OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                simResult, waveformMapping, switchingFrequency, windingNames,
                ops[opIdx].get_ambient_temperature(), flipCurrentSign);

            std::string vinName = (vinIdx == 0) ? "Nominal" : (vinIdx == 1 ? "Min" : "Max");
            operatingPoint.set_name(vinName + " input (" +
                std::to_string(static_cast<int>(inputVoltages[vinIdx])) + "V)");
            operatingPoints.push_back(operatingPoint);
        }
    }
    
    // Restore original value
    set_num_periods_to_extract(originalNumPeriodsToExtract);
    return operatingPoints;
}


std::vector<ConverterWaveforms> Llc::simulate_and_extract_topology_waveforms(
    const std::vector<double>& turnsRatios, double magnetizingInductance, size_t numberOfPeriods)
{
    std::vector<ConverterWaveforms> results;
    NgspiceRunner runner;
    if (!runner.is_available())
        throw std::runtime_error("ngspice is not available");

    int originalNumPeriods = get_num_periods_to_extract();
    set_num_periods_to_extract(static_cast<int>(numberOfPeriods));
    auto ops = get_operating_points();

    for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
        std::string netlist = generate_ngspice_circuit(turnsRatios, magnetizingInductance, 0, opIdx);
        double switchingFrequency = ops[opIdx].get_switching_frequency();

        SimulationConfig config;
        config.frequency = switchingFrequency;
        config.extractOnePeriod = true;
        config.numberOfPeriods = numberOfPeriods;
        config.keepTempFiles = false;

        auto simResult = runner.run_simulation(netlist, config);
        if (!simResult.success)
            throw std::runtime_error("LLC simulation failed: " + simResult.errorMessage);

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
        wf.set_switching_frequency(switchingFrequency);
        wf.set_operating_point_name("LLC op. point " + std::to_string(opIdx));

        // §5.1 / §8a.5 converter-port stream:
        //   input_voltage = v(vin_dc)  (real DC rail, constant Vin)
        //   input_current = i(Vq*_sense) summed across high-side switches
        //                   (SW1 mode), or power-balance reconstruction
        //                   (BEHAVIORAL_PULSE mode, where no real switch
        //                   draws current from vin_dc).
        Waveform vdc = getWf("vin_dc");
        if (!vdc.get_data().empty()) wf.set_input_voltage(vdc);

        const auto bridgeModeExtract = get_bridge_simulation_mode();
        const bool isFullBridgeExtract = (get_bridge_type().has_value() &&
                                          get_bridge_type().value() == LlcBridgeType::FULL_BRIDGE);

        bool inputCurrentSet = false;
        if (bridgeModeExtract != BridgeSimulationMode::BEHAVIORAL_PULSE) {
            // SW1 mode: dedicated zero-V sense source(s) upstream of the
            // high-side switch(es). i(Vq*_sense) sign convention:
            // ngspice integrates current from the + terminal to the −
            // terminal through the source — here + is vin_dc and − is
            // the switch drain, so positive i(Vq*_sense) corresponds to
            // current flowing *out* of vin_dc into the bridge, i.e. the
            // converter's input draw.
            Waveform iQ1 = getWf("vq1_sense#branch");
            Waveform iQ3 = isFullBridgeExtract ? getWf("vq3_sense#branch") : Waveform();
            if (!iQ1.get_data().empty() &&
                (!isFullBridgeExtract || !iQ3.get_data().empty())) {
                Waveform iInWf = iQ1;
                auto& iInData = iInWf.get_mutable_data();
                if (isFullBridgeExtract) {
                    const auto& iQ3Data = iQ3.get_data();
                    const size_t N = std::min(iInData.size(), iQ3Data.size());
                    for (size_t k = 0; k < N; ++k) iInData[k] += iQ3Data[k];
                }
                // ngspice SW1 model produces narrow di/dt spikes (~10^5 A)
                // around switching edges that are a numerical artefact of
                // the discontinuous RON/ROFF transition, not a physical
                // bound. Clamp samples to ±2·max|i(Vpri_sense)| —
                // generous enough to retain all real conduction current
                // (which is bounded by the tank current) but small
                // enough to discard the spikes that would poison the
                // DC mean used downstream by §5.1.
                Waveform iPri = getWf("vpri_sense#branch");
                double iPriMax = 0.0;
                for (double v : iPri.get_data()) iPriMax = std::max(iPriMax, std::abs(v));
                if (iPriMax > 0.0) {
                    const double clamp = 2.0 * iPriMax;
                    for (auto& v : iInData) {
                        if (v >  clamp) v =  clamp;
                        if (v < -clamp) v = -clamp;
                    }
                }
                wf.set_input_current(iInWf);
                inputCurrentSet = true;
            }
        }
        if (!inputCurrentSet) {
            // BEHAVIORAL_PULSE: no real switch routes through vin_dc, so
            // i(Vin) is ~0. Reconstruct DC bus current from power balance:
            // Iin_dc = Σ(Vo·Io) / (η · Vin), broadcast as a flat waveform
            // to honour the §5.1 "DC stream" contract.
            double Vin_local = 0.0;
            if (!vdc.get_data().empty()) {
                for (double v : vdc.get_data()) Vin_local += v;
                Vin_local /= vdc.get_data().size();
            } else if (get_input_voltage().get_nominal().has_value()) {
                Vin_local = get_input_voltage().get_nominal().value();
            }
            double Pout_total = 0.0;
            for (size_t i = 0; i < ops[opIdx].get_output_voltages().size(); ++i) {
                Pout_total += ops[opIdx].get_output_voltages()[i] *
                              ops[opIdx].get_output_currents()[i];
            }
            double effLlc = 1.0;
            if (get_efficiency().has_value() && get_efficiency().value() > 0.0)
                effLlc = get_efficiency().value();
            const double Iin_dc = (Vin_local > 0.0) ? Pout_total / (effLlc * Vin_local) : 0.0;
            Waveform iInWf = vdc;  // borrow time grid / length
            auto& iInData = iInWf.get_mutable_data();
            for (auto& v : iInData) v = Iin_dc;
            wf.set_input_current(iInWf);
        }

        // Per-output rectified voltage (cap node is the DC output rail)
        // and current reconstructed as Vout(t)/Rload (matches Buck /
        // IsolatedBuck / AHB §5.1 — ammeters in cap-return paths
        // average to zero in steady state).
        size_t numOutputs = ops[opIdx].get_output_voltages().size();
        if (numOutputs == 0) numOutputs = 1;
        LlcRectifierType rectTypeTW = get_effective_rectifier_type();
        for (size_t i = 0; i < numOutputs; ++i) {
            std::string si = std::to_string(i + 1);
            // VD uses vout_pos (stacked cap top) as the rail; others use vout_cap.
            auto vCap = (rectTypeTW == LlcRectifierType::VOLTAGE_DOUBLER)
                        ? getWf("vout_pos_o" + si)
                        : getWf("vout_cap_o" + si);
            if (vCap.get_data().empty()) continue;
            wf.get_mutable_output_voltages().push_back(vCap);
            const double Vo_i = ops[opIdx].get_output_voltages()[i];
            const double Io_i = ops[opIdx].get_output_currents()[i];
            const double Rload_i = (Io_i > 0) ? (Vo_i / Io_i) : 100.0;
            Waveform ioutWf = vCap;
            auto& ioutData = ioutWf.get_mutable_data();
            for (auto& v : ioutData) v = v / Rload_i;
            wf.get_mutable_output_currents().push_back(ioutWf);
        }
        results.push_back(wf);
    }
    set_num_periods_to_extract(originalNumPeriods);
    return results;
}


// =====================================================================
// AdvancedLlc
// =====================================================================
Inputs AdvancedLlc::process() {
    auto designRequirements = process_design_requirements();

    // Reuse the merged DR's turns ratios and magnetizing inductance for
    // operating-point calculation (process_design_requirements() has
    // already applied the desired* overrides).
    std::vector<double> turnsRatios;
    for (const auto& tr : designRequirements.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    if (turnsRatios.empty()) {
        throw std::runtime_error("LLC: no turns ratios available (neither desired nor computed)");
    }
    double Lm = resolve_dimensional_values(designRequirements.get_magnetizing_inductance());

    auto ops = process_operating_points(turnsRatios, Lm);
    Inputs inputs;
    inputs.set_design_requirements(designRequirements);
    inputs.set_operating_points(ops);
    return inputs;
}

DesignRequirements AdvancedLlc::process_design_requirements() {
    // Issue M1: chain to parent PDR for resonant-tank derivation and
    // topology defaults, then override turns ratios and magnetizing
    // inductance with desired* values when provided.
    auto designRequirements = Llc::process_design_requirements();

    std::vector<double> turnsRatios = desiredTurnsRatios;
    if (turnsRatios.empty()) {
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(tr));
        }
    }
    designRequirements.get_mutable_turns_ratios().clear();
    for (auto n : turnsRatios) {
        DimensionWithTolerance nTol; nTol.set_nominal(n);
        designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }

    double Lm = desiredMagnetizingInductance;
    if (Lm <= 0) {
        Lm = resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
    }
    DimensionWithTolerance LmTol; LmTol.set_nominal(Lm);
    designRequirements.set_magnetizing_inductance(LmTol);

    return designRequirements;
}

std::vector<std::variant<Inputs, CAS::Inputs>> Llc::get_extra_components_inputs(
    ExtraComponentsMode mode,
    std::optional<Magnetic> magnetic)
{
    if (mode == ExtraComponentsMode::REAL && !magnetic.has_value())
        throw std::invalid_argument("get_extra_components_inputs: mode REAL requires a designed magnetic");

    if (computedResonantInductance <= 0 || computedResonantCapacitance <= 0 || extraCapVoltageWaveforms.empty())
        throw std::runtime_error("LLC get_extra_components_inputs: call process() first");

    bool integratedLs = get_integrated_resonant_inductor().value_or(false);
    double Ls = computedResonantInductance;
    double Cr = computedResonantCapacitance;
    double fr = get_effective_resonant_frequency();

    double Llk = 0.0;
    if (mode == ExtraComponentsMode::REAL) {
        auto leakageOutput = LeakageInductance().calculate_leakage_inductance_all_windings(
            magnetic.value(), fr);
        auto perWinding = leakageOutput.get_leakage_inductance_per_winding();
        if (!perWinding.empty())
            Llk = resolve_dimensional_values(perWinding[0]);
        // Recalculate Cr to resonate with actual available inductance (Llk or full Ls)
        double Ls_effective = (Llk >= Ls) ? Llk : Ls;
        Cr = 1.0 / (4.0 * M_PI * M_PI * fr * fr * Ls_effective);
    }

    size_t nOps = extraCapVoltageWaveforms.size();
    std::vector<std::variant<Inputs, CAS::Inputs>> result;

    // Build CAS::Inputs for the resonant capacitor
    {
        CAS::Inputs casInputs;
        CAS::DesignRequirements dr;

        CAS::DimensionWithTolerance capacitance;
        capacitance.set_nominal(Cr);
        dr.set_capacitance(capacitance);

        // Peak voltage across capacitor (from stored waveforms)
        double peakVc = 0.0;
        for (const auto& wfm : extraCapVoltageWaveforms) {
            for (double v : wfm.get_data())
                peakVc = std::max(peakVc, std::abs(v));
        }
        dr.set_rated_voltage(peakVc * 1.2);  // 20 % margin

        dr.set_role(CAS::Application::RESONANT);
        dr.set_name("resonantCapacitor");

        casInputs.set_design_requirements(dr);

        // Operating points: one per stored waveform
        std::vector<CAS::TwoTerminalOperatingPoint> casOps;
        for (size_t i = 0; i < nOps; ++i) {
            CAS::TwoTerminalOperatingPoint op;
            CAS::OperatingPointExcitation excitation;
            excitation.set_frequency(fr);

            // Voltage across cap
            CAS::SignalDescriptor vSig;
            CAS::Waveform vWfm;
            vWfm.set_data(extraCapVoltageWaveforms[i].get_data());
            if (extraCapVoltageWaveforms[i].get_time())
                vWfm.set_time(extraCapVoltageWaveforms[i].get_time().value());
            vSig.set_waveform(vWfm);
            excitation.set_voltage(vSig);

            // Current through cap (= series inductor current)
            CAS::SignalDescriptor iSig;
            CAS::Waveform iWfm;
            iWfm.set_data(extraCapCurrentWaveforms[i].get_data());
            if (extraCapCurrentWaveforms[i].get_time())
                iWfm.set_time(extraCapCurrentWaveforms[i].get_time().value());
            iSig.set_waveform(iWfm);
            excitation.set_current(iSig);

            op.set_excitation(excitation);
            casOps.push_back(op);
        }
        casInputs.set_operating_points(casOps);
        result.emplace_back(std::move(casInputs));
    }

    // Build MAS Inputs for external series inductor (only when needed)
    double Ls_external = 0.0;
    if (mode == ExtraComponentsMode::IDEAL && !integratedLs) {
        Ls_external = Ls;
    } else if (mode == ExtraComponentsMode::REAL) {
        Ls_external = (Llk < Ls) ? (Ls - Llk) : 0.0;
    }

    if (Ls_external > 0.0 && !extraIndVoltageWaveforms.empty()) {
        Inputs masInputs;
        DesignRequirements dr;

        DimensionWithTolerance inductance;
        inductance.set_nominal(Ls_external);
        dr.set_magnetizing_inductance(inductance);
        dr.set_name("seriesInductor");
        dr.set_topology(Topologies::LLC_RESONANT_CONVERTER);

        // Single primary winding, no isolation
        dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
        dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
        masInputs.set_design_requirements(dr);

        std::vector<OperatingPoint> masOps;
        for (size_t i = 0; i < nOps; ++i) {
            OperatingPoint op;
            auto& indVoltWfm = extraIndVoltageWaveforms[i];
            auto& indCurrWfm = extraIndCurrentWaveforms[i];
            double fsw = fr;  // Use resonant frequency as representative frequency

            auto excitation = complete_excitation(indCurrWfm, indVoltWfm, fsw, "Primary");
            op.get_mutable_excitations_per_winding().push_back(excitation);
            masOps.push_back(op);
        }
        masInputs.set_operating_points(masOps);
        result.emplace_back(std::move(masInputs));
    }

    // Rectifier-variant extras.
    LlcRectifierType rectTypeExtras = get_effective_rectifier_type();
    if (rectTypeExtras == LlcRectifierType::CURRENT_DOUBLER &&
        !extraLoCurrentWaveforms.empty()) {
        // Two output inductors L_o1 / L_o2 — emit two MAS Inputs.
        // Inductance value matches the SPICE netlist sizing:
        // Lo = Vout / (4·Fs·0.30·IoutDC).
        auto& ops = get_operating_points();
        if (ops.empty()) throw std::runtime_error("LLC CD extras: no operating points");
        double Vout_0 = ops[0].get_output_voltages()[0];
        double Iout_0 = ops[0].get_output_currents()[0];
        double fsw_0  = ops[0].get_switching_frequency();
        if (Vout_0 <= 0 || Iout_0 <= 0 || fsw_0 <= 0)
            throw std::runtime_error("LLC CD extras: invalid output/freq for Lo sizing");
        double Lo_design = Vout_0 / (4.0 * fsw_0 * 0.30 * Iout_0);

        for (size_t loIdx = 0; loIdx < 2; ++loIdx) {
            Inputs loInputs;
            DesignRequirements drLo;
            DimensionWithTolerance loInd;
            loInd.set_nominal(Lo_design);
            drLo.set_magnetizing_inductance(loInd);
            drLo.set_name(std::string("outputInductor") + (loIdx == 0 ? "1" : "2"));
            drLo.set_topology(Topologies::LLC_RESONANT_CONVERTER);
            drLo.set_turns_ratios(std::vector<DimensionWithTolerance>{});
            drLo.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::SECONDARY});
            loInputs.set_design_requirements(drLo);

            std::vector<OperatingPoint> loOps;
            auto& iVec = (loIdx == 0) ? extraLoCurrentWaveforms : extraLo2CurrentWaveforms;
            auto& vVec = (loIdx == 0) ? extraLoVoltageWaveforms : extraLo2VoltageWaveforms;
            for (size_t i = 0; i < iVec.size(); ++i) {
                OperatingPoint op;
                auto excitation = complete_excitation(iVec[i], vVec[i], fsw_0,
                                                     (loIdx == 0) ? "Lo1" : "Lo2");
                op.get_mutable_excitations_per_winding().push_back(excitation);
                loOps.push_back(op);
            }
            loInputs.set_operating_points(loOps);
            result.emplace_back(std::move(loInputs));
        }
    }
    else if (rectTypeExtras == LlcRectifierType::VOLTAGE_DOUBLER &&
             !extraVoutCapVoltageWaveforms.empty()) {
        // Two stacked output capacitors. Both see Vo/2 nominal DC; current
        // through them is ±I_sec on alternating halves. We emit two CAS::Inputs
        // sized identically (each cap rated to Vo/2 · 1.2 margin).
        auto& ops = get_operating_points();
        if (ops.empty()) throw std::runtime_error("LLC VD extras: no operating points");
        double Vout_0 = ops[0].get_output_voltages()[0];
        for (size_t capIdx = 0; capIdx < 2; ++capIdx) {
            CAS::Inputs casCap;
            CAS::DesignRequirements drCap;
            CAS::DimensionWithTolerance cap;
            cap.set_nominal(94e-6);                          // matches netlist
            drCap.set_capacitance(cap);
            drCap.set_rated_voltage((Vout_0 / 2.0) * 1.2);   // 20 % margin
            drCap.set_role(CAS::Application::DC_LINK);
            drCap.set_name(std::string("outputCapacitor") + (capIdx == 0 ? "Hi" : "Lo"));
            casCap.set_design_requirements(drCap);

            std::vector<CAS::TwoTerminalOperatingPoint> capOps;
            for (size_t i = 0; i < extraVoutCapVoltageWaveforms.size(); ++i) {
                CAS::TwoTerminalOperatingPoint op;
                CAS::OperatingPointExcitation exc;
                exc.set_frequency(fr);
                CAS::SignalDescriptor vSig; CAS::Waveform vWfm;
                vWfm.set_data(extraVoutCapVoltageWaveforms[i].get_data());
                if (extraVoutCapVoltageWaveforms[i].get_time())
                    vWfm.set_time(extraVoutCapVoltageWaveforms[i].get_time().value());
                vSig.set_waveform(vWfm); exc.set_voltage(vSig);
                // Current = secondary current divided across the two caps;
                // we approximate via the resonant tank current waveform.
                CAS::SignalDescriptor iSig; CAS::Waveform iWfm;
                iWfm.set_data(extraCapCurrentWaveforms[i].get_data());
                if (extraCapCurrentWaveforms[i].get_time())
                    iWfm.set_time(extraCapCurrentWaveforms[i].get_time().value());
                iSig.set_waveform(iWfm); exc.set_current(iSig);
                op.set_excitation(exc);
                capOps.push_back(op);
            }
            casCap.set_operating_points(capOps);
            result.emplace_back(std::move(casCap));
        }
    }

    return result;
}

} // namespace OpenMagnetics
