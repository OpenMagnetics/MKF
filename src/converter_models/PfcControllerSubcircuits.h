// ─────────────────────────────────────────────────────────────────────────────
// PfcControllerSubcircuits.h
//
// Reusable ngspice `.subckt` definitions for the analog blocks used by the
// boost-PFC average-current-mode controller (and potentially other future
// converter controllers). Embedded as raw C++ string constants so the
// netlist builder doesn't depend on a runtime filesystem lookup.
//
// Both subcircuits are validated by `tests/TestSpiceSubcircuits.cpp`.
//
// ─── Design rationale (see PfcControllerDesign.h for the full PFC context) ──
//
// We deliberately avoid vendor-supplied opamp and comparator macromodels
// (LTspice / PSpice / SIMetrix encrypted libraries are not portable to
// ngspice, and unencrypted ones often use simulator-specific extensions).
// These two subckts are the minimum analog primitives needed to build an
// average-current-mode PFC: a single-pole opamp (for the voltage and current
// error amplifiers) and a hard-threshold comparator (reserved for variants
// where the inline B-source PWM comparator isn't sufficient).
//
// ─────────────────────────────────────────────────────────────────────────────

#pragma once

#include <string>

namespace OpenMagnetics {

/**
 * @brief OPAMP_IDEAL — single-pole opamp with output rail clamp
 *        AND anti-windup on the internal state.
 *
 * Pins: vp vn vo
 * Params: A0 (DC open-loop gain, V/V), GBW (unity-gain bandwidth, Hz),
 *         VSSPOS (positive rail, V), VSSNEG (negative rail, V),
 *         IC_FILT (initial voltage on the internal `filt` state node, V —
 *                  use this for warm-start; .ic on the external `vo` pin
 *                  is silently ignored because `vo` is a B-source output,
 *                  not a state variable).
 *
 * Topology:
 *   diff = A0·(V(vp)−V(vn))                          (B-source, no clamp)
 *   filt = LP(diff, fp = GBW/A0)                     (RC pole)
 *   vo   = clamp(V(filt), VSSNEG, VSSPOS)            (B-source rail clamp)
 *   anti-windup: a high-gain B-source pulls filt back inside
 *                [VSSNEG, VSSPOS] whenever it tries to wind past a rail.
 *                Linear-region behavior (when filt is between rails) is
 *                untouched — the AW source contributes zero current.
 *
 * Why anti-windup matters: when the opamp is wrapped in a control loop
 * where an external constraint (comparator saturation, load nonlinearity,
 * boost current sat at line zero crossing) prevents the loop from
 * converging, an unbounded integrator state can wind to ±∞ and pin
 * the output at one rail for many time constants after the constraint
 * releases. Empirically (PFC NCP1654 startup), without AW the current EA
 * filt wound to −5 V within 0.2 ms and never recovered — gate stayed
 * permanently at 0, no boost action. With AW, filt is hard-bounded to
 * the linear-control range and recovers in one EA time constant.
 *
 * Validated against (TestSpiceSubcircuits.cpp):
 *   - DC accuracy: open-loop gain |Δvo/Δvd| ≈ A0 for vd in linear range
 *   - Voltage-follower bandwidth: closed-loop −3 dB at fc ≈ GBW
 *   - Rail clamp: vo saturates at VSSPOS / VSSNEG when input drives past
 *   - Warm-start IC: `IC_FILT={value}` instance param sets vo at t=0
 */
inline constexpr const char* SPICE_SUBCKT_OPAMP_IDEAL = R"NGSPICE(
.subckt OPAMP_IDEAL vp vn vo params: A0=1e5 GBW=1e6 VSSPOS=5 VSSNEG=0 IC_FILT=0
B_diff diff 0 V = A0 * (V(vp) - V(vn))
R_pole diff filt 1k
C_pole filt 0 {A0/(6.283185307*1k*GBW)} IC={IC_FILT}
* Anti-windup rail clamps — near-ideal diodes with Vt=25mV·0.05=1.25 mV.
* When filt would exceed VSSPOS, D_aw_p conducts strongly to rail_p,
* draining whatever current the upstream R_pole pushes through. Same on
* the negative rail with D_aw_n. In the linear region (VSSNEG < filt <
* VSSPOS), both diodes are reverse-biased and contribute zero current.
.model DAWIDEAL D (IS=1e-12 RS=1m N=0.05)
V_rail_p_aw rail_p_aw 0 {VSSPOS}
V_rail_n_aw rail_n_aw 0 {VSSNEG}
D_aw_p filt        rail_p_aw DAWIDEAL
D_aw_n rail_n_aw   filt      DAWIDEAL
B_clamp vo 0 V = min(VSSPOS, max(VSSNEG, V(filt)))
.ends OPAMP_IDEAL
)NGSPICE";

/**
 * @brief COMPARATOR_IDEAL — zero-hysteresis B-source comparator.
 *
 * Pins: vp vn vo
 * Params: VHIGH (output high level, V), VLOW (output low level, V).
 *
 * Behaviour:
 *   vo = (V(vp) > V(vn)) ? VHIGH : VLOW
 *
 * Reserved for future controller variants (the in-line PWM comparator in
 * `generate_ngspice_switching_circuit` uses an equivalent inline B-source
 * to avoid the subckt-call overhead in the hottest loop of the netlist).
 */
inline constexpr const char* SPICE_SUBCKT_COMPARATOR_IDEAL = R"NGSPICE(
.subckt COMPARATOR_IDEAL vp vn vo params: VHIGH=5 VLOW=0
B_cmp vo 0 V = (V(vp) > V(vn)) ? VHIGH : VLOW
.ends COMPARATOR_IDEAL
)NGSPICE";

/**
 * @brief Concatenated prelude block — paste both subckts into a netlist
 *        with one `<<` (or one `+=`). Convenience helper.
 */
inline std::string spice_subckt_prelude_pfc_controller() {
    std::string s;
    s += SPICE_SUBCKT_OPAMP_IDEAL;
    s += SPICE_SUBCKT_COMPARATOR_IDEAL;
    return s;
}

} // namespace OpenMagnetics
