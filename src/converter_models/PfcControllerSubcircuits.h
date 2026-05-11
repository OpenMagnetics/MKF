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
 * @brief OPAMP_IDEAL — single-pole opamp with output rail clamp.
 *
 * Pins: vp vn vo
 * Params: A0 (DC open-loop gain, V/V), GBW (unity-gain bandwidth, Hz),
 *         VSSPOS (positive rail, V), VSSNEG (negative rail, V).
 *
 * Topology:
 *   diff = A0 · (V(vp) − V(vn))                    (B-source, no filtering)
 *   filt = LP(diff, fp = GBW/A0)                   (RC pole; unity-crossover at GBW)
 *   vo   = clamp(V(filt), VSSNEG, VSSPOS)          (B-source rail limiter)
 *
 * Validated against:
 *   - DC accuracy: open-loop gain |Δvo/Δvd| ≈ A0 for vd in linear range
 *   - Voltage-follower bandwidth: closed-loop −3 dB at fc ≈ GBW
 *   - Rail clamp: vo saturates at VSSPOS / VSSNEG when input drives past
 */
inline constexpr const char* SPICE_SUBCKT_OPAMP_IDEAL = R"NGSPICE(
.subckt OPAMP_IDEAL vp vn vo params: A0=1e5 GBW=1e6 VSSPOS=5 VSSNEG=0
B_diff diff 0 V = A0 * (V(vp) - V(vn))
R_pole diff filt 1k
C_pole filt 0 {A0/(6.283185307*1k*GBW)}
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
