// ─────────────────────────────────────────────────────────────────────────────
// PfcControllerDesign.h
//
// Design document + closed-form tuning helper for the OpenMagnetics native
// boost-PFC average-current-mode controller (no vendor SPICE IP). This header
// is purely a design reference — the netlist builder lives in
// PowerFactorCorrection.{h,cpp} (`generate_ngspice_switching_circuit`) and
// the runner in `simulate_with_ngspice_switching`.
//
// ─── Reference design (published, non-proprietary) ──────────────────────────
//   * Texas Instruments SNVA408B  — "Power Factor Correction Using the
//                                    UCC28019" (control law in math form,
//                                    no encrypted IP).
//   * Erickson & Maksimović §18  — averaged-switch boost PFC small-signal
//                                    model; multiplier-based current shaping.
//   * Christophe Basso, "Switch-Mode Power Supplies: SPICE Simulations and
//     Practical Designs", 2nd ed., Ch. 7 — UC3854-class average-current
//     compensator design walkthrough (used for cross-checking only; no
//     subcircuits copied).
//
// Rebuilt from first principles in plain ngspice primitives so we own the IP.
//
// ─── Topology view (BOOST CCM) ──────────────────────────────────────────────
//
//       + vin_rect ─[Vl_sense]─[ L1 ]─┬───[ D1 ]──┬─── vbus ──┬──┐
//                                     │           │            │  │
//                                  ┌──┴──┐     ┌──┴──┐         │  │
//                                  │ S1  │     │ Cbus│      Rload│
//                                  └──┬──┘     └──┬──┘         │  │
//                                     │           │            │  │
//        gnd ───────────────────────  └───────────┴────────────┴──┘
//                                     ▲
//                                  gate│ (PWM, 0 / VPWM)
//
// Controller (six analog blocks, all built from B/E/G sources + the
// `OPAMP_IDEAL` and `COMPARATOR_IDEAL` subckts in PfcControllerSubcircuits.h):
//
//   ┌─────────────────────────────────────────────────────────────────────┐
//   │  vbus ─[K_div]─►(−)Gv                                              │
//   │                  │                                                  │
//   │             Vref►(+)── vea ────►──────────┐                         │
//   │                                            │  multiplier            │
//   │  vin_rect ──────────────────────►──────────┤                        │
//   │              │                             │ i_ref = K_mul · vea ·  │
//   │              ▼                             │         vin_rect /     │
//   │           ┌──FF──┐    ┌──FF──┐  ┌─sq─┐     │         vrms_ff        │
//   │           │ LP1  │──► │ LP2  │─►│ x² │── vrms_ff ───►───────────────┘
//   │           └──────┘    └──────┘  └────┘                              │
//   │                                                                     │
//   │  i_L ──[Rs]──►(−)Gi                                                │
//   │                │                                                    │
//   │           i_ref►(+)── vc_i ──►(+)comp────► gate                    │
//   │                                  │                                  │
//   │                       saw ─PULSE►(−)                                │
//   └─────────────────────────────────────────────────────────────────────┘
//
// ─── Block 1: Voltage error amp (Gv) ───────────────────────────────────────
//
// Sense: K_div = Vref / Vbus_nom, where Vref = 2.5 V (industry-standard
// internal reference in UCC28019 / NCP1654 / L4981 / UC3854).
// → divider output = Vref when vbus = Vbus_nom (zero error).
//
// Compensator: type-II (integrator + zero + pole) for zero DC error and a
// single phase boost.
//
//   Gv(s) = (1 + s·Rz·Cz) / (s·Rz·(Cz+Cp)·(1 + s·Rz·Cz·Cp/(Cz+Cp)))
//
// Crossover frequency fc_v: must be WELL below the second harmonic of the
// line frequency (2·fline = 100 Hz at 50 Hz mains) to avoid the loop reacting
// to the natural 100 Hz Cbus ripple — that reaction would distort the input
// current shape (3rd-harmonic injection). Industry rule: fc_v ≤ fline / 5.
//
//   fc_v = max(5 Hz, fline / 5)        // 10 Hz @ 50 Hz, 12 Hz @ 60 Hz
//
// Phase boost: 60° → zero at fc_v / 3, pole at 3·fc_v.
//
//   Rz = 100 kΩ                         // arbitrary impedance scale
//   Cz = 1 / (2π · (fc_v/3) · Rz)
//   Cp = 1 / (2π · (3·fc_v) · Rz)
//
// Steady-state output `vea_nom`: opamp saturates between (vssNeg, vssPos).
// We design G_mul so vea = vea_nom = (vssPos + vssNeg)/2 = 2.5 V at rated
// power and nominal Vrms (mid-scale). Off-nominal operating points pull vea
// up (lower line voltage / higher load) or down (higher line / lighter load).
//
// ─── Block 2: RMS² feed-forward (vrms_ff) ──────────────────────────────────
//
// Two cascaded RC low-pass filters on vin_rect, each at fc_ff << 2·fline,
// then a B-source squarer. Two stages get us −40 dB/dec roll-off at the
// 100 Hz ripple frequency, sufficient to keep vrms_ff effectively constant
// over a line cycle (< 1 % residual ripple).
//
//   fc_ff = fline / 5                   // 10 Hz @ 50 Hz
//   tau   = 1 / (2π · fc_ff)
//
// Steady-state output (after both filters):
//   v_ff_filt ≈ mean(|vin_rect|) = (2/π) · Vpk = (2√2/π) · Vrms ≈ 0.9 · Vrms
//
// After squaring:
//   vrms_ff   ≈ ((2√2/π))² · Vrms²    = (8/π²) · Vrms²    ≈ 0.811 · Vrms²
//
// We capture this constant as `K_ff² = (2√2/π)² ≈ 0.811`.
//
// ─── Block 3: Multiplier ───────────────────────────────────────────────────
//
//   i_ref(t) = G_mul · vea · vin_rect(t) / vrms_ff
//
// Goal: at rated operating point, i_ref(t) tracks the desired line-current
// envelope, which by power balance (η = 1 in the ideal model) is
//
//   i_L_avg(t) = (Pout / Vrms²) · vin_rect(t)
//              = (√2 · Pout / Vrms) · |sin(ωt)|         [envelope peak √2·Pout/Vrms]
//
// At rated condition (vea = vea_nom, vin_rect = vpk·|sin|, vrms_ff = K_ff²·Vrms²):
//
//   i_ref_env_peak = G_mul · vea_nom · Vpk / (K_ff² · Vrms²)
//                  = G_mul · vea_nom · √2 / (K_ff² · Vrms)
//
// Set equal to the desired envelope peak √2 · Pout / Vrms:
//
//   G_mul = K_ff² · Pout / vea_nom            [units: A/(V·V/V) = A]
//
// (With K_ff² ≈ 0.811, vea_nom = 2.5 V, this gives G_mul ≈ 0.324·Pout in
// amps — see the worked numbers in derive_pfc_controller_tuning().)
//
// Numerical guard at line zero crossings: vrms_ff is approximately constant
// across the cycle so it never approaches zero, but we still clamp the
// divisor inside the B-source to MAX(vrms_ff, 0.01 · K_ff² · Vrms_nom²) to
// be defensive against startup transients before the FF filter has settled.
//
// ─── Block 4: Current error amp (Gi) ──────────────────────────────────────
//
// Plant (boost-CCM small-signal duty-to-inductor-current at envelope peak):
//
//   Gid(s) = Vbus / (s · L)
//
// Sense gain Rs is set to 1 (treat i_L directly as a voltage signal in the
// controller — physically equivalent to a 1 Ω shunt; the magnetic-component
// model doesn't care about the absolute voltage scale of the controller
// internals so long as the loop closes correctly at fc_i).
//
// PWM gain: 1 / V_pk_saw (V_pk_saw = peak of sawtooth = 1 V by convention).
//
// Loop gain (proportional only):
//   Tu(s) = Kp · Vbus / (s · L · V_pk_saw)
//
// Target crossover fc_i = fsw / 10 (industry rule: well below fsw to avoid
// switching-ripple coupling, well above 2·fline so the current loop tracks
// the line-frequency envelope cleanly).
//
//   Kp = 2π · fc_i · L · V_pk_saw / Vbus
//
// Add a type-II boost (zero at fc_i/3, pole at 3·fc_i) for ~60° phase margin
// — same topology as Gv:
//
//   Rz_i = 100 kΩ
//   Cz_i = 1 / (2π · (fc_i/3) · Rz_i)
//   Cp_i = 1 / (2π · (3·fc_i) · Rz_i)
//
// The integrator capacitor sets DC gain; in the type-II topology, total DC
// gain is set by the input resistor R_in_i = Rz_i / Kp (so the proportional
// gain at mid-band equals Kp).
//
// ─── Block 5: Oscillator ───────────────────────────────────────────────────
//
//   V_saw saw 0  PULSE(0  V_pk_saw  0  Tsw-1n  1n  1n  Tsw)
//
// Linear ramp 0 → V_pk_saw over each Tsw period, instantaneous reset.
// Identical to the existing Boost / Buck PWM oscillator.
//
// ─── Block 6: PWM comparator ──────────────────────────────────────────────
//
//   B_gate gate 0 V=(V(vc_i) > V(saw)) ? VPWM : 0
//
// Hard threshold (no hysteresis); the analog comparator subckt
// `COMPARATOR_IDEAL` is reserved for future variants where comparator
// dynamics matter.
//
// ─── Power stage steady-state initial conditions ──────────────────────────
//
// Cold-start of a switching PFC takes many line cycles for the voltage loop
// to charge Cbus from 0 to Vbus_nom. We bypass this entirely by warm-starting
// every state variable to its analytical steady-state value:
//
//   v(vbus) = Vbus_nom                 .ic
//   v(vea)  = vea_nom = 2.5 V          .ic
//   v(vff_lp1) = K_ff·Vrms             .ic   // mean(|vin|) = (2√2/π)·Vrms
//   v(vff_lp2) = K_ff·Vrms             .ic
//
// A `uic` flag on `.tran` tells ngspice to honour these and skip its OP
// solve (which often fails for nonlinear PFC topologies anyway).
//
// ─── ngspice solver settings ──────────────────────────────────────────────
//
// Hard-switching events at fsw with a line-cycle envelope require a balance
// between fidelity (small enough timestep to resolve switching edges) and
// throughput (3 line cycles × 50–100 kHz = 6 000 – 30 000 switching events).
//
//   .options METHOD=GEAR  TRTOL=7  RELTOL=1e-3  ABSTOL=1e-9  VNTOL=1e-6
//   .options ITL1=500     ITL4=200
//   .tran   stepMax=Tsw/40  tStop=N·Tline   tStart=0  stepMax  uic
//
// GEAR (vs default TRAP) damps switching-edge ringing that otherwise causes
// false convergence failures. TRTOL=7 (vs default 7) keeps the local
// truncation error loose enough to cross switching events without
// step-rejection storms. ITL1/4 raised to absorb startup convergence churn
// even with warm-start ICs.
//
// ─────────────────────────────────────────────────────────────────────────────

#pragma once

#include <string>

namespace OpenMagnetics {

/**
 * @brief Derived component values for the boost-PFC average-current-mode
 *        controller. Computed in closed form from the converter spec by
 *        `derive_pfc_controller_tuning()`; consumed by
 *        `PowerFactorCorrection::generate_ngspice_switching_circuit()`.
 *
 * All resistances in Ω, all capacitances in F, all frequencies in Hz, all
 * voltages in V. Naming convention follows the design note above:
 *   - prefix `gv_*` → voltage error amp (Block 1)
 *   - prefix `ff_*` → RMS² feed-forward (Block 2)
 *   - field  `g_mul` → multiplier gain (Block 3)
 *   - prefix `gi_*` → current error amp (Block 4)
 *   - prefix `pwm_*` → oscillator + comparator (Blocks 5–6)
 *   - prefix `ic_*` → initial conditions (warm start)
 */
struct PfcControllerTuning {
    // ── Block 1: voltage error amp (Gv) ─────────────────────────────────
    double vref       = 2.5;        ///< Internal reference (industry standard)
    double k_div      = 0.0;        ///< Vbus divider gain = vref / Vbus_nom
    double gv_fc_hz   = 0.0;        ///< Voltage-loop crossover [Hz]
    double gv_rz      = 100e3;      ///< Type-II Rz [Ω]
    double gv_cz      = 0.0;        ///< Type-II Cz [F] — sets zero at fc_v/3
    double gv_cp      = 0.0;        ///< Type-II Cp [F] — sets pole at 3·fc_v
    double vea_nom    = 2.5;        ///< Steady-state opamp output at rated cond.
    double vea_min    = 0.1;        ///< Opamp negative rail (just above 0)
    double vea_max    = 4.9;        ///< Opamp positive rail

    // ── Block 2: RMS² feed-forward (vrms_ff) ────────────────────────────
    double ff_fc_hz   = 0.0;        ///< Each LP stage corner [Hz]
    double ff_r       = 100e3;      ///< RC R [Ω]
    double ff_c       = 0.0;        ///< RC C [F] — set so 2π·R·C = 1/fc_ff
    double k_ff_sq    = 0.0;        ///< (2√2/π)² ≈ 0.811 — populated by
                                    ///< derive_pfc_controller_tuning() to keep
                                    ///< the value in exactly one place (the
                                    ///< canonical k2Sqrt2OverPi constant in
                                    ///< PfcControllerDesign.cpp).

    // ── Block 3: multiplier ─────────────────────────────────────────────
    double g_mul      = 0.0;        ///< Multiplier gain [A] — see design note
    double vrms_ff_floor = 0.0;     ///< Numerical clamp on divisor [V²]

    // ── Block 4: current error amp (Gi) ─────────────────────────────────
    double gi_fc_hz   = 0.0;        ///< Current-loop crossover [Hz]
    double gi_kp      = 0.0;        ///< Mid-band proportional gain [V/V]
    double gi_rz      = 100e3;      ///< Type-II Rz [Ω]
    double gi_rin     = 0.0;        ///< Input R [Ω] — sets Kp = Rz/Rin
    double gi_cz      = 0.0;        ///< Type-II Cz [F] — zero at fc_i/3
    double gi_cp      = 0.0;        ///< Type-II Cp [F] — pole at 3·fc_i
    double rs_sense   = 1.0;        ///< Current-sense gain [V/A]

    // ── Block 5–6: oscillator + comparator ──────────────────────────────
    double pwm_v_pk_saw = 1.0;      ///< Sawtooth peak [V]
    double pwm_v_high   = 5.0;      ///< Gate logic high [V]
    double pwm_t_rise   = 1e-9;     ///< Sawtooth reset edge [s]

    // ── Initial conditions (warm start) ─────────────────────────────────
    double ic_vbus    = 0.0;        ///< = Vpk (rectified peak — bridge clamp)
    double ic_vea     = 0.0;        ///< = vea_nom
    double ic_vff     = 0.0;        ///< = (2√2/π) · Vrms_nom — mean(|vin|)
    double ic_vc_i    = 0.0;        ///< = D_target · V_pk_saw
                                    ///<   D_target = 1 − Vpk/Vbus_nom
                                    ///<   (NOT D_avg — D_avg over-pumps Cbus
                                    ///<   to Vpk/(1−D_avg) at start; D_target
                                    ///<   = duty needed to hold Vbus = Vbus_nom
                                    ///<   exactly at the line peak).

    // ── Soft-start (NEW: boost-action establishment) ────────────────────
    // During the first t_ss_release seconds, the PWM control voltage seen
    // by the comparator is clamped from BELOW by a linearly-decaying floor
    //
    //   ss_floor(t) = ic_vc_i · max(0, 1 − t/t_ss_release)
    //
    // This forces the gate to fire at duty ≈ D_target initially (boosting
    // Cbus from Vpk up to Vbus_nom without overshoot), and gracefully
    // hands control back to the current EA as the floor decays. Without
    // this, the inrush charge pulse via L1+D1 each line cycle drives the
    // current EA's integrator capacitor C_fb_zi irreversibly negative
    // (i_sense > i_ref always when no boost action), deadlocking the
    // gate at OFF for the entire simulation.
    //
    // 30 ms = 1.5 line cycles at 50 Hz — long enough to span 3 inrush
    // events (rising-edge of each rectified half-cycle) so steady-state
    // operating point is established before the floor releases.
    double t_ss_release = 30e-3;    ///< Soft-start release time [s]
};

/**
 * @brief Compute the controller tuning in closed form from the PFC spec.
 *
 * All formulas are derived in the design note at the top of this header.
 * Throws `std::invalid_argument` if any input is non-positive or non-finite,
 * or `std::runtime_error` if the resulting design is internally inconsistent
 * (e.g. fc_i would be required to exceed fsw/2). NEVER returns silently with
 * default values — per project no-fallbacks rule.
 *
 * @param vrms_nom  Nominal line RMS [V]
 * @param vbus_nom  Regulated bus voltage [V]
 * @param pout      Rated output power [W]
 * @param fsw       Switching frequency [Hz]
 * @param fline     Line frequency [Hz]
 * @param inductance Boost inductance [H]
 */
PfcControllerTuning derive_pfc_controller_tuning(double vrms_nom,
                                                   double vbus_nom,
                                                   double pout,
                                                   double fsw,
                                                   double fline,
                                                   double inductance);

} // namespace OpenMagnetics
