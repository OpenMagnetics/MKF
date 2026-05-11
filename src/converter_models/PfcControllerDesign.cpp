// ─────────────────────────────────────────────────────────────────────────────
// PfcControllerDesign.cpp
//
// Implementation of derive_pfc_controller_tuning(). All formulas are derived
// in the design note at the top of PfcControllerDesign.h — that file is the
// canonical reference.
// ─────────────────────────────────────────────────────────────────────────────

#include "converter_models/PfcControllerDesign.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace OpenMagnetics {

namespace {
constexpr double kTwoPi = 6.283185307179586;
constexpr double kSqrt2 = 1.4142135623730951;
constexpr double k2OverPi = 0.6366197723675814;     // 2/π = mean(|sin|)
constexpr double k2Sqrt2OverPi = kSqrt2 * k2OverPi; // (2√2)/π ≈ 0.9003

void require_positive(const char* name, double v) {
    if (!std::isfinite(v) || v <= 0.0) {
        throw std::invalid_argument(
            std::string("derive_pfc_controller_tuning: ") + name +
            " must be positive and finite (got " + std::to_string(v) + ")");
    }
}
} // namespace

PfcControllerTuning derive_pfc_controller_tuning(double vrms_nom,
                                                   double vbus_nom,
                                                   double pout,
                                                   double fsw,
                                                   double fline,
                                                   double inductance) {
    require_positive("vrms_nom",   vrms_nom);
    require_positive("vbus_nom",   vbus_nom);
    require_positive("pout",       pout);
    require_positive("fsw",        fsw);
    require_positive("fline",      fline);
    require_positive("inductance", inductance);

    // Boost-PFC requires Vbus > Vpk (otherwise no boost action is possible).
    const double vpk = kSqrt2 * vrms_nom;
    if (vbus_nom <= vpk) {
        throw std::runtime_error(
            "derive_pfc_controller_tuning: Vbus_nom (" + std::to_string(vbus_nom) +
            " V) must exceed Vin_peak (" + std::to_string(vpk) +
            " V = √2·Vrms) for boost-PFC operation");
    }

    PfcControllerTuning t;

    // ── Block 1: voltage error amp (Gv) ─────────────────────────────────
    // Industry rule: fc_v ≤ fline / 5 to keep loop bandwidth well below
    // the natural 2·fline ripple on Cbus (avoids 3rd-harmonic input
    // current distortion). 5 Hz floor handles non-standard low-fline tests.
    t.vref     = 2.5;
    t.k_div    = t.vref / vbus_nom;
    t.gv_fc_hz = std::max(5.0, fline / 5.0);
    t.gv_rz    = 100e3;
    // 60° phase boost: zero at fc/3, pole at 3·fc.
    t.gv_cz    = 1.0 / (kTwoPi * (t.gv_fc_hz / 3.0) * t.gv_rz);
    t.gv_cp    = 1.0 / (kTwoPi * (3.0      * t.gv_fc_hz) * t.gv_rz);
    t.vea_nom  = 2.5;
    t.vea_min  = 0.1;
    t.vea_max  = 4.9;

    // ── Block 2: RMS² feed-forward (vrms_ff) ────────────────────────────
    t.ff_fc_hz = std::max(5.0, fline / 5.0);
    t.ff_r     = 100e3;
    t.ff_c     = 1.0 / (kTwoPi * t.ff_fc_hz * t.ff_r);
    t.k_ff_sq  = k2Sqrt2OverPi * k2Sqrt2OverPi;            // ≈ 0.811

    // ── Block 3: multiplier ─────────────────────────────────────────────
    // i_ref(t) = G_mul · vea · vin_rect(t) / vrms_ff
    // At rated conditions: i_ref_env_peak = G_mul · vea_nom · √2 / (k_ff² · Vrms_nom)
    // Set this equal to desired envelope peak √2 · Pout / Vrms_nom →
    //   G_mul = k_ff² · Pout / vea_nom
    t.g_mul          = t.k_ff_sq * pout / t.vea_nom;
    // Numerical floor on divisor (1 % of nominal vrms_ff) — defensive only,
    // vrms_ff stays approximately constant across the line cycle so this
    // clamp should never engage in steady state.
    t.vrms_ff_floor  = 0.01 * t.k_ff_sq * vrms_nom * vrms_nom;

    // ── Block 4: current error amp (Gi) ─────────────────────────────────
    // Plant: Gid(s) = Vbus / (s·L); PWM gain = 1/V_pk_saw.
    // Loop-crossover Kp: |Tu(j·2π·fc_i)| = Kp · Vbus / (2π·fc_i·L·V_pk_saw) = 1
    // → Kp = 2π·fc_i·L·V_pk_saw / Vbus
    t.pwm_v_pk_saw = 1.0;
    t.gi_fc_hz     = fsw / 10.0;
    if (t.gi_fc_hz >= fsw / 2.0) {
        throw std::runtime_error(
            "derive_pfc_controller_tuning: current-loop crossover ("+
            std::to_string(t.gi_fc_hz)+
            " Hz) violates Nyquist (fsw/2 = "+std::to_string(fsw/2.0)+" Hz)");
    }
    t.gi_kp        = kTwoPi * t.gi_fc_hz * inductance * t.pwm_v_pk_saw / vbus_nom;
    t.gi_rz        = 100e3;
    t.gi_rin       = t.gi_rz / t.gi_kp;
    t.gi_cz        = 1.0 / (kTwoPi * (t.gi_fc_hz / 3.0) * t.gi_rz);
    t.gi_cp        = 1.0 / (kTwoPi * (3.0     * t.gi_fc_hz) * t.gi_rz);
    t.rs_sense     = 1.0;

    // ── Block 5–6: oscillator + comparator ──────────────────────────────
    t.pwm_v_high   = 5.0;
    t.pwm_t_rise   = 1e-9;

    // ── Initial conditions (warm start) ─────────────────────────────────
    // ic_vbus = Vbus_nom: emulates the inrush-bypass-relay + soft-start
    // sequence that real PFC ICs perform in hardware. Initialising Cbus
    // to Vbus_nom keeps the boost diode D1 reverse-biased from t=0+
    // (since vsw < Vbus_nom always until the gate switches), which
    // eliminates the uncontrolled bridge-inrush current pulses that
    // otherwise flow through L1→D1→Cbus each time vin_rect > vbus
    // during the first few line cycles. Those inrush pulses (3–4 A
    // peak) drive the current EA's integrator capacitor C_fb_zi
    // irreversibly negative within the first half-cycle, dead-locking
    // vc_i at 0 V — gate stays OFF for the entire simulation.
    //
    // With ic_vbus = Vbus_nom, the controller has roughly Cbus·Vbus²/2P
    // = ~(Cbus·Vbus²)/(2·Pload) ≈ 30 ms (for the 100 W design) of
    // "grace period" during which Cbus discharges from Vbus_nom toward
    // Vpk through Rload alone — long enough for the boost loop to
    // engage and start sourcing the load power before D1 ever conducts.
    //
    // Using ic_vbus = Vpk was the previous "physical" choice, but the
    // bridge cannot sustain Vbus > Vpk on its own — the moment
    // Cbus < Vpk, an uncontrolled inrush surge flows through L1+D1
    // each line peak. That inrush deadlocks the EA. The relay-emulated
    // pre-charge is exactly how UCC28019 / NCP1654 / L4981 production
    // boards solve this in hardware (NTC + bypass relay + soft-start
    // pin); we just collapse it to a single .ic line.
    t.ic_vbus = vbus_nom;
    t.ic_vea  = t.vea_nom;
    t.ic_vff  = k2Sqrt2OverPi * vrms_nom;       // mean(|vin|)
    // D_target = duty required to hold Vbus exactly at Vbus_nom at the
    // line-voltage peak (worst-case, lowest duty). This is the duty the
    // soft-start floor will hold initially: Vbus = Vpk/(1−D) → D = 1 − Vpk/Vbus.
    // Using D_avg here would over-pump Cbus to Vpk/(1−D_avg) — observed
    // empirically to drive Vbus to ~630 V for a 400 V design.
    const double d_target = 1.0 - vpk / vbus_nom;
    t.ic_vc_i = d_target * t.pwm_v_pk_saw;
    t.t_ss_release = 30e-3;

    return t;
}

} // namespace OpenMagnetics
