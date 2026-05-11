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
    t.ic_vbus = vbus_nom;
    t.ic_vea  = t.vea_nom;
    t.ic_vff  = k2Sqrt2OverPi * vrms_nom;       // mean(|vin|)

    return t;
}

} // namespace OpenMagnetics
