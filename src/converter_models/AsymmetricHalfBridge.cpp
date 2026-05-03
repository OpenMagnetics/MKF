#include "converter_models/AsymmetricHalfBridge.h"
#include "support/Exceptions.h"
#include <cmath>
#include <stdexcept>

namespace OpenMagnetics {

// =========================================================================
// Constructors
// =========================================================================
AsymmetricHalfBridge::AsymmetricHalfBridge(const json& j) {
    from_json(j, *static_cast<MAS::AsymmetricHalfBridge*>(this));
}

AdvancedAsymmetricHalfBridge::AdvancedAsymmetricHalfBridge(const json& j) {
    from_json(j, *this);
}


// =========================================================================
// Input validation (P1 deliverable)
//
// Mirrors PSHB / PSFB run_checks: reject empty operating-point lists,
// per-OP missing voltages/currents, out-of-range duty cycle, non-positive
// switching frequency, invalid rectifier type. Per the project rule "no
// fallbacks, no defaults, no silent shortcuts": surface every malformed
// field as a std::runtime_error rather than substituting a default.
// =========================================================================
bool AsymmetricHalfBridge::run_checks(bool assertErrors) {
    _assertErrors = assertErrors;
    bool ok = true;

    auto& ops = get_operating_points();
    if (ops.empty()) {
        if (assertErrors) throw std::runtime_error(
            "AsymmetricHalfBridge: no operating points");
        return false;
    }

    // Input voltage must be present with at least one of nominal/min/max.
    auto& vIn = get_input_voltage();
    if (!vIn.get_nominal().has_value() &&
        !vIn.get_minimum().has_value() &&
        !vIn.get_maximum().has_value()) {
        if (assertErrors) throw std::runtime_error(
            "AsymmetricHalfBridge: inputVoltage must specify nominal, "
            "minimum, or maximum");
        ok = false;
    }

    // maximumDutyCycle, when present, must lie in (0, 1). The non-monotonic
    // gain peaks at D=0.5; conventional control uses D in [0, 0.5], but
    // the schema allows D up to 1.0 to support the upper-branch researchers.
    if (auto dMax = get_maximum_duty_cycle()) {
        if (*dMax <= 0.0 || *dMax >= 1.0) {
            if (assertErrors) throw std::runtime_error(
                "AsymmetricHalfBridge: maximumDutyCycle must lie in (0, 1)");
            ok = false;
        }
    }

    for (size_t i = 0; i < ops.size(); ++i) {
        auto& op = ops[i];
        const std::string opTag = "AsymmetricHalfBridge: OP[" +
                                  std::to_string(i) + "]";

        if (op.get_output_voltages().empty() ||
            op.get_output_currents().empty()) {
            if (assertErrors) throw std::runtime_error(
                opTag + " missing output voltages/currents");
            ok = false;
            continue;
        }
        if (op.get_output_voltages().size() != op.get_output_currents().size()) {
            if (assertErrors) throw std::runtime_error(
                opTag + " output voltage / current vectors length mismatch");
            ok = false;
        }

        if (op.get_switching_frequency() <= 0) {
            if (assertErrors) throw std::runtime_error(
                opTag + " switching frequency must be > 0");
            ok = false;
        }

        // dutyCycle is REQUIRED by the AHB schema and bounded [0, 1].
        // The schema-level constraint is enforced by quicktype, but we
        // double-check here so a hand-constructed object is also caught.
        double D = op.get_duty_cycle();
        if (D < 0.0 || D > 1.0) {
            if (assertErrors) throw std::runtime_error(
                opTag + " dutyCycle must lie in [0, 1]");
            ok = false;
        }

        // Multi-output cross-regulation in AHB is materially worse than
        // forward (asymmetric secondary voltage). Flag it loud now so the
        // P11 warn-once message is not the first time the user sees it.
        if (op.get_output_voltages().size() > 1) {
            // Soft warning only - actual warn-once lives in P11.
            // Do NOT downgrade ok.
        }
    }

    // rectifierType, when present, must be one of the four enum values.
    // (Trivially true since AhbRectifierType is closed; included for
    // symmetry with PSFB / PSHB run_checks.)
    if (auto rt = get_rectifier_type()) {
        switch (*rt) {
            case AhbRectifierType::AHB_FLYBACK:
            case AhbRectifierType::CENTER_TAPPED:
            case AhbRectifierType::CURRENT_DOUBLER:
            case AhbRectifierType::FULL_BRIDGE:
                break;
            default:
                if (assertErrors) throw std::runtime_error(
                    "AsymmetricHalfBridge: invalid rectifierType");
                ok = false;
        }
    }

    return ok;
}


// =========================================================================
// EVERYTHING BELOW IS A P2-P12 STUB. Throws std::runtime_error with a
// pointer to the implementation plan so callers do not silently get
// zero-initialised garbage. Per project rule: no defaults, no silent
// shortcuts, throw loud.
// =========================================================================

namespace {
[[noreturn]] void not_implemented(const char* method, const char* phase) {
    throw std::runtime_error(
        std::string("AsymmetricHalfBridge: ") + method +
        " not implemented yet (see ASYMMETRIC_HALF_BRIDGE_PLAN.md " +
        phase + ")");
}
}

DesignRequirements AsymmetricHalfBridge::process_design_requirements() {
    not_implemented("process_design_requirements", "P6");
}

std::vector<OperatingPoint> AsymmetricHalfBridge::process_operating_points(
    const std::vector<double>& /*turnsRatios*/,
    double /*magnetizingInductance*/) {
    not_implemented("process_operating_points", "P3");
}

std::vector<OperatingPoint> AsymmetricHalfBridge::process_operating_points(
    Magnetic /*magnetic*/) {
    not_implemented("process_operating_points(Magnetic)", "P3");
}

OperatingPoint AsymmetricHalfBridge::process_operating_point_for_input_voltage(
    double /*inputVoltage*/,
    const AhbOperatingPoint& /*opPoint*/,
    const std::vector<double>& /*turnsRatios*/,
    double /*magnetizingInductance*/) {
    not_implemented("process_operating_point_for_input_voltage", "P3");
}

// =========================================================================
// Static analytical helpers (P2 deliverable)
//
// All equation references in §-form are to ASYMMETRIC_HALF_BRIDGE_PLAN.md
// §4 ("Equations (canonical block for the new header docstring)").
// Notation:
//   Vin = input bus voltage [V]
//   Vo  = output voltage [V]
//   Vd  = rectifier voltage drop [V]
//   D   = duty cycle of the high-side switch Q1 (Q2 runs at 1-D) [-]
//   n   = primary-to-secondary turns ratio Np/Ns [-]
//          (half-secondary winding for CT, full secondary for FB / CD)
//   Tsw = 1/fsw [s], fsw = switching frequency [Hz]
//   Lo  = output inductor [H]
//   Lm  = primary magnetizing inductance [H]
//   Llk = primary leakage inductance [H]
//   Cb  = DC-blocking capacitance [F]
//   Coss = per-MOSFET output capacitance [F]
//   Np  = primary turn count [-]
//   Ae  = effective core cross-section [m^2]
//
// Per project rule "no fallbacks, no defaults, no silent shortcuts": every
// out-of-domain or sign-violating input throws std::invalid_argument. No
// clamping, no value_or, no "physically plausible minimum" substitution.
// =========================================================================

namespace {
[[noreturn]] void invalid(const char* helper, const std::string& msg) {
    throw std::invalid_argument(
        std::string("AsymmetricHalfBridge::") + helper + ": " + msg);
}
inline void require_positive(const char* helper, const char* name, double v) {
    if (!(v > 0.0) || !std::isfinite(v))
        invalid(helper, std::string(name) + " must be > 0 (got " +
                std::to_string(v) + ")");
}
inline void require_in_unit_interval(const char* helper, const char* name,
                                     double v) {
    if (!(v >= 0.0 && v <= 1.0) || !std::isfinite(v))
        invalid(helper, std::string(name) + " must lie in [0, 1] (got " +
                std::to_string(v) + ")");
}
} // namespace


// -------------------------------------------------------------------------
// (E1) DC-blocking cap voltage in steady state.
//        V_Cb = (1 - D) * Vin
// Reference: Imbertson & Mohan IEEE TIA 29(1) 1993, eq. 3.
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_dc_blocking_cap_voltage(double Vin,
                                                             double D) {
    require_positive("compute_dc_blocking_cap_voltage", "Vin", Vin);
    require_in_unit_interval("compute_dc_blocking_cap_voltage", "D", D);
    return (1.0 - D) * Vin;
}


// -------------------------------------------------------------------------
// (E5)/(E6) Conversion ratio Vo/Vin (lossless, ideal rectifier).
//   CT, FB:        M = 2 * D * (1 - D) / n
//   CD:            M =     D * (1 - D) / n     (factor of 2 absorbed by doubler)
//   AHB-Flyback:   M = D / ((1 - D) * n)        (Imbertson-Mohan flyback variant)
// The CT/FB/CD curves are NON-MONOTONIC in D, peaking at D=0.5
// with M_max = 1/(2n) (CT/FB) or 1/(4n) (CD). Conventional control
// uses D in [0, 0.5] (the rising branch).
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_conversion_ratio(double D, double n,
                                                      AhbRectifierType rect) {
    require_in_unit_interval("compute_conversion_ratio", "D", D);
    require_positive("compute_conversion_ratio", "n", n);
    const double DDc = D * (1.0 - D);
    switch (rect) {
        case AhbRectifierType::CENTER_TAPPED:
        case AhbRectifierType::FULL_BRIDGE:
            return 2.0 * DDc / n;
        case AhbRectifierType::CURRENT_DOUBLER:
            return DDc / n;
        case AhbRectifierType::AHB_FLYBACK:
            // Pseudo-flyback gain: Vo/Vin = D / ((1-D) * n). Singular at
            // D=1; bounded for D in [0, 1).
            if (D >= 1.0)
                invalid("compute_conversion_ratio",
                        "AHB_FLYBACK gain singular at D=1");
            return D / ((1.0 - D) * n);
    }
    invalid("compute_conversion_ratio", "unknown AhbRectifierType");
}


// -------------------------------------------------------------------------
// (S1) Turns ratio sized at Vin_min, D_max so M(D_max)*Vin_min = Vo + Vd.
//        n = M_topology(D_max) * Vin_min / (Vo + Vd_total)
// where Vd_total accounts for diode count: 1 diode in CT path, 2 diodes
// in FB path (each half-cycle goes through 2 diodes), 1 diode in CD path
// (each inductor sees 1 conducting diode at a time), 1 diode in flyback.
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_turns_ratio(double Vin_min, double Vo,
                                                 double D_max, double Vd,
                                                 AhbRectifierType rect) {
    require_positive("compute_turns_ratio", "Vin_min", Vin_min);
    require_positive("compute_turns_ratio", "Vo", Vo);
    require_in_unit_interval("compute_turns_ratio", "D_max", D_max);
    if (Vd < 0.0)
        invalid("compute_turns_ratio", "Vd must be >= 0");

    const double DDc = D_max * (1.0 - D_max);
    if (DDc <= 0.0)
        invalid("compute_turns_ratio",
                "D_max in {0,1} yields zero conversion ratio");

    double Vd_total = 0;
    double gain_factor = 0;
    switch (rect) {
        case AhbRectifierType::CENTER_TAPPED:
            Vd_total    = Vd;
            gain_factor = 2.0 * DDc;
            break;
        case AhbRectifierType::FULL_BRIDGE:
            Vd_total    = 2.0 * Vd;
            gain_factor = 2.0 * DDc;
            break;
        case AhbRectifierType::CURRENT_DOUBLER:
            Vd_total    = Vd;
            gain_factor = DDc;
            break;
        case AhbRectifierType::AHB_FLYBACK:
            // Vo + Vd = (D / ((1-D)*n)) * Vin_min
            //   => n  = (D / (1-D)) * Vin_min / (Vo + Vd)
            if (D_max >= 1.0)
                invalid("compute_turns_ratio",
                        "AHB_FLYBACK: D_max=1 is singular");
            return (D_max / (1.0 - D_max)) * Vin_min / (Vo + Vd);
    }
    return gain_factor * Vin_min / (Vo + Vd_total);
}


// -------------------------------------------------------------------------
// Inverse of (E5): solve for the duty cycle D in [0, 0.5] that yields a
// target conversion ratio M = Vo / Vin (lossy form: M_eff = (Vo+Vd) /
// (Vin * eta)). For CT/FB:
//        2 D (1-D) / n  =  M_eff
//   =>   D^2 - D + M_eff*n/2 = 0
//   =>   D = 0.5 - sqrt(0.25 - M_eff*n/2)            (lower branch)
// For CD: same shape with M_eff -> 2*M_eff. For AHB-Flyback: closed form
//        M_eff = D / ((1-D) * n)  =>  D = M_eff * n / (1 + M_eff * n)
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_duty_cycle(double Vin, double Vo,
                                                double n, double Vd,
                                                double efficiency,
                                                AhbRectifierType rect) {
    require_positive("compute_duty_cycle", "Vin", Vin);
    require_positive("compute_duty_cycle", "Vo", Vo);
    require_positive("compute_duty_cycle", "n", n);
    if (Vd < 0.0)
        invalid("compute_duty_cycle", "Vd must be >= 0");
    if (!(efficiency > 0.0 && efficiency <= 1.0))
        invalid("compute_duty_cycle", "efficiency must lie in (0, 1]");

    // Effective gain target accounting for diode drops + efficiency. The
    // efficiency factor reduces deliverable Vo, equivalently raises the
    // required gain.
    const double Vd_total =
        (rect == AhbRectifierType::FULL_BRIDGE) ? 2.0 * Vd : Vd;
    const double M_eff = (Vo + Vd_total) / (Vin * efficiency);

    if (rect == AhbRectifierType::AHB_FLYBACK) {
        const double Mn = M_eff * n;
        return Mn / (1.0 + Mn);
    }

    // For CT/FB:  D^2 - D + M_eff*n/2 = 0
    // For CD:     D^2 - D + M_eff*n   = 0
    const double c =
        (rect == AhbRectifierType::CURRENT_DOUBLER) ? M_eff * n
                                                    : M_eff * n * 0.5;
    const double disc = 0.25 - c;
    if (disc < 0.0)
        invalid("compute_duty_cycle",
                "target gain " + std::to_string(M_eff) +
                " exceeds the AHB peak 1/(2n) — increase n or decrease Vo");
    return 0.5 - std::sqrt(disc);  // lower / "rising" branch
}


// -------------------------------------------------------------------------
// (S2) Output inductor sized to bound peak-to-peak inductor ripple dILo.
// During the rectifier-OFF freewheel sub-interval, V_Lo = -(Vo+Vd) and
// the freewheel duration is (1 - 2D(1-D))*Tsw (the fraction of the
// period during which BOTH secondary half-cycles are off; the secondary
// "sees" 2D(1-D) duty after the centre-tap rectification).
//        Lo_min = Vo * (1 - 2 * D * (1-D)) / (dILo * fsw)
// At D=0.5 this collapses to Lo_min = 0.5 * Vo / (dILo * fsw) (the
// freewheel duration is half the period). dILo is the absolute peak-to-
// peak current ripple in [A] (NOT a ratio).
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_lo_min(double Vo, double D,
                                            double dILo, double fsw) {
    require_positive("compute_lo_min", "Vo",  Vo);
    require_in_unit_interval("compute_lo_min", "D", D);
    require_positive("compute_lo_min", "dILo", dILo);
    require_positive("compute_lo_min", "fsw",  fsw);
    return Vo * (1.0 - 2.0 * D * (1.0 - D)) / (dILo * fsw);
}


// -------------------------------------------------------------------------
// (S4) Magnetizing-inductance upper bound to deliver a target magnetizing-
// current ripple Im_target sufficient for ZVS at minimum load.
//   v_pri = (1-D)*Vin during D*Tsw, di_Lm/dt = (1-D)*Vin/Lm
//   Im_pp,target = ((1-D)*Vin*D*Tsw) / Lm
//   => Lm <= (1-D) * Vin * D * Tsw / (2 * Im_target)        (peak = pp/2)
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_lm_min_for_zvs(double Vin, double D,
                                                    double Tsw,
                                                    double Im_target) {
    require_positive("compute_lm_min_for_zvs", "Vin", Vin);
    require_in_unit_interval("compute_lm_min_for_zvs", "D", D);
    require_positive("compute_lm_min_for_zvs", "Tsw", Tsw);
    require_positive("compute_lm_min_for_zvs", "Im_target", Im_target);
    return ((1.0 - D) * Vin * D * Tsw) / (2.0 * Im_target);
}


// -------------------------------------------------------------------------
// (S3) DC-blocking capacitor sized to bound steady-state ripple dV_Cb.
// Charge balance: Q_Cb = I_pri,pk * D * Tsw flows in/out per half-cycle.
//        Cb_min = I_pri,pk * D / (fsw * dV_Cb)
// ON Semi AN-4153 eq. 16 recommends dV_Cb <= 5% of V_Cb steady state.
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_cb_min(double Iprimary_pk, double D,
                                            double dVCb, double fsw) {
    require_positive("compute_cb_min", "Iprimary_pk", Iprimary_pk);
    require_in_unit_interval("compute_cb_min", "D", D);
    require_positive("compute_cb_min", "dVCb", dVCb);
    require_positive("compute_cb_min", "fsw",  fsw);
    return (Iprimary_pk * D) / (fsw * dVCb);
}


// -------------------------------------------------------------------------
// (Z1) ZVS energy-balance margin. ZVS occurs iff:
//        0.5 * (Llk + Lm_refl) * Ipri^2  >=  2 * Coss * Vin^2
// Returns the LHS - RHS difference in joules. A non-negative result means
// ZVS is achievable; negative means hard switching. The "2" on the RHS
// accounts for both MOSFETs' Coss being charged/discharged through the
// resonant tank during the dead-time. Lm_refl is the Lm value reflected
// to the primary side (already in primary-side units; pass 0 to ignore
// magnetizing-current contribution).
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_zvs_energy_balance(double Llk,
                                                        double Lm_refl,
                                                        double Ipri,
                                                        double Coss,
                                                        double Vin) {
    require_positive("compute_zvs_energy_balance", "Llk",  Llk);
    if (Lm_refl < 0.0)
        invalid("compute_zvs_energy_balance", "Lm_refl must be >= 0");
    require_positive("compute_zvs_energy_balance", "Coss", Coss);
    require_positive("compute_zvs_energy_balance", "Vin",  Vin);
    if (Ipri < 0.0)
        invalid("compute_zvs_energy_balance", "Ipri must be >= 0");
    const double E_stored = 0.5 * (Llk + Lm_refl) * Ipri * Ipri;
    const double E_charge = 2.0 * Coss * Vin * Vin;
    return E_stored - E_charge;
}


// -------------------------------------------------------------------------
// (Z2) Dead-time matched to the resonant ZVS transition period.
//        td = pi * sqrt(Llk * 2 * Coss)
// (Quarter-period of the LC tank formed by Llk and the parallel Coss
// of both switches.) Mirrors the PSFB dead-time sizing.
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_dead_time(double Llk, double Coss) {
    require_positive("compute_dead_time", "Llk",  Llk);
    require_positive("compute_dead_time", "Coss", Coss);
    return M_PI * std::sqrt(Llk * 2.0 * Coss);
}


// -------------------------------------------------------------------------
// (F4) Steady-state peak flux density.
// Volt-second product per half-cycle: lambda+ = (1-D)*Vin * D*Tsw.
// Peak flux: Bpk = lambda+ / (2 * Np * Ae)   (BIPOLAR core utilisation)
//          = D * (1-D) * Vin * Tsw / (2 * Np * Ae)
// Peaks at D=0.5; degenerates as D -> 0 or D -> 1. Returns Bpk [T].
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_steady_state_flux_excursion(
    double Vin, double D, double Tsw, double Np, double Ae) {
    require_positive("compute_steady_state_flux_excursion", "Vin", Vin);
    require_in_unit_interval("compute_steady_state_flux_excursion", "D", D);
    require_positive("compute_steady_state_flux_excursion", "Tsw", Tsw);
    require_positive("compute_steady_state_flux_excursion", "Np",  Np);
    require_positive("compute_steady_state_flux_excursion", "Ae",  Ae);
    return D * (1.0 - D) * Vin * Tsw / (2.0 * Np * Ae);
}


// -------------------------------------------------------------------------
// (T1)/(T2) Worst-case TRANSIENT peak flux during a Vin step.
// V_Cb cannot follow Vin instantaneously (its time-constant is set by
// Lm/Rdcr); during the transient, volt-second balance (E4) is violated
// and an asymmetric volt-second adder appears on the primary:
//        delta_lambda_transient = D * dVin * Tsw    (worst case)
// Peak flux adder:
//        delta_B_transient      = D * dVin * Tsw / (2 * Np * Ae)
// Returned value is the TOTAL worst-case Bpk (steady-state + transient
// adder), so a magnetic adviser can compare directly against Bsat.
//
// Simplification: for a positive Vin step we use the worst-half-cycle
// adder (no exponential envelope). The Lm/Rdcr time-constant (Lm * Rdcr
// units of seconds) is reported to the caller via lastTransientFlux* so
// that controller designers can tune brown-out recovery; here we assume
// the transient lasts at least one switching cycle (the relevant case
// for core saturation).
// -------------------------------------------------------------------------
double AsymmetricHalfBridge::compute_transient_flux_excursion(
    double Vin, double dVin, double D, double Tsw, double Np, double Ae,
    double Lm, double Rdcr) {
    require_positive("compute_transient_flux_excursion", "Vin",  Vin);
    if (dVin < 0.0)
        invalid("compute_transient_flux_excursion", "dVin must be >= 0");
    require_in_unit_interval("compute_transient_flux_excursion", "D", D);
    require_positive("compute_transient_flux_excursion", "Tsw",  Tsw);
    require_positive("compute_transient_flux_excursion", "Np",   Np);
    require_positive("compute_transient_flux_excursion", "Ae",   Ae);
    require_positive("compute_transient_flux_excursion", "Lm",   Lm);
    require_positive("compute_transient_flux_excursion", "Rdcr", Rdcr);

    // Time-constant of the V_Cb relaxation (informational; not used in
    // the worst-case bound). Surfacing this in lastTransient* is a P5
    // deliverable.
    (void)(Lm / Rdcr);

    const double Bpk_ss     = D * (1.0 - D) * Vin * Tsw / (2.0 * Np * Ae);
    const double dB_trans   = D * dVin * Tsw          / (2.0 * Np * Ae);
    return Bpk_ss + dB_trans;
}

std::vector<std::variant<Inputs, CAS::Inputs>>
AsymmetricHalfBridge::get_extra_components_inputs(
    ExtraComponentsMode /*mode*/, std::optional<Magnetic> /*magnetic*/) {
    not_implemented("get_extra_components_inputs", "P9");
}

std::string AsymmetricHalfBridge::generate_ngspice_circuit(
    const std::vector<double>& /*turnsRatios*/,
    double /*magnetizingInductance*/,
    size_t /*inputVoltageIndex*/,
    size_t /*operatingPointIndex*/) {
    not_implemented("generate_ngspice_circuit", "P7");
}

std::vector<OperatingPoint>
AsymmetricHalfBridge::simulate_and_extract_operating_points(
    const std::vector<double>& /*turnsRatios*/,
    double /*magnetizingInductance*/) {
    not_implemented("simulate_and_extract_operating_points", "P7");
}

std::vector<ConverterWaveforms>
AsymmetricHalfBridge::simulate_and_extract_topology_waveforms(
    const std::vector<double>& /*turnsRatios*/,
    double /*magnetizingInductance*/,
    size_t /*numberOfPeriods*/) {
    not_implemented("simulate_and_extract_topology_waveforms", "P7");
}

Inputs AdvancedAsymmetricHalfBridge::process() {
    not_implemented("AdvancedAsymmetricHalfBridge::process", "P12");
}

} // namespace OpenMagnetics
