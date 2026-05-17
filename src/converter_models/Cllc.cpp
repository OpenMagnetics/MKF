#include "converter_models/Cllc.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/LeakageInductance.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <iostream>
#include <complex>
#include <limits>
#include "support/Exceptions.h"

namespace OpenMagnetics {

// Helper to extract value from either double or std::optional<double>
template<typename T>
double get_value_or(T&& val, double default_val) {
    if constexpr (std::is_same_v<std::decay_t<T>, std::optional<double>>) {
        return val.value_or(default_val);
    } else {
        return val;
    }
}

    // =========================================================================
    // Construction
    // =========================================================================

    CllcConverter::CllcConverter(const json& j) {
        // Parse base CllcResonant fields (v1)
        set_bidirectional(get_stack_optional<bool>(j, "bidirectional"));
        set_efficiency(get_stack_optional<double>(j, "efficiency"));
        set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
        set_max_switching_frequency(j.at("maxSwitchingFrequency").get<double>());
        set_min_switching_frequency(j.at("minSwitchingFrequency").get<double>());
        set_operating_points(j.at("operatingPoints").get<std::vector<CllcOperatingPoint>>());
        set_quality_factor(get_stack_optional<double>(j, "qualityFactor"));
        set_symmetric_design(get_stack_optional<bool>(j, "symmetricDesign"));

        // v2 schema additions (all optional; absence preserves defaults).
        // bridgeType: shares LlcBridgeType (deduplicated by quicktype since
        // the JSON enum values fullBridge/halfBridge are identical).
        set_bridge_type(get_stack_optional<LlcBridgeType>(j, "bridgeType"));
        set_integrated_resonant_inductor1(get_stack_optional<bool>(j, "integratedResonantInductor1"));
        set_integrated_resonant_inductor2(get_stack_optional<bool>(j, "integratedResonantInductor2"));
        set_resonant_inductor_ratio(get_stack_optional<double>(j, "resonantInductorRatio"));
        set_resonant_capacitor_ratio(get_stack_optional<double>(j, "resonantCapacitorRatio"));

        // v1 → v2 derived alias: when symmetricDesign was set explicitly to
        // false in v1 input but neither ratio is supplied, leave both at
        // their schema defaults (1.0). New code should set the ratios
        // directly. (No-op block; documented for grep.)
    }

    AdvancedCllcConverter::AdvancedCllcConverter(const json& j) {
        from_json(j, *this);
    }

    bool CllcConverter::run_checks(bool assert) {
        // Schema-level checks first (frequency window, OP non-empty, etc.).
        bool ok = Topology::run_checks_common(this, assert);

        // v2 field validation
        auto a = get_resonant_inductor_ratio();
        if (a.has_value() && (*a <= 0.0 || *a > 10.0)) {
            ok = false;
            if (assert) {
                throw std::invalid_argument(
                    "CLLC resonantInductorRatio (a) out of range (0, 10]: "
                    + std::to_string(*a));
            }
        }
        auto b = get_resonant_capacitor_ratio();
        if (b.has_value() && (*b <= 0.0 || *b > 10.0)) {
            ok = false;
            if (assert) {
                throw std::invalid_argument(
                    "CLLC resonantCapacitorRatio (b) out of range (0, 10]: "
                    + std::to_string(*b));
            }
        }
        // a·b ≈ 1 keeps fr1 = fr2 (Min IEEE TPE 2021). Warn if far off, but
        // don't fail — some unusual designs intentionally split the
        // resonant frequencies.
        if (a.has_value() && b.has_value()) {
            double product = (*a) * (*b);
            if (product < 0.5 || product > 2.0) {
                if (assert) {
                    throw std::invalid_argument(
                        "CLLC asymmetric ratios a*b should be ~1 (got "
                        + std::to_string(product)
                        + "); fr1 ≠ fr2 will result.");
                }
            }
        }
        return ok;
    }

    double CllcConverter::get_bridge_voltage_factor() const {
        auto bt = get_bridge_type();
        if (bt.has_value() && *bt == LlcBridgeType::HALF_BRIDGE) {
            return 0.5;
        }
        return 1.0; // FULL_BRIDGE (default)
    }

    void CllcConverter::resolve_resonant_ratios(double& a, double& b) const {
        auto aOpt = get_resonant_inductor_ratio();
        auto bOpt = get_resonant_capacitor_ratio();

        // v1 fallback when neither v2 ratio is set: symmetricDesign drives
        // the (a, b) pair. symmetricDesign==true (default) → (1, 1).
        // symmetricDesign==false → Infineon AN-2024-06 §2.3 typical
        // asymmetric values (0.95, 1.052).
        bool symmetric = get_symmetric_design().value_or(true);
        a = aOpt.has_value() ? *aOpt : (symmetric ? 1.0 : 0.95);
        b = bOpt.has_value() ? *bOpt : (symmetric ? 1.0 : 1.052);
    }

    // =========================================================================
    // Resonant Tank Parameter Calculation
    // =========================================================================

    /**
     * Implements the design procedure from Infineon AN-2024-06, Section 2.3.
     *
     * The CLLC resonant tank is defined by:
     *   Primary side:   C1 (series capacitor) + L1 (series inductor, = leakage inductance)
     *   Transformer:    Lm (magnetizing inductance)
     *   Secondary side: L2 (series inductor) + C2 (series capacitor)
     *
     * The resonant frequency is: fr = 1/(2π√(L1·C1)) = 1/(2π√(L2·C2))
     *
     * For a symmetric tank: L2 = L1/n², C2 = n²·C1  (so fr1 = fr2)
     */
    CllcResonantParameters CllcConverter::calculate_resonant_parameters() {
        CllcResonantParameters params;

        double nominalInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::NOMINAL);
        if (nominalInputVoltage == 0) {
            nominalInputVoltage = (resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM) +
                                   resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM)) / 2.0;
        }

        // Find representative output voltage and power from operating points
        // Use the first forward-mode operating point, or the first one available
        double nominalOutputVoltage = 0;
        double nominalOutputPower = 0;
        double representativeFrequency = 0;

        for (const auto& op : get_operating_points()) {
            if (op.get_output_voltages().size() > 0) {
                nominalOutputVoltage = op.get_output_voltages()[0];
                double current = op.get_output_currents().size() > 0 ? op.get_output_currents()[0] : 0;
                nominalOutputPower = nominalOutputVoltage * current;
                representativeFrequency = get_value_or(op.get_switching_frequency(), 0.0);
                break;
            }
        }

        if (nominalOutputVoltage == 0 || nominalOutputPower == 0) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS,
                "CLLC: Cannot determine output voltage/power from operating points");
        }

        // Step 1: Transformer turns ratio (Infineon AN Step 1)
        // n = Vin_nominal / Vout_nominal
        params.turnsRatio = nominalInputVoltage / nominalOutputVoltage;

        // Use a representative resonant frequency:
        // If operating point frequency is within min/max range, use it;
        // otherwise use the geometric mean of min and max switching frequency
        double fr;
        if (representativeFrequency >= get_min_switching_frequency() &&
            representativeFrequency <= get_max_switching_frequency()) {
            fr = representativeFrequency;
        } else {
            fr = sqrt(get_min_switching_frequency() * get_max_switching_frequency());
        }
        params.resonantFrequency = fr;

        // Step 4: Effective AC resistance (FHA) (Infineon AN Step 4)
        // Ro = (8·n²/π²) · (Vout²/Pout) = (8·n²/π²) · R_load
        double n = params.turnsRatio;
        double Rload = nominalOutputVoltage * nominalOutputVoltage / nominalOutputPower;
        double Ro = (8.0 * n * n / (M_PI * M_PI)) * Rload;
        params.equivalentAcResistance = Ro;

        // Step 3: Quality factor Q (use user-specified or default)
        double Q;
        if (get_quality_factor()) {
            Q = get_quality_factor().value();
        } else {
            Q = 0.3;  // Default from Infineon AN recommendation (between 0.2 and 0.4)
        }
        params.qualityFactor = Q;

        // Step 5: Primary resonant capacitor C1 (Infineon AN Step 5)
        // C1 = 1 / (2π·Q·fr·Ro)
        double C1 = 1.0 / (2.0 * M_PI * Q * fr * Ro);
        params.primaryResonantCapacitance = C1;

        // Step 6: Primary resonant inductor L1 (Infineon AN Step 6)
        // L1 = 1 / ((2πfr)²·C1)
        double omega_r = 2.0 * M_PI * fr;
        double L1 = 1.0 / (omega_r * omega_r * C1);
        params.primaryResonantInductance = L1;

        // Step 7: Magnetizing inductance Lm (Infineon AN Step 7)
        // Lm = k · L1
        double k = defaultInductanceRatio;
        double Lm = k * L1;
        params.magnetizingInductance = Lm;
        params.inductanceRatio = k;

        // Steps 8-9: Secondary resonant components — v2 schema ratios
        // (resonantInductorRatio, resonantCapacitorRatio) override the
        // v1 symmetricDesign alias. See resolve_resonant_ratios.
        double a, b;
        resolve_resonant_ratios(a, b);
        params.resonantInductorRatio = a;
        params.resonantCapacitorRatio = b;

        // L2 = a·L1/n² (Infineon AN Step 9)
        params.secondaryResonantInductance = a * L1 / (n * n);

        // C2 = n²·b·C1 (Infineon AN Step 9)
        params.secondaryResonantCapacitance = n * n * b * C1;

        return params;
    }

    // =========================================================================
    // FHA Voltage Gain Calculation
    // =========================================================================

    /**
     * Compute the FHA voltage gain |nVout/Vin| at a given switching frequency.
     *
     * From Infineon AN Equation 2 (for ωr1 = ωr2, i.e., symmetric or
     * when primary and secondary resonant frequencies match):
     *
     *   nVout/Vin = 1 / |Denominator|
     *
     * where the denominator in terms of s = jω is given by Equation 2.
     *
     * This implementation uses the impedance-based form:
     *   H(jω) = Zm·Ro / (Z1·Zm + Z1·Z2 + Z1·Ro + Zm·Z2 + Zm·Ro)
     *
     * where:
     *   Z1 = jωL1 + 1/(jωC1)        (primary resonant impedance)
     *   Z2 = n²·(jωL2 + 1/(jωC2))   (secondary resonant impedance referred to primary)
     *   Zm = jωLm                     (magnetizing impedance)
     *   Ro = 8n²/(π²)·R_load         (FHA equivalent AC load resistance)
     */
    double CllcConverter::get_voltage_gain(double switchingFrequency,
                                            const CllcResonantParameters& params) {
        double omega_s = 2.0 * M_PI * switchingFrequency;
        std::complex<double> j(0.0, 1.0);
        std::complex<double> s = j * omega_s;

        double n = params.turnsRatio;
        double L1 = params.primaryResonantInductance;
        double C1 = params.primaryResonantCapacitance;
        double L2 = params.secondaryResonantInductance;
        double C2 = params.secondaryResonantCapacitance;
        double Lm = params.magnetizingInductance;
        double Ro = params.equivalentAcResistance;

        // Primary resonant impedance: Z1 = sL1 + 1/(sC1)
        std::complex<double> Z1 = s * L1 + 1.0 / (s * C1);

        // Secondary resonant impedance referred to primary: Z2 = n²·(sL2 + 1/(sC2))
        std::complex<double> Z2 = n * n * (s * L2 + 1.0 / (s * C2));

        // Magnetizing impedance: Zm = sLm
        std::complex<double> Zm = s * Lm;

        // Transfer function: H = Zm·Ro / (Z1·Zm + Z1·Z2 + Z1·Ro + Zm·Z2 + Zm·Ro)
        std::complex<double> numerator = Zm * Ro;
        std::complex<double> denominator = Z1 * Zm + Z1 * Z2 + Z1 * Ro + Zm * Z2 + Zm * Ro;

        return std::abs(numerator / denominator);
    }

    // =========================================================================
    // Design Requirements
    // =========================================================================

    DesignRequirements CllcConverter::process_design_requirements() {
        auto params = calculate_resonant_parameters();

        DesignRequirements designRequirements;

        // Turns ratio
        designRequirements.get_mutable_turns_ratios().clear();
        DimensionWithTolerance turnsRatioWithTolerance;
        turnsRatioWithTolerance.set_nominal(roundFloat(params.turnsRatio, 2));
        designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);

        // Magnetizing inductance
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(params.magnetizingInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

        // Isolation sides: primary and one secondary
        designRequirements.set_isolation_sides(
            Topology::create_isolation_sides(1, false));

        designRequirements.set_topology(Topologies::CLLC_RESONANT_CONVERTER);

        return designRequirements;
    }

    // =========================================================================
    // Analytical Waveform Generation
    // =========================================================================

    /**
     * Generate analytical operating point for a single input voltage and operating condition.
     *
     * The CLLC converter waveforms are modeled using FHA:
     *
     * PRIMARY WINDING:
     *   Voltage: Bipolar rectangular ±Vin with dead time
     *            +Vin for [0, T/2 - td], 0 for dead time, -Vin for [T/2, T - td], 0
     *   Current: Sum of resonant current and magnetizing current
     *            - Magnetizing current: triangular, peak = Vin/(4·Lm·fs)
     *            - Resonant current: sinusoidal at frequency ≈ fr
     *            - Total Ip(t) ≈ I_res·sin(2π·fs·t) + Im(t)
     *
     * SECONDARY WINDING:
     *   Voltage: Bipolar rectangular ±Vout*n (reflected to primary, then /n to secondary)
     *   Current: Resonant current minus magnetizing current, scaled by n
     *            Isec(t) = n · (Ip(t) - Im(t))
     *
     * At resonance (fs = fr):
     *   - Primary current is quasi-sinusoidal
     *   - Secondary current completes exactly one half-sine per half period
     *   - Peak current: Ip_peak ≈ Pout/(Vin·η) · π/2 (fundamental approximation)
     */
    // =========================================================================
    // P2 — Time-Domain Analysis (TDA) solver
    //
    // For the symmetric tank (a = b = 1) the CLLC 5-element resonant network
    //   Cr1 — Lr1 — Lm — Lr2 — Cr2
    // collapses (per CLLC_REWRITE_PLAN.md §5.1) to an LLC-equivalent with
    //   Lr_eq = Lr1 + Lr2/n²              = 2·Lr1     (sym)
    //   Cr_eq = Cr1·(n²·Cr2)/(Cr1+n²·Cr2) = Cr1/2     (sym)
    //   Lm    unchanged
    // and ω_P = 1/√(Lr_eq·Cr_eq) = 1/√(Lr1·Cr1) ≡ ω_r exactly. This means
    // the symmetric CLLC and the equivalent LLC produce identical state
    // trajectories on the (i_Lr_eq, i_Lm, v_Cr_eq) phase space; the
    // physical primary current is i_Lr_eq, and the primary cap voltage is
    // v_Cr_eq · (Cr1+n²·Cr2)/(n²·Cr2) = 2·v_Cr_eq for symmetric.
    //
    // For asymmetric tanks (a≠1 or b≠1, plan §7) the collapsed form is no
    // longer exact. The asymmetric 5-state TDA is deferred to P7.
    //
    // The solver mirrors the LLC TDA in Llc.cpp (same author, same shapes):
    //   3-state vector  x = (i_Lr_eq, i_Lm, v_Cr_eq)   [ "iLs ", "iL", "vC" ]
    //   4 sub-states    P_POS, P_NEG  (power delivery, ±n·Vout clamp)
    //                   F             (freewheel, Lm in series; sign-symm.)
    //   2D Newton on    F(x0) = x(Ts/2) + x(0)        (half-wave antisym)
    //   Closed-form propagator per sub-state — never Euler.
    // =========================================================================

    namespace {  // file-static TDA solver internals

    // Sub-states observed during ONE positive half cycle of the primary
    // bridge (V_bridge = +Vi for the entire half). For the next half cycle
    // we exploit half-wave antisymmetry. Mirrors LLC convention precisely
    // (single freewheel sub-state F; P_POS / P_NEG distinguish the
    // direction of the secondary-rectifier conduction).
    enum class CllcSubState { P_POS, P_NEG, F };

    struct CllcStateVector {
        double iLs;   ///< collapsed series-inductor current (= primary i_Lr1 at sym)
        double iL;    ///< magnetizing current i_Lm
        double vC;    ///< collapsed series-cap voltage
    };

    struct CllcSubStateSegment {
        CllcSubState state;
        double t_start;
        double t_end;
        CllcStateVector x_start;
        CllcStateVector x_end;
    };

    /// Sub-state-specific natural frequency ω and characteristic Z =√(L/C).
    void cllc_substate_freq(CllcSubState s, double Lr_eq, double Lm, double Cr_eq,
                             double& w, double& Z) {
        if (s == CllcSubState::F) {
            // Freewheel: Lm in series with Lr_eq, secondary off.
            double L_eff = Lr_eq + Lm;
            w = 1.0 / std::sqrt(L_eff * Cr_eq);
            Z = std::sqrt(L_eff / Cr_eq);
        } else {
            // Power delivery: Lm clamped to ±Vo, only Lr_eq+Cr_eq oscillate.
            w = 1.0 / std::sqrt(Lr_eq * Cr_eq);
            Z = std::sqrt(Lr_eq / Cr_eq);
        }
    }

    /// Closed-form propagator for one sub-state. Computes x(t_start+dt)
    /// EXACTLY given x(t_start) — never numerical integration.
    /// Bridge voltage during this half cycle is +Vi (positive half).
    CllcStateVector cllc_propagate_substate(CllcSubState s, CllcStateVector x_in,
                                             double dt, double Vi, double Vo,
                                             double Lr_eq, double Lm, double Cr_eq) {
        CllcStateVector out{};
        if (dt <= 0) return x_in;
        double w, Z;
        cllc_substate_freq(s, Lr_eq, Lm, Cr_eq, w, Z);
        double cs = std::cos(w * dt);
        double sn = std::sin(w * dt);

        if (s == CllcSubState::F) {
            // Freewheel: i_Lr ≡ i_Lm. Drive is the bridge voltage +Vi.
            double V_eq = Vi;
            double dV = x_in.vC - V_eq;
            double iLs_new = x_in.iLs * cs - (dV / Z) * sn;
            double vC_new  = V_eq + dV * cs + x_in.iLs * Z * sn;
            out.iLs = iLs_new;
            out.iL  = iLs_new;
            out.vC  = vC_new;
            return out;
        }

        // Power-delivery sub-states (bridge = +Vi throughout):
        //   P_POS (Id = iLs - iL > 0): V_drive = Vi - Vo, dIL/dt = +Vo/Lm
        //   P_NEG (Id < 0):            V_drive = Vi + Vo, dIL/dt = -Vo/Lm
        double V_drive = (s == CllcSubState::P_POS) ? (Vi - Vo) : (Vi + Vo);
        double dIL_dt  = (s == CllcSubState::P_POS) ? (+Vo / Lm) : (-Vo / Lm);

        double dV = x_in.vC - V_drive;
        out.iLs = x_in.iLs * cs - (dV / Z) * sn;
        out.vC  = V_drive + dV * cs + x_in.iLs * Z * sn;
        out.iL  = x_in.iL + dIL_dt * dt;
        return out;
    }

    /// Trigger value at the END of a propagation step.
    ///   P_POS → F : g = iL - iLs (rises to 0 as Id = iLs - iL falls to 0)
    ///   P_NEG → F : g = iLs - iL (rises to 0 as Id rises to 0)
    ///   F → P    : VLm exits the ±Vo clamp window
    double cllc_trigger_value(CllcSubState s, CllcStateVector x_in, double dt,
                               double Vi, double Vo, double Lr_eq, double Lm, double Cr_eq) {
        CllcStateVector x = cllc_propagate_substate(s, x_in, dt, Vi, Vo, Lr_eq, Lm, Cr_eq);
        if (s == CllcSubState::P_POS) return x.iL  - x.iLs;
        if (s == CllcSubState::P_NEG) return x.iLs - x.iL;
        // F: VLm exits ±Vo window.  VLm = (Lm/(Lr_eq+Lm))·(Vi - vC)
        double VLm = (Lm / (Lr_eq + Lm)) * (Vi - x.vC);
        double pos_violation =  VLm - Vo;
        double neg_violation = -VLm - Vo;
        return std::max(pos_violation, neg_violation);
    }

    /// Coarse-grid bisection event finder (mirrors LLC pattern).
    double cllc_find_next_event(CllcSubState s, CllcStateVector x_in, double t_max,
                                 double Vi, double Vo, double Lr_eq, double Lm, double Cr_eq) {
        constexpr int COARSE_STEPS = 64;
        double dt_coarse = t_max / COARSE_STEPS;
        double prev_g = cllc_trigger_value(s, x_in, 1e-12, Vi, Vo, Lr_eq, Lm, Cr_eq);
        if (prev_g >= 0) return 0.0;
        for (int k = 1; k <= COARSE_STEPS; ++k) {
            double t = k * dt_coarse;
            double g = cllc_trigger_value(s, x_in, t, Vi, Vo, Lr_eq, Lm, Cr_eq);
            if (g >= 0 && std::isfinite(g)) {
                double lo = t - dt_coarse, hi = t, g_lo = prev_g;
                for (int it = 0; it < 50; ++it) {
                    double mid = 0.5 * (lo + hi);
                    double g_mid = cllc_trigger_value(s, x_in, mid, Vi, Vo, Lr_eq, Lm, Cr_eq);
                    if (g_mid * g_lo < 0) { hi = mid; }
                    else { lo = mid; g_lo = g_mid; }
                    if ((hi - lo) < 1e-12) break;
                }
                return 0.5 * (lo + hi);
            }
            prev_g = g;
        }
        return t_max;
    }

    /// Determine the next P sub-state after a freewheel ends (sign of VLm).
    CllcSubState cllc_next_state_after_F(CllcStateVector x_at_event,
                                          double Vi, double Lr_eq, double Lm) {
        double VLm = (Lm / (Lr_eq + Lm)) * (Vi - x_at_event.vC);
        return (VLm > 0) ? CllcSubState::P_POS : CllcSubState::P_NEG;
    }

    /// Initial sub-state at t=0+ of the half cycle (bridge just switched to +Vi).
    CllcSubState cllc_initial_substate(CllcStateVector x0, double Vi, double Vo,
                                        double Lr_eq, double Lm) {
        double Id = x0.iLs - x0.iL;
        if (std::abs(Id) > 1e-9) {
            return (Id > 0) ? CllcSubState::P_POS : CllcSubState::P_NEG;
        }
        double VLm = (Lm / (Lr_eq + Lm)) * (Vi - x0.vC);
        if (VLm >  Vo) return CllcSubState::P_POS;
        if (VLm < -Vo) return CllcSubState::P_NEG;
        return CllcSubState::F;
    }

    /// Drive the event loop over [0, Thalf] for one half cycle.
    std::vector<CllcSubStateSegment> cllc_propagate_half_cycle(
        CllcStateVector x0, double Thalf,
        double Vi, double Vo, double Lr_eq, double Lm, double Cr_eq)
    {
        std::vector<CllcSubStateSegment> segments;
        segments.reserve(8);
        CllcSubState current = cllc_initial_substate(x0, Vi, Vo, Lr_eq, Lm);
        CllcStateVector x = x0;
        double t = 0.0;
        constexpr int MAX_SEGMENTS = 16;
        for (int k = 0; k < MAX_SEGMENTS; ++k) {
            double remaining = Thalf - t;
            if (remaining <= 1e-15) break;
            double t_event = cllc_find_next_event(current, x, remaining,
                                                   Vi, Vo, Lr_eq, Lm, Cr_eq);
            double dt = std::min(t_event, remaining);
            if (dt < 1e-15 && k > 0) {
                CllcSubState next = (current == CllcSubState::F)
                    ? cllc_next_state_after_F(x, Vi, Lr_eq, Lm)
                    : CllcSubState::F;
                current = next;
                continue;
            }
            CllcStateVector x_end = cllc_propagate_substate(current, x, dt,
                                                              Vi, Vo, Lr_eq, Lm, Cr_eq);
            segments.push_back({current, t, t + dt, x, x_end});
            t += dt; x = x_end;
            if (t >= Thalf - 1e-15) break;
            if (current == CllcSubState::F) {
                current = cllc_next_state_after_F(x, Vi, Lr_eq, Lm);
            } else {
                current = CllcSubState::F;
            }
        }
        return segments;
    }

    /// 2D damped-Newton on F(x0) = propagate_half(x0).end + x0  (antisymmetry).
    /// Falls back to a Picard step on stagnation.
    CllcStateVector cllc_solve_steady_state(
        CllcStateVector x0_seed, double Thalf,
        double Vi, double Vo, double Lr_eq, double Lm, double Cr_eq,
        std::vector<CllcSubStateSegment>& outSegments,
        double& outResidual)
    {
        auto eval_F = [&](CllcStateVector x0) -> std::array<double, 3> {
            auto segs = cllc_propagate_half_cycle(x0, Thalf, Vi, Vo, Lr_eq, Lm, Cr_eq);
            CllcStateVector xe = segs.empty() ? x0 : segs.back().x_end;
            return {xe.iLs + x0.iLs, xe.iL + x0.iL, xe.vC + x0.vC};
        };
        auto norm = [](const std::array<double, 3>& f) {
            return std::sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
        };

        CllcStateVector x0 = x0_seed;
        auto F = eval_F(x0);
        double r = norm(F);
        constexpr int MAX_ITERS = 25;
        constexpr double TOL = 1e-7;
        double scale_i = std::max(1e-3, 0.01 * std::abs(x0.iLs) + 1e-3);
        double scale_v = std::max(1e-2, 0.01 * std::abs(x0.vC)  + 1e-2);
        double damping = 1.0;
        double prev_r = r;
        int stagnant = 0;

        for (int iter = 0; iter < MAX_ITERS && r > TOL; ++iter) {
            double J[3][3];
            CllcStateVector xp, xm;
            std::array<double, 3> Fp, Fm;
            xp = x0; xp.iLs += scale_i; xm = x0; xm.iLs -= scale_i;
            Fp = eval_F(xp); Fm = eval_F(xm);
            for (int i = 0; i < 3; ++i) J[i][0] = (Fp[i] - Fm[i]) / (2 * scale_i);
            xp = x0; xp.iL  += scale_i; xm = x0; xm.iL  -= scale_i;
            Fp = eval_F(xp); Fm = eval_F(xm);
            for (int i = 0; i < 3; ++i) J[i][1] = (Fp[i] - Fm[i]) / (2 * scale_i);
            xp = x0; xp.vC  += scale_v; xm = x0; xm.vC  -= scale_v;
            Fp = eval_F(xp); Fm = eval_F(xm);
            for (int i = 0; i < 3; ++i) J[i][2] = (Fp[i] - Fm[i]) / (2 * scale_v);

            double A[3][4];
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) A[i][j] = J[i][j];
                A[i][3] = -F[i];
            }
            bool singular = false;
            for (int col = 0; col < 3; ++col) {
                int piv = col; double maxv = std::abs(A[col][col]);
                for (int row = col + 1; row < 3; ++row) {
                    if (std::abs(A[row][col]) > maxv) { maxv = std::abs(A[row][col]); piv = row; }
                }
                if (maxv < 1e-14) { singular = true; break; }
                if (piv != col) std::swap(A[col], A[piv]);
                for (int row = col + 1; row < 3; ++row) {
                    double f = A[row][col] / A[col][col];
                    for (int k = col; k < 4; ++k) A[row][k] -= f * A[col][k];
                }
            }
            if (singular) break;
            double dx[3];
            for (int i = 2; i >= 0; --i) {
                double sum = A[i][3];
                for (int j = i + 1; j < 3; ++j) sum -= A[i][j] * dx[j];
                dx[i] = sum / A[i][i];
            }

            double try_d = damping;
            CllcStateVector x_new{}; std::array<double, 3> F_new{};
            double r_new = r; bool accepted = false;
            for (int ls = 0; ls < 6; ++ls) {
                x_new.iLs = x0.iLs + try_d * dx[0];
                x_new.iL  = x0.iL  + try_d * dx[1];
                x_new.vC  = x0.vC  + try_d * dx[2];
                F_new = eval_F(x_new);
                r_new = norm(F_new);
                if (std::isfinite(r_new) && r_new < r) { accepted = true; break; }
                try_d *= 0.5;
            }
            if (!accepted) {
                auto segs = cllc_propagate_half_cycle(x0, Thalf, Vi, Vo, Lr_eq, Lm, Cr_eq);
                if (!segs.empty()) {
                    CllcStateVector xe = segs.back().x_end;
                    x_new.iLs = -xe.iLs; x_new.iL = -xe.iL; x_new.vC = -xe.vC;
                    F_new = eval_F(x_new);
                    r_new = norm(F_new);
                }
            }
            x0 = x_new; F = F_new;
            if (r_new >= prev_r * 0.999) { stagnant++; damping *= 0.7; }
            else { stagnant = 0; damping = std::min(1.0, damping * 1.2); }
            prev_r = r_new; r = r_new;
            if (stagnant >= 4) break;
        }
        outSegments = cllc_propagate_half_cycle(x0, Thalf, Vi, Vo, Lr_eq, Lm, Cr_eq);
        outResidual = r;
        return x0;
    }

    /// Classify the segment chain into one of the three CLLC modes:
    ///   1 = sub-resonant boost  (chain CONTAINS a freewheel segment F:
    ///                            resonant half-sine completed before the
    ///                            bridge edge, leaving a freewheel tail)
    ///   2 = at resonance        (chain is all-P; Thalf ≈ natural P-state
    ///                            half-period π·√(Lr_eq·Cr_eq); P-segment
    ///                            naturally completes at the bridge edge)
    ///   3 = super-resonant buck (chain is all-P; Thalf < natural P-state
    ///                            half-period; bridge edge clips P-segment
    ///                            before its natural zero-crossing)
    /// The duration-vs-natural-half-period criterion is robust at high
    /// gain demand where the current amplitude shrinks (which would make
    /// the |Id_end|/|iL_start| ratio test misleading).
    /// Returns 0 if unrecognised.
    int cllc_classify_mode(const std::vector<CllcSubStateSegment>& segs,
                           double Thalf, double Lr_eq, double Cr_eq) {
        if (segs.empty()) return 0;
        // Any freewheel segment in the chain indicates that the rectifier
        // current naturally returned to zero before the next bridge edge —
        // characteristic of sub-resonant boost operation (mode 1).
        for (const auto& s : segs) {
            if (s.state == CllcSubState::F) return 1;
        }
        // Chain is all P-segments. Compare Thalf to the natural P-state
        // half-period Tp_half = π·√(Lr_eq·Cr_eq). At fsw = fr_eq this ratio
        // is exactly 1.0 (mode 2). At fsw > fr_eq (super-resonant buck) the
        // bridge clips the half-sine and the ratio drops below 1.
        double Tp_half = M_PI * std::sqrt(Lr_eq * Cr_eq);
        if (Tp_half <= 0) return 0;
        double ratio = Thalf / Tp_half;
        if (ratio < 0.9) return 3;   // bridge clipped → super-resonant
        return 2;                    // at resonance (or very close to it)
    }

    /// Convert a segment chain into uniformly-sampled (i_Lr, i_Lm, v_Cr)
    /// vectors of length N+1 over [0, Thalf].
    void cllc_sample_segments(const std::vector<CllcSubStateSegment>& segs,
                               double Thalf, int N,
                               double Vi, double Vo,
                               double Lr_eq, double Lm, double Cr_eq,
                               std::vector<double>& iLs_out,
                               std::vector<double>& iL_out,
                               std::vector<double>& vC_out)
    {
        double dt = Thalf / N;
        size_t segIdx = 0;
        for (int k = 0; k <= N; ++k) {
            double t = k * dt;
            if (t > Thalf) t = Thalf;
            while (segIdx + 1 < segs.size() && t > segs[segIdx].t_end + 1e-15) ++segIdx;
            if (segs.empty()) {
                iLs_out[k] = iL_out[k] = vC_out[k] = 0.0;
                continue;
            }
            const auto& seg = segs[segIdx];
            double t_local = t - seg.t_start;
            if (t_local < 0) t_local = 0;
            double seg_dt = seg.t_end - seg.t_start;
            if (t_local > seg_dt) t_local = seg_dt;
            CllcStateVector x = cllc_propagate_substate(seg.state, seg.x_start,
                                                          t_local, Vi, Vo,
                                                          Lr_eq, Lm, Cr_eq);
            iLs_out[k] = x.iLs;
            iL_out[k]  = x.iL;
            vC_out[k]  = x.vC;
        }
    }

    } // anonymous namespace

    // -------------------------------------------------------------------------

    OperatingPoint CllcConverter::process_operating_point_for_input_voltage(
        double inputVoltage,
        const CllcOperatingPoint& cllcOpPoint,
        double turnsRatio,
        double magnetizingInductance,
        const CllcResonantParameters& params) {

        OperatingPoint operatingPoint;

        if (cllcOpPoint.get_output_voltages().empty() || cllcOpPoint.get_output_currents().empty())
            throw std::runtime_error("CLLC: operating point missing output voltages/currents");

        double switchingFrequency = get_value_or(cllcOpPoint.get_switching_frequency(), 0.0);
        if (switchingFrequency <= 0)
            throw std::runtime_error("CLLC: operating point has invalid switching frequency");

        double outputVoltage = cllcOpPoint.get_output_voltages()[0];
        double outputCurrent = cllcOpPoint.get_output_currents()[0];
        double n  = turnsRatio;
        double Lm = magnetizingInductance;
        if (n <= 0 || Lm <= 0)
            throw std::runtime_error("CLLC: invalid turns ratio or magnetizing inductance");

        // -------------------------------------------------------------------
        // Collapse 5-element CLLC tank to 3-state LLC-equivalent (Lr_eq, Cr_eq, Lm)
        // referred to the primary. Lr2 is stored on the secondary side; refer
        // to primary by ×n². Cr2 is stored on the secondary side; refer to
        // primary by ÷n². Series combinations:
        //   Lr_eq = Lr1 + n²·Lr2_sec
        //   Cr_eq = (Cr1 · Cr2_pri) / (Cr1 + Cr2_pri)   with Cr2_pri = Cr2_sec/n²
        // For symmetric (a=b=1):  Lr_eq = 2·Lr1,  Cr_eq = Cr1/2.
        // -------------------------------------------------------------------
        double Lr1 = params.primaryResonantInductance;
        double Cr1 = params.primaryResonantCapacitance;
        double Lr2_sec = params.secondaryResonantInductance;
        double Cr2_sec = params.secondaryResonantCapacitance;
        if (Lr1 <= 0 || Cr1 <= 0 || Lr2_sec <= 0 || Cr2_sec <= 0)
            throw std::runtime_error("CLLC: resonant tank values invalid (Lr1/Cr1/Lr2/Cr2 must be > 0)");

        double Lr2_pri = Lr2_sec * n * n;
        double Cr2_pri = Cr2_sec / (n * n);
        double Lr_eq  = Lr1 + Lr2_pri;
        double Cr_eq  = (Cr1 * Cr2_pri) / (Cr1 + Cr2_pri);

        double k_bridge = get_bridge_voltage_factor();

        // P8b — REVERSE mode: per plan §5.5, the same TDA propagator runs
        // with Vi ↔ n·Vout swapped. The active inverter is now on the
        // secondary; reflected to primary terms its drive is k_bridge·n·Vout
        // (a transformer step-up of the secondary EMF). The rectifier is on
        // the primary; its load voltage IS the primary rail (Vin), already
        // in primary-referred units (no n multiplication). Tank state-space
        // matrices (Lr_eq, Cr_eq, Lm) stay primary-referred; the solver
        // returns iLs as the physical primary winding tank current — same
        // signal SPICE measures via Vpri_sense — so reverse PtP NRMSE
        // compares apples to apples without further conversion.
        bool isReverse = (cllcOpPoint.get_power_flow() == CllcPowerFlow::REVERSE);
        double Vi, Vo;
        if (isReverse) {
            Vi = k_bridge * n * outputVoltage;   // sec inverter, pri-referred
            Vo = inputVoltage;                    // pri rectifier load
        }
        else {
            Vi = k_bridge * inputVoltage;
            Vo = n * outputVoltage;               // reflected output (pri-referred)
        }

        double period = 1.0 / switchingFrequency;
        double Thalf  = period / 2.0;
        // Dead time is not represented by the TDA model (would require a 5th
        // sub-state with both bridge legs floating). Mirrors LLC convention.
        double Thalf_eff = Thalf;

        // -------------------------------------------------------------------
        // Seed the Newton solver with FHA-style estimates. The DC load
        // current referred to the primary tank differs by direction:
        //   FORWARD: load on secondary, I_load_pri = outputCurrent / n
        //   REVERSE: load on primary  , I_load_pri = (Vo·Io)/Vin
        //                              ≈ outputCurrent·outputVoltage/inputVoltage
        // -------------------------------------------------------------------
        double Im_pk_est = Vo * Thalf_eff / (2.0 * Lm);
        double Iload_reflected;
        if (isReverse) {
            double powerThru = outputVoltage * outputCurrent;
            Iload_reflected = powerThru / std::max(std::abs(inputVoltage), 1e-9);
        }
        else {
            Iload_reflected = outputCurrent / n;
        }
        double Ires_est = std::max(std::abs(Im_pk_est) + std::abs(Iload_reflected),
                                   std::abs(Iload_reflected) * 1.5);
        CllcStateVector x0_seed{ -Ires_est, -Im_pk_est, 0.0 };

        // LIP perturbation (mirrors LLC) to avoid the singular-Jacobian ridge
        // when Vi ≈ Vo and fsw ≈ fr.
        double Vi_solver = Vi;
        double denom_vo = std::max(std::abs(Vi), std::abs(Vo));
        if (denom_vo > 0 && std::abs(Vi - Vo) / denom_vo < 0.005) {
            Vi_solver = Vi * 1.005;
        }

        std::vector<CllcSubStateSegment> segments;
        double residual = 0.0;
        CllcStateVector x0 = cllc_solve_steady_state(x0_seed, Thalf_eff,
                                                       Vi_solver, Vo, Lr_eq, Lm, Cr_eq,
                                                       segments, residual);

        // FIXME (asymmetric 5-state TDA): per plan §7 and the comment at
        // line 357, asymmetric tanks (a≠1 or b≠1) are currently run through
        // the SYMMETRIC collapsed 3-state form (Lr_eq, Cr_eq, Lm). On the
        // KIT-20kW-Asymmetric reference (a=0.95, b=1.052) the Newton over
        // this collapsed form has multiple local minima:
        //   • OP[0] (Vi=804 nominal+LIP): converges to vC≈-4545 V (5.7×Vi)
        //     with residual ≈ 8.25; the cap discharge through the tank
        //     impedance Zr=√(Lr_eq/Cr_eq)≈15Ω yields a spurious peak
        //     current ≈ 4545/15 ≈ 300 A, matching the observed analytical
        //     i_range = ±301 A vs SPICE ±54 A.
        //   • OP[2] (Vi=900): false convergence (residual ≈ 1.7e-4) to a
        //     non-physical iLs=5.4e4 A, vC=-1.6 MV — caught by the sanity
        //     block below; the fallback seed is then re-propagated but is
        //     itself not a steady state, so the emitted waveform contains
        //     a transient.
        // The shape-only ptp_nrmse (PtpHelpers.h:213) normalises each
        // waveform by its own AC-RMS, so this 5.6× amplitude bug still
        // passes the 16 % NRMSE gate. The proper fix is implementing the
        // asymmetric 5-state TDA (4-D Newton on (iLr1, iLm, vCr1, vCr2)
        // with separate Lr1/Cr1 and Lr2/Cr2 sub-state ODEs). Sanity tweaks
        // here only shift the failure mode — fallback is non-steady-state.
        // Sanity-check against null-space blow-up; fall back to seed if violated.
        double sanity_iLs = std::max(10.0 * Ires_est, 20.0);
        double sanity_vC  = std::max(10.0 * std::abs(Vi), 10.0 * std::abs(Vo));
        if (sanity_vC < 200.0) sanity_vC = 200.0;
        if (!std::isfinite(x0.iLs) || !std::isfinite(x0.iL) || !std::isfinite(x0.vC) ||
            std::abs(x0.iLs) > sanity_iLs ||
            std::abs(x0.vC)  > sanity_vC) {
            x0 = x0_seed;
            residual = -1.0;
        }
        // Re-propagate with the authoritative Vi for waveform emission.
        segments = cllc_propagate_half_cycle(x0, Thalf_eff, Vi, Vo, Lr_eq, Lm, Cr_eq);

        // -------------------------------------------------------------------
        // Sample the segment chain across the positive half cycle.
        // -------------------------------------------------------------------
        const int N = 200;
        double dt = Thalf_eff / N;
        std::vector<double> ILs_pos(N + 1, 0.0), IL_pos(N + 1, 0.0), Vc_pos(N + 1, 0.0);
        cllc_sample_segments(segments, Thalf_eff, N,
                              Vi, Vo, Lr_eq, Lm, Cr_eq,
                              ILs_pos, IL_pos, Vc_pos);

        // VLm at each sample (closed-form per sub-state).
        std::vector<double> VLm_pos(N + 1, 0.0);
        size_t segIdx = 0;
        for (int k = 0; k <= N; ++k) {
            double t = std::min<double>(k * dt, Thalf_eff);
            while (segIdx + 1 < segments.size() && t > segments[segIdx].t_end + 1e-15) ++segIdx;
            if (segments.empty()) { VLm_pos[k] = 0.0; continue; }
            const auto& seg = segments[segIdx];
            switch (seg.state) {
                case CllcSubState::P_POS: VLm_pos[k] = +Vo; break;
                case CllcSubState::P_NEG: VLm_pos[k] = -Vo; break;
                case CllcSubState::F: VLm_pos[k] = (Lm/(Lr_eq+Lm)) * (Vi - Vc_pos[k]); break;
            }
        }

        // -------------------------------------------------------------------
        // Diagnostics (per CLLC_REWRITE_PLAN §5.4)
        // -------------------------------------------------------------------
        lastMode = cllc_classify_mode(segments, Thalf_eff, Lr_eq, Cr_eq);
        lastSubStateSequence.clear();
        for (const auto& seg : segments)
            lastSubStateSequence.push_back(static_cast<int>(seg.state));
        lastSteadyStateResidual = residual;
        // LIP frequency: 1/(2π·√(Lr_eq·Cr_eq)). Match Llc convention.
        computedLipFrequency = 1.0 / (2.0 * M_PI * std::sqrt(Lr_eq * Cr_eq));

        // Peak primary current and peak resonant-cap voltage (across both Cr1
        // and Cr2). Cr1 sees a fraction Cr2_pri/(Cr1+Cr2_pri) of the collapsed
        // vC; Cr2_pri sees Cr1/(Cr1+Cr2_pri). Voltage divider — capacitor-cap
        // voltage divider for series caps: V_each = vC_total · C_other/(C1+C2).
        double Ip_peak = 0.0, Vc_peak = 0.0;
        for (int k = 0; k <= N; ++k) {
            Ip_peak = std::max(Ip_peak, std::abs(ILs_pos[k]));
            Vc_peak = std::max(Vc_peak, std::abs(Vc_pos[k]));
        }
        lastPrimaryPeakCurrent = Ip_peak;
        // Cr1 peak voltage (referred to primary, since Cr1 is on primary):
        // V_Cr1 = vC · Cr2_pri/(Cr1 + Cr2_pri)
        double V_Cr1_peak = Vc_peak * Cr2_pri / (Cr1 + Cr2_pri);
        lastResonantCapPeakVoltage = V_Cr1_peak;

        // -------------------------------------------------------------------
        // ZVS diagnostics (Infineon AN-2024-06 Eq 28; CLLC_REWRITE_PLAN §5.4).
        //
        // At the bridge edge (t = Thalf_eff) the magnetizing current iLm must
        // be large enough to discharge the bridge-leg output capacitance
        // Coss_eq within the dead time td. The ZVS charge-balance criterion
        // for a full-bridge primary leg (two FETs share the swing through Vi
        // and the leg-pair reservoir is 2·Coss) is
        //
        //     |iLm(t_switch)| · td  ≥  2·Coss·Vi
        //   ⇒ iLm_threshold       =  2·Coss·Vi / td
        //   ⇒ ZVS margin          =  |iLm(t_switch)| − iLm_threshold
        //
        // For half-bridge the swing only crosses one Coss but uses the same
        // reservoir-charge formula scaled by k_bridge (= 0.5 HB / 1.0 FB).
        // The resonant transition time (slew of the bridge midpoint) is
        //
        //     t_resonant = 2·Coss·Vi / |iLm(t_switch)|
        //
        // and must satisfy t_resonant ≤ td for ZVS.
        //
        // Forward direction populates lastZvsMarginPrimary; reverse direction
        // populates lastZvsMarginSecondary using the secondary-referred
        // current iLs/n at the switch edge and the input voltage as seen
        // from the secondary winding (n·Vo for reverse). Reverse-mode
        // physical SR drive is P8; here we record the analytical figure
        // either way for the active side.
        // -------------------------------------------------------------------
        {
            double Coss = mosfetCoss;
            double td   = deadTime;
            // The plan's CllcConverter::process_operating_point_for_input_voltage
            // throws if td ≤ 0 in the SPICE path; here we guard the divisor
            // explicitly so a missing dead-time still yields an honest "no
            // ZVS data" reading rather than NaN/Inf in the diagnostic.
            if (td <= 0 || Coss <= 0) {
                lastZvsMarginPrimary = std::numeric_limits<double>::quiet_NaN();
                lastZvsMarginSecondary = std::numeric_limits<double>::quiet_NaN();
                lastResonantTransitionTime = std::numeric_limits<double>::quiet_NaN();
            }
            else {
                double iLm_at_switch = std::abs(IL_pos[N]);
                double Vi_for_zvs = std::abs(Vi);  // Vi already includes k_bridge
                double iLm_threshold = (Vi_for_zvs > 0)
                    ? (2.0 * Coss * Vi_for_zvs / td) : 0.0;
                double zvs_margin = iLm_at_switch - iLm_threshold;
                double t_resonant = (iLm_at_switch > 0)
                    ? (2.0 * Coss * Vi_for_zvs / iLm_at_switch)
                    : std::numeric_limits<double>::infinity();

                if (isReverse) {
                    // Active bridge in reverse mode is on the secondary side.
                    // Secondary-referred magnetizing current is iLm/n; the
                    // secondary bridge sees ±k_bridge·outputVoltage. After
                    // the P8b Vi/Vo swap, Vi is already k_bridge·n·Vout in
                    // primary terms, so Vsec_actual = |Vi|/n.
                    double iLm_sec = iLm_at_switch / n;
                    double Vsec    = std::abs(Vi) / n;
                    double iLm_thr_sec = (Vsec > 0)
                        ? (2.0 * Coss * Vsec / td) : 0.0;
                    lastZvsMarginSecondary = iLm_sec - iLm_thr_sec;
                    lastZvsMarginPrimary   = 0.0;
                }
                else {
                    lastZvsMarginPrimary   = zvs_margin;
                    lastZvsMarginSecondary = 0.0;
                }
                lastResonantTransitionTime = t_resonant;
            }
        }

        // -------------------------------------------------------------------
        // Build full-period waveforms by half-wave antisymmetry.
        // -------------------------------------------------------------------
        int totalSamples = 2 * N + 1;
        std::vector<double> time_full(totalSamples);
        std::vector<double> ILs_full(totalSamples);
        std::vector<double> IL_full(totalSamples);
        std::vector<double> VLm_full(totalSamples);
        std::vector<double> Vc_full(totalSamples);

        for (int k = 0; k <= N; ++k) {
            time_full[k] = k * dt;
            ILs_full[k] = std::isfinite(ILs_pos[k]) ? ILs_pos[k] : 0.0;
            IL_full[k]  = std::isfinite(IL_pos[k])  ? IL_pos[k]  : 0.0;
            VLm_full[k] = std::isfinite(VLm_pos[k]) ? VLm_pos[k] : 0.0;
            Vc_full[k]  = std::isfinite(Vc_pos[k])  ? Vc_pos[k]  : 0.0;
        }
        for (int k = 1; k <= N; ++k) {
            time_full[N + k] = Thalf_eff + k * dt;
            ILs_full[N + k] = -ILs_full[k];
            IL_full[N + k]  = -IL_full[k];
            VLm_full[N + k] = -VLm_full[k];
            Vc_full[N + k]  = -Vc_full[k];
        }

        // -------------------------------------------------------------------
        // Extra-component waveforms (Cr1, Cr2, Lr1, Lr2) for round-trip via
        // get_extra_components_inputs in P5. Voltages split via series cap
        // divider; currents are the same series tank current ILs.
        // -------------------------------------------------------------------
        {
            std::vector<double> Vcr1_full(totalSamples), Vcr2_full(totalSamples);
            // Capacitor voltage divider: vC_total = V_Cr1 + V_Cr2_pri
            //   V_Cr1   = vC * Cr2_pri/(Cr1 + Cr2_pri)
            //   V_Cr2_pri = vC * Cr1/(Cr1 + Cr2_pri)
            // V_Cr2 (secondary side) = V_Cr2_pri * n
            double k1 = Cr2_pri / (Cr1 + Cr2_pri);
            double k2 = Cr1     / (Cr1 + Cr2_pri);
            for (int k = 0; k < totalSamples; ++k) {
                Vcr1_full[k] = Vc_full[k] * k1;
                Vcr2_full[k] = Vc_full[k] * k2 * n;   // referred back to secondary
            }

            Waveform t_wf, vCr1_wf, vCr2_wf, iCr_wf;
            iCr_wf.set_data(ILs_full); iCr_wf.set_time(time_full);
            iCr_wf.set_ancillary_label(WaveformLabel::CUSTOM);
            vCr1_wf.set_data(Vcr1_full); vCr1_wf.set_time(time_full);
            vCr1_wf.set_ancillary_label(WaveformLabel::CUSTOM);
            vCr2_wf.set_data(Vcr2_full); vCr2_wf.set_time(time_full);
            vCr2_wf.set_ancillary_label(WaveformLabel::CUSTOM);

            extraCr1VoltageWaveforms.push_back(vCr1_wf);
            extraCr1CurrentWaveforms.push_back(iCr_wf);
            extraCr2VoltageWaveforms.push_back(vCr2_wf);
            // Cr2 on secondary side carries n×ILs (ampere-turn balance).
            std::vector<double> ICr2_sec(totalSamples);
            for (int k = 0; k < totalSamples; ++k) ICr2_sec[k] = n * (ILs_full[k] - IL_full[k]);
            Waveform iCr2_wf;
            iCr2_wf.set_data(ICr2_sec); iCr2_wf.set_time(time_full);
            iCr2_wf.set_ancillary_label(WaveformLabel::CUSTOM);
            extraCr2CurrentWaveforms.push_back(iCr2_wf);

            // Lr1 voltage = Vi_sq − V_Cr1 − VLm  (KVL on primary loop in P-mode)
            // Lr2 voltage (secondary side) = (Vi_sq − VLm)·... actually for full
            // generality V_Lr1 = (Lr1/Lr_eq)·(Vi_sq − vC − VLm). Likewise for Lr2.
            std::vector<double> VLr1_full(totalSamples), VLr2_full(totalSamples);
            for (int k = 0; k < totalSamples; ++k) {
                double Vi_sq = (k * dt < Thalf_eff) ? Vi : -Vi;  // bipolar bridge
                if (k > N) Vi_sq = -Vi;
                double V_loop = Vi_sq - Vc_full[k] - VLm_full[k];
                VLr1_full[k] = (Lr1     / Lr_eq) * V_loop;
                VLr2_full[k] = (Lr2_pri / Lr_eq) * V_loop / n;   // back to secondary side
            }
            Waveform vLr1_wf, vLr2_wf;
            vLr1_wf.set_data(VLr1_full); vLr1_wf.set_time(time_full);
            vLr1_wf.set_ancillary_label(WaveformLabel::CUSTOM);
            vLr2_wf.set_data(VLr2_full); vLr2_wf.set_time(time_full);
            vLr2_wf.set_ancillary_label(WaveformLabel::CUSTOM);
            extraLr1VoltageWaveforms.push_back(vLr1_wf);
            extraLr1CurrentWaveforms.push_back(iCr_wf);   // same as ILs
            extraLr2VoltageWaveforms.push_back(vLr2_wf);
            extraLr2CurrentWaveforms.push_back(iCr2_wf);  // same as ICr2_sec

            extraTimeVectors.push_back(time_full);
        }

        // -------------------------------------------------------------------
        // PRIMARY WINDING excitation: current = ILs (series tank current),
        // voltage = VLm (primary magnetizing voltage; this is what the core
        // sees, *not* the bridge voltage).
        // -------------------------------------------------------------------
        {
            Waveform currentWaveform;
            currentWaveform.set_data(ILs_full);
            currentWaveform.set_time(time_full);
            currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);

            Waveform voltageWaveform;
            voltageWaveform.set_data(VLm_full);
            voltageWaveform.set_time(time_full);
            voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);

            auto excitation = complete_excitation(currentWaveform, voltageWaveform,
                                                   switchingFrequency, "Primary");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // -------------------------------------------------------------------
        // SECONDARY WINDING excitation (single full-wave secondary, not
        // center-tapped). I_sec = n·(ILs − IL), V_sec = VLm/n.
        // -------------------------------------------------------------------
        {
            std::vector<double> iSec(totalSamples), vSec(totalSamples);
            for (int k = 0; k < totalSamples; ++k) {
                iSec[k] = n * (ILs_full[k] - IL_full[k]);
                vSec[k] = VLm_full[k] / n;
            }
            Waveform currentWaveform;
            currentWaveform.set_data(iSec);
            currentWaveform.set_time(time_full);
            currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);

            Waveform voltageWaveform;
            voltageWaveform.set_data(vSec);
            voltageWaveform.set_time(time_full);
            voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);

            auto excitation = complete_excitation(currentWaveform, voltageWaveform,
                                                   switchingFrequency, "Secondary 0");
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(cllcOpPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    // =========================================================================
    // Process Operating Points (all input voltages × all operating points)
    // =========================================================================

    std::vector<OperatingPoint> CllcConverter::process_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {

        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        // Clear extra-component waveform buffers; they accumulate per OP.
        extraCr1VoltageWaveforms.clear();
        extraCr1CurrentWaveforms.clear();
        extraCr2VoltageWaveforms.clear();
        extraCr2CurrentWaveforms.clear();
        extraLr1VoltageWaveforms.clear();
        extraLr1CurrentWaveforms.clear();
        extraLr2VoltageWaveforms.clear();
        extraLr2CurrentWaveforms.clear();
        extraTimeVectors.clear();

        double n = turnsRatios[0];

        // Calculate resonant parameters for ngspice and waveform generation
        CllcResonantParameters params = calculate_resonant_parameters();
        // Override with actual values provided
        params.turnsRatio = n;
        params.magnetizingInductance = magnetizingInductance;
        // Recalculate L1 from Lm and k
        params.primaryResonantInductance = magnetizingInductance / params.inductanceRatio;
        double L1 = params.primaryResonantInductance;
        double omega_r = 2.0 * M_PI * params.resonantFrequency;
        params.primaryResonantCapacitance = 1.0 / (omega_r * omega_r * L1);
        double a, b;
        resolve_resonant_ratios(a, b);
        params.secondaryResonantInductance = a * L1 / (n * n);
        params.secondaryResonantCapacitance = n * n * b * params.primaryResonantCapacitance;

        // Cache for get_extra_components_inputs.
        lastResonantParameters = params;

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto cllcOpPoint = get_operating_points()[opIndex];

                auto operatingPoint = process_operating_point_for_input_voltage(
                    inputVoltage, cllcOpPoint, n, magnetizingInductance, params);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(opIndex);
                }
                std::string flowStr = (cllcOpPoint.get_power_flow() == CllcPowerFlow::FORWARD) ?
                    " (forward)" : " (reverse)";
                name += flowStr;
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    std::vector<OperatingPoint> CllcConverter::process_operating_points(Magnetic magnetic) {
        CllcConverter::run_checks(_assertErrors);

        auto& settings = Settings::GetInstance();
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(settings.get_reluctance_model());
        double magnetizingInductance = magnetizingInductanceModel
            .calculate_inductance_from_number_turns_and_gapping(
                magnetic.get_mutable_core(), magnetic.get_mutable_coil())
            .get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();

        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    // =========================================================================
    // Ngspice Circuit Generation
    // =========================================================================

    /**
     * Generate a complete ngspice netlist for the CLLC converter.
     *
     * Circuit topology (forward mode):
     *
     *  Vin ─┬─ S1 ─┬─ S3 ─┬─ GND          GND ─┬─ Ds1 ─┬─ Ds3 ─┬─ Vout+
     *       │      A      │                     │       C       │
     *       │      │      │                     │       │       │
     *       └─ S2 ─┘─ S4 ─┘                    └─ Ds2 ─┘─ Ds4 ─┘
     *              B                                    D
     *
     *  A ── C1 ── L1 ── Lpri ──── B
     *                    ║ Lm
     *  C ── C2 ── L2 ── Lsec ──── D
     *
     * S1,S4 are controlled together (diagonal pair 1)
     * S2,S3 are controlled together (diagonal pair 2)
     * 50% duty cycle with dead time between transitions
     */
    std::string CllcConverter::generate_ngspice_circuit(
        double turnsRatio,
        const CllcResonantParameters& params,
        size_t inputVoltageIndex,
        size_t operatingPointIndex) {

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames_;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames_);

        if (inputVoltageIndex >= inputVoltages.size()) {
            throw std::invalid_argument("inputVoltageIndex out of range");
        }
        if (operatingPointIndex >= get_operating_points().size()) {
            throw std::invalid_argument("operatingPointIndex out of range");
        }

        double inputVoltage = inputVoltages[inputVoltageIndex];
        auto opPoint = get_operating_points()[operatingPointIndex];
        double switchingFrequency = get_value_or(opPoint.get_switching_frequency(), 0.0);
        if (switchingFrequency <= 0) {
            throw std::invalid_argument(
                "CLLC SPICE: operating point has invalid switching frequency");
        }
        double outputVoltage = opPoint.get_output_voltages()[0];
        double outputCurrent = opPoint.get_output_currents()[0];
        double n = turnsRatio;

        // P8 — bidirectional SPICE: in REVERSE the secondary side hosts the
        // EMF source and the primary side hosts the resistive load. Both
        // bridges are gated symmetrically off pwm1/pwm2 (same as forward); the
        // direction of net energy transfer flips because the source/load roles
        // swap. The TDA propagator and analytical waveforms remain forward-only
        // (FIXME-P8b: implement physically reverse propagator + reverse PtP
        // NRMSE fixture). This branch only re-physics the SPICE netlist so a
        // power-balance sanity test can confirm energy actually flows back.
        bool isReverse = (opPoint.get_power_flow() == CllcPowerFlow::REVERSE);

        double L1 = params.primaryResonantInductance;
        double C1 = params.primaryResonantCapacitance;
        double L2 = params.secondaryResonantInductance;
        double C2 = params.secondaryResonantCapacitance;
        double Lm = params.magnetizingInductance;
        double fr = params.resonantFrequency;

        double period = 1.0 / switchingFrequency;
        double halfPeriod = period / 2.0;
        double td = deadTime;

        // Ensure dead time doesn't exceed reasonable portion of half period
        if (td > halfPeriod * 0.05) {
            td = halfPeriod * 0.05;
        }

        double tOn = halfPeriod - td;

        // Simulation timing
        int periodsToExtract = get_num_periods_to_extract();
        int steadyPeriods = get_num_steady_state_periods();
        int totalPeriods = steadyPeriods + periodsToExtract;
        double simTime = totalPeriods * period;
        double startTime = steadyPeriods * period;
        // Time step: period/200 baseline, period/500 when fs > 1.5·fr
        // (super-resonant ringing demands tighter resolution).
        // CLLC_REWRITE_PLAN §6.2.
        double stepDivisor = (fr > 0 && switchingFrequency > 1.5 * fr) ? 500.0 : 200.0;
        double stepTime = period / stepDivisor;

        // Load resistance — divide-by-zero guard per plan §6.2 (R_load floor).
        double Rload = std::max(outputVoltage / std::max(outputCurrent, 1e-9), 1e-3);

        // P8 — primary-side load resistance for REVERSE mode. Sized so that
        // (Vin)^2 / R_pri matches the design-point throughput Vo·Io. Same
        // 1 mΩ floor as the secondary load.
        double powerTransferred = std::max(outputVoltage * outputCurrent, 1e-9);
        double Rload_pri = std::max(inputVoltage * inputVoltage / powerTransferred, 1e-3);

        // Secondary inductance for coupled inductor model
        double Lsec = Lm / (n * n);  // Lm referred to secondary for coupling

        std::ostringstream circuit;
        circuit << std::scientific;

        circuit << "* CLLC Bidirectional Resonant Converter - Generated by OpenMagnetics\n";
        circuit << "* Vin=" << inputVoltage << "V, Vout=" << outputVoltage
                << "V, Pout=" << (outputVoltage*outputCurrent) << "W, fs="
                << (switchingFrequency/1e3) << "kHz, fr=" << (fr/1e3) << "kHz\n";
        circuit << "* n=" << n << ", L1=" << (L1*1e6) << "uH, C1=" << (C1*1e9) << "nF\n";
        circuit << "* Lm=" << (Lm*1e6) << "uH, L2=" << (L2*1e6) << "uH, C2=" << (C2*1e9) << "nF\n\n";

        // DC Input — FORWARD: hard DC source on the primary rail. REVERSE:
        // primary rail is now a load; hold it with a bulk cap (UIC pre-charged
        // to inputVoltage) and dissipate returned power through Rload_pri,
        // sensed via a 0 V source so we can probe i(Vin_sense).
        if (isReverse) {
            circuit << "* DC Input (REVERSE — primary is load)\n";
            circuit << "Cin_pri vin_p 0 100u IC=" << inputVoltage << "\n";
            circuit << "Vin_sense vin_p vin_load 0\n";
            circuit << "Rload_pri vin_load 0 " << std::fixed << Rload_pri
                    << std::scientific << "\n\n";
        }
        else {
            circuit << "* DC Input\n";
            circuit << "Vin vin_p 0 " << inputVoltage << "\n\n";
        }

        // Switch model: VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg per plan §6.2
        circuit << "* Switch model (CLLC_REWRITE_PLAN §6.2)\n";
        circuit << ".model SW1 SW(VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg)\n\n";

        // PWM control signals for full bridge
        // Pair 1: S1 (high-side leg A), S4 (low-side leg B) - first half cycle
        // Pair 2: S2 (low-side leg A), S3 (high-side leg B) - second half cycle
        circuit << "* PWM control signals (complementary pairs with dead time)\n";
        circuit << "Vpwm1 pwm1 0 PULSE(0 5 0 10n 10n " << tOn
                << " " << period << ")\n";
        circuit << "Vpwm2 pwm2 0 PULSE(0 5 " << halfPeriod
                << " 10n 10n " << tOn << " " << period << ")\n\n";

        // Primary full bridge — 4 active switches.
        circuit << "* Primary Full Bridge (4 active switches)\n";
        circuit << "S1 vin_p node_a pwm1 0 SW1\n";
        circuit << "S2 node_a 0     pwm2 0 SW1\n";
        circuit << "S3 vin_p node_b pwm2 0 SW1\n";
        circuit << "S4 node_b 0     pwm1 0 SW1\n";
        // Per-switch snubbers (1k + 1nF) — 4 primary, 4 secondary = 8 total
        circuit << "* Primary snubbers (1k + 1nF across each switch)\n";
        circuit << "Rsn_S1 vin_p ns1 1k\n   Csn_S1 ns1 node_a 1n\n";
        circuit << "Rsn_S2 node_a ns2 1k\n  Csn_S2 ns2 0 1n\n";
        circuit << "Rsn_S3 vin_p ns3 1k\n   Csn_S3 ns3 node_b 1n\n";
        circuit << "Rsn_S4 node_b ns4 1k\n  Csn_S4 ns4 0 1n\n\n";

        // Bridge midpoint differential probe.
        circuit << "* Bridge differential voltage probe\n";
        circuit << "Evab vab 0 VALUE={V(node_a) - V(node_b)}\n\n";

        // Primary current sense
        circuit << "* Primary current sense\n";
        circuit << "Vpri_sense node_a pri_c1_in 0\n\n";

        // Primary resonant tank: C1 in series with L1
        circuit << "* Primary Resonant Tank (Cr1 series with Lr1)\n";
        circuit << "C_res1 pri_c1_in pri_l1_in " << C1 << "\n";
        circuit << "L_res1 pri_l1_in pri_trafo_in " << L1 << "\n\n";

        // Transformer: coupled inductors (Lpri = Lm, Lsec = Lm/n²). K = 0.9999
        // (never 1.0) so SPICE doesn't choke on a singular coupling matrix.
        circuit << "* Transformer (coupled inductors, K=0.9999 — never 1.0)\n";
        circuit << "Lpri pri_trafo_in node_b " << Lm << "\n";
        circuit << "Lsec sec_trafo_p sec_trafo_n " << Lsec << "\n";
        circuit << "Kpri_sec Lpri Lsec 0.9999\n\n";

        // Secondary resonant tank: L2 in series with C2
        circuit << "* Secondary Resonant Tank (Lr2 series with Cr2) + sense\n";
        circuit << "Vsec_sense sec_trafo_p sec_l2_in 0\n";
        circuit << "L_res2 sec_l2_in sec_c2_in " << L2 << "\n";
        circuit << "C_res2 sec_c2_in node_c " << C2 << "\n\n";

        // node_d = sec_trafo_n; high-Z DC path for solver convergence
        circuit << "* Secondary bridge reference node\n";
        circuit << "Vd_ref sec_trafo_n node_d 0\n";
        circuit << "Rdc_sec sec_trafo_n 0 1G\n\n";

        // ---- Active synchronous-rectifier full bridge on secondary ----
        // Per plan §6.1: in forward mode the SR is gated synchronously with
        // the primary PWM (pwm_s1 = pwm_p1, pwm_s2 = pwm_p2). This avoids
        // the chatter that an ideal-diode B-source emulator suffers at the
        // V(D,S) ≈ 0 transition (timestep collapse). Reverse mode (P8)
        // will drive the secondary as inverter and re-purpose the primary
        // as the synchronous rectifier with its own gating.
        //
        // Conduction map (forward, full-bridge):
        //   pwm1 high  → S1 (high-A) + S4 (low-B) ON  → node_a > node_b
        //              → secondary node_c > node_d
        //              → SR Sa (node_c → vout_p) and Sd (vout_n → node_d) ON
        //   pwm2 high  → opposite polarity  → SR Sb (vout_n → node_c) and
        //                                       Sc (node_d → vout_p) ON
        circuit << "* Secondary Active Synchronous Rectifier (4 SR switches)\n";
        circuit << "* pwm1 -> Sa, Sd    pwm2 -> Sb, Sc   (forward mode SR)\n";
        circuit << "Sa node_c vout_p pwm1 0 SW1\n";
        circuit << "Sb vout_n node_c pwm2 0 SW1\n";
        circuit << "Sc node_d vout_p pwm2 0 SW1\n";
        circuit << "Sd vout_n node_d pwm1 0 SW1\n";
        // Per-switch snubbers (1k + 1nF) on the SR side
        circuit << "* Secondary snubbers (1k + 1nF across each SR switch)\n";
        circuit << "Rsn_Sa node_c nsa 1k\n  Csn_Sa nsa vout_p 1n\n";
        circuit << "Rsn_Sb vout_n nsb 1k\n  Csn_Sb nsb node_c 1n\n";
        circuit << "Rsn_Sc node_d nsc 1k\n  Csn_Sc nsc vout_p 1n\n";
        circuit << "Rsn_Sd vout_n nsd 1k\n  Csn_Sd nsd node_d 1n\n\n";

        // Ground reference for secondary
        circuit << "* Secondary ground reference\n";
        circuit << "Vgnd_sec vout_n 0 0\n\n";

        // Output filter, sense source, and load — FORWARD path. In REVERSE
        // the secondary terminals host the EMF source instead (no Cout, no
        // Rload here). i(Vsec_src) is the current the secondary battery is
        // delivering into the converter.
        if (isReverse) {
            circuit << "* Secondary EMF (REVERSE — secondary is source)\n";
            circuit << "Vsec_src vout_p vout_n " << outputVoltage << "\n\n";
        }
        else {
            circuit << "* Output filter, current sense, and load\n";
            circuit << "Cout vout_p vout_n 100u IC=" << outputVoltage << "\n";
            circuit << "Vout_sense vout_p vout_load 0\n";
            circuit << "Rload vout_load vout_n " << std::fixed << Rload
                    << std::scientific << "\n\n";
        }

        // Transient analysis with UIC honoring the .ic statements below.
        circuit << "* Transient Analysis (UIC honors .ic pre-charge)\n";
        circuit << ".tran " << stepTime << " " << simTime
                << " " << startTime << " UIC\n\n";

        // Save signals — node-level and current branches that the
        // extractor maps into Primary/Secondary windings and ports. In
        // REVERSE the input/output current senses swap names: i(Vin_sense)
        // becomes the primary load current, i(Vsec_src) becomes the secondary
        // source current. Vin / Vout_sense don't exist in that branch.
        circuit << "* Save signals\n";
        circuit << ".save v(vab) v(vin_p) v(pri_trafo_in) v(node_a) v(node_b)\n";
        circuit << "+ v(node_c) v(node_d) v(sec_trafo_p) v(sec_trafo_n)\n";
        circuit << "+ v(vout_p) v(vout_n) v(pri_c1_in) v(sec_c2_in)\n";
        if (isReverse) {
            circuit << "+ i(Vin_sense) i(Vpri_sense) i(Vsec_sense) i(Vsec_src)\n\n";
        }
        else {
            circuit << "+ i(Vin) i(Vpri_sense) i(Vsec_sense) i(Vout_sense)\n\n";
        }

        // Solver options — verbatim from plan §6.2
        circuit << "* Solver options (CLLC_REWRITE_PLAN §6.2)\n";
        circuit << ".options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500\n";
        circuit << ".options METHOD=GEAR TRTOL=7\n\n";

        // Initial conditions: pre-charge output cap to Vo, leave resonant
        // caps at zero. SPICE node names for cap voltages: the node BETWEEN
        // C1 (pri_c1_in -- pri_l1_in) and BETWEEN C2 (sec_c2_in -- node_c).
        // In REVERSE there's no Cout (vout_p is held by Vsec_src), so pre-
        // charge the primary bulk cap to inputVoltage instead.
        circuit << "* Initial conditions (UIC pre-charge)\n";
        if (isReverse) {
            circuit << ".ic v(vin_p)=" << std::fixed << inputVoltage
                    << std::scientific << "\n";
        }
        else {
            circuit << ".ic v(vout_p)=" << std::fixed << outputVoltage
                    << std::scientific << "\n";
        }
        circuit << ".ic v(pri_l1_in)=0 v(sec_c2_in)=0\n\n";

        circuit << ".end\n";

        return circuit.str();
    }

    // =========================================================================
    // Simulation-based Operating Point Extraction
    // =========================================================================

    std::vector<OperatingPoint> CllcConverter::simulate_and_extract_operating_points(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance) {

        std::vector<OperatingPoint> operatingPoints;

        NgspiceRunner runner;
        if (!runner.is_available()) {
            std::cerr << "[CLLC] ngspice not available — falling back to analytical "
                         "operating-point generation" << std::endl;
            return process_operating_points(turnsRatios, magnetizingInductance);
        }

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        double n = turnsRatios[0];
        CllcResonantParameters params = calculate_resonant_parameters();
        params.turnsRatio = n;
        params.magnetizingInductance = magnetizingInductance;
        params.primaryResonantInductance = magnetizingInductance / params.inductanceRatio;
        double L1 = params.primaryResonantInductance;
        double omega_r = 2.0 * M_PI * params.resonantFrequency;
        params.primaryResonantCapacitance = 1.0 / (omega_r * omega_r * L1);
        double a, b;
        resolve_resonant_ratios(a, b);
        params.secondaryResonantInductance = a * L1 / (n * n);
        params.secondaryResonantCapacitance = n * n * b * params.primaryResonantCapacitance;

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto cllcOpPoint = get_operating_points()[opIndex];

                std::string netlist = generate_ngspice_circuit(n, params, inputVoltageIndex, opIndex);
                double switchingFrequency = get_value_or(cllcOpPoint.get_switching_frequency(), 0.0);

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = get_num_periods_to_extract();
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);

                if (!simResult.success) {
                    std::cerr << "[CLLC] ngspice run failed: " << simResult.errorMessage
                              << " — falling back to analytical for this op." << std::endl;
                    auto analytical = process_operating_point_for_input_voltage(
                        inputVoltages[inputVoltageIndex], cllcOpPoint, n,
                        magnetizingInductance, params);
                    std::string aname = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                    if (get_operating_points().size() > 1) {
                        aname += " op. point " + std::to_string(opIndex);
                    }
                    aname += " [analytical]";
                    analytical.set_name(aname);
                    operatingPoints.push_back(analytical);
                    continue;
                }

                // Map waveform names to winding excitations
                NgspiceRunner::WaveformNameMapping waveformMapping;

                // Primary winding: voltage across transformer primary = v(pri_trafo_in) - v(node_b)
                // Current through primary = i(Vpri_sense)
                waveformMapping.push_back({
                    {"voltage", "pri_trafo_in"},
                    {"current", "vpri_sense#branch"}
                });

                // Secondary winding: voltage = v(sec_trafo_p) - v(sec_trafo_n)
                // Current = i(Vsec_sense)
                waveformMapping.push_back({
                    {"voltage", "sec_trafo_p"},
                    {"current", "vsec_sense#branch"}
                });

                std::vector<std::string> windingNames = {"Primary", "Secondary"};
                std::vector<bool> flipCurrentSign = {false, false};

                OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                    simResult, waveformMapping, switchingFrequency,
                    windingNames, cllcOpPoint.get_ambient_temperature(),
                    flipCurrentSign);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt. (simulated)";
                if (get_operating_points().size() > 1) {
                    name += " op. point " + std::to_string(opIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }

        return operatingPoints;
    }

    // =========================================================================
    // Simulation-based Topology Waveform Extraction
    // =========================================================================

    std::vector<ConverterWaveforms> CllcConverter::simulate_and_extract_topology_waveforms(
        const std::vector<double>& turnsRatios,
        double magnetizingInductance,
        size_t numberOfPeriods) {

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

        double n = turnsRatios[0];
        CllcResonantParameters params = calculate_resonant_parameters();
        params.turnsRatio = n;
        params.magnetizingInductance = magnetizingInductance;
        params.primaryResonantInductance = magnetizingInductance / params.inductanceRatio;
        double L1 = params.primaryResonantInductance;
        double omega_r = 2.0 * M_PI * params.resonantFrequency;
        params.primaryResonantCapacitance = 1.0 / (omega_r * omega_r * L1);
        double a, b;
        resolve_resonant_ratios(a, b);
        params.secondaryResonantInductance = a * L1 / (n * n);
        params.secondaryResonantCapacitance = n * n * b * params.primaryResonantCapacitance;

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto opPoint = get_operating_points()[opIndex];

                std::string netlist = generate_ngspice_circuit(n, params, inputVoltageIndex, opIndex);
                double switchingFrequency = get_value_or(opPoint.get_switching_frequency(), 0.0);

                SimulationConfig config;
                config.frequency = switchingFrequency;
                config.extractOnePeriod = true;
                config.numberOfPeriods = numberOfPeriods;
                config.keepTempFiles = false;

                auto simResult = runner.run_simulation(netlist, config);
                if (!simResult.success) {
                    throw std::runtime_error("CLLC Simulation failed: " + simResult.errorMessage);
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

                // §5.1 converter-port stream. The "input" port is the
                // primary rail (vin_p) regardless of direction:
                //   FORWARD: vin_p is held by Vin source. i(Vin) flows
                //     OUT of source's + terminal → into the converter.
                //     Convention: converter draw = -i(Vin).
                //   REVERSE (P8b): vin_p is held by Cin_pri+Rload_pri.
                //     i(Vin_sense) flows from vin_p → vin_load (positive
                //     when load absorbs, i.e., reverse power flow). To
                //     keep the SAME sign convention "input current =
                //     current flowing INTO the converter from the input
                //     port", reverse mode emits -i(Vin_sense). The mean
                //     becomes negative, signalling power LEAVING the
                //     primary side — which is what reverse means.
                bool isReverse = (opPoint.get_power_flow() == CllcPowerFlow::REVERSE);
                wf.set_input_voltage(getWaveform("vin_p"));
                Waveform iInWf = isReverse
                    ? getWaveform("vin_sense#branch")
                    : getWaveform("vin#branch");
                auto& iInData = iInWf.get_mutable_data();
                for (auto& v : iInData) v = -v;
                wf.set_input_current(iInWf);

                wf.get_mutable_output_voltages().push_back(getWaveform("vout_p"));
                // Reconstruct Iout(t).
                //   FORWARD: Iout = Vout/Rload (DC by design, Cout filtered).
                //   REVERSE: secondary side hosts the EMF source Vsec_src;
                //     "output current" reflects current SOURCED by it,
                //     i.e., flowing OUT of Vsec_src + terminal. ngspice
                //     i(Vsec_src) is positive INTO + terminal, so source
                //     current = -i(Vsec_src).
                if (isReverse) {
                    Waveform ioutWf = getWaveform("vsec_src#branch");
                    auto& ioutData = ioutWf.get_mutable_data();
                    for (auto& v : ioutData) v = -v;
                    wf.get_mutable_output_currents().push_back(ioutWf);
                }
                else {
                    const double Vo = opPoint.get_output_voltages()[0];
                    const double Io = opPoint.get_output_currents()[0];
                    const double Rload = std::max(Vo / Io, 1e-3);
                    Waveform ioutWf = getWaveform("vout_p");
                    auto& ioutData = ioutWf.get_mutable_data();
                    for (auto& v : ioutData) v = v / Rload;
                    wf.get_mutable_output_currents().push_back(ioutWf);
                }

                results.push_back(wf);
            }
        }

        // Restore original value
        set_num_periods_to_extract(originalNumPeriodsToExtract);

        return results;
    }

    // =========================================================================
    // Extra components: Cr1, Cr2, Lr1, Lr2 ancillary MAS/CAS Inputs
    // =========================================================================
    //
    // Mirrors Llc::get_extra_components_inputs and Clllc::get_extra_components_inputs.
    // Cr1/Cr2 are always emitted (CAS::Inputs, RESONANT role). Lr1/Lr2 are
    // emitted as MAS::Inputs only when the corresponding integratedResonantInductor*
    // flag is false (or, in REAL mode, when the actual transformer leakage is
    // smaller than the required Lr — in which case the residual external
    // inductance is sized as Lr - Llk).
    std::vector<std::variant<Inputs, CAS::Inputs>>
    CllcConverter::get_extra_components_inputs(ExtraComponentsMode mode,
                                               std::optional<Magnetic> magnetic) {
        if (mode == ExtraComponentsMode::REAL && !magnetic.has_value()) {
            throw std::invalid_argument(
                "CllcConverter::get_extra_components_inputs: mode REAL requires a designed magnetic");
        }

        if (extraCr1VoltageWaveforms.empty()) {
            throw std::runtime_error(
                "CllcConverter::get_extra_components_inputs: call process_operating_points() first");
        }

        const double Lr1 = lastResonantParameters.primaryResonantInductance;
        const double Lr2 = lastResonantParameters.secondaryResonantInductance;
        double       Cr1 = lastResonantParameters.primaryResonantCapacitance;
        double       Cr2 = lastResonantParameters.secondaryResonantCapacitance;
        const double n   = lastResonantParameters.turnsRatio;
        const double fr  = lastResonantParameters.resonantFrequency;
        if (Lr1 <= 0 || Cr1 <= 0 || fr <= 0) {
            throw std::runtime_error(
                "CllcConverter::get_extra_components_inputs: cached resonant parameters invalid "
                "(call process_operating_points first)");
        }

        bool integratedLr1 = get_integrated_resonant_inductor1().value_or(true);
        bool integratedLr2 = get_integrated_resonant_inductor2().value_or(true);

        // REAL mode: rebalance Cr against the actual primary-side leakage so
        // the tank still resonates at fr when Lr1 is integrated. Mirrors
        // Llc.cpp:1631-1640.
        double Llk_pri = 0.0;
        if (mode == ExtraComponentsMode::REAL) {
            auto leakageOutput = LeakageInductance().calculate_leakage_inductance_all_windings(
                magnetic.value(), fr);
            auto perWinding = leakageOutput.get_leakage_inductance_per_winding();
            if (!perWinding.empty()) {
                Llk_pri = resolve_dimensional_values(perWinding[0]);
            }
            double Lr1_eff = (Llk_pri >= Lr1) ? Llk_pri : Lr1;
            Cr1 = 1.0 / (4.0 * M_PI * M_PI * fr * fr * Lr1_eff);
            // Maintain the symmetric / asymmetric ratio b = Cr2 / (n²·Cr1).
            double a_unused, b;
            resolve_resonant_ratios(a_unused, b);
            (void)a_unused;
            Cr2 = n * n * b * Cr1;
        }

        auto& ops = get_operating_points();
        std::vector<std::variant<Inputs, CAS::Inputs>> result;

        auto buildCapInputs = [&](const std::string& name,
                                  double C,
                                  const std::vector<Waveform>& vWfms,
                                  const std::vector<Waveform>& iWfms) {
            if (vWfms.size() != iWfms.size()) {
                throw std::runtime_error(
                    "CllcConverter::get_extra_components_inputs: " + name +
                    " voltage/current waveform count mismatch");
            }
            CAS::Inputs casInputs;
            CAS::DesignRequirements dr;

            CAS::DimensionWithTolerance capacitance;
            capacitance.set_nominal(C);
            dr.set_capacitance(capacitance);

            double peakV = 0.0;
            for (const auto& w : vWfms) {
                for (double v : w.get_data()) {
                    peakV = std::max(peakV, std::abs(v));
                }
            }
            // 1.5x peak rating — bidirectional CLLC sees ringing on direction
            // reversal (matches Clllc convention).
            dr.set_rated_voltage(peakV * 1.5);
            dr.set_role(CAS::Application::RESONANT);
            dr.set_name(name);
            casInputs.set_design_requirements(dr);

            std::vector<CAS::TwoTerminalOperatingPoint> casOps;
            for (size_t i = 0; i < vWfms.size(); ++i) {
                CAS::TwoTerminalOperatingPoint op;
                CAS::OperatingPointExcitation exc;
                double fsw = (i < ops.size())
                                 ? get_value_or(ops[i].get_switching_frequency(), fr)
                                 : fr;
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

        buildCapInputs("Cr1_resonantCapacitor_primary", Cr1,
                       extraCr1VoltageWaveforms, extraCr1CurrentWaveforms);
        buildCapInputs("Cr2_resonantCapacitor_secondary", Cr2,
                       extraCr2VoltageWaveforms, extraCr2CurrentWaveforms);

        // Discrete external resonant inductors. In REAL mode the magnetic's
        // own leakage absorbs part (or all) of Lr.
        auto buildInductorInputs = [&](const std::string& name,
                                       double L,
                                       const std::vector<Waveform>& iWfms,
                                       const std::vector<Waveform>& vWfms) {
            if (iWfms.size() != vWfms.size()) {
                throw std::runtime_error(
                    "CllcConverter::get_extra_components_inputs: " + name +
                    " current/voltage waveform count mismatch");
            }
            Inputs masInputs;
            DesignRequirements dr;
            DimensionWithTolerance inductance;
            inductance.set_nominal(L);
            dr.set_magnetizing_inductance(inductance);
            dr.set_name(name);
            dr.set_topology(Topologies::CLLC_RESONANT_CONVERTER);
            dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
            dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
            masInputs.set_design_requirements(dr);

            std::vector<OperatingPoint> masOps;
            for (size_t i = 0; i < iWfms.size(); ++i) {
                OperatingPoint op;
                double fsw = (i < ops.size())
                                 ? get_value_or(ops[i].get_switching_frequency(), fr)
                                 : fr;
                auto exc = complete_excitation(iWfms[i], vWfms[i], fsw, "Primary");
                op.get_mutable_excitations_per_winding().push_back(exc);
                masOps.push_back(op);
            }
            masInputs.set_operating_points(masOps);
            result.emplace_back(std::move(masInputs));
        };

        // Lr1 sizing: IDEAL → full Lr1 when not integrated; REAL → residual
        // beyond the magnetic's own primary leakage.
        double Lr1_external = 0.0;
        if (mode == ExtraComponentsMode::IDEAL && !integratedLr1) {
            Lr1_external = Lr1;
        }
        else if (mode == ExtraComponentsMode::REAL) {
            Lr1_external = (Llk_pri < Lr1) ? (Lr1 - Llk_pri) : 0.0;
        }
        if (Lr1_external > 0.0 && !extraLr1VoltageWaveforms.empty()) {
            buildInductorInputs("Lr1_seriesInductor_primary", Lr1_external,
                                extraLr1CurrentWaveforms, extraLr1VoltageWaveforms);
        }

        // Lr2 sizing: same logic on the secondary side. Lr2 is referred to
        // the secondary, so the leakage to subtract would be the secondary
        // winding's leakage; for now mirror the primary heuristic when REAL
        // and absorb fully when IDEAL+integrated. (Tightening to true
        // secondary-leakage subtraction is deferred to P7.)
        double Lr2_external = 0.0;
        if (mode == ExtraComponentsMode::IDEAL && !integratedLr2) {
            Lr2_external = Lr2;
        }
        else if (mode == ExtraComponentsMode::REAL) {
            Lr2_external = integratedLr2 ? 0.0 : Lr2;
        }
        if (Lr2_external > 0.0 && !extraLr2VoltageWaveforms.empty()) {
            buildInductorInputs("Lr2_seriesInductor_secondary", Lr2_external,
                                extraLr2CurrentWaveforms, extraLr2VoltageWaveforms);
        }

        return result;
    }

    // =========================================================================
    // Advanced CLLC Converter
    // =========================================================================

    Inputs AdvancedCllcConverter::process() {
        CllcConverter::run_checks(_assertErrors);

        Inputs inputs;

        double magnetizingInductance = get_desired_magnetizing_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();

        if (turnsRatios.empty()) {
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS,
                "CLLC: desiredTurnsRatios must not be empty");
        }

        // Build design requirements
        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();

        for (auto tr : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(tr, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(magnetizingInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        designRequirements.set_isolation_sides(
            Topology::create_isolation_sides(1, false));
        designRequirements.set_topology(Topologies::CLLC_RESONANT_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        // Generate operating points
        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;
        collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltagesNames);

        CllcResonantParameters params = calculate_resonant_parameters();
        params.turnsRatio = turnsRatios[0];
        params.magnetizingInductance = magnetizingInductance;

        // Override resonant components if user specified them
        if (get_desired_resonant_inductance_primary()) {
            params.primaryResonantInductance = get_desired_resonant_inductance_primary().value();
        } else {
            params.primaryResonantInductance = magnetizingInductance / params.inductanceRatio;
        }

        double L1 = params.primaryResonantInductance;

        if (get_desired_resonant_capacitance_primary()) {
            params.primaryResonantCapacitance = get_desired_resonant_capacitance_primary().value();
        } else {
            double omega_r = 2.0 * M_PI * params.resonantFrequency;
            params.primaryResonantCapacitance = 1.0 / (omega_r * omega_r * L1);
        }

        double n = turnsRatios[0];
        double a, b;
        resolve_resonant_ratios(a, b);

        if (get_desired_resonant_inductance_secondary()) {
            params.secondaryResonantInductance = get_desired_resonant_inductance_secondary().value();
        } else {
            params.secondaryResonantInductance = a * L1 / (n * n);
        }

        if (get_desired_resonant_capacitance_secondary()) {
            params.secondaryResonantCapacitance = get_desired_resonant_capacitance_secondary().value();
        } else {
            params.secondaryResonantCapacitance = n * n * b * params.primaryResonantCapacitance;
        }

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t opIndex = 0; opIndex < get_operating_points().size(); ++opIndex) {
                auto cllcOpPoint = get_operating_points()[opIndex];

                auto operatingPoint = process_operating_point_for_input_voltage(
                    inputVoltage, cllcOpPoint, n, magnetizingInductance, params);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(opIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

} // namespace OpenMagnetics
