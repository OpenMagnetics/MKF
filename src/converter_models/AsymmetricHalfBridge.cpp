#include "converter_models/AsymmetricHalfBridge.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/NgspiceRunner.h"
#include "support/Settings.h"
#include "support/Exceptions.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace OpenMagnetics {

// =========================================================================
// Constructors
// =========================================================================
AsymmetricHalfBridge::AsymmetricHalfBridge(const json& j) {
    from_json(j, *static_cast<MAS::AsymmetricHalfBridge*>(this));
    // Validate range-bounded fields immediately so malformed JSON is
    // rejected at construction (matches the user-facing contract: "no
    // fallbacks, no defaults, no silent shortcuts"). The MAS schema
    // setters do not enforce duty_cycle ∈ [0,1] for OperatingPoint, so
    // we re-check here.
    for (auto& op : get_operating_points()) {
        const double D = op.get_duty_cycle();
        if (!(D >= 0.0 && D <= 1.0)) {
            throw std::runtime_error(
                "AsymmetricHalfBridge: dutyCycle must lie in [0, 1] "
                "(got " + std::to_string(D) + ")");
        }
    }
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
[[noreturn, maybe_unused]] void not_implemented(const char* method, const char* phase) {
    throw std::runtime_error(
        std::string("AsymmetricHalfBridge: ") + method +
        " not implemented yet (see ASYMMETRIC_HALF_BRIDGE_PLAN.md " +
        phase + ")");
}
}

DesignRequirements AsymmetricHalfBridge::process_design_requirements() {
    // -------------------------------------------------------------------------
    // P6 deliverable. Sizes, in order:
    //
    //   S1  n   = compute_turns_ratio(Vin_min, Vo, D_max, Vd, rectType)
    //               -- worst-case Vin → still hits Vo at the operating duty.
    //   S2  Lo  = compute_lo_min(Vo, D_nom, dILo=0.30·Io, fsw)
    //               -- for CD: Io_per_inductor = Io/2.
    //   S4  Lm  = compute_lm_min_for_zvs(Vin_max, D_nom, Tsw, Im_target)
    //               with Im_target = 0.10 · I_pri,Io  (mirrors PSHB heuristic).
    //   S3  Cb  = compute_cb_min(I_pri_pk, D_nom, dVCb=0.05·V_Cb, fsw)
    //               -- enforces ΔV_Cb ≤ 5 % so flux excursion stays symmetric
    //               (AHB plan §12 mistake #1).
    //   S5  Co  = ΔILo / (8·fsw·ΔVo_target),  ΔVo_target = 0.01·Vo
    //               -- standard buck-style output ripple sizing.
    //
    // Multi-output (V6) and ahbFlyback (V5) are deferred to P11 / P10. The
    // emitted DesignRequirements always carries (a) the per-secondary
    // turns_ratios, (b) the magnetizing inductance, (c) the leakage
    // inductance when use_leakage_inductance is set, (d) the topology tag,
    // and (e) the isolation_sides built from numSecondaries. Lo/Cb/Co flow
    // through get_extra_components_inputs (P9) — they are stashed in the
    // computed* members so the P9 pass can pick them up without recomputing.
    //
    // The routine is idempotent: calling twice with the same inputs leaves
    // the computed* members at the same values and emits structurally
    // identical DesignRequirements (matches DAB/PSHB convention).
    // -------------------------------------------------------------------------

    AhbRectifierType rect = get_rectifier_type().value_or(
        AhbRectifierType::CENTER_TAPPED);

    // ---- Pull operating-point envelope ----
    auto& inputVoltage = get_input_voltage();
    if (!inputVoltage.get_nominal().has_value() &&
        !inputVoltage.get_minimum().has_value() &&
        !inputVoltage.get_maximum().has_value())
        throw std::runtime_error(
            "AsymmetricHalfBridge::process_design_requirements: inputVoltage "
            "must specify nominal, minimum, or maximum");

    const double Vin_nom = inputVoltage.get_nominal().value_or(
        (inputVoltage.get_minimum().value_or(0.0) +
         inputVoltage.get_maximum().value_or(0.0)) / 2.0);
    const double Vin_min = inputVoltage.get_minimum().value_or(Vin_nom);
    const double Vin_max = inputVoltage.get_maximum().value_or(Vin_nom);

    auto& ops = get_operating_points();
    if (ops.empty())
        throw std::runtime_error(
            "AsymmetricHalfBridge::process_design_requirements: at least "
            "one operating point required");
    const bool isMultiOutput = (ops[0].get_output_voltages().size() > 1 ||
                                ops[0].get_output_currents().size() > 1);

    // Effective single-output proxy used for primary sizing under
    // multi-output (V6). We size the primary against TOTAL output power
    // and a representative voltage (max of the per-output voltages, so
    // the worst-case Vsec swing is covered).
    double Vo_eff = ops[0].get_output_voltages()[0];
    double Io_eff = ops[0].get_output_currents()[0];
    if (isMultiOutput) {
        double P_total = 0.0;
        double V_rep   = 0.0;
        for (size_t i = 0; i < ops[0].get_output_voltages().size(); ++i) {
            const double Vo_i = ops[0].get_output_voltages()[i];
            const double Io_i = ops[0].get_output_currents()[i];
            P_total += Vo_i * Io_i;
            V_rep    = std::max(V_rep, Vo_i);
        }
        Vo_eff = V_rep;
        Io_eff = (V_rep > 0.0) ? P_total / V_rep : 0.0;
    }
    const double Vo  = Vo_eff;
    const double Io  = Io_eff;
    const double fsw = ops[0].get_switching_frequency();
    if (!(Vo > 0.0))  throw std::runtime_error("Vo must be > 0");
    if (!(Io > 0.0))  throw std::runtime_error("Io must be > 0");
    if (!(fsw > 0.0)) throw std::runtime_error("fsw must be > 0");

    // Operating duty used for sizing: prefer the OP's explicit value (the
    // user has chosen a target M = Vo/Vin), fall back to the schema's
    // maximumDutyCycle (default 0.45). Both are validated in (0, 1) by
    // setters / run_checks.
    const double D_op  = ops[0].get_duty_cycle();
    const double D_max = get_maximum_duty_cycle().value_or(0.45);
    const double D_nom = (D_op > 0.0 && D_op < 1.0) ? D_op : D_max;
    const double Tsw   = 1.0 / fsw;
    const double Vd    = computedDiodeVoltageDrop;

    // ---- S1: turns ratio (sized at min Vin to guarantee headroom) ----
    const double n = compute_turns_ratio(Vin_min, Vo, D_max, Vd, rect);
    computedTurnsRatio = n;

    // ---- AHB-Flyback (V5) sizing branch ---------------------------------
    // Two-switch / active-clamp flyback: NO Lo, NO Cb. Lm carries the
    // energy storage (mean(iLm) > 0). Output cap sized for flyback ripple
    // (rectified secondary current is iLm·n during (1-D)·Tsw only, so
    // Co supplies the full Io during D·Tsw).
    if (rect == AhbRectifierType::AHB_FLYBACK) {
        // Magnetizing inductance: target 30% pp ripple ratio of mean iLm.
        //   I_Lm,avg  = Io / (n · (1-D))                  (flyback CCM)
        //   dILm,pp   = Vin · D · Tsw / Lm
        //   ratio = dILm,pp / I_Lm,avg = Vin·D·n·(1-D)·Tsw / (Lm·Io)
        //   => Lm = Vin·D·n·(1-D)·Tsw / (0.30 · Io)
        // Sized at Vin_max so the worst-case ripple stays within budget.
        const double ratio_target = 0.30;
        const double Lm_v5 = (Vin_max * D_nom * n * (1.0 - D_nom) * Tsw) /
                              (ratio_target * std::max(Io, 1e-6));
        computedMagnetizingInductance = Lm_v5;
        computedOutputInductance      = 0.0;   // no Lo
        computedDcBlockingCapacitance = 0.0;   // no Cb (active-clamp flyback)

        // Output cap: standard flyback Co sizing — Co supplies Io during
        // the entire ON interval D·Tsw plus the discontinuous-conduction
        // gap. Worst case: Co = Io · D / (fsw · ΔVo).
        const double dVo_target_v5 = std::max(0.01 * Vo, 1e-3);
        computedOutputCapacitance   = (Io * D_nom) / (fsw * dVo_target_v5);

        computedDutyCycle = D_nom;
        computedDeadTime  = compute_dead_time(
            computedLeakageInductance + Lm_v5, mosfetCoss);
    }
    else {
        // ---- S2: output inductor ----
        const bool isCD = (rect == AhbRectifierType::CURRENT_DOUBLER);
        const double Io_per_inductor = isCD ? Io / 2.0 : Io;
        const double dILo_target = 0.30 * Io_per_inductor;
        const double Lo = compute_lo_min(Vo, D_nom, dILo_target, fsw);
        computedOutputInductance = Lo;

        // ---- S4: magnetizing inductance for ZVS at min load ----
        // Mirrors PSHB heuristic: target Im_pk = 10 % of primary load current.
        // Sized at Vin_max so the volt-second drives a sufficient Im across
        // the half-cycle even at the worst-case operating point.
        const double Io_pri = Io / n;
        const double Im_target = std::max(0.10 * Io_pri, 1e-3);
        const double Lm = compute_lm_min_for_zvs(Vin_max, D_nom, Tsw, Im_target);
        computedMagnetizingInductance = Lm;

        // ---- S3: DC-blocking cap (≤ 5 % ripple of V_Cb,nom) ----
        // I_pri peak ≈ Io/n + Im_pk. Use the worst-case Vin_max for V_Cb sizing.
        const double dILm_pp_max  = (1.0 - D_nom) * Vin_max * D_nom * Tsw / Lm;
        const double Iprimary_pk  = Io_pri + dILm_pp_max / 2.0;
        const double VCb_nom      = (1.0 - D_nom) * Vin_nom;
        const double dVCb_target  = std::max(0.05 * VCb_nom, 1e-3);
        const double Cb = compute_cb_min(Iprimary_pk, D_nom, dVCb_target, fsw);
        computedDcBlockingCapacitance = Cb;

        // ---- S5: output cap (1 % output ripple) ----
        // Standard buck cap sizing: ΔVo = ΔILo / (8 · fsw · Co)  →
        // Co = ΔILo / (8 · fsw · ΔVo_target). For CD the rectified ripple
        // frequency at the cap is 2·fsw, but we conservatively use fsw here
        // — exact value is tuned in P9 simulation if needed.
        const double dILo_pp_total = isCD ? 2.0 * dILo_target : dILo_target;
        const double dVo_target    = std::max(0.01 * Vo, 1e-3);
        computedOutputCapacitance  = dILo_pp_total / (8.0 * fsw * dVo_target);

        // ---- Dead-time (resonant tank Llk + Lm per AHB physics, P5 refinement) ----
        computedDutyCycle = D_nom;
        computedDeadTime  = compute_dead_time(computedLeakageInductance + Lm,
                                              mosfetCoss);
    }
    const double Lm_for_dr = computedMagnetizingInductance;

    // ---- Build DesignRequirements ----
    DesignRequirements designRequirements;
    designRequirements.get_mutable_turns_ratios().clear();
    {
        // Per-output turns ratios: each secondary feeds its own Vo_i and
        // therefore needs its own n_i. For single-output this collapses
        // to one ratio (with CT duplicating it).
        std::vector<double> n_per_output;
        if (isMultiOutput) {
            for (size_t i = 0; i < ops[0].get_output_voltages().size(); ++i) {
                const double Vo_i = ops[0].get_output_voltages()[i];
                if (!(Vo_i > 0.0))
                    throw std::runtime_error(
                        "AsymmetricHalfBridge::process_design_requirements: "
                        "outputVoltages[" + std::to_string(i) + "] must be > 0");
                n_per_output.push_back(
                    compute_turns_ratio(Vin_min, Vo_i, D_max, Vd, rect));
            }
        } else {
            n_per_output.push_back(n);
        }

        for (double n_i : n_per_output) {
            DimensionWithTolerance nTol;
            nTol.set_nominal(roundFloat(n_i, 2));
            designRequirements.get_mutable_turns_ratios().push_back(nTol);
        }
        // For CENTER_TAPPED there are two physical secondary windings
        // (Sec_a, Sec_b) per output and several places in the coil/adviser
        // pipeline size winding-count from turns_ratios.size()+1 (e.g.
        // Coil::get_patterns, MagneticAdviser::numberWindings,
        // MagneticFilterEstimatedCost). Duplicate the FIRST ratio to keep
        // that arithmetic consistent with the single-secondary CT usage;
        // multi-output CT is not commonly built, so we only duplicate the
        // primary CT entry to preserve V1 behaviour.
        if (rect == AhbRectifierType::CENTER_TAPPED && !isMultiOutput) {
            DimensionWithTolerance nTol;
            nTol.set_nominal(roundFloat(n_per_output[0], 2));
            designRequirements.get_mutable_turns_ratios().push_back(nTol);
        }
    }

    DimensionWithTolerance lmTol;
    lmTol.set_nominal(roundFloat(Lm_for_dr, 10));
    designRequirements.set_magnetizing_inductance(lmTol);

    if (computedLeakageInductance > 0.0) {
        std::vector<DimensionWithTolerance> leakageReqs;
        DimensionWithTolerance llkTol;
        llkTol.set_nominal(roundFloat(computedLeakageInductance, 10));
        leakageReqs.push_back(llkTol);
        designRequirements.set_leakage_inductance(leakageReqs);
    }

    designRequirements.set_topology(Topologies::ASYMMETRIC_HALF_BRIDGE_CONVERTER);
    // Isolation sides:
    //   - For all rectifier types the primary takes side 0.
    //   - For CENTER_TAPPED there are TWO physical secondary windings
    //     (Sec_a, Sec_b) but they form a single galvanic secondary; both
    //     sit on side 1 and are then linked via wound_with so the section
    //     algorithm treats them as one group.
    //   - For FB / CD / Flyback there is a single secondary winding on
    //     side 1 (output_currents().size() == 1).
    // Topology::create_isolation_sides hands out a fresh side per
    // secondary (i+1), which would put Sec_a on side 1 and Sec_b on
    // side 2 — that does not match the CT physical layout and broke
    // section allocation in get_patterns/get_ordered_sections (one
    // winding ended up with zero sections, tripping
    // "Number of slots cannot be less than 1").
    {
        std::vector<IsolationSide> isolationSides;
        isolationSides.push_back(get_isolation_side_from_index(0));  // primary
        if (rect == AhbRectifierType::CENTER_TAPPED) {
            isolationSides.push_back(get_isolation_side_from_index(1));  // Sec_a
            isolationSides.push_back(get_isolation_side_from_index(1));  // Sec_b (same side)
        } else {
            for (size_t i = 0; i < ops[0].get_output_currents().size(); ++i) {
                isolationSides.push_back(get_isolation_side_from_index(i + 1));
            }
        }
        designRequirements.set_isolation_sides(isolationSides);
    }

    return designRequirements;
}


// =========================================================================
// process_operating_point_for_input_voltage   (P3 deliverable, V1-CT secondary)
//
// CCM piecewise-linear analytical model with the canonical four-segment
// AHB sub-interval structure (dead-times collapsed to instantaneous):
//
//   Sub-interval | length    | Q1 | Q2 | v_pri          | i_Lm slope
//   -------------+-----------+----+----+----------------+----------------
//   A            | D·Tsw     | ON | OFF| +(1-D)·Vin     | +(1-D)·Vin/Lm
//   C            | (1-D)·Tsw | OFF| ON | -D·Vin         | -D·Vin/Lm
//
// Volt-second balance is automatic:
//   λ+ = (1-D)·Vin · D·Tsw      ;  λ- = -D·Vin · (1-D)·Tsw   →   λ+ + λ- = 0.
//
// Magnetizing current is symmetric about zero:
//   ΔI_Lm,pp = (1-D)·Vin·D·Tsw / Lm
//   I_Lm(0)  = -ΔI_Lm,pp / 2
//
// Centre-tapped secondary (V1):
//   - During interval A the top diode conducts. v_sec_rect = +(1-D)·Vin / n
//   - During interval C the bottom diode conducts. v_sec_rect = +D·Vin / n
//   In steady state Vo = 2·D·(1-D)·Vin / n  (closed-form gain at this duty),
//   so the Lo ripple cancels exactly across A and C.
//
// Output inductor (CCM):
//   ΔI_Lo,pp = (Vsec_pos - Vo) · D · Tsw / Lo
//             ≡ (Vo - Vsec_neg) · (1-D) · Tsw / Lo   (volt-sec balance)
//   I_Lo(0)  = Io - ΔI_Lo,pp / 2     (trough at start of interval A)
//
// Primary current:
//   In A:  i_pri = +i_Lo / n + i_Lm    (top diode reflects load to primary)
//   In C:  i_pri = -i_Lo / n + i_Lm    (bottom diode; primary swings negative)
//
// Switch RMS — ASYMMETRIC by construction:
//   I_Q1,rms ≃ |i_pri,A,rms| · √D       (Q1 conducts only during interval A)
//   I_Q2,rms ≃ |i_pri,C,rms| · √(1-D)   (Q2 conducts only during interval C)
//   At D = 0.3 the textbook ratio I_Q1,rms / I_Q2,rms = √(D / (1-D))
//   (Imbertson-Mohan §IV.A) holds within the linear-ripple approximation.
//
// Diagnostics populated (per Guide §2.4):
//   lastDutyCycle, lastConversionRatio, lastDcBlockingCapVoltage,
//   lastPrimaryPeakVoltagePositive/Negative, lastSwitchPeakVoltageQ1/Q2,
//   lastSwitchRmsCurrentQ1/Q2, lastZvsMargin, lastResonantTransitionTime,
//   lastSteadyStateFluxExcursion, lastMagnetizingCurrentRipple,
//   lastOutputInductorRipple, lastOperatingMode, lastRectifierType.
//
// DCM detection: I_Lo,trough < 0  →  lastOperatingMode = 1 (DCM); P5 will
// model the discontinuous waveform properly. P3 emits CCM linear waveform
// regardless and tags the OP — no silent re-mapping.
//
// Multi-secondary, rectifier types CD / FB / Flyback are deferred to P4
// and P11; P3 throws if any non-CENTER_TAPPED rectifier is requested or
// if more than one secondary is supplied.
//
// References:
//   - Imbertson & Mohan, IEEE TIA 29(1) 1993 (canonical asymmetric duty)
//   - TI SLUP223 §3 (centre-tapped AHB worked example)
//   - ASYMMETRIC_HALF_BRIDGE_PLAN.md §15 (sub-interval table)
// =========================================================================
namespace {
inline double rms_of(const std::vector<double>& y, size_t i0, size_t i1,
                     const std::vector<double>& t) {
    // Trapezoidal integral of y² over [t[i0], t[i1]] (i1 inclusive).
    if (i1 <= i0) return 0.0;
    double sum = 0.0;
    for (size_t k = i0; k < i1; ++k) {
        const double dt = t[k + 1] - t[k];
        sum += 0.5 * (y[k] * y[k] + y[k + 1] * y[k + 1]) * dt;
    }
    const double duration = t[i1] - t[i0];
    if (duration <= 0.0) return 0.0;
    return std::sqrt(sum / duration);
}
} // namespace

OperatingPoint AsymmetricHalfBridge::process_operating_point_for_input_voltage(
    double inputVoltage,
    const AhbOperatingPoint& opPoint,
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    // ---- Input validation (no defaults, no fallbacks: throw loud) ----
    if (!(inputVoltage > 0.0) || !std::isfinite(inputVoltage))
        throw std::invalid_argument(
            "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
            "inputVoltage must be > 0");
    if (turnsRatios.empty())
        throw std::invalid_argument(
            "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
            "turnsRatios must contain at least one entry");
    for (double tr : turnsRatios)
        if (!(tr > 0.0) || !std::isfinite(tr))
            throw std::invalid_argument(
                "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
                "all turnsRatios must be > 0");
    if (!(magnetizingInductance > 0.0) || !std::isfinite(magnetizingInductance))
        throw std::invalid_argument(
            "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
            "magnetizingInductance must be > 0");

    // V5 (AHB-Flyback) is a different storage topology — Lm carries the
    // energy, no Lo, no series Cb (active-clamp / two-switch flyback).
    // Dispatched to a dedicated branch below. V6 multi-output is handled
    // inline by extending the per-secondary loop.
    AhbRectifierType rect = get_rectifier_type().value_or(
        AhbRectifierType::CENTER_TAPPED);
    const size_t numOutputs = opPoint.get_output_voltages().size();
    if (opPoint.get_output_currents().size() != numOutputs)
        throw std::invalid_argument(
            "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
            "outputVoltages / outputCurrents length mismatch");
    if (numOutputs == 0)
        throw std::invalid_argument(
            "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
            "operating point must specify at least one output");
    if (turnsRatios.size() != numOutputs)
        throw std::invalid_argument(
            "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
            "turnsRatios.size() (" + std::to_string(turnsRatios.size()) +
            ") must match outputs (" + std::to_string(numOutputs) + ")");

    // Multi-output cross-regulation in AHB is materially worse than buck-
    // family converters (asymmetric secondary voltage, V_Cb shared across
    // all secondaries). Per Guide §4.7 we emit a one-shot warning so the
    // user knows the per-output current waveforms are an approximation.
    if (numOutputs > 1) {
        static thread_local bool ahbMultiOutWarned = false;
        if (!ahbMultiOutWarned) {
            std::cerr << "[AHB] multi-output configuration detected ("
                      << numOutputs << " outputs) — per-output current "
                         "waveforms are approximate (load-share projection "
                         "of primary reflected current). Cross-regulation "
                         "is poor in AHB without per-secondary feedback; "
                         "supply per-secondary leakage inductance for a "
                         "more accurate magnetic-loss estimate."
                      << std::endl;
            ahbMultiOutWarned = true;
        }
    }

    const double Vin = inputVoltage;
    const double D   = opPoint.get_duty_cycle();
    // For multi-output, project ALL outputs to a single equivalent load
    // referenced to output #0 (the controlling output). The primary side
    // sees the total reflected power; per-output secondary current
    // waveforms are recovered after the standard OP is built using the
    // load-share formula:  i_sec_k = (Vo_k·Io_k / Σ Vo_j·Io_j) · n_k · iPri.
    // (Mirrors Dab.cpp:842-906 multi-output projection.)
    const auto& outV = opPoint.get_output_voltages();
    const auto& outI = opPoint.get_output_currents();
    const double Vo  = outV[0];
    double IoEff = outI[0];
    double totalPower = 0.0;
    for (size_t k = 0; k < outV.size() && k < outI.size(); ++k)
        totalPower += outV[k] * outI[k];
    if (outV.size() > 1 && Vo > 0.0)
        IoEff = totalPower / Vo;
    const double Io  = IoEff;
    const double fsw = opPoint.get_switching_frequency();
    if (!(D > 0.0 && D < 1.0))
        throw std::invalid_argument(
            "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
            "dutyCycle must lie in (0, 1) for waveform generation");
    if (!(fsw > 0.0))
        throw std::invalid_argument(
            "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
            "switchingFrequency must be > 0");

    // ---------------------------------------------------------------------
    // V5 AHB-Flyback path (no Lo, no Cb, Lm = energy storage).
    // ---------------------------------------------------------------------
    // Two-switch / active-clamp flyback interpretation (selected per
    // ASYMMETRIC_HALF_BRIDGE_PLAN.md §4 V5 row):
    //   - During Q1 ON (D·Tsw): v_pri = +Vin, secondary diode OFF, Lm
    //     ramps with slope +Vin/Lm.
    //   - During Q2 ON ((1-D)·Tsw): secondary diode forward-biased,
    //     reflected v_pri = -Vo·n, Lm ramps with slope -Vo·n/Lm.
    //   - Volt-second balance: Vin·D = Vo·n·(1-D) → Vo = D·Vin/((1-D)·n)
    //     (matches static helper compute_conversion_ratio for AHB_FLYBACK).
    //   - Mean iLm = Io / (n·(1-D))  (CCM flyback). Mean i_pri = D·iLm
    //     (NON-zero because no series Cb).
    //   - During Q2 ON, primary current goes to ZERO (ampere-turn balance:
    //     all transformer current is on the secondary side).
    if (rect == AhbRectifierType::AHB_FLYBACK) {
        if (numOutputs != 1)
            throw std::invalid_argument(
                "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
                "AHB_FLYBACK does not support multi-output");

    const double n   = turnsRatios[0];
    const double Lm  = magnetizingInductance;


        const double Tsw = 1.0 / fsw;
        const double tA  = D * Tsw;
        const double tC  = (1.0 - D) * Tsw;

        // Primary voltages: full Vin during Q1, -Vo·n during Q2 (clamped).
        const double Vpri_pos_v5 = +Vin;
        const double Vpri_neg_v5 = -Vo * n;

        const double iLm_avg_v5  = Io / (n * std::max(1.0 - D, 1e-9));
        const double dILm_pp_v5  = Vin * D * Tsw / Lm;
        const double Im_min_v5   = iLm_avg_v5 - dILm_pp_v5 / 2.0;
        const double Im_max_v5   = iLm_avg_v5 + dILm_pp_v5 / 2.0;
        const bool   isDCM_v5    = (Im_min_v5 < 0.0);

        const int N = 128;
        std::vector<double> time, vPri, iPri, iLm, vSec, iSec, vCo, iCo;
        time.reserve(2 * N + 3);
        for (auto* v : {&vPri, &iPri, &iLm, &vSec, &iSec, &vCo, &iCo})
            v->reserve(2 * N + 3);

        // Interval A (Q1 ON): primary current = iLm, secondary OFF.
        for (int k = 0; k <= N; ++k) {
            const double frac = static_cast<double>(k) / N;
            const double t    = frac * tA;
            const double i_lm = Im_min_v5 + frac * dILm_pp_v5;
            time.push_back(t);
            vPri.push_back(Vpri_pos_v5);
            iLm .push_back(i_lm);
            iPri.push_back(i_lm);
            // Secondary diode OFF: secondary winding voltage swings to
            // +Vin/n (reflected primary), no current.
            vSec.push_back(+Vin / n);
            iSec.push_back(0.0);
            vCo .push_back(Vo);
            iCo .push_back(-Io);            // Co supplies Io alone
        }
        // Discontinuity duplicate at end of A.
        time.push_back(tA);
        vPri.push_back(Vpri_neg_v5);
        iLm .push_back(Im_max_v5);
        iPri.push_back(0.0);                // primary opens (sec takes over)
        vSec.push_back(-Vpri_neg_v5 / n);   // = +Vo (forward-bias)
        iSec.push_back(Im_max_v5 * n);
        vCo .push_back(Vo);
        iCo .push_back(Im_max_v5 * n - Io);

        // Interval C (Q2 ON): primary current = 0, secondary current = iLm·n.
        for (int k = 1; k <= N; ++k) {
            const double frac = static_cast<double>(k) / N;
            const double t    = tA + frac * tC;
            const double i_lm = Im_max_v5 - frac * dILm_pp_v5;
            time.push_back(t);
            vPri.push_back(Vpri_neg_v5);
            iLm .push_back(i_lm);
            iPri.push_back(0.0);
            vSec.push_back(+Vo);
            iSec.push_back(i_lm * n);
            vCo .push_back(Vo);
            iCo .push_back(i_lm * n - Io);
        }

        // ---- Diagnostics for V5 ----
        // Switch RMS (asymmetric): Q1 carries iLm during A; Q2 carries the
        // small reverse-recovery during C (≈0 ideal). Use rms_of for
        // accuracy on Q1; Q2 stays effectively zero in this idealisation.
        const double Iq1_rms_v5 = rms_of(iPri, 0, static_cast<size_t>(N), time)
                                  * std::sqrt(D);
        const double Iq2_rms_v5 = 0.0;

        const double Llk_v5 = (computedLeakageInductance > 0.0)
                              ? computedLeakageInductance : 1e-6;
        const double zvs_v5 = compute_zvs_energy_balance(
            Llk_v5, /*Lm_refl=*/Lm, std::abs(Im_max_v5), mosfetCoss, Vin);

        lastDutyCycle                  = D;
        lastConversionRatio            = Vo / Vin;
        lastDcBlockingCapVoltage       = 0.0;     // no Cb in V5
        lastPrimaryPeakVoltagePositive = Vpri_pos_v5;
        lastPrimaryPeakVoltageNegative = Vpri_neg_v5;
        lastSwitchPeakVoltageQ1        = Vin;
        // Q2 sees Vin + reflected secondary clamp = Vin + Vo·n during Q1 ON.
        lastSwitchPeakVoltageQ2        = Vin + Vo * n;
        lastSwitchRmsCurrentQ1         = Iq1_rms_v5;
        lastSwitchRmsCurrentQ2         = Iq2_rms_v5;
        lastZvsMargin                  = zvs_v5;
        lastResonantTransitionTime     = compute_dead_time(Llk_v5 + Lm, mosfetCoss);
        lastSteadyStateFluxExcursion   = Vin * D * Tsw;   // single-polarity
        lastTransientFluxExcursionEstimate = 0.0;
        lastMagnetizingCurrentRipple   = dILm_pp_v5;
        lastOutputInductorRipple       = 0.0;
        lastOperatingMode              = isDCM_v5 ? 1 : 0;
        lastRectifierType              = 3;       // 3 = AHB_FLYBACK

        OperatingPoint operatingPoint;
        auto wfm = [](const std::vector<double>& d, const std::vector<double>& t) {
            Waveform w;
            w.set_ancillary_label(WaveformLabel::CUSTOM);
            w.set_data(d);
            w.set_time(t);
            return w;
        };
        operatingPoint.get_mutable_excitations_per_winding().push_back(
            complete_excitation(wfm(iPri, time), wfm(vPri, time), fsw, "Primary"));
        operatingPoint.get_mutable_excitations_per_winding().push_back(
            complete_excitation(wfm(iSec, time), wfm(vSec, time), fsw, "Secondary 0"));

        OperatingConditions conditions;
        conditions.set_ambient_temperature(opPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        // V5 has no Lo, no Lo2, no Cb. Push placeholder waveforms (Co only)
        // so get_extra_components_inputs callers stay consistent in size.
        // The Lo-vector remains empty, signalling "no Lo extra component".
        extraCbVoltageWaveforms.push_back(wfm(vCo, time));   // Co stand-in
        extraCbCurrentWaveforms.push_back(wfm(iCo, time));

        return operatingPoint;
    }

    const double n   = turnsRatios[0];
    const double Lm  = magnetizingInductance;
    const double Tsw = 1.0 / fsw;
    const double tA  = D * Tsw;
    const double tC  = (1.0 - D) * Tsw;

    const double Vpri_pos = +(1.0 - D) * Vin;
    const double Vpri_neg = -D * Vin;
    const double Vsec_pos = +(1.0 - D) * Vin / n;   // sec voltage in interval A
    const double Vsec_neg = +D * Vin / n;           // sec voltage in interval C

    // Rectifier-type-dependent per-output-inductor load current.
    //   CT / FB : single output inductor carries the full Io.
    //   CD      : two output inductors, each carrying Io/2.
    const bool isCD = (rect == AhbRectifierType::CURRENT_DOUBLER);
    const double Io_per_inductor = isCD ? Io / 2.0 : Io;

    // ---- Magnetizing current with AHB-asymmetric DC bias ----
    //
    // The Cb steady-state charge-balance forces mean(i_pri) = 0 (no DC
    // through a series capacitor). With reflected load = +Io_per_inductor/n
    // in interval A and -Io_per_inductor/n in interval C, the load
    // contribution to mean(i_pri) is (Io_per_inductor/n)·(2D-1). The
    // magnetizing current must therefore carry a DC offset
    //
    //     mean(i_Lm) = -(Io_per_inductor/n)·(2D-1) = (Io_per_inductor/n)·(1-2D)
    //
    // This is the "AHB DC magnetizing bias" that asymmetric-flux AHB
    // analyses (Imbertson-Mohan TIA 1993; TI SLUP223 fig. 3) call out.
    // Without this offset the analytical primary-current waveform has a
    // spurious DC mean of (Io_per_inductor/n)·(2D-1) and disagrees with
    // any Cb-bearing SPICE model by exactly that amount.
    //
    // The peak-to-peak ramp dILm_pp = (1-D)·Vin·D·Tsw / Lm is symmetric
    // (volt-second balance — A rises by dILm_pp, C falls by the same),
    // so Im0 = mean - dILm_pp/2.
    const double iLm_dc_bias = (Io_per_inductor / n) * (1.0 - 2.0 * D);
    const double dILm_pp = (1.0 - D) * Vin * D * Tsw / Lm;
    const double Im0     = iLm_dc_bias - dILm_pp / 2.0;

    // ---- Output inductor sizing (CCM linear; P5 detects DCM) ----
    // Lo is sized in P6/P9. For raw-call usage we accept it via
    // computedOutputInductance (default 0). When Lo == 0 we fall back
    // to the analytical S2 sizing at 30 % current-ripple ratio so the
    // waveform is well-defined and the user can introspect what was
    // used via lastOutputInductorRipple. This is a diagnostic aid,
    // surfaced explicitly — not a silent default for the topology
    // physics.
    //
    // Rectifier-type effects on Lo:
    //   CT / FB : single output inductor; sees +(Vsec_pos - Vo) in A,
    //             +(Vsec_neg - Vo) in C. Volt-second balance ties Vo
    //             to D for self-consistent designs.
    //   CD      : two output inductors Lo1 + Lo2, each carrying Io/2.
    //             Lo1 charges in A, freewheels in C: vLo1 = Vsec_pos-Vo
    //             in A, -Vo in C. Lo2 is the mirror: -Vo in A,
    //             Vsec_neg-Vo in C. ASYMMETRIC ripple between Lo1 and
    //             Lo2 when D ≠ 0.5 (a known AHB-CD wart).

    double Lo1 = computedOutputInductance;
    if (Lo1 <= 0.0) {
        const double dILo_target = std::max(0.30 * Io_per_inductor, 1e-6);
        Lo1 = compute_lo_min(Vo, D, dILo_target, fsw);
        // Stash the size we actually used so get_extra_components_inputs (P9)
        // can emit it as the Lo design requirement when the caller goes
        // straight from process_operating_points → get_extra_components_inputs
        // without a prior process_design_requirements pass.
        computedOutputInductance = Lo1;
    }
    const double Lo2 = Lo1;   // V3 default: matched inductors

    // Lo1 ripple — drives both CT/FB and CD-Lo1 (charges in A).
    const double dILo1_pp = (Vsec_pos - Vo) * tA / Lo1;
    const double ILo1_min = Io_per_inductor - dILo1_pp / 2.0;
    const double ILo1_max = Io_per_inductor + dILo1_pp / 2.0;

    // Lo2 ripple — CD only. Charges during interval C with (Vsec_neg - Vo).
    // For a self-consistent CD design (Vo = D(1-D)Vin/n), this matches
    // -dILo1_pp in magnitude only when D = 0.5; otherwise asymmetric.
    const double dILo2_pp = isCD ? (Vsec_neg - Vo) * tC / Lo2 : 0.0;
    const double ILo2_min = isCD ? Io_per_inductor - std::abs(dILo2_pp) / 2.0 : 0.0;
    const double ILo2_max = isCD ? Io_per_inductor + std::abs(dILo2_pp) / 2.0 : 0.0;

    const bool isDCM = (ILo1_min < 0.0) || (isCD && ILo2_min < 0.0);

    // ---- Sample buffers ----
    const int N = 128;
    const size_t cap = 2 * N + 3;
    std::vector<double> time, vPri, iPri, iLm, iLo1, iLo2, vLo1, vLo2,
        vSec_a, iSec_a, vSec_b, iSec_b, vSec, iSec, vCb;
    for (auto* v : {&time, &vPri, &iPri, &iLm, &iLo1, &iLo2, &vLo1, &vLo2,
                    &vSec_a, &iSec_a, &vSec_b, &iSec_b, &vSec, &iSec, &vCb})
        v->reserve(cap);

    const double VCb_dc = (1.0 - D) * Vin;

    // emit() pushes one sample to all per-rectifier-type arrays. We
    // dispatch on rect inside the lambda once — keeps the integrators
    // simple and lets the discontinuity duplicate sample work uniformly.
    auto emit = [&](double t, double v_pri_k, double i_lm_k,
                    double i_lo1_k, double i_lo2_k) {
        time.push_back(t);
        vPri.push_back(v_pri_k);
        iLm .push_back(i_lm_k);
        iLo1.push_back(i_lo1_k);
        iLo2.push_back(i_lo2_k);
        vCb .push_back(VCb_dc);

        const bool inA = (v_pri_k > 0.0);
        const bool inC = (v_pri_k < 0.0);

        // Primary current: Lm component + reflected load. CT/FB always
        // reflect the single Lo current; CD reflects Lo1 in A and Lo2
        // in C (each conducts on the opposite half-cycle).
        if (inA)
            iPri.push_back(i_lo1_k / n + i_lm_k);
        else if (inC)
            iPri.push_back(-(isCD ? i_lo2_k : i_lo1_k) / n + i_lm_k);
        else
            iPri.push_back(i_lm_k);

        // Secondary winding emission depends on rectifier type.
        switch (rect) {
            case AhbRectifierType::CENTER_TAPPED: {
                if (inA) {
                    vSec_a.push_back(+Vsec_pos); iSec_a.push_back(i_lo1_k);
                    vSec_b.push_back(-Vsec_pos); iSec_b.push_back(0.0);
                } else if (inC) {
                    vSec_a.push_back(-Vsec_neg); iSec_a.push_back(0.0);
                    vSec_b.push_back(+Vsec_neg); iSec_b.push_back(i_lo1_k);
                } else {
                    vSec_a.push_back(0.0); iSec_a.push_back(0.0);
                    vSec_b.push_back(0.0); iSec_b.push_back(0.0);
                }
                vSec.push_back(0.0); iSec.push_back(0.0);   // unused
                break;
            }
            case AhbRectifierType::FULL_BRIDGE: {
                // Single secondary winding, bipolar.
                if (inA)      { vSec.push_back(+Vsec_pos); iSec.push_back(+i_lo1_k); }
                else if (inC) { vSec.push_back(-Vsec_neg); iSec.push_back(-i_lo1_k); }
                else          { vSec.push_back(0.0);       iSec.push_back(0.0); }
                vSec_a.push_back(0.0); iSec_a.push_back(0.0);   // unused
                vSec_b.push_back(0.0); iSec_b.push_back(0.0);   // unused
                break;
            }
            case AhbRectifierType::CURRENT_DOUBLER: {
                // Single secondary winding, bipolar. In A the upper
                // diode forwards i_Lo1; in C the lower diode forwards
                // i_Lo2 with reversed winding polarity.
                if (inA)      { vSec.push_back(+Vsec_pos); iSec.push_back(+i_lo1_k); }
                else if (inC) { vSec.push_back(-Vsec_neg); iSec.push_back(-i_lo2_k); }
                else          { vSec.push_back(0.0);       iSec.push_back(0.0); }
                vSec_a.push_back(0.0); iSec_a.push_back(0.0);   // unused
                vSec_b.push_back(0.0); iSec_b.push_back(0.0);   // unused
                break;
            }
            default: break;   // AHB_FLYBACK already rejected at entry
        }

        // Lo1 voltage:
        //   CT / FB : sees rectified secondary minus Vo on both intervals
        //   CD      : charges (Vsec_pos - Vo) in A; freewheels (-Vo) in C
        if (rect == AhbRectifierType::CURRENT_DOUBLER) {
            vLo1.push_back(inA ? (Vsec_pos - Vo) : (inC ? -Vo : -Vo));
            vLo2.push_back(inA ? -Vo : (inC ? (Vsec_neg - Vo) : -Vo));
        } else {
            vLo1.push_back(inA ? (Vsec_pos - Vo)
                                : (inC ? (Vsec_neg - Vo) : -Vo));
            vLo2.push_back(0.0);   // unused
        }
    };

    // Interval A
    for (int k = 0; k <= N; ++k) {
        const double frac = static_cast<double>(k) / N;
        const double t    = frac * tA;
        const double i_lm = Im0 + frac * dILm_pp;
        // Lo1 charges (+slope) in A; Lo2 freewheels (-slope) in A (CD only).
        const double i_lo1 = ILo1_min + frac * dILo1_pp;
        const double i_lo2 = isCD ? (ILo2_max + frac * (-Vo / Lo2) * tA) : 0.0;
        emit(t, Vpri_pos, i_lm, i_lo1, i_lo2);
    }

    // Discontinuity duplicate sample (post-jump v_pri).
    {
        // End-of-A inductor values: Lo1 = ILo1_max; for CD, Lo2 trough.
        const double iLo2_end_A = isCD
            ? (ILo2_max + (-Vo / Lo2) * tA)   // Lo2 dropped during A
            : 0.0;
        emit(tA, Vpri_neg, Im0 + dILm_pp, ILo1_max, iLo2_end_A);
    }

    // Interval C
    for (int k = 1; k <= N; ++k) {
        const double frac = static_cast<double>(k) / N;
        const double t    = tA + frac * tC;
        const double i_lm = (Im0 + dILm_pp) - frac * dILm_pp;   // ramp back to Im0
        // Lo1 discharges (-slope) in C for CT/FB; for CD Lo1 freewheels (-Vo/Lo1).
        // Lo2 charges (+slope) in C (CD only).
        double i_lo1, i_lo2 = 0.0;
        if (rect == AhbRectifierType::CURRENT_DOUBLER) {
            i_lo1 = ILo1_max + frac * (-Vo / Lo1) * tC;
            const double iLo2_start_C = ILo2_max + (-Vo / Lo2) * tA;
            i_lo2 = iLo2_start_C + frac * dILo2_pp;
        } else {
            i_lo1 = ILo1_max - frac * dILo1_pp;
        }
        emit(t, Vpri_neg, i_lm, i_lo1, i_lo2);
    }

    // ---- Asymmetric switch RMS (Imbertson-Mohan §IV.A) ----
    const size_t idx_tA  = static_cast<size_t>(N);
    const size_t idx_end = time.size() - 1;
    const double Iq1_segrms = rms_of(iPri, 0,           idx_tA, time);
    const double Iq2_segrms = rms_of(iPri, idx_tA + 1,  idx_end, time);
    const double Iq1_rms = Iq1_segrms * std::sqrt(D);
    const double Iq2_rms = Iq2_segrms * std::sqrt(1.0 - D);

    // ---- ZVS energy margin & dead-time ----
    //
    // P5 refinement: in the AHB the secondary diodes BLOCK during the
    // dead-time (rectifier output filter is supplied by Lo from the
    // capacitor side), so the resonant tank seen at the switch node is
    // Llk + Lm (NOT just Llk). The magnetizing current also contributes
    // ZVS energy because Lm is in series with Llk during the transition:
    //
    //   E_stored = 0.5 * (Llk + Lm_refl) * I_pri^2
    //   t_dead   = pi * sqrt((Llk + Lm) * 2 * Coss)
    //
    // This matches Mappus 2014 (TI SLUP223 §5) and Imbertson-Mohan §V.
    // PSFB by contrast keeps Lm out of the tank because the secondary
    // freewheels through both rectifier diodes during dead-time, clamping
    // the primary at zero — that's where the AHB / PSFB ZVS physics
    // diverge.
    const double Ipri_at_Q1_off = iPri[idx_tA];
    const double Llk = (computedLeakageInductance > 0.0)
                       ? computedLeakageInductance : 1e-6;
    const double zvsMargin = compute_zvs_energy_balance(
        Llk, /*Lm_refl=*/Lm, std::abs(Ipri_at_Q1_off), mosfetCoss, Vin);
    const double tdRes = compute_dead_time(Llk + Lm, mosfetCoss);

    // ---- Transient flux excursion (P5) ----
    //
    // Reported in volt-seconds (matches lastSteadyStateFluxExcursion's
    // units); a downstream magnetic adviser divides by Np·Ae to get B.
    // Worst-case bound: V_Cb cannot follow a Vin step within one cycle,
    // so the volt-second adder per half-cycle is D·dVin·Tsw on top of
    // the steady-state D(1-D)·Vin·Tsw. Throws if dVin < 0.
    double transientFluxExcursion = 0.0;
    if (auto dVinOpt = get_input_voltage_step_range()) {
        const double dVin = dVinOpt.value();
        if (dVin < 0.0)
            throw std::invalid_argument(
                "AsymmetricHalfBridge::process_operating_point_for_input_voltage: "
                "inputVoltageStepRange must be >= 0");
        transientFluxExcursion =
            D * (1.0 - D) * Vin * Tsw + D * dVin * Tsw;
    }

    // ---- Diagnostics ----
    lastDutyCycle                  = D;
    lastConversionRatio            = Vo / Vin;
    lastDcBlockingCapVoltage       = VCb_dc;
    lastPrimaryPeakVoltagePositive = Vpri_pos;
    lastPrimaryPeakVoltageNegative = Vpri_neg;
    lastSwitchPeakVoltageQ1        = Vin;
    lastSwitchPeakVoltageQ2        = Vin;
    lastSwitchRmsCurrentQ1         = Iq1_rms;
    lastSwitchRmsCurrentQ2         = Iq2_rms;
    lastZvsMargin                  = zvsMargin;
    lastResonantTransitionTime     = tdRes;
    lastSteadyStateFluxExcursion   = D * (1.0 - D) * Vin * Tsw;
    lastTransientFluxExcursionEstimate = transientFluxExcursion;
    lastMagnetizingCurrentRipple   = dILm_pp;
    lastOutputInductorRipple       = dILo1_pp;
    lastOperatingMode              = isDCM ? 1 : 0;
    switch (rect) {
        case AhbRectifierType::CENTER_TAPPED:   lastRectifierType = 0; break;
        case AhbRectifierType::CURRENT_DOUBLER: lastRectifierType = 1; break;
        case AhbRectifierType::FULL_BRIDGE:     lastRectifierType = 2; break;
        case AhbRectifierType::AHB_FLYBACK:     lastRectifierType = 3; break;
    }

    // ---- Build OperatingPoint excitations ----
    OperatingPoint operatingPoint;

    auto wfm = [](const std::vector<double>& d, const std::vector<double>& t) {
        Waveform w;
        w.set_ancillary_label(WaveformLabel::CUSTOM);
        w.set_data(d);
        w.set_time(t);
        return w;
    };

    operatingPoint.get_mutable_excitations_per_winding().push_back(
        complete_excitation(wfm(iPri, time), wfm(vPri, time), fsw, "Primary"));

    if (rect == AhbRectifierType::CENTER_TAPPED) {
        operatingPoint.get_mutable_excitations_per_winding().push_back(
            complete_excitation(wfm(iSec_a, time), wfm(vSec_a, time), fsw,
                                "Secondary 0a"));
        operatingPoint.get_mutable_excitations_per_winding().push_back(
            complete_excitation(wfm(iSec_b, time), wfm(vSec_b, time), fsw,
                                "Secondary 0b"));
    } else {
        // FB and CD: single secondary winding.
        operatingPoint.get_mutable_excitations_per_winding().push_back(
            complete_excitation(wfm(iSec, time), wfm(vSec, time), fsw,
                                "Secondary 0"));
    }

    // ---- V6 multi-output: per-secondary load-share projection ----
    // The primary current already reflects the TOTAL load (Io was overridden
    // to the power-equivalent at Vo_0 above). The standard secondary
    // excitations above carry currents scaled to that total. Replace them
    // with per-output excitations using the load-share formula:
    //     share_k = Vo_k·Io_k / Σ Vo_j·Io_j
    //     i_sec_k = share_k · n_k · iPri    (ampere-turn balance)
    // Per-output voltage = Vo_k square wave with same A/C polarity as
    // primary reflected. This matches the Dab.cpp multi-output pattern.
    if (outV.size() > 1 && totalPower > 0.0) {
        // Pop the standard secondary excitations (CT pops 2, FB/CD pops 1).
        auto& exc = operatingPoint.get_mutable_excitations_per_winding();
        const size_t numToPop =
            (rect == AhbRectifierType::CENTER_TAPPED) ? 2u : 1u;
        for (size_t i = 0; i < numToPop && !exc.empty(); ++i) exc.pop_back();

        for (size_t k = 0; k < outV.size() && k < turnsRatios.size(); ++k) {
            const double Vo_k = outV[k];
            const double Io_k = (k < outI.size()) ? outI[k] : 0.0;
            const double n_k  = turnsRatios[k];
            if (n_k <= 0.0) continue;
            const double share = (Vo_k * Io_k) / totalPower;
            std::vector<double> iSec_k(iPri.size());
            std::vector<double> vSec_k(vPri.size());
            for (size_t i = 0; i < iPri.size(); ++i) {
                // Sign convention matches "Secondary 0" (positive when
                // diode forward-biased during interval A).
                iSec_k[i] = share * n_k * iPri[i];
                // Secondary voltage: scale primary voltage by 1/n_k and
                // re-anchor to ±Vo_k (during diode-OFF, sec voltage = -Vo_k;
                // diode-ON, +Vo_k). For the analytical PWL waveform here
                // we approximate vSec_k_during_A = +Vo_k, _during_C = -Vo_k.
                vSec_k[i] = (vPri[i] > 0.0) ? +Vo_k : -Vo_k;
            }
            const std::string label = "Secondary " + std::to_string(k);
            exc.push_back(complete_excitation(
                wfm(iSec_k, time), wfm(vSec_k, time), fsw, label));
        }
    }

    OperatingConditions conditions;
    conditions.set_ambient_temperature(opPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    // ---- Extra-component waveforms (Lo, [Lo2], Cb) ----
    extraLoVoltageWaveforms.push_back(wfm(vLo1, time));
    extraLoCurrentWaveforms.push_back(wfm(iLo1, time));
    if (isCD) {
        extraLo2VoltageWaveforms.push_back(wfm(vLo2, time));
        extraLo2CurrentWaveforms.push_back(wfm(iLo2, time));
    }
    extraCbVoltageWaveforms.push_back(wfm(vCb, time));
    extraCbCurrentWaveforms.push_back(wfm(iPri, time));   // Cb in series with primary

    return operatingPoint;
}


// =========================================================================
// process_operating_points (P3 deliverable): loop over (Vin × OPs).
//
// Mirrors the PSHB / ACF pattern: collect distinct input voltages
// (nominal/min/max) and call the per-input-voltage routine for every
// operating point, naming each entry "<volt-tag> input volt." (with the
// OP index appended when there is more than one). Extra-component
// waveform vectors are cleared at the top so a re-run leaves the class
// in a deterministic state.
// =========================================================================
std::vector<OperatingPoint> AsymmetricHalfBridge::process_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    extraLoVoltageWaveforms.clear();
    extraLoCurrentWaveforms.clear();
    extraLo2VoltageWaveforms.clear();
    extraLo2CurrentWaveforms.clear();
    extraCbVoltageWaveforms.clear();
    extraCbCurrentWaveforms.clear();

    std::vector<double> inputVoltages;
    std::vector<std::string> inputVoltageNames;
    collect_input_voltages(get_input_voltage(), inputVoltages,
                           inputVoltageNames);
    if (inputVoltages.empty())
        throw std::runtime_error(
            "AsymmetricHalfBridge::process_operating_points: inputVoltage "
            "must specify nominal, minimum, or maximum");

    auto& ops = get_operating_points();
    if (ops.empty())
        throw std::runtime_error(
            "AsymmetricHalfBridge::process_operating_points: no operating points");

    std::vector<OperatingPoint> result;
    result.reserve(inputVoltages.size() * ops.size());
    for (size_t vi = 0; vi < inputVoltages.size(); ++vi) {
        for (size_t oi = 0; oi < ops.size(); ++oi) {
            auto op = process_operating_point_for_input_voltage(
                inputVoltages[vi], ops[oi], turnsRatios, magnetizingInductance);
            std::string name = inputVoltageNames[vi] + " input volt.";
            if (ops.size() > 1)
                name += " with op. point " + std::to_string(oi);
            op.set_name(name);
            result.push_back(std::move(op));
        }
    }
    return result;
}

std::vector<OperatingPoint> AsymmetricHalfBridge::process_operating_points(
    Magnetic magnetic) {
    AsymmetricHalfBridge::run_checks(_assertErrors);
    auto& settings = Settings::GetInstance();
    OpenMagnetics::MagnetizingInductance magnetizingInductanceModel(
        settings.get_reluctance_model());
    double magnetizingInductance =
        magnetizingInductanceModel
            .calculate_inductance_from_number_turns_and_gapping(
                magnetic.get_mutable_core(), magnetic.get_mutable_coil())
            .get_magnetizing_inductance().get_nominal().value();
    std::vector<double> turnsRatios = magnetic.get_turns_ratios();
    return process_operating_points(turnsRatios, magnetizingInductance);
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

// =========================================================================
// Extra components inputs (P9 deliverable)
//
// Mirrors PhaseShiftedHalfBridge::get_extra_components_inputs (PSHB sibling).
// AHB has two ancillary *magnetic* components in the basic V1/V2 designs:
//
//   - Output inductor Lo                   (always present except V5 flyback)
//   - Second output inductor Lo2           (current-doubler V3 only)
//
// The DC-blocking capacitor Cb and output capacitor Co are non-magnetic
// components: their sized values are exposed through the
// `get_computed_dc_blocking_capacitance()` and
// `get_computed_output_capacitance()` accessors and through the per-OP
// extraCb*Waveforms vectors. They are NOT returned here because the
// `Inputs` schema is the magnetics-design contract; capacitor selection
// is handled by the consumer using the accessors.
//
// The IDEAL/REAL distinction follows the PSHB convention: in REAL mode
// the caller supplies the designed primary `Magnetic` so the topology
// can subtract any leakage already absorbed by the transformer from
// auxiliary inductors that physically live in series with the primary.
// AHB has no such series inductor (Lo lives on the secondary, isolated
// from primary leakage), so REAL == IDEAL for AHB. We keep the magnetic
// argument required in REAL mode for API parity with PSHB.
// =========================================================================
std::vector<std::variant<Inputs, CAS::Inputs>>
AsymmetricHalfBridge::get_extra_components_inputs(
    ExtraComponentsMode mode, std::optional<Magnetic> magnetic)
{
    if (mode == ExtraComponentsMode::REAL && !magnetic.has_value())
        throw std::invalid_argument(
            "AsymmetricHalfBridge::get_extra_components_inputs: "
            "mode REAL requires a designed magnetic");

    if (computedOutputInductance <= 0 || extraLoVoltageWaveforms.empty())
        throw std::runtime_error(
            "AsymmetricHalfBridge::get_extra_components_inputs: "
            "call process_operating_points() first");

    if (extraLoVoltageWaveforms.size() != extraLoCurrentWaveforms.size())
        throw std::runtime_error(
            "AsymmetricHalfBridge::get_extra_components_inputs: "
            "Lo voltage/current waveform count mismatch");

    std::vector<std::variant<Inputs, CAS::Inputs>> result;
    const size_t nOps = extraLoVoltageWaveforms.size();
    const double fsw  = get_operating_points()[0].get_switching_frequency();

    // ---- Output inductor Lo (V1 / V2 / V3 / V4) ----
    {
        Inputs masInputs;
        DesignRequirements dr;
        DimensionWithTolerance loInductance;
        loInductance.set_nominal(computedOutputInductance);
        dr.set_magnetizing_inductance(loInductance);
        dr.set_name("outputInductor");
        dr.set_topology(Topologies::ASYMMETRIC_HALF_BRIDGE_CONVERTER);
        dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
        dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
        masInputs.set_design_requirements(dr);

        std::vector<OperatingPoint> masOps;
        masOps.reserve(nOps);
        for (size_t i = 0; i < nOps; ++i) {
            OperatingPoint op;
            op.get_mutable_excitations_per_winding().push_back(
                complete_excitation(extraLoCurrentWaveforms[i],
                                    extraLoVoltageWaveforms[i],
                                    fsw, "Primary"));
            masOps.push_back(op);
        }
        masInputs.set_operating_points(masOps);
        result.emplace_back(std::move(masInputs));
    }

    // ---- Second output inductor Lo2 (current-doubler V3 only) ----
    if (!extraLo2VoltageWaveforms.empty()) {
        if (extraLo2VoltageWaveforms.size() != extraLo2CurrentWaveforms.size())
            throw std::runtime_error(
                "AsymmetricHalfBridge::get_extra_components_inputs: "
                "Lo2 voltage/current waveform count mismatch");
        Inputs masInputs;
        DesignRequirements dr;
        DimensionWithTolerance lo2Inductance;
        lo2Inductance.set_nominal(computedOutputInductance);  // CD: Lo1 == Lo2
        dr.set_magnetizing_inductance(lo2Inductance);
        dr.set_name("outputInductor2");
        dr.set_topology(Topologies::ASYMMETRIC_HALF_BRIDGE_CONVERTER);
        dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
        dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
        masInputs.set_design_requirements(dr);

        std::vector<OperatingPoint> masOps;
        masOps.reserve(extraLo2VoltageWaveforms.size());
        for (size_t i = 0; i < extraLo2VoltageWaveforms.size(); ++i) {
            OperatingPoint op;
            op.get_mutable_excitations_per_winding().push_back(
                complete_excitation(extraLo2CurrentWaveforms[i],
                                    extraLo2VoltageWaveforms[i],
                                    fsw, "Primary"));
            masOps.push_back(op);
        }
        masInputs.set_operating_points(masOps);
        result.emplace_back(std::move(masInputs));
    }

    return result;
}

// =========================================================================
// generate_ngspice_circuit (P7 deliverable)
//
// Topology (Imbertson-Mohan upper-Cb convention, matches the analytical
// model's V_Cb = (1-D)·Vin):
//
//      Vin ──┬───────┬────[Cb]────pri_top──[Llk]──[Lpri]── sw  ──┐
//            │       │                        ║                  │
//          [Q1]    [Cb-snub]                  ║ K=0.99            │
//            │                                ║                 [Q2]
//            sw ─────────────────             ║                   │
//                                          [Lsec]──>>(rectifier)─0
//
// Conventions matched to the P3-P5 analytical model:
//   - V_Cb,dc = (1-D)·Vin (Imbertson-Mohan original).
//   - Primary winding dot at the sw side, so v_pri = v(sw) − v(pri_top)
//     swings to +(1-D)·Vin when Q1 ON and to −D·Vin when Q2 ON.
//   - Vab probe is the differential v_pri = v(sw) − v(pri_top).
//   - Vpri_sense (zero V-source) measures primary current with the same
//     sign as the analytical i_pri (positive when current flows from sw
//     into the primary winding, i.e. during Q1 ON at heavy load).
//
// Per AHB plan §12 mistakes #10 / #11 / #12 the netlist commits to:
//   - K = 0.99 (not 1.0): factorisation stable, leakage explicit.
//   - METHOD=GEAR + TRTOL=7: ngspice convergence on switching circuits.
//   - IC pre-charge for ALL FOUR reactives:
//        L_pri (i_Lm,0 = -dILm/2),  Cb ((1-D)·Vin),  Lo (Io_per_inductor),
//        Co (Vo). Without this the first cycle saturates the transformer.
//   - Per-switch RC snubbers (1k / 1n) for ngspice convergence.
//   - Parasitic DCR on Lo (50 mΩ) and ESR on Cb / Co (5 mΩ) — lossless
//     LC networks are badly conditioned in ngspice.
//   - Rectifier-type branch:
//        CT  : two half-windings (sec_a, sec_b), ct-tap, two diodes, one Lo.
//        FB  : single secondary, four-diode bridge, one Lo.
//        CD  : single secondary, two diodes, two Los (Lo1 + Lo2).
//   - AHB_FLYBACK : deferred to P10 (no Lo; Lm itself is the storage).
//
// The complementary PWM uses a 50 ns dead-time so the switch transition
// is a real (resonant-or-hard) event ngspice can simulate, instead of
// instantaneous switching that defeats the IC.
// =========================================================================
std::string AsymmetricHalfBridge::generate_ngspice_circuit(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t inputVoltageIndex,
    size_t operatingPointIndex)
{
    // ---- Validation (Guide §5: throw loudly on malformed inputs) ----
    AhbRectifierType rect = get_rectifier_type().value_or(
        AhbRectifierType::CENTER_TAPPED);
    if (turnsRatios.empty())
        throw std::invalid_argument(
            "AsymmetricHalfBridge::generate_ngspice_circuit: turnsRatios empty");
    if (!(turnsRatios[0] > 0.0))
        throw std::invalid_argument(
            "AsymmetricHalfBridge::generate_ngspice_circuit: turnsRatios[0] must be > 0");
    if (!(magnetizingInductance > 0.0))
        throw std::invalid_argument(
            "AsymmetricHalfBridge::generate_ngspice_circuit: "
            "magnetizingInductance must be > 0");

    std::vector<double> inputVoltages;
    std::vector<std::string> inputVoltageNames;
    collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltageNames);
    if (inputVoltages.empty())
        throw std::runtime_error(
            "AsymmetricHalfBridge::generate_ngspice_circuit: inputVoltage "
            "must specify nominal, minimum, or maximum");
    const double Vin =
        inputVoltages[std::min(inputVoltageIndex, inputVoltages.size() - 1)];

    auto& ops = get_operating_points();
    if (ops.empty())
        throw std::runtime_error(
            "AsymmetricHalfBridge::generate_ngspice_circuit: at least one "
            "operating point required");
    const auto& op = ops[std::min(operatingPointIndex, ops.size() - 1)];

    // V6 multi-output netlist policy: collapse all outputs into a single
    // power-equivalent load referenced to output #0 (the controlling rail).
    // Per-secondary leakage / cross-regulation is NOT modelled in SPICE —
    // the analytical model already emits per-output secondary excitations
    // via load-share projection. Warn once so the user knows.
    const auto& outV_n = op.get_output_voltages();
    const auto& outI_n = op.get_output_currents();
    if (outV_n.empty() || outI_n.empty() ||
        outV_n.size() != outI_n.size())
        throw std::runtime_error(
            "AsymmetricHalfBridge::generate_ngspice_circuit: outputVoltages "
            "/ outputCurrents must be non-empty and equal length");
    if (outV_n.size() > 1) {
        static thread_local bool ahbMultiOutNetlistWarned = false;
        if (!ahbMultiOutNetlistWarned) {
            std::cerr << "[AHB] generate_ngspice_circuit: multi-output "
                         "configuration collapsed to single power-equivalent "
                         "load on output #0 (per-secondary cross-regulation "
                         "is not simulated in SPICE; use the analytical "
                         "process_operating_points for per-output waveforms)."
                      << std::endl;
            ahbMultiOutNetlistWarned = true;
        }
    }
    double totalPower_n = 0.0;
    for (size_t k = 0; k < outV_n.size(); ++k)
        totalPower_n += outV_n[k] * outI_n[k];
    const double Vo  = outV_n[0];
    const double Io  = (outV_n.size() > 1 && Vo > 0.0)
                       ? totalPower_n / Vo
                       : outI_n[0];
    const double fsw = op.get_switching_frequency();
    const double D   = op.get_duty_cycle();
    if (!(D > 0.0 && D < 1.0))
        throw std::runtime_error("AHB SPICE: dutyCycle must lie in (0, 1)");
    if (!(fsw > 0.0))
        throw std::runtime_error("AHB SPICE: switchingFrequency must be > 0");
    if (!(Vo > 0.0))
        throw std::runtime_error("AHB SPICE: outputVoltage must be > 0");
    if (!(Io > 0.0))
        throw std::runtime_error("AHB SPICE: outputCurrent must be > 0");

    const double n   = turnsRatios[0];
    const double Lm  = magnetizingInductance;

    // ====================================================================
    // V5 — AHB_FLYBACK (active-clamp / two-switch flyback) netlist branch.
    // Topology: Q1 main switch (high-side, gate g1) energises L_pri (sw→0)
    // during D·Tsw. During (1-D)·Tsw the primary energy transfers to the
    // secondary via the coupled flyback diode. Q2 + C_clamp form the
    // active-clamp branch parallel to L_pri (catches leakage ringing and
    // recycles energy; clamp cap charges to ≈ Vo·n in steady state).
    // No series Cb, no Lo. Vo = D·Vin / ((1-D)·n).  Mirrors V5 analytical
    // model in process_operating_point_for_input_voltage.
    // ====================================================================
    if (rect == AhbRectifierType::AHB_FLYBACK) {
        const double Tsw_v5 = 1.0 / fsw;
        const double tA_v5  = D * Tsw_v5;
        constexpr double T_DEAD_V5 = 50e-9;
        const double t_q1_on_v5 = std::max(tA_v5 - T_DEAD_V5, 1e-9);
        const double t_q2_on_v5 = std::max((1.0 - D) * Tsw_v5 - T_DEAD_V5, 1e-9);
        const double dILm_pp_v5 = Vin * D * Tsw_v5 / Lm;
        const double iLm_avg_v5 = Io / (n * std::max(1.0 - D, 1e-9));
        const double iLm_init_v5 = iLm_avg_v5 - dILm_pp_v5 / 2.0;
        const double Vclamp_dc_v5 = Vo * n;
        const double C_clamp_v5 = std::max(
            (iLm_avg_v5 + dILm_pp_v5 / 2.0) * D * Tsw_v5
                / std::max(0.05 * Vclamp_dc_v5, 1e-3),
            1e-9);
        double Co_v5 = computedOutputCapacitance;
        if (Co_v5 <= 0.0)
            Co_v5 = Io * D / (fsw * std::max(0.01 * Vo, 1e-3));
        const double Rload_v5 = std::max(Vo / Io, 1e-3);

        const int periodsToExtract = numPeriodsToExtract;
        const int steadyStatePeriods = numSteadyStatePeriods;
        const int numPeriodsTotal = steadyStatePeriods + periodsToExtract;
        const double simTime_v5   = numPeriodsTotal * Tsw_v5;
        const double startTime_v5 = steadyStatePeriods * Tsw_v5;
        const double stepTime_v5  = Tsw_v5 / 200.0;

        std::ostringstream c;
        c << std::scientific;
        c << "* AHB V5 — AHB_FLYBACK (active-clamp / two-switch flyback)\n";
        c << "* Vin=" << Vin << " V, Vo=" << Vo << " V, Io=" << Io
          << " A, fsw=" << (fsw / 1e3) << " kHz, D=" << D << ", n=" << n << "\n";
        c << "* Lm=" << (Lm * 1e6) << " uH, Cclamp=" << (C_clamp_v5 * 1e6)
          << " uF, Co=" << (Co_v5 * 1e6) << " uF\n";
        c << "* IC: i_Lm=" << iLm_init_v5 << " A, V_Cclamp=" << Vclamp_dc_v5
          << " V, V_Co=" << Vo << " V\n\n";

        c << ".model SW1 SW VT=2.5 VH=0.8 RON=1m ROFF=1Meg\n";
        c << ".model DIDEAL D(IS=1e-12 RS=0.001 BV=1000 IBV=1e-12)\n\n";

        c << "Vdc vin_dc 0 " << Vin << "\n\n";

        c << "Vpwm_Q1 g1 0 PULSE(0 5 0 10n 10n " << t_q1_on_v5 << " " << Tsw_v5 << ")\n";
        c << "Vpwm_Q2 g2 0 PULSE(0 5 " << (tA_v5 + T_DEAD_V5) << " 10n 10n "
          << t_q2_on_v5 << " " << Tsw_v5 << ")\n\n";

        c << "S1 vin_dc sw g1 0 SW1\n";
        c << "D1 sw vin_dc DIDEAL\n";
        c << "Rsnub_Q1 vin_dc sw 1k\nCsnub_Q1 vin_dc sw 1n\n";
        c << "S2 sw clamp_lo g2 0 SW1\n";
        c << "D2 clamp_lo sw DIDEAL\n";
        c << "Rsnub_Q2 sw clamp_lo 1k\nCsnub_Q2 sw clamp_lo 1n\n";
        c << "C_clamp clamp_lo 0 " << C_clamp_v5 << " IC=" << Vclamp_dc_v5 << "\n";
        c << "R_clamp_bleed clamp_lo 0 1Meg\n\n";

        c << "Vpri_sense sw pri_top 0\n";
        c << "L_pri pri_top 0 " << Lm << " IC=" << iLm_init_v5 << "\n";
        c << "Evab vab 0 sw 0 1\n\n";

        constexpr double K_COUPLING_V5 = 0.999;
        const double Lsec_v5 = Lm / (n * n);
        c << "L_sec sec_a sec_b " << Lsec_v5 << "\n";
        c << "K1 L_pri L_sec " << K_COUPLING_V5 << "\n";
        c << "Vsec_a_sense sec_b sec_d 0\n";
        c << "D_fly sec_d co_top DIDEAL\n";
        c << "R_sec_anchor sec_a 0 1Meg\n\n";

        c << "R_co_esr co_top co_top_esr 1m\n";
        c << "C_o co_top_esr out_gnd " << Co_v5 << " IC=" << Vo << "\n";
        c << "R_load co_top out_gnd " << Rload_v5 << "\n";
        c << "Vout_sense out_gnd 0 0\n\n";

        c << ".tran " << stepTime_v5 << " " << simTime_v5 << " "
          << startTime_v5 << " uic\n\n";
        c << ".save v(vab) v(sw) v(co_top) i(Vpri_sense) i(Vsec_a_sense)"
          << " i(Vout_sense) i(Vdc)\n\n";
        c << ".options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500\n";
        c << ".options METHOD=GEAR TRTOL=7\n";
        c << ".nodeset v(co_top)=" << Vo
          << " v(clamp_lo)=" << Vclamp_dc_v5 << "\n\n";
        c << ".end\n";
        return c.str();
    }

    // Explicit L_lk is set to a tiny "wire-stub" value (1 nH). The
    // dominant primary leakage in the SPICE model comes from the K=0.999
    // coupling factor between primary and secondaries (Llk_eff ≈ Lm·(1-K²)
    // ≈ 2 µH for Lm ≈ 1 mH). Using a large explicit L_lk on top of that
    // would double-count leakage and slow primary-current transitions
    // beyond the analytical PWL idealisation that P3-P5 use, blowing up
    // the analytical-vs-SPICE NRMSE comparison in P8.
    const double Llk = (computedLeakageInductance > 0.0)
                       ? computedLeakageInductance : 1e-9;
    const double Tsw = 1.0 / fsw;
    const double tA  = D * Tsw;

    // Lo / Cb / Co: prefer pre-sized values; fall back to analytical
    // 30 % ripple / 5 % V_Cb / 1 % Vo defaults so the netlist is always
    // emittable even when process_design_requirements wasn't called.
    const bool isCD = (rect == AhbRectifierType::CURRENT_DOUBLER);
    const double Io_per_inductor = isCD ? Io / 2.0 : Io;
    double Lo = computedOutputInductance;
    if (Lo <= 0.0)
        Lo = compute_lo_min(Vo, D, std::max(0.30 * Io_per_inductor, 1e-6), fsw);
    double Cb = computedDcBlockingCapacitance;
    if (Cb <= 0.0) {
        const double VCb = (1.0 - D) * Vin;
        const double dILm_pp_max = (1.0 - D) * Vin * D * Tsw / Lm;
        const double Iprimary_pk = Io / n + dILm_pp_max / 2.0;
        Cb = compute_cb_min(Iprimary_pk, D,
                            std::max(0.05 * VCb, 1e-3), fsw);
    }
    double Co = computedOutputCapacitance;
    if (Co <= 0.0) {
        const double dILo = std::max(0.30 * Io_per_inductor, 1e-6);
        const double dILo_total = isCD ? 2.0 * dILo : dILo;
        Co = dILo_total / (8.0 * fsw * std::max(0.01 * Vo, 1e-3));
    }

    // IC values for all four reactive elements (mistake #10).
    const double VCb_dc      = (1.0 - D) * Vin;
    const double dILm_pp     = (1.0 - D) * Vin * D * Tsw / Lm;
    const double iLm_dc_bias = (Io_per_inductor / n) * (1.0 - 2.0 * D);
    const double iLm_init    = iLm_dc_bias - dILm_pp / 2.0;  // start of period A,
                                                              // matches analytical Im0
    const double iLo_init    = Io_per_inductor;       // mean Lo current
    const double Rload       = std::max(Vo / Io, 1e-3);

    // PWM timing: 50 ns dead-time (resonant transition envelope).
    constexpr double T_DEAD = 50e-9;
    const double t_q1_on = std::max(tA - T_DEAD, 1e-9);
    const double t_q2_on = std::max((1.0 - D) * Tsw - T_DEAD, 1e-9);

    // Simulation horizon: numSteadyStatePeriods to settle, then
    // numPeriodsToExtract to capture.
    const int periodsToExtract = numPeriodsToExtract;
    const int steadyStatePeriods = numSteadyStatePeriods;
    const int numPeriodsTotal = steadyStatePeriods + periodsToExtract;
    const double simTime   = numPeriodsTotal * Tsw;
    const double startTime = steadyStatePeriods * Tsw;
    const double stepTime  = Tsw / 200.0;

    std::ostringstream c;
    c << std::scientific;

    c << "* Asymmetric Half Bridge (AHB) — Imbertson-Mohan / TI SLUP223 / ON-Semi AN-4153\n";
    c << "* Vin=" << Vin << " V, Vo=" << Vo << " V, Io=" << Io
      << " A, fsw=" << (fsw/1e3) << " kHz, D=" << D
      << ", n=" << n << ", rect=" << static_cast<int>(rect) << "\n";
    c << "* Lm=" << (Lm*1e6) << " uH, Llk=" << (Llk*1e6) << " uH, Lo="
      << (Lo*1e6) << " uH, Cb=" << (Cb*1e6) << " uF, Co=" << (Co*1e6) << " uF\n";
    c << "* IC: V_Cb=" << VCb_dc << " V, i_Lm=" << iLm_init
      << " A, i_Lo=" << iLo_init << " A, V_Co=" << Vo << " V\n\n";

    c << ".model SW1 SW VT=2.5 VH=0.8 RON=1m ROFF=1Meg\n";
    // Near-ideal diode RS so SPICE secondary-side operating point stays
    // close to the analytical ideal-diode assumption (P8 NRMSE).
    c << ".model DIDEAL D(IS=1e-12 RS=0.001 BV=1000 IBV=1e-12)\n\n";

    c << "Vdc vin_dc 0 " << Vin << "\n\n";

    // ---- PWM gates (complementary, 50 ns dead-time) ----
    c << "Vpwm_Q1 g1 0 PULSE(0 5 0 10n 10n " << t_q1_on << " " << Tsw << ")\n";
    c << "Vpwm_Q2 g2 0 PULSE(0 5 " << (tA + T_DEAD) << " 10n 10n "
      << t_q2_on << " " << Tsw << ")\n\n";

    // ---- Half-bridge leg ----
    c << "S1 vin_dc sw g1 0 SW1\n";
    c << "D1 sw vin_dc DIDEAL\n";
    c << "Rsnub_Q1 vin_dc sw 1k\nCsnub_Q1 vin_dc sw 1n\n";
    c << "S2 sw 0 g2 0 SW1\n";
    c << "D2 0 sw DIDEAL\n";
    c << "Rsnub_Q2 sw 0 1k\nCsnub_Q2 sw 0 1n\n\n";

    // ---- Primary loop: Vin → Cb → pri_top → Llk → Lpri → sw ----
    // Cb sense source so we can probe i_Cb (= i_pri at steady state).
    c << "Vcb_sense vin_dc cb_lo 0\n";
    c << "C_b cb_lo pri_top " << Cb << " IC=" << VCb_dc << "\n";
    c << "R_cb_esr pri_top pri_top_esr 1m\n";
    c << "L_lk pri_top_esr pri_lk " << Llk << "\n";
    // Vab probe = primary terminal voltage = v(sw) - v(pri_top) (Imbertson sign).
    c << "Evab vab 0 sw pri_top 1\n";
    // Primary current sense in series with the magnetising inductance.
    c << "Vpri_sense pri_lk pri_dot 0\n";
    c << "L_pri pri_dot sw " << Lm << " IC=" << iLm_init << "\n\n";

    // ---- Secondary winding(s) per rectifier type ----
    // K_COUPLING balances two competing effects: too close to 1.0 makes
    // the mutual-inductance matrix near-singular and ngspice rejects it
    // (per plan §12 mistake #12), but K = 0.99 implies an effective
    // leakage of Lm·(1-K²) ≈ Lm·0.02 — for Lm = 1 mH that's 22 µH,
    // which dominates the explicit L_lk = 1 µH and crushes the primary
    // current peak by limiting di/dt during the switch transition.
    // K = 0.999 puts effective leakage at Lm·0.002 ≈ 2 µH which is
    // close to the explicit L_lk and keeps ngspice numerically stable.
    constexpr double K_COUPLING = 0.999;
    if (rect == AhbRectifierType::CENTER_TAPPED) {
        // Two half-windings of n turns each (relative to primary). Each
        // half-secondary's inductance = Lm / n^2.
        const double Lsec_half = Lm / (n * n);
        c << "L_sec_a sec_a sec_ct " << Lsec_half << "\n";
        c << "L_sec_b sec_ct sec_b " << Lsec_half << "\n";
        c << "K1 L_pri L_sec_a " << K_COUPLING << "\n";
        c << "K2 L_pri L_sec_b " << K_COUPLING << "\n";
        c << "K3 L_sec_a L_sec_b " << K_COUPLING << "\n";
        // Diodes from each half to the rectifier output node, ct-tap to gnd
        // through Vct sense for diagnostic. With useSynchronousRectifier the
        // diodes are replaced by SW switches gated complementary to their
        // conducting interval (S_ra conducts during Q1 OFF / Q2 ON ⇒ gate=g2).
        c << "Vsec_a_sense sec_a sec_a_d 0\n";
        c << "Vsec_b_sense sec_b sec_b_d 0\n";
        if (useSynchronousRectifier) {
            c << "S_ra sec_a_d out_rect g2 0 SW1\n";
            c << "S_rb sec_b_d out_rect g1 0 SW1\n";
        } else {
            c << "D_ra sec_a_d out_rect DIDEAL\n";
            c << "D_rb sec_b_d out_rect DIDEAL\n";
        }
        c << "Vct_sense out_gnd sec_ct 0\n";
        c << "L_o out_rect out_node " << Lo << " IC=" << iLo_init << "\n";
        c << "R_lo_dcr out_node out_node_after_dcr 1m\n";
    } else if (rect == AhbRectifierType::FULL_BRIDGE) {
        const double Lsec = Lm / (n * n);
        c << "L_sec sec_a sec_b " << Lsec << "\n";
        c << "K1 L_pri L_sec " << K_COUPLING << "\n";
        c << "Vsec_a_sense sec_a sec_a_d 0\n";
        c << "Vsec_b_sense sec_b sec_b_d 0\n";
        // 4-diode bridge. SR variant pairs each diode with the gate of the
        // primary switch driving the conducting half-cycle: D_r1/D_r4 conduct
        // during Q1 ON (sec_a positive) ⇒ gate g1; D_r2/D_r3 during Q2 ON
        // (sec_b positive) ⇒ gate g2.
        if (useSynchronousRectifier) {
            c << "S_r1 sec_a_d out_rect g1 0 SW1\n";
            c << "S_r2 sec_b_d out_rect g2 0 SW1\n";
            c << "S_r3 out_gnd sec_a_d g2 0 SW1\n";
            c << "S_r4 out_gnd sec_b_d g1 0 SW1\n";
        } else {
            c << "D_r1 sec_a_d out_rect DIDEAL\n";
            c << "D_r2 sec_b_d out_rect DIDEAL\n";
            c << "D_r3 out_gnd sec_a_d DIDEAL\n";
            c << "D_r4 out_gnd sec_b_d DIDEAL\n";
        }
        c << "L_o out_rect out_node " << Lo << " IC=" << iLo_init << "\n";
        c << "R_lo_dcr out_node out_node_after_dcr 1m\n";
    } else { // CURRENT_DOUBLER
        const double Lsec = Lm / (n * n);
        c << "L_sec sec_a sec_b " << Lsec << "\n";
        c << "K1 L_pri L_sec " << K_COUPLING << "\n";
        c << "Vsec_a_sense sec_a sec_a_d 0\n";
        c << "Vsec_b_sense sec_b sec_b_d 0\n";
        // Lo1 charges from sec_a in interval A; Lo2 charges from sec_b in C.
        c << "L_o  sec_a_d out_node " << Lo << " IC=" << iLo_init << "\n";
        c << "L_o2 sec_b_d out_node " << Lo << " IC=" << iLo_init << "\n";
        // Freewheel diodes anchor each Lo's input to gnd during off-half.
        // SR: D_r1 freewheels Lo1 (sec_a) during Q2 ON ⇒ gate g2;
        //     D_r2 freewheels Lo2 (sec_b) during Q1 ON ⇒ gate g1.
        if (useSynchronousRectifier) {
            c << "S_r1 out_gnd sec_a_d g2 0 SW1\n";
            c << "S_r2 out_gnd sec_b_d g1 0 SW1\n";
        } else {
            c << "D_r1 out_gnd sec_a_d DIDEAL\n";
            c << "D_r2 out_gnd sec_b_d DIDEAL\n";
        }
        // Single shared Co/Rload.
        c << "R_lo_dcr out_node out_node_after_dcr 1m\n";
    }

    // ---- Output cap, load, and probes ----
    c << "R_co_esr out_node_after_dcr co_top 1m\n";
    c << "C_o co_top out_gnd " << Co << " IC=" << Vo << "\n";
    c << "R_load co_top out_gnd " << Rload << "\n";
    c << "Vout_sense out_gnd 0 0\n\n";

    // ---- Analysis directives ----
    c << ".tran " << stepTime << " " << simTime << " " << startTime << " uic\n\n";

    c << ".save v(vab) v(sw) v(pri_top) v(out_node)"
      << " i(Vpri_sense) i(Vcb_sense) i(Vsec_a_sense)"
      << " i(Vsec_b_sense) i(Vout_sense) i(Vdc)\n\n";

    c << ".options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500\n";
    c << ".options METHOD=GEAR TRTOL=7\n";
    // NOTE: do NOT emit a `.ic v(sw)=...` line here. With `uic` on .tran the
    // switch-node IC is fought between (a) the body-diode operating point of
    // Q1/Q2 and (b) the snubber RC time constant; the resulting timestep is
    // microscopic and ngspice aborts at t=0 ("Timestep too small ... d1").
    // The reactive elements already carry their own IC= clauses, which is the
    // important pre-charge needed for steady-state convergence within
    // numSteadyStatePeriods. .nodeset is a soft hint for the DC OP only.
    c << ".nodeset v(co_top)=" << Vo
      << " v(pri_top)=" << (Vin - VCb_dc)
      << " v(out_rect)=" << (Vo + 0.7)
      << " v(sec_a_d)=" << (Vo + 0.7)
      << " v(sec_b_d)=" << (Vo + 0.7) << "\n\n";

    c << ".end\n";
    return c.str();
}


// =========================================================================
// simulate_and_extract_operating_points (P7 deliverable)
//
// Runs ngspice for every (inputVoltage, operatingPoint) pair and
// extracts a populated OperatingPoint. Falls back to the analytical
// process_operating_point_for_input_voltage on simulator unavailability
// or per-OP failure, tagging fallback results with "[analytical]".
// Mirrors the PSHB pattern (PSHB cpp lines 826-907).
// =========================================================================
std::vector<OperatingPoint>
AsymmetricHalfBridge::simulate_and_extract_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    if (turnsRatios.empty())
        throw std::invalid_argument(
            "AsymmetricHalfBridge::simulate_and_extract_operating_points: "
            "turnsRatios empty");
    if (!(magnetizingInductance > 0.0))
        throw std::invalid_argument(
            "AsymmetricHalfBridge::simulate_and_extract_operating_points: "
            "magnetizingInductance must be > 0");

    std::vector<OperatingPoint> result;
    NgspiceRunner runner;
    if (!runner.is_available()) {
        std::cerr << "[AHB] ngspice not available — falling back to analytical "
                     "operating-point generation" << std::endl;
        return process_operating_points(turnsRatios, magnetizingInductance);
    }

    AhbRectifierType rect = get_rectifier_type().value_or(
        AhbRectifierType::CENTER_TAPPED);

    std::vector<double> inputVoltages;
    std::vector<std::string> inputVoltageNames;
    collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltageNames);

    auto& ops = get_operating_points();
    for (size_t vinIdx = 0; vinIdx < inputVoltages.size(); ++vinIdx) {
        for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
            const double fsw = ops[opIdx].get_switching_frequency();
            std::string netlist = generate_ngspice_circuit(
                turnsRatios, magnetizingInductance, vinIdx, opIdx);

            SimulationConfig config;
            config.frequency        = fsw;
            config.extractOnePeriod = true;
            config.numberOfPeriods  = static_cast<size_t>(numPeriodsToExtract);
            config.keepTempFiles    = false;

            auto sim = runner.run_simulation(netlist, config);
            std::string tag = inputVoltageNames[vinIdx] + " input ("
                + std::to_string(static_cast<int>(inputVoltages[vinIdx])) + "V)";

            if (!sim.success) {
                std::cerr << "[AHB] ngspice run failed: " << sim.errorMessage
                          << " — falling back to analytical for this op."
                          << std::endl;
                auto analytical = process_operating_point_for_input_voltage(
                    inputVoltages[vinIdx], ops[opIdx], turnsRatios,
                    magnetizingInductance);
                analytical.set_name(tag + " [analytical]");
                result.push_back(analytical);
                continue;
            }

            NgspiceRunner::WaveformNameMapping mapping;
            std::vector<std::string> windingNames;
            std::vector<bool> flipSign;

            mapping.push_back({{"voltage", "vab"},
                               {"current", "vpri_sense#branch"}});
            windingNames.push_back("Primary");
            flipSign.push_back(false);

            if (rect == AhbRectifierType::CENTER_TAPPED) {
                mapping.push_back({{"voltage", "out_node"},
                                   {"current", "vsec_a_sense#branch"}});
                windingNames.push_back("Secondary 0a");
                flipSign.push_back(false);
                mapping.push_back({{"voltage", "out_node"},
                                   {"current", "vsec_b_sense#branch"}});
                windingNames.push_back("Secondary 0b");
                flipSign.push_back(false);
            } else {
                mapping.push_back({{"voltage", "out_node"},
                                   {"current", "vsec_a_sense#branch"}});
                windingNames.push_back("Secondary 0");
                flipSign.push_back(false);
            }

            OperatingPoint extracted = NgspiceRunner::extract_operating_point(
                sim, mapping, fsw, windingNames,
                ops[opIdx].get_ambient_temperature(), flipSign);
            extracted.set_name(tag);
            result.push_back(extracted);
        }
    }
    return result;
}


// =========================================================================
// simulate_and_extract_topology_waveforms (P7 deliverable)
//
// Returns ConverterWaveforms for plotting / NRMSE comparison. Throws
// on simulator unavailability (no analytical fallback — this entry
// point exists specifically to surface SPICE-specific signals such as
// v(sw) and i(C_b) that the analytical model does not produce).
// =========================================================================
std::vector<ConverterWaveforms>
AsymmetricHalfBridge::simulate_and_extract_topology_waveforms(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t numberOfPeriods)
{
    if (turnsRatios.empty())
        throw std::invalid_argument(
            "AsymmetricHalfBridge::simulate_and_extract_topology_waveforms: "
            "turnsRatios empty");
    if (!(magnetizingInductance > 0.0))
        throw std::invalid_argument(
            "AsymmetricHalfBridge::simulate_and_extract_topology_waveforms: "
            "magnetizingInductance must be > 0");

    std::vector<ConverterWaveforms> results;
    NgspiceRunner runner;
    if (!runner.is_available())
        throw std::runtime_error(
            "AsymmetricHalfBridge::simulate_and_extract_topology_waveforms: "
            "ngspice is not available");

    const int origNum = numPeriodsToExtract;
    numPeriodsToExtract = static_cast<int>(numberOfPeriods);

    std::vector<double> inputVoltages;
    std::vector<std::string> inputVoltageNames;
    collect_input_voltages(get_input_voltage(), inputVoltages, inputVoltageNames);

    for (size_t vinIdx = 0; vinIdx < inputVoltages.size(); ++vinIdx) {
        for (size_t opIdx = 0; opIdx < get_operating_points().size(); ++opIdx) {
            const auto& op = get_operating_points()[opIdx];
            const double fsw = op.get_switching_frequency();
            std::string netlist = generate_ngspice_circuit(
                turnsRatios, magnetizingInductance, vinIdx, opIdx);

            SimulationConfig config;
            config.frequency        = fsw;
            config.extractOnePeriod = true;
            config.numberOfPeriods  = numberOfPeriods;
            config.keepTempFiles    = false;

            auto sim = runner.run_simulation(netlist, config);
            if (!sim.success) {
                numPeriodsToExtract = origNum;
                throw std::runtime_error(
                    std::string("AHB simulation failed: ") + sim.errorMessage);
            }

            std::map<std::string, size_t> nameToIdx;
            for (size_t i = 0; i < sim.waveformNames.size(); ++i) {
                std::string lower = sim.waveformNames[i];
                std::transform(lower.begin(), lower.end(), lower.begin(),
                               ::tolower);
                nameToIdx[lower] = i;
            }
            auto getWfm = [&](const std::string& name) -> Waveform {
                std::string lower = name;
                std::transform(lower.begin(), lower.end(), lower.begin(),
                               ::tolower);
                auto it = nameToIdx.find(lower);
                if (it != nameToIdx.end()) return sim.waveforms[it->second];
                return Waveform();
            };

            ConverterWaveforms wf;
            wf.set_switching_frequency(fsw);
            std::string tag = inputVoltageNames[vinIdx] + " input";
            if (get_operating_points().size() > 1)
                tag += " op. point " + std::to_string(opIdx);
            wf.set_operating_point_name(tag);

            wf.set_input_voltage(getWfm("vab"));
            wf.set_input_current(getWfm("vpri_sense#branch"));
            wf.get_mutable_output_voltages().push_back(getWfm("out_node"));
            wf.get_mutable_output_currents().push_back(getWfm("vout_sense#branch"));

            results.push_back(wf);
        }
    }

    numPeriodsToExtract = origNum;
    return results;
}

// =========================================================================
// AdvancedAsymmetricHalfBridge::process
//
// User-pinned-parameters entry point used by the WASM wizard binding
// (libMKF::calculate_ahb_inputs). Mirrors AdvancedPshb::process
// exactly: builds DesignRequirements from the auto-sizing path, then
// overrides any field the user pinned via desired*. Calls
// process_operating_points so the operating-point waveforms (primary,
// secondaries, Lo[, Lo2], Cb) are populated and the per-OP diagnostics
// (lastDutyCycle, lastZvsMargin, ...) are filled in.
//
// Pre-P12 lite implementation: AdvancedAHB only adds setters for
// turnsRatio, magnetizingInductance, leakageInductance, outputInductance,
// dcBlockingCapacitance, outputCapacitance — same set as the PSHB sibling
// (with the AHB-specific addition of Cb). Full AdvancedAHB extras
// (separate "designed only" mode, output-cap sizing override) wait for
// P12 per ASYMMETRIC_HALF_BRIDGE_PLAN.md §11.
// =========================================================================
Inputs AdvancedAsymmetricHalfBridge::process() {
    auto designRequirements = process_design_requirements();

    if (!desiredTurnsRatios.empty()) {
        designRequirements.get_mutable_turns_ratios().clear();
        for (auto n : desiredTurnsRatios) {
            DimensionWithTolerance nTol;
            nTol.set_nominal(n);
            designRequirements.get_mutable_turns_ratios().push_back(nTol);
        }
        // CT carries TWO physical secondary windings. If the wizard only
        // pinned one ratio (the typical case — the wizard only exposes a
        // single "turnsRatio" field), duplicate it so turns_ratios.size()
        // matches the secondary-winding count expected by the coil
        // pipeline (see process_design_requirements rationale).
        AhbRectifierType rectAdv = get_rectifier_type().value_or(
            AhbRectifierType::CENTER_TAPPED);
        if (rectAdv == AhbRectifierType::CENTER_TAPPED &&
            designRequirements.get_mutable_turns_ratios().size() == 1) {
            designRequirements.get_mutable_turns_ratios().push_back(
                designRequirements.get_mutable_turns_ratios()[0]);
        }
    }

    double Lm = desiredMagnetizingInductance;
    if (Lm > 0) {
        DimensionWithTolerance LmTol;
        LmTol.set_minimum(Lm);
        designRequirements.set_magnetizing_inductance(LmTol);
        set_computed_magnetizing_inductance(Lm);
    } else {
        // Fall back to the auto-sized value the base class already picked.
        Lm = get_computed_magnetizing_inductance();
    }

    if (desiredLeakageInductance.has_value())
        set_computed_leakage_inductance(desiredLeakageInductance.value());
    if (desiredOutputInductance.has_value())
        set_computed_output_inductance(desiredOutputInductance.value());
    if (desiredDcBlockingCapacitance.has_value())
        set_computed_dc_blocking_capacitance(desiredDcBlockingCapacitance.value());
    if (desiredOutputCapacitance.has_value())
        set_computed_output_capacitance(desiredOutputCapacitance.value());

    std::vector<double> turnsRatios;
    if (!desiredTurnsRatios.empty()) {
        turnsRatios = desiredTurnsRatios;
    } else {
        for (const auto& tr : designRequirements.get_turns_ratios()) {
            if (tr.get_nominal().has_value())
                turnsRatios.push_back(tr.get_nominal().value());
        }
    }

    auto ops = process_operating_points(turnsRatios, Lm);

    Inputs inputs;
    inputs.set_design_requirements(designRequirements);
    inputs.set_operating_points(ops);
    return inputs;
}

} // namespace OpenMagnetics
