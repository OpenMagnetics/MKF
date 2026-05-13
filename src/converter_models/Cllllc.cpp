// Version: 1.0 (skeleton) — CLLLC bidirectional symmetric resonant converter.
//
// STATUS: This is a v1 SKELETON. The schema, MAS class, header API, and SPICE
// defaults are complete and integrated; the analytical solver, SPICE wrapper,
// and tests are still to be written (see CLLLC_PLAN.md sections A.5–A.10 +
// the per-OP diagnostics table in §A.3).
//
// Implementations marked with "TODO(cllllc-v1)" throw std::logic_error so any
// caller that touches an unimplemented surface fails loudly rather than
// silently returning bogus data — per OpenMagnetics' "no fallbacks, no
// defaults, throw" rule.

#include "converter_models/Cllllc.h"
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace OpenMagnetics {

// =====================================================================
// Construction
// =====================================================================
Cllllc::Cllllc(const json& j) {
    from_json(j, *static_cast<CllllcResonant*>(this));
}

AdvancedCllllc::AdvancedCllllc(const json& j) {
    from_json(j, *this);
}

// =====================================================================
// Bridge / frequency / ratio accessors
// =====================================================================
double Cllllc::get_primary_bridge_voltage_factor() const {
    auto bt = get_bridge_type_primary();
    if (bt.has_value() && bt.value() == LlcBridgeType::FULL_BRIDGE) return 1.0;
    return 0.5;
}

double Cllllc::get_secondary_bridge_voltage_factor() const {
    auto bt = get_bridge_type_secondary();
    if (bt.has_value() && bt.value() == LlcBridgeType::FULL_BRIDGE) return 1.0;
    return 0.5;
}

double Cllllc::get_effective_resonant_frequency() const {
    if (get_primary_resonant_frequency().has_value())
        return get_primary_resonant_frequency().value();
    return std::sqrt(get_min_switching_frequency() * get_max_switching_frequency());
}

double Cllllc::get_effective_inductance_ratio_k() const {
    auto kBase = MAS::CllllcResonant::get_inductance_ratio_k();
    if (kBase.has_value() && kBase.value() > 0) return kBase.value();
    return computedInductanceRatioK > 0 ? computedInductanceRatioK : 5.0;
}

double Cllllc::get_effective_tank_symmetry_ratio() const {
    auto s = get_tank_symmetry_ratio();
    if (s.has_value()) return s.value();
    return 1.0;  // symmetric by default — the headline CLLLC property.
}

// =====================================================================
// Static analytical helpers (these CAN be unit-tested in isolation
// before the solver is finished — they're pure math.)
// =====================================================================
double Cllllc::compute_primary_resonant_frequency(double Lr1, double Cr1) {
    if (Lr1 <= 0 || Cr1 <= 0)
        throw std::invalid_argument(
            "Cllllc::compute_primary_resonant_frequency: Lr1 and Cr1 must be > 0");
    return 1.0 / (2.0 * M_PI * std::sqrt(Lr1 * Cr1));
}

double Cllllc::compute_magnetizing_resonant_frequency_fwd(
        double Lr1, double Lm, double Cr1) {
    if (Lr1 <= 0 || Lm <= 0 || Cr1 <= 0)
        throw std::invalid_argument(
            "Cllllc::compute_magnetizing_resonant_frequency_fwd: positive args required");
    return 1.0 / (2.0 * M_PI * std::sqrt((Lr1 + Lm) * Cr1));
}

double Cllllc::compute_magnetizing_resonant_frequency_rev(
        double Lr2, double Lm, double n, double Cr2) {
    if (Lr2 <= 0 || Lm <= 0 || n <= 0 || Cr2 <= 0)
        throw std::invalid_argument(
            "Cllllc::compute_magnetizing_resonant_frequency_rev: positive args required");
    // From the LV-side reference frame: Lm reflected as Lm/n^2.
    return 1.0 / (2.0 * M_PI * std::sqrt((Lr2 + Lm / (n * n)) * Cr2));
}

double Cllllc::compute_fha_gain(double f_n, double K, double Q,
                                double a, double b) {
    if (f_n <= 0 || K <= 0)
        throw std::invalid_argument("Cllllc::compute_fha_gain: f_n and K must be > 0");
    // Symmetric (a = b = 1) reduces to LLC FHA shape; asymmetric per
    // IEEE 10362332 / MDPI Electronics 12(7):1605 §3 closed form.
    // v1 uses the symmetric reduction and applies a small (a+b)/2 weighting
    // for the asymmetric branch — refine when the asymmetric solver lands.
    double ab = 0.5 * (a + b);
    double term1 = 1.0 + ab * K * (1.0 - 1.0 / (f_n * f_n));
    double term2 = Q * (f_n - 1.0 / f_n);
    double denom = std::sqrt(term1 * term1 + term2 * term2);
    if (denom <= 0)
        throw std::runtime_error("Cllllc::compute_fha_gain: denominator collapsed");
    return 1.0 / denom;
}

double Cllllc::compute_inductance_ratio_K(double Lm, double Lr1) {
    if (Lr1 <= 0 || Lm <= 0)
        throw std::invalid_argument("Cllllc::compute_inductance_ratio_K: Lm and Lr1 must be > 0");
    return Lm / Lr1;
}

double Cllllc::compute_quality_factor(double Lr1, double Cr1, double R_load_reflected) {
    if (Lr1 <= 0 || Cr1 <= 0 || R_load_reflected <= 0)
        throw std::invalid_argument("Cllllc::compute_quality_factor: positive args required");
    return std::sqrt(Lr1 / Cr1) / R_load_reflected;
}

double Cllllc::compute_zvs_lm_max(double T_dead, double C_oss, double Fs) {
    if (T_dead <= 0 || C_oss <= 0 || Fs <= 0)
        throw std::invalid_argument("Cllllc::compute_zvs_lm_max: positive args required");
    // Huang SLUP263 §6 (extended for CLLLC mirror): Lm <= T_dead / (16*Coss*Fs)
    return T_dead / (16.0 * C_oss * Fs);
}

double Cllllc::compute_turns_ratio(double V_HV, double V_LV) {
    if (V_HV <= 0 || V_LV <= 0)
        throw std::invalid_argument("Cllllc::compute_turns_ratio: positive args required");
    return V_HV / V_LV;
}

double Cllllc::compute_symmetric_lr2(double Lr1, double n) {
    if (Lr1 <= 0 || n <= 0)
        throw std::invalid_argument("Cllllc::compute_symmetric_lr2: positive args required");
    return Lr1 / (n * n);
}

double Cllllc::compute_symmetric_cr2(double Cr1, double n) {
    if (Cr1 <= 0 || n <= 0)
        throw std::invalid_argument("Cllllc::compute_symmetric_cr2: positive args required");
    return Cr1 * (n * n);
}

// =====================================================================
// run_checks — validates the input deck (analogue of Llc::run_checks)
// =====================================================================
bool Cllllc::run_checks(bool assertErrors) {
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
DesignRequirements Cllllc::process_design_requirements() {
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

    // Turns ratio sized at nominal so M = 1 at f_r1 (Cllllc plan §A.5).
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

    double Q  = get_quality_factor().value_or(0.4);
    double K  = MAS::CllllcResonant::get_inductance_ratio_k().value_or(5.0);
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
    designRequirements.set_topology(Topologies::CLLLLC_RESONANT_CONVERTER);

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
// process_operating_points — TODO(cllllc-v1): full 5-state Nielsen
// solver per CLLLC_PLAN.md §A.5. Currently throws so callers fail loudly
// rather than getting bogus data.
// =====================================================================
std::vector<OperatingPoint> Cllllc::process_operating_points(
        const std::vector<double>& /*turnsRatios*/,
        double /*magnetizingInductance*/) {
    throw std::logic_error(
        "Cllllc::process_operating_points not yet implemented. "
        "v1 ships schema + design requirements + static helpers only; "
        "the 5-state Nielsen solver (CLLLC_PLAN.md §A.5) lands in v2.");
}

OperatingPoint Cllllc::process_operating_point_for_input_voltage(
        double /*inputVoltage*/,
        const CllllcOperatingPoint& /*op*/,
        const std::vector<double>& /*turnsRatios*/,
        double /*magnetizingInductance*/) {
    throw std::logic_error(
        "Cllllc::process_operating_point_for_input_voltage not yet implemented. "
        "Pending 5-state Nielsen solver in CLLLC_PLAN.md §A.5.");
}

// =====================================================================
// SPICE — TODO(cllllc-v1): direction-aware 8-switch netlist with three
// K-couplings (K1, K2, K_pri_sec) per CLLLC_PLAN.md §A.6.
// =====================================================================
std::string Cllllc::generate_ngspice_circuit(
        const std::vector<double>& /*turnsRatios*/,
        double /*magnetizingInductance*/,
        size_t /*inputVoltageIndex*/,
        size_t /*operatingPointIndex*/) {
    throw std::logic_error(
        "Cllllc::generate_ngspice_circuit not yet implemented. "
        "Pending direction-aware 8-switch netlist (CLLLC_PLAN.md §A.6).");
}

std::vector<OperatingPoint> Cllllc::simulate_and_extract_operating_points(
        const std::vector<double>& /*turnsRatios*/,
        double /*magnetizingInductance*/) {
    throw std::logic_error(
        "Cllllc::simulate_and_extract_operating_points not yet implemented. "
        "Depends on generate_ngspice_circuit (CLLLC_PLAN.md §A.7).");
}

std::vector<ConverterWaveforms> Cllllc::simulate_and_extract_topology_waveforms(
        const std::vector<double>& /*turnsRatios*/,
        double /*magnetizingInductance*/,
        size_t /*numberOfPeriods*/) {
    throw std::logic_error(
        "Cllllc::simulate_and_extract_topology_waveforms not yet implemented. "
        "Depends on generate_ngspice_circuit (CLLLC_PLAN.md §A.8).");
}

// =====================================================================
// Extra components — TODO(cllllc-v1) per CLLLC_PLAN.md §A.9 (transformer
// with two integrated leakages OR two discrete inductors, plus Cr1, Cr2,
// HV bus cap, LV bus cap).
// =====================================================================
std::vector<std::variant<Inputs, CAS::Inputs>> Cllllc::get_extra_components_inputs(
        ExtraComponentsMode /*mode*/,
        std::optional<Magnetic> /*magnetic*/) {
    throw std::logic_error(
        "Cllllc::get_extra_components_inputs not yet implemented. "
        "Pending solver (needs operating-point waveforms to size Cr1/Cr2/bus caps). "
        "See CLLLC_PLAN.md §A.9.");
}

// =====================================================================
// AdvancedCllllc::process — TODO(cllllc-v1)
// =====================================================================
Inputs AdvancedCllllc::process() {
    throw std::logic_error(
        "AdvancedCllllc::process not yet implemented. "
        "Depends on Cllllc::process_operating_points.");
}

} // namespace OpenMagnetics
