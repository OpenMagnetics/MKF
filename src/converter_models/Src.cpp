// =====================================================================
// Series Resonant Converter (SRC) — Phase 1+2 implementation.
//
// Phase 1: skeleton (constructor, run_checks, design-requirements stub).
// Phase 2: FHA (First-Harmonic-Approximation) steady-state solver.
//
// FHA assumes the resonant tank current is purely sinusoidal at the
// switching frequency. This is exact at fsw = fr and accurate to within
// a few percent for above-resonance CCM (Λ ≥ 1) at moderate Q (1–4).
//
// Below-resonance CCM (capacitive region) and DCM are not solved here —
// process_operating_point_for_input_voltage throws on Λ < 1 (deferred to
// future phases per SRC_PLAN.md §11).
//
// REFERENCES
//   [1] R. L. Steigerwald, "A Comparison of Half-Bridge Resonant Converter
//       Topologies", IEEE TPEL 1988.
//   [2] M. K. Kazimierczuk, "Resonant Power Converters", Wiley 2nd ed., Ch.4.
//   [3] src/converter_models/SRC_PLAN.md (Phase 2 scope).
// =====================================================================
#include "converter_models/Src.h"
#include "converter_models/Topology.h"
#include "physical_models/LeakageInductance.h"
#include "support/Utils.h"
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

namespace OpenMagnetics {

// =====================================================================
// Constructors
// =====================================================================
Src::Src(const json& j) {
    from_json(j, *static_cast<SeriesResonant*>(this));
}

AdvancedSrc::AdvancedSrc(const json& j) {
    from_json(j, *this);
}

// =====================================================================
// Helpers
// =====================================================================
double Src::get_bridge_voltage_factor() const {
    auto bt = get_bridge_type();
    if (bt.has_value() &&
        (bt.value() == SrcBridgeType::FULL_BRIDGE ||
         bt.value() == SrcBridgeType::FULL_BRIDGE_PHASE_SHIFT)) {
        return 1.0;
    }
    return 0.5;  // half-bridge (default per SRC_PLAN.md)
}

double Src::get_effective_resonant_frequency() const {
    if (get_resonant_frequency().has_value() && get_resonant_frequency().value() > 0)
        return get_resonant_frequency().value();
    if (computedResonantFrequency > 0) return computedResonantFrequency;

    // Seed: geometric mean of fsw window. We require both bounds be positive.
    double fmin = get_min_switching_frequency();
    double fmax = get_max_switching_frequency();
    if (fmin <= 0 || fmax <= 0)
        throw std::runtime_error(
            "SRC: cannot derive resonant frequency — no resonantFrequency set, "
            "no computed tank yet, and switching-frequency bounds are non-positive");
    return std::sqrt(fmin * fmax);
}

// =====================================================================
// run_checks
// =====================================================================
bool Src::run_checks(bool assertErrors) {
    _assertErrors = assertErrors;
    bool ok = true;
    auto& ops = get_operating_points();
    if (ops.empty()) {
        if (assertErrors) throw std::runtime_error("SRC: no operating points");
        return false;
    }
    if (get_min_switching_frequency() <= 0 || get_max_switching_frequency() <= 0) {
        if (assertErrors) throw std::runtime_error("SRC: switching-frequency bounds must be > 0");
        return false;
    }
    if (get_min_switching_frequency() > get_max_switching_frequency()) {
        if (assertErrors) throw std::runtime_error("SRC: minSwitchingFrequency > maxSwitchingFrequency");
        return false;
    }
    for (size_t i = 0; i < ops.size(); ++i) {
        const auto& op = ops[i];
        if (op.get_output_voltages().empty() || op.get_output_currents().empty()) {
            if (assertErrors) throw std::runtime_error("SRC: OP missing output voltages/currents");
            ok = false;
            continue;
        }
        if (op.get_output_voltages().size() != op.get_output_currents().size()) {
            if (assertErrors) throw std::runtime_error("SRC: voltage/current count mismatch");
            ok = false;
        }
        double fsw = op.get_switching_frequency();
        if (fsw < get_min_switching_frequency() * 0.99 ||
            fsw > get_max_switching_frequency() * 1.01) {
            if (assertErrors) throw std::runtime_error("SRC: fsw out of [min,max] range");
            ok = false;
        }
    }
    return ok;
}

// =====================================================================
// process_design_requirements
//
// SRC has TWO independent design parameters (one fewer than LLC because
// there is no Lm branch):
//   1. Q  (quality factor)     → controls Lr, Cr (tank impedance)
//   2. fr (resonant frequency) → sets the operating-region centre
//
// From these:
//   Vbridge_pk_fund = (4/π) · k_bridge · Vin             (FHA fundamental)
//   At fsw = fr (M=1 from FHA at HB → Vout_max = Vin/2 / n via n = Vin/(2·Vout)):
//     n   = k_bridge · Vin / Vout                         (turns ratio)
//     Rac = (8 · n²) / π² · Rload                         (FHA-reflected AC load)
//     Zr  = Q · Rac                                       (characteristic impedance)
//     Lr  = Zr / (2π·fr)
//     Cr  = 1 / (2π·fr·Zr)
// =====================================================================
DesignRequirements Src::process_design_requirements() {
    double k_bridge = get_bridge_voltage_factor();
    auto& inputVoltage = get_input_voltage();
    auto& ops = get_operating_points();
    if (ops.empty())
        throw std::runtime_error("SRC: at least one operating point is required");
    if (ops[0].get_output_voltages().empty() || ops[0].get_output_currents().empty())
        throw std::runtime_error("SRC: operating point is missing output voltages or currents");
    if (!inputVoltage.get_nominal().has_value() &&
        !(inputVoltage.get_minimum().has_value() && inputVoltage.get_maximum().has_value()))
        throw std::runtime_error("SRC: input voltage must specify nominal (or both min and max)");

    double Vin_nom = inputVoltage.get_nominal().value_or(
        (inputVoltage.get_minimum().value_or(0) + inputVoltage.get_maximum().value_or(0)) / 2.0);
    if (Vin_nom <= 0)
        throw std::runtime_error("SRC: nominal input voltage must be > 0");

    double mainOutputVoltage = ops[0].get_output_voltages()[0];
    double mainOutputCurrent = ops[0].get_output_currents()[0];
    if (mainOutputVoltage <= 0 || mainOutputCurrent <= 0)
        throw std::runtime_error("SRC: main output voltage and current must be > 0");

    // Reflected (FHA-equivalent) bridge voltage referred to secondary side.
    // SRC has |M| ≤ 1 (step-down). At resonance (Λ=1), peak fundamental of
    // bridge midpoint square-wave = (4/π) · k_bridge · Vin, and FHA says
    // Vout (averaged DC) = (2/π) · n · Vbridge_pk_fund · M_fha = k_bridge · Vin / n
    // when M=1. So n = k_bridge · Vin / Vout for FB, CT, and CD alike.
    // (CD's actual gain is ~3× lower than this FHA prediction — see
    // FIXME-src-3 in TestSrcReferenceDesignsPtp.cpp — but n_CD = n_FB
    // remains the canonical published-design choice; the gain shortfall
    // is normally compensated by operating closer to resonance.)
    double Vbridge_dc_equiv = k_bridge * Vin_nom;
    double mainTurnsRatio   = Vbridge_dc_equiv / mainOutputVoltage;

    size_t nOutputs = ops[0].get_output_voltages().size();
    std::vector<double> outputTurnsRatios;
    outputTurnsRatios.reserve(nOutputs);
    for (size_t i = 0; i < nOutputs; ++i) {
        double Vout_i = ops[0].get_output_voltages()[i];
        if (Vout_i <= 0)
            throw std::runtime_error("SRC: output voltage must be > 0");
        outputTurnsRatios.push_back(Vbridge_dc_equiv / Vout_i);
    }

    double Rload = mainOutputVoltage / mainOutputCurrent;
    // Reflected AC resistance — canonical FHA: Rac = (8/π²)·n²·Rload for
    // FB, CT, and CD alike. CD's secondary doesn't fit cleanly in FHA (Vsec
    // ≈ ±Vout, not ±2·Vout; non-sinusoidal harmonics carry significant
    // power), so absolute Vout predictions are 2–3× off, but shape/peak
    // metrics still match between analytical and SPICE because both share
    // the same FHA assumption. See FIXME-src-3.
    double Rac   = (8.0 * mainTurnsRatio * mainTurnsRatio) / (M_PI * M_PI) * Rload;
    double Q     = get_quality_factor().value_or(2.0);  // schema default
    if (Q <= 0)
        throw std::runtime_error("SRC: qualityFactor must be > 0");
    double fr = get_effective_resonant_frequency();
    if (fr <= 0)
        throw std::runtime_error("SRC: resonant frequency must be > 0");

    double Zr = Q * Rac;
    double Lr = Zr / (2.0 * M_PI * fr);
    double Cr = 1.0 / (2.0 * M_PI * fr * Zr);

    // User overrides take precedence (validation against published designs).
    if (userResonantInductance > 0)        Lr = userResonantInductance;
    else if (get_series_inductance().has_value() && get_series_inductance().value() > 0)
        Lr = get_series_inductance().value();
    if (userResonantCapacitance > 0)       Cr = userResonantCapacitance;
    else if (get_resonant_capacitance().has_value() && get_resonant_capacitance().value() > 0)
        Cr = get_resonant_capacitance().value();

    if (Lr <= 0 || Cr <= 0)
        throw std::runtime_error("SRC: derived Lr/Cr non-positive — check Q, fr, load");

    computedResonantInductance  = Lr;
    computedResonantCapacitance = Cr;
    computedResonantFrequency   = 1.0 / (2.0 * M_PI * std::sqrt(Lr * Cr));

    DesignRequirements designRequirements;
    designRequirements.get_mutable_turns_ratios().clear();
    size_t wpo = windings_per_output();
    for (auto n : outputTurnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(roundFloat(n, 2));
        for (size_t w = 0; w < wpo; ++w)
            designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }

    // Magnetizing inductance: SRC has no Lm in the resonant tank — the
    // transformer Lm is sized purely from the requirement that magnetizing
    // current stays well below resonant tank current, NOT from a resonance
    // constraint. We pick Lm = 10·Lr as a reasonable starting point (Im_pk
    // ≈ ILr_pk / 10). Designers can override via AdvancedSrc.
    double Lm_default = 10.0 * Lr;
    DimensionWithTolerance lmTol;
    lmTol.set_nominal(roundFloat(Lm_default, 10));
    designRequirements.set_magnetizing_inductance(lmTol);

    // Treat Lr as a separate (or integrated) leakage requirement — analogous
    // to LLC's integrated_resonant_inductor flag. SRC default is "external Lr",
    // so we always emit a leakage spec.
    std::vector<DimensionWithTolerance> leakageReqs;
    DimensionWithTolerance lrTol;
    lrTol.set_nominal(roundFloat(Lr, 10));
    leakageReqs.push_back(lrTol);
    designRequirements.set_leakage_inductance(leakageReqs);

    // Isolation: primary + windings_per_output() entries per output (CT splits
    // each output into 2 half-windings, FB/CD use 1 winding per output).
    std::vector<IsolationSide> isolationSides;
    isolationSides.push_back(get_isolation_side_from_index(0));
    for (size_t i = 0; i < nOutputs; ++i) {
        for (size_t w = 0; w < wpo; ++w)
            isolationSides.push_back(get_isolation_side_from_index(i + 1));
    }
    designRequirements.set_isolation_sides(isolationSides);
    designRequirements.set_topology(Topologies::SERIES_RESONANT_CONVERTER);

    return designRequirements;
}

// =====================================================================
// process_operating_points
// =====================================================================
std::vector<OperatingPoint> Src::process_operating_points(
    const std::vector<double>& turnsRatios, double /*magnetizingInductance*/)
{
    if (computedResonantInductance <= 0 || computedResonantCapacitance <= 0)
        process_design_requirements();

    // Clear extra-component waveforms — repopulated by per-OP solver.
    extraCapVoltageWaveforms.clear();
    extraCapCurrentWaveforms.clear();
    extraIndVoltageWaveforms.clear();
    extraIndCurrentWaveforms.clear();
    extraLoCurrentWaveforms.clear();
    extraLoVoltageWaveforms.clear();
    extraLo2CurrentWaveforms.clear();
    extraLo2VoltageWaveforms.clear();

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
        auto op = process_operating_point_for_input_voltage(Vin, ops[0], turnsRatios);
        op.set_name(name + " input (" + std::to_string(static_cast<int>(Vin)) + "V)");
        result.push_back(op);
    }
    return result;
}

std::vector<OperatingPoint> Src::process_operating_points(Magnetic /*magnetic*/) {
    auto req = process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios())
        turnsRatios.push_back(resolve_dimensional_values(tr));
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    return process_operating_points(turnsRatios, Lm);
}

// =====================================================================
// process_operating_point_for_input_voltage  (FHA, above-resonance CCM)
//
// Steady-state model (HB; FB scales by 2):
//   Vbridge_pk_fund = (4/π) · k_bridge · Vin                       [V]
//   M(Λ,Q)          = 1 / sqrt(1 + Q² (Λ - 1/Λ)²)                  (HB FHA)
//                     (Steigerwald [1], simplified by FHA absorbing
//                      the 4/π and 8n²/π² factors into M_dc-to-dc).
//   Vout_pred       = M · Vbridge_dc_equiv / n                     [V]
//   ILr_pk          = Vbridge_pk_fund / Z_in(Λ)                    [A]
//     where Z_in(Λ) = sqrt( Rac² + (ωLr - 1/(ωCr))² ),  Rac at primary
//     side = (8/π²) · n² · Rload.
//   φ               = atan2( ωLr - 1/(ωCr), Rac )                  (lagging > 0 above res)
//   VCr_pk          = ILr_pk / (ωCr)                               [V]
//
// Above-resonance ZVS holds when φ > 0 (current lags bridge → switch
// turns on with body diode conducting).
// =====================================================================
OperatingPoint Src::process_operating_point_for_input_voltage(
    double inputVoltage,
    const TopologyExcitation& srcOpPoint,
    const std::vector<double>& turnsRatios)
{
    OperatingPoint operatingPoint;

    if (srcOpPoint.get_output_voltages().empty() || srcOpPoint.get_output_currents().empty())
        throw std::runtime_error("SRC: operating point missing output voltages/currents");

    double fsw = srcOpPoint.get_switching_frequency();
    if (fsw <= 0)
        throw std::runtime_error("SRC: operating-point switching frequency must be > 0");

    double k_bridge = get_bridge_voltage_factor();
    double Lr = computedResonantInductance;
    double Cr = computedResonantCapacitance;
    if (Lr <= 0 || Cr <= 0)
        throw std::runtime_error("SRC: tank not initialized — call process_design_requirements first");

    double fr     = computedResonantFrequency;
    double Lambda = fsw / fr;
    if (Lambda < 1.0) {
        throw std::runtime_error(
            "SRC: below-resonance operation (Λ=" + std::to_string(Lambda) +
            " < 1) is not supported in Phase 2 (capacitive region / hard-switching). "
            "Increase switching frequency above the resonant frequency, or wait "
            "for Phase 3 (variant + below-res support).");
    }

    size_t nOutputs = srcOpPoint.get_output_voltages().size();
    size_t wpo = windings_per_output();
    auto get_n_for_output = [&](size_t outputIdx) -> double {
        if (turnsRatios.empty())
            return (k_bridge * inputVoltage) / srcOpPoint.get_output_voltages()[outputIdx];
        // turnsRatios is laid out per non-primary winding (wpo entries / output).
        // All wpo entries for a given output carry the same n_design (LLC CT
        // convention), so we read the first entry of the output's block.
        size_t flatIdx = outputIdx * wpo;
        if (flatIdx < turnsRatios.size()) return turnsRatios[flatIdx];
        return turnsRatios[0];
    };

    double n_main = get_n_for_output(0);
    if (n_main <= 0)
        throw std::runtime_error("SRC: turns ratio must be > 0");

    double Vout_main = srcOpPoint.get_output_voltages()[0];
    double Iout_main = srcOpPoint.get_output_currents()[0];
    double Rload    = Vout_main / Iout_main;
    // Canonical FHA Rac = (8/π²)·n²·Rload for FB, CT, and CD. CD's
    // secondary doesn't fit FHA cleanly so absolute Vout predictions are
    // 2–3× off (see FIXME-src-3), but shape/peak metrics are consistent
    // between analytical and SPICE.
    double Rac      = (8.0 * n_main * n_main) / (M_PI * M_PI) * Rload;

    double w   = 2.0 * M_PI * fsw;
    double XLr = w * Lr;
    double XCr = 1.0 / (w * Cr);
    double X   = XLr - XCr;            // > 0 above resonance (inductive)
    double Zin = std::sqrt(Rac * Rac + X * X);

    double Vbridge_pk_fund = (4.0 / M_PI) * k_bridge * inputVoltage;
    double ILr_pk          = Vbridge_pk_fund / Zin;
    double VCr_pk          = ILr_pk * XCr;
    double phi             = std::atan2(X, Rac);                  // rad (> 0 above res)
    double M_fha           = 1.0 / std::sqrt(1.0 + (Rac > 0 ? (X / Rac) * (X / Rac) : 0.0));

    lastGainM            = M_fha;
    lastNormalizedFsw    = Lambda;
    lastIrPeak           = ILr_pk;
    lastVcrPeak          = VCr_pk;
    lastIsAboveResonance = (Lambda >= 1.0);

    // ── Build full-period time grid ──
    const int N         = 256;                  // samples per period
    double    period    = 1.0 / fsw;
    double    dt        = period / N;
    int       totalSamples = N + 1;             // closed period (last == first)
    std::vector<double> time_full(totalSamples);
    for (int k = 0; k < totalSamples; ++k)
        time_full[k] = k * dt;

    // ── Primary winding waveforms (sinusoidal current; square bridge voltage) ──
    std::vector<double> ILr_full(totalSamples);
    std::vector<double> Vpri_full(totalSamples);
    std::vector<double> VCr_full(totalSamples);   // FHA: voltage across Cr (lags ILr by 90°)
    std::vector<double> VLr_full(totalSamples);   // FHA: voltage across Lr (leads ILr by 90°)
    double Vbridge_lvl = k_bridge * inputVoltage;   // half-bridge: ±Vin/2; full-bridge: ±Vin
    for (int k = 0; k < totalSamples; ++k) {
        double t = time_full[k];
        double theta = w * t;                       // bridge half-cycle ref (high half: 0..π)
        ILr_full[k]  = ILr_pk * std::sin(theta - phi);
        // VCr = (1/C)·∫i dt = -(ILr_pk/(ωC))·cos(θ-φ)   (FHA, zero-mean)
        VCr_full[k]  = -VCr_pk * std::cos(theta - phi);
        // VLr = Lr·di/dt = ω·Lr·ILr_pk·cos(θ-φ) = XLr·ILr_pk·cos(θ-φ)
        VLr_full[k]  =  XLr * ILr_pk * std::cos(theta - phi);
        // Bridge midpoint voltage referred to bridge-zero-DC frame:
        //   half 0 (0 ≤ θ mod 2π < π) → +Vbridge_lvl
        //   half 1 (π ≤ θ mod 2π < 2π) → -Vbridge_lvl
        double thetaMod = std::fmod(theta, 2.0 * M_PI);
        if (thetaMod < 0) thetaMod += 2.0 * M_PI;
        Vpri_full[k] = (thetaMod < M_PI) ? +Vbridge_lvl : -Vbridge_lvl;
    }

    // Stash Cr / Lr waveforms for get_extra_components_inputs.
    {
        Waveform vcW; vcW.set_ancillary_label(WaveformLabel::CUSTOM);
        vcW.set_data(VCr_full); vcW.set_time(time_full);
        extraCapVoltageWaveforms.push_back(vcW);

        Waveform icW; icW.set_ancillary_label(WaveformLabel::CUSTOM);
        icW.set_data(ILr_full); icW.set_time(time_full);     // same current through Cr and Lr (series)
        extraCapCurrentWaveforms.push_back(icW);

        Waveform vlW; vlW.set_ancillary_label(WaveformLabel::CUSTOM);
        vlW.set_data(VLr_full); vlW.set_time(time_full);
        extraIndVoltageWaveforms.push_back(vlW);

        Waveform ilW; ilW.set_ancillary_label(WaveformLabel::CUSTOM);
        ilW.set_data(ILr_full); ilW.set_time(time_full);
        extraIndCurrentWaveforms.push_back(ilW);
    }

    {
        Waveform currentWaveform;
        currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        currentWaveform.set_data(ILr_full);
        currentWaveform.set_time(time_full);

        Waveform voltageWaveform;
        voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        voltageWaveform.set_data(Vpri_full);
        voltageWaveform.set_time(time_full);

        auto excitation = complete_excitation(currentWaveform, voltageWaveform, fsw, "Primary");
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    // ── Secondary windings (per-output, per-rectifier-type) ──
    //
    // Phase 2+: supports FULL_BRIDGE_DIODE (1 winding/output), CENTER_TAPPED_DIODE
    // (2 half-windings/output, each carrying only its conducting half-cycle of
    // the primary tank current), and CURRENT_DOUBLER (1 secondary winding + 2
    // output filter inductors Lo1/Lo2 — emitted via get_extra_components_inputs).
    auto rt = get_effective_rectifier_type();

    // For current sharing across multi-output secondaries, we weight by load
    // conductance, matching the LLC convention.
    double total_g = 0.0;
    for (size_t i = 0; i < nOutputs; ++i) {
        double Vout_i = srcOpPoint.get_output_voltages()[i];
        double Iout_i = srcOpPoint.get_output_currents()[i];
        if (Vout_i > 0 && Iout_i > 0) total_g += Iout_i / Vout_i;
    }
    if (total_g <= 0) total_g = 1.0;

    for (size_t outputIdx = 0; outputIdx < nOutputs; ++outputIdx) {
        double n_i = get_n_for_output(outputIdx);
        if (n_i <= 0) n_i = 1.0;

        double Vout_i = srcOpPoint.get_output_voltages()[outputIdx];
        double Iout_i = srcOpPoint.get_output_currents()[outputIdx];
        double share  = (Vout_i > 0 && Iout_i > 0)
                         ? (Iout_i / Vout_i) / total_g
                         : (outputIdx == 0 ? 1.0 : 0.0);

        if (rt == SrcRectifierType::CENTER_TAPPED_DIODE) {
            // CT: two half-windings. Each conducts on alternating half-cycles.
            //   halfIdx=0 → conducts when primary tank current > 0
            //              ⇒ i_top = max(0, +ILr·n·share), v_top = +Vout (D1 on)
            //              ⇒ when D2 is on instead, v_top = -Vout (reflected)
            //   halfIdx=1 → symmetric, swapped polarity.
            // n_i is primary:half-secondary (same convention as LLC CT).
            for (size_t halfIdx = 0; halfIdx < 2; ++halfIdx) {
                std::vector<double> iSecData(totalSamples, 0.0);
                std::vector<double> vSecData(totalSamples, 0.0);
                for (int k = 0; k < totalSamples; ++k) {
                    double iPri = ILr_full[k];
                    double i_share = iPri * n_i * share;
                    double i_half = (halfIdx == 0)
                                    ? std::max(0.0, +i_share)
                                    : std::max(0.0, -i_share);
                    // Each half-winding voltage flips sign when the opposite
                    // diode is conducting. halfIdx=0: +Vout when iPri>0,
                    // -Vout when iPri<0. halfIdx=1: opposite.
                    double v_half = (halfIdx == 0)
                                    ? ((iPri >= 0.0) ? +Vout_i : -Vout_i)
                                    : ((iPri >= 0.0) ? -Vout_i : +Vout_i);
                    iSecData[k] = std::isfinite(i_half) ? i_half : 0.0;
                    vSecData[k] = std::isfinite(v_half) ? v_half : 0.0;
                }
                Waveform secCurrentWfm;
                secCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
                secCurrentWfm.set_data(iSecData);
                secCurrentWfm.set_time(time_full);

                Waveform secVoltageWfm;
                secVoltageWfm.set_ancillary_label(WaveformLabel::CUSTOM);
                secVoltageWfm.set_data(vSecData);
                secVoltageWfm.set_time(time_full);

                std::string windingName = "Secondary " + std::to_string(outputIdx)
                                        + " Half " + std::to_string(halfIdx + 1);
                auto excitation = complete_excitation(secCurrentWfm, secVoltageWfm, fsw, windingName);
                operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
            }
        }
        else if (rt == SrcRectifierType::CURRENT_DOUBLER) {
            // CD: single bipolar secondary winding feeding two diodes and two
            // output filter inductors Lo1/Lo2. Under huge-Lo CCM, the tank
            // sees an electrically FB-equivalent rectifier (Vsec swings ±Vo
            // square because one terminal is always clamped to Vo by an
            // active diode and the other is held near 0V by a freewheeling
            // Lo). The CD-specific output behaviour (Iout split 50/50 across
            // Lo1/Lo2 with ripple cancellation at 2·fs) is captured in the
            // Lo waveforms stashed for get_extra_components_inputs.
            std::vector<double> iSecData(totalSamples, 0.0);
            std::vector<double> vSecData(totalSamples, 0.0);
            std::vector<double> iLo1Data(totalSamples, 0.0);
            std::vector<double> iLo2Data(totalSamples, 0.0);
            std::vector<double> vLo1Data(totalSamples, 0.0);
            std::vector<double> vLo2Data(totalSamples, 0.0);
            const double Iout_dc = (Vout_i > 0 && Iout_i > 0) ? Iout_i : 0.0;
            for (int k = 0; k < totalSamples; ++k) {
                double iPri = ILr_full[k];
                double I_sec = iPri * n_i * share;             // bipolar; matches FB sign
                double V_sec = (iPri >= 0.0) ? +Vout_i : -Vout_i;
                iSecData[k] = std::isfinite(I_sec) ? I_sec : 0.0;
                vSecData[k] = std::isfinite(V_sec) ? V_sec : 0.0;
                // Lo1/Lo2 — KCL at output node: ILo1 + ILo2 = IoutDC. Each
                // carries IoutDC/2 mean with opposite-sign ripple sourced from
                // the secondary current (which averages 0 over a period).
                double ripple = I_sec / 4.0;                   // peak excursion ≈ |I_sec|/4
                iLo1Data[k] = Iout_dc / 2.0 + ripple;
                iLo2Data[k] = Iout_dc / 2.0 - ripple;
                // Each Lo sees ±Vo: +Vo when its diode is conducting (Lo
                // charging from the rectified secondary), −Vo when
                // freewheeling. Volt-second balance is exact.
                vLo1Data[k] = (I_sec > 0) ? (+Vout_i) : (-Vout_i);
                vLo2Data[k] = (I_sec < 0) ? (+Vout_i) : (-Vout_i);
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

            // Stash Lo waveforms for get_extra_components_inputs (one entry per
            // output appended in outputIdx order; consumers must respect that).
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
        }
        else {
            // FULL_BRIDGE_DIODE: single secondary winding.
            std::vector<double> iSecData(totalSamples, 0.0);
            std::vector<double> vSecData(totalSamples, 0.0);
            for (int k = 0; k < totalSamples; ++k) {
                // Secondary current (winding) = primary current × n × share
                // (full-wave bridge: secondary winding sees full sinusoid; the
                // diode bridge handles polarity inversion to feed Vout.)
                double iSec = ILr_full[k] * n_i * share;
                // Secondary winding voltage = +Vout_i during positive half of
                // primary current (rectifier conducts D1/D4), -Vout_i otherwise.
                double vSec = (ILr_full[k] >= 0.0) ? +Vout_i : -Vout_i;
                iSecData[k] = std::isfinite(iSec) ? iSec : 0.0;
                vSecData[k] = std::isfinite(vSec) ? vSec : 0.0;
            }

            Waveform secCurrentWfm;
            secCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            secCurrentWfm.set_data(iSecData);
            secCurrentWfm.set_time(time_full);

            Waveform secVoltageWfm;
            secVoltageWfm.set_ancillary_label(WaveformLabel::CUSTOM);
            secVoltageWfm.set_data(vSecData);
            secVoltageWfm.set_time(time_full);

            std::string windingName = "Secondary " + std::to_string(outputIdx);
            auto excitation = complete_excitation(secCurrentWfm, secVoltageWfm, fsw, windingName);
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }
    }

    OperatingConditions conditions;
    conditions.set_ambient_temperature(srcOpPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    return operatingPoint;
}


// =====================================================================
// generate_ngspice_circuit
//
// SRC SPICE topology (HB / FB, separate Lr, FB-diode rectifier):
//
//   Vdc_supply ─┬── (bridge) ─── Cr ─── Lr ─── Lpri ─┐
//               │                                    │  (K=0.999)
//               └────────────────────── Lpri_bot ────┘  ── Lsec_o<i> ──
//
// Bridge model: behavioural Vbridge PULSE source (no SW1, no body
// diodes, no snubber) — matches Llc's BEHAVIORAL_PULSE branch and is
// adequate for FHA/PtP regression where we are interested in the tank
// shape rather than switching detail. The DC bus Vdc_supply is a
// separate 1 MΩ-loaded probe so callers can read v(vdc_supply).
//
// Secondary: one winding per output → 4-diode bridge → Cout || Rload.
// Mirrors Llc::FULL_BRIDGE rectifier; SRC has no Lm-branch dynamics
// so a single Lpri couples to each Lsec_o<i> with K=0.999 (loose enough
// to keep the K matrix positive-definite, tight enough that leakage
// is negligible vs. the explicit Lr).
// =====================================================================
std::string Src::generate_ngspice_circuit(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t inputVoltageIndex,
    size_t operatingPointIndex)
{
    // SPICE codegen covers FULL_BRIDGE_DIODE, CENTER_TAPPED_DIODE, and
    // CURRENT_DOUBLER. All three rectifier types have analytical solver
    // support (process_operating_point_for_input_voltage) since the CD
    // analytical extension. The SPICE branch here therefore no longer
    // throws on any documented rectifierType.
    SrcRectifierType rectType = get_effective_rectifier_type();
    const bool ct = (rectType == SrcRectifierType::CENTER_TAPPED_DIODE);
    const bool cd = (rectType == SrcRectifierType::CURRENT_DOUBLER);

    auto& inputVoltageSpec = get_input_voltage();
    auto& ops              = get_operating_points();

    std::vector<double> inputVoltages;
    if (inputVoltageSpec.get_nominal().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_nominal().value());
    if (inputVoltageSpec.get_minimum().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_minimum().value());
    if (inputVoltageSpec.get_maximum().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_maximum().value());
    if (inputVoltages.empty())
        throw std::runtime_error("SRC generate_ngspice_circuit: no input voltage available");

    double inputVoltage = inputVoltages[std::min(inputVoltageIndex, inputVoltages.size() - 1)];
    auto& srcOp = ops[std::min(operatingPointIndex, ops.size() - 1)];

    double fsw        = srcOp.get_switching_frequency();
    if (fsw <= 0)
        throw std::runtime_error("SRC generate_ngspice_circuit: switching frequency must be > 0");
    double period     = 1.0 / fsw;
    double halfPeriod = period / 2.0;
    // Dead-time: SRC has no MAS field for it; pick a reasonable fraction of
    // the half-period (1 %, capped at 50 ns) to give the behavioural pulse
    // a small ramp so ngspice convergence is well-behaved.
    double deadTime   = std::min(halfPeriod * 0.01, 50e-9);
    double tOn        = halfPeriod - deadTime;
    if (tOn <= 0)
        throw std::runtime_error("SRC generate_ngspice_circuit: tOn computed non-positive");

    size_t numOutputs = srcOp.get_output_voltages().size();
    if (numOutputs == 0) numOutputs = 1;

    std::vector<double> nPerOutput(numOutputs);
    for (size_t i = 0; i < numOutputs; ++i) {
        if (i < turnsRatios.size())            nPerOutput[i] = turnsRatios[i];
        else if (!turnsRatios.empty())         nPerOutput[i] = turnsRatios[0];
        else                                   nPerOutput[i] = inputVoltage / srcOp.get_output_voltages()[i];
        if (nPerOutput[i] <= 0)                nPerOutput[i] = 1.0;
    }

    double Lr = computedResonantInductance;
    double Cr = computedResonantCapacitance;
    if (Lr <= 0 || Cr <= 0)
        throw std::runtime_error("SRC generate_ngspice_circuit: tank not initialized — call process_design_requirements first");

    // Transformer Lpri: SRC has no Lm-in-the-tank, so the only role of the
    // primary inductance is to set the secondary reflection ratio. Use the
    // caller-supplied magnetizingInductance (typically 10·Lr per the design-
    // requirements default) if positive; otherwise default to 10·Lr.
    double Lpri = (magnetizingInductance > 0) ? magnetizingInductance : 10.0 * Lr;
    double k    = 0.999;

    bool isFullBridge = (get_bridge_type().has_value() &&
                         (get_bridge_type().value() == SrcBridgeType::FULL_BRIDGE ||
                          get_bridge_type().value() == SrcBridgeType::FULL_BRIDGE_PHASE_SHIFT));

    int    numPeriodsTotal = get_num_steady_state_periods() + get_num_periods_to_extract();
    double simTime         = numPeriodsTotal * period;
    double startTime       = get_num_steady_state_periods() * period;
    double maxStep         = period / 200.0;

    std::ostringstream c;
    c << "* SRC Series Resonant Converter v0.1 (Behavioural PULSE bridge, separate Lr)\n";
    c << "* " << (isFullBridge ? "Full" : "Half") << "-Bridge, full-bridge diode rectifier\n";
    c << "* Vin=" << inputVoltage << "V  fsw=" << (fsw/1e3) << "kHz  outputs=" << numOutputs << "\n";
    c << "* Lr=" << (Lr*1e6) << "uH  Cr=" << (Cr*1e9) << "nF  Lpri=" << (Lpri*1e6) << "uH\n";
    for (size_t i = 0; i < numOutputs; ++i) {
        c << "*   output " << i << ": Vout=" << srcOp.get_output_voltages()[i]
          << "V Iout=" << srcOp.get_output_currents()[i] << "A n=" << nPerOutput[i] << "\n";
    }
    c << "\n";

    // DC bus probe (constant, lets extractors read v(vdc_supply)).
    c << "Vdc_supply vdc_supply 0 " << inputVoltage << "\n";
    c << "Rdc_supply_dummy vdc_supply 0 1Meg\n\n";

    const auto bridgeMode = get_bridge_simulation_mode();

    if (bridgeMode == BridgeSimulationMode::BEHAVIORAL_PULSE) {
        // Behavioural bridge: ±Vin for FB, ±Vin/2 for HB. Fast path; no
        // SW1/body diodes/snubbers. Adequate for FHA/PtP regression where
        // the tank shape matters, not switching detail.
        if (isFullBridge) {
            c << "Vbridge bridge_a bridge_b PULSE("
              << -inputVoltage << " " << inputVoltage << " 0 "
              << std::scientific << deadTime << " " << deadTime << " "
              << tOn << " " << period << std::fixed << ")\n";
            c << "Vpri_sense bridge_a cr_in 0\n";
            c << "Vbus_gnd  bridge_b 0 0\n\n";
        }
        else {
            c << "Vbridge sw_node mid_point PULSE("
              << -(inputVoltage / 2.0) << " " << (inputVoltage / 2.0) << " 0 "
              << std::scientific << deadTime << " " << deadTime << " "
              << tOn << " " << period << std::fixed << ")\n";
            c << "Vmid       mid_point 0 0\n";
            c << "Vpri_sense sw_node   cr_in 0\n\n";
        }
    }
    else {
        // Voltage-controlled-switch bridge — high-fidelity path. Models
        // MOSFETs as SW1 + antiparallel ideal-diode + RC snubber, mirrors
        // Llc.cpp's VOLTAGE_CONTROLLED_SWITCH branch. Required by downstream
        // tooling that needs real MOSFETs to anchor (e.g. Heaviside TAS
        // decomposer: Q1..Q4 sized from the deck).
        if (deadTime <= 0.0) {
            throw std::runtime_error(
                "SRC VOLTAGE_CONTROLLED_SWITCH bridge: deadTime must be > 0. "
                "Got " + std::to_string(deadTime));
        }
        if (2.0 * deadTime >= halfPeriod) {
            throw std::runtime_error(
                "SRC VOLTAGE_CONTROLLED_SWITCH bridge: deadTime ("
                + std::to_string(deadTime) + " s) is too large for halfPeriod ("
                + std::to_string(halfPeriod) + " s) — tOn would be non-positive.");
        }

        c << ".model SW1 SW VT=2.5 VH=0.8 RON=0.01 ROFF=1Meg\n";
        c << ".model DIDEAL D(IS=1e-12 RS=0.05)\n\n";

        if (isFullBridge) {
            // Diagonal QA+QD on first half, QB+QC on second half.
            c << "Vpwm_QA pwm_QA 0 PULSE(0 5 0 10n 10n "
              << std::scientific << tOn << " " << period << std::fixed << ")\n";
            c << "Vpwm_QB pwm_QB 0 PULSE(0 5 "
              << std::scientific << halfPeriod << " 10n 10n "
              << tOn << " " << period << std::fixed << ")\n";
            c << "Vpwm_QC pwm_QC 0 PULSE(0 5 "
              << std::scientific << halfPeriod << " 10n 10n "
              << tOn << " " << period << std::fixed << ")\n";
            c << "Vpwm_QD pwm_QD 0 PULSE(0 5 0 10n 10n "
              << std::scientific << tOn << " " << period << std::fixed << ")\n";
            c << "SQA vdc_supply bridge_a pwm_QA 0 SW1\n";
            c << "DQA 0 bridge_a DIDEAL\n";
            c << "SQB bridge_a 0 pwm_QB 0 SW1\n";
            c << "DQB bridge_a vdc_supply DIDEAL\n";
            c << "Rsnub_QA vdc_supply bridge_a 1k\nCsnub_QA vdc_supply bridge_a 1n\n";
            c << "Rsnub_QB bridge_a 0 1k\nCsnub_QB bridge_a 0 1n\n";
            c << "SQC vdc_supply bridge_b pwm_QC 0 SW1\n";
            c << "DQC 0 bridge_b DIDEAL\n";
            c << "SQD bridge_b 0 pwm_QD 0 SW1\n";
            c << "DQD bridge_b vdc_supply DIDEAL\n";
            c << "Rsnub_QC vdc_supply bridge_b 1k\nCsnub_QC vdc_supply bridge_b 1n\n";
            c << "Rsnub_QD bridge_b 0 1k\nCsnub_QD bridge_b 0 1n\n";
            c << "Vpri_sense bridge_a cr_in 0\n";
            c << "Vbus_gnd  bridge_b 0 0\n\n";
        } else {
            // Half-bridge: split caps form mid_point at Vin/2.
            c << "Cbus_hi vdc_supply mid_point 1u IC=" << (inputVoltage / 2.0) << "\n";
            c << "Cbus_lo mid_point 0 1u IC=" << (inputVoltage / 2.0) << "\n";
            c << "Rbal_hi vdc_supply mid_point 100k\n";
            c << "Rbal_lo mid_point 0 100k\n";
            c << "Vpwm_HI pwm_HI 0 PULSE(0 5 0 10n 10n "
              << std::scientific << tOn << " " << period << std::fixed << ")\n";
            c << "Vpwm_LO pwm_LO 0 PULSE(0 5 "
              << std::scientific << halfPeriod << " 10n 10n "
              << tOn << " " << period << std::fixed << ")\n";
            c << "SHI vdc_supply sw_node pwm_HI 0 SW1\n";
            c << "DHI 0 sw_node DIDEAL\n";
            c << "SLO sw_node 0 pwm_LO 0 SW1\n";
            c << "DLO sw_node vdc_supply DIDEAL\n";
            c << "Rsnub_HI vdc_supply sw_node 1k\nCsnub_HI vdc_supply sw_node 1n\n";
            c << "Rsnub_LO sw_node 0 1k\nCsnub_LO sw_node 0 1n\n";
            c << "Vpri_sense sw_node cr_in 0\n\n";
        }
    }

    // Resonant tank: separate Cr then Lr, feeding Lpri.
    c << "* Resonant tank (separate Lr)\n";
    c << "Cr cr_in cr_lr  " << std::scientific << Cr << "\n";
    c << "Lr cr_lr pri_top " << std::scientific << Lr << "\n\n";

    // Transformer primary (loosely coupled to per-output secondaries).
    // Secondary windings per output depend on rectifier type:
    //   FB → 1 secondary (Lsec)              — feeds 4-diode bridge
    //   CT → 2 half-windings (Lsec1/Lsec2)   — share a center tap, feed 2 diodes
    c << "Lpri pri_top pri_bot " << std::scientific << Lpri << "\n";
    for (size_t i = 0; i < numOutputs; ++i) {
        double n_i  = nPerOutput[i];
        double Lsec = Lpri / (n_i * n_i);
        std::string si = std::to_string(i + 1);
        if (ct) {
            c << "Lsec1_o" << si << " sec_top_sec_o" << si << " sec_ct_o" << si
              << " " << std::scientific << Lsec << "\n";
            c << "Lsec2_o" << si << " sec_ct_o" << si << " sec_bot_sec_o" << si
              << " " << std::scientific << Lsec << "\n";
        }
        else {
            c << "Lsec_o" << si << " sec_pos_sec_o" << si << " sec_neg_sec_o" << si
              << " " << std::scientific << Lsec << "\n";
        }
    }
    // Pairwise K matrix (positive-definite requirement of ngspice).
    std::vector<std::string> coupled = {"Lpri"};
    coupled.reserve(1 + windings_per_output() * numOutputs);
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        if (ct) {
            coupled.push_back("Lsec1_o" + si);
            coupled.push_back("Lsec2_o" + si);
        }
        else {
            coupled.push_back("Lsec_o" + si);
        }
    }
    int kIdx = 0;
    for (size_t a = 0; a < coupled.size(); ++a)
        for (size_t b = a + 1; b < coupled.size(); ++b)
            c << "K" << ++kIdx << " " << coupled[a] << " " << coupled[b] << " " << k << "\n";
    c << "\n";

    // Primary return path
    if (isFullBridge) c << "Rpri_ret pri_bot bridge_b 0.001\n\n";
    else              c << "Rpri_ret pri_bot mid_point 0.001\n\n";

    // Secondary rectifier per output (branches by rectifier type).
    c << ".model DRECT D(Is=1e-8 N=0.01 RS=0.01)\n";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        double Vout_i = srcOp.get_output_voltages()[i];
        double Iout_i = srcOp.get_output_currents()[i];
        double Rload_i = (Iout_i > 0) ? (Vout_i / Iout_i) : 100.0;
        c << "* Output " << si << " (" << (ct ? "CT" : cd ? "CD" : "FB")
          << " diode rectifier, Vout=" << Vout_i << "V Rload=" << Rload_i << ")\n";
        if (ct) {
            // 2 diodes (anodes on each half-secondary endpoint, cathodes
            // commoned on vout_pos); center-tap returns through Vsec_sense.
            c << "D1_o" << si << " sec_top_o" << si << " vout_pos_o" << si << " DRECT\n";
            c << "D2_o" << si << " sec_bot_o" << si << " vout_pos_o" << si << " DRECT\n";
            c << "Rsn1_o" << si << " sec_top_o" << si << " vout_pos_o" << si << " 100\n";
            c << "Csn1_o" << si << " sec_top_o" << si << " vout_pos_o" << si << " 100p\n";
            c << "Rsn2_o" << si << " sec_bot_o" << si << " vout_pos_o" << si << " 100\n";
            c << "Csn2_o" << si << " sec_bot_o" << si << " vout_pos_o" << si << " 100p\n";
            // Per-half-winding current sense + common center-tap return sense.
            c << "Vsec1_sense_o" << si << " sec_top_sec_o" << si << " sec_top_o" << si << " 0\n";
            c << "Vsec2_sense_o" << si << " sec_bot_sec_o" << si << " sec_bot_o" << si << " 0\n";
            c << "Vsec_sense_o"  << si << " sec_ct_o"      << si << " vout_neg_o" << si << " 0\n";
        }
        else if (cd) {
            // CURRENT_DOUBLER (canonical topology):
            //   D1: sec_pos → vout_pos
            //   D2: sec_neg → vout_pos
            //   Lo1: sec_pos → vout_neg   (return inductor for D1's half)
            //   Lo2: sec_neg → vout_neg   (return inductor for D2's half)
            //   Cout/Rload between vout_pos and vout_neg
            //
            // Conduction trace (V(sec_pos) > V(sec_neg)):
            //   sec_pos → D1 → vout_pos → Rload → vout_neg → Lo2 → sec_neg
            //   (winding closes the loop). Lo1 freewheels through D1+Rload+Lo2.
            // Each Lo carries Iout/2 DC; ripple cancels at 2·fs at vout_pos.
            //
            // NOTE: this differs from the LLC CD SPICE pattern
            // (Llc.cpp:1567-1596) which routes Lo between diode cathodes
            // and vout_pos and adds freewheel diodes — that topology lacks
            // a working return path to sec_neg. The canonical form below
            // is the textbook current-doubler used in Telecom 3 kW briks.
            //
            // Lo sizing for SPICE: must be large enough that the inductors
            // truly clamp their currents to ±Iout/2 over a switching period
            // (CCM with small ripple). The 30 % ripple target used by
            // get_extra_components_inputs is a designer-facing recommendation
            // for the BUILT magnetic; for simulation we want headroom so the
            // current-doubler operating mode actually obtains (otherwise the
            // rectifier collapses to a plain bridge and Iout = 2·Iref/π).
            // Use ~5 % ripple per Lo here: Lo = Vo / (4·fs·0.05·IoutDC).
            const double Iout_dc = (Iout_i > 0) ? Iout_i : 1.0;
            const double fsw_local = 1.0 / period;
            const double Lo = (Vout_i > 0)
                ? (Vout_i / (4.0 * fsw_local * 0.30 * Iout_dc))
                : 10e-6;
            c << "D1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " DRECT\n";
            c << "D2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " DRECT\n";
            // Snubbers across each diode (RC).
            c << "Rsn1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100\n";
            c << "Csn1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100p\n";
            c << "Rsn2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " 100\n";
            c << "Csn2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " 100p\n";
            // Output inductors with per-Lo current sense.
            c << "VLo1_sense_o" << si << " sec_pos_o" << si << " lo1_a_o" << si << " 0\n";
            c << "Lo1_o" << si << " lo1_a_o" << si << " vout_neg_o" << si << " " << std::scientific << Lo << "\n";
            c << "VLo2_sense_o" << si << " sec_neg_o" << si << " lo2_a_o" << si << " 0\n";
            c << "Lo2_o" << si << " lo2_a_o" << si << " vout_neg_o" << si << " " << std::scientific << Lo << "\n";
            // Single bipolar secondary winding (same naming as FB).
            c << "Vsec_sense_o" << si << " sec_pos_sec_o" << si << " sec_pos_o" << si << " 0\n";
            c << "Vsec_ret_o"   << si << " sec_neg_o"    << si << " sec_neg_sec_o" << si << " 0\n";
        }
        else {
            // 4-diode bridge
            c << "Dh1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " DRECT\n";
            c << "Dh2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " DRECT\n";
            c << "Dl1_o" << si << " vout_neg_o" << si << " sec_pos_o" << si << " DRECT\n";
            c << "Dl2_o" << si << " vout_neg_o" << si << " sec_neg_o" << si << " DRECT\n";
            c << "Rsnh1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100\n";
            c << "Csnh1_o" << si << " sec_pos_o" << si << " vout_pos_o" << si << " 100p\n";
            c << "Rsnh2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " 100\n";
            c << "Csnh2_o" << si << " sec_neg_o" << si << " vout_pos_o" << si << " 100p\n";
            c << "Vsec_sense_o" << si << " sec_pos_sec_o" << si << " sec_pos_o" << si << " 0\n";
            c << "Vsec_ret_o"   << si << " sec_neg_o"    << si << " sec_neg_sec_o" << si << " 0\n";
        }
        c << "Vgnd_o"  << si << " vout_neg_o" << si << " 0 0\n";
        c << "Resr_o"  << si << " vout_pos_o" << si << " vout_cap_o" << si << " 0.05\n";
        c << "Cout_o"  << si << " vout_cap_o" << si << " vout_neg_o" << si << " " << std::scientific << 47e-6 << "\n";
        c << "Rload_o" << si << " vout_cap_o" << si << " vout_neg_o" << si << " " << Rload_i << "\n\n";
    }

    // Initial conditions
    c << ".ic";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        c << " v(vout_cap_o" << si << ")=" << srcOp.get_output_voltages()[i];
    }
    c << " v(cr_lr)=0\n\n";

    c << ".options RELTOL=0.01 ABSTOL=1e-7 VNTOL=1e-4 ITL1=500 ITL4=500\n";
    c << ".options METHOD=GEAR TRTOL=7\n\n";
    c << ".tran " << std::scientific << maxStep << " " << simTime
      << " " << startTime << " " << maxStep << " UIC\n\n";

    c << ".save v(vdc_supply) i(Vdc_supply) v(pri_top) v(pri_bot) i(Vpri_sense)";
    if (isFullBridge) c << " v(bridge_a) v(bridge_b)";
    else              c << " v(sw_node) v(mid_point)";
    for (size_t i = 0; i < numOutputs; ++i) {
        std::string si = std::to_string(i + 1);
        if (ct) {
            c << " v(sec_top_o" << si << ") v(sec_bot_o" << si << ")"
              << " v(vout_pos_o" << si << ") v(vout_cap_o" << si << ")"
              << " i(Vsec1_sense_o" << si << ") i(Vsec2_sense_o" << si << ")"
              << " i(Vsec_sense_o" << si << ")";
        }
        else if (cd) {
            c << " v(sec_pos_o" << si << ") v(sec_neg_o" << si << ")"
              << " v(vout_pos_o" << si << ") v(vout_cap_o" << si << ")"
              << " i(Vsec_sense_o" << si << ")"
              << " i(VLo1_sense_o" << si << ") i(VLo2_sense_o" << si << ")";
        }
        else {
            c << " v(sec_pos_o" << si << ") v(sec_neg_o" << si << ")"
              << " v(vout_pos_o" << si << ") v(vout_cap_o" << si << ")"
              << " i(Vsec_sense_o" << si << ")";
        }
    }
    c << "\n\n.end\n";
    return c.str();
}


// =====================================================================
// simulate_and_extract_operating_points
//
// Mirrors Llc::simulate_and_extract_operating_points. For each (Vin × OP)
// pair, generate a netlist, run ngspice, and build an OperatingPoint with
// one excitation per winding using NgspiceRunner::extract_operating_point.
// =====================================================================
std::vector<OperatingPoint> Src::simulate_and_extract_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t numberOfPeriods)
{
    NgspiceRunner runner;
    if (!runner.is_available())
        return process_operating_points(turnsRatios, magnetizingInductance);

    int numPeriods = (numberOfPeriods > 0) ? static_cast<int>(numberOfPeriods)
                                           : get_num_periods_to_extract();
    int originalNumPeriodsToExtract = get_num_periods_to_extract();
    set_num_periods_to_extract(numPeriods);

    std::vector<double> inputVoltages;
    if (get_input_voltage().get_nominal().has_value())
        inputVoltages.push_back(get_input_voltage().get_nominal().value());
    if (get_input_voltage().get_minimum().has_value())
        inputVoltages.push_back(get_input_voltage().get_minimum().value());
    if (get_input_voltage().get_maximum().has_value())
        inputVoltages.push_back(get_input_voltage().get_maximum().value());

    std::vector<OperatingPoint> result;
    auto& ops = get_operating_points();
    for (size_t vinIdx = 0; vinIdx < inputVoltages.size(); ++vinIdx) {
        for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
            std::string netlist = generate_ngspice_circuit(
                turnsRatios, magnetizingInductance, vinIdx, opIdx);
            double fsw = ops[opIdx].get_switching_frequency();

            SimulationConfig config;
            config.frequency        = fsw;
            config.extractOnePeriod = true;
            config.numberOfPeriods  = get_num_periods_to_extract();
            config.keepTempFiles    = false;

            auto simResult = runner.run_simulation(netlist, config);
            if (!simResult.success) {
                std::cerr << "SRC sim failed: " << simResult.errorMessage
                          << ". Falling back to analytical." << std::endl;
                set_num_periods_to_extract(originalNumPeriodsToExtract);
                return process_operating_points(turnsRatios, magnetizingInductance);
            }

            NgspiceRunner::WaveformNameMapping waveformMapping;
            waveformMapping.push_back({{"voltage", "pri_top"},
                                       {"current", "vpri_sense#branch"}});
            std::vector<std::string> windingNames     = {"Primary"};
            std::vector<bool>        flipCurrentSign  = {false};

            size_t numOutputs = ops[opIdx].get_output_voltages().size();
            if (numOutputs == 0) numOutputs = 1;
            SrcRectifierType rectTypeSim = get_effective_rectifier_type();
            for (size_t outIdx = 0; outIdx < numOutputs; ++outIdx) {
                std::string si = std::to_string(outIdx + 1);
                if (rectTypeSim == SrcRectifierType::CENTER_TAPPED_DIODE) {
                    waveformMapping.push_back({{"voltage", "sec_top_o" + si},
                                               {"current", "vsec1_sense_o" + si + "#branch"}});
                    windingNames.push_back("Secondary " + std::to_string(outIdx) + " Half 1");
                    flipCurrentSign.push_back(true);
                    waveformMapping.push_back({{"voltage", "sec_bot_o" + si},
                                               {"current", "vsec2_sense_o" + si + "#branch"}});
                    windingNames.push_back("Secondary " + std::to_string(outIdx) + " Half 2");
                    flipCurrentSign.push_back(true);
                }
                else {
                    waveformMapping.push_back({{"voltage", "sec_pos_o" + si},
                                               {"current", "vsec_sense_o" + si + "#branch"}});
                    windingNames.push_back("Secondary " + std::to_string(outIdx));
                    flipCurrentSign.push_back(true);
                }
            }

            OperatingPoint operatingPoint = NgspiceRunner::extract_operating_point(
                simResult, waveformMapping, fsw, windingNames,
                ops[opIdx].get_ambient_temperature(), flipCurrentSign);

            std::string vinName = (vinIdx == 0) ? "Nominal"
                                                : (vinIdx == 1 ? "Min" : "Max");
            operatingPoint.set_name(vinName + " input (" +
                std::to_string(static_cast<int>(inputVoltages[vinIdx])) + "V)");
            result.push_back(operatingPoint);
        }
    }

    set_num_periods_to_extract(originalNumPeriodsToExtract);
    return result;
}


// =====================================================================
// simulate_and_extract_topology_waveforms
//
// Mirrors Llc::simulate_and_extract_topology_waveforms. Emits the §5.1
// converter-port stream (input_voltage, input_current — reconstructed
// from power balance as the behavioural bridge does not source physical
// DC current; output_voltages from vout_cap_o<i>; output_currents from
// Vout/Rload).
// =====================================================================
std::vector<ConverterWaveforms> Src::simulate_and_extract_topology_waveforms(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t numberOfPeriods)
{
    NgspiceRunner runner;
    if (!runner.is_available())
        throw std::runtime_error("SRC simulate_and_extract_topology_waveforms: ngspice not available");

    int originalNumPeriods = get_num_periods_to_extract();
    set_num_periods_to_extract(static_cast<int>(numberOfPeriods));

    std::vector<ConverterWaveforms> results;
    auto& ops = get_operating_points();
    for (size_t opIdx = 0; opIdx < ops.size(); ++opIdx) {
        std::string netlist = generate_ngspice_circuit(
            turnsRatios, magnetizingInductance, 0, opIdx);
        double fsw = ops[opIdx].get_switching_frequency();

        SimulationConfig config;
        config.frequency        = fsw;
        config.extractOnePeriod = true;
        config.numberOfPeriods  = numberOfPeriods;
        config.keepTempFiles    = false;

        auto simResult = runner.run_simulation(netlist, config);
        if (!simResult.success)
            throw std::runtime_error("SRC simulation failed: " + simResult.errorMessage);

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
        wf.set_switching_frequency(fsw);
        wf.set_operating_point_name("SRC op. point " + std::to_string(opIdx));

        Waveform vdc = getWf("vdc_supply");
        if (!vdc.get_data().empty()) wf.set_input_voltage(vdc);

        // Power-balance reconstruction of input current (see Llc.cpp note).
        double Vin_local = 0.0;
        if (!vdc.get_data().empty()) {
            for (double v : vdc.get_data()) Vin_local += v;
            Vin_local /= vdc.get_data().size();
        } else if (get_input_voltage().get_nominal().has_value()) {
            Vin_local = get_input_voltage().get_nominal().value();
        }
        double Pout_total = 0.0;
        for (size_t i = 0; i < ops[opIdx].get_output_voltages().size(); ++i)
            Pout_total += ops[opIdx].get_output_voltages()[i] *
                          ops[opIdx].get_output_currents()[i];
        double effLocal = 1.0;
        if (get_efficiency().has_value() && get_efficiency().value() > 0.0)
            effLocal = get_efficiency().value();
        const double Iin_dc = (Vin_local > 0.0) ? Pout_total / (effLocal * Vin_local) : 0.0;
        Waveform iInWf = vdc;
        auto& iInData = iInWf.get_mutable_data();
        for (auto& v : iInData) v = Iin_dc;
        wf.set_input_current(iInWf);

        size_t numOutputs = ops[opIdx].get_output_voltages().size();
        if (numOutputs == 0) numOutputs = 1;
        for (size_t i = 0; i < numOutputs; ++i) {
            std::string si = std::to_string(i + 1);
            auto vCap = getWf("vout_cap_o" + si);
            if (vCap.get_data().empty()) continue;
            wf.get_mutable_output_voltages().push_back(vCap);
            const double Vo_i    = ops[opIdx].get_output_voltages()[i];
            const double Io_i    = ops[opIdx].get_output_currents()[i];
            const double Rload_i = (Io_i > 0) ? (Vo_i / Io_i) : 100.0;
            Waveform ioutWf      = vCap;
            auto& ioutData       = ioutWf.get_mutable_data();
            for (auto& v : ioutData) v = v / Rload_i;
            wf.get_mutable_output_currents().push_back(ioutWf);
        }
        results.push_back(wf);
    }
    set_num_periods_to_extract(originalNumPeriods);
    return results;
}


// =====================================================================
// get_extra_components_inputs
//
// SRC has TWO discrete components alongside the transformer:
//   1. resonant capacitor Cr  — emitted as a CAS::Inputs
//   2. series inductor Lr     — emitted as MAS Inputs (always external;
//                               SRC has no "integrated Lr in transformer
//                               leakage" branch — that's the LLC pattern).
//
// mode=IDEAL  → series-inductor inductance = computedResonantInductance
// mode=REAL   → recompute Cr to resonate at fr with the actual available
//               inductance (leakage if Llk >= Lr, otherwise Lr), and emit
//               the external residual Lr_ext = max(Lr - Llk, 0).
//
// Mirrors Llc::get_extra_components_inputs (separate-Ls branch) with the
// LLC's CT/CD/VD output-side extras stripped — SRC Phase-1+2 supports
// only FB-diode rectifier.
// =====================================================================
std::vector<std::variant<Inputs, CAS::Inputs>> Src::get_extra_components_inputs(
    ExtraComponentsMode mode,
    std::optional<Magnetic> magnetic)
{
    if (mode == ExtraComponentsMode::REAL && !magnetic.has_value())
        throw std::invalid_argument(
            "SRC get_extra_components_inputs: mode REAL requires a designed magnetic");

    if (computedResonantInductance <= 0 || computedResonantCapacitance <= 0 ||
        extraCapVoltageWaveforms.empty())
        throw std::runtime_error(
            "SRC get_extra_components_inputs: call process_operating_points first");

    double Lr = computedResonantInductance;
    double Cr = computedResonantCapacitance;
    double fr = computedResonantFrequency;
    if (fr <= 0)
        throw std::runtime_error("SRC get_extra_components_inputs: resonant frequency not set");

    double Llk = 0.0;
    if (mode == ExtraComponentsMode::REAL) {
        auto leakageOutput = LeakageInductance().calculate_leakage_inductance_all_windings(
            magnetic.value(), fr);
        auto perWinding = leakageOutput.get_leakage_inductance_per_winding();
        if (!perWinding.empty())
            Llk = resolve_dimensional_values(perWinding[0]);
        // Recalculate Cr to resonate with the actual available inductance.
        double Lr_eff = (Llk >= Lr) ? Llk : Lr;
        Cr = 1.0 / (4.0 * M_PI * M_PI * fr * fr * Lr_eff);
    }

    size_t nOps = extraCapVoltageWaveforms.size();
    std::vector<std::variant<Inputs, CAS::Inputs>> result;

    // -------- (1) Resonant capacitor Cr ---------
    {
        CAS::Inputs casInputs;
        CAS::DesignRequirements dr;
        CAS::DimensionWithTolerance capacitance;
        capacitance.set_nominal(Cr);
        dr.set_capacitance(capacitance);

        double peakVc = 0.0;
        for (const auto& wfm : extraCapVoltageWaveforms)
            for (double v : wfm.get_data())
                peakVc = std::max(peakVc, std::abs(v));
        dr.set_rated_voltage(peakVc * 1.2);
        dr.set_role(CAS::Application::RESONANT);
        dr.set_name("resonantCapacitor");
        casInputs.set_design_requirements(dr);

        std::vector<CAS::TwoTerminalOperatingPoint> casOps;
        for (size_t i = 0; i < nOps; ++i) {
            CAS::TwoTerminalOperatingPoint op;
            CAS::OperatingPointExcitation exc;
            exc.set_frequency(fr);

            CAS::SignalDescriptor vSig; CAS::Waveform vWfm;
            vWfm.set_data(extraCapVoltageWaveforms[i].get_data());
            if (extraCapVoltageWaveforms[i].get_time())
                vWfm.set_time(extraCapVoltageWaveforms[i].get_time().value());
            vSig.set_waveform(vWfm); exc.set_voltage(vSig);

            CAS::SignalDescriptor iSig; CAS::Waveform iWfm;
            iWfm.set_data(extraCapCurrentWaveforms[i].get_data());
            if (extraCapCurrentWaveforms[i].get_time())
                iWfm.set_time(extraCapCurrentWaveforms[i].get_time().value());
            iSig.set_waveform(iWfm); exc.set_current(iSig);

            op.set_excitation(exc);
            casOps.push_back(op);
        }
        casInputs.set_operating_points(casOps);
        result.emplace_back(std::move(casInputs));
    }

    // -------- (2) Series inductor Lr (always external in SRC) ---------
    double Lr_external = (mode == ExtraComponentsMode::IDEAL)
        ? Lr
        : ((Llk < Lr) ? (Lr - Llk) : 0.0);

    if (Lr_external > 0.0 && !extraIndVoltageWaveforms.empty()) {
        Inputs masInputs;
        DesignRequirements dr;
        DimensionWithTolerance inductance;
        inductance.set_nominal(Lr_external);
        dr.set_magnetizing_inductance(inductance);
        dr.set_name("seriesInductor");
        dr.set_topology(Topologies::SERIES_RESONANT_CONVERTER);
        dr.set_turns_ratios(std::vector<DimensionWithTolerance>{});
        dr.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY});
        masInputs.set_design_requirements(dr);

        std::vector<OperatingPoint> masOps;
        for (size_t i = 0; i < nOps; ++i) {
            OperatingPoint op;
            auto& vWfm = extraIndVoltageWaveforms[i];
            auto& iWfm = extraIndCurrentWaveforms[i];
            auto excitation = complete_excitation(iWfm, vWfm, fr, "Primary");
            op.get_mutable_excitations_per_winding().push_back(excitation);
            masOps.push_back(op);
        }
        masInputs.set_operating_points(masOps);
        result.emplace_back(std::move(masInputs));
    }

    // -------- (3) CD output filter inductors Lo1 / Lo2 ---------
    // Only emitted when rectifier is CURRENT_DOUBLER and the per-OP solver
    // populated the Lo stash. Lo sizing mirrors the LLC CD SPICE convention
    // (Llc.cpp:2070): Lo = Vout / (4·Fs·0.30·IoutDC) — 30 % ripple of IoutDC.
    if (get_effective_rectifier_type() == SrcRectifierType::CURRENT_DOUBLER &&
        !extraLoCurrentWaveforms.empty()) {
        auto& ops = get_operating_points();
        if (ops.empty())
            throw std::runtime_error("SRC CD extras: no operating points");
        double Vout_0 = ops[0].get_output_voltages()[0];
        double Iout_0 = ops[0].get_output_currents()[0];
        double fsw_0  = ops[0].get_switching_frequency();
        if (Vout_0 <= 0 || Iout_0 <= 0 || fsw_0 <= 0)
            throw std::runtime_error("SRC CD extras: invalid output/freq for Lo sizing");
        double Lo_design = Vout_0 / (4.0 * fsw_0 * 0.30 * Iout_0);

        for (size_t loIdx = 0; loIdx < 2; ++loIdx) {
            Inputs loInputs;
            DesignRequirements drLo;
            DimensionWithTolerance loInd;
            loInd.set_nominal(Lo_design);
            drLo.set_magnetizing_inductance(loInd);
            drLo.set_name(std::string("outputInductor") + (loIdx == 0 ? "1" : "2"));
            drLo.set_topology(Topologies::SERIES_RESONANT_CONVERTER);
            drLo.set_turns_ratios(std::vector<DimensionWithTolerance>{});
            drLo.set_isolation_sides(std::vector<IsolationSide>{IsolationSide::SECONDARY});
            loInputs.set_design_requirements(drLo);

            std::vector<OperatingPoint> loOps;
            auto& iVec = (loIdx == 0) ? extraLoCurrentWaveforms  : extraLo2CurrentWaveforms;
            auto& vVec = (loIdx == 0) ? extraLoVoltageWaveforms  : extraLo2VoltageWaveforms;
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

    return result;
}


// =====================================================================
// AdvancedSrc — user supplies turns ratios / Lr / Cr directly.
// =====================================================================
DesignRequirements AdvancedSrc::process_design_requirements() {
    if (desiredResonantInductance.has_value() && desiredResonantInductance.value() > 0)
        set_user_resonant_inductance(desiredResonantInductance.value());
    if (desiredResonantCapacitance.has_value() && desiredResonantCapacitance.value() > 0)
        set_user_resonant_capacitance(desiredResonantCapacitance.value());

    auto designRequirements = Src::process_design_requirements();

    if (!desiredTurnsRatios.empty()) {
        designRequirements.get_mutable_turns_ratios().clear();
        for (auto n : desiredTurnsRatios) {
            DimensionWithTolerance nTol;
            nTol.set_nominal(roundFloat(n, 2));
            designRequirements.get_mutable_turns_ratios().push_back(nTol);
        }
    }
    return designRequirements;
}

} // namespace OpenMagnetics
