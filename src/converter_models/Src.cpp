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
#include "support/Utils.h"
#include <cmath>
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
    // when M=1. So n = k_bridge · Vin / Vout.
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
    for (auto n : outputTurnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(roundFloat(n, 2));
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

    // Isolation: primary + one entry per output (full-bridge rectifier model).
    std::vector<IsolationSide> isolationSides;
    isolationSides.push_back(get_isolation_side_from_index(0));
    for (size_t i = 0; i < nOutputs; ++i)
        isolationSides.push_back(get_isolation_side_from_index(i + 1));
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
    auto get_n_for_output = [&](size_t outputIdx) -> double {
        if (turnsRatios.empty())
            return (k_bridge * inputVoltage) / srcOpPoint.get_output_voltages()[outputIdx];
        if (outputIdx < turnsRatios.size()) return turnsRatios[outputIdx];
        return turnsRatios[0];
    };

    double n_main = get_n_for_output(0);
    if (n_main <= 0)
        throw std::runtime_error("SRC: turns ratio must be > 0");

    double Vout_main = srcOpPoint.get_output_voltages()[0];
    double Iout_main = srcOpPoint.get_output_currents()[0];
    double Rload    = Vout_main / Iout_main;
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
    double Vbridge_lvl = k_bridge * inputVoltage;   // half-bridge: ±Vin/2; full-bridge: ±Vin
    for (int k = 0; k < totalSamples; ++k) {
        double t = time_full[k];
        double theta = w * t;                       // bridge half-cycle ref (high half: 0..π)
        ILr_full[k]  = ILr_pk * std::sin(theta - phi);
        // Bridge midpoint voltage referred to bridge-zero-DC frame:
        //   half 0 (0 ≤ θ mod 2π < π) → +Vbridge_lvl
        //   half 1 (π ≤ θ mod 2π < 2π) → -Vbridge_lvl
        double thetaMod = std::fmod(theta, 2.0 * M_PI);
        if (thetaMod < 0) thetaMod += 2.0 * M_PI;
        Vpri_full[k] = (thetaMod < M_PI) ? +Vbridge_lvl : -Vbridge_lvl;
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

    // ── Secondary windings (one per output; full-bridge rectifier model) ──
    //
    // Phase 2 supports only full-bridge diode rectifier on the secondary
    // (rectifierType=fullBridgeDiode, the schema default). Center-tapped and
    // current-doubler are deferred to Phase 3.
    auto rt = get_rectifier_type();
    if (rt.has_value() && rt.value() != SrcRectifierType::FULL_BRIDGE_DIODE) {
        throw std::runtime_error(
            "SRC: only fullBridgeDiode rectifier is supported in Phase 2 "
            "(centerTappedDiode / currentDoubler deferred to Phase 3).");
    }

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

    OperatingConditions conditions;
    conditions.set_ambient_temperature(srcOpPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    return operatingPoint;
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
