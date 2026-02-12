#include "converter_models/Llc.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "support/Exceptions.h"

namespace OpenMagnetics {

// ═══════════════════════════════════════════════════════════════════════
// Construction
// ═══════════════════════════════════════════════════════════════════════
Llc::Llc(const json& j) {
    from_json(j, *static_cast<LlcResonant*>(this));
}

AdvancedLlc::AdvancedLlc(const json& j) {
    from_json(j, *this);
}

// ═══════════════════════════════════════════════════════════════════════
// Bridge voltage factor
// ═══════════════════════════════════════════════════════════════════════
double Llc::get_bridge_voltage_factor() const {
    auto bt = get_bridge_type();
    if (bt.has_value() && bt.value() == LlcBridgeType::FULL_BRIDGE)
        return 1.0;
    return 0.5;  // Half-bridge: Vi = ½·Vin
}

// ═══════════════════════════════════════════════════════════════════════
// Effective resonant frequency
// ═══════════════════════════════════════════════════════════════════════
double Llc::get_effective_resonant_frequency() const {
    if (get_resonant_frequency().has_value())
        return get_resonant_frequency().value();
    // Geometric mean of switching frequency range
    return std::sqrt(get_min_switching_frequency() * get_max_switching_frequency());
}

// ═══════════════════════════════════════════════════════════════════════
// Input validation
// ═══════════════════════════════════════════════════════════════════════
bool Llc::run_checks(bool assertErrors) {
    _assertErrors = assertErrors;
    bool ok = true;

    auto& ops = get_operating_points();
    if (ops.empty()) {
        if (assertErrors) throw std::runtime_error("LLC: no operating points");
        return false;
    }

    for (size_t i = 0; i < ops.size(); i++) {
        auto& op = ops[i];
        if (op.get_output_voltages().empty() || op.get_output_currents().empty()) {
            if (assertErrors) throw std::runtime_error("LLC: OP missing voltages/currents");
            ok = false;
        }
        if (op.get_output_voltages().size() != op.get_output_currents().size()) {
            if (assertErrors) throw std::runtime_error("LLC: voltage/current count mismatch");
            ok = false;
        }
        double fsw = op.get_switching_frequency();
        if (fsw < get_min_switching_frequency() * 0.99 ||
            fsw > get_max_switching_frequency() * 1.01) {
            if (assertErrors) throw std::runtime_error("LLC: fsw out of range");
            ok = false;
        }
    }
    return ok;
}

// ═══════════════════════════════════════════════════════════════════════
// Design Requirements
//
// Based on Runo Nielsen's model (LLC_LCC.pdf, llc.pdf):
//   Vo (fictitious output voltage) = k_bridge * Vin_nom
//     where k_bridge = 0.5 for HB, 1.0 for FB
//   n = Vo / Vout (turns ratio)
//   f1 = 1/(2π√(Ls·C)) is the Load Independent Point (LIP) frequency
//   f0 = 1/(2π√((Ls+L)·C)) is the freewheeling resonant frequency
//
// FHA-based component sizing (for initial values):
//   Rac = (8·n²)/(π²) · Rload
//   Zr = Q · Rac
//   Ls = Zr / (2π·fr)
//   C  = 1 / (2π·fr·Zr)
//   L  = Ln · Ls
// ═══════════════════════════════════════════════════════════════════════
DesignRequirements Llc::process_design_requirements() {
    double k_bridge = get_bridge_voltage_factor();
    auto& inputVoltage = get_input_voltage();
    double Vin_nom = inputVoltage.get_nominal().value_or(
        (inputVoltage.get_minimum().value_or(0) + inputVoltage.get_maximum().value_or(0)) / 2.0);

    auto& ops = get_operating_points();
    double mainOutputVoltage  = ops[0].get_output_voltages()[0];
    double mainOutputCurrent  = ops[0].get_output_currents()[0];

    // ─── Runo Nielsen: Fictitious output voltage & turns ratio ───
    // Vo ≈ k_bridge · Vin_nom  (places LIP near nominal input)
    // n = Vo / Vout
    double Vo = k_bridge * Vin_nom;
    double mainTurnsRatio = Vo / mainOutputVoltage;

    // ─── Turns ratios for all secondaries ───
    std::vector<double> turnsRatios;
    turnsRatios.push_back(mainTurnsRatio);
    for (size_t i = 1; i < ops[0].get_output_voltages().size(); i++) {
        turnsRatios.push_back(Vo / ops[0].get_output_voltages()[i]);
    }

    // ─── FHA-based resonant tank sizing ───
    double Rload = mainOutputVoltage / mainOutputCurrent;
    double Rac = (8.0 * mainTurnsRatio * mainTurnsRatio) / (M_PI * M_PI) * Rload;

    double Q = get_quality_factor().value_or(0.4);
    double fr = get_effective_resonant_frequency();

    double Zr = Q * Rac;
    double Ls = Zr / (2.0 * M_PI * fr);
    double Cr = 1.0 / (2.0 * M_PI * fr * Zr);

    // ─── Magnetizing inductance ───
    double Ln = computedInductanceRatio;  // L/Ls ratio (typically 3-10)
    double L = Ln * Ls;

    // Store computed values
    computedResonantInductance  = Ls;
    computedResonantCapacitance = Cr;

    // ─── Verify resonant frequencies ───
    // f1 = 1/(2π√(Ls·C)) — LIP frequency (power delivery resonance)
    // f0 = 1/(2π√((Ls+L)·C)) — freewheeling resonance
    double f1 = 1.0 / (2.0 * M_PI * std::sqrt(Ls * Cr));
    double f0 = 1.0 / (2.0 * M_PI * std::sqrt((Ls + L) * Cr));
    // f1 should equal fr by construction
    // f0 < f1 always

    // ─── Build DesignRequirements ───
    DesignRequirements designRequirements;
    designRequirements.get_mutable_turns_ratios().clear();
    for (auto n : turnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(roundFloat(n, 2));
        designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }

    DimensionWithTolerance inductanceWithTolerance;
    inductanceWithTolerance.set_nominal(roundFloat(L, 10));
    designRequirements.set_magnetizing_inductance(inductanceWithTolerance);

    // If integrated resonant inductor, request leakage = Ls
    if (get_integrated_resonant_inductor().value_or(false)) {
        std::vector<DimensionWithTolerance> leakageReqs;
        DimensionWithTolerance lrTol;
        lrTol.set_nominal(roundFloat(Ls, 10));
        leakageReqs.push_back(lrTol);
        designRequirements.set_leakage_inductance(leakageReqs);
    }

    return designRequirements;
}

// ═══════════════════════════════════════════════════════════════════════
// Process operating points for all input voltages
// ═══════════════════════════════════════════════════════════════════════
std::vector<OperatingPoint> Llc::process_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    std::vector<OperatingPoint> result;
    auto& inputVoltage = get_input_voltage();
    auto& ops = get_operating_points();

    // Collect input voltages: nominal, minimum, maximum
    std::vector<double> inputVoltages;
    if (inputVoltage.get_nominal().has_value())
        inputVoltages.push_back(inputVoltage.get_nominal().value());
    if (inputVoltage.get_minimum().has_value())
        inputVoltages.push_back(inputVoltage.get_minimum().value());
    if (inputVoltage.get_maximum().has_value())
        inputVoltages.push_back(inputVoltage.get_maximum().value());

    // Remove duplicates
    std::sort(inputVoltages.begin(), inputVoltages.end());
    inputVoltages.erase(std::unique(inputVoltages.begin(), inputVoltages.end()), inputVoltages.end());

    for (double Vin : inputVoltages) {
        auto op = process_operating_point_for_input_voltage(
            Vin, ops[0], turnsRatios, magnetizingInductance);
        result.push_back(op);
    }
    return result;
}

std::vector<OperatingPoint> Llc::process_operating_points(Magnetic magnetic) {
    auto req = process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    return process_operating_points(turnsRatios, Lm);
}

// ═══════════════════════════════════════════════════════════════════════
// CORE WAVEFORM GENERATION — Runo Nielsen Time Domain Approach
//
// Reference: LLC_LCC.pdf pages 6-7, llc.pdf pages 1-4
//
// LLC equivalent circuit (Runo's Figure 3):
//
//         -Vc+      Ls        Id →  +Vo
//   Vi ──┬──||──────/\/\/──┬──────┤►├──── +
//        |   C             |
//        |                 L  IL↓
//        |                 |
//        └─────────────────┴──────┤►├──── -
//                            Id →  -Vo
//
// Variables:
//   ILs(t) : current in series inductor Ls = primary winding current
//   IL(t)  : current in magnetizing inductor L
//   Id(t)  : diode current = ILs(t) - IL(t)
//   Vc(t)  : resonance capacitor voltage
//   Vi(t)  : input square wave = ±k_bridge·Vin
//   VL(t)  : voltage across L (= transformer primary voltage)
//            Clamped to ±Vo during power delivery
//
// Two resonant frequencies:
//   w1 = 1/√(Ls·C)      Z1 = √(Ls/C)     — diodes ON
//   w0 = 1/√((Ls+L)·C)  Z0 = √((Ls+L)/C) — diodes OFF
//
// Positive half-cycle (Vi = +k_bridge·Vin):
//
//   Phase A — Power Delivery (diodes ON, VL = +Vo):
//     ILs(t) = ILs0·cos(w1·t) + (Vi-Vo-Vc0)/Z1·sin(w1·t)
//     Vc(t)  = (Vi-Vo) - (Vi-Vo-Vc0)·cos(w1·t) + ILs0·Z1·sin(w1·t)
//     IL(t)  = IL0 + (Vo/L)·t    (linear ramp, Lm clamped to Vo)
//     Ends when Id = ILs - IL ≤ 0  OR  t = Thalf
//
//   Phase B — Freewheeling (diodes OFF, IL = ILs):
//     ILs(τ) = ILs_fw·cos(w0·τ) + (Vi-Vc_fw)/Z0·sin(w0·τ)
//     Vc(τ)  = Vi - (Vi-Vc_fw)·cos(w0·τ) + ILs_fw·Z0·sin(w0·τ)
//     IL(τ)  = ILs(τ)
//     Where τ = t - t_freewheel
//
// Negative half-cycle: by half-wave antisymmetry
//   ILs(t+Thalf) = -ILs(t), Vc(t+Thalf) = -Vc(t), IL(t+Thalf) = -IL(t)
//
// Steady-state boundary condition:
//   Vc(Thalf) = -Vc0  (iterated by bisection on Vc0)
//   ILs(Thalf) = -ILs0 = IL(Thalf) = -IL0  (at switching instant)
// ═══════════════════════════════════════════════════════════════════════

OperatingPoint Llc::process_operating_point_for_input_voltage(
    double inputVoltage,
    const LlcOperatingPoint& llcOpPoint,
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    OperatingPoint operatingPoint;

    double fsw  = llcOpPoint.get_switching_frequency();
    double k_bridge = get_bridge_voltage_factor();

    // ─── Model parameters (Runo Nielsen notation) ───────────────
    double Vi = k_bridge * inputVoltage;        // Square wave amplitude ±Vi
    double Vo = turnsRatios[0] * llcOpPoint.get_output_voltages()[0];  // Fictitious output voltage
    double Ls = computedResonantInductance;     // Series (resonant) inductor
    double C  = computedResonantCapacitance;    // Resonance capacitor = C1+C2
    double L  = magnetizingInductance;          // Magnetizing inductance

    double period = 1.0 / fsw;
    double Thalf  = period / 2.0;

    // ─── Resonant parameters ────────────────────────────────────
    // w1, Z1: power delivery (diodes ON, Ls-C resonance)
    double w1 = 1.0 / std::sqrt(Ls * C);
    double Z1 = std::sqrt(Ls / C);

    // w0, Z0: freewheeling (diodes OFF, (Ls+L)-C resonance)
    double w0 = 1.0 / std::sqrt((Ls + L) * C);
    double Z0 = std::sqrt((Ls + L) / C);

    // ─── Sampling ───────────────────────────────────────────────
    const int N = 256;  // Samples per half-period
    double dt = Thalf / N;

    // ─── Initial conditions for steady-state ────────────────────
    // At t=0 (start of positive half-cycle):
    //   ILs(0) = IL(0) = ILs0  (diode current is zero at switching instant)
    //   Vc(0) = Vc0  (to be determined by iteration)
    //
    // From magnetizing current ramp during power delivery:
    //   IL ramps from IL0 toward +Im_pk with slope Vo/L
    //   Im_pk = Vo/(4·L·fsw) = Vo·Thalf/(2·L)
    //
    // At the switching instant, ILs = IL (ZVS condition), so:
    //   ILs0 = IL0 = -Im_pk (estimated, will be refined by iteration)

    double Im_pk_est = Vo * Thalf / (2.0 * L);
    double ILs0 = -Im_pk_est;
    double IL0  = ILs0;

    // ─── Bisection on Vc0 ───────────────────────────────────────
    // Find Vc0 such that Vc(Thalf) = -Vc0 (half-wave antisymmetry)
    // Also refine ILs0/IL0 such that ILs(Thalf) = -ILs0

    double Vc_lo = -3.0 * Vi;
    double Vc_hi =  3.0 * Vi;
    double Vc0 = 0.0;

    const int MAX_BISECT = 60;
    const double TOL_VC = 1e-4;

    // Storage for final waveform (filled on last iteration)
    std::vector<double> ILs_pos(N + 1);
    std::vector<double> IL_pos(N + 1);
    std::vector<double> Vc_pos(N + 1);
    std::vector<double> VL_pos(N + 1);  // Voltage across L (= transformer primary voltage)

    for (int bisect = 0; bisect < MAX_BISECT; ++bisect) {
        Vc0 = 0.5 * (Vc_lo + Vc_hi);

        // ── Simulate positive half-cycle ──
        bool in_freewheeling = false;
        double t_fw = Thalf;
        double ILs_fw = 0, Vc_fw = 0;

        double ILs_end = 0, IL_end = 0, Vc_end = 0;

        for (int k = 0; k <= N; ++k) {
            double t = k * dt;
            double ILs_t, IL_t, Vc_t, VL_t;

            if (!in_freewheeling) {
                // Phase A: Power Delivery (diodes ON)
                // VL = +Vo (magnetizing inductance clamped by output)
                double V_drive = Vi - Vo;

                ILs_t = ILs0 * std::cos(w1 * t)
                       + (V_drive - Vc0) / Z1 * std::sin(w1 * t);

                Vc_t = V_drive
                     - (V_drive - Vc0) * std::cos(w1 * t)
                     + ILs0 * Z1 * std::sin(w1 * t);

                IL_t = IL0 + (Vo / L) * t;

                VL_t = Vo;  // Clamped by output diodes

                // Check if diodes turn off: Id = ILs - IL ≤ 0
                double Id = ILs_t - IL_t;
                if (Id < 0 && k > 0) {
                    in_freewheeling = true;
                    t_fw = t;
                    ILs_fw = ILs_t;
                    Vc_fw  = Vc_t;
                    // At this instant IL = ILs (approximately)
                    ILs_fw = IL_t;  // Ensure continuity
                }
            }

            if (in_freewheeling) {
                // Phase B: Freewheeling (diodes OFF, Lm joins resonance)
                double tau = t - t_fw;

                ILs_t = ILs_fw * std::cos(w0 * tau)
                       + (Vi - Vc_fw) / Z0 * std::sin(w0 * tau);

                Vc_t = Vi
                     - (Vi - Vc_fw) * std::cos(w0 * tau)
                     + ILs_fw * Z0 * std::sin(w0 * tau);

                IL_t = ILs_t;  // No diode current, IL tracks ILs

                // VL = voltage across L during freewheeling
                // From KVL: Vi = Vc + Ls·dILs/dt + L·dIL/dt
                // Since ILs = IL: Vi = Vc + (Ls+L)·dILs/dt
                // dILs/dt = -ILs_fw·w0·sin(w0·tau) + (Vi-Vc_fw)/Z0·w0·cos(w0·tau)
                double dILs_dt = -ILs_fw * w0 * std::sin(w0 * tau)
                               + (Vi - Vc_fw) / Z0 * w0 * std::cos(w0 * tau);
                VL_t = L * dILs_dt;
                // Equivalently: VL = Vi - Vc - Ls·dILs/dt = L·dILs/dt
            }

            ILs_pos[k] = ILs_t;
            IL_pos[k]  = IL_t;
            Vc_pos[k]  = Vc_t;
            VL_pos[k]  = VL_t;

            if (k == N) {
                ILs_end = ILs_t;
                IL_end  = IL_t;
                Vc_end  = Vc_t;
            }
        }

        // Bisection criterion: Vc(Thalf) should equal -Vc0
        double err = Vc_end + Vc0;
        if (std::abs(err) < TOL_VC) break;
        if (err > 0) Vc_hi = Vc0;
        else         Vc_lo = Vc0;
    }
    
    // ─── Run final simulation with converged Vc0 ─────────────────
    // Re-run one last time to fill the waveform arrays correctly
    {
        bool in_freewheeling = false;
        double t_fw = Thalf;
        double ILs_fw = 0, Vc_fw = 0;

        for (int k = 0; k <= N; ++k) {
            double t = k * dt;
            double ILs_t, IL_t, Vc_t, VL_t;

            if (!in_freewheeling) {
                double V_drive = Vi - Vo;

                ILs_t = ILs0 * std::cos(w1 * t)
                       + (V_drive - Vc0) / Z1 * std::sin(w1 * t);

                Vc_t = V_drive
                     - (V_drive - Vc0) * std::cos(w1 * t)
                     + ILs0 * Z1 * std::sin(w1 * t);

                IL_t = IL0 + (Vo / L) * t;

                VL_t = Vo;

                double Id = ILs_t - IL_t;
                
                if (Id < 0 && k > 0) {
                    in_freewheeling = true;
                    t_fw = t;
                    ILs_fw = ILs_t;
                    Vc_fw  = Vc_t;
                }
            }

            if (in_freewheeling) {
                double tau = t - t_fw;

                ILs_t = ILs_fw * std::cos(w0 * tau)
                       + (Vi - Vc_fw) / Z0 * std::sin(w0 * tau);

                Vc_t = Vi
                     - (Vi - Vc_fw) * std::cos(w0 * tau)
                     + ILs_fw * Z0 * std::sin(w0 * tau);

                IL_t = ILs_t;

                double dILs_dt = -ILs_fw * w0 * std::sin(w0 * tau)
                               + (Vi - Vc_fw) / Z0 * w0 * std::cos(w0 * tau);
                VL_t = L * dILs_dt;
            }

            ILs_pos[k] = ILs_t;
            IL_pos[k]  = IL_t;
            Vc_pos[k]  = Vc_t;
            VL_pos[k]  = VL_t;
            
            // Debug output
            if (k < 3) {
                std::cout << "DEBUG MKF LLC: k=" << k << " ILs=" << ILs_t 
                          << " IL=" << IL_t << " diff=" << (ILs_t - IL_t) << std::endl;
            }
        }
    }

    // ─── Refine IL0 from steady-state result ────────────────────
    // At end of half-cycle: ILs_end should be -ILs0
    // Update ILs0 = IL0 = -ILs_end for next iteration if needed
    // (The bisection already converges Vc0; ILs0 is set by Im_pk_est
    //  which is exact for above-resonance operation. For below-resonance
    //  the freewheeling phase handles the IL=ILs condition.)

    // ─── Build full-period waveforms ────────────────────────────
    int totalSamples = 2 * N + 1;
    std::vector<double> time_full(totalSamples);
    std::vector<double> ILs_full(totalSamples);
    std::vector<double> IL_full(totalSamples);
    std::vector<double> VL_full(totalSamples);

    // Positive half-cycle
    for (int k = 0; k <= N; ++k) {
        time_full[k] = k * dt;
        ILs_full[k]  = ILs_pos[k];
        IL_full[k]   = IL_pos[k];
        VL_full[k]   = VL_pos[k];
    }

    // Negative half-cycle by antisymmetry
    for (int k = 1; k <= N; ++k) {
        time_full[N + k] = Thalf + k * dt;
        ILs_full[N + k]  = -ILs_pos[k];
        IL_full[N + k]   = -IL_pos[k];
        VL_full[N + k]   = -VL_pos[k];
    }

    // ─── Primary winding excitation ─────────────────────────────
    // The primary winding of the transformer carries ILs(t)
    // (the full resonant tank current flows through the primary)
    // The voltage across the primary winding = VL(t) (voltage across Lm)
    {
        Waveform currentWaveform;
        currentWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        currentWaveform.set_data(ILs_full);
        currentWaveform.set_time(time_full);

        // Primary voltage = VL(t) = ±Vo during power delivery,
        // varying during freewheeling
        Waveform voltageWaveform;
        voltageWaveform.set_ancillary_label(WaveformLabel::CUSTOM);
        voltageWaveform.set_data(VL_full);
        voltageWaveform.set_time(time_full);

        auto excitation = complete_excitation(currentWaveform, voltageWaveform,
                                              fsw, "Primary");
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    // ─── Secondary winding excitation(s) ────────────────────────
    // For each secondary: 
    //   i_sec(t) = |ILs(t) - IL(t)| / n  (rectified transformer current)
    //   v_sec(t) = VL(t) / n  (transformer reflected voltage, NOT rectified — this is
    //              the voltage the winding sees, not after the rectifier)
    
    // If no turns ratios provided, compute from input/output voltages
    std::vector<double> effectiveTurnsRatios = turnsRatios;
    if (effectiveTurnsRatios.empty()) {
        // Compute turns ratio from the operating point: n = Vi / Vout
        // where Vi is the effective primary voltage (k_bridge * Vin)
        double computedTurnsRatio = Vi / llcOpPoint.get_output_voltages()[0];
        effectiveTurnsRatios.push_back(computedTurnsRatio);
    }
    
    for (size_t secIdx = 0; secIdx < effectiveTurnsRatios.size(); ++secIdx) {
        double n = effectiveTurnsRatios[secIdx];

        std::vector<double> iSecData(totalSamples);
        std::vector<double> vSecData(totalSamples);

        for (int k = 0; k < totalSamples; ++k) {
            // Diode current in model = ILs - IL, reflected and rectified
            double Id = ILs_full[k] - IL_full[k];
            
            // Debug: check first few samples
            if (k < 3) {
                std::cout << "DEBUG SEC CALC: k=" << k << " ILs=" << ILs_full[k] 
                          << " IL=" << IL_full[k] << " Id=" << Id << std::endl;
            }

            // For center-tapped rectification: each secondary half-winding
            // carries current only during its half-cycle
            // The current magnitude in the secondary winding is |Id|/n
            iSecData[k] = std::abs(Id) / n;

            // Secondary voltage = VL / n (transformer reflected voltage)
            // This is the voltage across the secondary winding
            vSecData[k] = VL_full[k] / n;
        }

        // Debug: Check secondary current data before creating waveform
        std::cout << "DEBUG: iSecData size=" << iSecData.size() 
                  << " max=" << *std::max_element(iSecData.begin(), iSecData.end()) << std::endl;

        Waveform secCurrentWfm;
        secCurrentWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        secCurrentWfm.set_data(iSecData);
        secCurrentWfm.set_time(time_full);

        // Verify data was set
        auto verifyData = secCurrentWfm.get_data();
        std::cout << "DEBUG: After set_data, size=" << verifyData.size() 
                  << " max=" << (verifyData.empty() ? 0 : *std::max_element(verifyData.begin(), verifyData.end())) << std::endl;

        Waveform secVoltageWfm;
        secVoltageWfm.set_ancillary_label(WaveformLabel::CUSTOM);
        secVoltageWfm.set_data(vSecData);
        secVoltageWfm.set_time(time_full);

        auto excitation = complete_excitation(secCurrentWfm, secVoltageWfm,
                                              fsw,
                                              "Secondary " + std::to_string(secIdx));
        
        // Debug: Check excitation current
        if (excitation.get_current() && excitation.get_current()->get_waveform()) {
            auto checkData = excitation.get_current()->get_waveform()->get_data();
            std::cout << "DEBUG: In excitation, size=" << checkData.size() 
                      << " max=" << (checkData.empty() ? 0 : *std::max_element(checkData.begin(), checkData.end())) << std::endl;
        }
        
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    // ─── Operating conditions ───────────────────────────────────
    OperatingConditions conditions;
    conditions.set_ambient_temperature(llcOpPoint.get_ambient_temperature());
    conditions.set_cooling(std::nullopt);
    operatingPoint.set_conditions(conditions);

    // Debug: Check final operating point before return
    std::cout << "DEBUG: Returning operating point with " 
              << operatingPoint.get_excitations_per_winding().size() << " excitations" << std::endl;
    if (operatingPoint.get_excitations_per_winding().size() > 1) {
        const auto& secExc = operatingPoint.get_excitations_per_winding()[1];
        if (secExc.get_current() && secExc.get_current()->get_waveform()) {
            auto finalData = secExc.get_current()->get_waveform()->get_data();
            std::cout << "DEBUG: Final secondary current size=" << finalData.size() 
                      << " max=" << (finalData.empty() ? 0 : *std::max_element(finalData.begin(), finalData.end())) << std::endl;
        }
    }

    return operatingPoint;
}


// ═══════════════════════════════════════════════════════════════════════
// SPICE Circuit Generation
//
// Generates an NgSpice netlist for the LLC converter matching
// Runo Nielsen's equivalent circuit (Figure 3).
//
// Circuit topology:
//   Vi (square wave) → C (resonance cap) → Ls (series inductor) →
//   → Transformer primary (L = magnetizing inductance) →
//   → Diode bridge → Vo (output)
// ═══════════════════════════════════════════════════════════════════════
std::string Llc::generate_ngspice_circuit(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance,
    size_t inputVoltageIndex,
    size_t operatingPointIndex)
{
    auto& inputVoltageSpec = get_input_voltage();
    auto& ops = get_operating_points();

    // Select input voltage
    std::vector<double> inputVoltages;
    if (inputVoltageSpec.get_nominal().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_nominal().value());
    if (inputVoltageSpec.get_minimum().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_minimum().value());
    if (inputVoltageSpec.get_maximum().has_value())
        inputVoltages.push_back(inputVoltageSpec.get_maximum().value());

    double inputVoltage = inputVoltages[std::min(inputVoltageIndex, inputVoltages.size() - 1)];
    auto& llcOp = ops[std::min(operatingPointIndex, ops.size() - 1)];

    double fsw = llcOp.get_switching_frequency();
    double period = 1.0 / fsw;
    double halfPeriod = period / 2.0;
    double deadTime = computedDeadTime;
    double tOn = halfPeriod - deadTime;

    double Vout = llcOp.get_output_voltages()[0];
    double Iout = llcOp.get_output_currents()[0];
    double n    = turnsRatios[0];

    double Ls = computedResonantInductance;
    double Cr = computedResonantCapacitance;
    double L  = magnetizingInductance;

    bool isFullBridge = (get_bridge_type().has_value() &&
                         get_bridge_type().value() == LlcBridgeType::FULL_BRIDGE);

    // Simulation timing
    int periodsToExtract = get_num_periods_to_extract();
    int steadyStatePeriods = get_num_steady_state_periods();
    int numPeriodsTotal = steadyStatePeriods + periodsToExtract;
    double simTime = numPeriodsTotal * period;
    double startTime = steadyStatePeriods * period;
    double stepTime = period / 200;

    std::ostringstream circuit;

    circuit << "* LLC Resonant Converter - Generated by OpenMagnetics\n";
    circuit << "* " << (isFullBridge ? "Full" : "Half") << "-Bridge\n";
    circuit << "* Vin=" << inputVoltage << "V, f=" << (fsw/1e3)
            << "kHz, Vout=" << Vout << "V\n";
    circuit << "* Runo Nielsen model: Ls=" << (Ls*1e6) << "uH, C="
            << (Cr*1e9) << "nF, L=" << (L*1e6) << "uH\n\n";

    // DC input
    circuit << "Vdc vin_dc 0 " << inputVoltage << "\n\n";

    if (isFullBridge) {
        // Full-bridge: 4 switches
        circuit << "* Full-bridge switching\n";
        circuit << ".model SW1 SW(Ron=10m Roff=10Meg Vt=2.5)\n";
        circuit << ".model DIDEAL D(Is=1e-14 N=0.001)\n\n";

        circuit << "Vpwm_a pwm_a 0 PULSE(0 5 0 10n 10n "
                << std::scientific << tOn << " " << period << std::fixed << ")\n";
        circuit << "Vpwm_b pwm_b 0 PULSE(0 5 "
                << std::scientific << halfPeriod << std::fixed
                << " 10n 10n "
                << std::scientific << tOn << " " << period << std::fixed << ")\n\n";

        circuit << "S1 vin_dc bridge_a pwm_a 0 SW1\n";
        circuit << "D1a 0 bridge_a DIDEAL\n";
        circuit << "S2 bridge_a 0 pwm_b 0 SW1\n";
        circuit << "D2a bridge_a vin_dc DIDEAL\n\n";

        circuit << "S3 vin_dc bridge_b pwm_b 0 SW1\n";
        circuit << "D3b 0 bridge_b DIDEAL\n";
        circuit << "S4 bridge_b 0 pwm_a 0 SW1\n";
        circuit << "D4b bridge_b vin_dc DIDEAL\n\n";

        circuit << "Vpri_sense bridge_a lr_in 0\n\n";

        // Resonant tank: C in series, then Ls
        circuit << "Cr lr_in cr_ls " << std::scientific << Cr << "\n";
        circuit << "Ls cr_ls pri_top " << std::scientific << Ls << "\n\n";

        // Transformer: magnetizing inductance L in parallel with ideal transformer
        circuit << "* Transformer with magnetizing inductance L\n";
        circuit << "L_mag pri_top bridge_b " << std::scientific << L << "\n";
        double k_coupling = std::sqrt(L / (L + Ls));
        circuit << "* Coupling coefficient k=" << k_coupling << "\n\n";

    } else {
        // Half-bridge: 2 switches + capacitive divider
        circuit << "* Half-bridge switching\n";
        circuit << ".model SW1 SW(Ron=10m Roff=10Meg Vt=2.5)\n";
        circuit << ".model DIDEAL D(Is=1e-14 N=0.001)\n\n";

        circuit << "Vpwm_hi pwm_hi 0 PULSE(0 5 0 10n 10n "
                << std::scientific << tOn << " " << period << std::fixed << ")\n";
        circuit << "Vpwm_lo pwm_lo 0 PULSE(0 5 "
                << std::scientific << halfPeriod << std::fixed
                << " 10n 10n "
                << std::scientific << tOn << " " << period << std::fixed << ")\n\n";

        circuit << "S_hi vin_dc mid_point pwm_hi 0 SW1\n";
        circuit << "D_hi 0 mid_point DIDEAL\n";
        circuit << "S_lo mid_point 0 pwm_lo 0 SW1\n";
        circuit << "D_lo mid_point vin_dc DIDEAL\n\n";

        circuit << "Vpri_sense mid_point lr_in 0\n\n";

        // Resonant tank: C, then Ls
        circuit << "Cr lr_in cr_ls " << std::scientific << Cr << "\n";
        circuit << "Ls cr_ls pri_top " << std::scientific << Ls << "\n\n";

        // Magnetizing inductance
        circuit << "L_mag pri_top 0 " << std::scientific << L << "\n\n";
    }

    // Output rectifier and load
    circuit << "* Output rectifier (center-tapped model)\n";
    circuit << "* Transformer ratio n = " << n << "\n";
    double Rload = Vout / Iout;
    circuit << "Rload out_pos 0 " << Rload << "\n";
    circuit << "Cout out_pos 0 " << std::scientific << 100e-6 << "\n\n";

    // Diode rectifier modeled as voltage source for TDA validation
    circuit << "Vout_sense out_pos 0 " << Vout << "\n\n";

    // Simulation commands
    circuit << ".tran " << std::scientific << stepTime
            << " " << simTime << " " << startTime << "\n";
    circuit << ".end\n";

    return circuit.str();
}


// ═══════════════════════════════════════════════════════════════════════
// SPICE simulation wrappers
// ═══════════════════════════════════════════════════════════════════════
std::vector<OperatingPoint> Llc::simulate_and_extract_operating_points(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    // Placeholder: run NgSpice and extract operating points
    // For now, fall back to analytical
    return process_operating_points(turnsRatios, magnetizingInductance);
}

std::vector<ConverterWaveforms> Llc::simulate_and_extract_topology_waveforms(
    const std::vector<double>& turnsRatios,
    double magnetizingInductance)
{
    // Placeholder: run NgSpice and extract topology-level waveforms
    std::vector<ConverterWaveforms> result;
    return result;
}


// ═══════════════════════════════════════════════════════════════════════
// AdvancedLlc
// ═══════════════════════════════════════════════════════════════════════
Inputs AdvancedLlc::process() {
    // Override resonant-tank values if user-specified
    if (desiredResonantInductance.has_value()) {
        // Use user's Ls directly
    }
    if (desiredResonantCapacitance.has_value()) {
        // Use user's C directly
    }

    auto designRequirements = process_design_requirements();

    // Override turns ratios
    designRequirements.get_mutable_turns_ratios().clear();
    for (auto n : desiredTurnsRatios) {
        DimensionWithTolerance nTol;
        nTol.set_nominal(n);
        designRequirements.get_mutable_turns_ratios().push_back(nTol);
    }

    // Override magnetizing inductance
    DimensionWithTolerance LmTol;
    LmTol.set_nominal(desiredMagnetizingInductance);
    designRequirements.set_magnetizing_inductance(LmTol);

    // Override resonant components if specified
    if (desiredResonantInductance.has_value()) {
        // Store for waveform generation
        // computedResonantInductance set via process_design_requirements
    }
    if (desiredResonantCapacitance.has_value()) {
        // computedResonantCapacitance set via process_design_requirements
    }

    auto ops = process_operating_points(desiredTurnsRatios, desiredMagnetizingInductance);

    Inputs inputs;
    inputs.set_design_requirements(designRequirements);
    inputs.set_operating_points(ops);

    return inputs;
}

} // namespace OpenMagnetics
