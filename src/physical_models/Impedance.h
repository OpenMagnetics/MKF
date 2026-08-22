#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "physical_models/MagnetizingInductance.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/CoreLosses.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include <MAS.hpp>
#include <complex>
#include <optional>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

// Frequency-independent building blocks of the differential-mode impedance:
// the leakage inductance, the winding DC resistance and the inter-winding
// (through-core) capacitance. Computed once per magnetic so a frequency sweep
// does not recompute the (expensive) stray-capacitance model at every point.
struct DifferentialModeParameters {
    double leakageInductance;
    double windingResistance;
    double interWindingCapacitance;
};

// One arm of the wideband terminal-impedance model: a parallel-RLC tank. The
// terminal impedance is the SERIES sum of these tanks (a Foster ladder), each
// contributing one resonance. usesCorePermeability distinguishes the two kinds:
//   - magnetizing tank: the inductance is the air-cored inductance to be
//     multiplied by the core complex permeability µ(f) at each frequency, and
//     its loss comes from µ''(f) (no series resistor). This is the first/main
//     resonance, present for every magnetic.
//   - leakage tank: a flat air-cored leakage inductance in series with the
//     (frequency-dependent) winding resistance, shunted by the inter-winding
//     capacitance. Each additional winding adds one — these are the higher
//     self-resonances (the common-mode-choke "second resonance", generalized to
//     any multi-winding magnetic).
// capacitance <= 0 means the tank has no shunt branch (a bare series arm).
struct ImpedanceTank {
    double inductance;
    double capacitance;
    bool usesCorePermeability;
    // Leakage tanks only: the leakage loop runs through winding 0 and this
    // secondary winding, so the series damping resistance, referred to winding 0,
    // is R_0(f) + turnsRatioSquared·R_secondary(f) with turnsRatioSquared = (N_0/N_j)².
    size_t secondaryWindingIndex = 0;
    double turnsRatioSquared = 0.0;
};

// Precomputed, frequency-independent building blocks of the full wideband
// terminal impedance. Built once per magnetic (the heavy models — reluctance,
// stray capacitance, leakage, DC resistance — run here, not per frequency) and
// evaluated cheaply at each frequency by impedance_from_model().
struct WidebandImpedanceModel {
    std::optional<CoreMaterial> coreMaterial;   // complex permeability µ(f) for the magnetizing tank
    // ABT #848: the core's cross-section (column width, depth) for the dimensional / eddy-
    // dielectric attenuation of the magnetizing permeability. Only acts when the material
    // carries permittivity data; empty or missing data means no correction, never a guess.
    std::vector<double> coreCrossSectionDimensions;
    double permeabilityScaling = 1.0;           // DC-bias rolloff factor µ(Hdc)/µ(0)
    std::vector<ImpedanceTank> tanks;           // [0] magnetizing; the rest leakage resonances
    // Per-winding frequency-dependent resistance data, indexed by winding, used to
    // damp the leakage tanks (each leakage loop referred to winding 0 as
    // R_0(f) + n²·R_j(f)). Empty for single-winding magnetics.
    std::vector<double> windingResistanceDc;          // DC resistance per winding (Ohm)
    std::vector<Wire> windingWire;                    // per winding, for the analytic skin factor
    std::vector<double> windingDcResistancePerMeter;  // per winding, denominator of the skin factor
    bool hasWindingResistanceModel = false;           // fast path: DC + analytic skin factor
    bool hasProximityModel = false;                   // slow path: field-based effective resistance
    std::optional<Magnetic> magnetic;                 // only needed by the slow proximity path
    bool fast = true;
    double temperature = Defaults().ambientTemperature;
};

class Impedance {
    private:
        bool _fastCapacitance;
        size_t _maximumNumberTurns = 200;
        // The magnetizing tank (air-cored inductance ∥ winding self-capacitance),
        // the first resonance shared by calculate_impedance and the wideband model.
        ImpedanceTank build_magnetizing_tank(Core& core, Coil& coil);
    protected:
    public:
    // The intended flagship is the FULL energy-based StrayCapacitance model (per-pair
    // Albach statics, floating-core self term, mirrored-winding SRF factor — ABT #848).
    // Controlled A/B on 107 measured WE common-mode chokes, same materials and complex-
    // permeability model on both arms (2026-08-21): full model 92% peak / 45% frequency
    // within 2x; fast path 99% / 74%. The whole remaining gap is the cased (nylon)
    // SC-CMC family, whose measured peaks scale with turn count (genuinely LC-set) at a
    // capacitance no honest electrostatic element delivers through a 1.25 mm case — the
    // fast path lands them via its near-contact acosh term, accuracy not yet explained
    // by physics. Until #848 resolves that, the fast path stays the default so shipped
    // results are the best available; flip deliberately, with a rescore.
    Impedance(bool fastCapacitance=true) {
        _fastCapacitance = fastCapacitance;
    } 

    virtual ~Impedance() = default;

        // ABT #848: complex factor mu_eff/mu for the EM wave crossing the core cross-section,
        // tan(k d/2)/(k d/2) per dimension (Snelling, Soft Ferrites, dimensional resonance), with
        // k = w*sqrt(mu0*mu*eps0*eps) from the material's complex permittivity. Returns 1 when the
        // material has no permittivity data.
        static std::complex<double> core_dimensional_attenuation(const CoreMaterial& material, double frequency, std::complex<double> complexPermeability, const std::vector<double>& crossSectionDimensions);

    std::complex<double> calculate_impedance(Magnetic magnetic, double frequency, double temperature = Defaults().ambientTemperature);
    std::complex<double> calculate_impedance(Core core, Coil coil, double frequency, double temperature = Defaults().ambientTemperature);
    std::complex<double> calculate_impedance(Core core, Coil coil, double frequency, double magneticFieldDcBias, double temperature);
    // Differential-mode impedance of a coupled (e.g. common-mode-choke) magnetic:
    // the opposing winding currents cancel the core flux, so the inductance seen
    // is the leakage inductance, which parallel-resonates with the inter-winding
    // capacitance (the second self-resonance).
    std::complex<double> calculate_differential_mode_impedance(Magnetic magnetic, double frequency, double temperature = Defaults().ambientTemperature);
    std::complex<double> calculate_differential_mode_impedance(Core core, Coil coil, double frequency, double temperature = Defaults().ambientTemperature);
    // Compute the frequency-independent DM parameters once (the leakage is taken at
    // referenceFrequency — it is essentially air-cored and flat). Then evaluate the
    // impedance cheaply at each frequency. A sweep should call the first once and the
    // second per point rather than calling calculate_differential_mode_impedance N times.
    DifferentialModeParameters calculate_differential_mode_parameters(Core core, Coil coil, double referenceFrequency, double temperature = Defaults().ambientTemperature);
    std::complex<double> differential_mode_impedance_from_parameters(const DifferentialModeParameters& parameters, double frequency);
    // Wideband terminal impedance as a series cascade of resonant tanks (a Foster
    // ladder): the magnetizing tank (first resonance) plus one leakage tank per
    // additional winding (the higher resonances — the CMC "second resonance",
    // generalized to any multi-winding magnetic; a single-winding inductor keeps
    // its single resonance). Build the model once per magnetic, then evaluate it
    // cheaply at each frequency. A sweep must call the first once and the second
    // per point. referenceFrequency is where the (frequency-flat) leakage is taken.
    WidebandImpedanceModel build_wideband_impedance_model(Magnetic magnetic, double referenceFrequency, double temperature = Defaults().ambientTemperature, bool fast = true);
    // Common-mode impedance of a coupled magnetic (e.g. a common-mode choke): the
    // datasheet/REDEXPERT CM measurement drives all windings in parallel, so the
    // leakage between them carries no net ampere-turns and NO leakage resonance
    // appears at the terminals — the CM impedance is the magnetizing tank alone
    // (core µ(f), damped by µ''(f), in parallel with the winding self-capacitance).
    // Contrast with build_wideband_impedance_model (one winding driven, the others
    // floating), whose leakage tanks model the higher self-resonances of THAT
    // measurement; serving that curve as "CM impedance" adds a near-undamped
    // leakage spike that measured CM curves do not show (ABT #167).
    WidebandImpedanceModel build_common_mode_impedance_model(Magnetic magnetic, double temperature = Defaults().ambientTemperature);
    std::complex<double> impedance_from_model(const WidebandImpedanceModel& model, double frequency);
    double calculate_q_factor(Magnetic magnetic, double frequency, double temperature = Defaults().ambientTemperature);
    double calculate_q_factor(Core core, Coil coil, double frequency, double temperature = Defaults().ambientTemperature);
    double calculate_self_resonant_frequency(Magnetic magnetic, double temperature = Defaults().ambientTemperature);
    double calculate_self_resonant_frequency(Core core, Coil coil, double temperature = Defaults().ambientTemperature);
    int64_t calculate_minimum_number_turns(Magnetic magnetic, Inputs inputs);

};

} // namespace OpenMagnetics