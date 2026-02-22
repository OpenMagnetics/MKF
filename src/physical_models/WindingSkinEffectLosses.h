#pragma once
#include "Constants.h"
#include "MAS.hpp"
#include "physical_models/Resistivity.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"
#include "Models.h"
#include "Defaults.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

using namespace MAS; // QUAL-001 TODO: qualify types and remove

namespace OpenMagnetics {


class WindingSkinEffectLossesModel {
 private:
 protected:
    std::map<size_t, std::map<double, std::map<double, double>>> _skinFactorPerWirePerFrequencyPerTemperature;
    std::optional<double> try_get_skin_factor(Wire wire, double frequency, double temperature);
    void set_skin_factor(Wire wire, double frequency, double temperature, double skinFactor);

 public:
    std::string methodName = "Default";
    virtual double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0) = 0;
    static std::shared_ptr<WindingSkinEffectLossesModel> factory(WindingSkinEffectLossesModels modelName);
};


class WindingSkinEffectLosses {
 private:
 protected:
 public:
    static std::shared_ptr<WindingSkinEffectLossesModel> get_model(WireType wireType, std::optional<WindingSkinEffectLossesModels> modelOverride = std::nullopt);
    static double calculate_skin_depth(WireMaterialDataOrNameUnion material, double frequency, double temperature);
    static double calculate_skin_depth(Wire wire, double frequency, double temperature);
    static WindingLossesOutput calculate_skin_effect_losses(Coil coil, double temperature, WindingLossesOutput windingLossesOutput, double windingLossesHarmonicAmplitudeThreshold = Defaults().harmonicAmplitudeThreshold, std::optional<WindingSkinEffectLossesModels> modelOverride = std::nullopt);
    static std::pair<double, std::vector<std::pair<double, double>>> calculate_skin_effect_losses_per_meter(Wire wire, SignalDescriptor current, double temperature, double currentDivider = 1, std::optional<WindingSkinEffectLossesModels> modelOverride = std::nullopt);
};


// ============================================================================
// DOWELL MODEL (1966)
// "Effects of eddy currents in transformer windings"
// P.L. Dowell, Proc. IEE, Vol. 113, No. 8, pp. 1387-1394, August 1966
//
// SKIN EFFECT ONLY (m=1 case of the general Dowell equation).
//
// Full Dowell equation for m layers (Eq. 10):
//   FR = M' + (m² − 1)/3 · D'
//   M' = Re[α·h·coth(α·h)]   ← SKIN term (implemented here)
//   D' = Re[2·α·h·tanh(α·h/2)]  ← PROXIMITY term (in proximity model)
//   α  = (1+j)/δ
//
// For m=1 (skin only):
//   FR_skin = M' = ζ·(sinh 2ζ + sin 2ζ) / (cosh 2ζ − cos 2ζ)
//   ζ = h/δ  (conductor height / skin depth)
//
// Round wire: equivalent square h_eq = d·√(π/4)  [same cross-sectional area]
// Limits: ζ→0: FR→1 (DC);  ζ→∞: FR→ζ (thin skin layer)
// ============================================================================
class WindingSkinEffectLossesDowellModel : public WindingSkinEffectLossesModel {
public:
    std::string methodName = "Dowell";
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// WOJDA MODEL (2018)
// "Winding Resistance and Power Loss of Inductors With Litz and Solid-Round
//  Wires" by R.P. Wojda, M.K. Kazimierczuk
// IEEE Trans. Industry Applications, Vol. 54, No. 5, 2018
//
// SKIN EFFECT ONLY.
// Penetration ratio (Eq. 7):
//   A = (π/4)^(3/4) · d/δ · √(d/d_outer)
// Skin factor:
//   Fs = (A/2)·(sinh A + sin A) / (cosh A − cos A)
// Litz wire: uses strand diameter and strand outer diameter.
// ============================================================================
class WindingSkinEffectLossesWojdaModel : public WindingSkinEffectLossesModel {
public:
    std::string methodName = "Wojda";
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_penetration_ratio(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// ALBACH MODEL (2017)
// "Induktivitäten in der Leistungselektronik" by M. Albach
// Springer Vieweg, 2017, Chapter 4
//
// SKIN EFFECT + internal bundle proximity for litz wire.
// Modified Bessel form (Eq. 4.7):
//   Ks = (ξ/2) · Re[I₀(ξ) / I₁(ξ)]    ξ = (1+j)·r/δ
//
// For bundled conductors:
//   Ks += n(n−1)·(rD/rO)² · Re[I₁(ξ)/I₀(ξ)]
//
// SEPARATION: This model includes BOTH self-induced skin effect AND
// internal bundle proximity (within litz strand bundle).
// External proximity is handled by the proximity model.
// ============================================================================
class WindingSkinEffectLossesAlbachModel : public WindingSkinEffectLossesModel {
public:
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// PAYNE MODEL (2021)
// "The AC Resistance of Rectangular Conductors" by A. Payne
// Application Note AP101, Issue 4, 2021
//
// SKIN EFFECT ONLY. Empirical, rectangular conductors only.
// No proximity effect modeled.
//
//   p  = √A / (1.26·δ·1000)      [A in mm²]
//   Ff = 1 − exp(−0.026·p)
//   Kc = 1 + Ff·(1.2/exp(2.1·b/a) + 1.2/exp(2.1·a/b))
//   x  = [2δ/b·(1+b/a) + 8(δ/b)³/(a/b)] / [(a/b)^0.33·exp(−3.5b/δ) + 1]
//   Rac/Rdc = Kc / (1 − exp(−x))
//   a = thin dim, b = thick dim
// ============================================================================
class WindingSkinEffectLossesPayneModel : public WindingSkinEffectLossesModel {
public:
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// NAN MODEL (2003)
// "An Improved Calculation of Proximity-Effect Loss in High-Frequency
//  Windings of Round Conductors"
// Xi Nan, C.R. Sullivan, IEEE PESC 2003, pp. 853–860
//
// SKIN EFFECT ONLY (paper is primarily about an improved proximity model).
// Uses Dowell's equivalent-foil skin expression (Eq. 2):
//   FR = η·(sinh 2η + sin 2η) / (cosh 2η − cos 2η)
//   η = d_eq / (2δ)    [half-thickness convention]
//   d_eq = d·√(π/4)   [Dowell equivalent square for round wire]
//
// SEPARATION: Nan's improved proximity (Eq. 7) with porosity correction:
//   FR_prox = η·(sinh 2η − sin 2η)/(cosh 2η + cos 2η) · [2(2m²+1)/3]
// is handled by WindingProximityEffectLossesNanModel.
// ============================================================================
class WindingSkinEffectLossesNanModel : public WindingSkinEffectLossesModel {
public:
    std::string methodName = "Nan";
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// KAZIMIERCZUK MODEL (2019)
// "An Analytical Correction to Dowell's Equation for Inductor and
//  Transformer Winding Losses Using Cylindrical Coordinates"
// D.J. Whitman, M.K. Kazimierczuk
// IEEE Trans. Power Electron., Vol. 35, No. 2, Feb 2020
// DOI: 10.1109/TPEL.2019.2904582
//
// SKIN EFFECT ONLY.
// Corrected Dowell (Eq. 25):
//   FR_cyl = FR_Dowell · C_F(R_in, R_out, h, δ)
//   C_F ≈ 1 + h²/(12·R_in·R_out)   [first-order curvature correction]
//
// Since Wire object does not carry winding radii, we implement the base
// Dowell skin factor. Apply C_F at winding level when radii are available.
//
// SEPARATION: Proximity term (m²−1)/3·D' also needs C_F correction,
// handled by WindingProximityEffectLossesKazimierczukModel.
// ============================================================================
class WindingSkinEffectLossesKazimierczukModel : public WindingSkinEffectLossesModel {
public:
    std::string methodName = "Kazimierczuk";
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// KUTKUT MODEL (1998)
// "A Simple Technique to Evaluate Winding Losses Including Two-Dimensional
//  Edge Effects" by N.H. Kutkut
// IEEE Trans. Power Electronics, Vol. 13, No. 5, September 1998
//
// SKIN EFFECT with 2D edge-effect correction for rectangular conductors.
// Combines low/high frequency asymptotes (Eq. 30):
//   Fr = [1 + (f/fl)^α + (f/fh)^β]^(1/γ)
//   fl = 3.22ρ/(8μ₀·b'·a')       (Eq. 31)
//   fh = π²ρ/(4μ₀·a'²)·K(c/b)⁻²  (Eq. 32, K = elliptic integral)
//   α=2, β=5.5, γ=11
//
// SEPARATION: Isolated conductor only. No proximity effect.
// ============================================================================
class WindingSkinEffectLossesKutkutModel : public WindingSkinEffectLossesModel {
public:
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// FERREIRA MODEL (1994)
// "Improved Analytical Modeling of Conductive Losses in Magnetic Components"
// J.A. Ferreira, IEEE Trans. Power Electron., Vol. 9, No. 1, January 1994
//
// SKIN EFFECT ONLY.
// Ferreira proves orthogonality of skin and proximity — they can be
// calculated independently and summed.
//
// Skin factor for round conductors (Eq. 8):
//   Fs = (ξ/4)·(sinh ξ + sin ξ) / (cosh ξ − cos ξ)
//   ξ = d/δ
//
// SEPARATION: Full FR (Eq. 6) = Fs + (4(m²−1)/3)·Fp
// Only Fs is here. Fp (Kelvin functions) is in the proximity model.
// ============================================================================
class WindingSkinEffectLossesFerreiraModel : public WindingSkinEffectLossesModel {
public:
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// DIMITRAKAKIS MODEL (2007)
// "A New Model for the Determination of Copper Losses in Transformer
//  Windings with Arbitrary Conductor Distribution"
// G.S. Dimitrakakis, E.C. Tatakis, E.J. Rikos, EPE 2007
//
// SKIN EFFECT ONLY (paper covers both skin + proximity).
// Full FR (Eq. 5): FR = (q/2)·r_skin + [4(m²−1)/3+1]·(q/2)·r_prox
//
// Skin portion only (Eq. 6):
//   FR_skin = (q/2)·r_skin
//   r_skin  = [ber(q)·bei'(q) − bei(q)·ber'(q)] / [ber'(q)²+bei'(q)²]
//   q = √2·r/δ
//
// Equivalent Bessel form (used in implementation):
//   FR_skin = Re[α·I₀(α)/I₁(α)] / 2,   α = (1+j)·r/δ
//
// SEPARATION: Proximity r_prox (ber₂/bei₂) is in the proximity model.
// Rectangular/foil fallback: Dowell 1D skin factor.
// ============================================================================
class WindingSkinEffectLossesDimitrakakisModel : public WindingSkinEffectLossesModel {
public:
    std::string methodName = "Dimitrakakis";
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// MUEHLETHALER MODEL (2012)
// "Modeling And Multi-objective Optimization Of Inductive Power Components"
// Jonas Mühlethaler, ETH Zurich, DISS. ETH NO. 20217
//
// SKIN EFFECT ONLY.
//
// Round conductors (Section 4.3.1, Eq. 4.6-4.7, Appendix A.8.1):
//   Bessel ODE: d²Jz/dr² + (1/r)·dJz/dr = jωμσ·Jz   (A.70)
//   Solution:   Jz = C·J₀(j^(3/2)·r/δ)               (A.71)
//   Kelvin form: FR = (ξ/(4√2))·2·[ber₀(ξ)bei₁(ξ)−bei₀(ξ)ber₁(ξ)]
//                                  / [ber₁²(ξ)+bei₁²(ξ)]
//   ξ = d/(√2·δ)
//   Bessel form: FR = Re[(1+j)·(r/δ)·I₀((1+j)r/δ)/I₁((1+j)r/δ)] / 2
//
// Foil conductors (Section 4.4.1, Eq. 4.20, Appendix A.7.1):
//   FF = (Δ/4)·(sinh Δ + sin Δ) / (cosh Δ − cos Δ),   Δ = h/δ
//
// SEPARATION:
//   Round prox: PP = RDC·GR(f)·Ĥe²  (GR uses ber₂/bei₂, Eq. 4.8-4.9)
//   Foil prox:  PP = RDC·GF(f)·Ĥe²  (GF = b²·Δ·(sinhΔ−sinΔ)/(coshΔ+cosΔ))
// Both handled by WindingProximityEffectLossesMuehlethalerModel.
// ============================================================================
class WindingSkinEffectLossesMuehlethalerModel : public WindingSkinEffectLossesModel {
public:
    std::string methodName = "Muehlethaler";
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// WANG MODEL (2018)
// "Improved Analytical Calculation of High Frequency Winding Losses in
//  Planar Inductors"
// X. Wang, L. Wang, L. Mao, Y. Zhang, IEEE APEC 2018
//
// SKIN EFFECT with 2D edge-effect correction for rectangular conductors.
// 2D model (Eq. 2):
//   P_2D = P_x(Hx) + P_y(Hy)  [orthogonal field components]
// 1D skin function: f(z) = z·(sinh 2z+sin 2z)/(cosh 2z−cos 2z)
//
// Boundary values for isolated conductor (Eq. 3-5):
//   Hx1 = I/(2c)
//   Hy1: linear current model,  λ = 0.01·(c/h) + 0.66
//   Combined: FR = f(h/δ) + (Hy/Hx)²·(h/c)·f(c/δ)
//
// Accuracy: 2–5% error vs FEM (vs 45–50% for 1D Dowell).
// SEPARATION: Multi-layer proximity (Section IV, Eq. 6-7) → proximity model.
// Round wire: fallback to Dowell skin factor.
// ============================================================================
class WindingSkinEffectLossesWangModel : public WindingSkinEffectLossesModel {
public:
    std::string methodName = "Wang";
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// HOLGUIN MODEL (2014)
// "Power Losses Calculations in Windings of Gapped Magnetic Components"
// F.A. Holguin, R. Prieto, R. Asensi, J.A. Cobos, IEEE APEC 2014
//
// SKIN EFFECT ONLY. Uses exact Bessel solution (from Lammeraner [13]):
//   FR = Re[(1+j)·(r₀/δ)·I₀((1+j)r₀/δ) / I₁((1+j)r₀/δ)] / 2
//
// Asymptotic forms (Eq. 9):
//   r₀/δ → 0: FR ≈ 1 + (1/48)·(r₀/δ)⁴
//   r₀/δ → ∞: FR ≈ (r₀/δ)/(2√2) + 3/(32·(r₀/δ))
//
// SEPARATION: Gap-fringing proximity from Eq. (1)-(5):
//   P_prox = (4πr₀²H₀²/σ)·[ber·bei'−ber'·bei]/[ber²+bei²]
// is handled by WindingProximityEffectLossesHolguinModel.
// Rectangular/foil fallback: Dowell 1D skin factor.
// ============================================================================
class WindingSkinEffectLossesHolguinModel : public WindingSkinEffectLossesModel {
public:
    std::string methodName = "Holguin";
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// LOTFI MODEL
// Uses complete elliptic integrals for rectangular/elliptical cross-sections.
// SKIN EFFECT ONLY for isolated conductor.
// ============================================================================
class WindingSkinEffectLossesLotfiModel : public WindingSkinEffectLossesModel {
public:
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


// ============================================================================
// PERRY MODEL (1979)
// "Multiple Layer Series Connected Winding Design For Minimum Losses"
// M.P. Perry, IEEE Trans. PAS, Vol. PAS-98, No. 1, Jan/Feb 1979
//
// SKIN EFFECT ONLY (single-layer case of a multi-layer paper).
// From Eq. (1) and Eq. (13a):
//   FR = Δ·(sinh 2Δ + sin 2Δ) / (cosh 2Δ − cos 2Δ),   Δ = T/δ
//
// Mathematically identical to Dowell's M'.
// Perry's contributions: cylindrical coordinate solution (Section III)
// and per-layer optimized thickness design (Eq. 18-19).
//
// SEPARATION: Multi-layer proximity (Eq. 12):
//   Q_n uses F₁'(Δ) [=skin term] and F₂'(Δ) [=proximity coupling term]
// Only F₁' is implemented here.
// ============================================================================
class WindingSkinEffectLossesPerryModel : public WindingSkinEffectLossesModel {
public:
    std::string methodName = "Perry";
    double calculate_skin_factor(const Wire& wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, double currentRms = 0);
};


} // namespace OpenMagnetics
