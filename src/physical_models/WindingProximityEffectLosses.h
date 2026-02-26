#pragma once
#include "Constants.h"
#include "MAS.hpp"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"
#include "Models.h"

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


class WindingProximityEffectLossesModel {
  private:
  protected:
    std::map<size_t, std::map<double, std::map<double, double>>> _proximityFactorPerWirePerFrequencyPerTemperature;
  public:
    std::string methodName = "Default";
    virtual double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) = 0;
    std::optional<double> try_get_proximity_factor(Wire wire, double frequency, double temperature);
    void set_proximity_factor(Wire wire, double frequency, double temperature, double proximityFactor);
    static std::shared_ptr<WindingProximityEffectLossesModel> factory(WindingProximityEffectLossesModels modelName);
};

class WindingProximityEffectLosses {
  private:
  protected:
  public:
    static std::shared_ptr<WindingProximityEffectLossesModel> get_model(WireType wireType, std::optional<WindingProximityEffectLossesModels> modelOverride = std::nullopt);
    static WindingLossesOutput calculate_proximity_effect_losses(Coil coil, double temperature, WindingLossesOutput windingLossesOutput, WindingWindowMagneticStrengthFieldOutput windingWindowMagneticStrengthFieldOutput, std::optional<WindingProximityEffectLossesModels> modelOverride = std::nullopt);
    static std::pair<double, std::vector<std::pair<double, double>>> calculate_proximity_effect_losses_per_meter(Wire wire, double temperature, std::vector<ComplexField> fields, std::optional<WindingProximityEffectLossesModels> modelOverride = std::nullopt);

};

// Based on Measurement and Characterization of High Frequency Losses in Nonideal Litz Wires by Hans Rossmanith
// https://sci-hub.st/10.1109/tpel.2011.2143729
class WindingProximityEffectLossesRossmanithModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Rossmanith";
    static double calculate_proximity_factor(Wire wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);

};

// Based on Improved Analytical Calculation of High Frequency Winding Losses in Planar Inductors by Xiaohui Wang
// https://sci-hub.st/10.1109/ECCE.2018.8558397
class WindingProximityEffectLossesWangModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Wang";
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// Based on A New Approach to Analyse Conduction Losses in High Frequency Magnetic Components by J.A. Ferreira
// https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9485268
class WindingProximityEffectLossesFerreiraModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Ferreira";
    static double calculate_proximity_factor(Wire wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};


// Based on Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften by Manfred Albach
// https://libgen.rocks/get.php?md5=94b7f2906f53602f19892d7f1dabd929&key=YMKCEJOWB653PYLL
class WindingProximityEffectLossesAlbachModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Albach";
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// Based on Eddy currents by Jiří Lammeraner
// https://archive.org/details/eddycurrents0000lamm
class WindingProximityEffectLossesLammeranerModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Lammeraner";
    static double calculate_proximity_factor(Wire wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// Based on Effects of eddy currents in transformer windings by P. L. Dowell
// https://ieeexplore.ieee.org/document/5247417
class WindingProximityEffectLossesDowellModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Dowell";
    static double calculate_proximity_factor(Wire wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// ============================================================================
// XI NAN MODEL (2003)
// "An Improved Calculation of Proximity-Effect Loss in High-Frequency
//  Windings of Round Conductors" by Xi Nan and Charles R. Sullivan
// IEEE PESC 2003, pp. 853-860
//
// Modified Dowell function fitted to 2-D FEM results for round conductors.
// Proximity factor G is unitless: P/l = G * H² / σ  (Eq. 1)
// Converted to our architecture: factor [Ω·m] = G · ρ
//
// Core functional form (Eq. 15/17/18):
//   G' = k1·√k2·X · [sinh(√k2·X) - sin(√k2·X)] / [cosh(√k2·X) + cos(√k2·X)]
//   g(X) = K·X / (X^(-3n) + b^(3n))^(1/3n)
//   G = (1-w)·G' + w·g(X)
// Default coefficients for square packing (v/d≈1, h/d≈1), Table I:
//   k1=1.45, k2=0.33, K=0.29, b=1.1, n=3, w=0.54
//
// Valid for: ROUND, LITZ. Accuracy <2% vs FEM.
// ============================================================================
class WindingProximityEffectLossesNanModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Nan";
    static double calculate_proximity_factor(Wire wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// ============================================================================
// WOJDA MODEL (2012)
// "Proximity-effect winding loss in different conductors using magnetic
//  field averaging" by R.P. Wojda and M.K. Kazimierczuk
// COMPEL Vol. 31 No. 6, 2012, pp. 1793-1814. DOI: 10.1108/03321641211267128
//
// Low/medium-frequency closed-form proximity model for multiple wire shapes.
// Derived via field averaging (Ampère's law), not Bessel functions.
//
//   Foil/Strip (Eq. 42): R_pe/l = ηh²·μ₀²·ω²·h³ / (12·ρ·b)
//   Square     (Eq. 58): R_pe   = ηs²·μ₀²·ω²·Nl²·lT·h² / (36·ρ)
//   Round      (Eq. 70): R_pe   = ηb²·π²·μ₀²·ω²·Nl²·lT·d² / (576·ρ)
//
// Proximity factor per single conductor per meter in external field He:
//   factor = μ₀²·ω²·A_geom / (K_shape·ρ)
//
// Valid for: FOIL, PLANAR, RECTANGULAR, ROUND, LITZ. Low-medium freq (d/δ < 2).
// NOTE: Skin effect is EXCLUDED (orthogonality principle). Pure proximity only.
// NOTE: Gap fringing NOT included. H is taken as external input.
// ============================================================================
class WindingProximityEffectLossesWojdaModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Wojda";
    static double calculate_proximity_factor(Wire wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// ============================================================================
// SULLIVAN SFD MODEL (2001)
// "Computationally Efficient Winding Loss Calculation with Multiple Windings,
//  Arbitrary Waveforms, and Two-Dimensional or Three-Dimensional Field Geometry"
// by Charles R. Sullivan
// IEEE Trans. Power Electronics, Vol. 16, No. 1, January 2001. DOI: 10.1109/63.903999
//
// Squared-Field-Derivative (SFD) method. Low-frequency approximation (d << δ).
//
// From Appendix A, Eq. (1):
//   P_inst/l = (π·d⁴)/(128·ρ) · |dB/dt|²
// For sinusoidal B = μ₀·H·sin(ωt):
//   factor = π·μ₀²·ω²·d⁴ / (128·ρ)  [Ω·m per (A/m)²]
//
// Valid for: ROUND, LITZ. Accuracy good only when d << δ.
// NOTE: Gap fringing effects belong in the MagneticField class, NOT here.
//       The H field is taken as external input from the MagneticField module.
// ============================================================================
class WindingProximityEffectLossesSullivanModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Sullivan";
    static double calculate_proximity_factor(Wire wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// ============================================================================
// BARTOLI MODEL (1996)
// "Modeling Litz-Wire Winding Losses in High-Frequency Power Inductors"
// by M. Bartoli, N. Noferi, A. Reatti, M.K. Kazimierczuk
// IEEE PESC 1996, pp. 1690-1696. DOI: 10.1109/PESC.1996.548808
//
// Separates EXTERNAL proximity (bundle in external field) and INTERNAL proximity
// (strand-to-strand within bundle) using Kelvin-Bessel functions.
//
// Bessel proximity factor (Eq. 10/13/15/16):
//   K2(y) = (y/2) · |ber₂y·ber'y - bei₂y·bei'y| / (ber²y + bei²y)
//   y_s = d_s · √2 / δ
//   External per strand: factor_ext = 2π·ρ·K2(y_s)
//   Internal per strand: factor_int = 2π·ρ·K2(y_s)·k_s·0.5  (×0.5 for twisted litz)
//   Total: factor = factor_ext + factor_int
//
// NOTE: Skin effect term in Eq. 19 is IDENTICAL to Ferreira's Bessel skin model.
//       No new skin model needed; use WindingSkinEffectLossesFerreiraModel.
// NOTE: Gap fringing NOT included.
//
// Valid for: LITZ (primary use), ROUND (ns=1 reduces to Ferreira proximity).
// ============================================================================
class WindingProximityEffectLossesBartoliModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Bartoli";
    static double calculate_proximity_factor(Wire wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// ============================================================================
// VANDELAC MODEL (1988)
// "A Novel Approach for Minimizing High-Frequency Transformer Copper Losses"
// by J.-P. Vandelac and P.D. Ziogas
// IEEE Trans. Power Electronics, Vol. 3, No. 3, July 1988. DOI: 10.1109/63.17953
//
// Layer-level eddy current model using real-valued F1/F2 functions (Eq. 25/29/30).
// Implemented here as PROXIMITY ONLY (α=1, β=0 — uniform external field, no self-current).
//
//   Q_prox = H²_e / (2·σ·δ) · [3·F1(p) - 4·F2(p)]  [W/m²]
//   F1(p) = [sinh(2p) + sin(2p)] / [cosh(2p) - cos(2p)]   (Eq. 29)
//   F2(p) = [sinh(p)·cos(p) + cosh(p)·sin(p)] / [cosh(2p) - cos(2p)]  (Eq. 30)
//   p = h / (δ·√2)
//   factor = a · ρ · [3·F1(p) - 4·F2(p)] / (2·δ)  [Ω·m per (A/m)²]
//
// Equivalence to Dowell: F1 = Re[αh·coth(αh)], F2 = Re[αh·csch(αh)]
// Skin-only case (α=0): Q_skin = H²·F1(p)/(2·σ·δ) → identical to Dowell M'
// → NO new skin model needed; Dowell skin already covers this.
//
// NOTE: The FULL Vandelac model with arbitrary α = H_z(0)/H_z(h) (for interleaved
//       flyback/push-pull winding optimization) requires layer-level boundary fields
//       and belongs in a separate winding design/optimization class.
// NOTE: Gap fringing NOT included. H is external input from MagneticField module.
//
// Valid for: FOIL, PLANAR, RECTANGULAR, ROUND (round→square equiv.), LITZ. All freqs.
// ============================================================================
class WindingProximityEffectLossesVandelacModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Vandelac";
    static double calculate_proximity_factor(Wire wire, double frequency, double temperature);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

} // namespace OpenMagnetics