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
    // Turn losses of a PHASOR field (MAS excitation convention, 2026-09-24): inPhaseData and
    // quadratureData are the in-phase and quadrature (Hx, Hy) at the same points, in the same
    // order (see WindingWindowMagneticStrengthFieldPhasorOutput); the loss must be that of
    // |H|^2 = Hx_i^2 + Hy_i^2 + Hx_q^2 + Hy_q^2. For a model whose loss is a quadratic form of
    // the field (every model that sums |H|^2 or squared component averages with fixed
    // coefficients) that is exactly loss(inPhase) + loss(quadrature), which is this default.
    // A model that bridges two DIFFERENT quadratic forms (a low/high-frequency harmonic mean)
    // is not additive and must override it and bridge the summed forms.
    virtual double calculate_turn_losses_from_phasors(Wire wire, double frequency, std::vector<ComplexFieldPoint> inPhaseData, std::vector<ComplexFieldPoint> quadratureData, double temperature);
    // Whether the model consumes the "widthsample" mesh points (width-resolved
    // perpendicular-field samples for flat conductors). Models that average over
    // the lumped surface points must not see them, so the dispatcher strips them
    // unless this returns true.
    virtual bool consumes_width_samples() const { return false; }
    virtual ~WindingProximityEffectLossesModel() = default;
    std::optional<double> try_get_proximity_factor(Wire wire, double frequency, double temperature);
    void set_proximity_factor(Wire wire, double frequency, double temperature, double proximityFactor);
    static std::shared_ptr<WindingProximityEffectLossesModel> factory(WindingProximityEffectLossesModels modelName);
};

class WindingProximityEffectLosses {
  private:
  protected:
  public:
    static std::shared_ptr<WindingProximityEffectLossesModel> get_model(WireType wireType, std::optional<WindingProximityEffectLossesModels> modelOverride = std::nullopt);
    // Phasor field (MagneticField::calculate_magnetic_field_strength_field): the loss of every
    // harmonic is that of the in-phase plus the quadrature field.
    static WindingLossesOutput calculate_proximity_effect_losses(Coil coil, double temperature, WindingLossesOutput windingLossesOutput, WindingWindowMagneticStrengthFieldPhasorOutput windingWindowMagneticStrengthFieldOutput, std::optional<WindingProximityEffectLossesModels> modelOverride = std::nullopt);
    // A plain (real, single-phase) field: there is no quadrature component.
    static WindingLossesOutput calculate_proximity_effect_losses(Coil coil, double temperature, WindingLossesOutput windingLossesOutput, WindingWindowMagneticStrengthFieldOutput windingWindowMagneticStrengthFieldOutput, std::optional<WindingProximityEffectLossesModels> modelOverride = std::nullopt);
    // quadratureFields, when given, must match fields harmonic by harmonic and point by point.
    static std::pair<double, std::vector<std::pair<double, double>>> calculate_proximity_effect_losses_per_meter(Wire wire, double temperature, std::vector<ComplexField> fields, std::optional<WindingProximityEffectLossesModels> modelOverride = std::nullopt, std::optional<std::vector<ComplexField>> quadratureFields = std::nullopt);
  private:
    static WindingLossesOutput calculate_proximity_effect_losses_impl(Coil coil, double temperature, WindingLossesOutput windingLossesOutput, const std::vector<ComplexField>& fieldPerFrequency, const std::optional<std::vector<ComplexField>>& quadratureFieldPerFrequency, std::optional<WindingProximityEffectLossesModels> modelOverride);

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
    // Its width-sample term bridges two different quadratic forms, so it is not additive.
    double calculate_turn_losses_from_phasors(Wire wire, double frequency, std::vector<ComplexFieldPoint> inPhaseData, std::vector<ComplexFieldPoint> quadratureData, double temperature) override;
    bool consumes_width_samples() const override { return true; }
  private:
    double calculate_turn_losses_impl(Wire& wire, double frequency, const std::vector<ComplexFieldPoint>& data, const std::vector<ComplexFieldPoint>* quadratureData, double temperature);
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

// Based on "Eddy currents in rectangular conductors: Analytical 2D loss model in the context of
// magnetic component design", Thomas Ewald and Jürgen Biela, ETH Zurich, EPE'23 ECCE Europe.
// https://doi.org/10.23919/epe23ecceeurope58414.2023.10264392
//
// Solves the Helmholtz equation INSIDE the rectangular conductor with a separable formulation
// whose four tangential surface fields H1..H4 are matched by hyperbolic profiles, takes the
// electric field from Ampere's law and the loss from Poynting's theorem. The result is the
// closed-form complex power of the paper's equation (7); the external (proximity) half of it is
// what this model returns:
//
//     P_ext = rho * Re[ psi_x * H_ext,x^2 + psi_y * H_ext,y^2 ]
//     psi_x = gamma * w * tanh(gamma * h / 2)      psi_y = gamma * h * tanh(gamma * w / 2)
//     gamma = (1 + i) / delta
//
// ITS STATED SIMPLIFICATION, which is the reason MKF carries it as a comparison rather than as
// the answer: section II-C assumes "the field on the conductor's surface is spatially
// homogeneous". That removes edge crowding by construction, and the paper's own Table I records
// the consequence — the model reads -12 % at 1 MHz on a square conductor and -19 % / -23 % at
// 100 kHz / 1 MHz on a mixed-field case, with the authors noting it "underestimates the losses,
// which might be problematic in terms of optimizations".
class WindingProximityEffectLossesEwaldModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Ewald";
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
};

// MARTINEZ MODEL (ABT #1188) - proximity loss with FIELD EXCLUSION.
//
// Every other proximity model in this file computes the loss a conductor dissipates in a field
// that passes through it unchanged. A real conductor thicker than about a skin depth does not
// let the field through: the eddy currents it induces produce their own field, which cancels the
// incident one inside the metal and crowds the remainder onto the surface, concentrated at the
// edges. That exclusion is the whole of this model.
//
// RECTANGULAR-SECTION WIRES ONLY. The functional below would evaluate for any cross-section, but
// round and litz wire keep Ferreira, which is the EXACT Bessel solution for an isolated cylinder
// and therefore already carries this exclusion without bridging between two asymptotes.
// calculate_turn_losses throws for them rather than silently give the worse answer.
//
// The geometry enters through the conductor's own cross-section, not a per-wire-type formula:
//
//   HIGH FREQUENCY (delta << thickness): the conductor excludes the field, so the exterior
//   problem is a perfectly conducting cross-section in a transverse field. Taking that section
//   as the ellipse of the conductor's semi-axes (the strip's edge-singularity-free equivalent -
//   a rectangle's surface current diverges at the corners, an ellipse's does not), the surface
//   field is H_s(v) = H0 (a + b) |cos v| / sqrt(a^2 sin^2 v + b^2 cos^2 v): zero at the middle of
//   the wide faces, peaking at the edges with the classical (1 + a/b) enhancement that the
//   standard demagnetising factor N = a/(a + b) also gives. The loss is then the surface-
//   resistance integral P = 1/2 (rho/delta) * CONTOUR-INTEGRAL H_s^2 dl, which is exact for a
//   circle (it returns 2*pi, the known perfectly-conducting-cylinder result) and grows like
//   ln(wide/thin) for a strip.
//
//   LOW FREQUENCY (delta >> thickness): the field penetrates fully, the induced currents are
//   resistance-limited, and the loss is the vector-potential integral over the section, growing
//   as omega^2.
//
//   The two are joined by the same harmonic-mean bridge the fringing kernel uses, applied to
//   EACH field component separately and then summed - NOT to the field magnitude, which
//   over-predicts (158 % on a representative winding) because the two components have different
//   exclusion factors, different characteristic dimensions and so different corner frequencies.
//
// The driving field is the TOTAL local field - the gap-fringing field (Roshen / Albach) and the
// transport-current field of every other turn, cross term included - sampled across the
// conductor's wide face where those samples exist.
//
// VALIDATION (ABT #1188): against 2D FEM on the stadium cross-section that real rectangular wire
// actually has, at 2:1 aspect and both field orientations, this model reads 0.94-1.08 of FEM
// (median 0.99) at 10 kHz, 100 kHz and 1 MHz. The Ewald model above reads 0.17-1.15 (median
// 0.44) over the same six points.
//
// WHAT IT DOES NOT YET DO, stated so nobody assumes otherwise: the MUTUAL reaction. Each
// conductor excludes the field from its own volume, but the field its induced currents throw
// back at its neighbours is not fed round a second time. That is a coupled solve over the
// induced dipole moments, and it is the natural next step; without it, closely packed turns are
// still charged more field than they really see.
class WindingProximityEffectLossesMartinezModel : public WindingProximityEffectLossesModel {
  public:
    std::string methodName = "Martinez";
    // The exclusion factor C for a field along one of the conductor's axes. A conductor excludes
    // BOTH field components, not only the one normal to its wide face: the same ellipse solution
    // serves each, with the roles of the two axes swapped.
    //   fieldAlongWide = false -> field normal to the wide face. The surface current runs round the
    //     section and crowds at the edges with the classical (1 + a/b) peak. 5.67 for a 2:1 wire.
    //   fieldAlongWide = true  -> field along the wide face. The conductor still excludes it, and
    //     the peak is (1 + b/a). 4.03 for the same wire -- NOT the 2.0 that "two flat faces, no
    //     crowding" would give, which is where this model was a factor 2.016 low until 2D FEM
    //     measured it (ABT #1188).
    static double calculate_exclusion_factor(double wideDimension, double thinDimension, bool fieldAlongWide);
    // Named alias for the perpendicular case, which is what "edge crowding" means at the call site.
    static double calculate_edge_crowding_factor(double wideDimension, double thinDimension);
    double calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature);
    // Its perpendicular term bridges two different quadratic forms, so it is not additive.
    double calculate_turn_losses_from_phasors(Wire wire, double frequency, std::vector<ComplexFieldPoint> inPhaseData, std::vector<ComplexFieldPoint> quadratureData, double temperature) override;
    bool consumes_width_samples() const override { return true; }
  private:
    double calculate_turn_losses_impl(Wire& wire, double frequency, const std::vector<ComplexFieldPoint>& data, const std::vector<ComplexFieldPoint>* quadratureData, double temperature);
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