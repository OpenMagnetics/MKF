#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "Models.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Wire.h"

#include "processors/Inputs.h"
#include "support/CoilMesher.h"
#include <MAS.hpp>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

using namespace MAS; // QUAL-001 TODO: qualify types and remove

namespace OpenMagnetics {

// ABT #1645: thin-strip transmission. A foil turn is a strip of width W (its wide conducting dimension) and thickness
// t << W. In a uniform field H0 normal to its wide face it carries, to first order, the eddy sheet current
// K(x) = j omega mu0 sigma t H0 x (Faraday, curl E = -j omega mu0 H, with x across the width from the centre and no
// net current; t << delta, so the current is uniform through the thickness). Line currents K dx at +-x give the
// centre a normal field -K dx / (2 pi x) each, so the sheet adds -(dK/dx) W / (2 pi) = -j omega mu0 sigma t W / (2 pi)
// H0 = -j W t / (pi delta^2) H0 (omega mu0 sigma = 2 / delta^2), opposing H0. Solved self-consistently (the eddy current is driven by the total normal field H, not
// by H0): H = H0 - j (W t / (pi delta^2)) H, so the strip transmits the normal component with
// T = 1 / (1 + j W t / (pi delta^2)) and leaves the tangential one alone. It applies to the field of any source
// reaching a conductor (foil, round or litz) through the strip: the strips counted are the FOIL turns strictly
// between the source and the induced point whose plane the straight source-to-point line crosses within the strip's
// width. T is complex and per harmonic; it multiplies the normal component of the phasor (in-phase and quadrature).
struct ScreeningStrip {
    double centreX;
    double centreY;
    double tangentX;  // along the wide dimension
    double tangentY;
    double normalX;   // normal to the wide face
    double normalY;
    double width;     // the wide conducting dimension, W
    double thickness; // the thin conducting dimension, t
    size_t windingIndex;
    int64_t turnIndex;
};

// A unit-current field contribution (Hx, Hy) that crosses `strips` (indices into the strip list, ordered from the
// source to the induced point) on its way.
struct ScreenedFieldSum {
    size_t windingIndex;
    std::vector<uint32_t> strips;
    double fieldX;
    double fieldY;
};

// The strips strictly between the source (sourceX, sourceY) and the induced point, ordered from the source; the
// source's own strip (excludedSourceTurn, for a real turn) and the induced point's own strip (excludedInducedTurn)
// never count.
std::vector<uint32_t> strips_between(const std::vector<ScreeningStrip>& strips, double sourceX, double sourceY, double pointX, double pointY,
                                     std::optional<int64_t> excludedSourceTurn, std::optional<int64_t> excludedInducedTurn);

class MagneticFieldStrengthModel {
    public:
        std::vector<Wire> _wirePerWinding;
        std::vector<double> _wireMaxOuterWidth;  // Precomputed for performance
        std::vector<double> _wireMaxOuterHeight; // Precomputed for performance
        // ABT #376: winding breadth b — the window dimension PARALLEL to the layers, over which
        // Dowell's one-dimensional field is assumed uniform (H = enclosed MMF / b). Only the
        // Dowell model reads it; the point-to-point models derive everything from the two points.
        double _windingWindowBreadth = 0;
        virtual ComplexFieldPoint get_magnetic_field_strength_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt) = 0;
        // (Hx, Hy) only: the point-pair sum calls this for every inducing/induced pair, and building a
        // ComplexFieldPoint (a heap-allocated point copy) per pair dominated its cost.
        virtual std::pair<double, double> get_magnetic_field_strength_components_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt) {
            auto complexFieldPoint = get_magnetic_field_strength_between_two_points(inducingFieldPoint, inducedFieldPoint, inducingWireIndex);
            return {complexFieldPoint.get_real(), complexFieldPoint.get_imaginary()};
        }
};

class MagneticFieldStrengthFringingEffectModel {
    public:
        virtual ~MagneticFieldStrengthFringingEffectModel() = default;
        virtual FieldPoint get_equivalent_inducing_point_for_gap(CoreGap gap, double magneticFieldStrengthGap) = 0;
        virtual ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point(CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint) = 0;

        // A WINDOW-BOUNDARY model (ABT #1304) places the magnetic circuit's MMF drops on the winding-window boundary
        // instead of evaluating gaps one by one. Its sources are only consistent with a conductor field built the same
        // way -- one point per turn with the complete-cell image lattice, a turn's own images counted -- and
        // MagneticField builds it so when this returns true.
        virtual bool is_window_boundary_model() const { return false; }
        virtual void setup_window_boundary([[maybe_unused]] Magnetic magnetic, [[maybe_unused]] std::optional<OperatingPoint> operatingPoint) {
            throw ModelNotAvailableException("This fringing effect model has no window-boundary sources");
        }
        // Field at (x, y) of the boundary sources for a net winding MMF of netAmpereTurns, Binns-Lawrenson sign convention.
        virtual std::pair<double, double> get_window_boundary_field([[maybe_unused]] double pointX, [[maybe_unused]] double pointY,
                                                                    [[maybe_unused]] double netAmpereTurns) const {
            throw ModelNotAvailableException("This fringing effect model has no window-boundary sources");
        }
        // The same field per unit net MMF, split by the foil strips it crosses on its way to (x, y) (ABT #1645): returns
        // the part that crosses none and appends the screened parts to `screened`, grouped by the strips they cross.
        virtual std::pair<double, double> get_window_boundary_field_screened([[maybe_unused]] double pointX, [[maybe_unused]] double pointY,
                                                                             [[maybe_unused]] const std::vector<ScreeningStrip>& strips,
                                                                             [[maybe_unused]] std::optional<int64_t> inducedTurn,
                                                                             [[maybe_unused]] std::vector<ScreenedFieldSum>& screened) const {
            throw ModelNotAvailableException("This fringing effect model has no window-boundary sources");
        }
};

// ABT #1314: the physical current split between the parallels of each winding. The parallels of a winding are strands
// connected at both ends: every strand of a winding sees the same terminal voltage and the strand currents sum to the
// winding's current. A strand runs through all the turns of its parallel in series, so its voltage is the sum over its
// turn crossings of (R_t I + j w l_t A_t), A_t the 2D vector potential at crossing t (per unit length, averaged over the
// conductor's section) from every conductor, its core images and the window-wall sources that return the net MMF.
// Solved per harmonic as one small dense complex system; DC splits by resistance (WindingOhmicLosses) as before.
struct ParallelCurrentSplit {
    // The harmonics solved, as frequencies k f_excitation (k >= 1, ascending), and per harmonic every frequency the
    // windings' harmonics lists give it (a list built from a waveform's own time span can differ from k f_excitation).
    std::vector<double> frequencies;
    std::vector<std::vector<double>> listedFrequenciesPerFrequency;
    // Per solved frequency, per turn: the complex peak current of the turn's strand, MAS current frame (into the dot),
    // phases referred to the gauge winding as CoilMesher::calculate_current_phase_per_winding does.
    std::vector<std::vector<std::complex<double>>> currentPerTurnPerFrequency;
    // Per solved frequency, per winding: the winding's complex peak current (amplitude x exp(j phase)).
    std::vector<std::vector<std::complex<double>>> currentPerWindingPerFrequency;
    // The strand (winding, parallel) of every turn, numbered 0 .. numberStrands - 1, and the winding of every strand.
    std::vector<size_t> strandPerTurn;
    std::vector<size_t> windingPerStrand;
    // The index in `frequencies` of a harmonic frequency, solved or listed; throws when it was not solved.
    size_t frequency_index(double frequency) const;
    // The same, empty when the harmonic was not solved.
    std::optional<size_t> find_frequency_index(double frequency) const;
};

class MagneticField {
    private:
    protected:
        MagneticFieldStrengthModels _magneticFieldStrengthModel;
        MagneticFieldStrengthFringingEffectModels _magneticFieldStrengthFringingEffectModel;
        std::shared_ptr<MagneticFieldStrengthModel>  _model;
        std::shared_ptr<MagneticFieldStrengthFringingEffectModel>  _fringingEffectModel;
    public:
        std::vector<Wire> _wirePerWinding;
        std::vector<double> _wireMaxOuterWidth;  // Precomputed for performance
        std::vector<double> _wireMaxOuterHeight; // Precomputed for performance

        MagneticField(MagneticFieldStrengthModels magneticFieldStrengthModel = Defaults().magneticFieldStrengthModelDefault,
                      MagneticFieldStrengthFringingEffectModels magneticFieldStrengthFringingEffectModel = Defaults().magneticFieldStrengthFringingEffectModelDefault) {
            _magneticFieldStrengthModel = magneticFieldStrengthModel;
            _magneticFieldStrengthFringingEffectModel = magneticFieldStrengthFringingEffectModel;
            _model = factory(_magneticFieldStrengthModel);
            _fringingEffectModel = factory(_magneticFieldStrengthFringingEffectModel);
        }

        // frequency: the excitation's, so a magnetizing current given on a closed axis ([t0, t0 + T])
        // is told from a sampled open one when it is compressed (ABT #1467).
        static SignalDescriptor calculate_magnetic_flux(SignalDescriptor magnetizingCurrent,
                                                        double reluctance,
                                                        double numberTurns,
                                                        double frequency);
        static SignalDescriptor calculate_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                                double area);
        static SignalDescriptor calculate_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                    double initialPermeability);

        // Returns the in-phase field (inherited field_per_frequency, what every earlier reader
        // used) AND the quadrature field (see WindingWindowMagneticStrengthFieldPhasorOutput).
        // Default current directions follow the MAS excitation convention by isolation side
        // (CoilMesher::calculate_current_direction_per_winding); a custom vector overrides
        // them. Phases always come from the windings' current waveforms.
        WindingWindowMagneticStrengthFieldPhasorOutput calculate_magnetic_field_strength_field(OperatingPoint operatingPoint, Magnetic magnetic, std::optional<Field> externalInducedField = std::nullopt, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt, std::optional<CoilMesherModels> coilMesherModel = std::nullopt, std::optional<ParallelCurrentSplit> parallelCurrentSplit = std::nullopt);

        static std::shared_ptr<MagneticFieldStrengthFringingEffectModel> factory(MagneticFieldStrengthFringingEffectModels modelName);
        static std::shared_ptr<MagneticFieldStrengthModel> factory(MagneticFieldStrengthModels modelName);
        static std::shared_ptr<MagneticFieldStrengthModel> factory();

        // The turn-field sums kept between calls (Settings::magnetic_field_turn_sums_cache_bytes):
        // how often a mesh was found, how often it was summed afresh, and what is held now.
        struct TurnSumsCacheStatistics {
            size_t hits = 0;
            size_t misses = 0;
            size_t entries = 0;
            size_t bytes = 0;
        };
        // ABT #1314: the split of every winding's harmonic currents between its parallels (see ParallelCurrentSplit), at the
        // conductors' temperature. Empty when no winding has more than one parallel.
        static std::optional<ParallelCurrentSplit> calculate_parallel_current_split(Magnetic magnetic, OperatingPoint operatingPoint, double temperature,
                                                                                    std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt);
        static TurnSumsCacheStatistics get_turn_sums_cache_statistics();
        static void clear_turn_sums_cache();
};



// Based on Effects of eddy currents in transformer windings by P. L. Dowell.
// https://sci-hub.st/10.1049/piee.1966.0236
//
// ABT #376: Dowell is a ONE-DIMENSIONAL model and does not fit the point-to-point interface the
// way Lammeraner or Binns-Lawrenson do — there is no "field of one conductor at a point" in his
// formulation. He assumes the layers span the whole winding breadth b, so the field is parallel
// to the layers, uniform along them, and set purely by the MMF enclosed between the point and
// the zero-field boundary: H = (enclosed ampere-turns) / b. Expressed per inducing conductor
// that is exactly the MMF staircase — each conductor adds its own I/b to every point on ONE
// side of it and nothing on the other — which is what this implementation returns, so summing
// over conductors reproduces Dowell's staircase exactly.
//
// The consequence is worth stating: unlike the 2-D models this one has NO distance dependence,
// so it does not resolve fringing near a gap, and its accuracy relies on the layers being long
// against the window (Dowell's own assumption). It is the right model for a classic layered
// transformer window and the wrong one for a gapped inductor.
class MagneticFieldStrengthDowellModel : public MagneticFieldStrengthModel {
    public:
        std::string methodName = "Dowell";
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt);
};



// Based on Improved Analytical Calculation of High Frequency Winding Losses in Planar Inductors by Xiaohui Wang
// https://sci-hub.st/10.1109/ECCE.2018.8558397
class MagneticFieldStrengthWangModel : public MagneticFieldStrengthModel {
    public:
        std::string methodName = "Wang";
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt);
};


// Based on Analysis and Computation of Electric and Magnetic Field Problems By Kenneth John Binns and P. J. Lawrenson
// https://library.lol/main/78A1F466FA9F3EFBCB6165283FC346B6
// Equation 3.34 and 5.4
class MagneticFieldStrengthBinnsLawrensonModel : public MagneticFieldStrengthModel {
    public:
        std::string methodName = "BinnsLawrenson";
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt);
        std::pair<double, double> get_magnetic_field_strength_components_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt) override;
};


// Based on Eddy currents by Jiří Lammeraner
// https://archive.org/details/eddycurrents0000lamm
class MagneticFieldStrengthLammeranerModel : public MagneticFieldStrengthModel {
    public:
        std::string methodName = "Lammeraner";
        ComplexFieldPoint get_magnetic_field_strength_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt);
        std::pair<double, double> get_magnetic_field_strength_components_between_two_points(const FieldPoint& inducingFieldPoint, const FieldPoint& inducedFieldPoint, std::optional<size_t> inducingWireIndex = std::nullopt) override;
};


// Based on Fringing Field Formulas and Winding Loss Due to an Air Gap by Waseem A. Roshen
// https://sci-hub.st/10.1109/tmag.2007.898908
// https://sci-hub.st/10.1109/tmag.2008.2002302
class MagneticFieldStrengthRoshenModel : public MagneticFieldStrengthFringingEffectModel {
    public:
        std::string methodName = "Roshen";
        ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point(CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint);
        FieldPoint get_equivalent_inducing_point_for_gap([[maybe_unused]]CoreGap gap, [[maybe_unused]]double magneticFieldStrengthGap) {
            throw std::runtime_error("Fringing field not implemented for this model");
        }
};


// Based on Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre parasitären Eigenschaften by Manfred Albach
// https://libgen.rocks/get.php?md5=94b7f2906f53602f19892d7f1dabd929&key=YMKCEJOWB653PYLL
// https://sci-hub.st/10.1109/tpel.2011.2143729
class MagneticFieldStrengthAlbachModel : public MagneticFieldStrengthFringingEffectModel {
    public:
        std::string methodName = "Albach";
        ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point([[maybe_unused]]CoreGap gap, [[maybe_unused]]double magneticFieldStrengthGap, [[maybe_unused]]FieldPoint inducedFieldPoint) {
            throw std::runtime_error("Fringing field not implemented for this model");
        }
        FieldPoint get_equivalent_inducing_point_for_gap(CoreGap gap, double magneticFieldStrengthGap);
        // Albach's fitted polynomials are only valid for small xi = gapLength / columnDiameter
        // (the fitted denominator crosses zero at xi ~ 0.2755). Gaps outside this range are
        // routed through the Roshen conformal model by the field calculation instead of throwing.
        static bool is_gap_within_validity_range(CoreGap gap);
};

// ============================================================================
// SULLIVAN Fringing Field Model (2D Image Method / Biot-Savart)
// ============================================================================
/**
 * @brief Fringing field model based on the method of images (Biot-Savart approach)
 *
 * Based on the "shapeopt" MATLAB tool by C.R. Sullivan et al. at Dartmouth College.
 * Reference: http://thayer.dartmouth.edu/inductor
 *
 * Key publications:
 * - J. Hu, C.R. Sullivan, "Optimization of shapes for round-wire high-frequency
 *   gapped-inductor windings", IEEE IAS Annual Meeting, 1998
 * - C.R. Sullivan, "Optimal Choice for Number of Strands in a Litz-Wire
 *   Transformer Winding", IEEE Trans. Power Electron., 14(2), March 1999
 *
 * This model uses the 2D method of images to compute fringing magnetic field
 * strength (H) from an air gap. The gap MMF is modeled as distributed current
 * filaments across the gap width. Each filament is placed at the gap position
 * along with its image (reflected about the centerpost face), and these pairs
 * are replicated across image units of the winding window to enforce boundary
 * conditions imposed by the high-permeability core walls.
 *
 * The field is computed via 2D Biot-Savart superposition:
 *   B = (mu_0 * I) / (2*pi) * sum_images( (R - R_source) / |R - R_source|^2 )
 * converted to H = B / mu_0
 *
 * Unlike Roshen (closed-form) or Albach (equivalent current loop), this model
 * captures the effect of multiple reflections from all core walls, making it
 * more accurate for narrow winding windows or points close to core walls.
 */
class MagneticFieldStrengthSullivanModel : public MagneticFieldStrengthFringingEffectModel {
 public:
    std::string methodName = "Sullivan";

    /**
     * @brief Compute fringing H-field at an arbitrary point due to one gap
     *
     * Uses 2D Biot-Savart with method of images. The gap is subdivided into
     * _gapDivisions filaments, each replicated across (2*_imageUnitsX+1) x
     * (2*_imageUnitsY+1) image cells.
     *
     * @param gap                      CoreGap with coordinates, length, section_dimensions
     * @param magneticFieldStrengthGap Peak H in the gap (A/m), i.e. B_gap / mu_0
     * @param inducedFieldPoint        Point where the field is evaluated
     * @return ComplexFieldPoint with real = Hx (radial), imaginary = Hy (axial)
     */
    ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point(
        CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint);

    /**
     * @brief Not used - Sullivan computes field directly, not via equivalent point
     */
    FieldPoint get_equivalent_inducing_point_for_gap(
        [[maybe_unused]] CoreGap gap, [[maybe_unused]] double magneticFieldStrengthGap) {
        throw std::runtime_error("Sullivan fringing model computes field directly, not via equivalent point");
    }

 private:
    /// Number of subdivisions of the gap length (higher = more accurate, slower)
    int _gapDivisions = 10;
    /// Number of image unit repetitions in the X direction (perpendicular to gap)
    int _imageUnitsX = 2;
    /// Number of image unit repetitions in the Y direction (parallel to gap)
    int _imageUnitsY = 2;
};

// MMF_BOUNDARY (ABT #1304): the winding-window field as a charge-neutral image problem.
//
// The image lattice of a single rectangular window mirrors every conductor in the four core walls, taken as mu -> inf.
// A winding with net ampere-turns NI then sits in a closed box that cannot hold net current (Ampere round the box would
// need H inside the walls), and a truncated lattice carries a spurious, nearly uniform field. In the real core the NI
// drops across the gaps and along the ferrite. This model puts exactly those drops on the window boundary as uniform
// current sheets returning -NI: each gap's share of the MMF (its reluctance over the core's total, MKF's own reluctance
// composition) across its slot mouth, and the ferrite's share along the four walls in proportion to length x H with
// H = B / mu. Sheets and conductors are imaged by the same complete-cell lattice, so every cell is neutral and the sum
// converges; the tangential field just inside each wall equals that wall's MMF density, as the wall boundary condition
// requires. No fitted constant: geometry and the reluctance model only.
//
// It generalises SULLIVAN (gap MMF as imaged filaments) to every MMF drop in the circuit, made consistent with the
// conductors' own lattice. Single rectangular windows with the pure 2-D BINNS_LAWRENSON kernel only: it throws for
// toroids, multi-window cores, LAMMERANER (its finite-length factor breaks the lattice's neutrality), DOWELL (1-D, no
// images), WANG (its own two-point mesh) and ALBACH. WindingLosses picks BINNS_LAWRENSON for it when the field model is
// chosen automatically.
class MagneticFieldStrengthMmfBoundaryModel : public MagneticFieldStrengthFringingEffectModel {
  private:
    std::vector<WindowWallSheet> _sheets;
    // ABT #1645: per sheet, its centre and the normalised moments <zeta^k> (k = 0, 2, ..., 12) of a uniform line source,
    // for the multipole form of calculate_sheet_field far from the sheet (see get_window_boundary_field).
    struct SheetMultipole {
        double centreX;
        double centreY;
        double farRadiusSquared;
        double momentRe[7];
        double momentIm[7];
    };
    std::vector<SheetMultipole> _sheetMultipoles;
    bool _isSetUp = false;
    WindowGapMmf _gapMmf = WindowGapMmf::RELUCTANCE_SHARE;
  public:
    std::string methodName = "MmfBoundary";
    MagneticFieldStrengthMmfBoundaryModel() = default;
    // ROSHEN_COUNT_ONCE: the same construction with Roshen's gap field at each gap mouth. Roshen's closed form is
    // roshenGapFieldFactor H_g times the field of the mouth as a uniform sheet of density H_g lying on the mu -> inf face
    // (sheet plus coincident image, which the lattice supplies), so the mouth is that sheet carrying 0.9 H_g g, counted
    // once: the ferrite walls take NI minus the gaps' Roshen MMF (see CoilMesher::generate_window_wall_sheets).
    explicit MagneticFieldStrengthMmfBoundaryModel(WindowGapMmf gapMmf) : _gapMmf(gapMmf) {
        if (gapMmf == WindowGapMmf::ROSHEN_GAP_FIELD) {
            methodName = "RoshenCountOnce";
        }
    }
    ComplexFieldPoint get_magnetic_field_strength_between_gap_and_point([[maybe_unused]] CoreGap gap, [[maybe_unused]] double magneticFieldStrengthGap,
                                                                        [[maybe_unused]] FieldPoint inducedFieldPoint) {
        throw ModelNotAvailableException("MMF_BOUNDARY is not a per-gap model: the ferrite carries part of the MMF; use get_window_boundary_field");
    }
    FieldPoint get_equivalent_inducing_point_for_gap([[maybe_unused]] CoreGap gap, [[maybe_unused]] double magneticFieldStrengthGap) {
        throw ModelNotAvailableException("MMF_BOUNDARY is not a per-gap model: the ferrite carries part of the MMF; use get_window_boundary_field");
    }
    bool is_window_boundary_model() const override { return true; }
    void setup_window_boundary(Magnetic magnetic, std::optional<OperatingPoint> operatingPoint) override;
    std::pair<double, double> get_window_boundary_field(double pointX, double pointY, double netAmpereTurns) const override;
    std::pair<double, double> get_window_boundary_field_screened(double pointX, double pointY, const std::vector<ScreeningStrip>& strips,
                                                                 std::optional<int64_t> inducedTurn, std::vector<ScreenedFieldSum>& screened) const override;
    // Closed-form field of one uniform sheet (see the definition), Binns-Lawrenson sign convention.
    static std::pair<double, double> calculate_sheet_field(const WindowWallSheet& sheet, double current, double pointX, double pointY);
};



// ============================================================================
// ALBACH H-Field Model (Air Coil / Biot-Savart approach)
// ============================================================================

/**
 * @brief Structure representing a turn's position for the ALBACH H-field model
 * 
 * For round wires: width and height should be 0 (point filament approximation)
 * For rectangular wires: width = radial extent, height = axial extent
 *   The current is assumed uniformly distributed across the cross-section
 *   and is modeled using filamentary subdivision (multiple circular filaments)
 */
struct AlbachTurnPosition {
    double r;           ///< Radial position of conductor center (m)
    double z;           ///< Axial position of conductor center (m)
    double current;     ///< Current amplitude (A)
    size_t turnIndex;   ///< Index in the coil's turn list
    
    // For rectangular wire support (set to 0 for round wires)
    double width = 0;   ///< Radial extent of conductor cross-section (m)
    double height = 0;  ///< Axial extent of conductor cross-section (m)
    
    // For frequency-dependent current distribution
    double skinDepth = 1e9;  ///< Skin depth at current frequency (m), large default = uniform distribution

    /// @brief Check if this turn represents a rectangular conductor
    bool isRectangular() const { return width > 1e-10 && height > 1e-10; }
};

/**
 * @brief Magnetic field strength model using air coil (Biot-Savart) approach
 * 
 * Based on M. Albach, "Two-dimensional calculation of winding losses in transformers",
 * PESC 2000.
 * 
 * This model computes the H-field contribution from TURNS only using the
 * analytical Biot-Savart formula with elliptic integrals for circular current
 * filaments. Gap fringing effects are handled separately by the fringing
 * effect models (ROSHEN or ALBACH fringing).
 */
class MagneticFieldStrengthAlbach2DModel : public MagneticFieldStrengthModel {
public:
    std::string methodName = "Albach";
    
    /**
     * @brief Calculate H field between inducing and induced points
     * 
     * This model calculates field from all turns at once via calculateTotalFieldAtPoint().
     * This per-turn-pair method should not be called and will throw an error.
     */
    ComplexFieldPoint get_magnetic_field_strength_between_two_points(
        const FieldPoint& inducingFieldPoint,
        const FieldPoint& inducedFieldPoint, 
        std::optional<size_t> inducingWireIndex = std::nullopt
    ) override;
    
    /**
     * @brief Set the turn positions for field calculation
     */
    void setTurns(const std::vector<AlbachTurnPosition>& turns) {
        _turns = turns;
    }
    
    /**
     * @brief Update turn currents for a specific harmonic
     * This allows reusing the setup with different current values
     */
    void updateTurnCurrents(const std::vector<double>& currents) {
        for (size_t i = 0; i < _turns.size() && i < currents.size(); ++i) {
            _turns[i].current = currents[i];
        }
    }

    /**
     * @brief Update skin depths for all turns based on frequency
     * 
     * At high frequency, current concentrates at conductor edges due to skin effect.
     * 
     * @param skinDepth Skin depth at current frequency (m)
     */
    void updateSkinDepths(double skinDepth) {
        for (auto& turn : _turns) {
            turn.skinDepth = skinDepth;
        }
    }
    
    /**
     * @brief Calculate total field at an induced point from all turns
     * This is more efficient than the per-turn-pair interface
     */
    ComplexFieldPoint calculateTotalFieldAtPoint(FieldPoint inducedFieldPoint);
    
    /**
     * @brief Setup the model from a Magnetic component
     * 
     * Extracts turn positions from the coil and sets up for field calculation.
     * 
     * @param magnetic The magnetic component (core + coil)
     * @param wirePerWinding Wire definitions for each winding
     */
    void setupFromMagnetic(Magnetic magnetic, const std::vector<Wire>& wirePerWinding);
    
private:
    std::vector<AlbachTurnPosition> _turns;
    
    /**
     * @brief Calculate H field at point (r, z) from all turns using Biot-Savart
     * 
     * @param r Radial coordinate of field point (m)
     * @param z Axial coordinate of field point (m)
     * @return pair<H_r, H_z> Magnetic field components in A/m
     */
    std::pair<double, double> calculateMagneticField(double r, double z) const;
};

// ============================================================================
// IMAGED_MMF_SHEETS: turns and gap MMF sheets imaged together (ABT #1409)
// ============================================================================
/**
 * @brief Window field of the turns and of every functional gap, imaged as one system.
 *
 * A gap in an otherwise high-permeability core is equivalent to a current sheet on the
 * leg face next to the window, carrying the MMF the gap absorbs and OPPOSING the winding
 * (Ampere around the window). Turns and gap sheets are placed in the same method-of-images
 * lattice the CoilMesher uses for the core walls, so the window carries ~no net current and
 * the lattice sum is the ideal-core window field. Away from the gap plane the two core
 * halves are equipotential and the field there collapses; it concentrates at the gap plane.
 * Free-space fringing kernels (Roshen) superposed on a separately computed turn field
 * cannot reproduce that, which is what read R_ac/R_dc 887 against 2D FEM's 268 on a stack
 * of flat turns beside 3.9 mm spacer gaps.
 *
 * Everything is closed form: a flat conductor (rectangular, foil, planar) is a uniform-
 * current rectangle, a round/litz turn a line current, a gap sheet a uniform line segment.
 * Gap MMF per unit magnetizing current: flux N/R_total (the reluctance model's total),
 * shared between lateral columns in proportion to their gap conductances, times the gap's
 * own reluctance. Only functional (subtractive/additive) gaps carry a sheet (ABT #832).
 *
 * The own turn is excluded at its own points (its field there is skin effect); its images
 * are kept, being the core's response. Scope: single-window, non-toroidal cores. Anything
 * else throws.
 */
class MagneticFieldStrengthImagedMmfSheetsModel : public MagneticFieldStrengthModel {
  public:
    std::string methodName = "ImagedMmfSheets";

    struct Source {
        double x1, x2, y1, y2;             // extent; x1 == x2 for a sheet, both extents collapse for a filament
        double weight;                     // image multiplier (times turn divider and direction for turns)
        std::optional<size_t> turnIndex;   // set on turn sources
        bool isImage;                      // false only for the real (m = n = 0) copy
        size_t windingIndex;
    };

    ComplexFieldPoint get_magnetic_field_strength_between_two_points(
        const FieldPoint& inducingFieldPoint,
        const FieldPoint& inducedFieldPoint,
        std::optional<size_t> inducingWireIndex = std::nullopt) override;

    void setup(Magnetic magnetic, const std::vector<Wire>& wirePerWinding,
               const std::vector<double>& currentDividerPerTurn,
               const std::vector<int8_t>& currentDirectionPerWinding,
               double frequency);

    // Field of all turns (per winding current amplitudes, signed), own real turn excluded.
    std::pair<double, double> turns_field(double x, double y, const std::vector<double>& currentPerWinding,
                                          std::optional<size_t> excludedTurn, const std::vector<double>* multiplierPerTurn = nullptr) const;
    // Field of all gap sheets per ampere of magnetizing current.
    std::pair<double, double> gaps_field(double x, double y) const;

    static std::pair<double, double> rectangle_field(double x, double y, double x1, double x2, double y1, double y2, double current);
    static std::pair<double, double> vertical_sheet_field(double x, double y, double xs, double y1, double y2, double current);
    static std::pair<double, double> filament_field(double x, double y, double xs, double ys, double current);

    const std::vector<Source>& get_turn_sources() const { return _turnSources; }
    const std::vector<Source>& get_gap_sources() const { return _gapSources; }

  private:
    std::vector<Source> _turnSources;
    std::vector<Source> _gapSources;
    static std::pair<double, double> source_field(const Source& source, double x, double y, double current);
};

} // namespace OpenMagnetics
