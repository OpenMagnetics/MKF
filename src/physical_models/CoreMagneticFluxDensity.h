#pragma once
#include "constructive_models/Magnetic.h"
#include "physical_models/InitialPermeability.h"
#include <MAS.hpp>
#include <array>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief One cell of the core's flux density map: a small piece of the core's cross-section, in the
 * Painter's 2D frame (metres; x across the core, y along the main column for a two-piece set, the ring
 * plane for a toroid).
 */
struct CoreFluxDensityCell {
    std::vector<std::array<double, 2>> polygon;
    std::string region;
    // The core volume the cell stands for (m3): its area times the core's depth, times the copies of the drawn part
    // the whole core holds (an E core is drawn from the main column's axis out, half of it).
    double volume;
    // The cell's material, an index into the map's (and the waveform's) materials: 0 is the core's own; a drum's ring
    // or shell of another grade, or a moulded body's regions, have theirs.
    size_t materialIndex = 0;
};

/**
 * @brief The flux density map of the core at one instant: |B| and the secant relative permeability
 * B / (mu0 H) of every cell, aligned with CoreMagneticFluxDensity::get_cells().
 */
struct CoreFluxDensityInstant {
    double magneticFlux;
    std::vector<double> fluxDensities;
    std::vector<double> relativePermeabilities;
    // The magnetomotive force across every region (A), aligned with CoreMagneticFluxDensity::get_regions(); for a
    // two-piece set the regions are in series around the loop, so their sum is the loop's ferrite MMF.
    std::vector<double> regionMagnetomotiveForces;
};

/**
 * @brief The core's flux density map over an operating point's waveform.
 * fluxDensities[sample][cell] is signed: positive while the flux in the wound column is positive.
 */
struct CoreFluxDensityWaveform {
    std::vector<CoreFluxDensityCell> cells;
    std::vector<double> time;
    std::vector<double> magneticFlux;
    std::vector<std::vector<double>> fluxDensities;
    // The sample whose |flux| is largest: the instant the core is nearest saturation.
    size_t peakSampleIndex;
    // |B| of every cell at that sample.
    std::vector<double> peakFluxDensities;
    // (max - min) / 2 of every cell's signed B over the waveform.
    std::vector<double> alternatingFluxDensityAmplitudes;
    double saturationFluxDensity;
    // The materials the cells' materialIndex points into (0: the core's own).
    std::vector<CoreMaterial> materials;
    // The operating point's frequency and core temperature the map was read at.
    double frequency;
    double temperature;
};

/**
 * @brief The core losses cell by cell: every cell's own B(t) priced by MKF's core-loss model (CoreLosses, the same model
 * cascade, span checks and DC-bias factor simulate uses), times the cell's volume.
 */
struct CoreLossDensityMap {
    // W/m3 of every cell, aligned with the waveform's cells.
    std::vector<double> volumetricLosses;
    // W: the sum of every cell's volumetric losses times its volume.
    double coreLosses;
    // m3: the sum of the cells' volumes.
    double volume;
    // Secondary materials with no loss data, priced at zero (CoreLosses' rule for a semi-shielded drum's shell and a
    // moulded body's regions).
    std::vector<std::string> unpricedMaterials;
};

/**
 * @class CoreMagneticFluxDensity
 * @brief Flux density inside the core, cell by cell, from the core's shape (ABT #2534). Analytical: no mesh solve.
 *
 * The core is cut into the regions of its magnetic path. A straight limb carries its flux uniformly. A
 * right-angle bend between two rectangular limbs is the Schwarz-Christoffel L channel already used for the
 * core's exact reluctance (right_angle_bend_squares, ABT #1766):
 *     dz/dzeta = (b/pi) sqrt((zeta - q) / (zeta + 1)) / zeta,  q = (a/b)^2,
 * integrated in closed form with s = sqrt((zeta - q)/(zeta + 1)):
 *     z(s) = (2/pi) [b atanh(s) - a atan(b s / a)],
 * the window's re-entrant corner at z = 0, the limb of width a along -y and the limb of width b along +x. Its
 * streamlines (arg zeta constant) and equipotentials (|zeta| constant) cut the bend into curvilinear squares:
 * equal flux per row, and |B| = flux per depth / (b |s|) in the linear limit, which crowds at the window corner
 * and vanishes at the core's outer corner. A toroid is cut into concentric rings, B ~ 1/r.
 *
 * Saturation is local: every cell reads its own permeability on the material's magnetisation curve
 * (InitialPermeability::get_magnetisation_curve, the datasheet B(H) of ABT #2470). A region is a series of slices
 * between equipotentials (the bend's lines of constant |zeta|, a ring's whole section); the cells of a slice are in
 * parallel and share its magnetomotive force, so a cell that saturates passes its flux to its neighbours in the same
 * slice: the flux steps around a saturating window corner and comes back past it. In the linear limit every
 * curvilinear square has the same reluctance and the map is the conformal field exactly.
 *
 * The flux in the wound column is imposed, the flux MKF itself computes for the operating point
 * (MagnetizingInductance::calculate_inductance_and_magnetic_flux_density: B(t) times Ae), so the map spreads the
 * flux MKF reports; it does not recompute it.
 *
 * Covered: every closed shape MAS describes, as the Painter draws it, and toroids.
 *  - Two-piece sets and pieces with their plate (E, planar E, EI, U, UI, UR, UT, C, ...): rectangular columns of one
 *    depth keep the core's own depth, its exact 2D geometry. Round or irregular posts (ETD, PQ, RM, P, EP, ...) take
 *    their area over their drawn width as depth, a bend between depths is mapped at the yoke's depth and stretched
 *    back, and the yoke's depth is the one whose linear map has MKF's effective le / Ae (an area-equivalent 2D map: the
 *    3D crowding where a round post meets a flat yoke is not resolved).
 *  - Moulded bodies, region by region in their own powders (no air region).
 *  - The drum family (drum, rod, shielded and semi-shielded drum, drum with a plate) as MKF's drum model has it: a
 *    uniform post, radial flanges with B ~ 1/rho, the ring or shell axial in its own grade; no axial fall-off.
 * Other cores throw.
 */
class CoreMagneticFluxDensity {
  public:
    // One parallel path of a slice: its length along the flux and width across it, the cells painted with its B, and,
    // for a gapped ring, the gap length in series with it. A region's end cut keeps lengthFraction of the path (0 when
    // it removes the path's cell entirely): the cut decides how much of the slice the region holds, not the field in
    // it, so the path's B comes from its whole length and only that fraction of the slice's MMF is the region's.
    struct Path {
        std::vector<size_t> cellIndexes;
        double length;
        double width;
        double gapLength = 0;
        double lengthFraction = 1;
    };
    struct Slice {
        std::vector<Path> paths;
    };
    struct Region {
        std::string name;
        std::vector<Slice> slices;
        double depth;
        // Share of the wound column's flux the region carries.
        double fluxShare;
        // The region's material (an index into the map's materials) and so its magnetisation curve.
        size_t materialIndex = 0;
    };

  private:
    std::vector<CoreFluxDensityCell> _cells;
    std::vector<Region> _regions;
    // The core's materials and their magnetisation curves at the map's temperature and frequency: 0 is the core's own;
    // regions built while _currentMaterialIndex is set take that one.
    std::vector<CoreMaterial> _materials;
    std::vector<std::shared_ptr<const InitialPermeability::MagnetisationCurve>> _curves;
    size_t _currentMaterialIndex = 0;
    double _temperature;
    std::optional<double> _frequency;
    size_t add_material(const CoreMaterial& material);
    size_t _numberRowsPerBend;
    // For every region, the first region with the very same slices, depth and flux share (a core's mirror-image corners):
    // its solution is the twin's, bit for bit, so calculate_flux_densities solves it once.
    std::vector<size_t> _regionTwins;
    // The canonical geometry of a bend (its curvilinear squares, clipped at the cuts), by (widthA, widthB, cutA, cutB):
    // a core's mirror-image corners are computed once and placed twice.
    struct BendGeometry;
    std::map<std::array<double, 4>, std::shared_ptr<const BendGeometry>> _bendGeometries;
    static BendGeometry compute_bend_geometry(double widthA, double widthB, double cutA, double cutB, size_t numberRows);
    void find_region_twins();

    // Every closed core of two limbs per side, as the Painter draws it: a two-piece set or a piece and its plate, of any
    // family (rectangular or round posts). yokeDepth: the plate's depth; nullopt for the core's own (rectangular columns)
    // or the one that makes the linear map's reluctance MKF's effective le / Ae (calibrate_closed_core).
    void build_closed_core(Core core, double yokeDepth);
    void calibrate_closed_core(Core core);
    bool closed_core_columns_rectangular(Core core);
    void build_toroid(Core core);
    void add_straight_region(const std::string& name, std::array<double, 2> origin, std::array<double, 2> alongFlux,
                             double length, std::array<double, 2> acrossFlux, double width, double depth, double fluxShare);
    void add_bend_region(const std::string& name, std::array<double, 2> corner, bool mirrorX, bool mirrorY, double widthA,
                         double widthB, double cutA, double cutB, double depth, double fluxShare, double limbDepth);

    using Curve = InitialPermeability::MagnetisationCurve;
    double field_strength(const Curve& curve, double fluxDensity) const;
    double flux_density(const Curve& curve, double fieldStrength) const;
    double incremental_permeability(const Curve& curve, double fieldStrength) const;
    double path_flux_density(const Curve& curve, const Path& path, double magnetomotiveForce) const;
    double solve_slice(const Curve& curve, const Slice& slice, double sliceFlux, double depth, std::vector<double>& pathFluxDensities) const;
    // The drum family, as paint_drum_family_core draws it (an axisymmetric half section): drum, rod, shielded drum,
    // semi-shielded drum, drum with a plate.
    void build_drum_family(Core core);

  public:
    /**
     * @param magnetic A magnetic whose core is processed and has its material.
     * @param temperature Core temperature (C) the magnetisation curve is read at.
     * @param frequency Frequency the magnetisation curve is read at, as MKF's magnetizing model reads it.
     * @param numberRowsPerBend Streamline rows each bend is cut into (resolution of the map).
     */
    CoreMagneticFluxDensity(const Magnetic& magnetic, double temperature, std::optional<double> frequency = std::nullopt,
                            size_t numberRowsPerBend = 12);
    /**
     * @brief The map of a core alone, its flux driven through the main column (the column a winding on the main
     * column links; the Magnetic constructor checks that every winding is there).
     */
    CoreMagneticFluxDensity(Core core, double temperature, std::optional<double> frequency = std::nullopt, size_t numberRowsPerBend = 12);

    const std::vector<CoreFluxDensityCell>& get_cells() const { return _cells; }
    const std::vector<Region>& get_regions() const { return _regions; }
    double get_saturation_flux_density() const { return _curves.front()->saturationFluxDensity; }
    const std::vector<CoreMaterial>& get_materials() const { return _materials; }

    /**
     * @brief |B| and mu_r of every cell when the wound column carries the flux magneticFlux (Wb, signed; the map is
     * odd in it).
     */
    CoreFluxDensityInstant calculate_flux_densities(double magneticFlux) const;

    /**
     * @brief |B| (T) at a point of the Painter frame for the wound-column flux magneticFlux. Throws outside the
     * core's ferrite.
     */
    double calculate_flux_density_at_point(double x, double y, double magneticFlux) const;

    /**
     * @brief The map over the operating point's flux waveform, as MKF's magnetizing model computes it.
     */
    static CoreFluxDensityWaveform calculate_waveform(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                      size_t numberRowsPerBend = 12);

    /**
     * @brief The core losses of the map: each cell's signed B(t) over the period is a flux density excitation of its
     * own, priced by CoreLosses::get_core_volumetric_losses at the waveform's frequency and temperature.
     */
    static CoreLossDensityMap calculate_core_losses(const Magnetic& magnetic, const CoreFluxDensityWaveform& waveform);
    static CoreLossDensityMap calculate_core_losses(const CoreFluxDensityWaveform& waveform);
    /**
     * @brief The map over a flux density waveform B(t) read over the core's reference area
     * (MagnetizingInductance::get_flux_density_reference_area), one period of the given frequency, at the map's
     * temperature.
     */
    static CoreFluxDensityWaveform calculate_waveform(const CoreMagneticFluxDensity& map, const Core& core,
                                                      const Waveform& fluxDensityWaveform, double frequency);
};

} // namespace OpenMagnetics
