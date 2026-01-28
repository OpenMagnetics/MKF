#pragma once

#include "Constants.h"
#include "Defaults.h"
#include "Models.h"
#include "MagneticField.h"
#include "constructive_models/Magnetic.h"
#include "processors/Inputs.h"
#include "support/CoilMesher.h"
#include <MAS.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <map>
#include <numbers>
#include <optional>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief Structure representing a single air gap in the core
 */
struct Albach2DGap {
    double length;      // Gap length (m)
    double positionZ;   // Z-position of gap center (m)
    double c_lower;     // z = c_l (lower edge of gap), derived
    double c_upper;     // z = c_u (upper edge of gap), derived
    
    void computeDerivedValues() {
        c_lower = positionZ - length / 2;
        c_upper = positionZ + length / 2;
    }
    
    double height() const { return c_upper - c_lower; }
};

/**
 * @brief Structure representing the geometry of a pot core or similar axisymmetric core
 * with a winding window and optional air gaps (supports multiple distributed gaps).
 * 
 * Based on the boundary value problem in Albach's paper:
 * "The influence of air gap size and winding position on the proximity losses 
 * in high frequency transformers" (PESC 2001)
 * 
 * Extended to support multiple gaps by solving a coupled BVP with separate
 * gap regions, each with its own eigenvalue expansion and boundary matching.
 * 
 * For toroidal cores:
 * - a = 0 (no center leg, winding is in the hole)
 * - b = inner radius of core (edge of winding area)
 * - c = height/2
 * - isToroidal = true
 * - Image currents are placed at r' = b²/r to satisfy ∂A/∂r = 0 at r = b
 */
struct Albach2DCoreGeometry {
    double a;           // Inner radius of winding window (center leg outer radius, 0 for toroidal)
    double b;           // Outer radius of winding window (inner core surface for toroidal)
    double c;           // Half-height of winding window (total height = 2c)
    double corePermeability; // Relative permeability of the core
    
    bool isToroidal = false;  // True for toroidal cores (no center leg)
    bool useImageCurrents = true; // Enable image currents for high-μ boundaries
    
    // Multiple gap support
    std::vector<Albach2DGap> gaps;  // Vector of air gaps (can be empty, single, or multiple)
    
    // Legacy single-gap interface (for backward compatibility)
    double gapLength = 0;    // Total air gap length (sum of all gaps)
    double gapPositionZ = 0; // Z-position of gap center (for single gap)
    double c_lower = 0;      // z = c_l (lower edge of gap)
    double c_upper = 0;      // z = c_u (upper edge of gap)
    
    /**
     * @brief Add a gap to the geometry
     */
    void addGap(double length, double positionZ) {
        Albach2DGap gap;
        gap.length = length;
        gap.positionZ = positionZ;
        gap.computeDerivedValues();
        gaps.push_back(gap);
    }
    
    /**
     * @brief Check if this geometry has any gaps
     */
    bool hasGaps() const { return !gaps.empty() && gaps[0].length > 1e-10; }
    
    /**
     * @brief Get total gap length (sum of all gaps)
     */
    double totalGapLength() const {
        double total = 0;
        for (const auto& gap : gaps) {
            total += gap.length;
        }
        return total;
    }
    
    void computeDerivedValues() {
        // Compute derived values for all gaps
        for (auto& gap : gaps) {
            gap.computeDerivedValues();
        }
        
        // For backward compatibility, also set single-gap values
        if (!gaps.empty()) {
            gapLength = totalGapLength();
            gapPositionZ = gaps[0].positionZ;
            c_lower = gaps[0].c_lower;
            c_upper = gaps[0].c_upper;
        } else {
            // Legacy single-gap mode
            c_lower = gapPositionZ - gapLength / 2;
            c_upper = gapPositionZ + gapLength / 2;
        }
    }
};

/**
 * @brief Structure representing a turn's position in cylindrical coordinates
 * 
 * For round wires: width and height should be 0 (point filament approximation)
 * For rectangular wires: width = radial extent, height = axial extent
 *   The current is assumed uniformly distributed across the cross-section
 *   and is modeled using filamentary subdivision (multiple circular filaments)
 */
struct Albach2DTurnPosition {
    double r;           ///< Radial position of conductor center (m)
    double z;           ///< Axial position of conductor center (m)
    double current;     ///< Current amplitude (A)
    size_t turnIndex;   ///< Index in the coil's turn list
    
    // For rectangular wire support (set to 0 for round wires)
    double width = 0;   ///< Radial extent of conductor cross-section (m)
    double height = 0;  ///< Axial extent of conductor cross-section (m)
    
    // For frequency-dependent current distribution (Wang 2018)
    double skinDepth = 1e9;  ///< Skin depth at current frequency (m), large default = uniform distribution

    /// @brief Check if this turn represents a rectangular conductor
    bool isRectangular() const { return width > 1e-10 && height > 1e-10; }
};

/**
 * @brief Coefficients for a single gap region in the multi-gap BVP solution
 */
struct Albach2DGapSolution {
    double C30;                     // DC term for this gap's A3
    Eigen::VectorXd C3m;            // Modified Bessel coefficients for this gap's A3
    Eigen::VectorXd p3m;            // Eigenvalues for this gap: m*pi/(c_u - c_l)
    double c_lower;                 // Lower z boundary of this gap
    double c_upper;                 // Upper z boundary of this gap
    
    double height() const { return c_upper - c_lower; }
};

/**
 * @brief Solution coefficients for the Albach 2D boundary value problem
 * 
 * The vector potential is expressed as:
 * 
 * Region 1 (winding area): A = A_aircoil + A1 + A2
 *   where A1 uses radial Bessel functions, A2 uses z-direction Fourier
 * 
 * Region 2, 3, ... (air gaps): A = A3_i for each gap i
 *   uses modified Bessel functions, with separate coefficients per gap
 * 
 * For multiple distributed gaps, each gap region has its own set of
 * coefficients and eigenvalues. The boundary conditions at each gap
 * boundary (z = c_l_i, z = c_u_i) couple the winding region solution
 * to each gap region solution.
 */
struct Albach2DSolution {
    // Air coil contribution (from superposition of current loops)
    // Stored separately, computed analytically
    
    // Region 1: Winding area coefficients
    double C10, D10;                // DC terms for A1
    Eigen::VectorXd C1n;            // Bessel function coefficients for A1
    Eigen::VectorXd D1n;            // Bessel function coefficients for A1
    
    double C20;                     // DC term for A2
    Eigen::VectorXd C2n;            // Fourier coefficients for A2
    Eigen::VectorXd D2n;            // Fourier coefficients for A2
    
    // Multi-gap support: each gap has its own solution coefficients
    std::vector<Albach2DGapSolution> gapSolutions;
    
    // Legacy single-gap interface (for backward compatibility)
    double C30;                     // DC term for A3 (first gap)
    Eigen::VectorXd C3m;            // Modified Bessel coefficients for A3 (first gap)
    
    // Eigenvalues
    Eigen::VectorXd p1n;            // Eigenvalues for radial expansion (from S1n roots)
    Eigen::VectorXd p2n;            // Eigenvalues for z expansion: n*pi/c
    Eigen::VectorXd p3m;            // Eigenvalues for gap region: m*pi/(c_u - c_l) (first gap)
    
    size_t nMax;                    // Number of terms in expansions
    size_t mMax;                    // Number of terms per gap expansion
    size_t numGaps = 0;             // Number of gap regions
    
    bool isValid = false;
};

/**
 * @brief Albach 2D Boundary Value Solver for magnetic field calculation
 * 
 * Implements the analytical solution from:
 * [1] M. Albach, H. Rossmanith, "The influence of air gap size and winding position
 *     on the proximity losses in high frequency transformers", PESC 2001
 * [2] M. Albach, "Two-dimensional calculation of winding losses in transformers",
 *     PESC 2000
 * 
 * Extended to support multiple distributed gaps by solving a coupled BVP
 * with separate gap regions, each with its own eigenvalue expansion.
 * 
 * For N gaps, the linear system couples:
 * - Region 1 (winding area) with N boundary condition sets
 * - N gap regions, each with its own coefficients
 * 
 * The boundary conditions at each gap boundary (z = c_l_i, z = c_u_i) are:
 * 1. Continuity of A (vector potential)
 * 2. Continuity of tangential H (related to ∂A/∂r)
 * 
 * The solver computes the vector potential A(r,z) in cylindrical coordinates
 * for a pot core geometry with an arbitrary air gap, then derives the magnetic
 * field strength H from H = (1/μ₀) * curl(A).
 */
class MagneticFieldAlbach2DBoundaryValueSolver {
public:
    MagneticFieldAlbach2DBoundaryValueSolver() = default;
    
    /**
     * @brief Solve the boundary value problem for the given geometry and turns
     * 
     * @param geometry Core geometry parameters
     * @param turns Vector of turn positions and currents
     * @param nMax Maximum number of terms for Bessel/Fourier expansions (default 50)
     * @param mMax Maximum number of terms for gap region expansion (default 50)
     * @return Solution coefficients
     */
    Albach2DSolution solve(
        const Albach2DCoreGeometry& geometry,
        const std::vector<Albach2DTurnPosition>& turns,
        size_t nMax = 50,
        size_t mMax = 50
    );
    
    /**
     * @brief Calculate the vector potential A at a point (r, z)
     * 
     * @param solution Previously computed solution
     * @param geometry Core geometry
     * @param turns Turn positions (for air coil contribution)
     * @param r Radial coordinate
     * @param z Axial coordinate
     * @return Vector potential A (only φ component, as A = φ̂ * A_φ)
     */
    double calculateVectorPotential(
        const Albach2DSolution& solution,
        const Albach2DCoreGeometry& geometry,
        const std::vector<Albach2DTurnPosition>& turns,
        double r, double z
    ) const;
    
    /**
     * @brief Calculate the magnetic field H at a point (r, z)
     * 
     * H_r = -(1/μ₀) * ∂A/∂z
     * H_z = (1/μ₀r) * ∂(rA)/∂r
     * 
     * @param solution Previously computed solution
     * @param geometry Core geometry
     * @param turns Turn positions
     * @param r Radial coordinate
     * @param z Axial coordinate
     * @return pair<H_r, H_z> Magnetic field components
     */
    std::pair<double, double> calculateMagneticField(
        const Albach2DSolution& solution,
        const Albach2DCoreGeometry& geometry,
        const std::vector<Albach2DTurnPosition>& turns,
        double r, double z
    ) const;
    
private:
    /**
     * @brief Calculate air coil vector potential from all current loops
     * Uses the formula from Eq. (5) in the paper
     */
    double calculateAirCoilPotential(
        const std::vector<Albach2DTurnPosition>& turns,
        double r, double z
    ) const;
    
    /**
     * @brief Calculate the S1n Bessel function combination
     * S1n(r) = J1(p1n*r)*Y0(p1n*a) - Y1(p1n*r)*J0(p1n*a)
     */
    double S1n(double p1n, double r, double a) const;
    
    /**
     * @brief Find eigenvalues p1n from the characteristic equation
     * S1n(p1n, b, a) = 0
     */
    Eigen::VectorXd findEigenvaluesP1n(double a, double b, size_t nMax) const;
    
    /**
     * @brief Calculate the orthogonality integral for Bessel functions
     */
    double besselOrthogonalityIntegral(double p1n, double a, double b) const;
    
    // Constants
    const double mu0 = Constants().vacuumPermeability;
};


/**
 * @brief Magnetic field strength model using the full 2D Albach boundary value solution
 * 
 * This is a complete 2D analytical solution that properly handles:
 * - Arbitrary air gap size and position
 * - Two-dimensional field distribution in the winding window
 * - Proper boundary conditions at ferrite surfaces
 * 
 * Unlike the simplified Albach fringing model (equivalent current point),
 * this solves the full Laplace equation with proper boundary conditions.
 */
class MagneticFieldStrengthAlbach2DModel : public MagneticFieldStrengthModel {
public:
    std::string methodName = "Albach2D";
    
    /**
     * @brief Calculate H field between inducing and induced points
     * 
     * For this 2D model, we need to solve the full boundary value problem
     * for all inducing turns, then evaluate at the induced point.
     * 
     * This is called per-point, but internally caches the solution per geometry.
     */
    ComplexFieldPoint get_magnetic_field_strength_between_two_points(
        FieldPoint inducingFieldPoint, 
        FieldPoint inducedFieldPoint, 
        std::optional<size_t> inducingWireIndex = std::nullopt
    ) override;
    
    /**
     * @brief Calculate the complete H field distribution for all turns
     * 
     * This is more efficient than calling the point-to-point method,
     * as it solves the BVP once and evaluates at all points.
     */
    WindingWindowMagneticStrengthFieldOutput calculateFieldDistribution(
        Magnetic magnetic,
        OperatingPoint operatingPoint,
        double harmonicAmplitudeThreshold
    );
    
    /**
     * @brief Set the core geometry for this model
     */
    void setCoreGeometry(const Albach2DCoreGeometry& geometry) {
        _geometry = geometry;
        _geometrySet = true;
        _solutionValid = false;
    }
    
    /**
     * @brief Set the turn positions
     */
    void setTurns(const std::vector<Albach2DTurnPosition>& turns) {
        _turns = turns;
        _solutionValid = false;
    }
    
    /**
     * @brief Update turn currents for a specific harmonic
     * This allows reusing the BVP solution with different current values
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
     * This updates the skin depth used for edge-weighted current distribution
     * in rectangular conductors (Wang 2018 model).
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
     * @brief Pre-solve the boundary value problem
     * Call this after setting geometry and turns to solve the BVP once
     * @param nMax Number of terms in Bessel/Fourier expansion (default 10 for speed)
     */
    void preSolve(size_t nMax = 10, size_t mMax = 10) {
        if (_geometrySet && !_turns.empty()) {
            _geometry.computeDerivedValues();
            _solution = _solver.solve(_geometry, _turns, nMax, mMax);
            _solutionValid = true;
        }
    }
    
    /**
     * @brief Setup the model from a Magnetic component
     * 
     * Extracts core geometry, sets up turn positions, and pre-solves the BVP.
     * This encapsulates all the setup logic in one place.
     * 
     * @param magnetic The magnetic component (core + coil)
     * @param wirePerWinding Wire definitions for each winding
     * @param nMax Number of terms for Bessel/Fourier expansion (default 10)
     * @param mMax Number of terms for gap expansion (default 10)
     */
    void setupFromMagnetic(Magnetic magnetic, 
                           const std::vector<Wire>& wirePerWinding,
                           size_t nMax = 10, 
                           size_t mMax = 10);
    
private:
    MagneticFieldAlbach2DBoundaryValueSolver _solver;
    Albach2DCoreGeometry _geometry;
    std::vector<Albach2DTurnPosition> _turns;
    Albach2DSolution _solution;
    bool _geometrySet = false;
    bool _solutionValid = false;
    
    void ensureSolutionValid();
};

} // namespace OpenMagnetics
