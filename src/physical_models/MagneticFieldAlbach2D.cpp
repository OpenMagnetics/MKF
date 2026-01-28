#include "MagneticFieldAlbach2D.h"
#include "support/Utils.h"

#include <cmath>
#include <iostream>
#include <numbers>

using namespace std::numbers;

namespace OpenMagnetics {

namespace {
    // Local constant for vacuum permeability
    const double mu0 = Constants().vacuumPermeability;
    
    /**
     * @brief Calculate air coil vector potential at a point
     * Uses the formula from Eq. (5) in the paper
     */
    double calculateAirCoilPotentialLocal(
        const std::vector<Albach2DTurnPosition>& turns,
        double r, double z
    ) {
        double A = 0.0;
        
        for (const auto& turn : turns) {
            double r0 = turn.r;
            double z0 = turn.z;
            double I = turn.current;
            
            double deltaZ = z - z0;
            double sumR = r + r0;
            
            double denom = sumR * sumR + deltaZ * deltaZ;
            if (denom < 1e-20) continue;
            
            double k2 = 4 * r * r0 / denom;
            double k = std::sqrt(k2);
            
            if (k > 0.999999) k = 0.999999;
            
            if (r > 1e-15 && k > 1e-10) {
                double K_k = std::comp_ellint_1(k);
                double E_k = std::comp_ellint_2(k);
                double prefactor = (mu0 * I / pi) * std::sqrt(r0 / r);
                A += prefactor * ((1 - k2 / 2) * K_k - E_k) / k;
            }
        }
        
        return A;
    }
    
    /**
     * @brief Calculate radial derivative of air coil potential at a point
     * Uses numerical differentiation
     */
    double calculateAirCoilPotentialDerivativeR(
        const std::vector<Albach2DTurnPosition>& turns,
        double r, double z
    ) {
        double dr = std::max(1e-6, 1e-4 * r);
        double A_plus = calculateAirCoilPotentialLocal(turns, r + dr, z);
        double A_minus = calculateAirCoilPotentialLocal(turns, r - dr, z);
        return (A_plus - A_minus) / (2 * dr);
    }
    
    /**
     * @brief Calculate axial derivative of air coil potential at a point
     * Uses numerical differentiation
     */
    double calculateAirCoilPotentialDerivativeZ(
        const std::vector<Albach2DTurnPosition>& turns,
        double r, double z, double c
    ) {
        double dz = std::max(1e-6, 1e-4 * c);
        double A_plus = calculateAirCoilPotentialLocal(turns, r, z + dz);
        double A_minus = calculateAirCoilPotentialLocal(turns, r, z - dz);
        return (A_plus - A_minus) / (2 * dz);
    }
    
    /**
     * @brief Compute A1 coefficients from boundary conditions at z=0 and z=c
     * Following Albach's Equations (10)-(13)
     * 
     * From Eq. (10) and (11):
     * C1n is determined by the BC ∂A/∂z = 0 at z = c (symmetric mode)
     * D1n is determined by the BC ∂A/∂z = 0 at z = 0 (antisymmetric mode)
     * 
     * These require projecting ∂A_aircoil/∂z onto the radial eigenfunctions S1n(r)
     */
    void computeA1Coefficients(
        Albach2DSolution& solution,
        const Albach2DCoreGeometry& geometry,
        const std::vector<Albach2DTurnPosition>& turns
    ) {
        size_t nMax = solution.nMax;
        double a = geometry.a;
        double b = geometry.b;
        double c = geometry.c;
        
        solution.C1n.resize(nMax);
        solution.D1n.resize(nMax);
        solution.C1n.setZero();
        solution.D1n.setZero();
        solution.C10 = 0;
        solution.D10 = 0;
        
        if (a < 1e-10 || turns.empty()) return;
        
        // For each eigenvalue p1n, compute the projection integral
        // Eq. (10): C1n involves integral of (∂A_aircoil/∂z)|_{z=c} * S1n(r) * r dr from a to b
        // Eq. (11): D1n involves integral of (∂A_aircoil/∂z)|_{z=0} * S1n(r) * r dr from a to b
        
        const int numIntegrationPoints = 32;
        double dr = (b - a) / numIntegrationPoints;
        
        for (size_t n = 0; n < nMax; ++n) {
            double p1n = solution.p1n(n);
            
            // Compute normalization integral: ∫ r * S1n(r)² dr from a to b
            double normIntegral = 0;
            for (int i = 0; i <= numIntegrationPoints; ++i) {
                double r = a + i * dr;
                double J1_pr = std::cyl_bessel_j(1, p1n * r);
                double Y0_pa = std::cyl_neumann(0, p1n * a);
                double Y1_pr = std::cyl_neumann(1, p1n * r);
                double J0_pa = std::cyl_bessel_j(0, p1n * a);
                double S1n = J1_pr * Y0_pa - Y1_pr * J0_pa;
                
                double weight = (i == 0 || i == numIntegrationPoints) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
                normIntegral += weight * r * S1n * S1n;
            }
            normIntegral *= dr / 3.0;
            
            if (std::abs(normIntegral) < 1e-20) continue;
            
            // Compute projection integral at z = c for C1n
            double projIntegral_C = 0;
            for (int i = 0; i <= numIntegrationPoints; ++i) {
                double r = a + i * dr;
                double dAdz = calculateAirCoilPotentialDerivativeZ(turns, r, c, c);
                
                double J1_pr = std::cyl_bessel_j(1, p1n * r);
                double Y0_pa = std::cyl_neumann(0, p1n * a);
                double Y1_pr = std::cyl_neumann(1, p1n * r);
                double J0_pa = std::cyl_bessel_j(0, p1n * a);
                double S1n = J1_pr * Y0_pa - Y1_pr * J0_pa;
                
                double weight = (i == 0 || i == numIntegrationPoints) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
                projIntegral_C += weight * r * dAdz * S1n;
            }
            projIntegral_C *= dr / 3.0;
            
            // C1n coefficient - normalized and scaled by hyperbolic factor
            double coshFactor = std::cosh(p1n * c);
            if (std::abs(coshFactor) > 1e-10) {
                solution.C1n(n) = -projIntegral_C / (normIntegral * p1n * std::sinh(p1n * c));
            }
            
            // Compute projection integral at z = 0 for D1n
            double projIntegral_D = 0;
            for (int i = 0; i <= numIntegrationPoints; ++i) {
                double r = a + i * dr;
                double dAdz = calculateAirCoilPotentialDerivativeZ(turns, r, 0, c);
                
                double J1_pr = std::cyl_bessel_j(1, p1n * r);
                double Y0_pa = std::cyl_neumann(0, p1n * a);
                double Y1_pr = std::cyl_neumann(1, p1n * r);
                double J0_pa = std::cyl_bessel_j(0, p1n * a);
                double S1n = J1_pr * Y0_pa - Y1_pr * J0_pa;
                
                double weight = (i == 0 || i == numIntegrationPoints) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
                projIntegral_D += weight * r * dAdz * S1n;
            }
            projIntegral_D *= dr / 3.0;
            
            // D1n coefficient
            double sinhFactor = std::sinh(p1n * c);
            if (std::abs(sinhFactor) > 1e-10) {
                solution.D1n(n) = -projIntegral_D / (normIntegral * p1n * std::cosh(p1n * c));
            }
        }
    }
    
    /**
     * @brief Compute C2n coefficients from boundary condition at r=b
     * Following Albach's Equation (15)
     * 
     * At r = b (outer core boundary), ∂A/∂r = 0 for high-μ core
     * This gives: C2n from projection of ∂A_aircoil/∂r onto cos(p2n * z)
     */
    void computeC2nCoefficients(
        Albach2DSolution& solution,
        const Albach2DCoreGeometry& geometry,
        const std::vector<Albach2DTurnPosition>& turns
    ) {
        size_t nMax = solution.nMax;
        double b = geometry.b;
        double c = geometry.c;
        
        solution.C2n.resize(nMax);
        solution.D2n.resize(nMax);
        solution.C2n.setZero();
        solution.D2n.setZero();  // D2n computed separately in coupled system
        solution.C20 = 0;
        
        if (turns.empty()) return;
        
        // Eq. (15): C2n from integral of (∂A_aircoil/∂r)|_{r=b} * cos(p2n * z) dz from -c to c
        // Normalization: ∫cos²(p2n * z) dz from -c to c = c (for n > 0)
        
        const int numIntegrationPoints = 32;
        double dz = (2 * c) / numIntegrationPoints;
        
        for (size_t n = 0; n < nMax; ++n) {
            double p2n = solution.p2n(n);
            
            // Compute projection integral
            double projIntegral = 0;
            for (int i = 0; i <= numIntegrationPoints; ++i) {
                double z = -c + i * dz;
                double dAdr = calculateAirCoilPotentialDerivativeR(turns, b, z);
                double cosFactor = std::cos(p2n * z);
                
                double weight = (i == 0 || i == numIntegrationPoints) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
                projIntegral += weight * dAdr * cosFactor;
            }
            projIntegral *= dz / 3.0;
            
            // C2n = -projIntegral / (normalization * radial_factor)
            // For high-μ core, the radial factor involves modified Bessel functions
            // Simplified: C2n ≈ -projIntegral / c
            solution.C2n(n) = -projIntegral / c;
        }
    }
    
    /**
     * @brief Compute C30 from integral of vector potential at r=a over gap
     * Following Albach's Equation (19a):
     * 
     * C30 * (c_u - c_l) = (1/μ₀) * ∫[A_aircoil + A1 + A2]_{r=a} dz from c_l to c_u
     * 
     * For a high-μ core, the integral of A at the gap boundary gives the
     * average flux density, which relates to the gap field.
     */
    double computeC30FromIntegral(
        const Albach2DSolution& solution,
        const Albach2DCoreGeometry& geometry,
        const std::vector<Albach2DTurnPosition>& turns,
        double c_l,
        double c_u
    ) {
        double a = geometry.a;
        double gapHeight = c_u - c_l;
        
        if (a < 1e-10 || gapHeight < 1e-10 || turns.empty()) {
            return 0;
        }
        
        // Integrate A at r = a over the gap region [c_l, c_u]
        const int numIntegrationPoints = 16;
        double dz = gapHeight / numIntegrationPoints;
        
        double integral = 0;
        for (int i = 0; i <= numIntegrationPoints; ++i) {
            double z = c_l + i * dz;
            
            // A_aircoil contribution
            double A = calculateAirCoilPotentialLocal(turns, a, z);
            
            // A1 contribution (if computed)
            // For now, A1 is typically small for high-μ cores
            // A += sum over n of C1n * S1n(a) * cosh(p1n*z)/cosh(p1n*c) + D1n * ...
            // Note: S1n(a) = 0 by definition of the eigenfunction, so A1 contribution at r=a is 0
            
            // A2 contribution would be added here if non-zero
            
            double weight = (i == 0 || i == numIntegrationPoints) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
            integral += weight * A;
        }
        integral *= dz / 3.0;
        
        // C30 from Eq. (19a): C30 = (1/μ₀) * integral / gapHeight
        double C30 = integral / (mu0 * gapHeight);
        
        return C30;
    }
    
    /**
     * @brief Compute C3m coefficients from integral formula Eq. (19b)
     * 
     * From Albach Eq. (19b):
     * C3m * a * I1(p3m*a) / I0(p3m*a) = (1/μ₀) * ∫[A + r*∂A/∂r]_{r=a} * cos(p3m*(z-c_l)) dz
     * 
     * where the integral is over the gap region [c_l, c_u].
     * This gives the Fourier coefficients for the fringing field profile.
     */
    void computeC3mFromIntegral(
        const Albach2DSolution& solution,
        const Albach2DCoreGeometry& geometry,
        const std::vector<Albach2DTurnPosition>& turns,
        Albach2DGapSolution& gapSol
    ) {
        size_t mMax = gapSol.C3m.size();
        double a = geometry.a;
        double c_l = gapSol.c_lower;
        [[maybe_unused]] double c_u = gapSol.c_upper;
        double gapHeight = gapSol.height();
        
        // For very thin gaps or toroidal cores, use uniform field approx
        if (a < 1e-10 || gapHeight < 1e-6 || turns.empty()) {
            gapSol.C3m.setZero();
            return;
        }
        
        const int numIntegrationPoints = 16;
        double dz = gapHeight / numIntegrationPoints;
        
        for (size_t m = 0; m < mMax; ++m) {
            double p3m = gapSol.p3m(m);
            
            // Compute Fourier integral
            // ∫[A + r*∂A/∂r]_{r=a} * cos(p3m*(z-c_l)) dz
            double integral = 0;
            for (int i = 0; i <= numIntegrationPoints; ++i) {
                double z = c_l + i * dz;
                
                double A = calculateAirCoilPotentialLocal(turns, a, z);
                double dAdr = calculateAirCoilPotentialDerivativeR(turns, a, z);
                double f = A + a * dAdr;  // The integrand from boundary matching
                
                double cosFactor = std::cos(p3m * (z - c_l));
                
                double weight = (i == 0 || i == numIntegrationPoints) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
                integral += weight * f * cosFactor;
            }
            integral *= dz / 3.0;
            
            // Normalization: ∫cos²(p3m*(z-c_l)) dz from c_l to c_u = gapHeight/2
            integral *= 2.0 / gapHeight;
            
            // Bessel function ratio: I1(p3m*a) / I0(p3m*a)
            double I0_val = std::cyl_bessel_i(0, p3m * a);
            double I1_val = std::cyl_bessel_i(1, p3m * a);
            
            if (std::abs(I1_val) > 1e-15 && std::abs(I0_val) > 1e-15) {
                double besselRatio = I1_val / I0_val;
                // C3m * a * besselRatio = integral / μ₀
                // C3m = integral / (μ₀ * a * besselRatio)
                gapSol.C3m(m) = integral / (mu0 * a * besselRatio);
            } else {
                gapSol.C3m(m) = 0;
            }
        }
    }
    
    /**
     * @brief Compute D2n from gap boundary condition Eq. (21)
     * 
     * From Albach Eq. (21), D2n couples to C3m through the boundary matching
     * at r = a. The radial derivative of A2 must match that of A3:
     * 
     * D2n = (1/c) * Σ_gaps Σ_m [C3m * coupling_integral]
     * 
     * where the coupling integral involves cos(p2n*z) * cos(p3m*(z-c_l)) over the gap.
     */
    void computeD2nFromGapBoundary(
        Albach2DSolution& solution,
        const Albach2DCoreGeometry& geometry
    ) {
        size_t nMax = solution.nMax;
        double a = geometry.a;
        double c = geometry.c;
        
        if (a < 1e-10 || solution.gapSolutions.empty()) {
            solution.D2n.setZero();
            return;
        }
        
        // Check if all C3m are zero (uniform field)
        bool allGapsUniform = true;
        for (const auto& gapSol : solution.gapSolutions) {
            for (Eigen::Index m = 0; m < gapSol.C3m.size(); ++m) {
                if (std::abs(gapSol.C3m(m)) > 1e-15) {
                    allGapsUniform = false;
                    break;
                }
            }
            if (!allGapsUniform) break;
        }
        
        if (allGapsUniform) {
            solution.D2n.setZero();
            return;
        }
        
        // Compute D2n from boundary matching
        for (size_t n = 0; n < nMax; ++n) {
            double p2n = solution.p2n(n);
            double D2n_sum = 0;
            
            for (const auto& gapSol : solution.gapSolutions) {
                double c_l = gapSol.c_lower;
                [[maybe_unused]] double c_u = gapSol.c_upper;
                double gapHeight = gapSol.height();
                
                // For each C3m, compute the coupling integral
                for (Eigen::Index m = 0; m < gapSol.C3m.size(); ++m) {
                    if (std::abs(gapSol.C3m(m)) < 1e-15) continue;
                    
                    double p3m = gapSol.p3m(m);
                    double C3m = gapSol.C3m(m);
                    
                    // Coupling integral: ∫cos(p2n*z) * cos(p3m*(z-c_l)) dz from c_l to c_u
                    // Using product-to-sum: cos(A)*cos(B) = 0.5*[cos(A-B) + cos(A+B)]
                    // Integral becomes sum of sin terms
                    
                    double integral = 0;
                    const int numPoints = 16;
                    double dz = gapHeight / numPoints;
                    
                    for (int i = 0; i <= numPoints; ++i) {
                        double z = c_l + i * dz;
                        double cos_p2n = std::cos(p2n * z);
                        double cos_p3m = std::cos(p3m * (z - c_l));
                        
                        double weight = (i == 0 || i == numPoints) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
                        integral += weight * cos_p2n * cos_p3m;
                    }
                    integral *= dz / 3.0;
                    
                    // Bessel function factor: p3m * I1(p3m*a) / I0(p3m*a)
                    double I0_val = std::cyl_bessel_i(0, p3m * a);
                    double I1_val = std::cyl_bessel_i(1, p3m * a);
                    
                    if (std::abs(I0_val) > 1e-15) {
                        double besselFactor = p3m * I1_val / I0_val;
                        D2n_sum += C3m * besselFactor * integral;
                    }
                }
            }
            
            // Normalize by the projection onto cos(p2n*z)
            // Normalization: ∫cos²(p2n*z) dz from -c to c = c
            solution.D2n(n) = D2n_sum / c;
        }
    }
} // anonymous namespace

// ============================================================================
// MagneticFieldAlbach2DBoundaryValueSolver Implementation
// ============================================================================

Albach2DSolution MagneticFieldAlbach2DBoundaryValueSolver::solve(
    const Albach2DCoreGeometry& geometry,
    const std::vector<Albach2DTurnPosition>& turns,
    size_t nMax,
    size_t mMax
) {
    Albach2DSolution solution;
    solution.nMax = nMax;
    solution.mMax = mMax;
    
    double a = geometry.a;
    double b = geometry.b;
    double c = geometry.c;
    
    // Step 1: Find eigenvalues p1n from S1n(p1n * b, a) = 0
    solution.p1n = findEigenvaluesP1n(a, b, nMax);
    
    // Step 2: Calculate p2n eigenvalues (straightforward)
    solution.p2n.resize(nMax);
    for (size_t n = 1; n <= nMax; ++n) {
        solution.p2n(n - 1) = n * pi / c;
    }
    
    // Step 3: Set up gap region - SINGLE GAP ONLY
    // The Albach 2D model is formulated for a single gap. Multiple gaps would
    // require extending the boundary matching to multiple z-regions.
    size_t numGaps = geometry.gaps.size();
    if (numGaps > 1) {
        throw std::runtime_error("ALBACH_2D model only supports single gap. Use a different model for multiple gaps.");
    }
    solution.numGaps = numGaps;
    solution.gapSolutions.resize(numGaps);
    
    // For each gap, compute eigenvalues and initialize coefficients
    for (size_t gapIdx = 0; gapIdx < numGaps; ++gapIdx) {
        const auto& gap = geometry.gaps[gapIdx];
        auto& gapSol = solution.gapSolutions[gapIdx];
        
        gapSol.c_lower = gap.c_lower;
        gapSol.c_upper = gap.c_upper;
        double gapHeight = gapSol.height();
        
        // Eigenvalues for this gap region: m*pi/(c_u - c_l)
        gapSol.p3m.resize(mMax);
        for (size_t m = 1; m <= mMax; ++m) {
            gapSol.p3m(m - 1) = m * pi / gapHeight;
        }
        
        gapSol.C3m.resize(mMax);
        gapSol.C3m.setZero();
        gapSol.C30 = 0;
    }
    
    // Legacy single-gap interface (for backward compatibility)
    if (numGaps > 0) {
        solution.p3m = solution.gapSolutions[0].p3m;
        solution.C3m = solution.gapSolutions[0].C3m;
        solution.C30 = solution.gapSolutions[0].C30;
    } else {
        // No gaps - create empty eigenvalue vector
        solution.p3m.resize(mMax);
        solution.p3m.setZero();
        solution.C3m.resize(mMax);
        solution.C3m.setZero();
        solution.C30 = 0;
    }
    
    // Initialize coefficient vectors for winding region
    solution.C1n.resize(nMax);
    solution.D1n.resize(nMax);
    solution.C2n.resize(nMax);
    solution.D2n.resize(nMax);
    
    solution.C1n.setZero();
    solution.D1n.setZero();
    solution.C2n.setZero();
    solution.D2n.setZero();
    solution.C10 = 0;
    solution.D10 = 0;
    solution.C20 = 0;
    
    // Step 4: Solve the boundary value problem following Albach's method
    // 
    // The method from Albach's paper:
    // 1. Compute A1 coefficients (C1n, D1n) from BCs at z=0, z=c (Eqs. 10-13)
    // 2. Compute C2n from BC at r=b (Eq. 15)
    // 3. Compute C30 and C3m from Eq. (19) - depends on A2 (which has D2n)
    // 4. Compute D2n from Eq. (21) - depends on C3m
    // 5. Steps 3-4 form a coupled system - iterate until convergence
    
    // First, compute A1 coefficients from boundary conditions at z=0 and z=c
    // Following Eq. (10)-(13): integrate the boundary conditions over r from a to b
    computeA1Coefficients(solution, geometry, turns);
    
    // Compute C2n coefficients from boundary condition at r=b
    // Following Eq. (15): integrate over z from 0 to c
    computeC2nCoefficients(solution, geometry, turns);
    
    // Calculate total MMF (N*I)
    double totalNI = 0;
    for (const auto& turn : turns) {
        totalNI += turn.current;
    }
    
    if (numGaps > 0 && geometry.totalGapLength() > 1e-10) {
        // Solve the coupled C3m - D2n system iteratively
        // 
        // From Albach's paper:
        // - Eq. (19): C30 and C3m depend on integrals of (A_aircoil + A1 + A2) at r=a
        // - Eq. (21): D2n depends on C3m through the gap boundary condition
        // - A2 contains D2n terms, so this creates a coupled system
        //
        // However, for practical magnetic components, the gap field is dominated by
        // the magnetic circuit relationship: H_gap ≈ N*I / l_gap
        // 
        // The C3m coefficients from Eq. (19) describe the FRINGING field profile,
        // while C30 from the magnetic circuit gives the AVERAGE gap field.
        //
        // IMPLEMENTATION APPROACH:
        // 1. Use magnetic circuit for C30 (dominant term)
        // 2. Compute C3m from Eq. (19) for fringing profile
        // 3. Compute D2n from Eq. (21) based on C3m
        // 4. Optionally iterate for better accuracy
        
        const int maxIterations = 3;  // Usually converges in 1-2 iterations
        
        for (int iter = 0; iter < maxIterations; ++iter) {
            // For each gap, compute gap coefficients
            for (size_t gapIdx = 0; gapIdx < numGaps; ++gapIdx) {
                auto& gapSol = solution.gapSolutions[gapIdx];
                double c_l = gapSol.c_lower;
                double c_u = gapSol.c_upper;
                double gapHeight = gapSol.height();
                
                // C30 from magnetic circuit (dominant contribution)
                // For distributed gaps, MMF is split among them
                double MMF_per_gap = totalNI / numGaps;
                double C30_magnetic_circuit = MMF_per_gap / (2 * gapHeight);
                
                // C30 from Eq. (19a): integral of A at r=a over gap region
                // Note: This typically underestimates because it doesn't account for
                // the flux focusing effect of the high-μ core
                [[maybe_unused]] double C30_from_integral = computeC30FromIntegral(solution, geometry, turns, c_l, c_u);
                
                // Use magnetic circuit value as the primary source (physically correct)
                gapSol.C30 = C30_magnetic_circuit;
                
                // Compute C3m from Eq. (19b) for fringing field profile
                computeC3mFromIntegral(solution, geometry, turns, gapSol);
            }
            
            // Compute D2n from Eq. (21) based on the current C3m values
            computeD2nFromGapBoundary(solution, geometry);
        }
        
        // Update legacy interface with first gap's solution
        solution.C3m = solution.gapSolutions[0].C3m;
        solution.C30 = solution.gapSolutions[0].C30;
    }
    
    solution.isValid = true;
    return solution;
}

double MagneticFieldAlbach2DBoundaryValueSolver::calculateVectorPotential(
    const Albach2DSolution& solution,
    const Albach2DCoreGeometry& geometry,
    const std::vector<Albach2DTurnPosition>& turns,
    double r, double z
) const {
    if (!solution.isValid) {
        throw std::runtime_error("Albach2DSolution is not valid");
    }
    
    double A = 0.0;
    
    // Add air coil contribution
    A += calculateAirCoilPotential(turns, r, z);
    
    double a = geometry.a;
    double b = geometry.b;
    double c = geometry.c;
    
    // Add A1 contribution (radial Bessel expansion)
    // From Albach Eq. (8a): A1 = Σ S1n(r) * [C1n*cosh(p1n*z) + D1n*sinh(p1n*z)]
    // where S1n satisfies ∂S1n/∂r = 0 at r=a and S1n = 0 at r=b
    if (r >= a && a > 1e-10) {
        for (size_t n = 0; n < solution.nMax; ++n) {
            double p = solution.p1n(n);
            double S = S1n(p, r, a);
            
            // Guard against overflow for large p*c
            double cosh_pc = std::cosh(p * c);
            double sinh_pc = std::sinh(p * c);
            
            if (std::abs(cosh_pc) > 1e-10 && std::abs(sinh_pc) > 1e-10) {
                // Symmetric terms (cosh for zero derivative at z = ±c)
                double zFactor_C = std::cosh(p * z) / cosh_pc;
                // Antisymmetric terms (sinh for zero derivative at z = ±c)
                double zFactor_D = std::sinh(p * z) / sinh_pc;
                
                A += solution.C1n(n) * S * zFactor_C;
                A += solution.D1n(n) * S * zFactor_D;
            }
        }
    }
    
    // Add A2 contribution (z-direction Fourier expansion)
    // From Albach Eq. (8b): A2 = Σ R2n(r) * [C2n + D2n] * cos(p2n*z)
    // where R2n involves modified Bessel functions I0, K0
    if (r >= a && a > 1e-10) {
        for (size_t n = 0; n < solution.nMax; ++n) {
            double p2 = solution.p2n(n);
            double cosFactor = std::cos(p2 * z);
            
            // Radial function using modified Bessel functions
            // Full form: R2n(r) = I0(p2*r)*K0(p2*b) - K0(p2*r)*I0(p2*b)
            // This satisfies R2n = 0 at r=b (BC for high-mu core)
            // Simplified using log ratio for numerical stability when p2*r is small
            double I0_r = std::cyl_bessel_i(0, p2 * r);
            double K0_r = (p2 * r > 1e-10) ? std::cyl_bessel_k(0, p2 * r) : 0;
            double I0_b = std::cyl_bessel_i(0, p2 * b);
            double K0_b = std::cyl_bessel_k(0, p2 * b);
            double I0_a = std::cyl_bessel_i(0, p2 * a);
            double K0_a = std::cyl_bessel_k(0, p2 * a);
            
            // R2n(r) normalized to 1 at r=a
            double R2n_r = I0_r * K0_b - K0_r * I0_b;
            double R2n_a = I0_a * K0_b - K0_a * I0_b;
            double radialFunc = (std::abs(R2n_a) > 1e-15) ? R2n_r / R2n_a : 0;
            
            A += (solution.C2n(n) + solution.D2n(n)) * radialFunc * cosFactor;
        }
    }
    
    // DC terms (typically small for high-mu cores)
    // A += solution.C10 * lnRatio + solution.D10;
    
    // IMPORTANT: According to Albach's paper (Eq. 6):
    // - For r > a (winding region): A = A_aircoil + A1 + A2
    // - For r < a (gap region only): A = A3
    //
    // The gap effect on the winding region comes through the BOUNDARY CONDITIONS
    // that couple D2n to C3m (Eq. 21 in the 2001 paper). The gap field does NOT
    // directly add to A in the winding region - it's only indirectly felt through
    // the modified coefficients.
    //
    // Currently A1/A2 are disabled, so the field is just from A_aircoil.
    // The gap contribution (A3) is only added for r < a (inside center leg).
    
    double totalCurrentActual = 0;
    for (const auto& turn : turns) {
        totalCurrentActual += turn.current;
    }
    double numTurns = static_cast<double>(turns.size());
    double currentScale = (numTurns > 0) ? (totalCurrentActual / numTurns) : 1.0;
    
    // Only add A3 (gap contribution) for r < a (inside center leg)
    // According to Albach Eq. (16): A3 = μ₀ C30 r + Σ μ₀ C3m I0(p3m r) cos(p3m(z-c_l))
    // Note: It's I0 (order 0), not I1!
    if (r < a && a > 1e-10) {
        for (size_t gapIdx = 0; gapIdx < solution.numGaps; ++gapIdx) {
            const auto& gapSol = solution.gapSolutions[gapIdx];
            
            double c_l = gapSol.c_lower;
            double c_u = gapSol.c_upper;
            
            // Only apply A3 within the gap z-range
            if (z >= c_l && z <= c_u) {
                double C30_scaled = gapSol.C30 * currentScale;
                
                // DC term: A3_DC = μ₀ * C30 * r
                double A3_dc = mu0 * C30_scaled * r;
                A += A3_dc;
                
                // Higher-order terms: A3_m = μ₀ * C3m * I0(p3m r) * cos(p3m(z-c_l))
                for (size_t m = 0; m < solution.mMax; ++m) {
                    double p3 = gapSol.p3m(m);
                    double cosFactor = std::cos(p3 * (z - c_l));
                    double I0_val = std::cyl_bessel_i(0, p3 * r);
                    
                    A += mu0 * gapSol.C3m(m) * currentScale * I0_val * cosFactor;
                }
            }
        }
    }
    
    return A;
}

std::pair<double, double> MagneticFieldAlbach2DBoundaryValueSolver::calculateMagneticField(
    const Albach2DSolution& solution,
    const Albach2DCoreGeometry& geometry,
    const std::vector<Albach2DTurnPosition>& turns,
    double r, double z
) const {
    // Calculate H from A using Albach's Eq. (2):
    // H_r = -(1/(μ₀ r)) * ∂A/∂z
    // H_z = (1/(μ₀ r)) * ∂(rA)/∂r = (1/μ₀) * (A/r + ∂A/∂r)
    
    // Handle r near zero to avoid division by zero
    if (r < 1e-10) {
        // At r=0, by symmetry H_r = 0 and H_z can be computed from limit
        return {0.0, 0.0};
    }
    
    // Use numerical differentiation with reasonable step sizes
    // Too small steps cause numerical noise, too large steps miss spatial features
    // Use ~0.1% of characteristic dimensions, minimum 1 micron
    double dr = std::max(1e-6, 1e-4 * r);
    double dz = std::max(1e-6, 1e-4 * geometry.c);
    
    double A_center = calculateVectorPotential(solution, geometry, turns, r, z);
    double A_r_plus = calculateVectorPotential(solution, geometry, turns, r + dr, z);
    double A_r_minus = calculateVectorPotential(solution, geometry, turns, r - dr, z);
    double A_z_plus = calculateVectorPotential(solution, geometry, turns, r, z + dz);
    double A_z_minus = calculateVectorPotential(solution, geometry, turns, r, z - dz);
    
    double dA_dr = (A_r_plus - A_r_minus) / (2 * dr);
    double dA_dz = (A_z_plus - A_z_minus) / (2 * dz);
    
    // From Albach Eq. (2): H_r = -(1/(μ₀ r)) * ∂A/∂z
    double H_r = -dA_dz / (mu0 * r);
    // From Albach Eq. (2): H_z = (1/μ₀) * (A/r + ∂A/∂r)
    double H_z = (A_center / r + dA_dr) / mu0;
    
    return {H_r, H_z};
}

double MagneticFieldAlbach2DBoundaryValueSolver::calculateAirCoilPotential(
    const std::vector<Albach2DTurnPosition>& turns,
    double r, double z
) const {
    // Sum contributions from all current loops
    // Using the exact formula for a circular current loop (Eq. 5 in paper):
    // A_φ = (μ₀ I / π) * sqrt(r₀/r) * [(1-k²/2)K(k) - E(k)] / k
    // where k² = 4r₀r / [(r + r₀)² + (z - z₀)²]
    
    double A = 0.0;
    
    for (const auto& turn : turns) {
        double r0 = turn.r;
        double z0 = turn.z;
        double I = turn.current;
        
        // Avoid singularity at the wire location
        double deltaZ = z - z0;
        double sumR = r + r0;
        
        double denom = sumR * sumR + deltaZ * deltaZ;
        if (denom < 1e-20) {
            continue; // Skip self-contribution
        }
        
        double k2 = 4 * r * r0 / denom;
        double k = std::sqrt(k2);
        
        if (k > 0.999999) {
            k = 0.999999; // Avoid singularity
        }
        
        // Complete elliptic integrals
        double K_k = comp_ellint_1(k);
        double E_k = comp_ellint_2(k);
        
        // Vector potential from this loop
        if (r > 1e-15 && k > 1e-10) {
            double prefactor = (mu0 * I / pi) * std::sqrt(r0 / r);
            A += prefactor * ((1 - k2 / 2) * K_k - E_k) / k;
        }
    }
    
    return A;
}

double MagneticFieldAlbach2DBoundaryValueSolver::S1n(double p1n, double r, double a) const {
    // S1n(r) = J1(p1n*r) * Y0(p1n*a) - Y1(p1n*r) * J0(p1n*a)
    // This combination satisfies the BC ∂A/∂r = 0 at r = a
    
    // Use fast real-valued Bessel functions for performance
    double z_r = p1n * r;
    double z_a = p1n * a;
    
    double J1_r = bessel_j1_fast(z_r);
    double Y0_a = bessel_y0_fast(z_a);
    double Y1_r = bessel_y1_fast(z_r);
    double J0_a = bessel_j0_fast(z_a);
    
    return J1_r * Y0_a - Y1_r * J0_a;
}

Eigen::VectorXd MagneticFieldAlbach2DBoundaryValueSolver::findEigenvaluesP1n(
    double a, double b, size_t nMax
) const {
    // Find roots of S1n(p * b, a) = 0
    // These are the eigenvalues where the radial function satisfies BC at both r = a and r = b
    
    Eigen::VectorXd eigenvalues(nMax);
    
    // For efficiency, use approximate eigenvalues directly
    // The roots are approximately spaced by π/(b-a) for the radial problem
    double dp = pi / (b - a);
    
    for (size_t n = 0; n < nMax; ++n) {
        // Approximate eigenvalue: slightly offset from n*π/(b-a)
        eigenvalues(n) = (n + 1) * dp;
    }
    
    return eigenvalues;
}

double MagneticFieldAlbach2DBoundaryValueSolver::besselOrthogonalityIntegral(
    double p1n, double a, double b
) const {
    // Integral of r * S1n(r)^2 from a to b
    // This is needed for normalization of the eigenfunction expansion
    
    // Use numerical integration (Simpson's rule for better accuracy with fewer points)
    const int numPoints = 32;  // Reduced from 100, using Simpson's for accuracy
    double dr = (b - a) / numPoints;
    double integral = 0.0;
    
    // Simpson's 1/3 rule: integral = (h/3) * [f(a) + 4*f(a+h) + 2*f(a+2h) + 4*f(a+3h) + ... + f(b)]
    for (int i = 0; i <= numPoints; ++i) {
        double r = a + i * dr;
        double S = S1n(p1n, r, a);
        double f = r * S * S;
        
        double weight;
        if (i == 0 || i == numPoints) {
            weight = 1.0;
        } else if (i % 2 == 1) {
            weight = 4.0;
        } else {
            weight = 2.0;
        }
        integral += weight * f;
    }
    integral *= dr / 3.0;
    
    return integral;
}

// ============================================================================
// MagneticFieldStrengthAlbach2DModel Implementation
// ============================================================================

ComplexFieldPoint MagneticFieldStrengthAlbach2DModel::get_magnetic_field_strength_between_two_points(
    FieldPoint inducingFieldPoint, 
    FieldPoint inducedFieldPoint, 
    std::optional<size_t> inducingWireIndex
) {
    // ALBACH_2D model calculates field from all turns at once via calculateTotalFieldAtPoint()
    // This per-turn-pair method should not be called
    throw std::runtime_error("ALBACH_2D model does not support per-turn-pair field calculation. Use calculateTotalFieldAtPoint() instead.");
}

/**
 * @brief Calculate magnetic field from a single circular filament (current loop)
 * 
 * Uses the analytical formula for the magnetic field of a circular current loop
 * based on complete elliptic integrals of the first and second kind.
 * 
 * @param r Radial position of field point (m)
 * @param z Axial position of field point (m)
 * @param r0 Radial position of the current loop (m)
 * @param z0 Axial position of the current loop (m)
 * @param I Current in the loop (A)
 * @return std::pair<double, double> (H_r, H_z) field components in A/m
 */
static std::pair<double, double> calculateCircularFilamentField(
    double r, double z, double r0, double z0, double I
) {
    double H_r = 0.0;
    double H_z = 0.0;
    
    if (std::abs(I) < 1e-15) return {0.0, 0.0}; // No current
    
    double deltaZ = z - z0;
    double sumR = r + r0;
    double diffR = r - r0;
    
    double denom = sumR * sumR + deltaZ * deltaZ;
    if (denom > 1e-20 && r > 1e-15 && r0 > 1e-15) {
        double k2 = 4 * r * r0 / denom;
        double k = std::sqrt(k2);
        
        if (k > 0.999999) {
            k = 0.999999; // Avoid singularity near the wire
        }
        
        if (k > 1e-10) {
            // Complete elliptic integrals
            double K_k = std::comp_ellint_1(k);
            double E_k = std::comp_ellint_2(k);
            
            double sqrtDenom = std::sqrt(denom);
            double denomDiffR = diffR * diffR + deltaZ * deltaZ;
            
            if (denomDiffR > 1e-20) {
                double prefactor = I / (2 * std::numbers::pi);
                
                H_r = prefactor * deltaZ / (r * sqrtDenom) * 
                      (-K_k + E_k * (r0*r0 + r*r + deltaZ*deltaZ) / denomDiffR);
                
                H_z = prefactor / sqrtDenom * 
                      (K_k + E_k * (r0*r0 - r*r - deltaZ*deltaZ) / denomDiffR);
            }
        }
    }
    
    return {H_r, H_z};
}

/**
 * @brief Calculate magnetic field from a rectangular conductor using filamentary subdivision
 * 
 * Divides the rectangular cross-section into a grid of circular filaments,
 * each carrying a proportional fraction of the total current. The field
 * contributions from all filaments are summed.
 * 
 * Based on the approach described in Binns & Lawrenson (1973) and consistent
 * with Albach's treatment of distributed current in rectangular conductors.
 * 
 * @param r Radial position of field point (m)
 * @param z Axial position of field point (m)
 * @param r0 Radial position of conductor center (m)
 * @param z0 Axial position of conductor center (m)
 * @param width Radial extent of conductor (m)
 * @param height Axial extent of conductor (m)
 * @param I Total current in conductor (A)
 * @param numSubR Number of subdivisions in radial direction (default 3)
 * @param numSubZ Number of subdivisions in axial direction (default 3)
 * @return std::pair<double, double> (H_r, H_z) field components in A/m
 */
static std::pair<double, double> calculateRectangularConductorField(
    double r, double z, double r0, double z0, double width, double height, double I,
    int numSubR = 3, int numSubZ = 3
) {
    double H_r_total = 0.0;
    double H_z_total = 0.0;
    
    // Current per filament (uniformly distributed)
    double dI = I / (numSubR * numSubZ);
    
    // Create a grid of filaments across the rectangular cross-section
    for (int ir = 0; ir < numSubR; ++ir) {
        for (int iz = 0; iz < numSubZ; ++iz) {
            // Position of this filament (centered in its sub-cell)
            double r_filament = r0 - width/2 + width * (ir + 0.5) / numSubR;
            double z_filament = z0 - height/2 + height * (iz + 0.5) / numSubZ;
            
            // Ensure positive radius (filament must be on positive r side)
            if (r_filament <= 0) r_filament = 1e-10;
            
            auto [dH_r, dH_z] = calculateCircularFilamentField(r, z, r_filament, z_filament, dI);
            H_r_total += dH_r;
            H_z_total += dH_z;
        }
    }
    
    return {H_r_total, H_z_total};
}

ComplexFieldPoint MagneticFieldStrengthAlbach2DModel::calculateTotalFieldAtPoint(FieldPoint inducedFieldPoint) {
    // Calculate total field at the induced point from ALL turns at once
    // using the BVP solution which includes gap effects
    
    if (!_geometrySet) {
        throw std::runtime_error("Albach2D model requires geometry to be set via setCoreGeometry() before use");
    }
    
    // Ensure BVP solution is computed
    ensureSolutionValid();
    
    // Extract induced point coordinates - in 2D cross section: [0] = x (radial), [1] = y (axial)
    double r = std::abs(inducedFieldPoint.get_point()[0]);
    double z = inducedFieldPoint.get_point()[1];
    
    // Use the BVP solver's calculateMagneticField which includes:
    // - Air coil contribution from all turns
    // - A1 contribution (radial Bessel expansion)
    // - A2 contribution (z-direction Fourier expansion)
    // - Gap contributions (A3) from the C30 and C3m coefficients
    auto [H_r, H_z] = _solver.calculateMagneticField(_solution, _geometry, _turns, r, z);
    
    // Convert to 2D Cartesian: real = radial (Hx), imaginary = axial (Hy)
    ComplexFieldPoint result;
    result.set_real(H_r);
    result.set_imaginary(H_z);
    result.set_point(inducedFieldPoint.get_point());
    if (inducedFieldPoint.get_turn_index()) {
        result.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_label()) {
        result.set_label(inducedFieldPoint.get_label().value());
    }
    
    return result;
}

WindingWindowMagneticStrengthFieldOutput MagneticFieldStrengthAlbach2DModel::calculateFieldDistribution(
    Magnetic magnetic,
    OperatingPoint operatingPoint,
    double harmonicAmplitudeThreshold
) {
    // This method efficiently calculates the field at all points
    // by solving the BVP once and evaluating at multiple locations
    
    WindingWindowMagneticStrengthFieldOutput output;
    output.set_field_per_frequency({});
    output.set_method_used("Albach2D");
    output.set_origin(ResultOrigin::SIMULATION);
    
    // Extract geometry from the magnetic component
    auto core = magnetic.get_mutable_core();
    auto coil = magnetic.get_mutable_coil();
    
    // For a pot core, extract the winding window dimensions
    // This would need to be adapted based on the actual core shape API
    auto processedDescription = core.get_processed_description();
    if (!processedDescription) {
        throw std::runtime_error("Core processed description not available");
    }
    
    auto columns = processedDescription->get_columns();
    auto windingWindows = processedDescription->get_winding_windows();
    
    if (windingWindows.size() == 0) {
        throw std::runtime_error("No winding window found in core");
    }
    
    auto& windingWindow = windingWindows[0];
    
    // Set up geometry
    _geometry.a = columns[0].get_width() / 2; // Center leg radius
    _geometry.b = windingWindow.get_radial_height().value() + _geometry.a; // Outer radius
    _geometry.c = windingWindow.get_height().value() / 2;
    
    // Get gap information
    auto gapping = core.get_functional_description().get_gapping();
    if (!gapping.empty()) {
        // Assuming center leg gap
        _geometry.gapLength = gapping[0].get_length();
        // Gap position - assuming center for now
        _geometry.gapPositionZ = 0; 
    } else {
        _geometry.gapLength = 0;
        _geometry.gapPositionZ = 0;
    }
    
    _geometry.corePermeability = 1e6; // Approximate ferrite as very high permeability
    _geometry.computeDerivedValues();
    _geometrySet = true;
    
    // Extract turn positions from the coil
    _turns.clear();
    auto turnsDescription = coil.get_turns_description();
    if (!turnsDescription) {
        throw std::runtime_error("Coil turns description not available");
    }
    
    size_t turnIndex = 0;
    for (const auto& turn : turnsDescription.value()) {
        auto coords = turn.get_coordinates();
        Albach2DTurnPosition turnPos;
        // In the 2D cross-section, x coordinate represents radial position
        turnPos.r = std::abs(coords[0]); // Radial distance from center
        turnPos.z = coords[1];            // Axial position (y in cross-section)
        turnPos.current = 1.0; // Normalized; actual current applied later
        turnPos.turnIndex = turnIndex++;
        _turns.push_back(turnPos);
    }
    
    // Solve the boundary value problem
    _solution = _solver.solve(_geometry, _turns);
    _solutionValid = true;
    
    // Now evaluate at all mesh points needed for loss calculation
    // This part would integrate with the existing mesher and loss calculation
    
    return output;
}

void MagneticFieldStrengthAlbach2DModel::setupFromMagnetic(
    Magnetic magnetic, 
    const std::vector<Wire>& wirePerWinding,
    size_t nMax, 
    size_t mMax
) {
    auto& core = magnetic.get_mutable_core();
    
    if (!core.is_gap_processed()) {
        core.process_gap();
    }
    auto gapping = core.get_functional_description().get_gapping();
    
    // Determine if toroidal
    bool isToroidalCore = (core.get_type() == CoreType::TOROIDAL);
    
    // Get winding windows
    auto windingWindows = core.get_winding_windows();
    if (windingWindows.empty()) {
        throw std::runtime_error("Core has no winding windows for ALBACH_2D model");
    }
    auto& windingWindow = windingWindows[0];
    
    // Extract geometry
    double windowWidth = 0;
    double windowHeight = 0;
    double centerLegRadius = 0;
    
    if (isToroidalCore) {
        // For toroidal cores:
        // - radial_height = inner radius of winding window (hole radius = B/2)
        // - The winding is in the hole, so inner boundary is at r=0 (center)
        // - Outer boundary is at r = radial_height (inner edge of core)
        // For the ALBACH_2D model:
        // - a = 0 (center of hole)
        // - b = radial_height (inner radius of core = edge of winding area)
        // - c = height/2
        if (windingWindow.get_radial_height()) {
            windowWidth = windingWindow.get_radial_height().value();
        } else {
            throw std::runtime_error("Winding Window is missing radial height");
        }
        windowHeight = windowWidth * 2; // Full diameter for height calculation
        centerLegRadius = 0; // No center leg - winding is in the hole
    } else {
        // For pot cores and other rectangular window shapes
        windowWidth = windingWindow.get_width().value();
        windowHeight = windingWindow.get_height().value();
        
        // Get center leg radius (inner radius of winding window)
        auto columns = core.get_columns();
        for (auto& column : columns) {
            if (column.get_type() == ColumnType::CENTRAL) {
                centerLegRadius = column.get_width() / 2;
                break;
            }
        }
        
        // If no central column found, use first column
        if (centerLegRadius == 0 && !columns.empty()) {
            centerLegRadius = columns[0].get_width() / 2;
        }
    }
    
    // Get gap information - support for multiple distributed gaps
    _geometry.gaps.clear();
    for (auto& gap : gapping) {
        if (gap.get_type() == GapType::SUBTRACTIVE || gap.get_type() == GapType::ADDITIVE) {
            double gapLength = gap.get_length();
            double gapPositionZ = 0;
            if (gap.get_coordinates()) {
                gapPositionZ = gap.get_coordinates().value()[1];
            }
            _geometry.addGap(gapLength, gapPositionZ);
        }
    }
    
    // Get core permeability
    double corePermeability = core.get_initial_permeability();
    
    // Set up geometry
    _geometry.a = centerLegRadius;
    _geometry.b = centerLegRadius + windowWidth;
    _geometry.c = windowHeight / 2;
    _geometry.corePermeability = corePermeability;
    _geometry.isToroidal = isToroidalCore;
    _geometry.useImageCurrents = true;
    _geometry.computeDerivedValues();
    _geometrySet = true;

    // Get turns description
    if (!magnetic.get_coil().get_turns_description()) {
        throw std::runtime_error("Missing turns description in coil");
    }
    auto turns = magnetic.get_coil().get_turns_description().value();
    
    // Set up all turns from the coil
    _turns.clear();
    for (size_t turnIdx = 0; turnIdx < turns.size(); ++turnIdx) {
        auto& turn = turns[turnIdx];
        Albach2DTurnPosition albachTurn;
        
        // In 2D cross-section, x = radial, y = axial
        albachTurn.r = std::abs(turn.get_coordinates()[0]);
        albachTurn.z = turn.get_coordinates()[1];
        albachTurn.current = 1.0; // Will be scaled per harmonic
        albachTurn.turnIndex = turnIdx;

        // Get wire info for this turn to set dimensions for rectangular wires
        auto windingIndex = magnetic.get_mutable_coil().get_winding_index_by_name(turn.get_winding());
        if (windingIndex < wirePerWinding.size()) {
            auto& wire = wirePerWinding[windingIndex];
            if (wire.get_type() != WireType::ROUND && wire.get_type() != WireType::LITZ) {
                // Rectangular, foil, or planar wire - set dimensions for subdivision
                if (wire.get_conducting_width()) {
                    albachTurn.width = resolve_dimensional_values(wire.get_conducting_width().value());
                }
                if (wire.get_conducting_height()) {
                    albachTurn.height = resolve_dimensional_values(wire.get_conducting_height().value());
                }
            }
            // For round/litz wires, width and height stay at 0 (point filament)
        }
        _turns.push_back(albachTurn);
    }
    
    // Pre-solve the boundary value problem
    preSolve(nMax, mMax);
}

void MagneticFieldStrengthAlbach2DModel::ensureSolutionValid() {
    if (!_solutionValid && _geometrySet && !_turns.empty()) {
        _geometry.computeDerivedValues();
        _solution = _solver.solve(_geometry, _turns);
        _solutionValid = true;
    }
}

} // namespace OpenMagnetics
