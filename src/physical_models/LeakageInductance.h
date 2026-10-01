#pragma once
#include "Defaults.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include "support/CoilMesher.h"
#include "physical_models/MagneticShunt.h"
#include <MAS.hpp>
#include <array>

using namespace MAS;

namespace OpenMagnetics {

class LeakageInductance{
    public:

        LeakageInductance(){
        };
        virtual ~LeakageInductance() = default;


    // Method selection: "ReluctanceNetwork" for windings on different columns; "Shunt" when the
    // magnetic carries a shunt placed inWindow or betweenSections (MAS-RFC 0015); "Energy" otherwise.
    // onColumn shunts only fill column gaps (they enter the magnetizing network, see
    // MagnetizingInductance's Magnetic overload); outsideWindow shunts throw.
    LeakageInductanceOutput calculate_leakage_inductance(Magnetic magnetic, double frequency, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1);

    // "Shunt" method in detail (ABT #1176): the Energy method's winding leakage with the shunts
    // removed, minus the air each sheet displaces, plus the energy of each sheet's reluctance network
    // (MagneticShuntModel). Also returns the per-shunt flux density per ampere and, where the material
    // has complex permeability over the frequency, the sheet losses per ampere squared.
    // relativePermeabilityPerShunt, when given, replaces the materials' permeability (one entry per
    // magnetic.shunts element) for design exploration or to reproduce a published case; the
    // materials are then not read and no losses are reported.
    MagneticShuntLeakageResult calculate_shunt_leakage(Magnetic magnetic, double frequency, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1,
                                                       std::optional<std::vector<double>> relativePermeabilityPerShunt = std::nullopt);
    ComplexField calculate_leakage_magnetic_field(Magnetic magnetic, double frequency, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1);
    LeakageInductanceOutput calculate_leakage_inductance_all_windings(Magnetic magnetic, double frequency, size_t sourceIndex = 0, size_t harmonicIndex = 1);
    // In-phase field only (see calculate_magnetic_field_phasor) and the grid cell area.
    std::pair<ComplexField, double> calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt);
    // The window field as a phasor (MAS excitation convention, 2026-09-24): in-phase and
    // quadrature components on the same grid (WindingWindowMagneticStrengthFieldPhasorOutput).
    // Stored energy is the in-phase energy plus the quadrature energy.
    struct PhasorField {
        ComplexField inPhase;
        ComplexField quadrature;
        double dA;
    };
    PhasorField calculate_magnetic_field_phasor(OperatingPoint operatingPoint, Magnetic magnetic, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt);
    std::pair<size_t, size_t> calculate_number_points_needed_for_leakage(Coil coil);

    // Leakage magnetic-field energy (Joules, peak stored energy) in the winding window for an
    // arbitrary signed per-winding current vector. The sign of each entry encodes current
    // direction; the magnitude is the PEAK of a sinusoid in amperes (an entry of 1.0 is the
    // 1 A peak unit reference). The parameter name keeps its historical "Rms" spelling.
    // This is the quadratic-form primitive used to assemble the multi-winding leakage
    // inductance matrix; fringing is disabled internally (leakage field only), matching
    // calculate_leakage_inductance.
    double calculate_leakage_field_energy(Magnetic magnetic, const std::vector<double>& currentsRmsSigned, double frequency, size_t harmonicIndex = 1);

    // Full symmetric N×N leakage inductance matrix Λ (henries), in the per-winding physical
    // current basis, assembled from the energy quadratic form via polarization:
    //   Λ_aa = 2·W(e_a)/a²,  Λ_ab = [W(e_a+e_b) − W(e_a) − W(e_b)]/a²  (a: unit reference peak)
    // It is well-conditioned (contains NO magnetizing term — the core flux is excluded from
    // the window-energy integral). The diagonal Λ_aa is the self-leakage of winding a, the
    // off-diagonal Λ_ab the mutual leakage. For an ampere-turn-balanced pair it reproduces
    // calculate_leakage_inductance exactly. Cost: N + N(N−1)/2 field solves.
    std::vector<std::vector<double>> calculate_leakage_inductance_matrix(Magnetic magnetic, double frequency, size_t harmonicIndex = 1);

    // ---- Toroidal cores (CoreShapeFamily::T) --------------------------------------------------
    // The leakage field of a toroid lives in the air of the bore and outside the core. With a high core
    // permeability the ferrite decouples the two air regions, and the model is:
    //
    //  1. Ring-plane (2-D) energy per unit length, exact. Each bore crossing of a turn is a line current in
    //     the bore, imaged ONLY in the bore wall (radius B/2); each outer crossing is a line current outside
    //     the core, imaged ONLY in the outer wall (radius A/2). A line current I at z in a region bounded by a
    //     high-permeability circle of radius R has the images k·I at R²/conj(z) and −k·I at the centre
    //     (k = (µ−1)/(µ+1)), which zero the tangential H on the wall. The energy is the closed-form sum
    //     ½·Σ I_j·A(z_j) with each conductor's geometric mean radius in its self term (in-wire energy
    //     included). No grid is involved, so there is no grid truncation, far-point drop or grid convergence.
    //     A region with a net current (a single winding, or ampere-turns that do not balance) links the core:
    //     the wall then carries the magnetizing circulation Σ I_bore / (2π R), the bore has no centre image
    //     and the outside has −(1 + k)·Σ I at the centre. The energy includes that circulation between the
    //     conductors and the wall, so it does not depend on the unit of length (see the .cpp).
    //  2. Extrusion by the half turn length ℓ∞ = h + w of the conductor envelope (the rectangle through the
    //     conductor centres): each ring-plane crossing stands for half of its turn.
    //  3. Three-dimensional correction for the azimuthal MMF. The winding MMF V(φ) on the envelope (a
    //     staircase rising by the turn current at each turn) drives a field that also leaves through the top
    //     and bottom faces and closes far from the core. For each azimuthal harmonic m the exterior Laplace
    //     problem around the envelope (a body of revolution) is solved by finite elements in the meridian
    //     plane; it gives an effective height h_eff(m), the 3-D energy of the mode over its 2-D per-length
    //     energy. h_eff(m) → ℓ∞ as m → ∞ (the field then hugs the surface), so the correction
    //     Σ_m w_m·(h_eff(m) − ℓ∞), with w_m = µ0·|S_m|²/(π·m) the 2-D per-length energy of harmonic m of the
    //     MMF sheet and S_m = Σ_k I_k·exp(−i·m·φ_k), converges fast. It vanishes for windings whose
    //     ampere-turns balance at every angle (interleaved or layered windings).
    struct ToroidalLineCurrent {
        double x;
        double y;
        double current;
        double geometricMeanRadius;
    };
    struct ToroidalLeakageEnergy {
        double ringPlaneEnergyPerLength;      // J/m, bore + outside (step 1)
        double extrusionLength;               // m, ℓ∞ (step 2)
        double threeDimensionalCorrection;    // J (step 3)
        double energy;                        // J = ringPlaneEnergyPerLength·extrusionLength + threeDimensionalCorrection
        double envelopeInnerRadius;           // m
        double envelopeOuterRadius;           // m
        double envelopeHeight;                // m
    };
    static constexpr size_t TOROIDAL_LEAKAGE_NUMBER_MODES = 16;
    // Per-length energy (J/m) of line currents in one air region bounded by a high-permeability circle.
    // interiorRegion: the conductors are inside the circle (bore); otherwise they are outside it.
    static double calculate_ring_plane_region_energy_per_length(const std::vector<ToroidalLineCurrent>& conductors, double wallRadius, double imageFactor, bool interiorRegion);
    // H field (A/m) of the same line currents and their images at (x, y). Inside a conductor the in-wire field
    // of a uniform current is used (conductor radius = geometricMeanRadius·e^{1/4}).
    static std::array<double, 2> calculate_ring_plane_region_field(const std::vector<ToroidalLineCurrent>& conductors, double wallRadius, double imageFactor, bool interiorRegion, double x, double y);
    // 2-D per-length energy (J/m) of azimuthal harmonic m of the MMF sheet: µ0·|Σ_k I_k·e^{−i·m·φ_k}|²/(π·m).
    static double calculate_sheet_mode_energy_per_length(const std::vector<double>& angles, const std::vector<double>& currents, size_t mode);
    // Effective height (m) of azimuthal harmonic m for a rectangular envelope [innerRadius, outerRadius] × [−height/2, height/2]:
    // the exterior energy of the potential cos(mφ) prescribed on the envelope, divided by 2πm (its 2-D per-length value).
    // refinement scales the mesh density; farFieldFactor puts the zero-potential boundary at farFieldFactor·outerRadius.
    static double calculate_body_of_revolution_effective_height(double innerRadius, double outerRadius, double height, size_t mode, double refinement = 1.0, double farFieldFactor = 30.0);
    // The ring-plane line currents of a toroidal magnetic for signed currents per winding (one per winding, A):
    // the bore crossing of every turn carries its winding current / parallels, the outer crossing the opposite.
    // Turns of a winding with zero current are left out. Shared by the leakage energy and by the Painter's field.
    struct ToroidalRingPlaneConductors {
        std::vector<ToroidalLineCurrent> bore;
        std::vector<ToroidalLineCurrent> outer;
        std::vector<double> turnAngles;       // azimuth of each bore crossing, same order as bore
        std::vector<double> turnCurrents;     // current of each turn, same order as bore
        double innerWallRadius;               // m, B/2
        double outerWallRadius;               // m, A/2
        double coreHeight;                    // m, C times the number of stacks
        double imageFactor;                   // (µ − 1)/(µ + 1), initial permeability at the ambient temperature
    };
    static ToroidalRingPlaneConductors calculate_toroidal_ring_plane_conductors(Magnetic magnetic, const std::vector<double>& currentPerWinding);
    // H (A/m) at (x, y) of the ring-plane model: in the bore (r < B/2) and outside the core (r > A/2) it is
    // calculate_ring_plane_region_field of that region's crossings. In the ferrite (B/2 <= r <= A/2) it is the
    // field each single-interface image solution transmits into the permeable side: a line current I in air at z
    // next to a permeable region continues there as (1 − k)·I at z (Binns, Lawrenson and Trowbridge, ch. 3), and
    // the bore crossings' remaining k·Σ I_bore sits at the centre, so that the circulation in the ferrite is
    // Σ I_bore (Ampère; it is the magnetizing field of the core).
    static std::array<double, 2> calculate_toroidal_ring_plane_field(const ToroidalRingPlaneConductors& conductors, double x, double y);
    // Leakage energy (J) of a toroidal magnetic for signed peak currents per winding.
    static ToroidalLeakageEnergy calculate_toroidal_leakage_energy(Magnetic magnetic, const std::vector<double>& currentPerWinding);

    private:
        static constexpr double NEGLIGIBLE_CURRENT = 1e-9;
        static constexpr double SINUSOIDAL_PEAK_TO_PEAK = 2.0;
        static constexpr double SINUSOIDAL_DUTY_CYCLE = 0.5;
        static constexpr double SINUSOIDAL_OFFSET = 0.0;
        static constexpr double DEGREES_IN_CIRCLE = 360.0;
        static constexpr double MINIMUM_PRECISION_LEVEL = 1.0;

        OperatingPoint create_leakage_operating_point(Magnetic& magnetic, size_t sourceIndex, size_t destinationIndex, double frequency);
        // Build an operating point that drives every winding with the requested signed RMS
        // current (sign encodes direction). Used for arbitrary multi-winding excitations.
        OperatingPoint create_excitation_operating_point(Magnetic& magnetic, const std::vector<double>& currentsRmsSigned, double frequency);
        // Integrate the leakage magnetic-field energy over the winding window from a computed
        // field. Extracted from calculate_leakage_inductance so the energy core is reusable.
        double integrate_leakage_energy(Magnetic& magnetic, ComplexField& field, double dA);
        // Peak amplitude of the excitation's current harmonic at the frequency the field was solved at:
        // the amplitude the field model drove the turns with (throws when there is none).
        static double harmonic_peak_current_at_field_frequency(const OperatingPointExcitation& excitation, double fieldFrequency);
        std::pair<size_t, size_t> calculate_grid_points(Magnetic& magnetic, double frequency);
        static double cached_body_of_revolution_effective_height(double innerRadius, double outerRadius, double height, size_t mode);

};
} // namespace OpenMagnetics
