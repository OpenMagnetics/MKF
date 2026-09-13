#pragma once
#include "Defaults.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include "support/CoilMesher.h"
#include "physical_models/MagneticShunt.h"
#include <MAS.hpp>

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
    std::pair<ComplexField, double> calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t sourceIndex = 0, size_t destinationIndex = 1, size_t harmonicIndex = 1, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt);
    std::pair<size_t, size_t> calculate_number_points_needed_for_leakage(Coil coil);

    // Leakage magnetic-field energy (Joules, per unit reference amplitude) stored in the
    // winding window for an arbitrary signed per-winding RMS current vector. The sign of
    // each entry encodes current direction; the magnitude is the RMS current in amperes.
    // This is the quadratic-form primitive used to assemble the multi-winding leakage
    // inductance matrix; fringing is disabled internally (leakage field only), matching
    // calculate_leakage_inductance.
    double calculate_leakage_field_energy(Magnetic magnetic, const std::vector<double>& currentsRmsSigned, double frequency, size_t harmonicIndex = 1);

    // Full symmetric N×N leakage inductance matrix Λ (henries), in the per-winding physical
    // current basis, assembled from the energy quadratic form via polarization:
    //   Λ_aa = 4·W(e_a),  Λ_ab = 2·[W(e_a+e_b) − W(e_a) − W(e_b)]
    // It is well-conditioned (contains NO magnetizing term — the core flux is excluded from
    // the window-energy integral). The diagonal Λ_aa is the self-leakage of winding a, the
    // off-diagonal Λ_ab the mutual leakage. For an ampere-turn-balanced pair it reproduces
    // calculate_leakage_inductance exactly. Cost: N + N(N−1)/2 field solves.
    std::vector<std::vector<double>> calculate_leakage_inductance_matrix(Magnetic magnetic, double frequency, size_t harmonicIndex = 1);

    private:
        static constexpr double NEGLIGIBLE_CURRENT = 1e-9;
        static constexpr double PLANAR_THICKNESS_RATIO_THRESHOLD = 10.0;
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
        CoilMesherModels select_mesh_model(Magnetic& magnetic);
        std::pair<size_t, size_t> calculate_grid_points(Magnetic& magnetic, double frequency);

};
} // namespace OpenMagnetics
