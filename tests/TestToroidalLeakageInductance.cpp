#include <source_location>
#include <fstream>
#include <numbers>
#include "physical_models/LeakageInductance.h"
#include "physical_models/Impedance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// Toroidal leakage model (LeakageInductance.h, "Toroidal cores"). The ring-plane field of a toroid lives in
// two air regions, the bore and the outside, which a high-permeability core decouples. Until 2026-09 MKF
// imaged every conductor in BOTH walls (an outer conductor's image landed in the bore and cancelled the bore
// field), dropped every grid point more than two wire diameters from a conductor (the bore centre and the
// gap between the windings, where a common-mode choke's differential field lives), stopped the grid at the
// outermost conductor, and extruded a z-invariant 2-D field. On a WE common-mode choke that gave 0.9 uH
// against 10 uH measured.

namespace {

std::vector<LeakageInductance::ToroidalLineCurrent> sector_conductors(size_t turnsPerWinding, double radius, double wireRadius, double firstSectorCentre, double sectorSpan, double sign) {
    // Two windings, +1 A and -1 A, each spread uniformly over sectorSpan, centred at firstSectorCentre and firstSectorCentre + pi.
    std::vector<LeakageInductance::ToroidalLineCurrent> conductors;
    for (size_t winding = 0; winding < 2; ++winding) {
        double centre = firstSectorCentre + winding * std::numbers::pi;
        double current = sign * (winding == 0 ? 1.0 : -1.0);
        for (size_t k = 0; k < turnsPerWinding; ++k) {
            double angle = centre - sectorSpan / 2 + sectorSpan * (k + 0.5) / turnsPerWinding;
            conductors.push_back({radius * std::cos(angle), radius * std::sin(angle), current, std::exp(-0.25) * wireRadius});
        }
    }
    return conductors;
}

// ∫ ½µ0|H|² dA over an annulus [innerRadius, outerRadius] on a polar midpoint grid (geometric radial spacing when
// the annulus is wide), from the image-set field.
double integrate_field_energy(const std::vector<LeakageInductance::ToroidalLineCurrent>& conductors, double wallRadius, bool interiorRegion, double innerRadius, double outerRadius,
                              size_t radialPoints, size_t angularPoints, bool geometric) {
    double vacuumPermeability = 4e-7 * std::numbers::pi;
    double energy = 0;
    double dphi = 2 * std::numbers::pi / angularPoints;
    for (size_t i = 0; i < radialPoints; ++i) {
        double r0, r1;
        if (geometric) {
            r0 = innerRadius * std::pow(outerRadius / innerRadius, double(i) / radialPoints);
            r1 = innerRadius * std::pow(outerRadius / innerRadius, double(i + 1) / radialPoints);
        }
        else {
            r0 = innerRadius + (outerRadius - innerRadius) * double(i) / radialPoints;
            r1 = innerRadius + (outerRadius - innerRadius) * double(i + 1) / radialPoints;
        }
        double r = (r0 + r1) / 2;
        double area = r * (r1 - r0) * dphi;
        for (size_t j = 0; j < angularPoints; ++j) {
            double phi = (j + 0.5) * dphi;
            auto h = LeakageInductance::calculate_ring_plane_region_field(conductors, wallRadius, 1.0, interiorRegion, r * std::cos(phi), r * std::sin(phi));
            energy += 0.5 * vacuumPermeability * (h[0] * h[0] + h[1] * h[1]) * area;
        }
    }
    return energy;
}

} // namespace

TEST_CASE("Toroidal leakage: per-region images leave no tangential field on the core walls", "[physical-model][leakage-inductance][toroidal][toroidal-leakage]") {
    double innerWall = 4e-3;
    double outerWall = 7e-3;
    auto bore = sector_conductors(5, 3.4e-3, 0.25e-3, 0.3, 1.9, 1.0);
    auto outside = sector_conductors(5, 7.6e-3, 0.25e-3, 0.3, 1.9, -1.0);
    double fieldScale = 1.0 / (2 * std::numbers::pi * innerWall);
    for (size_t k = 0; k < 90; ++k) {
        double phi = 2 * std::numbers::pi * k / 90;
        auto hBore = LeakageInductance::calculate_ring_plane_region_field(bore, innerWall, 1.0, true, innerWall * std::cos(phi), innerWall * std::sin(phi));
        double tangentialBore = -hBore[0] * std::sin(phi) + hBore[1] * std::cos(phi);
        CHECK(std::abs(tangentialBore) < 1e-9 * fieldScale);
        auto hOutside = LeakageInductance::calculate_ring_plane_region_field(outside, outerWall, 1.0, false, outerWall * std::cos(phi), outerWall * std::sin(phi));
        double tangentialOutside = -hOutside[0] * std::sin(phi) + hOutside[1] * std::cos(phi);
        CHECK(std::abs(tangentialOutside) < 1e-9 * fieldScale);
    }
}

TEST_CASE("Toroidal leakage: closed-form ring-plane energy equals the integrated field energy", "[physical-model][leakage-inductance][toroidal][toroidal-leakage]") {
    // A sectored winding pair (the common-mode-choke topology). The closed-form energy (with the in-wire self term)
    // must equal ½µ0∫|H|² over the WHOLE bore and the WHOLE outside. The grid levels show the integral converging;
    // the truncated outside and the near-conductor band show what the old toroid grid missed.
    double innerWall = 4e-3;
    double outerWall = 7e-3;
    double wireRadius = 0.25e-3;
    auto bore = sector_conductors(12, innerWall - 0.6e-3, wireRadius, std::numbers::pi / 2, 2.2, 1.0);
    auto outside = sector_conductors(12, outerWall + 0.6e-3, wireRadius, std::numbers::pi / 2, 2.2, -1.0);
    double closedBore = LeakageInductance::calculate_ring_plane_region_energy_per_length(bore, innerWall, 1.0, true);
    double closedOutside = LeakageInductance::calculate_ring_plane_region_energy_per_length(outside, outerWall, 1.0, false);

    std::vector<double> errors;
    for (size_t level : {1, 2, 4}) {
        double gridBore = integrate_field_energy(bore, innerWall, true, 0, innerWall, 100 * level, 360 * level, false);
        double gridOutside = integrate_field_energy(outside, outerWall, false, outerWall, 30 * outerWall, 150 * level, 360 * level, true);
        errors.push_back(std::abs((gridBore + gridOutside) / (closedBore + closedOutside) - 1));
    }
    CHECK(errors[1] < errors[0]);
    CHECK(errors[2] < errors[1]);
    CHECK(errors[2] < 0.005);

    // Outside region: stopping at 30·A/2 loses < 0.5 %; stopping two wire diameters past the conductors (the old
    // grid edge) loses a large share of the outside energy.
    double farTruncated = integrate_field_energy(outside, outerWall, false, outerWall, 30 * outerWall, 600, 1440, true);
    double farer = integrate_field_energy(outside, outerWall, false, outerWall, 60 * outerWall, 660, 1440, true);
    CHECK_THAT(farTruncated, WithinRel(farer, 0.005));
    double nearOnly = integrate_field_energy(outside, outerWall, false, outerWall, outerWall + 0.6e-3 + 4 * wireRadius, 200, 1440, false);
    CHECK(nearOnly < 0.8 * closedOutside);
}

TEST_CASE("Toroidal leakage: a net current keeps the core's circulation at the wall and a unit-free energy", "[physical-model][leakage-inductance][toroidal][toroidal-leakage]") {
    // One sectored winding alone (the self-leakage of the inductance matrix, or any excitation whose ampere-turns do
    // not balance): the bore carries +12 A net and the outside -12 A. That current links the core, so the wall carries
    // the core's magnetizing H = 12 A / (2 pi R), and the energy must be the field energy of that configuration. The
    // centre-image kernel alone (right for a balanced pair) makes the wall field-free and gives an energy that
    // changes sign when the geometry is written in millimetres instead of metres.
    double innerWall = 4e-3;
    double outerWall = 7e-3;
    double wireRadius = 0.25e-3;
    auto borePair = sector_conductors(12, innerWall - 0.6e-3, wireRadius, std::numbers::pi / 2, 2.2, 1.0);
    auto outsidePair = sector_conductors(12, outerWall + 0.6e-3, wireRadius, std::numbers::pi / 2, 2.2, -1.0);
    std::vector<LeakageInductance::ToroidalLineCurrent> bore(borePair.begin(), borePair.begin() + 12);
    std::vector<LeakageInductance::ToroidalLineCurrent> outside(outsidePair.begin(), outsidePair.begin() + 12);

    // Tangential field on each wall: the uniform magnetizing circulation of the 12 A linking the core.
    for (size_t k = 0; k < 90; ++k) {
        double phi = 2 * std::numbers::pi * k / 90;
        auto hBore = LeakageInductance::calculate_ring_plane_region_field(bore, innerWall, 1.0, true, innerWall * std::cos(phi), innerWall * std::sin(phi));
        CHECK_THAT(-hBore[0] * std::sin(phi) + hBore[1] * std::cos(phi), WithinRel(12.0 / (2 * std::numbers::pi * innerWall), 1e-9));
        auto hOutside = LeakageInductance::calculate_ring_plane_region_field(outside, outerWall, 1.0, false, outerWall * std::cos(phi), outerWall * std::sin(phi));
        CHECK_THAT(-hOutside[0] * std::sin(phi) + hOutside[1] * std::cos(phi), WithinRel(12.0 / (2 * std::numbers::pi * outerWall), 1e-9));
    }

    // Closed form = integrated field energy, each region.
    double closedBore = LeakageInductance::calculate_ring_plane_region_energy_per_length(bore, innerWall, 1.0, true);
    double closedOutside = LeakageInductance::calculate_ring_plane_region_energy_per_length(outside, outerWall, 1.0, false);
    CHECK_THAT(integrate_field_energy(bore, innerWall, true, 0, innerWall, 400, 1440, false), WithinRel(closedBore, 0.005));
    CHECK_THAT(integrate_field_energy(outside, outerWall, false, outerWall, 30 * outerWall, 600, 1440, true), WithinRel(closedOutside, 0.005));

    // 2-D energy per length is scale free: the same winding 1000x larger or smaller has the same energy per length,
    // for the realistic image factor too.
    auto scaled = [](std::vector<LeakageInductance::ToroidalLineCurrent> conductors, double factor) {
        for (auto& conductor : conductors) {
            conductor.x *= factor;
            conductor.y *= factor;
            conductor.geometricMeanRadius *= factor;
        }
        return conductors;
    };
    for (double imageFactor : {1.0, 0.998}) {
        double reference = LeakageInductance::calculate_ring_plane_region_energy_per_length(bore, innerWall, imageFactor, true) +
                           LeakageInductance::calculate_ring_plane_region_energy_per_length(outside, outerWall, imageFactor, false);
        CHECK(reference > 0);
        for (double factor : {1e-3, 1e3}) {
            double energy = LeakageInductance::calculate_ring_plane_region_energy_per_length(scaled(bore, factor), innerWall * factor, imageFactor, true) +
                            LeakageInductance::calculate_ring_plane_region_energy_per_length(scaled(outside, factor), outerWall * factor, imageFactor, false);
            CHECK_THAT(energy, WithinRel(reference, 1e-9));
        }
    }

    // A uniform ring of N conductors carrying I in total: the annulus energy mu0 I^2 / (4 pi) ln(R / rho) in the bore and
    // mu0 I^2 / (4 pi) ln(rho / R) outside, relative to the same ring against a wall 1.5x farther away.
    auto ring = [](double radius, double current) {
        std::vector<LeakageInductance::ToroidalLineCurrent> conductors;
        for (size_t k = 0; k < 400; ++k) {
            double angle = 2 * std::numbers::pi * (k + 0.5) / 400;
            conductors.push_back({radius * std::cos(angle), radius * std::sin(angle), current / 400, 1e-9});
        }
        return conductors;
    };
    double annulus = 4e-7 * std::numbers::pi / (4 * std::numbers::pi) * std::log(1.5);
    CHECK_THAT(LeakageInductance::calculate_ring_plane_region_energy_per_length(ring(3e-3, 1.0), 1.5 * innerWall, 1.0, true) -
               LeakageInductance::calculate_ring_plane_region_energy_per_length(ring(3e-3, 1.0), innerWall, 1.0, true), WithinRel(annulus, 1e-6));
    CHECK_THAT(LeakageInductance::calculate_ring_plane_region_energy_per_length(ring(8e-3, -1.0), outerWall / 1.5, 1.0, false) -
               LeakageInductance::calculate_ring_plane_region_energy_per_length(ring(8e-3, -1.0), outerWall, 1.0, false), WithinRel(annulus, 1e-6));
}

TEST_CASE("Toroidal leakage: sectored-winding ring-plane limit L' = 28 zeta(3) mu0 N^2 / pi^3", "[physical-model][leakage-inductance][toroidal][toroidal-leakage]") {
    // Two windings of N turns, each spread over a half circle, hugging a high-permeability wall. In the continuous
    // limit the MMF V(φ) is a triangle wave; its harmonics give L' = (2µ0/π)·Σ_m |S_m|²/m per unit length, which for
    // half-circle sectors sums to 28·ζ(3)·µ0·N²/π³ (both the bore and the outside counted). The discrete sum of
    // line currents with the per-region images converges to it as N grows.
    double vacuumPermeability = 4e-7 * std::numbers::pi;
    double zeta3 = 1.2020569031595942;
    double expected = 28 * zeta3 * vacuumPermeability / std::pow(std::numbers::pi, 3);
    double innerWall = 4e-3;
    double outerWall = 7e-3;
    size_t numberTurns = 800;
    double wireRadius = std::numbers::pi * innerWall / numberTurns / 4;
    auto bore = sector_conductors(numberTurns, innerWall - wireRadius, wireRadius, std::numbers::pi / 2, std::numbers::pi, 1.0);
    auto outside = sector_conductors(numberTurns, outerWall + wireRadius, wireRadius, std::numbers::pi / 2, std::numbers::pi, -1.0);
    double energyPerLength = LeakageInductance::calculate_ring_plane_region_energy_per_length(bore, innerWall, 1.0, true) +
                             LeakageInductance::calculate_ring_plane_region_energy_per_length(outside, outerWall, 1.0, false);
    double numberTurnsSquared = double(numberTurns) * double(numberTurns);
    CHECK_THAT(2 * energyPerLength / numberTurnsSquared, WithinRel(expected, 0.005));

    // The MMF-sheet harmonics used by the 3-D correction carry the same per-length energy.
    std::vector<double> angles;
    std::vector<double> currents;
    for (auto& conductor : bore) {
        angles.push_back(std::atan2(conductor.y, conductor.x));
        currents.push_back(conductor.current);
    }
    double sheetEnergy = 0;
    for (size_t mode = 1; mode <= 1000; ++mode) {
        sheetEnergy += LeakageInductance::calculate_sheet_mode_energy_per_length(angles, currents, mode);
    }
    CHECK_THAT(2 * sheetEnergy / numberTurnsSquared, WithinRel(expected, 0.005));
}

TEST_CASE("Toroidal leakage: body-of-revolution effective height converges and has the right limits", "[physical-model][leakage-inductance][toroidal][toroidal-leakage]") {
    double innerRadius = 3.13e-3;
    double outerRadius = 7.87e-3;
    double height = 10.7e-3;
    double halfTurnLength = height + outerRadius - innerRadius;

    // Mesh and far-boundary convergence.
    for (size_t mode : {1, 4}) {
        double coarse = LeakageInductance::calculate_body_of_revolution_effective_height(innerRadius, outerRadius, height, mode, 1.0, 30.0);
        double fine = LeakageInductance::calculate_body_of_revolution_effective_height(innerRadius, outerRadius, height, mode, 2.0, 30.0);
        double far = LeakageInductance::calculate_body_of_revolution_effective_height(innerRadius, outerRadius, height, mode, 1.0, 60.0);
        CHECK_THAT(coarse, WithinRel(fine, 0.003));
        CHECK_THAT(coarse, WithinRel(far, 0.001));
    }

    // Low harmonics reach beyond the envelope (h_eff > ℓ∞); high harmonics hug it (h_eff → ℓ∞ = h + w).
    double previous = std::numeric_limits<double>::max();
    for (size_t mode : {1, 2, 4, 8, 16, 32}) {
        double effectiveHeight = LeakageInductance::calculate_body_of_revolution_effective_height(innerRadius, outerRadius, height, mode);
        CHECK(effectiveHeight < previous);
        CHECK(effectiveHeight > halfTurnLength);
        previous = effectiveHeight;
    }
    CHECK_THAT(previous, WithinRel(halfTurnLength, 0.03));

    // Tall limit: every added unit of height adds exactly the 2-D per-length energy (h_eff grows by 1 per unit height).
    for (size_t mode : {1, 3}) {
        double shortCylinder = LeakageInductance::calculate_body_of_revolution_effective_height(innerRadius, outerRadius, 40e-3, mode);
        double tallCylinder = LeakageInductance::calculate_body_of_revolution_effective_height(innerRadius, outerRadius, 80e-3, mode);
        CHECK_THAT((tallCylinder - shortCylinder) / 40e-3, WithinRel(1.0, 0.01));
    }
}

TEST_CASE("Toroidal leakage of a WE common-mode choke against its measured differential-mode inductance", "[physical-model][leakage-inductance][toroidal][toroidal-leakage][cmc]") {
    // WE 744822222 (WE-CMB, T 14/8/9 MnZn ring, 0.6 mm epoxy coating, 2 x 18 turns of 0.5 mm wire). Measured
    // differential-mode (series-opposing loop) inductance: 9.9-10.7 uH from the WE .s4p between 100 kHz and 10 MHz
    // (100 ohm differential reference), 10.9 uH in the CISPR 17 release test at 10 kHz. The fixture carries the
    // functional description only, so MKF winds each winding over most of its half of the ring, as the part is.
    // Tolerance 20 %: the measurements spread by 10 %, and the model moves by about 10 % for every 20 degrees of
    // winding sector, which the datasheet does not give. The previous ring-plane model gave 0.85 uH (-92 %) and a
    // 2-D model without the 3-D correction gives 6 uH (-40 %).
    settings.reset();
    auto testDataPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_we_744822222_functional.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    OpenMagnetics::Magnetic magnetic(nlohmann::json::parse(file));
    magnetic = magnetic_autocomplete(magnetic);

    double measuredDifferentialModeInductance = 10.0e-6;
    auto leakage = LeakageInductance().calculate_leakage_inductance(magnetic, 100e3, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    CHECK_THAT(leakage, WithinRel(measuredDifferentialModeInductance, 0.2));

    // The 3-D correction is a large share of it: the ring-plane part alone stays well below the measurement.
    auto parts = LeakageInductance::calculate_toroidal_leakage_energy(magnetic, {1.0, -1.0});
    CHECK(2 * parts.ringPlaneEnergyPerLength * parts.extrusionLength < 0.75 * measuredDifferentialModeInductance);
    CHECK_THAT(2 * parts.energy, WithinRel(leakage, 1e-12));
}

TEST_CASE("Toroidal leakage matrix reproduces the pairwise leakage", "[physical-model][leakage-inductance][toroidal][toroidal-leakage]") {
    settings.reset();
    auto testDataPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_we_744822222_functional.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    OpenMagnetics::Magnetic magnetic(nlohmann::json::parse(file));
    magnetic = magnetic_autocomplete(magnetic);
    auto pairwise = LeakageInductance().calculate_leakage_inductance(magnetic, 100e3, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    auto matrix = LeakageInductance().calculate_leakage_inductance_matrix(magnetic, 100e3);
    // Equal turns: the balanced pair (1, -1) gives Λ00 + Λ11 − 2Λ01.
    CHECK_THAT(matrix[0][0] + matrix[1][1] - 2 * matrix[0][1], WithinRel(pairwise, 1e-6));
    CHECK_THAT(matrix[0][1], WithinRel(matrix[1][0], 1e-12));
}

TEST_CASE("Differential-mode impedance uses the loop resistance of both windings", "[physical-model][impedance][cmc][toroidal-leakage]") {
    // The differential-mode current flows through the two windings in series; at low frequency the DM impedance is
    // the loop resistance R1 + R2 (it used winding 0 only, half the loop for a symmetric choke).
    settings.reset();
    auto testDataPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_we_744822222_functional.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    OpenMagnetics::Magnetic magnetic(nlohmann::json::parse(file));
    magnetic = magnetic_autocomplete(magnetic);
    double temperature = 25;
    auto resistances = WindingOhmicLosses::calculate_dc_resistance_per_winding(magnetic.get_coil(), temperature);
    auto parameters = Impedance().calculate_differential_mode_parameters(magnetic.get_core(), magnetic.get_coil(), 100e3, temperature);
    CHECK_THAT(parameters.windingResistance, WithinRel(resistances[0] + resistances[1], 1e-12));
    auto lowFrequency = Impedance().differential_mode_impedance_from_parameters(parameters, 10.0);
    CHECK_THAT(lowFrequency.real(), WithinRel(resistances[0] + resistances[1], 1e-3));
}
