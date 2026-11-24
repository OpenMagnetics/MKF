#include "physical_models/CoreMagneticFluxDensity.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/CoreLosses.h"
#include "processors/Inputs.h"
#include "support/Exceptions.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "json.hpp"
#include <cmath>
#include <numbers>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;
using json = nlohmann::json;

// ABT #2534: the flux density inside the core, cell by cell, from the core's shape.

namespace {

double polygon_area(const std::vector<std::array<double, 2>>& polygon) {
    double area = 0;
    for (size_t index = 0; index < polygon.size(); ++index) {
        const auto& current = polygon[index];
        const auto& next = polygon[(index + 1) % polygon.size()];
        area += current[0] * next[1] - next[0] * current[1];
    }
    return std::abs(area) / 2;
}

OpenMagnetics::Magnetic magnetic_with_core(const std::string& shape, json gapping, const std::string& material = "N87") {
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagneticsTesting::get_quick_core(shape, gapping, 1, material));
    magnetic.set_coil(OpenMagneticsTesting::get_quick_coil({20}, {1}, shape));
    return magnetic;
}

// mu_i of the curve the map reads: the slope of its first segment.
double initial_relative_permeability(const CoreMagneticFluxDensity& map) {
    auto instant = map.calculate_flux_densities(1e-12);
    return instant.relativePermeabilities[0];
}

}  // namespace

TEST_CASE("Core flux density map of an E core tiles its ferrite and matches the exact core reluctance", "[physical-model][core-flux-density][smoke-test]") {
    auto magnetic = magnetic_with_core("E 42/21/20", json::parse("[]"));
    auto core = magnetic.get_core();
    CoreMagneticFluxDensity map(magnetic, 25, 100000);

    // The cells cover the drawn half core (x from the main column's axis out), window excluded, exactly once.
    auto processed = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    double halfWidth = processed.get_width() / 2;
    double windowWidth = core.get_winding_windows()[0].get_width().value();
    double expectedArea = halfWidth * processed.get_height() - windowWidth * mainColumn.get_height();
    double cellArea = 0;
    for (const auto& cell : map.get_cells()) {
        cellArea += polygon_area(cell.polygon);
    }
    CHECK_THAT(cellArea, Catch::Matchers::WithinRel(expectedArea, 1e-3));

    // Low flux (linear): the loop's MMF over the centre-leg flux is the core reluctance, which MKF computes
    // independently from the same corners (exact C1 = le / Ae, ABT #1766).
    double flux = 1e-9;
    auto instant = map.calculate_flux_densities(flux);
    double loopMagnetomotiveForce = 0;
    for (double magnetomotiveForce : instant.regionMagnetomotiveForces) {
        loopMagnetomotiveForce += magnetomotiveForce;
    }
    double mu = initial_relative_permeability(map);
    double coreFactor = core.get_effective_length() / core.get_effective_area();
    double expectedReluctance = coreFactor / (constants.vacuumPermeability * mu);
    CHECK_THAT(loopMagnetomotiveForce / flux, Catch::Matchers::WithinRel(expectedReluctance, 0.01));

    // The straight parts carry their share uniformly: the half centre leg half the flux, on half its area.
    for (size_t cell = 0; cell < map.get_cells().size(); ++cell) {
        if (map.get_cells()[cell].region == "central column") {
            CHECK_THAT(instant.fluxDensities[cell], Catch::Matchers::WithinRel(flux / mainColumn.get_area(), 1e-9));
        }
    }
}

TEST_CASE("Core flux density map crowds at the window corner and the flux leaves it as it saturates", "[physical-model][core-flux-density][smoke-test]") {
    auto magnetic = magnetic_with_core("E 42/21/20", json::parse("[]"));
    auto core = magnetic.get_core();
    CoreMagneticFluxDensity map(magnetic, 25, 100000);
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});

    auto bend_extremes = [&](const CoreFluxDensityInstant& instant) {
        double highest = 0;
        double lowest = std::numeric_limits<double>::max();
        for (size_t cell = 0; cell < map.get_cells().size(); ++cell) {
            if (map.get_cells()[cell].region == "central column / top yoke corner") {
                highest = std::max(highest, std::abs(instant.fluxDensities[cell]));
                lowest = std::min(lowest, std::abs(instant.fluxDensities[cell]));
            }
        }
        return std::pair<double, double>{highest, lowest};
    };

    // Linear: the window corner is far above the leg's B, the core's outer corner far below. The field vanishes at the
    // outer corner point itself; its cell, a finite square, averages 0.30 of the leg's B at 12 rows (OMFEM 2D averaged
    // over the same polygon: 0.25).
    double smallFlux = 1e-9;
    double legFluxDensity = smallFlux / mainColumn.get_area();
    auto [linearHighest, linearLowest] = bend_extremes(map.calculate_flux_densities(smallFlux));
    CHECK(linearHighest > 2 * legFluxDensity);
    CHECK(linearLowest < 0.4 * legFluxDensity);

    // The leg at 90 % of B_sat: the corner saturates, its rows carry less, and its peak falls below the linear map's.
    double saturatingFlux = 0.9 * map.get_saturation_flux_density() * mainColumn.get_area();
    auto saturated = map.calculate_flux_densities(saturatingFlux);
    auto [saturatedHighest, saturatedLowest] = bend_extremes(saturated);
    double scale = saturatingFlux / smallFlux;
    CHECK(saturatedHighest < linearHighest * scale);
    CHECK(saturatedLowest > linearLowest * scale);

    // The map is odd in the flux.
    auto negative = map.calculate_flux_densities(-saturatingFlux);
    for (size_t cell = 0; cell < map.get_cells().size(); ++cell) {
        CHECK(negative.fluxDensities[cell] == -saturated.fluxDensities[cell]);
    }
}

TEST_CASE("Core flux density map of a gapped E core leaves the gap out", "[physical-model][core-flux-density][smoke-test]") {
    auto magnetic = magnetic_with_core("E 42/21/20", OpenMagneticsTesting::get_ground_gap(0.001));
    auto core = magnetic.get_core();
    CoreMagneticFluxDensity map(magnetic, 25, 100000);
    auto processed = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    double halfWidth = processed.get_width() / 2;
    double windowWidth = core.get_winding_windows()[0].get_width().value();
    double gapArea = 0;
    for (const auto& gap : core.get_functional_description().get_gapping()) {
        auto column = core.find_closest_column_by_coordinates(gap.get_coordinates().value());
        double drawnWidth = column.get_coordinates()[0] == 0 ? column.get_width() / 2 : column.get_width();
        if (column.get_coordinates()[0] < 0) {
            continue;  // the left lateral leg is not in the drawn half
        }
        gapArea += drawnWidth * gap.get_length();
    }
    double expectedArea = halfWidth * processed.get_height() - windowWidth * mainColumn.get_height() - gapArea;
    double cellArea = 0;
    for (const auto& cell : map.get_cells()) {
        cellArea += polygon_area(cell.polygon);
    }
    CHECK_THAT(cellArea, Catch::Matchers::WithinRel(expectedArea, 1e-3));
    CHECK_THROWS_AS(map.calculate_flux_density_at_point(0.0001, 0, 1e-6), InvalidInputException);
}

TEST_CASE("Core flux density map of a toroid falls as 1/r", "[physical-model][core-flux-density][smoke-test]") {
    auto magnetic = magnetic_with_core("T 40/24/16", json::parse("[]"));
    auto core = magnetic.get_core();
    CoreMagneticFluxDensity map(magnetic, 25, 100000);
    auto processed = core.get_processed_description().value();
    double outerRadius = processed.get_width() / 2;
    double innerRadius = outerRadius - core.find_closest_column_by_coordinates({0, 0, 0}).get_width();
    double height = core.find_closest_column_by_coordinates({0, 0, 0}).get_depth();

    double flux = 1e-9;
    double mu = initial_relative_permeability(map);
    auto instant = map.calculate_flux_densities(flux);
    // Exact ring reluctance: 2 pi / (mu0 mu h ln(ro / ri)).
    double expectedReluctance = 2 * std::numbers::pi / (constants.vacuumPermeability * mu * height * std::log(outerRadius / innerRadius));
    CHECK_THAT(instant.regionMagnetomotiveForces[0] / flux, Catch::Matchers::WithinRel(expectedReluctance, 1e-3));
    double inner = map.calculate_flux_density_at_point(innerRadius * 1.01, 0, flux);
    double outer = map.calculate_flux_density_at_point(0, outerRadius * 0.99, flux);
    CHECK_THAT(inner / outer, Catch::Matchers::WithinRel(outerRadius / innerRadius, 0.1));
}

TEST_CASE("Core flux density map over an operating point peaks where MKF's flux peaks", "[physical-model][core-flux-density][smoke-test]") {
    auto magnetic = magnetic_with_core("E 42/21/20", OpenMagneticsTesting::get_ground_gap(0.001));
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 4, 0.5, 5);
    auto operatingPoint = inputs.get_operating_point(0);
    auto result = CoreMagneticFluxDensity::calculate_waveform(magnetic, operatingPoint);
    REQUIRE(!result.fluxDensities.empty());
    double peakFlux = std::abs(result.magneticFlux[result.peakSampleIndex]);
    for (double flux : result.magneticFlux) {
        CHECK(std::abs(flux) <= peakFlux);
    }
    // The centre leg at the peak carries the wound column's flux on its own area.
    auto core = magnetic.get_core();
    double centreArea = core.find_closest_column_by_coordinates({0, 0, 0}).get_area();
    CoreMagneticFluxDensity map(magnetic, 25, 100000);
    for (size_t cell = 0; cell < map.get_cells().size(); ++cell) {
        if (map.get_cells()[cell].region == "central column") {
            CHECK_THAT(result.peakFluxDensities[cell], Catch::Matchers::WithinRel(peakFlux / centreArea, 1e-9));
        }
        CHECK(result.alternatingFluxDensityAmplitudes[cell] <= result.peakFluxDensities[cell] * (1 + 1e-9));
    }
}

TEST_CASE("Core loss map prices each cell's own flux density with MKF's core-loss model", "[physical-model][core-flux-density][smoke-test]") {
    auto magnetic = magnetic_with_core("E 42/21/20", OpenMagneticsTesting::get_ground_gap(0.001));
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 4, 0.5, 5);
    auto waveform = CoreMagneticFluxDensity::calculate_waveform(magnetic, inputs.get_operating_point(0));
    auto lossMap = CoreMagneticFluxDensity::calculate_core_losses(magnetic, waveform);
    auto core = magnetic.get_core();

    // The cells stand for the whole core: the drawn half's area, gaps left out, times the depth, twice.
    auto processed = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    double drawnArea = 0;
    for (const auto& cell : waveform.cells) {
        drawnArea += polygon_area(cell.polygon);
    }
    CHECK_THAT(lossMap.volume, Catch::Matchers::WithinRel(2 * drawnArea * mainColumn.get_depth(), 1e-3));

    // A centre-leg cell carries B = Phi / A_centre all through the period: its losses are the kernel's at that B(t).
    double centreArea = mainColumn.get_area();
    Waveform centreWaveform;
    centreWaveform.set_time(waveform.time);
    std::vector<double> centreData;
    for (double flux : waveform.magneticFlux) {
        centreData.push_back(flux / centreArea);
    }
    centreWaveform.set_data(centreData);
    SignalDescriptor fluxDensity;
    fluxDensity.set_waveform(centreWaveform);
    fluxDensity.set_harmonics(OpenMagnetics::Inputs::calculate_harmonics_data(centreWaveform, waveform.frequency));
    fluxDensity.set_processed(OpenMagnetics::Inputs::calculate_processed_data(fluxDensity, centreWaveform, false));
    OperatingPointExcitation excitation;
    excitation.set_frequency(waveform.frequency);
    excitation.set_magnetic_flux_density(fluxDensity);
    CoreLosses coreLosses;
    double expected = coreLosses.get_core_volumetric_losses(core.resolve_material(), excitation, waveform.temperature);
    size_t centreCells = 0;
    double sum = 0;
    for (size_t cell = 0; cell < waveform.cells.size(); ++cell) {
        if (waveform.cells[cell].region == "central column") {
            ++centreCells;
            CHECK_THAT(lossMap.volumetricLosses[cell], Catch::Matchers::WithinRel(expected, 1e-6));
        }
        sum += lossMap.volumetricLosses[cell] * waveform.cells[cell].volume;
    }
    CHECK(centreCells > 0);
    CHECK_THAT(lossMap.coreLosses, Catch::Matchers::WithinRel(sum, 1e-12));
}

TEST_CASE("Core loss map of a toroid costs more than its mean flux density uniformly spread", "[physical-model][core-flux-density][smoke-test]") {
    auto magnetic = magnetic_with_core("T 40/24/16", json::parse("[]"));
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 0.1, 0.5, 0.05);
    auto waveform = CoreMagneticFluxDensity::calculate_waveform(magnetic, inputs.get_operating_point(0));
    auto lossMap = CoreMagneticFluxDensity::calculate_core_losses(magnetic, waveform);
    auto core = magnetic.get_core();
    auto processed = core.get_processed_description().value();
    double outerRadius = processed.get_width() / 2;
    auto column = core.find_closest_column_by_coordinates({0, 0, 0});
    double innerRadius = outerRadius - column.get_width();
    double height = column.get_depth();

    // The rings fill the ring exactly.
    CHECK_THAT(lossMap.volume, Catch::Matchers::WithinRel(std::numbers::pi * (outerRadius * outerRadius - innerRadius * innerRadius) * height, 1e-9));

    // Steinmetz losses grow faster than B (beta > 1): B ~ 1/r costs more than its volume-weighted mean B everywhere
    // (Jensen). A ring's volume grows with r, so that mean is B at the mean radius, (ro + ri) / 2, not the section's mean.
    double meanRadius = (outerRadius + innerRadius) / 2;
    Waveform meanWaveform;
    meanWaveform.set_time(waveform.time);
    std::vector<double> meanData;
    for (double flux : waveform.magneticFlux) {
        meanData.push_back(flux / (height * meanRadius * std::log(outerRadius / innerRadius)));
    }
    meanWaveform.set_data(meanData);
    SignalDescriptor fluxDensity;
    fluxDensity.set_waveform(meanWaveform);
    fluxDensity.set_harmonics(OpenMagnetics::Inputs::calculate_harmonics_data(meanWaveform, waveform.frequency));
    fluxDensity.set_processed(OpenMagnetics::Inputs::calculate_processed_data(fluxDensity, meanWaveform, false));
    OperatingPointExcitation excitation;
    excitation.set_frequency(waveform.frequency);
    excitation.set_magnetic_flux_density(fluxDensity);
    CoreLosses coreLosses;
    double uniform = coreLosses.get_core_volumetric_losses(core.resolve_material(), excitation, waveform.temperature) * lossMap.volume;
    CHECK(lossMap.coreLosses > uniform);
    // The density falls ring by ring from the bore out (the rings are built from the inner radius outwards).
    for (size_t ring = 1; ring < lossMap.volumetricLosses.size(); ++ring) {
        CHECK(lossMap.volumetricLosses[ring] < lossMap.volumetricLosses[ring - 1]);
    }
}
