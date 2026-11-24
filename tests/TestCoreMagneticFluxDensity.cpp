#include "physical_models/CoreMagneticFluxDensity.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/CoreLosses.h"
#include "processors/Inputs.h"
#include "processors/MagneticSimulator.h"
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

TEST_CASE("Distributed core-loss setting prices the core losses from the flux density map", "[physical-model][core-flux-density][smoke-test]") {
    // Settings are global: put them back however the test ends.
    struct ResetSettings {
        ~ResetSettings() { Settings::GetInstance().reset(); }
    } resetSettings;
    Settings::GetInstance().reset();
    CHECK(Settings::GetInstance().get_core_losses_flux_density_distribution() == CoreLossesFluxDensityDistribution::EFFECTIVE_PARAMETERS);

    auto magnetic = magnetic_with_core("E 42/21/20", OpenMagneticsTesting::get_ground_gap(0.001));
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 4, 0.5, 5);
    auto operatingPoint = inputs.get_operating_point(0);
    MagnetizingInductance magnetizingInductanceModel;
    auto fluxDensity = magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(magnetic, &operatingPoint).second;
    auto excitation = operatingPoint.get_excitations_per_winding()[0];
    excitation.set_magnetic_flux_density(fluxDensity);
    auto core = magnetic.get_core();
    const double temperature = 60;

    CoreLosses coreLosses;
    auto effective = coreLosses.calculate_core_losses(core, excitation, temperature);
    CHECK(effective.get_method_used() != "DistributedFluxDensityMap");

    // DISTRIBUTED: the same entry point returns the loss map's total, built here on its own from the same B(t).
    Settings::GetInstance().set_core_losses_flux_density_distribution(CoreLossesFluxDensityDistribution::DISTRIBUTED);
    auto distributed = coreLosses.calculate_core_losses(core, excitation, temperature);
    CHECK(distributed.get_method_used() == "DistributedFluxDensityMap");
    CoreMagneticFluxDensity map(core, temperature, excitation.get_frequency());
    auto waveform = CoreMagneticFluxDensity::calculate_waveform(map, core, fluxDensity.get_waveform().value(), excitation.get_frequency());
    auto lossMap = CoreMagneticFluxDensity::calculate_core_losses(waveform);
    CHECK_THAT(distributed.get_core_losses(), Catch::Matchers::WithinRel(lossMap.coreLosses, 1e-12));
    CHECK_THAT(distributed.get_volumetric_losses().value(), Catch::Matchers::WithinRel(lossMap.coreLosses / lossMap.volume, 1e-12));
    // A second call reads the cached map: the same number.
    CHECK(coreLosses.calculate_core_losses(core, excitation, temperature).get_core_losses() == distributed.get_core_losses());
    // The map prices the same B(t) where it actually is: the same order as the effective-parameter losses, not equal.
    CHECK(distributed.get_core_losses() != effective.get_core_losses());
    CHECK_THAT(distributed.get_core_losses(), Catch::Matchers::WithinRel(effective.get_core_losses(), 0.5));

    // simulate (and its thermal loop) goes through it too.
    MagneticSimulator simulator;
    auto mas = simulator.simulate(inputs, magnetic);
    CHECK(mas.get_outputs()[0].get_core_losses().value().get_method_used() == "DistributedFluxDensityMap");

    Settings::GetInstance().reset();
    CHECK(Settings::GetInstance().get_core_losses_flux_density_distribution() == CoreLossesFluxDensityDistribution::EFFECTIVE_PARAMETERS);
}

namespace {

// A sinusoidal flux density excitation of amplitude peak, read over the core's reference area, one period sampled.
OperatingPointExcitation sinusoidal_flux_density_excitation(double peak, double frequency) {
    std::vector<double> time, data;
    for (size_t sample = 0; sample <= 64; ++sample) {
        double phase = static_cast<double>(sample) / 64;
        time.push_back(phase / frequency);
        data.push_back(peak * std::sin(2 * std::numbers::pi * phase));
    }
    Waveform waveform;
    waveform.set_time(time);
    waveform.set_data(data);
    SignalDescriptor fluxDensity;
    fluxDensity.set_waveform(waveform);
    fluxDensity.set_harmonics(OpenMagnetics::Inputs::calculate_harmonics_data(waveform, frequency));
    fluxDensity.set_processed(OpenMagnetics::Inputs::calculate_processed_data(fluxDensity, waveform, false));
    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    excitation.set_magnetic_flux_density(fluxDensity);
    return excitation;
}

Core processed_core(json functionalDescription) {
    json coreJson;
    coreJson["functionalDescription"] = functionalDescription;
    Core core(coreJson);
    core.process_data();
    core.process_gap();
    return core;
}

}  // namespace

TEST_CASE("Core flux density map covers every closed and drum family, and prices their losses", "[physical-model][core-flux-density][smoke-test]") {
    struct ResetSettings {
        ~ResetSettings() { Settings::GetInstance().reset(); }
    } resetSettings;
    Settings::GetInstance().reset();
    const double frequency = 100000;
    const double temperature = 25;
    auto excitation = sinusoidal_flux_density_excitation(0.05, frequency);

    std::vector<std::pair<std::string, Core>> cores;
    // Round and irregular posts (the plate depth calibrated to MKF's le / Ae), a piece and plate of each kind, a drum.
    for (std::string shape : {"ETD 29/16/10", "PQ 32/20", "RM 8", "P 18/11", "EP 13", "ET 20", "UT 20", "DRH-14X20-4C"}) {
        cores.emplace_back(shape, OpenMagneticsTesting::get_quick_core(shape, json::array(), 1, "3C90"));
    }
    // A shielded drum, its ring of the core's own grade.
    cores.emplace_back("DR 2.3 + SRI 3.0", OpenMagneticsTesting::get_quick_core("DR 2.3 + SRI 3.0", json::array(), 1, "3C90"));
    // A semi-shielded drum, its shell of magnetic epoxy (the core geometry test's design).
    json semishieldedShape = {
        {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "drumSemishielded"},
        {"aliases", json::array()}, {"name", "LQS-like 4018"},
        {"dimensions", {
            {"A", {{"nominal", 0.0038}}}, {"B", {{"nominal", 0.0018}}}, {"C", {{"nominal", 0.0015}}},
            {"D", {{"nominal", 0.0004}}}, {"E", {{"nominal", 0.0010}}}, {"F", {{"nominal", 0.0004}}},
            {"J", {{"nominal", 0.0040}}}, {"K", {{"nominal", 0.0040}}}, {"L", {{"nominal", 0.0018}}}}}};
    cores.emplace_back("semi-shielded drum", processed_core({
        {"type", "pieceAndPlate"}, {"material", "3C90"}, {"shape", semishieldedShape},
        {"gapping", json::array()}, {"numberStacks", 1},
        {"coating", {{"type", "magneticEpoxy"}, {"thickness", {{"nominal", 0.0001}}}, {"material", "Kool Mµ 26"}}}}));
    // A moulded body of two powders (the per-region inductance test's design).
    json moldedShape = {
        {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
        {"aliases", json::array()}, {"name", "MAPI-like 4020"},
        {"dimensions", {
            {"A", {{"nominal", 0.0041}}}, {"B", {{"nominal", 0.0021}}}, {"C", {{"nominal", 0.0041}}},
            {"D", {{"nominal", 0.0014}}}, {"E", {{"nominal", 0.0030}}}, {"F", {{"nominal", 0.0012}}}}}};
    cores.emplace_back("moulded body", processed_core({
        {"type", "closedShape"}, {"material", json::array({"Kool Mµ 60", "Kool Mµ 26"})}, {"shape", moldedShape},
        {"gapping", json::array()}, {"numberStacks", 1}}));

    for (auto& [name, core] : cores) {
        INFO(name);
        CoreMagneticFluxDensity map(core, temperature, frequency);
        REQUIRE(!map.get_cells().empty());
        const auto family = core.get_shape_family();
        const bool calibrated = family == CoreShapeFamily::ETD || family == CoreShapeFamily::PQ || family == CoreShapeFamily::RM ||
                                family == CoreShapeFamily::P || family == CoreShapeFamily::EP;
        if (calibrated) {
            // Round or irregular posts: the linear map has MKF's own reluctance, le / (mu0 mu Ae).
            auto instant = map.calculate_flux_densities(1e-12);
            double magnetomotiveForce = 0;
            for (double regionMagnetomotiveForce : instant.regionMagnetomotiveForces) {
                magnetomotiveForce += regionMagnetomotiveForce;
            }
            double mu = initial_relative_permeability(map) * constants.vacuumPermeability;
            CHECK_THAT(magnetomotiveForce * mu / 1e-12, Catch::Matchers::WithinRel(core.get_effective_length() / core.get_effective_area(), 1e-6));
        }
        if (name == "semi-shielded drum" || name == "moulded body") {
            CHECK(map.get_materials().size() == 2);
        }

        CoreLosses coreLosses;
        Settings::GetInstance().set_core_losses_flux_density_distribution(CoreLossesFluxDensityDistribution::DISTRIBUTED);
        auto distributed = coreLosses.calculate_core_losses(core, excitation, temperature);
        Settings::GetInstance().set_core_losses_flux_density_distribution(CoreLossesFluxDensityDistribution::EFFECTIVE_PARAMETERS);
        CHECK(distributed.get_method_used().starts_with("DistributedFluxDensityMap"));
        CHECK(std::isfinite(distributed.get_core_losses()));
        CHECK(distributed.get_core_losses() > 0);
    }
}
