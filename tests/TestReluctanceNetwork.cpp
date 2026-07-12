#include "physical_models/ReluctanceNetwork.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/Reluctance.h"
#include "physical_models/Inductance.h"
#include "constructive_models/Core.h"
#include "constructive_models/Bobbin.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;
using json = nlohmann::json;

// Phase 3 of multi-column winding support: the per-column reluctance network with one
// MMF source per wound column, generalizing the lumped N^2/R magnetizing model.

namespace {

OpenMagnetics::Magnetic buildTwoWindingMagnetic(std::optional<int64_t> secondaryWindow) {
    auto& settings = Settings::GetInstance();
    settings.set_core_per_column_winding_windows(true);
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/20", json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, 0.001, 0.001);
    json coilJson;
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    coilJson["functionalDescription"].push_back(json{{"name", "Primary"}, {"numberTurns", 20}, {"numberParallels", 1},
                                                     {"isolationSide", "primary"}, {"wire", "Round 0.475 - Grade 1"}});
    coilJson["functionalDescription"].push_back(json{{"name", "Secondary"}, {"numberTurns", 10}, {"numberParallels", 1},
                                                     {"isolationSide", "secondary"}, {"wire", "Round 0.475 - Grade 1"}});
    OpenMagnetics::Coil coil(coilJson, false);
    if (secondaryWindow) {
        coil.get_mutable_functional_description()[1].set_winding_window(secondaryWindow.value());
    }
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    return magnetic;
}

ReluctanceNetwork buildCircuit(Core core) {
    auto reluctanceModel = ReluctanceModel::factory();
    auto reluctanceOutput = reluctanceModel->get_core_reluctance(core, std::optional<OperatingPoint>(std::nullopt));
    return ReluctanceNetwork(core, reluctanceOutput.get_ungapped_core_reluctance().value(),
                           reluctanceOutput.get_reluctance_per_gap().value_or(std::vector<AirGapReluctanceOutput>{}));
}

}  // namespace

TEST_CASE("ReluctanceNetwork_MainPlacement_ReproducesLumpedModel", "[physical-model][magnetic-circuit][multi-column][smoke-test]") {
    auto magnetic = buildTwoWindingMagnetic(std::nullopt);
    auto core = magnetic.get_core();
    CHECK(!ReluctanceNetwork::has_non_main_placement(magnetic));

    auto reluctanceModel = ReluctanceModel::factory();
    auto reluctanceOutput = reluctanceModel->get_core_reluctance(core, std::optional<OperatingPoint>(std::nullopt));
    double lumpedReluctance = reluctanceOutput.get_core_reluctance();

    auto circuit = buildCircuit(core);
    auto inductanceMatrix = circuit.calculate_magnetizing_inductance_matrix(magnetic);

    // With both windings on the main column, the network is exactly the lumped model.
    CHECK_THAT(inductanceMatrix[0][0], Catch::Matchers::WithinRel(20.0 * 20.0 / lumpedReluctance, 1e-9));
    CHECK_THAT(inductanceMatrix[1][1], Catch::Matchers::WithinRel(10.0 * 10.0 / lumpedReluctance, 1e-9));
    CHECK_THAT(inductanceMatrix[0][1], Catch::Matchers::WithinRel(20.0 * 10.0 / lumpedReluctance, 1e-9));
    CHECK(inductanceMatrix[0][1] == inductanceMatrix[1][0]);
}

TEST_CASE("ReluctanceNetwork_LateralPlacement_FluxDividerCoupling", "[physical-model][magnetic-circuit][multi-column][smoke-test]") {
    auto magnetic = buildTwoWindingMagnetic(2);
    auto core = magnetic.get_core();
    REQUIRE(ReluctanceNetwork::has_non_main_placement(magnetic));

    auto circuit = buildCircuit(core);
    auto columnIndexes = ReluctanceNetwork::resolve_winding_column_indexes(magnetic);
    size_t mainColumnIndex = circuit.get_main_column_index();
    REQUIRE(columnIndexes[0] == mainColumnIndex);
    REQUIRE(columnIndexes[1] != mainColumnIndex);

    auto columnReluctances = circuit.get_column_reluctances();
    REQUIRE(columnReluctances.size() == 3);

    // Closed-form driving-point inductances from the branch reluctances.
    auto parallelOf = [](double a, double b) { return 1 / (1 / a + 1 / b); };
    size_t secondaryColumn = columnIndexes[1];
    size_t otherLateralColumn = 3 - mainColumnIndex - secondaryColumn;
    double primaryDrivingReluctance =
        columnReluctances[mainColumnIndex] + parallelOf(columnReluctances[secondaryColumn], columnReluctances[otherLateralColumn]);
    double secondaryDrivingReluctance =
        columnReluctances[secondaryColumn] + parallelOf(columnReluctances[mainColumnIndex], columnReluctances[otherLateralColumn]);

    auto inductanceMatrix = circuit.calculate_magnetizing_inductance_matrix(magnetic);
    CHECK_THAT(inductanceMatrix[0][0], Catch::Matchers::WithinRel(20.0 * 20.0 / primaryDrivingReluctance, 1e-9));
    CHECK_THAT(inductanceMatrix[1][1], Catch::Matchers::WithinRel(10.0 * 10.0 / secondaryDrivingReluctance, 1e-9));

    // Coupling between leg-separated windings is imperfect and negative in the common
    // branch orientation (flux up one leg comes down the other).
    double couplingCoefficient = inductanceMatrix[0][1] / std::sqrt(inductanceMatrix[0][0] * inductanceMatrix[1][1]);
    CHECK(couplingCoefficient < 0);
    CHECK(std::abs(couplingCoefficient) < 1.0);
    CHECK(std::abs(couplingCoefficient) > 0.1);

    // Flux conservation: the column fluxes sum to zero, and the driven column carries
    // the largest magnitude.
    auto columnFluxes = circuit.calculate_column_magnetic_fluxes(magnetic, {1.0, 0.0});
    double fluxSum = 0;
    for (auto flux : columnFluxes) {
        fluxSum += flux;
    }
    CHECK_THAT(fluxSum, Catch::Matchers::WithinAbs(0.0, std::abs(columnFluxes[mainColumnIndex]) * 1e-9));
    CHECK(std::abs(columnFluxes[mainColumnIndex]) > std::abs(columnFluxes[secondaryColumn]));
    CHECK(std::abs(columnFluxes[mainColumnIndex]) > std::abs(columnFluxes[otherLateralColumn]));

    auto columnFluxDensities = circuit.calculate_column_magnetic_flux_densities(magnetic, {1.0, 0.0});
    REQUIRE(columnFluxDensities.size() == 3);
    CHECK(std::abs(columnFluxDensities[mainColumnIndex]) > 0);
}

TEST_CASE("ReluctanceNetwork_MagnetizingInductance_PlacementGate", "[physical-model][magnetic-circuit][multi-column][smoke-test]") {
    // Placing the SECONDARY on a lateral leg must not change the primary's
    // magnetizing inductance (its driving-point reluctance is unchanged), and the
    // full MagnetizingInductance pipeline must agree with the lumped value.
    auto baselineMagnetic = buildTwoWindingMagnetic(std::nullopt);
    auto lateralSecondaryMagnetic = buildTwoWindingMagnetic(2);

    MagnetizingInductance magnetizingInductanceModel;
    double baselineInductance = magnetizingInductanceModel
                                    .calculate_inductance_from_number_turns_and_gapping(
                                        baselineMagnetic.get_mutable_core(), baselineMagnetic.get_mutable_coil(), nullptr)
                                    .get_magnetizing_inductance()
                                    .get_nominal()
                                    .value();
    double lateralSecondaryInductance = magnetizingInductanceModel
                                            .calculate_inductance_from_number_turns_and_gapping(
                                                lateralSecondaryMagnetic.get_mutable_core(), lateralSecondaryMagnetic.get_mutable_coil(), nullptr)
                                            .get_magnetizing_inductance()
                                            .get_nominal()
                                            .value();
    CHECK_THAT(lateralSecondaryInductance, Catch::Matchers::WithinRel(baselineInductance, 1e-6));

    // Placing the PRIMARY on a lateral leg does change it: the lateral driving-point
    // reluctance is higher (smaller leg area, longer return path).
    auto lateralPrimaryMagnetic = buildTwoWindingMagnetic(std::nullopt);
    auto lateralPrimaryCoil = lateralPrimaryMagnetic.get_coil();
    lateralPrimaryCoil.get_mutable_functional_description()[0].set_winding_window(2);
    lateralPrimaryMagnetic.set_coil(lateralPrimaryCoil);
    double lateralPrimaryInductance = magnetizingInductanceModel
                                          .calculate_inductance_from_number_turns_and_gapping(
                                              lateralPrimaryMagnetic.get_mutable_core(), lateralPrimaryMagnetic.get_mutable_coil(), nullptr)
                                          .get_magnetizing_inductance()
                                          .get_nominal()
                                          .value();
    CHECK(lateralPrimaryInductance < baselineInductance);
}

TEST_CASE("ReluctanceNetwork_InductanceMatrix_PlacementGate", "[physical-model][magnetic-circuit][multi-column][smoke-test]") {
    auto magnetic = buildTwoWindingMagnetic(2);
    magnetic.get_mutable_coil().set_core_columns(magnetic.get_core().get_processed_description()->get_columns());
    REQUIRE(magnetic.get_mutable_coil().wind());

    Inductance inductanceModel;
    auto matrix = inductanceModel.calculate_inductance_matrix(magnetic, 100000);
    auto magnitude = matrix.get_magnitude();
    double selfPrimary = magnitude["Primary"]["Primary"].get_nominal().value();
    double selfSecondary = magnitude["Secondary"]["Secondary"].get_nominal().value();
    double mutual = magnitude["Primary"]["Secondary"].get_nominal().value();

    CHECK(selfPrimary > 0);
    CHECK(selfSecondary > 0);
    // Leg-separated windings must not report ideal coupling.
    double couplingCoefficient = std::abs(mutual) / std::sqrt(selfPrimary * selfSecondary);
    CHECK(couplingCoefficient < 0.999);
    CHECK(magnitude["Primary"]["Secondary"].get_nominal().value() == magnitude["Secondary"]["Primary"].get_nominal().value());
}
