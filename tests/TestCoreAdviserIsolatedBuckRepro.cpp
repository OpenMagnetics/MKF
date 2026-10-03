#include <fstream>
#include <source_location>
#include "advisers/MagneticAdviser.h"
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"
/**
 * Repro for the silent "no suitable core found" failure observed when the
 * WebFrontend Isolated Buck wizard reaches the Core Adviser with its
 * default inputs.
 *
 * Fixture: tests/testData/isolated_buck_coreadviser_inputs.json
 *   Captured directly from MagneticBuilder/src/stores/taskQueue.js::adviseCore
 *   right before calling mkf.calculate_advised_cores().
 *
 * Expected (the contract the WebFrontend test enforces): with default wizard
 * inputs, Core Adviser MUST return at least one candidate. A zero-result
 * outcome with empty log is a silent-failure bug we want to surface.
 */
#include "advisers/CoreAdviser.h"
#include "processors/Inputs.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <source_location>

using namespace MAS;
using namespace OpenMagnetics;

TEST_CASE("Test_CoreAdviser_IsolatedBuck_DefaultWizardInputs_HasResults",
          "[adviser][core-adviser][standard-cores][isolated-buck][bug-repro]") {
    clear_databases();

    // Load the captured fixture.
    auto path = OpenMagneticsTesting::get_test_data_path(
        std::source_location::current(), "isolated_buck_coreadviser_inputs.json");
    std::ifstream f(path);
    REQUIRE(f.is_open());
    json captured;
    f >> captured;

    // The capture stores {inputs, weights, adviserSettings, mkfResult, raw}.
    REQUIRE(captured.contains("inputs"));
    REQUIRE(captured.contains("weights"));

    OpenMagnetics::Inputs inputs(captured["inputs"]);

    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    double externalSum = 0.0;
    for (auto it = captured["weights"].begin(); it != captured["weights"].end(); ++it) {
        externalSum += it.value().get<double>();
    }
    REQUIRE(externalSum > 0.0);
    for (auto it = captured["weights"].begin(); it != captured["weights"].end(); ++it) {
        CoreAdviser::CoreAdviserFilters filter;
        from_json(it.key(), filter);
        weights[filter] = it.value().get<double>() / externalSum;
    }

    // Mirror the WebLibMKF settings (calculate_advised_cores in libMKF.cpp).
    // Wizard defaults: allow concentric cores, no toroidal, no stock filter.
    auto& settings = Settings::GetInstance();
    settings.set_core_adviser_include_distributed_gaps(false);
    settings.set_core_adviser_include_stacks(false);
    settings.set_use_toroidal_cores(false);
    settings.set_use_concentric_cores(true);
    settings.set_use_only_cores_in_stock(false);
    settings.set_core_adviser_include_margin(true);
    settings.set_coil_delimit_and_compact(true);

    CoreAdviser adviser;
    adviser.set_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);

    std::cout << "[REPRO] turnsRatios.size = "
              << inputs.get_design_requirements().get_turns_ratios().size() << std::endl;
    std::cout << "[REPRO] isolationSides = ";
    auto isolationSides = inputs.get_design_requirements().get_isolation_sides();
    if (isolationSides) {
        // MAS getters return by value: iterating .value() of the temporary directly binds to
        // storage destroyed at the end of the full expression (gcc warns -Wdangling-reference).
        auto isolationSidesValue = isolationSides.value();
        for (auto s : isolationSidesValue) {
            std::cout << static_cast<int>(s) << " ";
        }
    }
    std::cout << std::endl;
    std::cout << "[REPRO] operatingPoints[0].excitations.size = "
              << inputs.get_operating_points()[0].get_excitations_per_winding().size() << std::endl;
    if (inputs.get_design_requirements().get_magnetizing_inductance().get_minimum()) {
        std::cout << "[REPRO] L_min = "
                  << inputs.get_design_requirements().get_magnetizing_inductance().get_minimum().value()
                  << std::endl;
    }

    auto results = adviser.get_advised_core(inputs, weights, /*maximumNumberResults=*/1);

    auto log = read_log();
    std::cout << "[REPRO] MKF log size = " << log.size() << std::endl;
    if (!log.empty()) std::cout << "[REPRO] MKF log:\n" << log << std::endl;
    std::cout << "[REPRO] results.size = " << results.size() << std::endl;

    // The bug: this currently fails with results.size()==0 and an empty log.
    REQUIRE(results.size() > 0);
}

namespace WebAdviserRequests1542 {

// The web drops JSON nulls before handing the inputs to the engine.
void strip_nulls(json& value) {
    if (value.is_array()) {
        for (auto& element : value) {
            strip_nulls(element);
        }
    }
    else if (value.is_object()) {
        for (auto it = value.begin(); it != value.end();) {
            if (it->is_null() || (it->is_string() && it->get<std::string>() == "null")) {
                it = value.erase(it);
            }
            else {
                strip_nulls(*it);
                ++it;
            }
        }
    }
}

struct WebAdviserRequest {
    OpenMagnetics::Inputs inputs;
    std::map<MagneticFilters, double> weights;
    CoreAdviser::CoreAdviserModes coreMode;
    size_t maximumNumberResults;
};

// Applies the settings the web sets through set_settings (WebLibMKF) and returns the request.
WebAdviserRequest load_web_adviser_request(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open adviser request " + path.string());
    }
    json request = json::parse(file);
    auto webSettings = request.at("settings");
    auto adviserSettings = request.at("adviserSettings");
    auto magneticAdviserSettings = request.at("magneticAdviserSettings");

    settings.reset();
    settings.set_magnetizing_inductance_include_air_inductance(webSettings.at("magnetizingInductanceIncludeAirInductance"));
    settings.set_coil_allow_margin_tape(webSettings.at("coilAllowMarginTape"));
    settings.set_coil_allow_insulated_wire(webSettings.at("coilAllowInsulatedWire"));
    settings.set_coil_fill_sections_with_margin_tape(webSettings.at("coilFillSectionsWithMarginTape"));
    settings.set_coil_wind_even_if_not_fit(webSettings.at("coilWindEvenIfNotFit"));
    settings.set_coil_delimit_and_compact(webSettings.at("coilDelimitAndCompact"));
    settings.set_coil_only_one_turn_per_layer_in_contiguous_rectangular(webSettings.at("coilOnlyOneTurnPerLayerInContiguousRectangular"));
    settings.set_coil_try_rewind(webSettings.at("coilTryRewind"));
    settings.set_coil_maximum_layers_planar(webSettings.at("coilMaximumLayersPlanar"));
    settings.set_preferred_core_material_ferrite_manufacturer(webSettings.at("preferredCoreMaterialFerriteManufacturer"));
    settings.set_preferred_core_material_powder_manufacturer(webSettings.at("preferredCoreMaterialPowderManufacturer"));
    settings.set_preferred_wire_standard(WireStandard::IEC_60317);  // the web always sends IEC 60317
    settings.set_coil_include_additional_coordinates(webSettings.at("coilIncludeAdditionalCoordinates"));
    settings.set_coil_use_real_winding_geometry(webSettings.at("coilUseRealWindingGeometry"));
    settings.set_use_only_cores_in_stock(adviserSettings.at("useOnlyCoresInStock"));
    settings.set_coil_adviser_maximum_number_wires(webSettings.at("coilAdviserMaximumNumberWires"));
    settings.set_core_adviser_include_margin(webSettings.at("coreIncludeMargin"));
    settings.set_core_adviser_include_stacks(adviserSettings.at("allowStacks"));
    settings.set_core_adviser_include_distributed_gaps(adviserSettings.at("allowDistributedGaps"));
    settings.set_use_toroidal_cores(adviserSettings.at("allowToroidalCores"));
    settings.set_use_concentric_cores(webSettings.at("useConcentricCores"));

    auto inputsJson = request.at("inputs");
    strip_nulls(inputsJson);

    std::map<std::string, double> weightsByName = magneticAdviserSettings.at("weights");
    double weightSum = 0;
    for (auto const& [name, weight] : weightsByName) {
        weightSum += weight;
    }
    std::map<MagneticFilters, double> weights;
    for (auto const& [name, weight] : weightsByName) {
        MagneticFilters filter;
        from_json(name, filter);
        weights[filter] = weight / weightSum;
    }

    CoreAdviser::CoreAdviserModes coreMode;
    from_json(adviserSettings.at("coreAdviseMode"), coreMode);

    return {OpenMagnetics::Inputs(inputsJson), weights, coreMode, magneticAdviserSettings.at("maximumNumberResults").get<size_t>()};
}

std::vector<std::pair<OpenMagnetics::Mas, double>> run_web_adviser_request(WebAdviserRequest& request) {
    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(request.coreMode);
    return adviser.get_advised_magnetic(request.inputs, request.weights, request.maximumNumberResults);
}

}  // namespace WebAdviserRequests1542

// ABT #1542: the isolated buck-boost wizard (three windings of 13/11/8 turns, basic insulation
// for 400 V mains, 15 uH) got 0 designs. add_initial_turns_by_inductance sized N for 15 uH, then
// snapped it up to the next valid turns-ratio combination without re-solving the gap, so L grew
// by (N_snapped / N)^2 and the inductance filter rejected 546 of 557 candidates.
TEST_CASE("MagneticAdviser advises the isolated buck-boost web request", "[adviser][magnetic-adviser][standard-cores][abt1542]") {
    clear_databases();
    auto request = WebAdviserRequests1542::load_web_adviser_request(OpenMagneticsTesting::get_test_data_path(
        std::source_location::current(), "abt1542_isolated_buckboost_adviser_request.json"));
    auto results = WebAdviserRequests1542::run_web_adviser_request(request);
    settings.reset();
    REQUIRE(results.size() >= 1);
}

