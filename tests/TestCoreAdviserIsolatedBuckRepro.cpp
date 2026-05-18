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
    if (inputs.get_design_requirements().get_isolation_sides()) {
        for (auto s : inputs.get_design_requirements().get_isolation_sides().value()) {
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
