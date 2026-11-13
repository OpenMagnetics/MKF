// ABT #2376: the core adviser's per-filter cull (evaluate_and_cull) treated every OpenMagneticsException
// as "this candidate cannot work". An error in the caller's inputs, which no candidate can avoid, then
// culled the whole pool and the advise returned no cores and no reason. Such an error now stops the
// advise with its own message.
#include "advisers/CoreAdviser.h"
#include "physical_models/Temperature.h"
#include "processors/Inputs.h"
#include "support/Exceptions.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <fstream>
#include <map>
#include <source_location>
#include <string>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::ContainsSubstring;

namespace {

// The first concentric (two-piece) cores of the adviser tests' inventory.
std::vector<Core> first_concentric_test_cores(size_t count) {
    auto inventoryPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_cores.ndjson");
    std::ifstream ndjsonFile(inventoryPath);
    std::string jsonLine;
    std::vector<Core> cores;
    while (cores.size() < count && std::getline(ndjsonFile, jsonLine)) {
        json coreJson = json::parse(jsonLine);
        try {
            Core core(coreJson, false, true, false);
            if (core.get_functional_description().get_type() == CoreType::TWO_PIECE_SET) {
                cores.push_back(core);
            }
        }
        catch (const CoreShapeNotFoundException&) {
            continue;  // as load_test_data in TestCoreAdviser.cpp: shapes missing from the database are skipped
        }
    }
    return cores;
}

} // namespace

TEST_CASE("Core adviser: a heatsink without its interface resistance stops the advise with that message",
          "[adviser][core-adviser][available-cores][temperature][abt-2376]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    clear_databases();
    REQUIRE(settings.get_core_adviser_enable_temperature_filter());

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(
        100000, 100e-6, 25, WaveformLabel::SINUSOIDAL, 20, 0.5, 0, std::vector<double>{});
    MAS::Cooling heatsink;
    heatsink.set_thermal_resistance(2.5);
    { MAS::Dimensions d; d.set_width(0.05); d.set_height(0.05); d.set_depth(0.02); heatsink.set_dimensions(d); }
    for (auto& operatingPoint : inputs.get_mutable_operating_points()) {
        operatingPoint.get_mutable_conditions().set_cooling(heatsink);
    }

    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 1;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

    auto cores = first_concentric_test_cores(40);
    REQUIRE(cores.size() == 40);

    CoreAdviser coreAdviser;
    coreAdviser.set_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
    size_t returned = 0;
    try {
        returned = coreAdviser.get_advised_core(inputs, weights, &cores, 5).size();
        FAIL("the advise returned " << returned << " cores instead of refusing the heatsink");
    }
    catch (const InvalidInputException& e) {
        CHECK(e.code() == ErrorCode::MISSING_DATA);
        CHECK_THAT(std::string(e.what()), ContainsSubstring("interfaceThermalResistance"));
    }
    settings.reset();
}

// The heatsink model mounts on a top yoke; a toroid's thermal network has none. That is a missing model,
// said with a typed exception, not a std::runtime_error (which escaped the adviser's cull).
TEST_CASE("Temperature: a heatsink on a toroid is refused as a model that does not exist yet",
          "[temperature][cooling][abt-2376]") {
    auto jsonPath = OpenMagneticsTesting::get_test_data_path(
        std::source_location::current(), "toroidal_inductor_round_wire_multilayer.json");
    auto mas = OpenMagneticsTesting::mas_loader(jsonPath);
    auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
    REQUIRE(magnetic.get_core().get_shape_family() == CoreShapeFamily::T);

    TemperatureConfig config;
    config.ambientTemperature = 25;
    config.coreLosses = 1.0;
    config.windingLosses = 1.0;
    config.coreOnly = true;
    MAS::Cooling heatsink;
    heatsink.set_thermal_resistance(2.5);
    heatsink.set_interface_thermal_resistance(3e-4);
    { MAS::Dimensions d; d.set_width(0.05); d.set_height(0.05); d.set_depth(0.02); heatsink.set_dimensions(d); }
    config.masCooling = heatsink;

    REQUIRE_THROWS_AS(Temperature(magnetic, config).calculateTemperatures(), NotImplementedException);
    REQUIRE_THROWS_WITH(Temperature(magnetic, config).calculateTemperatures(), ContainsSubstring("toroid"));
}
