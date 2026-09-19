// Real-winding connection blocking: designs MVB++ builds as FEM-grade CAD must keep getting
// blocking applied by magnetic_autocomplete. MVB++ refuses a real-winding build whose blocking was
// not applied (magnetic_autocomplete_safe), so a regression here turns its fixtures red.

#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include "support/Settings.h"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <source_location>

using json = nlohmann::json;

namespace {

json load_test_data(const std::string& name) {
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("testData")
                           .append(name));
    REQUIRE(file.good());
    return json::parse(file);
}

}  // namespace

// ABT #1194: MVB++'s cm37 fixture (E16, Primary + Secondary, 22 turns x 2 parallels each,
// interleaved). Built exactly as MVB++'s magnetic_autocomplete_safe builds it: Core and Coil from
// the stored MAS (the coil NOT re-wound by its constructor), then magnetic_autocomplete with real
// winding geometry on and no fit relaxations.
TEST_CASE("Real winding: the interleaved N-filar E16 (MVB++ cm37) gets connection blocking",
          "[constructive-model][coil][real-winding][abt1194]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    const auto masJson = load_test_data("abt1194_interleaved_nfilar_e16.json");
    const auto& magneticJson = masJson.at("magnetic");

    settings.set_coil_use_real_winding_geometry(true);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    auto enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
    REQUIRE(enriched.get_coil().is_real_winding_blocking_applied());
    settings.reset();
}
