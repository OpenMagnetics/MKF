#include "constructive_models/Core.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

namespace {

// The "Dummy" core material is the synthetic placeholder a core carries before a real grade is
// chosen (CoreAdviser exploration, the web's shapes table). Core::resolve_material builds it from
// scratch and writes it back into the functional description, so every required member of
// MAS::CoreMaterial must be set there: quicktype leaves an enum member uninitialised, and
// serialising an indeterminate CoreMaterialType threw "Unexpected value in enumeration" for every
// shape in calculate_all_core_data_from_shapes, which then returned no cores at all.
TEST_CASE("Dummy core material states its material type", "[constructive-model][core][dummy-material][smoke-test]") {
    auto material = Core::resolve_material(std::string("Dummy"));
    REQUIRE(material.get_name() == "Dummy");
    REQUIRE(material.get_material() == MaterialType::FERRITE);
    REQUIRE(material.get_type() == CoreMaterialType::CUSTOM);

    json materialJson;
    REQUIRE_NOTHROW(to_json(materialJson, material));
    REQUIRE(materialJson.at("type") == "custom");
}

// Same construction as WebLibMKF's calculate_all_core_data_from_shapes: shape + "Dummy" material,
// processed, then serialised with its written-back material.
TEST_CASE("A core with the Dummy material processes and serialises for each core type", "[constructive-model][core][dummy-material][smoke-test]") {
    const std::vector<std::string> shapeNames = {"E 42/21/15", "PQ 32/30", "T 25/15/10"};
    for (const auto& shapeName : shapeNames) {
        INFO(shapeName);
        auto shape = find_core_shape_by_name(shapeName);
        Core core;
        CoreFunctionalDescription functionalDescription;
        functionalDescription.set_shape(shape);
        functionalDescription.set_material("Dummy");
        functionalDescription.set_number_stacks(1);
        if (shape.get_family() == CoreShapeFamily::T) {
            functionalDescription.set_type(CoreType::TOROIDAL);
        }
        else {
            functionalDescription.set_type(CoreType::TWO_PIECE_SET);
        }
        core.set_functional_description(functionalDescription);
        core.process_data();

        json coreJson;
        REQUIRE_NOTHROW(to_json(coreJson, core));
        REQUIRE(coreJson.at("functionalDescription").at("material").at("name") == "Dummy");
        REQUIRE(coreJson.at("functionalDescription").at("material").at("type") == "custom");
    }
}

}  // namespace
