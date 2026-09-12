// ABT #1170 (WP1 "Spacers end to end"): the spacer elements MKF emits into a core's
// geometricalDescription are the only description of the plastic shims that MVB++ draws,
// STEP exports and OMFEM meshes. Before this ticket they validated against NEITHER branch
// of MAS's `coreGeometricalDescriptionElement` oneOf:
//
//   * the material was written to the PIECE field `material` instead of the spacer field
//     `insulationMaterial` (which core/spacer.json REQUIRES), and
//   * a `rotation` was set, which core/spacer.json does not define (and it is
//     `additionalProperties: false`).
//
// and only the LATERAL columns got one, although `create_spacer_gapping` makes one ADDITIVE
// gap per column: the centre leg's gap had no solid.
//
// The validation here goes through MAS's own validator (MAS/scripts/validate-samples.py,
// jsonschema Draft 2020-12 with a $id registry over MAS + PEAS). MKF has no C++ JSON-Schema
// validator and this test does not invent one.
#include <cmath>
#include <cstdlib>
#include <limits>
#include <filesystem>
#include <fstream>
#include <source_location>
#include <string>
#include <vector>

#include "constructive_models/Core.h"
#include "Constants.h"
#include "Fixtures.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

namespace {

std::filesystem::path testsDir() {
    return std::filesystem::path{std::source_location::current().file_name()}.parent_path();
}

std::filesystem::path masDir() {
    return std::filesystem::canonical(testsDir() / ".." / "MAS");
}

// quicktype's generated `to_json` writes every optional member, absent ones as JSON `null`.
// MAS's schemas model "absent" as the key being absent (every object is
// `additionalProperties: false` and the optional members are typed, not nullable), so the
// generated document has to be normalised before it can be validated. This drops nothing a
// schema would have accepted: a REQUIRED field emitted as null is removed here and then
// fails validation as missing, which is the same verdict.
void strip_nulls(json& node) {
    if (node.is_object()) {
        for (auto it = node.begin(); it != node.end();) {
            if (it.value().is_null()) {
                it = node.erase(it);
            }
            else {
                strip_nulls(it.value());
                ++it;
            }
        }
    }
    else if (node.is_array()) {
        for (auto& element : node) {
            strip_nulls(element);
        }
    }
}

// Run MAS's validator over every *.json in `documentsDir` against `schemaRelativeToMas`.
// Returns the process exit status: 0 = every document valid.
int validate_with_mas(const std::filesystem::path& documentsDir,
                      const std::string& schemaRelativeToMas) {
    auto command = std::string("python3 ") + (masDir() / "scripts" / "validate-samples.py").string() +
                   " --root " + masDir().string() +
                   " --samples " + documentsDir.string() +
                   " --schema " + schemaRelativeToMas;
    return std::system(command.c_str());
}

}  // namespace

TEST_CASE("Core_Additive_Gap_Spacers_Are_Schema_Valid_And_One_Per_Column",
          "[constructive-model][core][geometrical-description][spacer][abt1170]") {
    const double spacerThickness = 0.5e-3;
    auto core = OpenMagneticsTesting::get_quick_core("E 42/21/15",
                                                     OpenMagneticsTesting::get_spacer_gap(spacerThickness));

    REQUIRE(core.get_processed_description());
    auto columns = core.get_processed_description().value().get_columns();
    REQUIRE(columns.size() == 3u);

    // ---------------------------------------------------------------- one spacer per column
    auto spacers = core.get_spacers();
    INFO("spacers emitted: " << spacers.size());
    REQUIRE(spacers.size() == columns.size());

    OpenMagnetics::Constants constants;
    const double protrusion = 1 + constants.spacerProtudingPercentage;

    for (size_t index = 0; index < columns.size(); ++index) {
        auto& column = columns[index];
        // The spacer for this column is the one whose x is nearest to it: the emission order
        // follows the column order, but the test must not depend on that.
        size_t closest = 0;
        double closestDistance = std::numeric_limits<double>::max();
        for (size_t candidate = 0; candidate < spacers.size(); ++candidate) {
            double distance = std::fabs(spacers[candidate].get_coordinates()[0] - column.get_coordinates()[0]);
            if (distance < closestDistance) {
                closestDistance = distance;
                closest = candidate;
            }
        }
        auto& spacer = spacers[closest];
        UNSCOPED_INFO("column " << index << " ("
                      << (column.get_type() == ColumnType::CENTRAL ? "central" : "lateral")
                      << ") width " << column.get_width() << " depth " << column.get_depth()
                      << " at x " << column.get_coordinates()[0]);

        REQUIRE(spacer.get_type() == CoreGeometricalDescriptionElementType::SPACER);
        REQUIRE(spacer.get_dimensions());
        auto dimensions = spacer.get_dimensions().value();
        REQUIRE(dimensions.size() == 3u);

        // [column width, thickness, column depth] — the footprint of the leg it separates,
        // grown by Constants::spacerProtudingPercentage so the shim can be seen and gripped
        // outside the set. The thickness is the additive gap itself, exactly.
        CHECK_THAT(dimensions[0], Catch::Matchers::WithinRel(column.get_width() * protrusion, 1e-9));
        CHECK_THAT(dimensions[1], Catch::Matchers::WithinAbs(spacerThickness, 1e-12));
        CHECK_THAT(dimensions[2], Catch::Matchers::WithinRel(column.get_depth() * protrusion, 1e-9));

        // Centred on the column it separates. The lateral shims are pushed outwards by the
        // protruding part (they overhang the outside of the set, not the winding window), so
        // only the centre leg's shim is exactly concentric with its column; the laterals stay
        // within the protrusion they were grown by.
        double centringToleranceX = column.get_type() == ColumnType::CENTRAL
                                        ? 1e-12
                                        : column.get_width() * constants.spacerProtudingPercentage;
        CHECK_THAT(spacer.get_coordinates()[0],
                   Catch::Matchers::WithinAbs(column.get_coordinates()[0], centringToleranceX));
        CHECK_THAT(spacer.get_coordinates()[1],
                   Catch::Matchers::WithinAbs(column.get_coordinates()[1], 1e-12));
        CHECK_THAT(spacer.get_coordinates()[2],
                   Catch::Matchers::WithinAbs(column.get_coordinates()[2], 1e-12));

        // The shim covers the whole leg in depth: nothing of the column rests on air.
        CHECK(dimensions[2] / 2 - std::fabs(spacer.get_coordinates()[2] - column.get_coordinates()[2])
              >= column.get_depth() / 2);

        // core/spacer.json REQUIRES insulationMaterial and defines no `material`/`rotation`.
        REQUIRE(spacer.get_insulation_material());
        CHECK(std::get<std::string>(spacer.get_insulation_material().value()) == "plastic");
        CHECK(!spacer.get_material());
        CHECK(!spacer.get_rotation());
        CHECK(!spacer.get_shape());
        CHECK(!spacer.get_machining());
    }

    // ------------------------------------------------------------------- schema validation
    auto outputDir = testsDir() / ".." / "output" / "abt1170_spacers";
    std::filesystem::remove_all(outputDir);
    std::filesystem::create_directories(outputDir);

    json coreJson;
    MAS::to_json(coreJson, static_cast<const MAS::MagneticCore&>(core));
    strip_nulls(coreJson);
    {
        std::ofstream file(outputDir / "core.json");
        file << coreJson.dump(2);
    }
    INFO("processed core written to " << (outputDir / "core.json").string());
    // The whole processed core, spacers included, against MAS's core schema. The spacers are
    // reached through `geometricalDescription`'s oneOf(piece, spacer): an element that is
    // neither fails here.
    CHECK(validate_with_mas(outputDir, "schemas/magnetic/core.json") == 0);

    // ...and each spacer element on its own, so a failure names the spacer schema directly
    // instead of "does not match any of the oneOf branches".
    auto spacersDir = testsDir() / ".." / "output" / "abt1170_spacer_elements";
    std::filesystem::remove_all(spacersDir);
    std::filesystem::create_directories(spacersDir);
    for (size_t index = 0; index < spacers.size(); ++index) {
        json spacerJson;
        MAS::to_json(spacerJson, spacers[index]);
        strip_nulls(spacerJson);
        std::ofstream file(spacersDir / ("spacer_" + std::to_string(index) + ".json"));
        file << spacerJson.dump(2);
    }
    CHECK(validate_with_mas(spacersDir, "schemas/magnetic/core/spacer.json") == 0);
}


// The fixture the downstream packages build on: the autocompleted MAS of the WP1 design
// (E 42/21/15, 0.5 mm ADDITIVE gap on all three columns). MVB++'s SpacerBuilder test and
// OMFEM's spacer-region validation both start from it, so it has to be schema-valid and it
// has to actually contain the three shims.
//
// Regenerating it: run this case and copy output/abt1170_fixture/magnetic.json into the
// "data" member of the `additive-gap-e-core` line of tests/fixtures/mas-additive-gap.ndjson.
TEST_CASE("Core_Additive_Gap_Fixture_Is_The_Autocompleted_Design_And_Schema_Valid",
          "[constructive-model][core][geometrical-description][spacer][abt1170]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 42/21/15",
                        OpenMagneticsTesting::get_spacer_gap(0.5e-3), {10});

    // NOT `MAS::to_json(..., static_cast<const MAS::Magnetic&>(magnetic))`: OpenMagnetics::
    // Magnetic SHADOWS the base's core/coil members (the ABT #611 family of bugs), so the base
    // serializes both as null and the document collapses to `{}` — which, since magnetic.json
    // has no required members, then "validates" perfectly. OpenMagnetics' own to_json writes
    // the derived members.
    json generated = magnetic;
    strip_nulls(generated);
    REQUIRE(generated.contains("core"));
    REQUIRE(generated.at("core").contains("geometricalDescription"));

    auto outputDir = testsDir() / ".." / "output" / "abt1170_fixture";
    std::filesystem::remove_all(outputDir);
    std::filesystem::create_directories(outputDir);
    {
        std::ofstream file(outputDir / "magnetic.json");
        file << generated.dump(2);
    }
    INFO("regenerated fixture at " << (outputDir / "magnetic.json").string());
    CHECK(validate_with_mas(outputDir, "schemas/magnetic.json") == 0);

    // ...and the committed fixture IS that design, byte for byte: a fixture that has drifted
    // from what MKF now emits is a fixture the downstream packages are testing the past with.
    auto fixture = OpenMagneticsTesting::fixtures::get_json("additive-gap-e-core");
    CHECK(fixture == generated);
    auto fixtureDir = testsDir() / ".." / "output" / "abt1170_fixture_committed";
    std::filesystem::remove_all(fixtureDir);
    std::filesystem::create_directories(fixtureDir);
    {
        std::ofstream file(fixtureDir / "magnetic.json");
        file << fixture.dump(2);
    }
    CHECK(validate_with_mas(fixtureDir, "schemas/magnetic.json") == 0);

    REQUIRE(fixture.contains("core"));
    REQUIRE(fixture.at("core").contains("geometricalDescription"));
    size_t spacersInFixture = 0;
    for (auto& element : fixture.at("core").at("geometricalDescription")) {
        if (element.at("type") == "spacer") {
            ++spacersInFixture;
            // The contract the downstream packages read: the spacer field, not the piece one.
            CHECK(element.contains("insulationMaterial"));
            CHECK(!element.contains("material"));
            CHECK(!element.contains("rotation"));
            REQUIRE(element.at("dimensions").size() == 3u);
            CHECK_THAT(element.at("dimensions")[1].get<double>(),
                       Catch::Matchers::WithinAbs(0.5e-3, 1e-12));
        }
    }
    CHECK(spacersInFixture == 3u);
}
