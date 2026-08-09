// Regression coverage for ABT #606: OpenMagnetics::compat::migrate_pre_1_0
// translates pre-1.0 MAS enum strings ("Industrial" -> "industrial",
// "Phase-Shifted Full-Bridge Converter" -> "phaseShiftedFullBridgeConverter",
// "Printed" -> "printed", "two-piece set" -> "twoPieceSet", etc. — 110
// entries in MasMigration.cpp) but was only wired into Magnetic.h's
// from_file(), a filesystem-path loader nothing in the product calls. Every
// real entry point (every PyOM/WASM binding) constructs a Magnetic from an
// in-memory json via Magnetic::from_json, which never ran the migration —
// so any pre-1.0 MAS export was rejected outright with the generic
// quicktype "Input JSON does not conform to schema!" before any physics
// ever ran. Reported symptom: "Error calculating impedance" on a magnetic
// saved months earlier, plus (very likely, same root cause) a silently
// blank 3D core render — Core3DVisualizer.vue catches any WASM exception
// and just doesn't draw.
//
// user_report_abt606_legacy_casing_magnetic.json is the literal magnetic
// sub-object from the user's shared export (custom_magnetic (3).json,
// PSFB / Custom E 38/8/25). It carries FOUR legacy-cased strings
// simultaneously (coil.groupsDescription[0].type = "Printed",
// core.functionalDescription.type = "two-piece set",
// core.geometricalDescription[0/1].type = "half set") — the exact document
// that failed sweep_impedance_over_frequency with the schema exception
// before this fix.
//
// Note: core.functionalDescription.type alone would NOT catch a regression
// here — OpenMagnetics::Core has its own constructor (Core.cpp:29) that
// already calls migrate_pre_1_0 independently of Magnetic::from_json. The
// coil.groupsDescription[0].type field is the one that actually pins this
// fix: OpenMagnetics::Coil's free from_json (Coil.h:726, as opposed to its
// constructors at Coil.cpp:276/291) never migrated on its own, and is only
// covered because Magnetic::from_json now migrates the WHOLE tree before
// descending into "coil"/"core".
//
// Calls the free function from_json(const json&, Magnetic&) directly rather
// than json::get<Magnetic>() — the latter's two-parameter deduction
// (get<ValueTypeCV, ValueType>()) misbehaves inside Catch2's expression-
// decomposing macros and inside a lambda return-type-deduction context on
// this GCC/nlohmann-json combination (reproducible: fails when called as
// `x = j.get<Magnetic>()` or inside `[]{ return j.get<Magnetic>(); }`,
// succeeds as a direct-initializer `Magnetic x = j.get<Magnetic>();`).
// from_json is the exact function .get<T>() would call either way.

#include "constructive_models/Magnetic.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <fstream>

using json = nlohmann::json;
using namespace OpenMagneticsTesting;

TEST_CASE("Test_MasMigration_Magnetic_From_Json_Accepts_Pre_1_0_Enum_Casing", "[constructive-model][magnetic][mas-migration][smoke-test]") {
    auto testDataPath = get_test_data_path(std::source_location::current(), "user_report_abt606_legacy_casing_magnetic.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    json magJson;
    file >> magJson;

    // Confirm the fixture still carries the legacy strings this test exists
    // for — if a future edit accidentally normalizes the fixture, this test
    // would otherwise keep passing for the wrong reason.
    REQUIRE(magJson.at("coil").at("groupsDescription").at(0).at("type").get<std::string>() == "Printed");
    REQUIRE(magJson.at("core").at("functionalDescription").at("type").get<std::string>() == "two-piece set");

    OpenMagnetics::Magnetic magnetic;
    REQUIRE_NOTHROW(OpenMagnetics::from_json(magJson, magnetic));

    double effectiveArea = magnetic.get_core().get_processed_description().value().get_effective_parameters().get_effective_area();
    REQUIRE(effectiveArea > 0);
    REQUIRE(magnetic.get_coil().get_functional_description().size() > 0);
}

TEST_CASE("Test_MasMigration_Magnetic_From_Json_Rejects_Genuinely_Unknown_Enum_Value", "[constructive-model][magnetic][mas-migration][smoke-test]") {
    // The fix must not turn from_json into something that silently accepts
    // ANY string — migrate_pre_1_0 only rewrites the 110 known pre-1.0 forms;
    // a value that was never valid at any point must still throw, not be
    // swallowed into a default. Guards against a fix that (e.g.) wraps the
    // whole function in a try/catch instead of translating known aliases.
    auto testDataPath = get_test_data_path(std::source_location::current(), "user_report_abt606_legacy_casing_magnetic.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    json magJson;
    file >> magJson;

    magJson["core"]["functionalDescription"]["type"] = "not-a-real-core-type-ever";
    OpenMagnetics::Magnetic rejected;
    REQUIRE_THROWS(OpenMagnetics::from_json(magJson, rejected));
}

TEST_CASE("Test_MasMigration_Magnetic_Vector_From_Json_Accepts_Pre_1_0_Enum_Casing", "[constructive-model][magnetic][mas-migration][smoke-test]") {
    // from_json(const json&, std::vector<Magnetic>&) is the untouched twin of
    // the singular overload fixed for ABT #606 — it iterates elements without
    // migrating each one. Same fixture, wrapped in a list, exercised through
    // the vector overload specifically.
    auto testDataPath = get_test_data_path(std::source_location::current(), "user_report_abt606_legacy_casing_magnetic.json");
    std::ifstream file(testDataPath);
    REQUIRE(file.good());
    json magJson;
    file >> magJson;
    REQUIRE(magJson.at("coil").at("groupsDescription").at(0).at("type").get<std::string>() == "Printed");

    json magneticsArrayJson = json::array({magJson});
    std::vector<OpenMagnetics::Magnetic> magnetics;
    REQUIRE_NOTHROW(OpenMagnetics::from_json(magneticsArrayJson, magnetics));

    REQUIRE(magnetics.size() == 1);
    double effectiveArea = magnetics[0].get_core().get_processed_description().value().get_effective_parameters().get_effective_area();
    REQUIRE(effectiveArea > 0);
    REQUIRE(magnetics[0].get_coil().get_functional_description().size() > 0);
}
