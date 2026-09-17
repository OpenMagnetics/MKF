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
#include "constructive_models/Mas.h"
#include "physical_models/StrayCapacitance.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

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

// ABT #1200: older MKF wrote its synthesised core spacer in two legacy forms, and saved designs
// still carry both. Before ABT #1170: {"type": "spacer", "material": "plastic"} -- the piece-only
// `material` key, which core/spacer.json forbids, and no insulationMaterial, which it requires.
// Between #1170 and #1200: {"type": "spacer", "insulationMaterial": "plastic"} -- a name the
// insulation database does not carry. Both load as the dielectric MKF now stamps
// (Defaults().defaultSpacerMaterial, PET). The match is on those exact legacy forms only.
TEST_CASE("Test_MasMigration_Legacy_Plastic_Spacer_Becomes_Valid_Spacer", "[constructive-model][core][mas-migration][abt1200]") {
    json elements = json::array({
        json{{"type", "halfSet"}, {"material", "N87"}, {"shape", "E 42/21/15"}},
        json{{"type", "spacer"}, {"material", "plastic"}, {"dimensions", {0.01, 0.001, 0.01}}, {"coordinates", {0.02, 0, 0}}},
    });
    OpenMagnetics::compat::migrate_pre_1_0(elements);

    const auto& spacer = elements.at(1);
    CHECK(spacer.at("insulationMaterial") == OpenMagnetics::Defaults().defaultSpacerMaterial);
    CHECK_FALSE(spacer.contains("material"));
    CHECK(spacer.at("dimensions") == json({0.01, 0.001, 0.01}));
    // A piece's `material` is its own field and is not touched.
    CHECK(elements.at(0).at("material") == "N87");
}

TEST_CASE("Test_MasMigration_Plastic_InsulationMaterial_Spacer_Becomes_PET", "[constructive-model][core][mas-migration][abt1200]") {
    json elements = json::array({
        json{{"type", "halfSet"}, {"material", "N87"}, {"shape", "E 42/21/15"}},
        json{{"type", "spacer"}, {"insulationMaterial", "plastic"}, {"dimensions", {0.01, 0.001, 0.01}}, {"coordinates", {0.02, 0, 0}}},
    });
    OpenMagnetics::compat::migrate_pre_1_0(elements);

    const auto& spacer = elements.at(1);
    CHECK(spacer.at("insulationMaterial") == OpenMagnetics::Defaults().defaultSpacerMaterial);
    CHECK_FALSE(spacer.contains("material"));
    CHECK(spacer.at("coordinates") == json({0.02, 0, 0}));
    CHECK(elements.at(0).at("material") == "N87");
}

TEST_CASE("Test_MasMigration_Plastic_Inline_Record_Spacer_Becomes_PET", "[constructive-model][core][mas-migration][abt1200]") {
    // The #1170-era value round-tripped as an inline record, with no dielectric of its own.
    json spacer = json{{"type", "spacer"}, {"insulationMaterial", json{{"name", "plastic"}}},
                       {"dimensions", {0.01, 0.001, 0.01}}, {"coordinates", {0.02, 0, 0}}};
    OpenMagnetics::compat::migrate_pre_1_0(spacer);
    CHECK(spacer.at("insulationMaterial") == OpenMagnetics::Defaults().defaultSpacerMaterial);
    CHECK_FALSE(spacer.contains("material"));

    // A record named "plastic" that DOES carry its own relativePermittivity is real data: kept.
    json withData = json{{"type", "spacer"}, {"insulationMaterial", json{{"name", "plastic"}, {"relativePermittivity", 2.7}}}};
    json withDataBefore = withData;
    OpenMagnetics::compat::migrate_pre_1_0(withData);
    CHECK(withData == withDataBefore);
}

TEST_CASE("Test_MasMigration_Does_Not_Guess_Other_Spacer_Materials", "[constructive-model][core][mas-migration][abt1200]") {
    // An unknown legacy value is a malformed record, not the form MKF wrote: left exactly as it is.
    json unknown = json{{"type", "spacer"}, {"material", "cardboard"}};
    json unknownBefore = unknown;
    OpenMagnetics::compat::migrate_pre_1_0(unknown);
    CHECK(unknown == unknownBefore);

    // An unknown insulationMaterial is not mapped either.
    json unknownInsulation = json{{"type", "spacer"}, {"insulationMaterial", "cardboard"}};
    json unknownInsulationBefore = unknownInsulation;
    OpenMagnetics::compat::migrate_pre_1_0(unknownInsulation);
    CHECK(unknownInsulation == unknownInsulationBefore);

    // "plastic" on anything that is not a spacer is not this migration's business.
    json notSpacer = json{{"type", "bobbin"}, {"insulationMaterial", "plastic"}};
    json notSpacerBefore = notSpacer;
    OpenMagnetics::compat::migrate_pre_1_0(notSpacer);
    CHECK(notSpacer == notSpacerBefore);

    // A spacer that already declares its dielectric keeps it.
    json declared = json{{"type", "spacer"}, {"material", "plastic"}, {"insulationMaterial", "Kapton HN"}};
    json declaredBefore = declared;
    OpenMagnetics::compat::migrate_pre_1_0(declared);
    CHECK(declared == declaredBefore);

    // A spacer with no material at all is not invented one either.
    json bare = json{{"type", "spacer"}, {"dimensions", {0.01, 0.001, 0.01}}};
    json bareBefore = bare;
    OpenMagnetics::compat::migrate_pre_1_0(bare);
    CHECK(bare == bareBefore);
}


// ABT #1200, backward compatibility (Alf: "make sure that old files with plastic still work").
// Old saved designs carry MKF's earlier spacer forms, and every MKF path that turns a MAS document
// into a Core or a Magnetic has to migrate them, or that path throws on an old file even though
// the others do not. Each entry point is its own SECTION so a path that does not migrate is named.
// For each: the spacers come out as PET with no `material`, and stray capacitance -- which reads
// the spacer dielectric through the gap classification -- runs end to end.
//
//   pre-#1170 form   {"type": "spacer", "material": "plastic"}   the fixtures as they are on disk
//   #1170-era forms  {"type": "spacer", "insulationMaterial": "plastic"}, as a bare name and as an
//                    inline record {"name": "plastic"} -- built from the same document in memory,
//                    so no fixture is edited or added
namespace {

json abt1200_load_document(const std::string& fixture) {
    std::ifstream file(get_test_data_path(std::source_location::current(), fixture));
    REQUIRE(file.good());
    json document;
    file >> document;
    return document;
}

// Rewrites every spacer of a (full-MAS or bare-magnetic) document into the #1170-era form, as a
// bare name string or (inlineRecord) as an inline record carrying only the name.
json abt1200_as_1170_era(json document, bool inlineRecord = false) {
    json& core = document.contains("magnetic") ? document["magnetic"]["core"] : document["core"];
    size_t rewritten = 0;
    for (auto& element : core["geometricalDescription"]) {
        if (element.value("type", "") == "spacer") {
            element.erase("material");
            element["insulationMaterial"] = inlineRecord ? json{{"name", "plastic"}} : json("plastic");
            rewritten += 1;
        }
    }
    REQUIRE(rewritten > 0);
    return document;
}

void abt1200_check_spacers_are_pet(const OpenMagnetics::Core& core) {
    REQUIRE(core.get_geometrical_description());
    auto elements = core.get_geometrical_description().value();
    size_t spacers = 0;
    for (const auto& element : elements) {
        if (element.get_type() != MAS::CoreGeometricalDescriptionElementType::SPACER) {
            continue;
        }
        spacers += 1;
        REQUIRE(element.get_insulation_material());
        auto materialUnion = element.get_insulation_material().value();
        REQUIRE(std::holds_alternative<std::string>(materialUnion));
        CHECK(std::get<std::string>(materialUnion) == OpenMagnetics::Defaults().defaultSpacerMaterial);
        CHECK_FALSE(element.get_material());
    }
    REQUIRE(spacers > 0);
}

void abt1200_check_loads_and_computes(OpenMagnetics::Magnetic magnetic) {
    abt1200_check_spacers_are_pet(magnetic.get_core());
    auto coil = magnetic.get_coil();
    auto core = magnetic.get_core();
    MAS::StrayCapacitanceOutput output;
    REQUIRE_NOTHROW(output = OpenMagnetics::StrayCapacitance().calculate_capacitance(coil, core));
    REQUIRE(output.get_capacitance_among_windings());
}

}  // namespace

TEST_CASE("Test_MasMigration_Old_Plastic_Spacer_Designs_Load_And_Compute_On_Every_Entry_Point",
          "[constructive-model][core][mas-migration][abt1200][stray-capacitance]") {
    for (std::string fixture : {std::string("custom_magnetic_flyback.json"), std::string("huge_losses.json")}) {
        for (int form = 0; form < 3; ++form) {
            const bool era1170 = form > 0;
            DYNAMIC_SECTION(fixture << (form == 0 ? " as the pre-#1170 material form"
                                        : form == 1 ? " as the #1170-era insulationMaterial name"
                                                    : " as the #1170-era insulationMaterial inline record")) {
                json document = abt1200_load_document(fixture);
                if (era1170) {
                    document = abt1200_as_1170_era(document, form == 2);
                }
                const json magneticJson = document.at("magnetic");
                const json coreJson = magneticJson.at("core");

                SECTION("OpenMagnetics::from_file(path, Mas)") {
                    auto path = std::filesystem::temp_directory_path() / ("abt1200_" + std::to_string(form) + "_" + fixture);
                    {
                        std::ofstream out(path);
                        out << document.dump();
                    }
                    OpenMagnetics::Mas mas;
                    OpenMagnetics::from_file(path, mas);
                    std::filesystem::remove(path);
                    abt1200_check_loads_and_computes(mas.get_magnetic());
                }
                SECTION("OpenMagneticsTesting::mas_loader(path)") {
                    auto path = std::filesystem::temp_directory_path() / ("abt1200_loader_" + std::to_string(form) + "_" + fixture);
                    {
                        std::ofstream out(path);
                        out << document.dump();
                    }
                    auto mas = OpenMagneticsTesting::mas_loader(path);
                    std::filesystem::remove(path);
                    abt1200_check_loads_and_computes(mas.get_magnetic());
                }
                SECTION("from_json(json, Mas)") {
                    OpenMagnetics::Mas mas;
                    OpenMagnetics::from_json(document, mas);
                    abt1200_check_loads_and_computes(mas.get_magnetic());
                }
                SECTION("from_json(json, Magnetic)") {
                    OpenMagnetics::Magnetic magnetic;
                    OpenMagnetics::from_json(magneticJson, magnetic);
                    abt1200_check_loads_and_computes(magnetic);
                }
                SECTION("from_json(json, std::vector<Magnetic>)") {
                    std::vector<OpenMagnetics::Magnetic> magnetics;
                    OpenMagnetics::from_json(json::array({magneticJson}), magnetics);
                    REQUIRE(magnetics.size() == 1);
                    abt1200_check_loads_and_computes(magnetics[0]);
                }
                SECTION("OpenMagnetics::Magnetic magnetic(json) -- the PyOpenMagnetics / WASM binding form") {
                    OpenMagnetics::Magnetic magnetic(magneticJson);
                    abt1200_check_loads_and_computes(magnetic);
                }
                SECTION("magnetic_autocomplete(Magnetic(json)) -- the autocomplete flow") {
                    OpenMagnetics::Magnetic magnetic(magneticJson);
                    auto completed = OpenMagnetics::magnetic_autocomplete(magnetic);
                    abt1200_check_loads_and_computes(completed);
                }
                SECTION("OpenMagnetics::Core(json)") {
                    OpenMagnetics::Core core(coreJson);
                    abt1200_check_spacers_are_pet(core);
                }
                SECTION("json.get<OpenMagnetics::Core>()") {
                    OpenMagnetics::Core core = coreJson.get<OpenMagnetics::Core>();
                    abt1200_check_spacers_are_pet(core);
                }
            }
        }
    }
}
