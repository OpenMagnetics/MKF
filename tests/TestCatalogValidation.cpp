// ABT #307: the shape catalog must never hand an adviser something it cannot build.
//
// Two different problems get two different answers at load time:
//   * a family whose geometry class does not exist yet (MAS ships UI/PQI records
//     ahead of MKF) is a KNOWN GAP — those shapes are left out of the database and
//     reported once, so a sweep over the catalog is unaffected;
//   * geometry that cannot physically exist (a toroid whose inner diameter is at
//     least as large as its outer one, a non-positive principal dimension) is
//     CORRUPT DATA — it is rejected outright, because a record like that used to
//     kill entire adviser runs with a message that named nothing useful
//     ("IEC 63182 effective parameters cannot be negative or 0", ABT #306).
#include <source_location>
#include <filesystem>
#include <fstream>
#include <magic_enum.hpp>
#include "constructive_models/Core.h"
#include "constructive_models/CorePiece.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "advisers/MagneticAdviser.h"
#include "processors/Inputs.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/MagneticSimulator.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

TEST_CASE("Test_Catalog_Every_Loaded_Shape_Is_Buildable", "[catalog][smoke-test]") {
    settings.reset();
    clear_databases();

    auto names = get_core_shape_names();
    REQUIRE(names.size() > 1000);

    std::vector<std::string> unbuildable;
    for (auto& name : names) {
        try {
            auto piece = CorePiece::factory(find_core_shape_by_name(name), true);
            auto effective = piece->get_partial_effective_parameters();
            if (!(effective.get_effective_area() > 0) || !(effective.get_effective_length() > 0)) {
                unbuildable.push_back(name + " (non-positive effective parameters)");
            }
        }
        catch (const std::exception& e) {
            unbuildable.push_back(name + " (" + e.what() + ")");
        }
    }

    INFO("shapes that loaded but cannot be built: " << [&] {
        std::string joined;
        for (auto& entry : unbuildable) {
            joined += "\n  " + entry;
        }
        return joined;
    }());
    CHECK(unbuildable.empty());
}

TEST_CASE("Test_Catalog_Unsupported_Families_Are_Not_Loaded", "[catalog][smoke-test]") {
    settings.reset();
    clear_databases();

    // MAS carries records whose families MKF has no CorePiece for; they must not reach the
    // database. The guard is kept and the LIST is updated as families are implemented, per its
    // original instruction.
    //
    // UI and PQI moved OFF this list (ABT #274 / #275, user-approved): both now have geometry
    // classes and load normally, so the expectation flips to the assertions below.
    // Declared in CoreShapeFamily but with no CorePiece geometry: these must stay unsupported
    // so load_core_shapes keeps skipping them instead of half-building a core.
    // DRUM moved off this list (ABT #331): it has a geometry class and an open-core model.
    // EI moved off this list (ABT #625): it has a geometry class now too, see below.
    // ROD moved off this list (ABT #933): CorePieceRod + the rod open-core model. Note that MAS
    // still carries ZERO rod records, so unlike UI/DRUM there is no catalogued name to resolve
    // below — a rod reaches MKF as an inline shape today.
    for (auto family : {CoreShapeFamily::BLOCK, CoreShapeFamily::H}) {
        REQUIRE_FALSE(CorePiece::is_family_supported(family));
    }
    REQUIRE(CorePiece::is_family_supported(CoreShapeFamily::ROD));

    // UI is now supported (ABT #274, user-approved) and its shapes must resolve.
    REQUIRE(CorePiece::is_family_supported(CoreShapeFamily::UI));
    CHECK_NOTHROW(find_core_shape_by_name(std::string("UI 93/76/20")));

    // DRUM is supported as an OPEN shape (ABT #331): its records must load.
    REQUIRE(CorePiece::is_family_supported(CoreShapeFamily::DRUM));
    CHECK_NOTHROW(find_core_shape_by_name(std::string("DRH-14X20-4C")));
    CHECK_NOTHROW(find_core_shape_by_name(std::string("Bobbin 9643001015")));

    // DRUM_RING (shielded drum, ABT #366) is supported as a piece-and-plate closed circuit:
    // the ACME DR + SRI pair records must load.
    REQUIRE(CorePiece::is_family_supported(CoreShapeFamily::DRUM_RING));
    CHECK_NOTHROW(find_core_shape_by_name(std::string("DR 2.3 + SRI 3.0")));
    CHECK_NOTHROW(find_core_shape_by_name(std::string("DR 2.0 + SRI 2.95")));

    // DRUM_SEMISHIELDED (ABT #362) and MOLDED (ABT #357) have piece classes too; they carry no
    // catalogue shapes yet (constructions are reconstructed per part, not sold as bare cores),
    // so support is asserted at the piece level only.
    REQUIRE(CorePiece::is_family_supported(CoreShapeFamily::DRUM_SEMISHIELDED));
    REQUIRE(CorePiece::is_family_supported(CoreShapeFamily::MOLDED));

    // PQI is supported too (ABT #275): the PQ clause of IEC 60205 covers the plate case, and the
    // geometry validates against TDK's published planar data to 0.3%.
    REQUIRE(CorePiece::is_family_supported(CoreShapeFamily::PQI));
    CHECK_NOTHROW(find_core_shape_by_name(std::string("PQI 16/7.8")));

    // EI is supported too (ABT #625): an E piece closed by a flat I plate rather than a mirrored
    // second E-half. No catalogue entries ship yet (see Test_ABT625_Ei_Matches_Published_Data for
    // the geometry validation against a real vendor datasheet), so only the family-level guard is
    // checked here -- there is no find_core_shape_by_name(...) to exercise.
    REQUIRE(CorePiece::is_family_supported(CoreShapeFamily::EI));
}

// ABT #625: MAS lists "ei" as a legal CoreShapeFamily (an E piece closed by a flat I plate, same
// pattern as UI/PQI above) but MKF's CorePiece factory did not implement it, throwing "Unknown shape
// family: ei" for every EI core. Added CorePieceEi (CorePiece.cpp), validated here against the
// Magnetics 2022 Ferrite Catalog (p.32-35, planar E/I dimensions and magnetic data), which publishes
// the E-piece-ALONE (mirrored-pair) le/Ae under the E-code row and the E+I ASSEMBLY's le/Ae under the
// mated I-code row.
TEST_CASE("Test_ABT625_Ei_Matches_Published_Data", "[constructive-model][core][bug]") {
    settings.reset();
    clear_databases();

    // E 40/8/10 (0_44008EC) alone: CorePieceE's own (single-piece) formula doubled for the
    // TWO_PIECE_SET mirrored pair matches the E-code row's published le=51.9mm / Ae=101mm2.
    {
        json eShape = json::parse(R"({"family":"e","type":"custom","name":"E 40/8/10 test",
            "dimensions":{"A":{"nominal":0.04065},"B":{"nominal":0.00851},"C":{"nominal":0.0107},
                          "D":{"nominal":0.00406},"E":{"nominal":0.03045},"F":{"nominal":0.01015}}})");
        auto piece = CorePiece::factory(eShape, true);
        auto [le, ae, minArea] = piece->get_shape_constants_iec63182();
        CHECK_THAT(le * 2000, Catch::Matchers::WithinRel(51.9, 0.02));
        CHECK_THAT(ae * 1e6, Catch::Matchers::WithinRel(101.0, 0.01));
    }
    // I 40/4/10 (0_44008IC): the same E piece mated with a B2=4.45mm plate (symmetric: equal to the
    // piece's own yoke thickness B-D=4.45mm). Published le=43.8mm / Ae=99.5mm2.
    {
        json eiShape = json::parse(R"({"family":"ei","type":"custom","name":"EI 40/8/10 test",
            "dimensions":{"A":{"nominal":0.04065},"B":{"nominal":0.00851},"B2":{"nominal":0.00445},
                          "C":{"nominal":0.0107},"D":{"nominal":0.00406},"E":{"nominal":0.03045},
                          "F":{"nominal":0.01015}}})");
        auto piece = CorePiece::factory(eiShape, true);
        auto [le, ae, minArea] = piece->get_shape_constants_iec63182();
        CHECK_THAT(le * 1000, Catch::Matchers::WithinRel(43.8, 0.02));
        CHECK_THAT(ae * 1e6, Catch::Matchers::WithinRel(99.5, 0.01));
    }
    // I 43/4/28 (0_44308IC): B2=4.1mm vs the piece's own yoke B-D=4.32mm -- an ASYMMETRIC
    // plate/yoke case. Published le=48.6mm / Ae=227mm2.
    {
        json eiShape2 = json::parse(R"({"family":"ei","type":"custom","name":"EI 43/8/28 test",
            "dimensions":{"A":{"nominal":0.0432},"B":{"nominal":0.00851},"B2":{"nominal":0.0041},
                          "C":{"nominal":0.0279},"D":{"nominal":0.00419},"E":{"nominal":0.0344},
                          "F":{"nominal":0.00813}}})");
        auto piece = CorePiece::factory(eiShape2, true);
        auto [le, ae, minArea] = piece->get_shape_constants_iec63182();
        CHECK_THAT(le * 1000, Catch::Matchers::WithinRel(48.6, 0.02));
        CHECK_THAT(ae * 1e6, Catch::Matchers::WithinRel(227.0, 0.04));
    }
    settings.reset();
}

// Regression for the ticket's exact repro: a pieceAndPlate core built from a custom "ei" shape must
// process end to end (not just at the CorePiece::factory level) and produce sane effective parameters.
TEST_CASE("Test_ABT625_PieceAndPlate_Core_Processes", "[constructive-model][core][bug]") {
    settings.reset();
    clear_databases();

    json coreJson = json::parse(R"({"name":"test","functionalDescription":{"type":"pieceAndPlate",
        "material":"3C95","gapping":[],"numberStacks":1,
        "shape":{"type":"custom","family":"ei","name":"test ei",
            "dimensions":{"A":{"nominal":0.02},"B":{"nominal":0.01},"B2":{"nominal":0.003},
                          "C":{"nominal":0.008},"D":{"nominal":0.005},"E":{"nominal":0.015},
                          "F":{"nominal":0.006}}}}})");

    Core core(coreJson);
    REQUIRE(core.get_effective_area() > 0);
    REQUIRE(core.get_effective_length() > 0);
    settings.reset();
}

TEST_CASE("Test_Catalog_Impossible_Toroid_Is_Rejected", "[catalog][smoke-test]") {
    settings.reset();
    clear_databases();

    json corrupt = json::parse(R"({"family":"t","type":"standard","magneticCircuit":"closed",
        "name":"T TEST 32.5/78.6/20.3",
        "dimensions":{"A":{"nominal":0.0325},"B":{"nominal":0.0786},"C":{"nominal":0.0203}}})");
    CHECK_THROWS_AS(load_core_shapes(true, corrupt.dump()), InvalidInputException);
    clear_databases();

    json negative = json::parse(R"({"family":"e","type":"standard","magneticCircuit":"open",
        "name":"E TEST negative",
        "dimensions":{"A":{"nominal":-0.02},"B":{"nominal":0.01},"C":{"nominal":0.005}}})");
    CHECK_THROWS_AS(load_core_shapes(true, negative.dump()), InvalidInputException);
    clear_databases();

    // A signed auxiliary dimension is legal and must still load (EFD's K is an offset).
    json signedOffset = json::parse(R"({"family":"efd","type":"standard","magneticCircuit":"open",
        "name":"EFD TEST signed",
        "dimensions":{"A":{"nominal":0.0105},"B":{"nominal":0.0053},"C":{"nominal":0.003},
                      "K":{"nominal":-0.0002}}})");
    CHECK_NOTHROW(load_core_shapes(true, signedOffset.dump()));
    clear_databases();
}

// ABT #370: shielded-drum CORES must reach the database, or the advisers can only ever evaluate
// a user-built drumRing and never PROPOSE one. Both pairs are ACME B45 (same NiZn grade on drum
// and ring, so MKF's single-material sectioned circuit is exact for them).
//
// NOTE the stock gate: load_cores() honours Settings::use_only_cores_in_stock, which DEFAULTS TO
// TRUE and then reads cores_stock.ndjson — a distributor-stock claim these bare ACME drum+ring
// pairs are not part of. A catalogue adviser aimed at this family must opt out of the stock
// filter; fabricating stock rows to paper over that would be a lie about purchasability. Both
// paths are pinned below.
TEST_CASE("Test_Catalog_Shielded_Drum_Cores_Are_Available", "[catalog][drum-ring][smoke-test]") {
    settings.reset();
    clear_databases();

    CHECK_NOTHROW(find_core_material_by_name(std::string("B45")));

    settings.set_use_only_cores_in_stock(false);
    clear_loaded_cores();
    load_cores();
    size_t drumRingCoresInDatabase = 0;
    for (auto& core : coreDatabase) {
        if (core.get_shape_family() == CoreShapeFamily::DRUM_RING) {
            drumRingCoresInDatabase++;
        }
    }
    CHECK(drumRingCoresInDatabase >= 2);
    auto shieldedCore = find_core_by_name(std::string("DR 2.3 + SRI 3.0 - B45 - Shielded"));
    CHECK(shieldedCore.get_functional_description().get_type() == CoreType::PIECE_AND_PLATE);

    // Stock-only (the default): the pairs are absent, and that is the honest answer.
    settings.set_use_only_cores_in_stock(true);
    clear_loaded_cores();
    load_cores();
    size_t drumRingCoresInStock = 0;
    for (auto& core : coreDatabase) {
        if (core.get_shape_family() == CoreShapeFamily::DRUM_RING) {
            drumRingCoresInStock++;
        }
    }
    CHECK(drumRingCoresInStock == 0);
    settings.reset();
}

// ABT #366/#370: the adviser must actually RUN with a shielded-drum catalogue. This is the
// integration nobody had exercised, and there was a concrete reason to worry: CoreAdviser tries
// to solve a GROUND gap for every candidate and then calls process_gap() (CoreAdviserDataset.cpp
// :653-691), while a drumRing REJECTS all user gapping — nothing can be ground on an assembled
// drum+ring. Those call sites happen to be inside catch(std::exception), so the throw degrades
// to "candidate keeps its structural gaps" instead of killing the run; this test pins that
// contract end-to-end rather than trusting the reading. A catalogue advise over these parts must
// return without throwing, and every returned magnetic must still carry its two structural
// annular clearances (the adviser must not have re-gapped or un-gapped them).
TEST_CASE("Test_Adviser_Runs_Over_Shielded_Drum_Catalogue", "[catalog][drum-ring][adviser]") {
    settings.reset();
    clear_databases();
    settings.set_use_only_cores_in_stock(false);
    load_cores();

    std::vector<OpenMagnetics::Magnetic> catalogueMagnetics;
    for (auto& core : coreDatabase) {
        if (core.get_shape_family() != CoreShapeFamily::DRUM_RING) continue;
        json coilJson;
        coilJson["bobbin"] = "Dummy";
        coilJson["functionalDescription"] = json::array({{
            {"name", "winding 0"}, {"numberTurns", 12}, {"numberParallels", 1},
            {"isolationSide", "primary"}, {"wire", "Round 0.1 - Grade 1"}}});
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
        catalogueMagnetics.push_back(OpenMagnetics::magnetic_autocomplete(magnetic));
    }
    REQUIRE(catalogueMagnetics.size() >= 2);

    // Every candidate entered the adviser with its structural clearances intact.
    for (auto& magnetic : catalogueMagnetics) {
        REQUIRE(magnetic.get_core().get_functional_description().get_gapping().size() == 2);
    }

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        1000000, 14e-6, 25, WaveformLabel::TRIANGULAR, 0.2, 0.5, 0.3);

    OpenMagnetics::MagneticAdviser adviser;
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    // strict=false: these tiny parts need not satisfy every requirement — the point is that the
    // flow completes over a family whose gapping cannot be solved for.
    REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(inputs, catalogueMagnetics, 2, false));

    for (auto& [mas, score] : results) {
        auto gapping = mas.get_magnetic().get_core().get_functional_description().get_gapping();
        CHECK(gapping.size() == 2);
        for (auto& gap : gapping) {
            CHECK(gap.get_type() == GapType::RESIDUAL);
            CHECK(gap.get_length() > 0);
        }
        CHECK(std::isfinite(score));
    }
    settings.reset();
}

// ABT #357: advisers must be able to PROPOSE molded parts, not merely evaluate one handed to
// them. The 183 WE-MAPI reconstructions are FITTED (public L/DCR/dims + measured curves through
// scripts/fit_we_mapi_internals.py), NOT vendor internals, so they are deliberately kept OUT of
// the canonical MAS catalogue — shipping a reconstruction as a catalogue part would dress a fit
// up as manufacturer data. They are instead usable as a caller-supplied catalogue, which is what
// this pins: build magnetics from the reconstructions and let MagneticAdviser rank them.
TEST_CASE("Test_Adviser_Proposes_Molded_Reconstructions", "[catalog][molded][adviser]") {
    settings.reset();
    clear_databases();
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path().append("testData").append("we_mapi_reconstructions.json");
    std::ifstream reconstructionsFile(path);
    REQUIRE(reconstructionsFile.good());
    json reconstructions = json::parse(reconstructionsFile);
    REQUIRE(reconstructions.size() == 183);

    // The fitted material must be REGISTERED, not inlined: MKF's simulate stage resolves the
    // core material by NAME from the database (CORE_MATERIAL_NOT_FOUND otherwise), even though
    // the material object is already inside the core. Registering it is also what a real
    // consumer does with a fitted grade, so this is the honest setup rather than a workaround.
    // The inline-material limitation is filed separately.
    std::string fittedMaterialName = "WE metal alloy (fitted mu_eff0)";
    std::vector<OpenMagnetics::Magnetic> catalogueMagnetics;
    for (auto& reconstruction : reconstructions) {
        if (reconstruction.at("caseCode").get<std::string>() != "4020") continue;  // one family
        json materialJson = {
            {"name", "WE metal alloy (fitted mu_eff0)"}, {"type", "custom"},
            {"material", "powder"}, {"materialComposition", "proprietary"},
            {"manufacturerInfo", {{"name", "FITTED — not vendor material data"}}},
            {"permeability", {{"initial", {{"value", reconstruction.at("mu_eff0").get<double>()}}}}},
            {"resistivity", json::array({{{"value", 1.0}}})},
            {"density", 5500}, {"curieTemperature", 500},
            {"saturation", json::array({{{"magneticFluxDensity", 1.2}, {"magneticField", 40000}, {"temperature", 25}}})},
            {"volumetricLosses", {{"default", json::array({{{"method", "lossFactor"},
                {"factors", json::array({{{"value", 1e-5}, {"frequency", 100000}}})}}})}}}
        };
        json shapeJson = {
            {"magneticCircuit", "closed"}, {"type", "custom"}, {"family", "molded"},
            {"aliases", json::array()},
            {"name", "WE-MAPI " + reconstruction.at("orderCode").get<std::string>()},
            {"dimensions", {
                {"A", {{"nominal", reconstruction.at("A").get<double>()}}},
                {"B", {{"nominal", reconstruction.at("B").get<double>()}}},
                {"C", {{"nominal", reconstruction.at("C").get<double>()}}},
                {"D", {{"nominal", reconstruction.at("D").get<double>()}}},
                {"E", {{"nominal", reconstruction.at("E").get<double>()}}},
                {"F", {{"nominal", reconstruction.at("F").get<double>()}}}}}
        };
        materialJson["name"] = fittedMaterialName;
        // Load the built-in materials FIRST, then add the fitted one on top: passing content to
        // load_core_materials seeds the database from that content alone, so registering the
        // fitted grade into an empty database would leave the process without the catalogue
        // materials and break every later test that looks one up (it did — "Kool Mu 26 not
        // found" in the next test case). The database is cleared again at the end of this test.
        if (coreMaterialDatabase.empty()) {
            load_core_materials();
        }
        load_core_materials(materialJson.dump());
        json coreJson;
        coreJson["functionalDescription"] = {
            {"type", "closedShape"}, {"material", fittedMaterialName}, {"shape", shapeJson},
            {"gapping", json::array()}, {"numberStacks", 1}};
        OpenMagnetics::Core core(coreJson);
        core.process_data();
        core.process_gap();

        json coilJson;
        coilJson["bobbin"] = "Dummy";
        coilJson["functionalDescription"] = json::array({{
            {"name", "winding 0"},
            {"numberTurns", reconstruction.at("numberTurns").get<int>()},
            {"numberParallels", 1}, {"isolationSide", "primary"},
            {"wire", reconstruction.at("wire").get<std::string>()}}});
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(OpenMagnetics::Coil(coilJson, false));
        try {
            catalogueMagnetics.push_back(OpenMagnetics::magnetic_autocomplete(magnetic));
        }
        catch (const std::exception&) {
            // A reconstruction whose fitted cavity cannot hold its own fitted winding is a
            // known limitation of fitting from public data (the fit reports packing usage per
            // part); skip it rather than failing the adviser contract.
        }
    }
    REQUIRE(catalogueMagnetics.size() >= 5);

    // Ask for something these parts can actually serve: take the target inductance from a real
    // candidate and a current its fine reconstructed wire can carry. (With 1 A DC and a flat
    // 1 uH request, DC_CURRENT_DENSITY and MAGNETIZING_INDUCTANCE both reject every candidate —
    // correct behaviour, but it tests the requirements rather than the family.)
    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    double targetInductance = magnetizingInductanceModel
        .calculate_inductance_from_number_turns_and_gapping(catalogueMagnetics[0].get_core(),
                                                            catalogueMagnetics[0].get_coil())
        .get_magnetizing_inductance().get_nominal().value();
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        1000000, targetInductance, 25, WaveformLabel::TRIANGULAR, 0.05, 0.5, 0.1);

    // WHY THE FLOW IS TRIMMED HERE — a real limitation, not a workaround. Run the default
    // catalogue flow and every molded candidate is rejected, because the IMPEDANCE filter throws
    // [MATERIAL_DATA_MISSING] for a fitted material: reconstructing mu_eff0 from public L/DCR
    // data gives a permeability, not the complex-permeability curves that impedance needs. A
    // control candidate (E 25/13/7, catalogue material) passes the same flow, so this is specific
    // to fitted materials, not to the molded family or the flow itself. Consequence for ABT #357
    // phase 2: to be adviser-ready against impedance requirements, a fitted molded material would
    // need complex permeability measured or fitted too — inventing it here would be fabrication.
    // This design has no impedance requirement, so the flow is trimmed to the filters its
    // requirements actually imply.
    OpenMagnetics::MagneticAdviser adviser;
    std::vector<MagneticFilterOperation> filterFlow;
    for (auto& operation : adviser._defaultCatalogueMagneticFilterFlow) {
        if (operation.get_filter() == MagneticFilters::IMPEDANCE) continue;
        filterFlow.push_back(operation);
    }

    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(inputs, catalogueMagnetics, filterFlow, 3, false));
    UNSCOPED_INFO("candidates " << catalogueMagnetics.size() << " -> results " << results.size());
    CHECK(results.size() > 0);
    for (auto& [mas, score] : results) {
        CHECK(mas.get_magnetic().get_core().get_shape_family() == CoreShapeFamily::MOLDED);
        CHECK(std::isfinite(score));
    }

    // The default flow (with IMPEDANCE) is pinned as the known limitation, so the day fitted
    // materials carry complex permeability this test fails and gets updated deliberately.
    auto resultsWithImpedance = adviser.get_advised_magnetic(inputs, catalogueMagnetics, 3, false);
    CHECK(resultsWithImpedance.empty());

    // The fitted grade must not outlive this test: it was injected into the process-global
    // material database, and leaving it there changes what every later test resolves.
    clear_databases();

    settings.reset();
}

// The families a UI may offer are the ENGINE's, not the database's. get_core_shape_families()
// (Utils.h) reports only families with a loaded catalogue shape, so DRUM_SEMISHIELDED and MOLDED
// — real, buildable families that ship no bare-core records because the construction is
// reconstructed per part — were absent from the shape-family dropdown entirely, and a user had no
// way to reach them. get_supported_core_shape_families() answers the capability question instead.
TEST_CASE("Test_Supported_Families_Are_The_Engine_Not_The_Database", "[catalog][smoke-test]") {
    settings.reset();
    clear_databases();

    auto supported = OpenMagnetics::get_supported_core_shape_families();
    REQUIRE(supported.size() > 30);

    // Every family the list names must actually be buildable, and every buildable family must be
    // named: the two answers come from one array precisely so they cannot drift apart.
    for (auto family : supported) {
        INFO("listed as supported: " << magic_enum::enum_name(family));
        CHECK(CorePiece::is_family_supported(family));
    }
    for (auto family : magic_enum::enum_values<CoreShapeFamily>()) {
        if (!CorePiece::is_family_supported(family)) {
            continue;
        }
        INFO("buildable but missing from the list: " << magic_enum::enum_name(family));
        CHECK(std::find(supported.begin(), supported.end(), family) != supported.end());
    }

    // The point of the function: it must be a STRICT superset of what the database happens to
    // carry. If these ever became equal the function would be redundant and the families with no
    // catalogue shape would have silently vanished again.
    auto inDatabase = get_core_shape_families();
    REQUIRE_FALSE(inDatabase.empty());
    for (auto family : inDatabase) {
        INFO("in the shape database but not buildable: " << magic_enum::enum_name(family));
        CHECK(std::find(supported.begin(), supported.end(), family) != supported.end());
    }
    CHECK(supported.size() > inDatabase.size());

    // Named explicitly, because these two are the whole reason the distinction exists: buildable,
    // and carrying no catalogue shape of their own.
    for (auto family : {CoreShapeFamily::DRUM_SEMISHIELDED, CoreShapeFamily::MOLDED}) {
        CHECK(std::find(supported.begin(), supported.end(), family) != supported.end());
        CHECK(std::find(inDatabase.begin(), inDatabase.end(), family) == inDatabase.end());
    }
}

// ABT #1007: get_shape_family_dimensions answered "what does this family need?" by scanning the
// shape database for published shapes of that family. A family that ships no bare-core record
// therefore reported ZERO dimensions, and the builder — which draws one input per returned
// dimension — offered the family and then no fields at all, so a custom shape in it could not be
// defined. Eleven of the engine's families were in that state: molded, drumSemishielded, rod, ei,
// el, ef, epc, epq, ept, epw and lep.
TEST_CASE("Test_Family_Dimensions_Are_Declared_For_Every_Supported_Family", "[catalog][smoke-test]") {
    settings.reset();
    clear_databases();

    auto supported = OpenMagnetics::get_supported_core_shape_families();
    REQUIRE(supported.size() > 30);

    // The bug, stated directly: every family the engine can build must be describable.
    for (auto family : supported) {
        INFO("family " << magic_enum::enum_name(family));
        std::vector<std::string> required;
        REQUIRE_NOTHROW(required = OpenMagnetics::get_core_shape_family_required_dimensions(family));
        CHECK_FALSE(required.empty());

        std::vector<std::string> offered;
        REQUIRE_NOTHROW(offered = get_shape_family_dimensions(family));
        CHECK_FALSE(offered.empty());

        // Whatever the geometry requires must reach the caller.
        for (auto& dimension : required) {
            INFO("required dimension " << dimension << " missing from the offered set");
            CHECK(std::find(offered.begin(), offered.end(), dimension) != offered.end());
        }
    }

    // The families that carry no catalogue shape are the whole point: before this, each returned
    // an empty list. Their dimensions can now come only from the declaration.
    for (auto family : {CoreShapeFamily::MOLDED, CoreShapeFamily::DRUM_SEMISHIELDED,
                        CoreShapeFamily::ROD, CoreShapeFamily::EI, CoreShapeFamily::LEP}) {
        INFO("shape-less family " << magic_enum::enum_name(family));
        CHECK(get_shape_family_dimensions(family).size() >= 3);
    }

    // Spot-checks against the geometry, so a wrong declaration is caught and not just a missing
    // one. A rod is a cylinder: outer diameter, length, bore. A molded block is an E-like set.
    CHECK(get_shape_family_dimensions(CoreShapeFamily::ROD) ==
          std::vector<std::string>{"A", "B", "H"});
    CHECK(get_shape_family_dimensions(CoreShapeFamily::MOLDED) ==
          std::vector<std::string>{"A", "B", "C", "D", "E", "F"});
    // The shielded drums add the shield's own J/K/L envelope on top of the drum's dimensions.
    for (auto& dimension : {"J", "K", "L"}) {
        auto offered = get_shape_family_dimensions(CoreShapeFamily::DRUM_SEMISHIELDED);
        CHECK(std::find(offered.begin(), offered.end(), dimension) != offered.end());
    }
}

// Backward compatibility: no family that worked before may lose an input field. The catalogue
// union is still folded in, because some of those keys are optional geometry this engine reads
// behind a guard (drum's A2, a toroid's R/r0) and some it never reads while another consumer
// does — MVB++ renders EC's T and PM's alpha off the shape dimensions.
TEST_CASE("Test_Family_Dimensions_Still_Cover_Every_Published_Shape", "[catalog][smoke-test]") {
    settings.reset();
    clear_databases();

    auto names = get_core_shape_names();
    REQUIRE(names.size() > 1000);

    size_t checked = 0;
    for (auto& name : names) {
        auto shape = find_core_shape_by_name(name);
        if (!shape.get_dimensions()) {
            continue;
        }
        auto offered = get_shape_family_dimensions(shape.get_family());
        // Materialise the optional: get_dimensions() returns by value, so binding the range-for
        // to .value() directly is a dangling reference (-Werror=dangling-reference).
        auto shapeDimensions = shape.get_dimensions().value();
        for (auto& [dimensionKey, dimensionValue] : shapeDimensions) {
            INFO(name << " (" << magic_enum::enum_name(shape.get_family())
                      << ") carries dimension " << dimensionKey << ", which is no longer offered");
            CHECK(std::find(offered.begin(), offered.end(), dimensionKey) != offered.end());
        }
        checked++;
    }
    CHECK(checked > 1000);
}

// A moulded body may be wound on an OVAL post. Several WE moulded families state both axes of the
// bore ('1.5/1.0'), and read as a single diameter the pole area is wrong — which makes the winding
// window (E - F)/2 wrong with it, since the window is measured from the post. F is the smaller
// axis and F2 the larger, matching the nomenclature the other oblong-column families already use.
TEST_CASE("Test_Molded_Oval_Post_Uses_Both_Axes", "[catalog][constructive-model][smoke-test]") {
    settings.reset();
    clear_databases();

    auto shape = [](double F, std::optional<double> F2) {
        json j = {{"type", "custom"}, {"family", "molded"}, {"name", "MLD oval probe"},
                  {"magneticCircuit", "closed"}, {"aliases", json::array()},
                  {"dimensions", {{"A", 0.004}, {"B", 0.002}, {"C", 0.004},
                                  {"D", 0.0009}, {"E", 0.0030}, {"F", F}}}};
        if (F2) {
            j["dimensions"]["F2"] = *F2;
        }
        return CoreShape(j);
    };

    // Round: no F2 at all. The post is a circle of F and nothing about this changes.
    auto roundPiece = CorePiece::factory(shape(0.0010, std::nullopt), true);
    auto roundColumn = roundPiece->get_columns()[0];
    REQUIRE(roundColumn.get_shape() == ColumnShape::ROUND);
    CHECK(std::abs(roundColumn.get_area() - std::numbers::pi / 4 * 0.0010 * 0.0010) < 1e-9);   // roundFloat quantises areas at 1e-9

    // Oval: F 1.0 mm x F2 1.5 mm. A stadium is a rectangle with a half-disc on each end, so the
    // area is pi/4*F^2 + (F2 - F)*F — strictly MORE than the circle, which is the whole point:
    // reading only F understates the pole and therefore overstates the reluctance.
    auto ovalPiece = CorePiece::factory(shape(0.0010, 0.0015), true);
    auto ovalColumn = ovalPiece->get_columns()[0];
    REQUIRE(ovalColumn.get_shape() == ColumnShape::OBLONG);
    CHECK(std::abs(ovalColumn.get_width() - 0.0010) < 1e-9);
    CHECK(std::abs(ovalColumn.get_depth() - 0.0015) < 1e-9);
    CHECK(std::abs(ovalColumn.get_area() -
                   (std::numbers::pi / 4 * 0.0010 * 0.0010 + (0.0015 - 0.0010) * 0.0010)) < 1e-9);   // roundFloat quantises areas at 1e-9
    CHECK(ovalColumn.get_area() > roundColumn.get_area());

    // A degenerate or inverted F2 must NOT silently produce an oblong post: F2 == F is a circle,
    // and F2 < F would mean the "larger" axis is smaller, which is a data error, not a shape.
    CHECK(CorePiece::factory(shape(0.0010, 0.0010), true)->get_columns()[0].get_shape()
          == ColumnShape::ROUND);
    CHECK(CorePiece::factory(shape(0.0010, 0.0008), true)->get_columns()[0].get_shape()
          == ColumnShape::ROUND);
}
