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
    for (auto family : {CoreShapeFamily::ROD, CoreShapeFamily::BLOCK,
                        CoreShapeFamily::EI, CoreShapeFamily::H}) {
        REQUIRE_FALSE(CorePiece::is_family_supported(family));
    }

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
