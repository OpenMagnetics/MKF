#include "support/LibraryContext.h"
#include "support/Utils.h"
#include "advisers/MagneticAdviser.h"
#include "advisers/WireAdviser.h"
#include "processors/Inputs.h"
#include "json.hpp"

#include <algorithm>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

namespace {

TEST_CASE("LibraryContext: TypeFilterSet predicate", "[library-context][smoke-test]") {
    TypeFilterSet f;
    REQUIRE(f.empty());
    REQUIRE(f.accepts("anything"));

    f.allowed = {"E", "ETD"};
    REQUIRE(f.accepts("E"));
    REQUIRE(f.accepts("ETD"));
    REQUIRE_FALSE(f.accepts("RM"));

    f.blocked = {"ETD"};
    REQUIRE(f.accepts("E"));
    REQUIRE_FALSE(f.accepts("ETD"));
}

TEST_CASE("LibraryContext: enum predicates case-insensitive", "[library-context][smoke-test]") {
    TypeFilterSet f;
    f.allowed = {"etd"};
    REQUIRE(acceptsCoreShapeFamily(f, CoreShapeFamily::ETD));
    REQUIRE_FALSE(acceptsCoreShapeFamily(f, CoreShapeFamily::E));

    TypeFilterSet w;
    w.blocked = {"LITZ"};
    REQUIRE(acceptsWireType(w, WireType::ROUND));
    REQUIRE_FALSE(acceptsWireType(w, WireType::LITZ));

    TypeFilterSet m;
    m.allowed = {"FERRITE"};
    REQUIRE(acceptsCoreMaterialType(m, MaterialType::FERRITE));
    REQUIRE_FALSE(acceptsCoreMaterialType(m, MaterialType::POWDER));
}

TEST_CASE("LibraryContext: bad JSON throws", "[library-context][smoke-test]") {
    LibraryContext ctx;
    REQUIRE_THROWS(ctx.loadFromString("{ this is not json", LibraryContext::LoadMode::Merge));
    REQUIRE_THROWS(ctx.loadFromString(R"({"cores":[{"name":"x","bogus":1}]})", LibraryContext::LoadMode::Replace));
}

TEST_CASE("LibraryContext: empty context Scope is a no-op", "[library-context][smoke-test]") {
    if (wireDatabase.empty()) load_wires();
    auto before = wireDatabase.size();
    {
        LibraryContext ctx;
        auto scope = ctx.applyScoped();
        REQUIRE(wireDatabase.size() == before);
    }
    REQUIRE(wireDatabase.size() == before);
}

TEST_CASE("LibraryContext: Replace mode replaces the WHOLE library", "[library-context][smoke-test]") {
    // Replace previously wiped only the kinds the context provided. That left
    // dangling cross-references: with an inventory of shapes+materials+wires
    // (no cores), the full core catalog kept resolving shape NAMES against a
    // replaced shape database and the advisers threw CORE_SHAPE_NOT_FOUND
    // ('E 100/60/28') — reproduced live on openmagnetics.com (ABT). Replace
    // now means: the context IS the library; unprovided kinds become empty.
    if (wireDatabase.empty()) load_wires();
    if (coreShapeDatabase.empty()) load_core_shapes();
    if (coreDatabase.empty()) load_cores();
    auto coreShapeCountBefore = coreShapeDatabase.size();
    auto wireCountBefore = wireDatabase.size();
    auto coreCountBefore = coreDatabase.size();
    REQUIRE(coreCountBefore > 0);

    LibraryContext ctx;
    ctx.loadFromString(R"({"wires":{}})", LibraryContext::LoadMode::Replace);
    {
        auto scope = ctx.applyScoped();
        REQUIRE(wireDatabase.empty());
        REQUIRE(coreShapeDatabase.empty());
        REQUIRE(coreDatabase.empty());
    }
    REQUIRE(wireDatabase.size() == wireCountBefore);
    REQUIRE(coreShapeDatabase.size() == coreShapeCountBefore);
    REQUIRE(coreDatabase.size() == coreCountBefore);
}

TEST_CASE("LibraryContext: filterCoresByConstraints respects shape family", "[library-context][smoke-test]") {
    if (coreDatabase.empty()) load_cores();
    if (coreShapeDatabase.empty()) load_core_shapes();

    AdviserConstraints c;
    c.shapeFamily.allowed = {"ETD"};
    auto out = filterCoresByConstraints(coreDatabase, c);
    REQUIRE(out.size() < coreDatabase.size());
    for (const auto& core : out) {
        const auto& s = core.get_functional_description().get_shape();
        CoreShapeFamily fam;
        if (std::holds_alternative<MAS::CoreShape>(s)) {
            fam = std::get<MAS::CoreShape>(s).get_family();
        } else {
            fam = coreShapeDatabase.at(std::get<std::string>(s)).get_family();
        }
        REQUIRE(fam == CoreShapeFamily::ETD);
    }
}

TEST_CASE("LibraryContext: filterWiresByConstraints respects wire type", "[library-context][smoke-test]") {
    if (wireDatabase.empty()) load_wires();
    AdviserConstraints c;
    c.wireType.allowed = {"litz"};  // case-insensitive
    auto out = filterWiresByConstraints(wireDatabase, c);
    REQUIRE(out.size() > 0);
    REQUIRE(out.size() < wireDatabase.size());
    for (const auto& [name, wire] : out) {
        REQUIRE(wire.get_type() == WireType::LITZ);
    }
}

} // namespace

TEST_CASE("LibraryContext: Replace ctx adviser returns results from the context only", "[library-context][adviser]") {
    // End-to-end regression for the live 'only my inventory' failure: a
    // shapes+materials+wires-only Replace context must let the adviser run
    // (no CORE_SHAPE_NOT_FOUND from dangling catalog cores) and every result
    // must come from the context. Also asserts ctx=nullptr equivalence and
    // catalog restoration.
    if (coreShapeDatabase.empty()) load_core_shapes();
    if (coreMaterialDatabase.empty()) load_core_materials();
    if (wireDatabase.empty()) load_wires();

    json ctxJson;
    ctxJson["coreShapes"] = json::object();
    ctxJson["coreMaterials"] = json::object();
    ctxJson["wires"] = json::object();
    std::vector<std::string> shapeNames = {"PQ 26/20", "PQ 26/25", "PQ 32/20", "PQ 32/30", "PQ 40/40"};
    std::vector<std::string> materialNames = {"N87", "N97", "3C95"};
    size_t wireCount = 0;
    for (auto& name : shapeNames) {
        INFO("shape: " << name);
        REQUIRE(coreShapeDatabase.count(name) == 1);
        json shapeJson;
        to_json(shapeJson, coreShapeDatabase.at(name));
        ctxJson["coreShapes"][name] = shapeJson;
    }
    for (auto& name : materialNames) {
        INFO("material: " << name);
        REQUIRE(coreMaterialDatabase.count(name) == 1);
        json materialJson;
        to_json(materialJson, coreMaterialDatabase.at(name));
        ctxJson["coreMaterials"][name] = materialJson;
    }
    std::vector<std::string> wireNames = {
        "Round 0.2 - Grade 1", "Round 0.3 - Grade 1", "Round 0.5 - Grade 1",
        "Round 0.80 - Grade 1", "Round 1.00 - Grade 1", "Round 0.2 - Grade 2",
        "Round 0.3 - Grade 2", "Round 0.5 - Grade 2", "Round 0.80 - Grade 2",
        "Round 1.00 - Grade 2"};
    for (auto& name : wireNames) {
        INFO("wire: " << name);
        REQUIRE(wireDatabase.count(name) == 1);
        json wireJson;
        to_json(wireJson, wireDatabase.at(name));
        ctxJson["wires"][name] = wireJson;
        ++wireCount;
    }

    LibraryContext ctx;
    ctx.loadFromString(ctxJson.dump(), LibraryContext::LoadMode::Replace);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(
        100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 10, 0.5, 5);

    std::map<OpenMagnetics::MagneticFilters, double> weights;
    weights[OpenMagnetics::MagneticFilters::COST] = 0.33;
    weights[OpenMagnetics::MagneticFilters::LOSSES] = 0.33;
    weights[OpenMagnetics::MagneticFilters::DIMENSIONS] = 0.33;

    // Match the web app's adviser setup: candidates are GENERATED from
    // shapes x materials, so the stock filter (default on) must be off, and
    // concentric shapes need the concentric path enabled.
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.set_use_only_cores_in_stock(false);
    settings.set_use_toroidal_cores(false);
    settings.set_use_concentric_cores(true);

    // Control: the same inputs against the FULL catalog must find magnetics,
    // otherwise the environment (settings/inputs), not the ctx, is broken.
    OpenMagnetics::MagneticAdviser adviserBaseline;
    adviserBaseline.set_core_mode(OpenMagnetics::CoreAdviser::CoreAdviserModes::STANDARD_CORES);
    auto baselineResults = adviserBaseline.get_advised_magnetic(inputs, weights, 2);
    CAPTURE(baselineResults.size());
    REQUIRE(baselineResults.size() > 0);

    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(OpenMagnetics::CoreAdviser::CoreAdviserModes::STANDARD_CORES);
    // Probe: inside the scope, a Core built from a round-tripped ctx shape
    // must process like a catalog one (this is exactly what the dataset
    // builder does per candidate).
    {
        auto scope = ctx.applyScoped();
        REQUIRE(coreShapeDatabase.size() == shapeNames.size());
        OpenMagnetics::Core probeCore(coreShapeDatabase.begin()->second);
        probeCore.process_data();
        REQUIRE(probeCore.process_gap());
    }

    OpenMagnetics::read_log();  // drain, so the dump below covers this run only
    auto ctxResults = adviser.get_advised_magnetic(inputs, weights, 2, &ctx, OpenMagnetics::AdviserConstraints{});
    UNSCOPED_INFO("MKF log:\n" << OpenMagnetics::read_log());
    REQUIRE(ctxResults.size() > 0);
    for (auto& [mas, score] : ctxResults) {
        auto shapeUnion = mas.get_magnetic().get_core().get_functional_description().get_shape();
        std::string shapeName = std::holds_alternative<CoreShape>(shapeUnion)
            ? std::get<CoreShape>(shapeUnion).get_name().value()
            : std::get<std::string>(shapeUnion);
        CAPTURE(shapeName);
        REQUIRE(std::find(shapeNames.begin(), shapeNames.end(), shapeName) != shapeNames.end());
    }

    // ctx=nullptr must be equivalent to the base overload (both non-empty).
    OpenMagnetics::MagneticAdviser adviserBase;
    auto baseResults = adviserBase.get_advised_magnetic(inputs, weights, 2);
    OpenMagnetics::MagneticAdviser adviserNull;
    auto nullCtxResults = adviserNull.get_advised_magnetic(inputs, weights, 2, nullptr, OpenMagnetics::AdviserConstraints{});
    REQUIRE(baseResults.size() > 0);
    REQUIRE(nullCtxResults.size() == baseResults.size());

    // The public catalog must be intact after the scoped calls.
    REQUIRE(coreShapeDatabase.size() > 100);
}

TEST_CASE("LibraryContext: wire advising inside a scope never lazily reloads the public catalog", "[library-context][adviser]") {
    // A Replace context WITHOUT wires: inside its scope the wire pool is
    // empty and must STAY empty. The lazy load_wires() fallbacks in
    // WireAdviser/CoilAdviser would otherwise silently refill the shared
    // catalog and un-restrict 'only my inventory' wire advising (same
    // defect class as the gated load_cores() in CoreAdviser). Without the
    // gate this test FAILS: the adviser returns wires from the full catalog.
    if (coreShapeDatabase.empty()) load_core_shapes();
    if (wireDatabase.empty()) load_wires();

    json ctxJson;
    ctxJson["coreShapes"] = json::object();
    REQUIRE(coreShapeDatabase.count("PQ 26/20") == 1);
    json shapeJson;
    to_json(shapeJson, coreShapeDatabase.at("PQ 26/20"));
    ctxJson["coreShapes"]["PQ 26/20"] = shapeJson;

    LibraryContext ctx;
    ctx.loadFromString(ctxJson.dump(), LibraryContext::LoadMode::Replace);

    OpenMagnetics::Winding winding;
    winding.set_isolation_side(IsolationSide::PRIMARY);
    winding.set_name("primary");
    winding.set_number_parallels(1);
    winding.set_number_turns(10);
    winding.set_wire("Dummy");

    Section section;
    section.set_dimensions({0.005, 0.015});
    section.set_coordinate_system(CoordinateSystem::CARTESIAN);

    SignalDescriptor current;
    ProcessedWaveform processed;
    processed.set_peak_to_peak(2 * 2 * 1.4142);
    processed.set_rms(2);
    processed.set_effective_frequency(13456);
    processed.set_offset(0);
    processed.set_label(MAS::WaveformLabel::TRIANGULAR);
    current.set_processed(processed);

    OpenMagnetics::WireAdviser wireAdviser;
    auto results = wireAdviser.get_advised_wire(winding, section, current, 25, 1, 5,
                                                &ctx, OpenMagnetics::AdviserConstraints{});
    REQUIRE(results.empty());

    // The public catalog is restored once the scoped call returns.
    REQUIRE(!wireDatabase.empty());
}
