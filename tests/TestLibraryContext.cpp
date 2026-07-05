#include "support/LibraryContext.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

namespace {

TEST_CASE("LibraryContext: TypeFilterSet predicate", "[library-context]") {
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

TEST_CASE("LibraryContext: enum predicates case-insensitive", "[library-context]") {
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

TEST_CASE("LibraryContext: bad JSON throws", "[library-context]") {
    LibraryContext ctx;
    REQUIRE_THROWS(ctx.loadFromString("{ this is not json", LibraryContext::LoadMode::Merge));
    REQUIRE_THROWS(ctx.loadFromString(R"({"cores":[{"name":"x","bogus":1}]})", LibraryContext::LoadMode::Replace));
}

TEST_CASE("LibraryContext: empty context Scope is a no-op", "[library-context]") {
    if (wireDatabase.empty()) load_wires();
    auto before = wireDatabase.size();
    {
        LibraryContext ctx;
        auto scope = ctx.applyScoped();
        REQUIRE(wireDatabase.size() == before);
    }
    REQUIRE(wireDatabase.size() == before);
}

TEST_CASE("LibraryContext: Replace mode wipes only provided kinds", "[library-context]") {
    if (wireDatabase.empty()) load_wires();
    if (coreShapeDatabase.empty()) load_core_shapes();
    auto coreShapeCountBefore = coreShapeDatabase.size();
    auto wireCountBefore = wireDatabase.size();

    LibraryContext ctx;
    // Provide an empty wires map in Replace mode -> wireDatabase becomes empty
    // inside the scope, but coreShapeDatabase is untouched.
    ctx.loadFromString(R"({"wires":{}})", LibraryContext::LoadMode::Replace);
    {
        auto scope = ctx.applyScoped();
        REQUIRE(wireDatabase.empty());
        REQUIRE(coreShapeDatabase.size() == coreShapeCountBefore);
    }
    REQUIRE(wireDatabase.size() == wireCountBefore);
}

TEST_CASE("LibraryContext: filterCoresByConstraints respects shape family", "[library-context]") {
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

TEST_CASE("LibraryContext: filterWiresByConstraints respects wire type", "[library-context]") {
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
