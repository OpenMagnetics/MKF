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
#include "constructive_models/Core.h"
#include "constructive_models/CorePiece.h"
#include "support/Utils.h"
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

    // MAS carries these records; MKF has no CorePiece for their families yet
    // (ABT #274 / #275). They must not reach the database. When a family IS
    // implemented, its shapes load again and this expectation flips — update the
    // list rather than deleting the guard.
    REQUIRE_FALSE(CorePiece::is_family_supported(CoreShapeFamily::UI));
    REQUIRE_FALSE(CorePiece::is_family_supported(CoreShapeFamily::PQI));
    for (auto& name : {std::string("UI 93/76/20"), std::string("PQI 16/7.8")}) {
        CHECK_THROWS_AS(find_core_shape_by_name(name), CoreShapeNotFoundException);
    }
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
