// Validation for the 2026-07 shape families added to CorePiece::factory
// (ABT #263-#272). Each new family is modelled as inheriting an existing piece's
// geometry; this test confirms (a) the factory recognises the family without
// throwing, (b) effective parameters compute, and (c) for identical A-F dimensions
// the child family yields the SAME Ae/Le as its parent piece (i.e. the inheritance
// is wired correctly).
#include "constructive_models/Core.h"
#include "TestingUtils.h"
#include "json.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;

static std::pair<double, double> ae_le(CoreShape shape, CoreShapeFamily family) {
    shape.set_family(family);
    auto piece = CorePiece::factory(shape, true);
    auto ep = piece->get_partial_effective_parameters();
    return {ep.get_effective_area(), ep.get_effective_length()};
}

TEST_CASE("New shape families inherit parent geometry", "[corepiece][newfamilies]") {
    // {parent shape name, parent family, child family}
    struct Case { std::string shapeName; CoreShapeFamily parent; CoreShapeFamily child; };
    std::vector<Case> cases = {
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::EPC},
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::EPQ},
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::EPW},
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::EPT},
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::LEP},
        {"ETD 29/16/10", CoreShapeFamily::ETD, CoreShapeFamily::EER},
        {"E 20/10/6",    CoreShapeFamily::E,   CoreShapeFamily::EF},
    };
    for (auto& c : cases) {
        auto shape = find_core_shape_by_name(c.shapeName);
        auto [aeP, leP] = ae_le(shape, c.parent);
        auto [aeC, leC] = ae_le(shape, c.child);
        INFO("shape=" << c.shapeName << " Ae parent=" << aeP << " child=" << aeC
                      << " Le parent=" << leP << " child=" << leC);
        REQUIRE(aeC > 0.0);
        REQUIRE(leC > 0.0);
        REQUIRE_THAT(aeC, Catch::Matchers::WithinRel(aeP, 1e-9));
        REQUIRE_THAT(leC, Catch::Matchers::WithinRel(leP, 1e-9));
    }
}
