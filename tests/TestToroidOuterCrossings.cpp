// ABT #2546: the outer crossings of a toroid (additionalCoordinates) are packed by the turn's real
// cross-section. They used to be packed and checked as circles of the wire's outer height (its
// tangential size), so a flattened litz turn, wider radially than tangentially, stepped later outer
// rings out by the height while the conductor is the width deep: on 57 turns of the design below the
// outer crossings of turns 38 and 52 overlapped by 1.22 mm. Its bore rings also pitched and stacked the
// oval by its size at its centre radius, so its turns overlapped in the bore too (0.38 mm).
//
// Design (catalogue data only, found by the Alvier SMC study design audit, 2026-11-24): T 34/21.8/15,
// High Flux 26, ungapped, bobbin "Basic", one winding of "Litz 15x0.5 - Grade 1 - Single Served", litz
// flattening on.
//
// Design change approved by Alf on 2026-11-24 ("Refuse 57 + wind 55"): this test used to wind 57 turns
// and check they do not overlap. 57 turns only ever fitted this window BECAUSE their turns overlapped:
// once the bore rings are pitched and stacked by the oval's real outline, its rings hold 19/15/11/7/3 =
// 55 turns at the deepest flattening and a sixth ring has no room, either way round. So 57 turns must now
// be refused (CoilException, COIL_WINDING_ERROR), and 55 turns, on the same core and wire, must wind with no
// two turns overlapping.

#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "support/Exceptions.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <optional>
#include <string>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
OpenMagnetics::Magnetic t34_litz15x05(int numberTurns) {
    json input = json::parse(R"({
        "core": {"functionalDescription": {"type": "toroidal", "shape": "T 34/21.8/15", "material": "High Flux 26",
                                           "gapping": [], "numberStacks": 1}},
        "coil": {"bobbin": "Basic", "functionalDescription": [{"name": "Primary", "numberTurns": 0, "numberParallels": 1,
                 "isolationSide": "primary", "wire": "Litz 15x0.5 - Grade 1 - Single Served"}]}
    })");
    input["coil"]["functionalDescription"][0]["numberTurns"] = numberTurns;
    return OpenMagnetics::Magnetic(input);
}
}

TEST_CASE("Test_Toroid_Flattened_Litz_That_Only_Fits_By_Overlapping_Is_Refused", "[constructive-model][coil][toroid][bug]") {
    settings.reset();
    settings.set_coil_litz_deformation(true);
    OpenMagnetics::Magnetic magnetic = t34_litz15x05(57);

    // Judged by the exception type and its ErrorCode, never by the message text (Alf, 2026-11-24:
    // answers between components are enums). MKF has no code of its own for "does not fit its window"
    // yet: autocomplete refuses with the coil-winding code, so that is what is asserted.
    bool refused = false;
    std::optional<ErrorCode> code;
    std::string message;
    try {
        magnetic_autocomplete(magnetic);
    }
    catch (const CoilException& e) {
        refused = true;
        code = e.code();
        message = e.what();
    }
    INFO("message: " << message);
    REQUIRE(refused);
    CHECK(code == ErrorCode::COIL_WINDING_ERROR);
    CHECK(!message.empty());
    settings.reset();
}

TEST_CASE("Test_Toroid_Outer_Crossings_Of_Flattened_Litz_Do_Not_Overlap", "[constructive-model][coil][toroid][bug]") {
    settings.reset();
    settings.set_coil_litz_deformation(true);
    OpenMagnetics::Magnetic magnetic = t34_litz15x05(55);
    OpenMagnetics::Magnetic wound = magnetic_autocomplete(magnetic);
    auto& coil = wound.get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    const auto turns = coil.get_turns_description().value();
    REQUIRE(turns.size() == 55);

    // The premise: the litz is flattened (an oval, either way round: no orientation asserted, per Alf
    // 2026-11-24) and the turns have outer crossings on more than one ring.
    size_t crossings = 0;
    for (const auto& turn : turns) {
        REQUIRE(turn.get_cross_sectional_shape());
        CHECK(turn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::OVAL);
        REQUIRE(turn.get_dimensions());
        const double longSide = std::max(turn.get_dimensions().value()[0], turn.get_dimensions().value()[1]);
        const double shortSide = std::min(turn.get_dimensions().value()[0], turn.get_dimensions().value()[1]);
        CHECK(longSide > shortSide * 1.2);
        if (turn.get_additional_coordinates()) {
            ++crossings;
        }
    }
    CHECK(crossings > 50);
    REQUIRE(coil.get_layers_description());
    size_t conductionLayers = 0;
    const auto layers = coil.get_layers_description().value();   // bound: never range over a temporary
    for (const auto& layer : layers) {
        conductionLayers += layer.get_type() == ElectricalType::CONDUCTION;
    }
    CHECK(conductionLayers >= 2);

    // No two turns' real cross-sections (inner stations and outer crossings) overlap by more than 1 um.
    double deepest = 0;
    std::string deepestPair = "none";
    for (const auto& overlap : coil.judge_copper_overlaps()) {
        if (overlap.kind != CopperOverlapKind::TURN_TURN) {
            continue;
        }
        if (overlap.overlap > deepest) {
            deepest = overlap.overlap;
            deepestPair = overlap.first + " x " + overlap.second;
        }
    }
    INFO("deepest turn-turn overlap " << deepest * 1e3 << " mm: " << deepestPair);
    CHECK(deepest <= 1e-6);
    settings.reset();
}
