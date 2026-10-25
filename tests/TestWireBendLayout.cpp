// ABT #1290 (Alf, 2026-09-20): under REAL WINDING, MKF must not draw a conductor corner tighter
// than the wire's own minimum bend radius, and the limit must come from the wire's standard --
// "no magic numbers, all decisions are physics, manufacturing logic and standards".
//
// THESE TESTS ASSERT THE RULE, NOT A NUMBER. Nothing here is pinned to a particular radius: each
// case asks the WIRE itself (WireBend, i.e. IEC 60317-0-1 Table 6 read through IEC 60851-3
// 5.1.1) for its own flexibility radius, measures the standoff the COIL actually placed its
// innermost turn at, and then builds the two formers that put the drawn corner just below and
// just above that limit. A layout change moves both together, so the test keeps asking the same
// question; a test pinned to "0,215 mm throws" would pass for the wrong reason the day the
// layout moves.
//
// The former used is a RECTANGULAR column (an E core): a round or oblong column is all corner,
// its radius is its own, and the rule has nothing to judge there.
#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "constructive_models/Wire.h"
#include "advisers/MagneticFilter.h"
#include "processors/ResultsAnalyzer.h"
#include "processors/MagneticSimulator.h"
#include "physical_models/WireBend.h"
#include "physical_models/MagneticField.h"
#include "physical_models/StrayCapacitance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/CoilMesher.h"
#include "support/Exceptions.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <set>
#include <source_location>
#include <limits>
#include <numbers>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

// A former whose corner radius is exactly what the caller asks for, so the drawn corner can be
// placed either side of the wire's limit without touching anything else in the layout: the
// corner radius feeds the bend only -- the column's faces, the window and hence every turn's
// radial position stay where they were.
void set_former_corner_radius(OpenMagnetics::Coil& coil, double cornerRadius) {
    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto processed = bobbin.get_processed_description().value();
    processed.set_column_corner_radius(cornerRadius);
    bobbin.set_processed_description(processed);
    coil.set_bobbin(bobbin);
}

// The radial gap between the former's face and the centre of the turn laid closest to it: the
// standoff the coil itself chose, read back off the wound geometry rather than assumed.
double innermost_turn_standoff(OpenMagnetics::Coil& coil) {
    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    const auto processed = bobbin.get_processed_description().value();
    const double columnWidth = processed.get_column_width().value();
    double closest = std::numeric_limits<double>::max();
    // get_turns_description() hands the optional back BY VALUE, so the vector has to be a named
    // local or the range-for binds into a temporary that is already gone (-Werror=dangling-reference).
    const auto turns = coil.get_turns_description().value();
    for (const auto& turn : turns) {
        closest = std::min(closest, std::abs(turn.get_coordinates()[0]));
    }
    return closest - columnWidth;
}

OpenMagnetics::Coil build_test_coil() {
    auto wires = std::vector<OpenMagnetics::Wire>({find_wire_by_name("Round 0.4 - Grade 1")});
    return OpenMagneticsTesting::get_quick_coil({6}, {1}, "E 25/13/7", 1,
                                               WindingOrientation::OVERLAPPING,
                                               WindingOrientation::OVERLAPPING,
                                               CoilAlignment::CENTERED,
                                               CoilAlignment::CENTERED, wires);
}

}  // namespace

TEST_CASE("Real winding refuses a corner tighter than the wire's own bend limit, and accepts one "
          "just above it",
          "[constructive-model][coil][real-geometry][wire-bend][abt1290]") {
    auto& bendSettings = Settings::GetInstance();
    bendSettings.reset();
    auto coil = build_test_coil();
    bendSettings.set_coil_use_real_winding_geometry(true);

    // The wire's own limit, from the wire's own standard. Nothing in this test knows the number.
    auto wire = coil.resolve_wire(0);
    const double flexibilityRadius =
        WireBend::get_minimum_bend_radius(wire, BendCriterion::FLEXIBILITY, BendAxis::ROUND);
    REQUIRE(flexibilityRadius > 0);

    // A former comfortably above the limit: this must wind, and it is also how the standoff the
    // coil chooses is measured.
    set_former_corner_radius(coil, flexibilityRadius);
    REQUIRE_NOTHROW(coil.wind());
    REQUIRE(coil.get_turns_description());
    const double standoff = innermost_turn_standoff(coil);
    REQUIRE(standoff > 0);

    // The two formers either side of the limit. The drawn corner is the former's corner plus the
    // turn's standoff, so these put it at 0,99 and 1,01 times the wire's flexibility radius.
    const double justUnder = 0.99 * flexibilityRadius - standoff;
    const double justOver = 1.01 * flexibilityRadius - standoff;
    REQUIRE(justUnder > 0);  // both formers must be real formers, not negative radii
    REQUIRE(justOver > justUnder);

    SECTION("a corner at 0,99 of the wire's flexibility radius is refused") {
        set_former_corner_radius(coil, justUnder);
        REQUIRE_THROWS_AS(coil.wind(), UnwindableBendException);
    }

    SECTION("the refusal names the winding, the drawn radius, the limit and the criterion") {
        set_former_corner_radius(coil, justUnder);
        using Catch::Matchers::ContainsSubstring;
        REQUIRE_THROWS_WITH(coil.wind(),
                            ContainsSubstring(coil.get_functional_description()[0].get_name()) &&
                                ContainsSubstring("bend radius") &&
                                ContainsSubstring("FLEXIBILITY") &&
                                ContainsSubstring("IEC 60317-0-1"));
    }

    SECTION("a corner at 1,01 of the same radius is not refused") {
        set_former_corner_radius(coil, justOver);
        REQUIRE_NOTHROW(coil.wind());
        REQUIRE(coil.get_turns_description());
    }

    bendSettings.reset();
}

TEST_CASE("The ideal layout is not judged against the wire's bend limit",
          "[constructive-model][coil][wire-bend][abt1290]") {
    // Alf's scope, explicit: the rule is enforced ONLY behind the real-winding flag. The same
    // sharp former that is refused above must wind exactly as it always has with the flag off --
    // the classic 2D layout treats the column as a mathematically sharp rectangle by design, and
    // every existing result depends on that.
    auto& bendSettings = Settings::GetInstance();
    bendSettings.reset();
    REQUIRE_FALSE(bendSettings.get_coil_use_real_winding_geometry());
    auto coil = build_test_coil();
    set_former_corner_radius(coil, 0.0);  // as sharp as a former can be asked to be
    REQUIRE_NOTHROW(coil.wind());
    REQUIRE(coil.get_turns_description());
}

// User-facing text states the fact, never a tracker number: exception messages that used to cite
// internal tickets ("Real winding (ABT #1487): ...", "Neutral window field (ABT #1304): ...") are
// read by users who cannot look those numbers up. Each case below triggers one such refusal and
// checks its what() carries no "ABT" while still saying what went wrong.
TEST_CASE("User-facing refusals state the fact and cite no tracker numbers",
          "[constructive-model][coil][wire-bend][no-abt-in-messages]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto messageOf = [](auto&& call) {
        std::string message;
        try {
            call();
        }
        catch (const std::exception& exception) {
            message = exception.what();
        }
        return message;
    };
    auto checkClean = [](const std::string& message, const std::string& fact) {
        INFO(message);
        REQUIRE_FALSE(message.empty());
        CHECK(message.find(fact) != std::string::npos);
        CHECK(message.find("ABT") == std::string::npos);
        CHECK(message.find("abt #") == std::string::npos);
    };

    SECTION("the real-winding bend refusal (Coil)") {
        auto coil = build_test_coil();
        settings.set_coil_use_real_winding_geometry(true);
        auto wire = coil.resolve_wire(0);
        const double flexibilityRadius =
            WireBend::get_minimum_bend_radius(wire, BendCriterion::FLEXIBILITY, BendAxis::ROUND);
        set_former_corner_radius(coil, flexibilityRadius);
        REQUIRE_NOTHROW(coil.wind());
        const double justUnder = 0.99 * flexibilityRadius - innermost_turn_standoff(coil);
        REQUIRE(justUnder > 0);
        set_former_corner_radius(coil, justUnder);
        checkClean(messageOf([&]() { coil.wind(); }), "Real winding only");
    }

    SECTION("the neutral window field on a toroid (CoilMesher)") {
        auto magnetic = OpenMagneticsTesting::get_quick_magnetic("T 36/23/15", json::array(), {12, 12}, 1, "3C90");
        checkClean(messageOf([&]() { CoilMesher::generate_window_wall_sheets(magnetic, std::nullopt); }),
                   "Neutral window field: toroids have no rectangular window boundary");
    }

    SECTION("the parallel current split losses with no per-turn losses (WindingOhmicLosses)") {
        auto coil = build_test_coil();
        REQUIRE(coil.get_turns_description());
        checkClean(messageOf([&]() {
                       WindingOhmicLosses::apply_parallel_current_split(coil, WindingLossesOutput(), ParallelCurrentSplit(), 25);
                   }),
                   "Parallel current split losses: the split needs the losses per turn and per winding");
    }

    SECTION("the gap capacitance of an additive gap with no geometrical description (StrayCapacitance)") {
        auto coil = OpenMagneticsTesting::get_quick_coil({20, 20}, {1, 1}, "ETD 49/25/16", 1,
                                                         WindingOrientation::CONTIGUOUS,
                                                         WindingOrientation::OVERLAPPING);
        auto spacerCore = OpenMagneticsTesting::get_quick_core("ETD 49/25/16", OpenMagneticsTesting::get_spacer_gap(0.001), 1, "N87");
        spacerCore.set_geometrical_description(std::nullopt);
        checkClean(messageOf([&]() { StrayCapacitance::core_gap_topology(spacerCore, coil, "winding 0", "winding 1"); }),
                   "capacitance depends directly on that dielectric");
    }

    SECTION("a catalogue mutation while the databases are frozen (Utils)") {
        set_databases_frozen(true);
        auto message = messageOf([&]() { clear_loaded_cores(); });
        set_databases_frozen(false);
        checkClean(message, "would mutate a shared catalog while databases are frozen");
    }

    settings.reset();
}

// ABT #2004 (Alf, 2026-10-22). IEC 60317-0-2 Table 6 (5 x width above 10 mm) is the mandrel-winding
// TEST of the coating (IEC 60851-3 5.1.1): it qualifies the enamel at that diameter. It is not the
// smallest radius copper strip can be formed to, and an edgewound coil is formed far tighter. So under
// real winding an EDGEWISE bend of rectangular wire below that mandrel is DRAWN, and analyze_results
// reports one warning per winding; the round (above) and FLATWISE (below) refusals of ABT #1290 stay.
namespace {

// Rectangular IEC 60317 grade-1 copper 4.50 x 1 mm from the catalogue, laid as the coil lays it:
// its 4.5 mm width radial, in the bend plane, so the former's corners bend it EDGEWISE.
OpenMagnetics::Wire edgewise_strip() {
    return find_wire_by_name("Rectangular 4.50x1 - Grade 1");
}

// A rectangular IEC 60317 wire stood on its edge: its conducting WIDTH (the dimension the coil lays
// radially, in the bend plane) is the SMALL one, so the corner bends it FLATWISE. The catalogue's
// 2 x 0.80 mm with width and height exchanged.
OpenMagnetics::Wire flatwise_strip_0_8_by_2() {
    auto wire = find_wire_by_name("Rectangular 2x0.80 - Grade 1");
    const auto conductingWidth = wire.get_conducting_width().value();
    const auto conductingHeight = wire.get_conducting_height().value();
    const auto outerWidth = wire.get_outer_width().value();
    const auto outerHeight = wire.get_outer_height().value();
    wire.set_conducting_width(conductingHeight);
    wire.set_conducting_height(conductingWidth);
    wire.set_outer_width(outerHeight);
    wire.set_outer_height(outerWidth);
    wire.set_name("Rectangular 0.80x2 - Grade 1");
    return wire;
}

// A moulded former's 0.5 mm column corner: the strip's edgewise bend around it is far tighter than
// the 4 x width mandrel of IEC 60317-0-2 Table 6.
const double edgewiseFormerCornerRadius = 0.0005;
const std::string edgewiseCoreShape = "E 70/33/32";

OpenMagnetics::Magnetic edgewise_magnetic(const OpenMagnetics::Coil& coil) {
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagneticsTesting::get_quick_core(edgewiseCoreShape, OpenMagneticsTesting::get_ground_gap(0.001), 1, "N87"));
    magnetic.set_coil(coil);
    return magnetic;
}

}  // namespace

TEST_CASE("Real winding draws an edgewise bend of rectangular wire below the IEC 60317-0-2 mandrel, "
          "and analyze_results warns once for its winding",
          "[constructive-model][coil][real-geometry][wire-bend][abt-2004]") {
    auto& bendSettings = Settings::GetInstance();
    bendSettings.reset();
    auto wires = std::vector<OpenMagnetics::Wire>({edgewise_strip()});
    // An edgewound inductor that fits its window and whose terminal leads route: 10 turns of the
    // strip on an E 70/33/32. (The reported 13.6 x 2.5 mm strip on a U 46/40/28 is, once its bends
    // are drawn, refused one step later by its terminal-lead route, which plans the lead at the same
    // Table 6 mandrel; that is a separate rule.)
    auto coil = OpenMagneticsTesting::get_quick_coil({10}, {1}, edgewiseCoreShape, 1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED, wires);
    bendSettings.set_coil_use_real_winding_geometry(true);
    set_former_corner_radius(coil, edgewiseFormerCornerRadius);

    // Wound, fitting and drawn: no refusal (the reported COIL_WINDING_ERROR). The bend rule judges
    // only a layout that fits, so the fit is part of what this proves. Real winding names each
    // station where a turn BEGINS and adds each layer's closing station "..._ending".
    bool fits = false;
    REQUIRE_NOTHROW(fits = coil.wind());
    REQUIRE(fits);
    REQUIRE(coil.get_turns_description());
    {
        const auto stations = coil.get_turns_description().value();
        size_t numberEndings = 0;
        for (const auto& station : stations) {
            if (station.get_name().ends_with("_ending")) {
                ++numberEndings;
            }
        }
        CHECK(stations.size() - numberEndings == 10);
        CHECK(numberEndings == coil.get_layers_description()->size());
    }

    // Not vacuous: the layout really bends this wire EDGEWISE below its own Table 6 limit. Nothing
    // here knows a number; the limit is the wire's own and the radius the coil's own.
    const auto judgement = coil.judge_corners_against_wire_flexibility();
    REQUIRE(judgement.numberCornersJudged > 0);
    REQUIRE(judgement.windingsBelowFlexibility.size() == 1);
    const auto& winding = judgement.windingsBelowFlexibility[0];
    CHECK(winding.axis == BendAxis::EDGEWISE);
    CHECK(winding.bendRadius < winding.asLaid.flexibilityRadius);
    CHECK(winding.asLaid.flexibilityRadius ==
          WireBend::get_minimum_bend_radius(coil.resolve_wire(0), BendCriterion::FLEXIBILITY, BendAxis::EDGEWISE));

    using Catch::Matchers::ContainsSubstring;
    auto magnetic = edgewise_magnetic(coil);

    SECTION("the bend is reported once, for its winding, with the drawn radius and the mandrel") {
        const auto findings = ResultsAnalyzer::analyze_winding_bends(magnetic);
        REQUIRE(findings.size() == 1);
        CHECK(findings[0].code == ResultFindingCode::EDGEWISE_BEND_BELOW_COATING_MANDREL);
        CHECK_FALSE(findings[0].operatingPointIndex);
        CHECK_THAT(findings[0].message,
                   ContainsSubstring(coil.get_functional_description()[0].get_name()) &&
                       ContainsSubstring("edgewise") && ContainsSubstring("centreline radius") &&
                       ContainsSubstring("mandrel") && ContainsSubstring("IEC 60317-0-2 Table 6"));
        json findingsJson;
        to_json(findingsJson, findings);
        CHECK_FALSE(findingsJson[0].contains("operatingPointIndex"));
        CHECK(findingsJson[0]["code"] == "EDGEWISE_BEND_BELOW_COATING_MANDREL");
    }

    SECTION("analyze_results of the simulated design carries that one finding") {
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
            100000, 50e-6, 25, WaveformLabel::TRIANGULAR, 5, 0.5, 0);
        MagneticSimulator magneticSimulator;
        OpenMagnetics::Mas simulated;
        REQUIRE_NOTHROW(simulated = magneticSimulator.simulate(inputs, magnetic));
        const auto findings = ResultsAnalyzer::analyze_results(simulated);
        size_t numberEdgewise = 0;
        for (const auto& finding : findings) {
            if (finding.code == ResultFindingCode::EDGEWISE_BEND_BELOW_COATING_MANDREL) {
                ++numberEdgewise;
                CHECK_FALSE(finding.operatingPointIndex);
            }
        }
        CHECK(numberEdgewise == 1);
    }

    SECTION("the windability filter agrees with real winding: the candidate stays valid") {
        MagneticFilterWindability filter;
        const auto [valid, scoring] = filter.evaluate_magnetic(&magnetic, nullptr);
        CHECK(valid);
        // IEC 60317-0-2 9 states the heat-shock bend flatwise only, so an edgewise bend has no
        // heat-shock radius to fall short of: WireBend gives none and the filter scores nothing.
        CHECK(scoring == 0);
    }

    bendSettings.reset();
}

TEST_CASE("Real winding still refuses a FLATWISE bend of rectangular wire below its IEC 60317-0-2 limit",
          "[constructive-model][coil][real-geometry][wire-bend][abt-2004][abt1290]") {
    auto& bendSettings = Settings::GetInstance();
    bendSettings.reset();
    auto wires = std::vector<OpenMagnetics::Wire>({flatwise_strip_0_8_by_2()});
    auto coil = OpenMagneticsTesting::get_quick_coil({6}, {1}, "E 42/21/15", 1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED, wires);
    bendSettings.set_coil_use_real_winding_geometry(true);

    auto wire = coil.resolve_wire(0);
    REQUIRE(WireBend::axis_from_bend_plane_dimension(wire, true) == BendAxis::FLATWISE);
    const double flexibilityRadius =
        WireBend::get_minimum_bend_radius(wire, BendCriterion::FLEXIBILITY, BendAxis::FLATWISE);
    REQUIRE(flexibilityRadius > 0);

    set_former_corner_radius(coil, flexibilityRadius);
    REQUIRE_NOTHROW(coil.wind());
    REQUIRE(coil.get_turns_description());
    const double standoff = innermost_turn_standoff(coil);
    REQUIRE(standoff > 0);
    const double justUnder = 0.99 * flexibilityRadius - standoff;
    REQUIRE(justUnder > 0);

    set_former_corner_radius(coil, justUnder);
    using Catch::Matchers::ContainsSubstring;
    REQUIRE_THROWS_AS(coil.wind(), UnwindableBendException);
    REQUIRE_THROWS_WITH(coil.wind(), ContainsSubstring(coil.get_functional_description()[0].get_name()) &&
                                         ContainsSubstring("FLEXIBILITY") && ContainsSubstring("IEC 60317-0-2"));

    // The adviser's filter refuses the same former: it judges the turn tangent to the former
    // (corner + half the outer width), which is no further out than the coil's own standoff.
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagneticsTesting::get_quick_core("E 42/21/15", OpenMagneticsTesting::get_ground_gap(0.001), 1, "N87"));
    magnetic.set_coil(coil);
    MagneticFilterWindability filter;
    CHECK_FALSE(filter.evaluate_magnetic(&magnetic, nullptr).first);

    bendSettings.reset();
}

// ABT #2012 (Alf, 2026-10-23, "same as #2004"): a rectangular wire's terminal lead is planned for its
// FLATWISE limit and its coated radius, not for the IEC 60317-0-2 Table 6 EDGEWISE mandrel (a coating
// test, not a forming limit). An edgewise lead bend below that mandrel is drawn, and the winding's
// EDGEWISE_BEND_BELOW_COATING_MANDREL finding names the lead. The design is the reported one: the
// 13.6 x 2.5 mm strip on a U 46/40/28 former with a 3.2 mm corner, with as many turns as fit its window.
namespace {

// Rectangular IEC 60317 grade-1 copper 13.6 x 2.5 mm, the strip of the reported edgewound inductor:
// the catalogue's neighbour 13.20 x 2.50 with its width moved to 13.6 mm; its coating build-up, edge
// radius, standard and grade are kept, and the conducting area follows the new width.
OpenMagnetics::Wire edgewise_strip_13_6_by_2_5() {
    auto wire = find_wire_by_name("Rectangular 13.20x2.50 - Grade 1");
    const double oldWidth = resolve_dimensional_values(wire.get_conducting_width().value());
    const double newWidth = 0.0136;
    const double height = resolve_dimensional_values(wire.get_conducting_height().value());
    const double edgeRadius = resolve_dimensional_values(wire.get_edge_radius().value());
    MAS::DimensionWithTolerance conductingWidth;
    conductingWidth.set_nominal(newWidth);
    wire.set_conducting_width(conductingWidth);
    auto outerWidth = wire.get_outer_width().value();
    outerWidth.set_minimum(outerWidth.get_minimum().value() - oldWidth + newWidth);
    outerWidth.set_maximum(outerWidth.get_maximum().value() - oldWidth + newWidth);
    wire.set_outer_width(outerWidth);
    MAS::DimensionWithTolerance conductingArea;
    conductingArea.set_nominal(newWidth * height - (4 - std::numbers::pi) * edgeRadius * edgeRadius);
    wire.set_conducting_area(conductingArea);
    wire.set_name("Rectangular 13.60x2.50 - Grade 1");
    return wire;
}

const std::string reportedCoreShape = "U 46/40/28";
const double reportedFormerCornerRadius = 0.0032;   // the reported former's ~3.2 mm corner

}  // namespace

TEST_CASE("Real winding draws the edgewise terminal leads of the reported strip below the IEC 60317-0-2 mandrel, "
          "and the winding's edgewise finding names them",
          "[constructive-model][coil][real-geometry][wire-bend][abt-2012]") {
    auto& bendSettings = Settings::GetInstance();
    bendSettings.reset();
    const auto strip = edgewise_strip_13_6_by_2_5();
    const double edgewise = WireBend::get_flexibility_bend_radius_if_standardised(strip, BendAxis::EDGEWISE).value();
    // 16 turns: the most of this strip that fit the U 46/40/28 window under real winding (18, the
    // reported count, overflow it).
    auto coil = OpenMagneticsTesting::get_quick_coil({16}, {1}, reportedCoreShape, 1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED, {strip});
    bendSettings.set_coil_use_real_winding_geometry(true);
    set_former_corner_radius(coil, reportedFormerCornerRadius);

    bool fits = false;
    REQUIRE_NOTHROW(fits = coil.wind());
    REQUIRE(fits);
    REQUIRE(coil.is_real_winding_blocking_applied());

    // Both leads are drawn, planned below the edgewise mandrel: at the strip's flatwise limit or its
    // coated radius, whichever is larger.
    const double flatwise = WireBend::get_flexibility_bend_radius_if_standardised(strip, BendAxis::FLATWISE).value();
    const double coatedRadius = 0.5 * std::max(strip.get_maximum_outer_width(), strip.get_maximum_outer_height());
    size_t terminals = 0;
    for (const auto& route : coil.get_connection_layout().routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) {
            continue;
        }
        ++terminals;
        CHECK(route.plannedBendRadius < edgewise);
        CHECK_THAT(route.plannedBendRadius, Catch::Matchers::WithinRel(std::max(coatedRadius, flatwise), 1e-12));
    }
    REQUIRE(terminals == 2);

    const auto leads = coil.judge_leads_against_wire_flexibility();
    REQUIRE(leads.size() == 1);
    REQUIRE(leads[0].leads.size() == 2);
    for (const auto& lead : leads[0].leads) {
        CHECK(lead.plannedBendRadius < lead.asLaid.flexibilityRadius);
        CHECK(lead.asLaid.flexibilityRadius == edgewise);
    }

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagneticsTesting::get_quick_core(reportedCoreShape, OpenMagneticsTesting::get_ground_gap(0.001), 1, "N87"));
    magnetic.set_coil(coil);
    const auto findings = ResultsAnalyzer::analyze_winding_bends(magnetic);
    REQUIRE(findings.size() == 1);   // one per winding: its turns and its leads together
    CHECK(findings[0].code == ResultFindingCode::EDGEWISE_BEND_BELOW_COATING_MANDREL);
    using Catch::Matchers::ContainsSubstring;
    CHECK_THAT(findings[0].message, ContainsSubstring("entrance lead of '") && ContainsSubstring("exit lead of '") &&
                                        ContainsSubstring("planned for edgewise bends of") &&
                                        ContainsSubstring("IEC 60317-0-2 Table 6"));
    bendSettings.reset();
}
// ABT #2061 (Alf, 2026-10-25: draw and report): a one-turn conductor has no next revolution and no
// pitch, so its tangent leads are drawn level at the turn's height, both shortening that one turn, and
// each is reported (LEAD_ROUTE_NOT_MODELLED). The design is the reported strip with a single turn,
// wound OVERLAPPING, whose 6.84 mm lead bend does not fit on the 6.67 mm connection-face half straight.
TEST_CASE("Real winding draws the tangent leads of a one-turn edgewise strip level and reports them",
          "[constructive-model][coil][real-geometry][wire-bend][tangent-departure][abt-2061]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto coil = OpenMagneticsTesting::get_quick_coil({1}, {1}, reportedCoreShape, 1, WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING, CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED, {edgewise_strip_13_6_by_2_5()});
    settings.set_coil_use_real_winding_geometry(true);
    set_former_corner_radius(coil, reportedFormerCornerRadius);
    bool fits = false;
    REQUIRE_NOTHROW(fits = coil.wind());
    REQUIRE(fits);
    const auto turns = coil.get_turns_description().value();
    REQUIRE(turns.size() == 1);
    size_t level = 0;
    for (const auto& route : coil.get_connection_layout().routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) {
            continue;
        }
        REQUIRE(route.tangentDeparture.has_value());
        const auto& departure = route.tangentDeparture.value();
        REQUIRE(departure.unmodelled.size() == 1);
        CHECK(departure.unmodelled[0].condition == TangentLeadCondition::SINGLE_REVOLUTION);
        // Level: it leaves at the station's own height and runs out without climbing.
        CHECK(departure.point[1] == turns[0].get_coordinates()[1]);
        CHECK(departure.direction[1] == 0.0);
        CHECK(departure.shortenedTurn == turns[0].get_name());
        ++level;
    }
    REQUIRE(level == 2);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagneticsTesting::get_quick_core(reportedCoreShape, OpenMagneticsTesting::get_ground_gap(0.001), 1, "N87"));
    magnetic.set_coil(coil);
    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    size_t reported = 0;
    for (const auto& finding : ResultsAnalyzer::analyze_results(mas)) {
        if (finding.code != ResultFindingCode::LEAD_ROUTE_NOT_MODELLED) {
            continue;
        }
        ++reported;
        using Catch::Matchers::ContainsSubstring;
        CHECK_THAT(finding.message, ContainsSubstring("lead of '" + turns[0].get_name() + "'") &&
                                        ContainsSubstring("its conductor has a single turn and no pitch"));
        CHECK(finding.message.find("ABT") == std::string::npos);
    }
    CHECK(reported == 2);
    settings.reset();
}

// ABT #2071 (Alf, 2026-10-25): a real-winding design whose ideal wind does not fit is still blocked and
// laid out as wound -- turns, connections and leads -- overflowing, and the copper it overlaps is
// reported and painted red, never squished or clamped. The design is the reported one: 18 turns of the
// strip on the U 46/40/28 former, layers CONTIGUOUS (.mkf-p2004 abt2012_magnetic_n18). Its 35 stations
// run radially past the window, and its radial exit lead is drawn across the turns of the layers it
// leaves through.
TEST_CASE("Real winding lays out a design whose ideal wind does not fit, and reports and paints its copper overlaps",
          "[constructive-model][coil][real-geometry][painter][abt-2071]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto coil = OpenMagneticsTesting::get_quick_coil({18}, {1}, reportedCoreShape, 1, WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::CONTIGUOUS, CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED, {edgewise_strip_13_6_by_2_5()});
    settings.set_coil_use_real_winding_geometry(true);
    set_former_corner_radius(coil, reportedFormerCornerRadius);
    bool fits = true;
    REQUIRE_NOTHROW(fits = coil.wind());
    CHECK_FALSE(fits);   // the verdict stays the ideal wind's: it does not fit
    REQUIRE(coil.is_real_winding_blocking_applied());   // ... and it is laid out for real all the same
    const auto turns = coil.get_turns_description().value();
    CHECK(turns.size() == 35);   // 18 turns over 17 layers: one closing station per layer
    CHECK_FALSE(coil.are_turns_inside_winding_window());   // overflowing, not squeezed in
    size_t terminals = 0;
    for (const auto& route : coil.get_connection_layout().routes) {
        terminals += route.kind == ConnectionKind::TERMINAL_ENTRANCE || route.kind == ConnectionKind::TERMINAL_EXIT ? 1 : 0;
    }
    CHECK(terminals == 2);

    const auto overlaps = coil.judge_copper_overlaps();
    REQUIRE_FALSE(overlaps.empty());
    std::set<std::string> crossedTurns;
    std::string crossingLead;
    for (const auto& overlap : overlaps) {
        CHECK(overlap.kind == CopperOverlapKind::LEAD_TURN);
        CHECK(overlap.overlap > 0);
        CHECK(overlap.region.size() >= 3);
        crossingLead = overlap.first;
        crossedTurns.insert(overlap.second);
    }
    CHECK(crossedTurns.size() == overlaps.size());   // one overlap per turn the lead crosses

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagneticsTesting::get_quick_core(reportedCoreShape, OpenMagneticsTesting::get_ground_gap(0.001), 1, "N87"));
    magnetic.set_coil(coil);
    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    size_t reported = 0;
    for (const auto& finding : ResultsAnalyzer::analyze_results(mas)) {
        if (finding.code != ResultFindingCode::LEAD_COLLIDES_WITH_TURNS) {
            continue;
        }
        ++reported;
        using Catch::Matchers::ContainsSubstring;
        CHECK_THAT(finding.message, ContainsSubstring(crossingLead) &&
                                        ContainsSubstring("the worst of " + std::to_string(overlaps.size()) + " turns it crosses"));
        CHECK(finding.message.find("ABT") == std::string::npos);
    }
    CHECK(reported == 1);


    // The stations past the window's outer side: TURNS_OUTSIDE_WINDOW names each with its overrun, and
    // their copper past the window is painted red too.
    const auto outside = coil.judge_turns_outside_window();
    REQUIRE_FALSE(outside.empty());
    auto bobbin = coil.resolve_bobbin();
    const auto window = bobbin.get_processed_description().value().get_winding_windows()[0];
    const double windowOuterSide = window.get_coordinates().value()[0] + window.get_width().value() / 2;
    double reach = 0;
    std::string farthest;
    for (const auto& turn : turns) {
        const double outerEdge = turn.get_coordinates()[0] + turn.get_dimensions().value()[0] / 2;
        if (outerEdge > reach) {
            reach = outerEdge;
            farthest = turn.get_name();
        }
    }
    size_t outsidePieces = 0;
    double worstOverrun = 0;
    const double windowInnerSide = window.get_coordinates().value()[0] - window.get_width().value() / 2;
    for (const auto& turn : outside) {
        // Each station's overrun is its own outline against the window's sides: the stations run out
        // past the outer side, and the first layer's closing station reaches in past the inner one.
        const auto station = std::find_if(turns.begin(), turns.end(), [&](const Turn& candidate) { return candidate.get_name() == turn.turn; });
        REQUIRE(station != turns.end());
        const double halfWidth = station->get_dimensions().value()[0] / 2;
        const double outward = station->get_coordinates()[0] + halfWidth - windowOuterSide;
        const double inward = windowInnerSide - (station->get_coordinates()[0] - halfWidth);
        CHECK(turn.side == (outward >= inward ? TurnOutsideWindow::Side::OUTER : TurnOutsideWindow::Side::INNER));
        CHECK_THAT(turn.overrun, Catch::Matchers::WithinAbs(std::max(outward, inward), 1e-12));
        CHECK(turn.overrun > 0);
        worstOverrun = std::max(worstOverrun, turn.overrun);
        outsidePieces += turn.regions.size();
    }
    CHECK_THAT(worstOverrun, Catch::Matchers::WithinAbs(reach - windowOuterSide, 1e-12));
    size_t outsideReported = 0;
    for (const auto& finding : ResultsAnalyzer::analyze_results(mas)) {
        if (finding.code != ResultFindingCode::TURNS_OUTSIDE_WINDOW) {
            continue;
        }
        ++outsideReported;
        using Catch::Matchers::ContainsSubstring;
        CHECK_THAT(finding.message, ContainsSubstring(std::to_string(outside.size()) + " turn stations outside its winding window") &&
                                        ContainsSubstring("the worst, '" + farthest + "'") &&
                                        ContainsSubstring("radially outward"));
        CHECK(finding.message.find("ABT") == std::string::npos);
    }
    CHECK(outsideReported == 1);
    settings.reset();
}

// ABT #2071 (Alf, 2026-10-25): magnetic_autocomplete winds a stored design with the layers orientation its
// sections state. The reported design (.mkf-p2004 abt2012_magnetic_n18: 18 turns of the strip on the
// U 46/40/28 former, its one section CONTIGUOUS) came back OVERLAPPING, 20 stations. It must be wound
// CONTIGUOUS as given, blocked under real winding, and report the turns it lays outside the window.
TEST_CASE("Autocomplete winds a stored CONTIGUOUS section as given, and reports the turns it lays outside the window",
          "[constructive-model][coil][real-geometry][autocomplete][abt-2071]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("testData")
                           .append("abt2071_edgewise_strip_U46_n18_contiguous.json")
                           .string());
    REQUIRE(file.good());
    const auto stored = json::parse(file);
    REQUIRE(stored.at("coil").at("sectionsDescription").size() == 1);
    REQUIRE(stored.at("coil").at("sectionsDescription")[0].at("layersOrientation") == "contiguous");
    settings.set_coil_use_real_winding_geometry(true);
    auto magnetic = OpenMagnetics::magnetic_autocomplete(OpenMagnetics::Magnetic(stored), json{});
    auto& coil = magnetic.get_mutable_coil();
    const auto sections = coil.get_sections_description().value();
    size_t conduction = 0;
    for (const auto& section : sections) {
        if (section.get_type() == ElectricalType::CONDUCTION) {
            ++conduction;
            CHECK(section.get_layers_orientation() == WindingOrientation::CONTIGUOUS);
        }
    }
    CHECK(conduction == 1);
    const auto layers = coil.get_layers_description().value();
    for (const auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            CHECK(layer.get_orientation() == WindingOrientation::CONTIGUOUS);
        }
    }
    CHECK(coil.get_turns_description().value().size() == 35);   // 18 turns over 17 layers, as stored
    CHECK(coil.is_real_winding_blocking_applied());
    CHECK_FALSE(coil.judge_turns_outside_window().empty());
    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    size_t outsideReported = 0;
    for (const auto& finding : ResultsAnalyzer::analyze_results(mas)) {
        outsideReported += finding.code == ResultFindingCode::TURNS_OUTSIDE_WINDOW ? 1 : 0;
    }
    CHECK(outsideReported == 1);


    // A caller that asks for an orientation itself is honoured over the stored one.
    auto asked = OpenMagnetics::magnetic_autocomplete(OpenMagnetics::Magnetic(stored), json{{"layersOrientation", "overlapping"}});
    const auto askedSections = asked.get_coil().get_sections_description().value();
    for (const auto& section : askedSections) {
        if (section.get_type() == ElectricalType::CONDUCTION) {
            CHECK(section.get_layers_orientation() == WindingOrientation::OVERLAPPING);
        }
    }
    settings.reset();
}

// ABT #2071: two turn stations whose outlines overlap are one TURNS_OVERLAP finding naming them and the
// depth. New design (no wind lays turns over each other: an overflowing one runs out of the window):
// two turns of Round 0.4 on an E 25/13/7, ideal winding, the second station moved down onto the first
// by half its coated diameter -- a stated layout drawn as given.
TEST_CASE("Turn stations laid over each other are reported as overlapping",
          "[constructive-model][coil][abt-2071]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto coil = OpenMagneticsTesting::get_quick_coil({2}, {1}, "E 25/13/7", 1, WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING, CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED, {find_wire_by_name("Round 0.4 - Grade 1")});
    auto turns = coil.get_turns_description().value();
    REQUIRE(turns.size() == 2);
    const double diameter = turns[0].get_dimensions().value()[0];
    turns[1].set_coordinates({turns[0].get_coordinates()[0], turns[0].get_coordinates()[1] + diameter / 2});
    coil.set_turns_description(turns);

    const auto overlaps = coil.judge_copper_overlaps();
    REQUIRE(overlaps.size() == 1);
    CHECK(overlaps[0].kind == CopperOverlapKind::TURN_TURN);
    CHECK(overlaps[0].first == turns[0].get_name());
    CHECK(overlaps[0].second == turns[1].get_name());
    CHECK_THAT(overlaps[0].overlap, Catch::Matchers::WithinAbs(diameter / 2, 1e-12));
    // The lens both circles cover, sampled on their arcs: inside both outlines.
    for (const auto& point : overlaps[0].region) {
        for (const auto& turn : turns) {
            CHECK(std::hypot(point[0] - turn.get_coordinates()[0], point[1] - turn.get_coordinates()[1]) <= diameter / 2 + 1e-12);
        }
    }

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagneticsTesting::get_quick_core("E 25/13/7", OpenMagneticsTesting::get_ground_gap(0.001), 1, "N87"));
    magnetic.set_coil(coil);
    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    size_t reported = 0;
    for (const auto& finding : ResultsAnalyzer::analyze_results(mas)) {
        if (finding.code != ResultFindingCode::TURNS_OVERLAP) {
            continue;
        }
        ++reported;
        using Catch::Matchers::ContainsSubstring;
        CHECK_THAT(finding.message, ContainsSubstring("2 turn stations of the coil overlap one another, in 1 place") &&
                                        ContainsSubstring("'" + turns[0].get_name() + "'") &&
                                        ContainsSubstring("'" + turns[1].get_name() + "'"));
        CHECK(finding.message.find("ABT") == std::string::npos);
    }
    CHECK(reported == 1);
    settings.reset();
}
