// ABT #1705: MKF gives every terminal lead its own lane on the connection face.
//
// Fixture: 14_dab_xfmr_pm8770_n97 (PM 87/70, N97), Primary and Secondary each 12 turns of 4 parallels of
// Round 0.80 Grade 1 (coated OD 0.855 mm), as the C2 corpus processed it. Under real winding each
// winding is laid as an inner layer (11 turns per parallel) and an outer layer (the 12th turn,
// descending to the bottom flange), so a winding's ENTRANCE leads (inner layer) and EXIT leads (outer
// layer) leave by the same flange, on rows 0.855 mm apart: in the window section the two bundles
// merely touch, so the old lane rule let them share the face, and in 3D -- where every lead turns off
// its helix and onto its radial run through bends of its planned radius -- the inner-layer entrance
// leads ran out through the outer layer exactly where the exit leads attach (MVB++: 7 violations, the
// worst 815 um inside the coated envelope).
//
// Every expected value is derived here from the routes' own data (coated diameters from the wires,
// waypoints, plannedBendRadius, rampLength), independently of terminal_exit_slots.
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <source_location>
#include "constructive_models/Wire.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/ResultsAnalyzer.h"
#include "support/Painter.h"
#include "TestingUtils.h"
#include "json.hpp"
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <algorithm>
#include <numbers>
#include <iterator>
#include <sstream>

using json = nlohmann::json;
using OpenMagnetics::Settings;

namespace {

OpenMagnetics::Coil wind_dab_fixture() {
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("testData")
                           .append("abt1705_dab_pm8770_4_parallels.json"));
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    // As MVB++ enriches it for the C2 corpus: real winding, corners planned at 1.05 coated radii
    // (MVB++ kRoundCornerBendFactor), no declared minimum bend, leads ending at the window border.
    auto& settings = Settings::GetInstance();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    auto enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
    auto wound = enriched.get_coil();
    settings.reset();
    REQUIRE(wound.is_real_winding_blocking_applied());
    return wound;
}

}  // namespace

namespace {

using OpenMagnetics::ConnectionKind;
using OpenMagnetics::ConnectionRoute;

bool is_terminal(const ConnectionRoute& route) {
    return route.kind == ConnectionKind::TERMINAL_ENTRANCE || route.kind == ConnectionKind::TERMINAL_EXIT;
}

std::string label(const ConnectionRoute& route) {
    return route.winding + " p" + std::to_string(route.parallel) +
           (route.kind == ConnectionKind::TERMINAL_ENTRANCE ? " entrance" : " exit");
}

// The lane plus, for a ramped lead, the ramp and its radial bend's level run on the side its turn lies
// (ConnectionRoute::rampLength): an entrance's turn leaves towards -x, an exit's arrives from +x.
std::pair<double, double> footprint(const ConnectionRoute& route) {
    const double x = route.exitSlot.value();
    const double span = route.rampLength ? *route.rampLength + route.plannedBendRadius : 0.0;
    return route.kind == ConnectionKind::TERMINAL_ENTRANCE ? std::pair{x - span, x} : std::pair{x, x + span};
}

double gap_between(const std::pair<double, double>& a, const std::pair<double, double>& b) {
    return std::max(a.first - b.second, b.first - a.second);
}

double segment_distance_2d(const std::vector<double>& p0, const std::vector<double>& p1, const std::vector<double>& q0,
                           const std::vector<double>& q1) {
    // Two segments in the plane: the minimum is at an endpoint of one against the other, or 0 when they cross.
    auto point_segment = [](const std::vector<double>& p, const std::vector<double>& a, const std::vector<double>& b) {
        const double dx = b[0] - a[0];
        const double dy = b[1] - a[1];
        const double length2 = dx * dx + dy * dy;
        double t = length2 > 0 ? ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length2 : 0.0;
        t = std::min(1.0, std::max(0.0, t));
        return std::hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy));
    };
    auto cross = [](const std::vector<double>& o, const std::vector<double>& a, const std::vector<double>& b) {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    };
    const double d1 = cross(q0, q1, p0);
    const double d2 = cross(q0, q1, p1);
    const double d3 = cross(p0, p1, q0);
    const double d4 = cross(p0, p1, q1);
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
        return 0.0;
    }
    return std::min({point_segment(p0, q0, q1), point_segment(p1, q0, q1), point_segment(q0, p0, p1), point_segment(q1, p0, p1)});
}

double section_distance(const ConnectionRoute& a, const ConnectionRoute& b) {
    double best = std::numeric_limits<double>::max();
    for (size_t i = 0; i + 1 < a.waypoints.size(); ++i) {
        for (size_t j = 0; j + 1 < b.waypoints.size(); ++j) {
            best = std::min(best, segment_distance_2d(a.waypoints[i], a.waypoints[i + 1], b.waypoints[j], b.waypoints[j + 1]));
        }
    }
    return best;
}

}  // namespace

TEST_CASE("14_dab: every terminal lead has its own lane, clear of every other bundle, with room for its bends (ABT #1705)",
          "[constructive-model][coil][real-winding][abt-1705]") {
    auto coil = wind_dab_fixture();
    const auto layout = coil.get_connection_layout();
    const auto wires = coil.get_wires();
    std::vector<ConnectionRoute> leads;
    for (const auto& route : layout.routes) {
        if (is_terminal(route)) {
            leads.push_back(route);
        }
    }
    // 2 windings x 4 parallels x (entrance, exit), every one with a lane (a round column: no tangential departures).
    REQUIRE(leads.size() == 16);
    auto diameter_of = [&](const ConnectionRoute& route) {
        auto wire = wires[coil.get_winding_index_by_name(route.winding)];
        return std::max({wire.get_maximum_outer_width(), wire.get_maximum_outer_height(), route.sleeveOuterDiameter.value_or(0.0)});
    };
    for (const auto& lead : leads) {
        INFO(label(lead));
        REQUIRE(lead.exitSlot.has_value());
        REQUIRE_FALSE(lead.tangentDeparture.has_value());
        REQUIRE(lead.plannedBendRadius >= diameter_of(lead) / 2);
    }

    size_t crossBundlePairs = 0;
    size_t primaryEntranceVsExit = 0;
    double minimumCrossBundleSurplus = std::numeric_limits<double>::max();
    double minimumCrossBundleSlack = std::numeric_limits<double>::max();   // gap beyond the required one
    for (size_t i = 0; i < leads.size(); ++i) {
        for (size_t j = i + 1; j < leads.size(); ++j) {
            const auto& a = leads[i];
            const auto& b = leads[j];
            if (a.side != b.side) {
                continue;   // another isolation side's face
            }
            INFO(label(a) << " at x " << *a.exitSlot * 1e3 << " mm vs " << label(b) << " at x " << *b.exitSlot * 1e3 << " mm");
            const double radii = (diameter_of(a) + diameter_of(b)) / 2;
            const double gap = gap_between(footprint(a), footprint(b));
            const double section = section_distance(a, b);
            // No two leads within their summed coated radii, siblings included.
            CHECK((gap >= radii - 1e-9 || section >= radii - 1e-9));
            const bool sameBundle = a.winding == b.winding && a.kind == b.kind;
            if (sameBundle) {
                continue;
            }
            // Two bundles whose routes come within reach of each other's corners (the coated radii plus
            // both bends) stand the coated radii plus the larger bend apart along the face: disjoint,
            // with a clearance gap, never an exact touch.
            if (section < radii + a.plannedBendRadius + b.plannedBendRadius) {
                ++crossBundlePairs;
                if (a.winding == "Primary" && b.winding == "Primary") {
                    ++primaryEntranceVsExit;
                }
                const double required = radii + std::max(a.plannedBendRadius, b.plannedBendRadius);
                CHECK(gap >= required - 1e-9);
                CHECK(gap > radii);
                minimumCrossBundleSurplus = std::min(minimumCrossBundleSurplus, gap - radii);
                minimumCrossBundleSlack = std::min(minimumCrossBundleSlack, gap - required);
            }
        }
    }
    // The pairs this test is about exist: the Primary's inner-layer entrances against its outer-layer exits.
    REQUIRE(primaryEntranceVsExit == 16);
    REQUIRE(crossBundlePairs >= primaryEntranceVsExit);
    INFO("smallest cross-bundle clearance beyond the coated radii: " << minimumCrossBundleSurplus * 1e3 << " mm");
    CHECK(minimumCrossBundleSurplus > 0);
    // The tightest bundle boundary stands EXACTLY the required gap apart: the anchor is solved, not
    // stepped along a pitch grid that would waste up to a pitch at every boundary (Alf, ABT #1705).
    INFO("smallest cross-bundle slack beyond the required gap: " << minimumCrossBundleSlack * 1e3 << " mm");
    CHECK(std::abs(minimumCrossBundleSlack) < 1e-9);

    // Each lead leaves its helix with room for its bends: a straight stub from its turn to its row
    // holds two bends of plannedBendRadius (one in the turn surface, one onto the radial run), so it is
    // at least 2R long, or it is published as a ramp on the turn surface (ABT #1336) whose span the
    // lanes above already keep clear.
    size_t stubs = 0;
    for (const auto& lead : leads) {
        // The route runs from the turn (entrance: last waypoint; exit: first) to the border.
        const auto& turnEnd = lead.kind == ConnectionKind::TERMINAL_ENTRANCE ? lead.waypoints.back() : lead.waypoints.front();
        const auto& next = lead.kind == ConnectionKind::TERMINAL_ENTRANCE ? lead.waypoints[lead.waypoints.size() - 2] : lead.waypoints[1];
        const bool axialStub = std::abs(turnEnd[0] - next[0]) < 1e-12 && std::abs(turnEnd[1] - next[1]) > 0;
        if (!axialStub) {
            continue;
        }
        ++stubs;
        const double height = std::abs(turnEnd[1] - next[1]);
        INFO(label(lead) << ": stub " << height * 1e3 << " mm, bend radius " << lead.plannedBendRadius * 1e3 << " mm");
        CHECK((height >= 2 * lead.plannedBendRadius - 1e-12 || lead.rampLength.has_value()));
    }
    REQUIRE(stubs > 0);
}

TEST_CASE("check_terminal_lanes refuses lanes that put two bundles' bent leads together (ABT #1705)",
          "[constructive-model][coil][real-winding][abt-1705]") {
    auto coil = wind_dab_fixture();
    const auto routes = coil.get_connection_layout().routes;
    const auto wires = coil.get_wires();
    std::vector<double> diameters;
    std::vector<std::optional<double>> slots;
    std::optional<size_t> entrance;
    std::optional<size_t> exit;
    for (size_t index = 0; index < routes.size(); ++index) {
        auto wire = wires[coil.get_winding_index_by_name(routes[index].winding)];
        diameters.push_back(std::max({wire.get_maximum_outer_width(), wire.get_maximum_outer_height()}));
        slots.push_back(routes[index].exitSlot);
        if (routes[index].winding == "Primary" && routes[index].parallel == 0 && routes[index].kind == ConnectionKind::TERMINAL_ENTRANCE) {
            entrance = index;
        }
        if (routes[index].winding == "Primary" && routes[index].parallel == 3 && routes[index].kind == ConnectionKind::TERMINAL_EXIT) {
            exit = index;
        }
    }
    REQUIRE(entrance.has_value());
    REQUIRE(exit.has_value());
    // MKF's own lanes pass.
    CHECK_NOTHROW(OpenMagnetics::Coil::check_terminal_lanes(routes, diameters, slots, true));
    // The lanes the old rule handed out put the Primary's parallel 3 exit on its parallel 0 entrance's lane:
    // their rows touch in the window section (0.855 mm), so the old test saw no clash.
    slots[exit.value()] = slots[entrance.value()];
    CHECK_THROWS_WITH(OpenMagnetics::Coil::check_terminal_lanes(routes, diameters, slots, true),
                      Catch::Matchers::ContainsSubstring("Terminal lanes collide") &&
                          Catch::Matchers::ContainsSubstring("Primary' parallel 3 exit") &&
                          Catch::Matchers::ContainsSubstring("different bundles"));
    // Two leads on one point collide even without the bend rule.
    auto coincident = routes;
    coincident[exit.value()].waypoints = coincident[entrance.value()].waypoints;
    CHECK_THROWS_WITH(OpenMagnetics::Coil::check_terminal_lanes(coincident, diameters, slots, false),
                      Catch::Matchers::ContainsSubstring("their coated radii need"));
}

TEST_CASE("terminal_exit_slots: two bundles that touch in the section take disjoint lanes only under real winding (ABT #1705)",
          "[constructive-model][coil][abt-1705]") {
    // 14_dab's Primary in miniature (mm -> m): the inner-layer entrance row and the outer-layer exit row,
    // one coated OD apart in the section; two parallels each, side by side on their rows.
    const double od = 0.855e-3;
    const double bend = 1.05 * od / 2;
    auto route = [&](int64_t parallel, ConnectionKind kind, std::vector<std::vector<double>> waypoints) {
        ConnectionRoute r;
        r.winding = "Primary";
        r.parallel = parallel;
        r.kind = kind;
        r.side = 0;
        r.waypoints = waypoints;
        r.plannedBendRadius = bend;
        (kind == ConnectionKind::TERMINAL_ENTRANCE ? r.toTurn : r.fromTurn) = "t";
        return r;
    };
    std::vector<ConnectionRoute> routes = {
        route(0, ConnectionKind::TERMINAL_ENTRANCE, {{34.5025e-3, -21.4223e-3}, {18.8475e-3, -21.4223e-3}}),
        route(1, ConnectionKind::TERMINAL_ENTRANCE, {{34.5025e-3, -21.4223e-3}, {18.8475e-3, -21.4223e-3}, {18.8475e-3, -19.6957e-3}}),
        route(0, ConnectionKind::TERMINAL_EXIT, {{19.7025e-3, -18.7212e-3}, {19.7025e-3, -20.5673e-3}, {34.5025e-3, -20.5673e-3}}),
        route(1, ConnectionKind::TERMINAL_EXIT, {{19.7025e-3, -17.8086e-3}, {19.7025e-3, -20.5673e-3}, {34.5025e-3, -20.5673e-3}}),
    };
    const std::vector<double> attach = {-21.4223e-3, -19.6957e-3, -18.7212e-3, -17.8086e-3};
    const std::vector<double> diameters(routes.size(), od);

    SECTION("ideal winding: the rows touch, so the old rule lets both bundles start at the plane") {
        auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attach, false);
        CHECK(slots[0].value() == 0.0);
        CHECK(slots[1].value() == od);
        CHECK(slots[3].value() == 0.0);   // exit order by span descending: parallel 1 first
        CHECK(slots[2].value() == od);
    }
    SECTION("real winding: the exit bundle starts exactly d + R past the entrance bundle") {
        auto slots = OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attach, true);
        CHECK(slots[0].value() == 0.0);
        CHECK(slots[1].value() == od);
        // The second entrance stands at od, so the first exit stands at od + (od + bend) = 2.159 mm -- not
        // rounded up to the next pitch (3 od = 2.565 mm).
        const double firstExit = od + od + bend;
        CHECK(std::abs(slots[3].value() - firstExit) < 1e-12);
        CHECK(std::abs(slots[2].value() - (firstExit + od)) < 1e-12);
    }
    SECTION("real winding without a planned bend throws") {
        routes[2].plannedBendRadius = 0;
        CHECK_THROWS_WITH(OpenMagnetics::Coil::terminal_exit_slots(routes, diameters, attach, true),
                          Catch::Matchers::ContainsSubstring("planned bend radius"));
    }
}

// ============================================================================================
// ABT #1383: the terminal-lead exit of CLOSED cores (drum + ring SEATED on a base plate; molded
// body with a lead frame), planned by Coil::plan_closed_core_lead_exits and carried on
// ConnectionLayout (routes' exitRun, leadExitSolids) for MVB++ (#1381) to draw.
//
// Fixtures (generic models; each carries its provenance in magnetic.manufacturerInfo.datasheetInfo,
// saying which values are published and which are inferred or stand-ins):
// - lead_exit_drumring_7332.json: a 7.3 x 7.3 x 3.2 mm drum + ring part, the ring SEATED on a square
//   base plate M = 7.3 mm (the footprint), N = 0 (no chamfer value published), B = the published
//   3.2 mm part height with E = B - D - F and L = B - F, pads 1.5 x 1.5 mm; NiZn ferrite (the
//   published material class; the grade is a stand-in); wire inferred from the published DCR.
// - lead_exit_molded_4030.json: a 4.0 x 4.0 x 3.1 mm molded part, pads 3.25 x 0.82 mm; wire inferred
//   from the published DCR; the powder material is a stand-in.
// - lead_exit_molded_overfull.json: test geometry, not a part: 7 turns of 0.355 mm wire in the same
//   cavity; the copper fits but the lead lane passes the body side, which keeps draw-and-report tested.
// Both bobbins 'Dummy': a drum is its own former and a molded coil has none.
// (Moved here from the lead-exit branch's own test file: terminal-lead tests live together.)
namespace lead_exit_suite {

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::StartsWith;

namespace {

json load_fixture(const std::string& name) {
    std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), name));
    REQUIRE(file.good());
    return json::parse(file);
}

double dimension(const json& fixture, const std::string& letter) {
    return fixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"][letter]["nominal"].get<double>();
}

OpenMagnetics::Magnetic autocomplete(json magneticJson) {
    OpenMagnetics::Magnetic magnetic(magneticJson);
    return magnetic_autocomplete(magnetic);
}

double polyline_length(const std::vector<std::vector<double>>& points, size_t from = 0) {
    double length = 0;
    for (size_t k = from; k + 1 < points.size(); ++k) {
        double sum = 0;
        for (size_t axis = 0; axis < points[k].size(); ++axis) {
            sum += std::pow(points[k + 1][axis] - points[k][axis], 2);
        }
        length += std::sqrt(sum);
    }
    return length;
}

bool is_terminal(const ConnectionRoute& route) {
    return route.kind == ConnectionKind::TERMINAL_ENTRANCE || route.kind == ConnectionKind::TERMINAL_EXIT;
}

size_t count_kind(const ConnectionLayout& layout, LeadExitSolid::Kind kind) {
    size_t count = 0;
    for (const auto& solid : layout.leadExitSolids) {
        count += solid.kind == kind;
    }
    return count;
}

// Every terminal route's routedLength is its in-window run plus its exit run (arc at arc length).
void check_routed_lengths(const ConnectionLayout& layout) {
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        REQUIRE(route.exitRunArc);
        const auto& arc = route.exitRunArc.value();
        const double expected = polyline_length(route.waypoints) + arc.radius * std::abs(arc.toAzimuth - arc.fromAzimuth) +
                                polyline_length(route.exitRun, 1);
        CHECK_THAT(route.routedLength, WithinAbs(expected, 1e-8));
    }
}

}  // namespace

TEST_CASE("A drum with its ring SEATED on a square base plate lays its leads out through ring notches to its electrodes",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::RING_NOTCH);

    const auto layout = coil.get_connection_layout();
    const auto spaces = coil.get_connection_reserved_spaces();
    const double diameter = coil.get_wires()[0].get_maximum_outer_width();
    const double plateTop = -dimension(fixture, "E") / 2;   // the window floor: groove E, centred
    const double ringOuterRadius = dimension(fixture, "J") / 2;

    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        // An edge-routed lead is on the plate-side (bottom) row: its in-window run's outer end is
        // below the window centre. (An outermost-layer end runs radially at its own turn's height.)
        const auto& border = entrance ? route.waypoints.front() : route.waypoints.back();
        for (const auto& space : spaces) {
            if (space.isTerminal && space.kind == route.kind && space.edgeDepth > 0) {
                CHECK(border[1] < 0);
            }
        }
        REQUIRE(route.exitRun.size() >= 3);
        // Every lead passes its notch on the plate-side rows: the exit run's point just outside the
        // ring (the start of its corner run) is below the window centre and on or above the plate.
        std::optional<std::vector<double>> found;
        for (const auto& point : route.exitRun) {
            if (std::abs(std::hypot(point[0], point[2]) - (ringOuterRadius + diameter / 2)) < 1e-9) {
                found = point;
            }
        }
        REQUIRE(found);
        const auto& throughNotch = found.value();
        CHECK(throughNotch[1] < 0);
        CHECK(throughNotch[1] - diameter / 2 >= plateTop - 1e-9);
        REQUIRE(route.exitRunEnd);
        CHECK(route.exitRunEnd.value() == ConnectionRoute::ExitRunEnd::ELECTRODE);
        // The run ends over its own side's terminal: entrance -X, exit +X, both on the -Z half.
        const auto& joint = route.exitRun.back();
        CHECK((entrance ? joint[0] < -ringOuterRadius : joint[0] > ringOuterRadius));
        CHECK(joint[2] < 0);
        // ... lying on the plate top.
        CHECK_THAT(joint[1], WithinAbs(plateTop + diameter / 2, 1e-9));
    }
    CHECK(terminals == 2);
    check_routed_lengths(layout);

    // Two notches, on the diagonals to the two connection-plane corners.
    REQUIRE(count_kind(layout, LeadExitSolid::Kind::CORE_CUT) == 2);
    for (const auto& solid : layout.leadExitSolids) {
        if (solid.kind != LeadExitSolid::Kind::CORE_CUT) {
            continue;
        }
        REQUIRE(solid.body.xz.size() == 4);
        double cx = 0, cz = 0;
        for (const auto& point : solid.body.xz) {
            cx += point[0] / 4;
            cz += point[1] / 4;
        }
        CHECK(cz < 0);
        CHECK(solid.entrance == (cx < 0));
        CHECK_THAT(std::abs(cx), WithinAbs(std::abs(cz), 1e-9));   // on the 45 degree diagonal
        // Width = parallels (1) x the largest coated OD: the inner edge xz[0] -> xz[3].
        const double width = std::hypot(solid.body.xz[3][0] - solid.body.xz[0][0], solid.body.xz[3][1] - solid.body.xz[0][1]);
        CHECK_THAT(width, WithinAbs(1 * diameter, 1e-9));
        // Through the ring wall: the outer edge at the ring's outer radius J/2.
        const double outerCentre = std::hypot((solid.body.xz[1][0] + solid.body.xz[2][0]) / 2, (solid.body.xz[1][1] + solid.body.xz[2][1]) / 2);
        CHECK_THAT(outerCentre, WithinAbs(ringOuterRadius, 1e-9));
        // Height = the lead's height above the plate: its edge row's edgeDepth, or one OD for an
        // outermost-layer end, which descends to lie on the plate.
        CHECK_THAT(solid.body.y0, WithinAbs(plateTop, 1e-12));
        std::optional<double> expectedHeight;
        for (const auto& space : spaces) {
            if (!space.isTerminal || space.kind != (solid.entrance ? ConnectionKind::TERMINAL_ENTRANCE : ConnectionKind::TERMINAL_EXIT)) {
                continue;
            }
            if (space.edgeDepth > 0) {
                expectedHeight = std::max(expectedHeight.value_or(0.0), space.edgeDepth);
            }
        }
        CAPTURE(solid.name, solid.body.y1 - solid.body.y0);
        CHECK_THAT(solid.body.y1 - solid.body.y0, WithinAbs(expectedHeight.value_or(diameter), 1e-9));
        CHECK(solid.material.empty());
    }
    // Each end's terminal: a side leg and a mounting pad (the board face), of the sourced stock.
    CHECK(count_kind(layout, LeadExitSolid::Kind::ELECTRODE) == 4);
    size_t boardFaces = 0;
    for (const auto& solid : layout.leadExitSolids) {
        if (solid.kind == LeadExitSolid::Kind::ELECTRODE) {
            CHECK(solid.material == "phosphor bronze");
            boardFaces += solid.boardFace.has_value();
        }
    }
    CHECK(boardFaces == 2);
    settings.reset();
}

TEST_CASE("A bifilar drum-ring end leaves through one notch, its parallels side by side, every row on the plate",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberTurns"] = 4;
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberParallels"] = 2;
    // Two leads per end side by side on the 1.5 mm terminal. With the fixture's wire both the ideal
    // and the real winding leave room for the two lanes inside the ring bore (the refusal when they
    // do not is covered in "A closed core's lead exit refuses what it cannot build").
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    const auto layout = coil.get_connection_layout();
    const auto spaces = coil.get_connection_reserved_spaces();
    const double diameter = coil.get_wires()[0].get_maximum_outer_width();
    // Every edge-routed terminal row is on the plate-side (bottom) edge.
    size_t edgeRows = 0;
    for (const auto& space : spaces) {
        if (space.isTerminal && space.edgeDepth > 0 &&
            (space.kind == ConnectionKind::TERMINAL_ENTRANCE || space.kind == ConnectionKind::TERMINAL_EXIT)) {
            ++edgeRows;
            CAPTURE(space.coordinates[1]);
            CHECK(space.coordinates[1] < 0);
        }
    }
    CHECK(edgeRows > 0);
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        terminals += is_terminal(route) && !route.exitRun.empty();
    }
    CHECK(terminals == 4);
    // One notch per END, two parallels wide.
    REQUIRE(count_kind(layout, LeadExitSolid::Kind::CORE_CUT) == 2);
    for (const auto& solid : layout.leadExitSolids) {
        if (solid.kind != LeadExitSolid::Kind::CORE_CUT) {
            continue;
        }
        CHECK_FALSE(solid.parallel.has_value());
        const double width = std::hypot(solid.body.xz[3][0] - solid.body.xz[0][0], solid.body.xz[3][1] - solid.body.xz[0][1]);
        CHECK_THAT(width, WithinAbs(2 * diameter, 1e-9));
    }
    check_routed_lengths(layout);
    settings.reset();
}

TEST_CASE("A closed core's terminal leads are laid on the plate-side edge even where the top edge is nearer",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    // Four parallels: the first layer holds the bundle's first turns in slots 1-4 of 6, so the
    // fourth parallel's entrance turn sits ABOVE the window centre and its nearest edge is the top.
    // A closed core's leads leave over its base, so that entrance must still be laid on the bottom.
    // The ring is widened (K 8.0, J = M 9.6 mm) and the pad made 2.0 mm wide so the four lanes and
    // the four leads side by side fit: test geometry, not a catalogue part.
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    auto& winding = fixture["magnetic"]["coil"]["functionalDescription"][0];
    winding["numberTurns"] = 3;
    winding["numberParallels"] = 4;
    for (auto& connection : winding["connections"]) {
        connection["padWidth"] = 0.002;
    }
    auto& dimensions = fixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"];
    dimensions["K"]["nominal"] = 0.008;
    dimensions["J"]["nominal"] = 0.0096;
    dimensions["M"]["nominal"] = 0.0096;
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    // The premise: some entrance turn is above the window centre.
    bool entranceAboveCentre = false;
    const auto layout = coil.get_connection_layout();
    const auto turns = coil.get_turns_description().value();
    for (const auto& route : layout.routes) {
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE) {
            continue;
        }
        for (const auto& turn : turns) {
            if (turn.get_name() == route.toTurn && turn.get_coordinates()[1] > 0) {
                entranceAboveCentre = true;
            }
        }
    }
    REQUIRE(entranceAboveCentre);
    size_t edgeRows = 0;
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        if (space.isTerminal && space.edgeDepth > 0 &&
            (space.kind == ConnectionKind::TERMINAL_ENTRANCE || space.kind == ConnectionKind::TERMINAL_EXIT)) {
            ++edgeRows;
            CAPTURE(space.coordinates[1]);
            CHECK(space.coordinates[1] < 0);
        }
    }
    CHECK(edgeRows >= 4);
    size_t exitRuns = 0;
    for (const auto& route : layout.routes) {
        exitRuns += is_terminal(route) && !route.exitRun.empty();
    }
    CHECK(exitRuns == 8);
    settings.reset();
}

TEST_CASE("A molded body's leads are welded to lead frames under the cavity floor",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_molded_4030.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::LEAD_FRAME);
    // The winding fits its cavity (ideal and real winding alike), so the lead exit is laid out.
    double outerRadius = 0;
    const auto layers = coil.get_layers_description().value();
    for (const auto& layer : layers) {
        outerRadius = std::max(outerRadius, layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2);
    }
    CAPTURE(outerRadius);
    CHECK(outerRadius <= dimension(fixture, "E") / 2);
    const auto layout = coil.get_connection_layout();
    const auto spaces = coil.get_connection_reserved_spaces();
    const double diameter = coil.get_wires()[0].get_maximum_outer_width();
    const double frameThickness = 0.000127;   // src/data/dfm_rules.json leadExitStock
    const double bodyWidth = dimension(fixture, "A"), bodyDepth = dimension(fixture, "C");
    const double bodyHeight = dimension(fixture, "B"), cavityHeight = dimension(fixture, "D");

    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        // An edge-routed lead is on the plate-side (bottom) row (an outermost-layer end runs
        // radially at its own turn's height and descends at its tab).
        const auto& border = entrance ? route.waypoints.front() : route.waypoints.back();
        for (const auto& space : spaces) {
            if (space.isTerminal && space.kind == route.kind && space.edgeDepth > 0) {
                CHECK(border[1] < 0);
            }
        }
        REQUIRE(route.exitRunEnd);
        CHECK(route.exitRunEnd.value() == ConnectionRoute::ExitRunEnd::LEAD_FRAME_WELD);
        // Ends lying on the cavity floor, over its weld.
        CHECK_THAT(route.exitRun.back()[1], WithinAbs(-cavityHeight / 2 + diameter / 2, 1e-9));
    }
    CHECK(terminals == 2);
    check_routed_lengths(layout);

    REQUIRE(count_kind(layout, LeadExitSolid::Kind::LEAD_FRAME) == 2);
    REQUIRE(count_kind(layout, LeadExitSolid::Kind::WELD) == 2);
    double entranceInner = std::numeric_limits<double>::lowest();
    double exitInner = std::numeric_limits<double>::max();
    for (const auto& solid : layout.leadExitSolids) {
        if (solid.kind == LeadExitSolid::Kind::LEAD_FRAME) {
            CHECK(solid.material == "phosphor bronze");
            // In the frame plane under the cavity floor, inside the base pressing and the footprint.
            CHECK_THAT(solid.body.y1, WithinAbs(-cavityHeight / 2, 1e-12));
            CHECK_THAT(solid.body.y1 - solid.body.y0, WithinAbs(frameThickness, 1e-12));
            CHECK(solid.body.y0 >= -bodyHeight / 2 - 1e-12);
            for (const auto& point : solid.body.xz) {
                CHECK(std::abs(point[0]) <= bodyWidth / 2 + 1e-12);
                CHECK(std::abs(point[1]) <= bodyDepth / 2 + 1e-12);
                if (solid.entrance) {
                    entranceInner = std::max(entranceInner, point[0]);
                }
                else {
                    exitInner = std::min(exitInner, point[0]);
                }
            }
        }
        if (solid.kind == LeadExitSolid::Kind::WELD) {
            // One coated OD long and wide, from the frame to the top of the wire lying on it.
            CHECK_THAT(solid.body.y0, WithinAbs(-cavityHeight / 2, 1e-12));
            CHECK_THAT(solid.body.y1 - solid.body.y0, WithinAbs(diameter, 1e-12));
            CHECK_THAT(solid.body.xz[1][0] - solid.body.xz[0][0], WithinAbs(diameter, 1e-12));
            CHECK_THAT(solid.body.xz[2][1] - solid.body.xz[1][1], WithinAbs(diameter, 1e-12));
        }
    }
    // The two frames stand at least one stock thickness apart.
    CAPTURE(entranceInner, exitInner);
    CHECK(exitInner - entranceInner >= frameThickness - 1e-12);
    settings.reset();
}

TEST_CASE("A molded body whose lead lane passes its side draws the lead and reports it",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_molded_overfull.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::LEAD_FRAME);
    // Real winding lays these 7 turns of 0.355 mm in TWO layers whose outer face (1.567 mm) is
    // inside the 1.75 mm cavity: the copper fits the body...
    double outerRadius = 0;
    const auto layers = coil.get_layers_description().value();
    REQUIRE(layers.size() == 2);
    for (const auto& layer : layers) {
        outerRadius = std::max(outerRadius, layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2);
        CHECK(std::abs(layer.get_coordinates()[1]) + layer.get_dimensions()[1] / 2 <= dimension(fixture, "D") / 2 + 1e-12);
    }
    CAPTURE(outerRadius);
    CHECK(outerRadius < dimension(fixture, "E") / 2);
    // ... but its lead lane does not: outer face + 0.384 mm ride-over + half the 0.384 mm pitch puts
    // the lane at 2.142 mm, 2.334 mm with its coated radius, past the body's side at 2.050 mm.
    // Drawn and reported, never refused: each lead is laid where the plan puts it, welded to its
    // lead frame, and the layout records one collision per lead.
    const auto layout = coil.get_connection_layout();
    size_t pastSide = 0;
    for (const auto& clash : layout.leadExitClashes) {
        if (clash.kind == LeadExitClash::Kind::LANE_PAST_BODY_SIDE) {
            ++pastSide;
            CHECK_THAT(clash.description, ContainsSubstring("past the body's side"));
        }
    }
    CHECK(pastSide == 2);
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (is_terminal(route)) {
            ++terminals;
            CHECK(route.exitRunEnd == ConnectionRoute::ExitRunEnd::LEAD_FRAME_WELD);
        }
    }
    CHECK(terminals == 2);
    // ... and analyze_results says so, once per clash.
    const auto findings = ResultsAnalyzer::analyze_lead_exit_clashes(magnetic);
    CHECK(findings.size() == layout.leadExitClashes.size());
    for (const auto& finding : findings) {
        CHECK(finding.code == ResultFindingCode::LEAD_EXIT_CLASHES);
        CHECK_THAT(finding.message, !ContainsSubstring("ABT"));
    }
    settings.reset();
}

// Alf, 2026-10-29: a closed core's lead-exit clashes are reported even when the turns' outlines cannot
// be judged (TURN_OUTLINES_NOT_EVALUABLE): the lead-exit plan records them without any turn outline.
// The molded body above, wound as it is, then its turns' stored shapes and its wire's outer diameter
// taken away, so that no turn has an outline.
TEST_CASE("A molded body's lead-exit clashes are reported when the turns' outlines are unknown",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_molded_overfull.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto coil = magnetic.get_coil();
    auto turns = coil.get_turns_description().value();
    for (auto& turn : turns) {
        turn.set_cross_sectional_shape(std::nullopt);
    }
    coil.set_turns_description(turns);
    auto wire = coil.resolve_wire(size_t(0));
    REQUIRE(wire.get_type() == WireType::ROUND);
    REQUIRE(wire.get_outer_diameter());
    wire.set_outer_diameter(std::nullopt);
    coil.get_mutable_functional_description()[0].set_wire(wire);
    magnetic.set_coil(coil);
    REQUIRE(magnetic.get_mutable_coil().turn_outlines_unknown_reason());
    const auto clashes = magnetic.get_mutable_coil().judge_lead_exit_clashes();
    REQUIRE_FALSE(clashes.empty());

    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    std::vector<ResultFinding> findings;
    REQUIRE_NOTHROW(findings = ResultsAnalyzer::analyze_results(mas));
    size_t notEvaluable = 0;
    size_t reported = 0;
    for (const auto& finding : findings) {
        CAPTURE(finding.message);
        notEvaluable += finding.code == ResultFindingCode::TURN_OUTLINES_NOT_EVALUABLE;
        if (finding.code == ResultFindingCode::LEAD_EXIT_CLASHES) {
            ++reported;
            CHECK_THAT(finding.message, StartsWith("The ") && ContainsSubstring(" lead of turn ") &&
                                        !ContainsSubstring("Warning") && !ContainsSubstring("ABT"));
        }
    }
    CHECK(notEvaluable == 1);
    CHECK(reported == clashes.size());
    settings.reset();
}

TEST_CASE("The fixtures' inferred windings reproduce their published DC resistance",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    // The wire of each fixture is INFERRED from its published DC resistance (see its provenance);
    // this pins that inference: under real winding (turns as laid plus the routed terminal leads)
    // the MKF DC resistance at 25 C lies within the published typical-to-maximum spread of the
    // typical value. The inductance is deliberately NOT checked against the published one.
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    const std::string name = GENERATE(std::string("lead_exit_molded_4030.json"), std::string("lead_exit_drumring_7332.json"));
    CAPTURE(name);
    auto fixture = load_fixture(name);
    const auto& published = fixture["magnetic"]["manufacturerInfo"]["datasheetInfo"]["electrical"][0]["dcResistance"];
    const double typical = published["nominal"].get<double>();
    const double maximum = published["maximum"].get<double>();
    auto magnetic = autocomplete(fixture["magnetic"]);
    const double resistance = WindingOhmicLosses::calculate_dc_resistance_per_winding(magnetic.get_mutable_coil(), 25)[0];
    CAPTURE(resistance, typical, maximum);
    CHECK_THAT(resistance, WithinAbs(typical, maximum - typical));
    settings.reset();
}

TEST_CASE("A catalogue drumRing (no base plate M) has no lead exit, and says so",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["core"]["functionalDescription"]["shape"] = "DR 5.3x3.6 + RIS 7.3x3.5x5.7";
    // Autocomplete and the magnetics stay as they were for a catalogue part...
    auto magnetic = autocomplete(fixture["magnetic"]);
    CHECK_FALSE(magnetic.get_mutable_coil().get_core_lead_exit());
    // ... and laying its leads out is refused, naming why.
    REQUIRE_THROWS_WITH(magnetic.get_mutable_coil().get_connection_layout(), ContainsSubstring("has no lead exit"));
    REQUIRE_THROWS_WITH(magnetic.get_mutable_core().get_lead_exit_geometry(), ContainsSubstring("ring-around-both-flanges"));
}

// Alf, 2026-10-26: laying out a catalogue drumRing's leads stays refused (get_connection_layout throws),
// but analyze_results checks the core for its base plate M BEFORE any lead judgement lays the leads out,
// and reports the part's leads as not modelled instead of throwing.
TEST_CASE("A catalogue drumRing under real winding: analyze_results reports its leads as not laid out, and does not throw",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["core"]["functionalDescription"]["shape"] = "DR 5.3x3.6 + RIS 7.3x3.5x5.7";
    auto magnetic = autocomplete(fixture["magnetic"]);
    REQUIRE(magnetic.get_mutable_coil().is_real_winding_blocking_applied());
    REQUIRE(magnetic.get_mutable_coil().get_turns_description());
    // Laying the leads out is still refused...
    REQUIRE_THROWS_WITH(magnetic.get_mutable_coil().get_connection_layout(), ContainsSubstring("has no lead exit"));
    // ... and analyze_results says so in one finding instead of throwing it.
    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    std::vector<ResultFinding> findings;
    REQUIRE_NOTHROW(findings = ResultsAnalyzer::analyze_results(mas));
    size_t notLaidOut = 0;
    for (const auto& finding : findings) {
        CAPTURE(finding.message);
        CHECK(finding.code != ResultFindingCode::LEAD_EXIT_CLASHES);
        if (finding.code == ResultFindingCode::LEAD_ROUTE_NOT_MODELLED) {
            ++notLaidOut;
            CHECK_THAT(finding.message, ContainsSubstring("states no base plate (letter M)"));
            CHECK_THAT(finding.message, !ContainsSubstring("ABT"));
        }
    }
    CHECK(notLaidOut == 1);
    settings.reset();
}

// Alf, 2026-10-26: painting a catalogue drumRing (no base plate M) under real winding paints its core,
// bobbin and turns, and no terminal leads: the part states no lead exit, so its leads have no route.
TEST_CASE("A catalogue drumRing under real winding paints its core, bobbin and turns without leads, and does not throw",
          "[lead-exit][constructive-model][coil][painter][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["core"]["functionalDescription"]["shape"] = "DR 5.3x3.6 + RIS 7.3x3.5x5.7";
    auto magnetic = autocomplete(fixture["magnetic"]);
    REQUIRE(magnetic.get_mutable_coil().is_real_winding_blocking_applied());
    REQUIRE(magnetic.get_mutable_coil().get_turns_description());
    const size_t turns = magnetic.get_mutable_coil().get_turns_description()->size();
    auto outFile = std::filesystem::path{std::source_location::current().file_name()}.parent_path().append("..").append("output");
    std::filesystem::create_directories(outFile);
    outFile.append("Test_Lead_Exit_Catalogue_DrumRing_Painted.svg");
    std::filesystem::remove(outFile);
    Painter painter(outFile);
    REQUIRE_NOTHROW(painter.paint_magnetic(magnetic));
    painter.export_svg();
    REQUIRE(std::filesystem::exists(outFile));
    std::ifstream svgFile(outFile);
    const std::string svg((std::istreambuf_iterator<char>(svgFile)), std::istreambuf_iterator<char>());
    CHECK_THAT(svg, ContainsSubstring("ferrite"));
    CHECK_THAT(svg, ContainsSubstring("bobbin"));
    CHECK_THAT(svg, ContainsSubstring("copper"));
    for (size_t turn = 0; turn < turns; ++turn) {
        CAPTURE(turn);
        CHECK_THAT(svg, ContainsSubstring(magnetic.get_mutable_coil().get_turns_description().value()[turn].get_name()));
    }
    CHECK_THAT(svg, !ContainsSubstring("class=\"connection_terminal\""));
    settings.reset();
}

TEST_CASE("A closed core's lead exit refuses what it cannot build",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    SECTION("no pad data on the winding's connections") {
        auto fixture = load_fixture("lead_exit_drumring_7332.json");
        fixture["magnetic"]["coil"]["functionalDescription"][0].erase("connections");
        REQUIRE_THROWS_WITH(autocomplete(fixture["magnetic"]).get_mutable_coil().get_connection_layout(),
                            ContainsSubstring("padWidth and padDepth"));
    }
    SECTION("the lead frame is thicker than the molded body's base") {
        auto fixture = load_fixture("lead_exit_molded_4030.json");
        // Base (B - D)/2 = 0.1 mm, under the 0.127 mm sourced frame stock.
        fixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"]["B"]["nominal"] = 0.0023;
        REQUIRE_THROWS_WITH(autocomplete(fixture["magnetic"]).get_mutable_coil().get_connection_layout(),
                            ContainsSubstring("thicker than the molded body's base"));
    }
    settings.reset();
}

namespace {

// The leads of a drum-ring part are drawn to their electrodes; returns how many terminal routes there are.
size_t check_leads_drawn_to_electrodes(const ConnectionLayout& layout) {
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        CHECK(route.exitRunEnd == ConnectionRoute::ExitRunEnd::ELECTRODE);
        CHECK(route.exitRun.size() >= 2);
        CHECK(route.exitRunArc);
    }
    return terminals;
}

size_t count_ring_bore_clashes(const ConnectionLayout& layout) {
    size_t pastBore = 0;
    for (const auto& clash : layout.leadExitClashes) {
        if (clash.kind == LeadExitClash::Kind::LANE_PAST_RING_BORE) {
            ++pastBore;
            CHECK_THAT(clash.description, ContainsSubstring("past the ring bore"));
        }
    }
    return pastBore;
}

}  // namespace

TEST_CASE("A drum-ring lead lane that passes the ring bore is drawn and reported",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    SECTION("the band between the winding and the ring bore is too narrow for the lane") {
        settings.set_coil_use_real_winding_geometry(false);
        auto fixture = load_fixture("lead_exit_drumring_7332.json");
        // 10 turns of 0.355 mm wire (0.384 mm pitch), ideal winding: two layers whose outer face is at
        // 2.017 mm. The lane sits at 2.592 mm (outer face + 0.384 mm ride-over + half its pitch) and
        // with its 0.192 mm bend reaches 2.784 mm, past a ring bore brought in to K/2 = 2.675 mm.
        auto& winding = fixture["magnetic"]["coil"]["functionalDescription"][0];
        winding["wire"] = "Round 0.355 - Grade 1";
        winding["numberTurns"] = 10;
        winding["numberParallels"] = 1;
        fixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"]["K"]["nominal"] = 0.00535;
        auto magnetic = autocomplete(fixture["magnetic"]);
        ConnectionLayout layout;
        REQUIRE_NOTHROW(layout = magnetic.get_mutable_coil().get_connection_layout());
        CHECK(check_leads_drawn_to_electrodes(layout) == 2);
        CHECK(count_ring_bore_clashes(layout) == 2);
    }
    SECTION("a bifilar end's two lanes and bend do not fit inside the ring bore under real winding") {
        settings.set_coil_use_real_winding_geometry(true);
        auto fixture = load_fixture("lead_exit_drumring_7332.json");
        // 4 bifilar turns of 0.315 mm wire (0.341 mm pitch): real winding lays the wires in three
        // layers with an outer face at 2.274 mm. The outer lane of an end sits at 2.787 mm and with
        // its 0.328 mm bend reaches 3.115 mm, past the 2.85 mm ring bore: one clash per lead.
        auto& winding = fixture["magnetic"]["coil"]["functionalDescription"][0];
        winding["wire"] = "Round 0.315 - Grade 1";
        winding["numberTurns"] = 4;
        winding["numberParallels"] = 2;
        auto magnetic = autocomplete(fixture["magnetic"]);
        REQUIRE(magnetic.get_mutable_coil().is_real_winding_blocking_applied());
        ConnectionLayout layout;
        REQUIRE_NOTHROW(layout = magnetic.get_mutable_coil().get_connection_layout());
        CHECK(check_leads_drawn_to_electrodes(layout) == 4);
        CHECK(count_ring_bore_clashes(layout) == 4);
    }
    settings.reset();
}

TEST_CASE("A seated drumRing closes its bottom through the seat joint, not a radial clearance",
          "[lead-exit][constructive-model][core][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    OpenMagnetics::Core seated(fixture["magnetic"]["core"]);
    auto gapping = seated.get_functional_description().get_gapping();
    REQUIRE(gapping.size() == 2);
    // Top: the radial clearance (K - A)/2 across the top flange; bottom: the residual seat joint.
    const auto& top = gapping[0].get_coordinates().value()[1] > 0 ? gapping[0] : gapping[1];
    const auto& seat = gapping[0].get_coordinates().value()[1] > 0 ? gapping[1] : gapping[0];
    CHECK_THAT(top.get_length(), WithinAbs((dimension(fixture, "K") - dimension(fixture, "A")) / 2, 1e-12));
    CHECK_THAT(seat.get_length(), WithinAbs(Constants().residualGap, 1e-15));
    CHECK_THAT(seat.get_area().value(), WithinAbs(std::numbers::pi / 4 * (std::pow(dimension(fixture, "J"), 2) - std::pow(dimension(fixture, "K"), 2)), 1e-9));
    // The envelope is the square plate.
    CHECK_THAT(seat.get_coordinates().value()[1], WithinAbs(-dimension(fixture, "B") / 2 + dimension(fixture, "F"), 1e-12));
    auto processed = seated.get_processed_description().value();
    CHECK_THAT(processed.get_width(), WithinAbs(dimension(fixture, "M"), 1e-12));
    CHECK_THAT(processed.get_depth(), WithinAbs(dimension(fixture, "M"), 1e-12));
    // Magnetically: the seat joint (Constants::residualGap across the ring annulus) replaces the
    // bottom flange's radial clearance, so the seated part's inductance is the higher one. Like for
    // like: the same drum and ring letters in the catalogue construction (no plate; the ring spans
    // both flanges, L = B).
    auto seatedMagnetic = autocomplete(fixture["magnetic"]);
    auto catalogueFixture = fixture;
    auto& catalogueDimensions = catalogueFixture["magnetic"]["core"]["functionalDescription"]["shape"]["dimensions"];
    catalogueDimensions.erase("M");
    catalogueDimensions.erase("N");
    catalogueDimensions["L"]["nominal"] = dimension(fixture, "B");
    auto catalogueMagnetic = autocomplete(catalogueFixture["magnetic"]);
    const double seatedInductance = MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(seatedMagnetic)
                                        .get_magnetizing_inductance().get_nominal().value();
    const double catalogueInductance = MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(catalogueMagnetic)
                                           .get_magnetizing_inductance().get_nominal().value();
    CHECK(std::isfinite(seatedInductance));
    CHECK(seatedInductance > catalogueInductance);
    // A plate that does not carry its ring is refused.
    auto small = fixture["magnetic"]["core"];
    small["functionalDescription"]["shape"]["dimensions"]["M"]["nominal"] = 0.0070;
    REQUIRE_THROWS_WITH(OpenMagnetics::Core(small), ContainsSubstring("must carry the ring"));
    // Chamfers that leave no corner outside the ring are refused.
    auto chamfered = fixture["magnetic"]["core"];
    chamfered["functionalDescription"]["shape"]["dimensions"]["N"]["nominal"] = 0.0022;
    REQUIRE_THROWS_WITH(OpenMagnetics::Core(chamfered), ContainsSubstring("no corner area"));
}

// The lane round the winding clears MKF's OWN ride-over. On the seated drumRing 7332 fixture a return
// is laid under the outer layer (a Z dragback: 7 + 2 turns), so the outer face as wound stands one
// ride (one OD) further out than the layers' outer face; a consumer lifting the attached turn's end
// by that ride must still find the lead one coated diameter outside it. Before the lane took the
// ride, exitRun[0] sat ~3 um outside the lifted attach.
TEST_CASE("A closed core's lead lane clears the ride-over MKF's own returns impose on the winding's outer face",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::RING_NOTCH);
    const auto layout = coil.get_connection_layout();
    const double diameter = coil.get_wires()[0].get_maximum_outer_width();
    double windingOuterRadius = std::numeric_limits<double>::lowest();
    const auto layers = coil.get_layers_description().value();
    for (const auto& layer : layers) {
        windingOuterRadius = std::max(windingOuterRadius, layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2);
    }
    const auto turns = coil.get_turns_description().value();
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        CAPTURE(entrance, route.side);
        const double ride = layout.ride_at(windingOuterRadius, route.side);
        CAPTURE(ride, windingOuterRadius, diameter);
        REQUIRE(ride > 0);   // anti-vacuous: this fixture's lanes DO have a ride to clear
        REQUIRE(route.exitRunRide);
        CHECK_THAT(route.exitRunRide.value(), WithinAbs(ride, 1e-12));
        REQUIRE(route.exitRunArc);
        CHECK(route.exitRunArc->radius >= windingOuterRadius + ride + diameter / 2 - 1e-9);
        const std::string turnName = entrance ? route.toTurn : route.fromTurn;
        std::optional<double> attachedRadius;
        for (const auto& turn : turns) {
            if (turn.get_name() == turnName) {
                attachedRadius = turn.get_coordinates()[0];
            }
        }
        REQUIRE(attachedRadius);
        REQUIRE(route.exitRun.size() >= 1);
        const double startRadius = std::hypot(route.exitRun[0][0], route.exitRun[0][2]);
        CAPTURE(attachedRadius.value(), startRadius);
        CHECK(startRadius - (attachedRadius.value() + ride) >= diameter - 1e-9);
    }
    CHECK(terminals == 2);
    settings.reset();
}

TEST_CASE("A molded body whose turns share one layer has no ride for its lead lane to clear",
          "[lead-exit][constructive-model][coil][abt1383][abt-1383]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    CAPTURE(realWinding);
    auto fixture = load_fixture("lead_exit_molded_4030.json");
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_core_lead_exit());
    REQUIRE(coil.get_core_lead_exit()->kind == CoreLeadExitGeometry::Kind::LEAD_FRAME);
    const auto layout = coil.get_connection_layout();
    // No return lays a lane (no Z dragback, no edge continuation), so nothing rides over anything.
    CHECK(layout.rideLevels.empty());
    size_t terminals = 0;
    for (const auto& route : layout.routes) {
        CHECK(route.kind != ConnectionKind::Z_DRAGBACK);
        CHECK(route.kind != ConnectionKind::EDGE_CONTINUATION);
        if (!is_terminal(route)) {
            continue;
        }
        ++terminals;
        REQUIRE(route.exitRunRide);
        CHECK(route.exitRunRide.value() == 0);
    }
    CHECK(terminals == 2);
    settings.reset();
}

// ABT #1424 / #1383: a bifilar Z winding's layer changes are DRAGBACKS -- each parallel finishes
// layer 0 at the top of the window and restarts layer 1 at the bottom -- and they are judged as
// what they are. The #1424 plane-link rule used to classify them by the slope between the two
// stations (dy <= dx + 2 OD K), which admitted this 1.3 mm return over a 0.27 mm radial step as a
// "link" and refused the two parallels 0.0577 mm apart in the connection plane, where nothing is
// drawn: a dragback steps out radially at the height its layer finished and runs back along the
// lane at the next layer's radius, and sibling dragbacks lie side by side along that lane.
//
// Everything checked here is re-derived from the stations, not read from the rule it guards:
//   - the return is a dragback because the conductor climbed through layer 0 and goes BACK the
//     other way to start layer 1 -- the definition, from the turn order alone;
//   - the straight station-to-station diagonals of the two parallels pass closer than one OD (the
//     number the slope-classified rule refused), while the drawn chains share the lane radius and
//     overlap along it, so the two need lanes of their own around the column;
//   - two lanes around a round column need 2 asin(OD / r) of azimuth at the innermost radius r
//     the chains occupy, and a revolution has 2 pi;
//   - MKF's ride model agrees they lie side by side: ONE ride level, one OD tall, at the lane.
TEST_CASE("A bifilar drum-ring's Z dragbacks are laid side by side along their lane, not judged as plane links",
          "[lead-exit][constructive-model][coil][abt1424-dragback][abt-1424]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberTurns"] = 4;
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberParallels"] = 2;
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.resolve_bobbin().get_processed_description().value().get_column_shape() == ColumnShape::ROUND);
    const auto turns = coil.get_turns_description().value();
    const double od = coil.get_wires()[0].get_maximum_outer_height();
    const double width = coil.get_wires()[0].get_maximum_outer_width();

    struct Return {
        std::string from, to;
        double xa, ya, xb, yb;
    };
    std::vector<Return> returns;
    for (int64_t parallel = 0; parallel < 2; ++parallel) {
        std::vector<const Turn*> sequence;
        for (const auto& turn : turns) {
            if (turn.get_parallel() == parallel) {
                sequence.push_back(&turn);
            }
        }
        for (size_t k = 1; k + 1 < sequence.size(); ++k) {
            const Turn& previous = *sequence[k - 1];
            const Turn& a = *sequence[k];
            const Turn& b = *sequence[k + 1];
            if (previous.get_layer() != a.get_layer() || a.get_layer() == b.get_layer()) {
                continue;
            }
            const double climb = a.get_coordinates()[1] - previous.get_coordinates()[1];
            const double back = b.get_coordinates()[1] - a.get_coordinates()[1];
            if (climb * back < 0) {
                returns.push_back({a.get_name(), b.get_name(), a.get_coordinates()[0], a.get_coordinates()[1],
                                   b.get_coordinates()[0], b.get_coordinates()[1]});
            }
        }
    }
    REQUIRE(returns.size() == 2);

    const auto layout = coil.get_connection_layout();
    for (const auto& r : returns) {
        size_t found = 0;
        for (const auto& route : layout.routes) {
            if (route.fromTurn == r.from && route.toTurn == r.to) {
                ++found;
                CHECK(route.kind == ConnectionKind::Z_DRAGBACK);
            }
        }
        CAPTURE(r.from, r.to);
        CHECK(found == 1);
    }

    const auto& p = returns[0];
    const auto& q = returns[1];
    const auto pointToSegment = [](double px, double py, double ax, double ay, double bx, double by) {
        const double vx = bx - ax, vy = by - ay;
        const double u = std::clamp(((px - ax) * vx + (py - ay) * vy) / (vx * vx + vy * vy), 0.0, 1.0);
        return std::hypot(px - (ax + u * vx), py - (ay + u * vy));
    };
    const double diagonals = std::min({pointToSegment(p.xa, p.ya, q.xa, q.ya, q.xb, q.yb),
                                       pointToSegment(p.xb, p.yb, q.xa, q.ya, q.xb, q.yb),
                                       pointToSegment(q.xa, q.ya, p.xa, p.ya, p.xb, p.yb),
                                       pointToSegment(q.xb, q.yb, p.xa, p.ya, p.xb, p.yb)});
    CAPTURE(diagonals, od);
    CHECK(diagonals < od);
    // The descents share the lane radius and overlap along it: one plane cannot hold both.
    CHECK_THAT(p.xb, WithinAbs(q.xb, 1e-9));
    CHECK(std::max(std::min(p.ya, p.yb), std::min(q.ya, q.yb)) < std::min(std::max(p.ya, p.yb), std::max(q.ya, q.yb)));
    const double innermost = std::min({p.xa, p.xb, q.xa, q.xb});
    CAPTURE(innermost);
    REQUIRE(od <= innermost);
    CHECK(2 * std::asin(od / innermost) <= 2 * std::numbers::pi);

    size_t levelsAtLane = 0;
    for (const auto& level : layout.rideLevels) {
        if (std::abs(level.radius - p.xb) <= 0.5 * width) {
            ++levelsAtLane;
            CHECK_THAT(level.height, WithinAbs(width, 1e-12));
        }
    }
    CHECK(levelsAtLane == 1);
    settings.reset();
}

// The dragback rule refuses a lane that cannot hold its siblings. The same wound bifilar layout,
// moved radially inward (every station by the same amount, so the layout's kinds -- decided by
// the winding order and the layers -- are unchanged) until the source layer sits at 0.99 OD from
// the axis: two radial run-outs there cannot be separated by one diameter at any azimuth. (Not
// further in: the terminal lead leaving that layer must still ramp off its turn, which needs the
// turn's radius above the lead's planned bend radius, ABT #1336 -- and that refusal would come
// first.)
TEST_CASE("Sibling Z dragbacks whose lane cannot hold them side by side are refused",
          "[constructive-model][coil][abt1424-dragback][abt-1424]") {
    auto& settings = Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto fixture = load_fixture("lead_exit_drumring_7332.json");
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberTurns"] = 4;
    fixture["magnetic"]["coil"]["functionalDescription"][0]["numberParallels"] = 2;
    auto magnetic = autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    const double od = coil.get_wires()[0].get_maximum_outer_height();
    // As wound: accepted.
    REQUIRE_NOTHROW(coil.refuse_sibling_layer_links_closer_than_the_wire());
    auto turns = coil.get_turns_description().value();
    double innermost = std::numeric_limits<double>::max();
    for (const auto& turn : turns) {
        innermost = std::min(innermost, turn.get_coordinates()[0]);
    }
    const double shift = innermost - 0.99 * od;
    REQUIRE(shift > 0);
    for (auto& turn : turns) {
        auto coordinates = turn.get_coordinates();
        coordinates[0] -= shift;
        turn.set_coordinates(coordinates);
    }
    coil.set_turns_description(turns);
    REQUIRE_THROWS_WITH(coil.refuse_sibling_layer_links_closer_than_the_wire(),
                        ContainsSubstring("Z dragbacks") && ContainsSubstring("side by side along the lane") &&
                            ContainsSubstring("no azimuth separates two radial run-outs by one diameter"));
    settings.reset();
}

}  // namespace lead_exit_suite

// ABT #1397: the seated drum-ring's coil ends are joined to their terminal legs by a BOND solid (the
// wire's own conductor, thermocompression-bonded): from the leg's top face (a shared face) up to the
// wire's centreline, so the bare conductor -- one coating thickness above the leg -- reaches the metal.
TEST_CASE("A seated drum-ring's coil ends are bonded to their terminal legs through the wire's coating gap",
          "[lead-exit][constructive-model][coil][abt-1397]") {
    using namespace lead_exit_suite;
    auto& settings = Settings::GetInstance();
    settings.reset();
    const bool realWinding = GENERATE(false, true);
    settings.set_coil_use_real_winding_geometry(realWinding);
    auto fixture = lead_exit_suite::load_fixture("lead_exit_drumring_7332.json");
    auto magnetic = lead_exit_suite::autocomplete(fixture["magnetic"]);
    auto& coil = magnetic.get_mutable_coil();
    const auto layout = coil.get_connection_layout();
    const double coated = coil.get_wires()[0].get_maximum_outer_width();
    const double bare = coil.get_wires()[0].get_maximum_conducting_width();
    REQUIRE(bare < coated);
    size_t bonds = 0;
    for (const auto& route : layout.routes) {
        if (!lead_exit_suite::is_terminal(route)) {
            continue;
        }
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        const LeadExitSolid* bond = nullptr;
        const LeadExitSolid* leg = nullptr;
        for (const auto& solid : layout.leadExitSolids) {
            if (solid.entrance != entrance) continue;
            if (solid.kind == LeadExitSolid::Kind::BOND && solid.parallel == route.parallel) bond = &solid;
            if (solid.kind == LeadExitSolid::Kind::ELECTRODE && solid.name.find("terminal leg") != std::string::npos) leg = &solid;
        }
        REQUIRE(bond != nullptr);
        REQUIRE(leg != nullptr);
        ++bonds;
        // Shared face: the bond stands on the leg's top face.
        CHECK(std::abs(bond->body.y0 - leg->body.y1) < 1e-12);
        // It reaches the wire's centreline (the run's last point) and so overlaps the bare conductor,
        // whose bottom is one coating thickness above the leg.
        const auto& end = route.exitRun.back();
        CHECK(std::abs(bond->body.y1 - end[1]) < 1e-12);
        CHECK(bond->body.y1 > leg->body.y1 + (coated - bare) / 2);
        // Volume > 0 (shoelace area x height), and the wire's own conductor.
        double area = 0;
        const auto& xz = bond->body.xz;
        for (size_t k = 0; k < xz.size(); ++k) {
            const auto& p = xz[k];
            const auto& q = xz[(k + 1) % xz.size()];
            area += p[0] * q[1] - q[0] * p[1];
        }
        CHECK(std::abs(area) / 2 * (bond->body.y1 - bond->body.y0) > 0);
        CHECK(bond->material == coil.resolve_wire(0).resolve_material().get_name());
    }
    CHECK(bonds == 2);
    settings.reset();
}
