// ABT #1172 (WP3 of docs/2026-09-12_manufacturing_features_proposal.md): lead-to-pin assignment,
// MAS RFC 0013 section 3. Kept out of Coil.cpp so the winding code and the pin plan can be read
// (and reverted) separately.
//
// What this file decides, and where each decision comes from:
//   * Rows are isolation groups (rule 1). One isolation side per row while there are enough
//     rows; more sides than rows share rows as contiguous blocks with the spare pins between
//     them. Every pin pair of different sides is then checked against the required creepage
//     over two paths: the straight distance between pin bases (a lower bound on any surface
//     path along the flange) and pin base -> core bounding box -> pin base.
//   * Inside a row or block, the first-wound winding takes the INNER pins and later windings
//     nest outward around it (Wuerth DFM 2024-04-29 08:42: "wind the first windings using the
//     internal bobbin pins first and then moving to the outer ones"). Nesting is what makes
//     rule 6 (no lead crossings) hold by construction: an outer winding's leads leave the coil
//     further out and land on pins outside every inner winding's pins.
//   * Every start sits towards the row's START end and every finish towards its FINISH end
//     (Wuerth, 'Transformers - what characteristics are important to you?', 2022-08-18, 15:22:
//     "horizontal packages: all the starts on one end of the rails and all the finishes on the
//     other end ... vertical package: starts in one diagonal corner"). Horizontal bobbins start
//     every row at -X; vertical bobbins start row 0 at -X and row 1 at +X, the diagonal. For
//     the innermost winding of a row this makes its start and finish adjacent pins (rule 2).
//   * Series sections put a tap pin between start and finish (rule 3); strands share a pin
//     while the wraps fit (rule 4, numbers in src/data/dfm_rules.json "pinAssignment"), else
//     one pin per strand, adjacent, with `parallel` set.
//   * TIW/insulated and margin-wound windings never take a corner pin (rule 7).
//   * Shield windings (sections or layers of type shielding) terminate both ends on pins when a
//     spare pin exists, else the finish is a buried FLYING_LEAD (rule 5 as amended by WP8 R14).
//
// Deliberate limits, each of which throws or is reported rather than guessed:
//   * Round winding windows (toroids on a family `t` base, WP4 / ABT #1173) follow the toroid rule in
//     assign_pins_on_toroid_base: each winding takes the base side its half of the ring faces.
//   * Role tags (switched-node end, primary return shared with a shield) do not exist in MAS, so
//     rule 5's "switched node is the start" cannot be applied; the first-wound end is the start.
//   * Removable pins are read but change nothing: unused pins are never counted as a creepage
//     bridge, so removing one cannot lengthen a path in this model.
//   * A tap pin is recorded as ONE connection (end: tap). RFC 0013 speaks of a tap connection on
//     both sections, but coil.json makes connections[] uniqueItems and the two would be identical.
//     The inter-section jumper itself is still routed inside the coil by the connection
//     reservation; routing it out through the tap pin is a follow-up.

#include "constructive_models/Coil.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Core.h"
#include "constructive_models/Insulation.h"
#include "physical_models/WireBend.h"
#include "support/Settings.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>
#include <numbers>
#include <cmrc/cmrc.hpp>
#include <iostream>
#include <magic_enum.hpp>
#include <iomanip>
#include <limits>
#include <map>
#include <array>
#include <unordered_map>
#include <bit>
#include <tuple>
#include <optional>
#include <set>
#include <mutex>
#include <sstream>

CMRC_DECLARE(dfmData);

namespace OpenMagnetics {

namespace {

constexpr double pinCoordinateTolerance = 1e-6;  // metres: two pins closer than this in their row coordinate share a row

std::string format_millimetres(double metres) {
    std::ostringstream stream;
    stream.setf(std::ios::fixed);
    stream.precision(3);
    stream << metres * 1000.0 << " mm";
    return stream.str();
}

const json& pin_assignment_rules() {
    static const json rules = [] {
        auto fs = cmrc::dfmData::get_filesystem();
        auto data = fs.open("src/data/dfm_rules.json");
        auto all = json::parse(std::string(data.begin(), data.end()));
        if (!all.contains("pinAssignment")) {
            throw std::runtime_error("src/data/dfm_rules.json has no pinAssignment block (RFC 0013 section 3 numbers)");
        }
        return all.at("pinAssignment");
    }();
    return rules;
}

double pin_assignment_number(const std::string& key) {
    const auto& rules = pin_assignment_rules();
    if (!rules.contains(key)) {
        throw std::runtime_error("src/data/dfm_rules.json pinAssignment has no number '" + key + "'");
    }
    return rules.at(key).get<double>();
}

bool is_pin_connection_type(const std::optional<ConnectionType>& type) {
    if (!type) {
        return true;  // untyped: the design did not say it is anything but a pin
    }
    return type.value() == ConnectionType::PIN || type.value() == ConnectionType::THT || type.value() == ConnectionType::SMT;
}

std::string end_name(End end) {
    switch (end) {
        case End::START: return "start";
        case End::FINISH: return "finish";
        case End::TAP: return "tap";
    }
    return "?";
}

double distance(const std::vector<double>& a, const std::vector<double>& b) {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

// Distance from a point to the core's bounding box, centred on the main column (0 inside it).
double distance_to_core_box(const std::vector<double>& point, double halfWidth, double halfHeight, double halfDepth) {
    const double dx = std::max(std::abs(point[0]) - halfWidth, 0.0);
    const double dy = std::max(std::abs(point[1]) - halfHeight, 0.0);
    const double dz = std::max(std::abs(point[2]) - halfDepth, 0.0);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

bool is_insulated_wire(Wire& wire) {
    auto coating = wire.resolve_coating();
    if (!coating || !coating->get_type()) {
        return false;
    }
    return coating->get_type().value() == InsulationWireCoatingType::INSULATED ||
           coating->get_type().value() == InsulationWireCoatingType::EXTRUDED;
}

// One connection a winding needs.
struct PinSlot {
    End end = End::START;
    std::optional<int64_t> parallel;  // absent: every strand of the winding terminates here
    size_t tapIndex = 0;              // 1..K-1 for taps (the junction after section k)
    std::optional<size_t> connectionIndex;  // the design's own connection this slot fills, if any
    std::optional<std::string> pinName;
    bool userPin = false;             // the design gave this pinName
    bool offPin = false;              // this end does not go to a pin (design type, or buried shield end)
    bool startHalf = true;            // allocated on the start half of the block (starts and taps)
};

struct WindingPinPlan {
    size_t windingIndex = 0;
    std::string name;
    IsolationSide side = IsolationSide::PRIMARY;
    size_t buildRank = 0;
    size_t numberSections = 0;
    int64_t numberParallels = 1;
    bool separateStrands = false;
    bool avoidCorners = false;
    bool isShield = false;
    bool legacy = false;              // connections carry pinNames but no `end`: validated, not planned
    double firstTurnY = 0;
    // Toroid base only (ABT #1173), MKF ring frame (x, y); the base's rows stand at z = -+rowDistance/2 in
    // the MVB++ toroid frame, where a turn's y is drawn as z. meanTurnY is the y-centroid of the winding's
    // turns (which half of the ring it occupies), first/lastTurnX the x of parallel 0's first and last turn
    // (which end of a row its start and finish leads leave towards).
    double meanTurnY = 0;
    double firstTurnX = 0;
    double lastTurnX = 0;
    std::vector<ConnectionElement> originalConnections;
    std::vector<PinSlot> slots;
};

std::string side_name(IsolationSide side) {
    json j;
    to_json(j, side);
    return j.get<std::string>();
}

}  // namespace

std::vector<PlacedPin> Coil::place_pins(const std::vector<MAS::Pin>& pins) {
    std::vector<PlacedPin> placed;
    for (const auto& pin : pins) {
        PlacedPin placedPin;
        if (!pin.get_name()) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "A bobbin pin has no name, so no winding end can be assigned to it.");
        }
        placedPin.name = pin.get_name().value();
        if (!pin.get_coordinates() || pin.get_coordinates()->size() < 3) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "Bobbin pin '" + placedPin.name + "' has no 3D coordinates; pin assignment needs where the pin is.");
        }
        if (pin.get_dimensions().size() < 3) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "Bobbin pin '" + placedPin.name + "' carries " + std::to_string(pin.get_dimensions().size()) +
                " dimensions; pin assignment needs [diameter, depth, length].");
        }
        placedPin.centre = pin.get_coordinates().value();
        placedPin.diameter = pin.get_dimensions()[0];
        placedPin.length = pin.get_dimensions()[2];
        placedPin.removable = pin.get_removable().value_or(false);  // bobbin.json default: false
        placedPin.type = pin.get_type();

        // Which way the pin leaves: Bobbin::expand_pinout writes no rotation for a vertical pin
        // (along -Y) and [90, 0, 0] for a horizontal one (along -Z). Nothing else is understood.
        bool alongZ = false;
        if (pin.get_rotation()) {
            auto rotation = pin.get_rotation().value();
            bool allZero = true;
            for (double angle : rotation) {
                if (std::abs(angle) > 1e-9) {
                    allZero = false;
                }
            }
            bool quarterAboutX = rotation.size() >= 1 && std::abs(rotation[0] - 90.0) < 1e-9 &&
                                 (rotation.size() < 2 || std::abs(rotation[1]) < 1e-9) &&
                                 (rotation.size() < 3 || std::abs(rotation[2]) < 1e-9);
            if (quarterAboutX) {
                alongZ = true;
            }
            else if (!allZero) {
                std::string text;
                for (double angle : rotation) {
                    text += (text.empty() ? "" : ", ") + std::to_string(angle);
                }
                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                    "Bobbin pin '" + placedPin.name + "' has rotation [" + text + "]. Pin assignment understands the "
                    "two placements Bobbin::expand_pinout produces: vertical (no rotation, along -Y) and horizontal "
                    "([90, 0, 0], along -Z).");
            }
        }
        placedPin.hangsAlongZ = alongZ;
        placedPin.base = placedPin.centre;
        if (alongZ) {
            placedPin.base[2] += placedPin.length / 2;
        }
        else {
            placedPin.base[1] += placedPin.length / 2;
        }
        placed.push_back(placedPin);
    }
    if (placed.empty()) {
        return placed;
    }
    for (const auto& placedPin : placed) {
        if (placedPin.hangsAlongZ != placed.front().hangsAlongZ) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "Bobbin pins '" + placed.front().name + "' and '" + placedPin.name + "' leave in different directions; "
                "a former whose pins are part vertical and part horizontal is not a footprint pin assignment understands.");
        }
    }

    // Rows: a vertical pin's row is its z, a horizontal pin's its y (Bobbin::expand_pinout).
    auto rowCoordinate = [](const PlacedPin& p) { return p.hangsAlongZ ? p.centre[1] : p.centre[2]; };
    std::vector<double> rowKeys;
    for (const auto& placedPin : placed) {
        double key = rowCoordinate(placedPin);
        bool known = false;
        for (double existing : rowKeys) {
            if (std::abs(existing - key) < pinCoordinateTolerance) {
                known = true;
            }
        }
        if (!known) {
            rowKeys.push_back(key);
        }
    }
    std::sort(rowKeys.begin(), rowKeys.end());
    for (auto& placedPin : placed) {
        double key = rowCoordinate(placedPin);
        for (size_t rowIndex = 0; rowIndex < rowKeys.size(); ++rowIndex) {
            if (std::abs(rowKeys[rowIndex] - key) < pinCoordinateTolerance) {
                placedPin.row = rowIndex;
            }
        }
    }
    for (size_t rowIndex = 0; rowIndex < rowKeys.size(); ++rowIndex) {
        std::vector<PlacedPin*> row;
        for (auto& placedPin : placed) {
            if (placedPin.row == rowIndex) {
                row.push_back(&placedPin);
            }
        }
        std::sort(row.begin(), row.end(), [](const PlacedPin* a, const PlacedPin* b) { return a->centre[0] < b->centre[0]; });
        for (size_t index = 1; index < row.size(); ++index) {
            // Two pins of one row closer than their radii overlap: that is not a footprint, and an
            // end assigned to either would share copper with the other.
            const double gap = row[index]->centre[0] - row[index - 1]->centre[0];
            if (gap < (row[index]->diameter + row[index - 1]->diameter) / 2) {
                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                    "Bobbin pins '" + row[index - 1]->name + "' and '" + row[index]->name + "' overlap: their centres are " +
                    format_millimetres(gap) + " apart along their row, less than their radii. A pinout with centralPitch 0 "
                    "puts the two middle pins on top of each other.");
            }
        }
        for (size_t index = 0; index < row.size(); ++index) {
            row[index]->indexAlongRow = index;
            row[index]->corner = (index == 0 || index + 1 == row.size());
        }
    }
    std::stable_sort(placed.begin(), placed.end(), [](const PlacedPin& a, const PlacedPin& b) {
        if (a.row != b.row) {
            return a.row < b.row;
        }
        return a.indexAlongRow < b.indexAlongRow;
    });
    return placed;
}

namespace {

// Closest distance between segments [p0, p1] and [q0, q1] (3D). Degenerate segments are points.
double segment_distance(const std::vector<double>& p0, const std::vector<double>& p1,
                        const std::vector<double>& q0, const std::vector<double>& q1) {
    double u[3], v[3], w[3];
    for (int k = 0; k < 3; ++k) {
        u[k] = p1[k] - p0[k];
        v[k] = q1[k] - q0[k];
        w[k] = p0[k] - q0[k];
    }
    auto dot = [](const double* x, const double* y) { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; };
    const double a = dot(u, u);
    const double b = dot(u, v);
    const double c = dot(v, v);
    const double d = dot(u, w);
    const double e = dot(v, w);
    const double denominator = a * c - b * b;
    double sN, sD = denominator, tN, tD = denominator;
    const double small = 1e-30;
    if (a < small && c < small) {
        return std::sqrt(dot(w, w));
    }
    if (a < small) {            // p is a point
        sN = 0; sD = 1; tN = e; tD = c;
    }
    else if (c < small) {       // q is a point
        tN = 0; tD = 1; sN = -d; sD = a;
    }
    else if (denominator < 1e-12 * a * c) {   // parallel
        sN = 0; sD = 1; tN = e; tD = c;
    }
    else {
        sN = b * e - c * d;
        tN = a * e - b * d;
    }
    if (sN < 0) {
        sN = 0; tN = e; tD = c;
    }
    else if (sN > sD) {
        sN = sD; tN = e + b; tD = c;
    }
    if (tN < 0) {
        tN = 0;
        if (-d < 0) { sN = 0; }
        else if (-d > a) { sN = sD; }
        else { sN = -d; sD = a; }
    }
    else if (tN > tD) {
        tN = tD;
        if ((-d + b) < 0) { sN = 0; }
        else if ((-d + b) > a) { sN = sD; }
        else { sN = -d + b; sD = a; }
    }
    const double sc = std::abs(sN) < small ? 0.0 : sN / sD;
    const double tc = std::abs(tN) < small ? 0.0 : tN / tD;
    double dp[3];
    for (int k = 0; k < 3; ++k) {
        dp[k] = w[k] + sc * u[k] - tc * v[k];
    }
    return std::sqrt(dot(dp, dp));
}

// Smallest distance between two polylines; `skipLastOfA` leaves A's last segment out.
double polyline_distance(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b,
                         bool skipLastOfA = false) {
    double best = std::numeric_limits<double>::max();
    const size_t segmentsA = a.size() < 2 ? 0 : a.size() - 1 - (skipLastOfA ? 1 : 0);
    for (size_t i = 0; i < segmentsA; ++i) {
        for (size_t j = 0; j + 1 < b.size(); ++j) {
            best = std::min(best, segment_distance(a[i], a[i + 1], b[j], b[j + 1]));
        }
    }
    return best;
}

bool is_terminal_route(const ConnectionRoute& route) {
    return route.kind == ConnectionKind::TERMINAL_ENTRANCE || route.kind == ConnectionKind::TERMINAL_EXIT;
}

std::string terminal_route_label(const ConnectionRoute& route) {
    return "winding '" + route.winding + "' parallel " + std::to_string(route.parallel) +
           (route.kind == ConnectionKind::TERMINAL_ENTRANCE ? " entrance" : " exit");
}

// A route's 2D polyline (radial, axial) lifted to 3D for the shared segment distance.
std::vector<std::vector<double>> route_polyline(const ConnectionRoute& route) {
    std::vector<std::vector<double>> points;
    for (const auto& waypoint : route.waypoints) {
        points.push_back({waypoint.at(0), waypoint.at(1), 0.0});
    }
    return points;
}

// What a copper piece occupies ALONG the connection face: [low, high] in x. A lead occupies its slot;
// a RAMPED lead (ABT #1336, ConnectionRoute::rampLength, real winding only) also occupies the stretch
// of its turn it climbs along -- an entrance towards lower x (the wrap leaves the crossing towards -x
// and the lead joins it at its slot, on x <= 0 under real winding, ABT #1705), an exit towards higher
// x (its last turn arrives at the slot from +x) -- the ramp plus the plannedBendRadius of level run its radial bend
// takes. Links and dragbacks sit at the crossing, x = 0.
std::pair<double, double> lane_footprint(const ConnectionRoute& route, double x) {
    if (!is_terminal_route(route)) {
        return {0.0, 0.0};
    }
    if (route.rampLength && route.laneStandoff) {
        throw std::logic_error("The terminal lead of " + terminal_route_label(route) + " has both a ramp and a lane standoff.");
    }
    // Spec section 18: an unramped lead's lane standoff is its own stretch the same way.
    const double span = route.rampLength ? route.rampLength.value() + route.plannedBendRadius
                                         : (route.laneStandoff ? route.laneStandoff.value() : 0.0);
    return route.kind == ConnectionKind::TERMINAL_ENTRANCE ? std::pair<double, double>{x - span, x}
                                                           : std::pair<double, double>{x, x + span};
}

double footprint_gap(const std::pair<double, double>& a, const std::pair<double, double>& b) {
    return std::max(a.first - b.second, b.first - a.second);
}

// Spec RW-size (Alf, 2026-10-25): the unit direction (radial, axial) of a terminal lead's FIRST leg in
// the winding window, from its attach waypoint (an entrance's last, an exit's first) to its neighbour.
std::pair<double, double> first_leg_direction(const ConnectionRoute& route) {
    if (route.waypoints.size() < 2) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "The terminal route of " + terminal_route_label(route) + " has no legs.");
    }
    const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
    // Spec section 21: an off-face exit's lift planned as a cosine S (ConnectionRoute::liftLength > 0) is
    // tangent to the turn axis where it leaves its turn, so its first leg is axial, towards the stub's
    // row, although waypoints 0 -> 1 run from the turn to the stub's x. A plain step (0) leaves radially.
    if (!entrance && route.liftLength && route.liftLength.value() > 0) {
        const double da = route.waypoints[1].at(1) - route.waypoints[0].at(1);
        if (!(std::abs(da) > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The off-face lift of the terminal route of " + terminal_route_label(route) + " climbs no height, so its direction is unknown.");
        }
        return {0.0, da > 0 ? 1.0 : -1.0};
    }
    const auto& attach = entrance ? route.waypoints.back() : route.waypoints.front();
    const auto& next = entrance ? route.waypoints[route.waypoints.size() - 2] : route.waypoints[1];
    const double dr = next.at(0) - attach.at(0);
    const double da = next.at(1) - attach.at(1);
    const double length = std::hypot(dr, da);
    if (!(length > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The first leg of the terminal route of " + terminal_route_label(route) + " has no length, so its direction is unknown.");
    }
    return {dr / length, da / length};
}

// Spec RW-size (Alf, 2026-10-25): what a terminal lead's copper occupies ALONG THE LANE (x on the
// connection face) and across its runs in the window section. The lead conductor's section is its
// radial width W by its axial height H on the turn (TerminalTurnPath::stubRadialExtent /
// stubAxialExtent; a round one a disc). A bend keeps the dimension perpendicular to its plane: a lead
// whose first leg is an axial STUB bends in the face plane, so H stands along x and W across every
// leg in the section; a lead that leaves RADIALLY at once, or tangentially off the turn, bends in the
// plane of the turn, so W stands along x and H across. A sleeve (ConnectionRoute::sleeveOuterDiameter)
// is a disc around the lead, so either extent is at least its outer diameter.
std::pair<double, double> lead_section_axes(const ConnectionRoute& route, const TerminalTurnPath& path, bool tangential) {
    const double width = path.stubRectangular ? path.stubRadialExtent : path.stubDiameter;
    const double height = path.stubRectangular ? path.stubAxialExtent : path.stubDiameter;
    if (!(width > 0) || !(height > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
            "The lead conductor of " + terminal_route_label(route) + " carries no coated size in its turn path (" +
            std::to_string(width * 1e3) + " x " + std::to_string(height * 1e3) + " mm).");
    }
    double alongLane = 0, across = 0;
    if (tangential) {
        alongLane = width;
        across = height;
    }
    else {
        const auto [dr, da] = first_leg_direction(route);
        alongLane = std::abs(dr) * width + std::abs(da) * height;
        across = std::abs(da) * width + std::abs(dr) * height;
    }
    return {alongLane, across};
}
// A lane is one straight run at its x from the attach turn to the window border, and a sleeved lead
// is sleeved from the border on (ABT #1269), so along the lane it holds its sleeve.
double lead_lane_extent(const ConnectionRoute& route, const TerminalTurnPath& path, bool tangential) {
    return std::max(lead_section_axes(route, path, tangential).first, route.sleeveOuterDiameter.value_or(0.0));
}
// Spec section 15 (Alf, 2026-10-26, "the terminal space should be reserved for the normal wire, without
// sleeve"): across its runs in the section the lead is its BARE coated section; the sleeve's own
// thickness applies only on the stretch the sleeve covers (sleeved_part, thick_lead_envelope).
double lead_bare_section_thickness(const ConnectionRoute& route, const TerminalTurnPath& path, bool tangential) {
    return lead_section_axes(route, path, tangential).second;
}
// Spec section 15: the sleeved stretch of a terminal lead's in-window polyline -- from its window-border
// end (an entrance's first point, an exit's last) for ConnectionRoute::sleeveSpan -- and the bare rest,
// each in path order. A span past the polyline's end sleeves all of it (the remainder lies on the attach
// turn). An unsleeved route is all bare. Throws on a sleeved route whose span was not planned (its
// sleeve states no overlapIntoWinding, so where it ends in the window is unknown).
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> split_at_sleeve_end(
    const ConnectionRoute& route, const std::vector<std::vector<double>>& points) {
    if (!route.sleeveOuterDiameter) {
        return {{}, points};
    }
    if (!route.sleeveSpan) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The terminal lead of " + terminal_route_label(route) + " is sleeved, but its sleeve states no overlapIntoWinding, "
            "so where the sleeve ends in the window is unknown.");
    }
    const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
    auto fromBorder = points;
    if (!entrance) {
        std::reverse(fromBorder.begin(), fromBorder.end());
    }
    std::vector<std::vector<double>> sleeved{fromBorder.front()};
    std::vector<std::vector<double>> bare;
    double left = route.sleeveSpan.value();
    size_t k = 0;
    for (; k + 1 < fromBorder.size(); ++k) {
        const auto& a = fromBorder[k];
        const auto& b = fromBorder[k + 1];
        const double length = std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
        if (length < left) {
            sleeved.push_back(b);
            left -= length;
            continue;
        }
        const double t = length > 0 ? left / length : 0.0;
        const std::vector<double> end{a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])};
        sleeved.push_back(end);
        bare.push_back(end);
        for (size_t m = k + 1; m < fromBorder.size(); ++m) {
            bare.push_back(fromBorder[m]);
        }
        break;
    }
    if (!entrance) {
        std::reverse(sleeved.begin(), sleeved.end());
        std::reverse(bare.begin(), bare.end());
    }
    return {sleeved, bare};
}
// The attach turn's own section in the window: half its radial width (the thickness of its run in the
// envelope) and how far its axial height reaches past that half width along the turn axis (0 for a
// round wire, whose disc the capsule covers).
std::pair<double, double> turn_section_half(const ConnectionRoute& route, const TerminalTurnPath& path) {
    const double width = path.turnRectangular ? path.turnRadialExtent : path.turnWireDiameter;
    const double height = path.turnRectangular ? path.turnAxialExtent : path.turnWireDiameter;
    if (!(width > 0) || !(height > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
            "The attach turn of " + terminal_route_label(route) + " carries no coated size in its turn path (" +
            std::to_string(width * 1e3) + " x " + std::to_string(height * 1e3) + " mm).");
    }
    return {width / 2, std::max(0.0, (height - width) / 2)};
}
// The support of a section centred on the origin along unit direction (ur, ua): a rectangle W x H
// (radial x axial) reaches |ur| W/2 + |ua| H/2, a disc its radius.
double section_support(bool rectangular, double width, double height, double diameter, double ur, double ua) {
    return rectangular ? (std::abs(ur) * width + std::abs(ua) * height) / 2 : diameter / 2;
}
// ABT #1705 (Alf, 2026-10-24, "lanes only where leads can meet"): the section envelope of a terminal
// lead's DRAWN path (Coil::terminal_leads_can_meet): its route polyline (radial, axial) from its attach
// turn to its window exit -- the stub its ramp climbs included -- plus, at every corner where the
// lead turns inside the section, the arc of its plannedBendRadius it is drawn with there (tangent to
// both legs; it cuts the corner, so it can come nearer to a lead inside that corner than the
// polyline does). The bend off the helix lies in the plane of the face and the lead's first leg, so
// in the section it falls on that leg. Points are 3D (z = 0) for the shared segment distance.
// Piece 0 is the route polyline, piece 1 the attach turn's run (extended by `turnReachPast` at both
// ends: a rectangular turn's axial height past its half width), every further piece a bend arc.
std::vector<std::vector<std::vector<double>>> terminal_lead_envelope(const ConnectionRoute& route,
                                                                     const std::pair<double, double>& turnSpan,
                                                                     double turnReachPast = 0.0) {
    const auto points = route_polyline(route);   // a copy: `pieces` grows below
    std::vector<std::vector<std::vector<double>>> pieces{points};
    // The attach turn's own run past the lead, at its radius, over the heights it climbs there.
    const auto& attach = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? route.waypoints.back() : route.waypoints.front();
    pieces.push_back({{attach.at(0), turnSpan.first - turnReachPast, 0.0}, {attach.at(0), turnSpan.second + turnReachPast, 0.0}});
    const double radius = route.plannedBendRadius;
    for (size_t k = 1; k + 1 < points.size(); ++k) {
        const double ax = points[k][0] - points[k - 1][0], ay = points[k][1] - points[k - 1][1];
        const double bx = points[k + 1][0] - points[k][0], by = points[k + 1][1] - points[k][1];
        const double la = std::hypot(ax, ay), lb = std::hypot(bx, by);
        if (!(la > 0) || !(lb > 0) || !(radius > 0)) {
            continue;
        }
        const double ux = ax / la, uy = ay / la, vx = bx / lb, vy = by / lb;
        const double turning = std::acos(std::clamp(ux * vx + uy * vy, -1.0, 1.0));   // 0: straight on
        if (turning < 1e-9 || turning > std::numbers::pi - 1e-9) {
            continue;
        }
        const double tangent = radius * std::tan(turning / 2);
        const double side = (ux * vy - uy * vx) > 0 ? 1.0 : -1.0;   // the corner turns left (+1) or right
        const double cx = points[k][0] - tangent * ux - side * radius * uy;
        const double cy = points[k][1] - tangent * uy + side * radius * ux;
        const double start = std::atan2(points[k][1] - tangent * uy - cy, points[k][0] - tangent * ux - cx);
        const size_t count = std::max<size_t>(8, size_t(std::ceil(turning / (std::numbers::pi / 64))));
        std::vector<std::vector<double>> arc;
        for (size_t step = 0; step <= count; ++step) {
            const double angle = start + side * turning * double(step) / double(count);
            arc.push_back({cx + radius * std::cos(angle), cy + radius * std::sin(angle), 0.0});
        }
        pieces.push_back(std::move(arc));
    }
    return pieces;
}
double envelope_distance(const std::vector<std::vector<std::vector<double>>>& a,
                         const std::vector<std::vector<std::vector<double>>>& b) {
    double best = std::numeric_limits<double>::max();
    for (const auto& pa : a) {
        for (const auto& pb : b) {
            best = std::min(best, polyline_distance(pa, pb));
        }
    }
    return best;
}
// One lead's envelope with the half thickness of each piece: `leadHalf` for its own copper (polyline
// and bend arcs), `turnHalf` for its attach turn's run (piece 1), and, spec section 15, `sleeveHalf` for
// the pieces its sleeve covers (`sleeved`, per piece).
struct ThickEnvelope {
    std::vector<std::vector<std::vector<double>>> pieces;
    double leadHalf = 0;
    double turnHalf = 0;
    double sleeveHalf = 0;
    std::vector<bool> sleeved;   // per piece; empty or short = bare
    double half(size_t piece) const {
        if (piece == 1) {
            return turnHalf;
        }
        return piece < sleeved.size() && sleeved[piece] ? sleeveHalf : leadHalf;
    }
};
// Two envelopes meet when some pair of their pieces comes within the two pieces' half thicknesses,
// an exact touch included.
bool envelopes_meet(const ThickEnvelope& a, const ThickEnvelope& b) {
    for (size_t i = 0; i < a.pieces.size(); ++i) {
        for (size_t j = 0; j < b.pieces.size(); ++j) {
            const double reach = a.half(i) + b.half(j);
            if (polyline_distance(a.pieces[i], b.pieces[j]) <= reach + 1e-9) {
                return true;
            }
        }
    }
    return false;
}
// Spec RW-size: the envelope of a terminal lead with its real sections (round leads: their discs).
// Spec section 15 (Alf, 2026-10-26): the lead is its bare section inside the winding and its sleeve
// only where the sleeve covers it (ConnectionRoute::sleeveSpan from the window border): piece 0 keeps
// the bare stretch of the polyline, and the sleeved stretch is appended as its own piece at the
// sleeve's outer radius. A bend arc is the lead's bare copper (the bends sit at the turn and in the
// window; a bend arc whose corner lies on the sleeved stretch is the sleeve's, else the bare copper).
ThickEnvelope thick_lead_envelope(const ConnectionRoute& route, const TerminalTurnPath& path, bool tangential,
                                  const std::pair<double, double>& turnSpan) {
    const auto [turnHalf, turnReachPast] = turn_section_half(route, path);
    ThickEnvelope envelope{terminal_lead_envelope(route, turnSpan, turnReachPast), lead_bare_section_thickness(route, path, tangential) / 2,
                           turnHalf};
    if (route.sleeveOuterDiameter) {
        const auto points = envelope.pieces.front();
        auto [sleevedPart, barePart] = split_at_sleeve_end(route, points);
        envelope.sleeved.assign(envelope.pieces.size(), false);
        // Arc pieces (2, 3, ...) follow the polyline's inner corners in order (terminal_lead_envelope
        // skips a straight-on corner, so each kept arc is matched to its corner by turning).
        const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        std::vector<double> fromBorder(points.size(), 0.0);   // arc length of each point from the border end
        for (size_t k = 1; k < points.size(); ++k) {
            const size_t a = entrance ? k - 1 : points.size() - k;
            const size_t b = entrance ? k : points.size() - k - 1;
            fromBorder[b] = fromBorder[a] + std::hypot(points[b][0] - points[a][0], points[b][1] - points[a][1]);
        }
        size_t arc = 2;
        for (size_t k = 1; k + 1 < points.size() && arc < envelope.pieces.size(); ++k) {
            const double ax = points[k][0] - points[k - 1][0], ay = points[k][1] - points[k - 1][1];
            const double bx = points[k + 1][0] - points[k][0], by = points[k + 1][1] - points[k][1];
            const double la = std::hypot(ax, ay), lb = std::hypot(bx, by);
            if (!(la > 0) || !(lb > 0) || !(route.plannedBendRadius > 0)) {
                continue;
            }
            const double turning = std::acos(std::clamp((ax * bx + ay * by) / (la * lb), -1.0, 1.0));
            if (turning < 1e-9 || turning > std::numbers::pi - 1e-9) {
                continue;
            }
            envelope.sleeved[arc++] = fromBorder[k] < route.sleeveSpan.value();
        }
        envelope.pieces.front() = barePart.size() >= 2 ? barePart : std::vector<std::vector<double>>{sleevedPart.back(), sleevedPart.back()};
        envelope.pieces.push_back(sleevedPart);
        envelope.sleeved.push_back(true);
        envelope.sleeveHalf = std::max(envelope.leadHalf, route.sleeveOuterDiameter.value() / 2);
    }
    return envelope;
}
bool same_bundle(const ConnectionRoute& a, const ConnectionRoute& b) {
    return a.winding == b.winding && a.kind == b.kind && a.side == b.side;
}

struct RoutedPinLead {
    size_t request = 0;
    size_t pinIndex = 0;
    double radius = 0;                              // coated radius of the lead
    double wrapRadius = 0;                          // pin radius + lead radius (wrap centreline)
    std::vector<std::vector<double>> points;        // window exit -> pin axis
    std::vector<std::vector<double>> wrapAxis;      // the wrap's extent along the pin axis, {top, bottom}
    double wrapTop = 0;                             // pin-axis coordinate of the wrap's first turn
    double wrapBottom = 0;                          // ... and of its last
    std::optional<std::string> collision;           // spec section 17: why no candidate cleared (PIN_RUN_COLLIDES)
    double dropX = 0;                               // x of the leg along the pin axis
};

}  // namespace

namespace {

// Distance from segment [a, b] to an axis-aligned box (0 when they touch or cross). The distance to a
// convex set is convex along the segment, so a golden-section search finds its minimum.
double segment_box_distance(const std::vector<double>& a, const std::vector<double>& b, const Bobbin::PinRailBlock& box) {
    auto at = [&](double t) {
        double d2 = 0;
        for (int k = 0; k < 3; ++k) {
            const double p = a[k] + t * (b[k] - a[k]);
            const double excess = std::max(std::abs(p - box.centre[k]) - box.halfExtents[k], 0.0);
            d2 += excess * excess;
        }
        return std::sqrt(d2);
    };
    const double golden = (std::sqrt(5.0) - 1) / 2;
    double lo = 0, hi = 1;
    double t1 = hi - golden * (hi - lo), t2 = lo + golden * (hi - lo);
    double f1 = at(t1), f2 = at(t2);
    for (int iteration = 0; iteration < 90; ++iteration) {
        if (f1 < f2) {
            hi = t2; t2 = t1; f2 = f1; t1 = hi - golden * (hi - lo); f1 = at(t1);
        }
        else {
            lo = t1; t1 = t2; f1 = f2; t2 = lo + golden * (hi - lo); f2 = at(t2);
        }
    }
    return std::min({at(0.0), at(1.0), f1, f2});
}

}  // namespace

bool Coil::terminal_leads_can_meet(const ConnectionRoute& a, double diameterA, const std::pair<double, double>& turnSpanA,
                                   const ConnectionRoute& b, double diameterB, const std::pair<double, double>& turnSpanB) {
    for (const auto* route : {&a, &b}) {
        if (!is_terminal_route(*route) || route->waypoints.size() < 2) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "terminal_leads_can_meet needs two terminal routes with legs, got " + terminal_route_label(*route) + ".");
        }
    }
    return envelopes_meet({terminal_lead_envelope(a, turnSpanA), diameterA / 2, diameterA / 2},
                          {terminal_lead_envelope(b, turnSpanB), diameterB / 2, diameterB / 2});
}

bool Coil::terminal_leads_meet_as_planned(const ConnectionRoute& a, const TerminalTurnPath& pathA,
                                         const ConnectionRoute& b, const TerminalTurnPath& pathB) {
    for (const auto* route : {&a, &b}) {
        if (!is_terminal_route(*route) || route->waypoints.size() < 2) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "terminal_leads_meet_as_planned needs two terminal routes with legs, got " + terminal_route_label(*route) + ".");
        }
        if (!route->laneTurnSpan || !route->laneTangential) {
            throw InvalidInputException(ErrorCode::MISSING_DATA,
                "The terminal lead of " + terminal_route_label(*route) + " carries no lane planner record (laneTurnSpan, "
                "laneTangential): it was not laid by the real-winding lane planner.");
        }
    }
    return envelopes_meet(thick_lead_envelope(a, pathA, a.laneTangential.value(), a.laneTurnSpan.value()),
                          thick_lead_envelope(b, pathB, b.laneTangential.value(), b.laneTurnSpan.value()));
}

double Coil::terminal_turn_height_at(const ConnectionRoute& route, double advance, const TerminalTurnPath& path, double x) {
    if (!is_terminal_route(route) || route.waypoints.empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "terminal_turn_height_at needs a terminal route with waypoints, not " + terminal_route_label(route) + ".");
    }
    if (!std::isfinite(advance) || !std::isfinite(x)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The attach turn of " + terminal_route_label(route) + " has no finite advance or lane.");
    }
    if (!(path.length > 0) || !(path.faceHalfStraight >= 0) || !(path.cornerRadius > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The attach turn of " + terminal_route_label(route) + " has no path along its connection face (length " +
            std::to_string(path.length * 1e3) + " mm, straight " + std::to_string(path.faceHalfStraight * 1e3) + " mm, arc " +
            std::to_string(path.cornerRadius * 1e3) + " mm).");
    }
    const double reach = std::abs(x);
    if (!(reach <= path.faceHalfStraight + path.cornerRadius)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The lane x = " + std::to_string(x * 1e3) + " mm of " + terminal_route_label(route) + " lies beyond its turn's face, " +
            std::to_string((path.faceHalfStraight + path.cornerRadius) * 1e3) + " mm from the crossing, so its turn does not pass there.");
    }
    return terminal_turn_height_along(route, advance, path, terminal_lane_along(path, x));
}

double Coil::terminal_turn_height_along(const ConnectionRoute& route, double advance, const TerminalTurnPath& path, double along) {
    if (!is_terminal_route(route) || route.waypoints.empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "terminal_turn_height_along needs a terminal route with waypoints, not " + terminal_route_label(route) + ".");
    }
    if (!std::isfinite(advance) || !std::isfinite(along) || !(path.length > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The attach turn of " + terminal_route_label(route) + " has no finite advance, position or path length.");
    }
    // Spec section 17: a tangential lead planned past the far side of its turn (plan_terminal_tangential_lead,
    // reported there as CROSSES_FAR_SIDE) leaves it more than half a turn from the crossing; its height
    // there follows the same uniform climb.
    const auto& attach = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? route.waypoints.back() : route.waypoints.front();
    // Every turn leaves the crossing towards -x and climbs its advance uniformly over its whole path
    // (MVB++ appendRoundWrap / appendRectWrap), so `along` upstream of the crossing (> 0) it stands that
    // share of its advance below its station, downstream (< 0) above it. Entrances and exits alike.
    return attach.at(1) - advance * along / path.length;
}

TerminalTangentialLead Coil::plan_terminal_tangential_lead(const ConnectionRoute& route, const TerminalTurnPath& path, double lane) {
    const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
    if (!is_terminal_route(route)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "plan_terminal_tangential_lead needs a terminal route, not " + terminal_route_label(route) + ".");
    }
    if (!std::isfinite(lane) || lane == 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The lead of " + terminal_route_label(route) + " has no tangential lane off the crossing (lane x = " + std::to_string(lane * 1e3) +
            " mm).");
    }
    // Spec S0 (Alf, 2026-10-25): an exit at x < 0 or an entrance at x > 0 is accepted; it lies on the
    // FAR half, where its turn does not run towards the back of the column (TerminalTangentialLead).
    const bool farHalf = entrance ? lane > 0 : lane < 0;
    const double bend = route.plannedBendRadius;
    const double reach = std::abs(lane);
    if (!(bend > 0) || !(path.cornerRadius > 0) || !(path.faceHalfStraight >= 0) || !(path.lateralHalfStraight >= 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The lead of " + terminal_route_label(route) + " has no planned bend radius or its turn no path across the column.");
    }
    if (reach < path.laneBound - 1e-12) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The lead of " + terminal_route_label(route) + " at lane x = " + std::to_string(lane * 1e3) +
            " mm lies on its turn's face (half extent " + std::to_string(path.laneBound * 1e3) +
            " mm): it leaves the turn there, not tangentially.");
    }
    if (farHalf) {
        // The turn runs along the face straight past the crossing; the lead leaves level along the face
        // line at the straight's end (or earlier, when the lane is nearer than that end plus one bend
        // radius) and one 90 degree corner, turning away from the turn, puts it on its radial run.
        const double sign = lane < 0 ? -1.0 : 1.0;
        const double faceLine = -(path.lateralHalfStraight + path.rideOwnFace + path.cornerRadius);
        if (!(reach - bend >= 0)) {
            std::ostringstream out;
            out << std::fixed << std::setprecision(4) << "The lead of " << terminal_route_label(route) << " at lane x = " << lane * 1e3
                << " mm lies on the far half of the crossing for its kind and beyond its turn's half extent (" << path.laneBound * 1e3
                << " mm), so it leaves along the face line; its corner of radius " << bend * 1e3 << " mm would have to begin "
                << (bend - reach) * 1e3 << " mm past the crossing, on the other half.";
            throw CoilException(ErrorCode::COIL_WINDING_ERROR, out.str());
        }
        const double tangentReach = std::min(path.faceHalfStraight, reach - bend);
        TerminalTangentialLead lead;
        lead.lane = lane;
        lead.farHalf = true;
        lead.tangentPoint = {sign * tangentReach, faceLine};
        lead.cornerStart = {sign * (reach - bend), faceLine};
        lead.cornerCentre = {sign * (reach - bend), faceLine - bend};
        lead.cornerEnd = {lane, faceLine - bend};
        lead.cornerAngle = std::numbers::pi / 2;
        lead.straightLength = reach - bend - tangentReach;
        lead.along = terminal_lane_along(path, sign * tangentReach);
        return lead;
    }
    // Canonical frame: the exit's half (+x), the connection face at -z, the wire running clockwise
    // (towards -x on the face). The turn's back corner circle and the lead's corner circle are wrapped
    // the same way, so the straight is their direct common tangent on the outer side. A dragback ride
    // pushes the connection face out by rideOwnFace and the back face by rideBackFace (ride_at), which
    // lengthens both lateral straights.
    const double h = path.faceHalfStraight;
    const double rho = path.cornerRadius;
    const double backCentre = path.lateralHalfStraight + path.rideBackFace;
    const double lateral = 2 * path.lateralHalfStraight + path.rideOwnFace + path.rideBackFace;
    const double centreX = reach - bend;
    auto describe = [&]() {
        std::ostringstream out;
        out << std::fixed << std::setprecision(4) << "The lead of " << terminal_route_label(route) << " at lane x = " << lane * 1e3
            << " mm (corner radius " << bend * 1e3 << " mm; turn face straight " << h * 1e3 << " mm, lateral straight "
            << lateral * 1e3 << " mm, corner " << rho * 1e3 << " mm)";
        return out.str();
    };
    // The back-right corner carries the tangent for a straight leaving at up to 90 degrees off the side;
    // a corner larger than the turn's own reaches round to the back-left one (the same circle on a round
    // or oblong turn, whose back face has no straight).
    double theta = 0;
    double cornerX = h;
    for (double candidateX : {h, -h}) {
        const double dx = centreX - candidateX;
        const double dz = -backCentre;
        const double distance = std::hypot(dx, dz);
        if (!(distance >= std::abs(rho - bend) - 1e-12)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                describe() + " cannot leave its turn tangentially: its corner's centre on the mid plane lies " +
                std::to_string(distance * 1e3) + " mm from the centre of the turn's back corner, inside the " +
                std::to_string(std::abs(rho - bend) * 1e3) + " mm the two radii differ by, so no straight is tangent to both.");
        }
        theta = std::atan2(dz, dx) + std::acos(std::clamp((rho - bend) / std::max(distance, 1e-300), -1.0, 1.0));
        cornerX = candidateX;
        if (theta <= std::numbers::pi / 2 + 1e-12 || h == 0) {
            break;
        }
    }
    const bool backLeft = cornerX < h;
    if (theta < -1e-9 || !(theta < std::numbers::pi) || (backLeft && theta < std::numbers::pi / 2 - 1e-9)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            describe() + " finds no straight tangent to the back of its turn that a single corner turns onto its radial "
            "run (the tangent would stand at " + std::to_string(theta * 180 / std::numbers::pi) + " deg).");
    }
    theta = std::max(theta, 0.0);   // the lane at exactly the half extent: rounding only, checked above
    const double nx = std::cos(theta);
    const double nz = std::sin(theta);
    const double sign = entrance ? -1.0 : 1.0;
    TerminalTangentialLead lead;
    lead.lane = lane;
    lead.tangentPoint = {sign * (cornerX + rho * nx), backCentre + rho * nz};
    lead.cornerCentre = {sign * centreX, 0.0};
    lead.cornerStart = {sign * (centreX + bend * nx), bend * nz};
    lead.cornerEnd = {lane, 0.0};
    lead.cornerAngle = theta;
    lead.straightLength = std::hypot(lead.cornerStart[0] - lead.tangentPoint[0], lead.cornerStart[1] - lead.tangentPoint[1]);
    // Path length from the crossing: the face straight, the front corner, the lateral straight, then
    // round the back corner (and, past 90 degrees, along the back face to the back-left corner).
    const double along = backLeft ? h + rho * std::numbers::pi / 2 + lateral + rho * std::numbers::pi / 2 + 2 * h +
                                        rho * (theta - std::numbers::pi / 2)
                                  : h + rho * std::numbers::pi / 2 + lateral + rho * theta;
    if (!(along <= path.length / 2 + 1e-12)) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(4) << describe() << " cannot leave its turn tangentially: the straight tangent to "
            << "its corner meets the turn " << along * 1e3 << " mm from the crossing, past the far side of the turn (half its "
            << path.length * 1e3 << " mm). The lane lies " << (reach - path.laneBound) * 1e3 << " mm beyond the turn's half extent ("
            << path.laneBound * 1e3 << " mm).";
        // Spec section 17 (Alf, 2026-10-26): the lead is planned anyway -- its straight crosses the far side
        // of its own turn, as it must, and the turn is still wound -- and reported (CROSSES_FAR_SIDE).
        lead.crossesFarSide = true;
        lead.farSideWitness = out.str();
    }
    lead.along = sign * along;
    return lead;
}

int Coil::terminal_wrap_direction(const Turn& turn) {
    if (!turn.get_orientation()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "The turn '" + turn.get_name() + "' states no orientation, so the direction it travels along its connection face is unknown.");
    }
    // terminal_turn_height_at: every clockwise turn leaves the crossing towards -x (MVB++ appendRoundWrap /
    // appendRectWrap) in its face's own frame; a counter-clockwise one is that wrap reversed.
    return turn.get_orientation().value() == TurnOrientation::CLOCKWISE ? -1 : 1;
}

double Coil::terminal_lead_along(const ConnectionRoute& route, const TerminalTurnPath& path, double lane) {
    if (std::abs(lane) < path.laneBound - 1e-12) {   // the tangential test plan_terminal_tangential_departures applies
        return terminal_lane_along(path, lane);
    }
    return plan_terminal_tangential_lead(route, path, lane).along;
}

double Coil::terminal_lane_along(const TerminalTurnPath& path, double x) {
    const double reach = std::abs(x);
    if (!(reach <= path.faceHalfStraight + path.cornerRadius)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The lane x = " + std::to_string(x * 1e3) + " mm lies beyond its turn's face, " +
            std::to_string((path.faceHalfStraight + path.cornerRadius) * 1e3) + " mm from the crossing.");
    }
    // Path length from the crossing: the face straight, then the arc beyond it (a round turn is all arc).
    const double along = reach <= path.faceHalfStraight
                             ? reach
                             : path.faceHalfStraight + path.cornerRadius * std::asin((reach - path.faceHalfStraight) / path.cornerRadius);
    return x < 0 ? -along : along;
}

TerminalTurnPath Coil::terminal_turn_path(const ConnectionRoute& route, const Turn& attachTurn,
                                          const std::vector<ConnectionRoute>& routes) {
    if (!attachTurn.get_section()) {
        throw CoilException(ErrorCode::COIL_WINDING_ERROR,
            "The turn '" + attachTurn.get_name() + "' of " + terminal_route_label(route) + " has no section, so its frame is unknown.");
    }
    const auto frame = get_wound_column_frame_for_section(attachTurn.get_section().value());
    const double turnX = std::abs(attachTurn.get_coordinates()[0]);
    const double radius = frame.axisX == 0 ? turnX : std::abs(turnX - frame.axisX);
    TerminalTurnPath path;
    path.shape = frame.shape;
    // The wires' own coated sizes: the attach turn's, and its stub's -- the conductor the lead leaves by
    // (its stated lead conductor, else the winding's wire; ABT #1496), without any sleeve.
    const auto turnWire = resolve_wire(get_winding_index_by_name(route.winding));
    path.turnWireDiameter = std::max(turnWire.get_maximum_outer_width(), turnWire.get_maximum_outer_height());
    const auto stubWire = resolve_lead_conductor(route.winding,
                                                 route.kind == ConnectionKind::TERMINAL_ENTRANCE ? End::START : End::FINISH,
                                                 route.parallel).value_or(turnWire);
    path.stubDiameter = std::max(stubWire.get_maximum_outer_width(), stubWire.get_maximum_outer_height());
    // Spec RW-size (Alf, 2026-10-25): the REAL sections, width along the layer axis (radial) by height
    // along the turn axis (axial), the frame MKF lays every wire in; a round conductor is its disc.
    auto rectangular_section = [](const Wire& wire) {
        return wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL || wire.get_type() == WireType::PLANAR;
    };
    path.turnRectangular = rectangular_section(turnWire);
    path.turnRadialExtent = path.turnRectangular ? turnWire.get_maximum_outer_width() : path.turnWireDiameter;
    path.turnAxialExtent = path.turnRectangular ? turnWire.get_maximum_outer_height() : path.turnWireDiameter;
    path.stubRectangular = rectangular_section(stubWire);
    path.stubRadialExtent = path.stubRectangular ? stubWire.get_maximum_outer_width() : path.stubDiameter;
    path.stubAxialExtent = path.stubRectangular ? stubWire.get_maximum_outer_height() : path.stubDiameter;
    if (!(path.turnRadialExtent > 0) || !(path.turnAxialExtent > 0) || !(path.stubRadialExtent > 0) || !(path.stubAxialExtent > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
            "The turn '" + attachTurn.get_name() + "' of " + terminal_route_label(route) + " or its lead conductor has no coated outer size (turn " +
            std::to_string(path.turnRadialExtent * 1e3) + " x " + std::to_string(path.turnAxialExtent * 1e3) + " mm, lead " +
            std::to_string(path.stubRadialExtent * 1e3) + " x " + std::to_string(path.stubAxialExtent * 1e3) + " mm).");
    }
    if (frame.shape == ColumnShape::ROUND) {
        const auto& attach = route.kind == ConnectionKind::TERMINAL_ENTRANCE ? route.waypoints.back() : route.waypoints.front();
        path.cornerRadius = attach.at(0);
        path.faceHalfStraight = 0;
        path.lateralHalfStraight = 0;
        path.length = 2 * std::numbers::pi * path.cornerRadius;
        path.laneBound = path.cornerRadius;
        return path;
    }
    // The dragback rides lengthen both lateral faces (MVB++ appendRectWrap: 2 ride0 + 2 rideBack0).
    ConnectionLayout rides;
    rides.rideLevels = compute_ride_levels(routes);
    const double rideLength = 2 * rides.ride_at(turnX, route.side) + 2 * rides.ride_at(turnX, 1 - route.side);
    std::optional<double> length;
    if (frame.shape == ColumnShape::OBLONG) {
        length = get_turn_length_in_frame(frame, turnX);
        path.cornerRadius = radius;
        path.faceHalfStraight = 0;
        // get_turn_length_in_frame: 2 pi r + 4 (depth - width), the stadium's two lateral straights.
        path.lateralHalfStraight = frame.columnDepth - frame.columnWidth;
        path.laneBound = radius;
    }
    else if (frame.shape == ColumnShape::RECTANGULAR || frame.shape == ColumnShape::IRREGULAR) {
        const auto charged = get_turn_bend_radius(attachTurn.get_name());
        if (!charged) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR,
                "The turn '" + attachTurn.get_name() + "' of " + terminal_route_label(route) +
                " was charged no corner radius, so its path along the connection face is unknown.");
        }
        length = get_turn_length_in_frame(frame, turnX, charged.value());
        path.cornerRadius = charged.value();
        path.faceHalfStraight = frame.columnWidth + (radius - frame.columnWidth) - charged.value();
        path.lateralHalfStraight = frame.columnDepth + (radius - frame.columnWidth) - charged.value();
        // The turn's half extent along the face: a lane on the corner arc leaves the arc as on a round turn.
        path.laneBound = path.faceHalfStraight + path.cornerRadius;
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
            "The column of " + terminal_route_label(route) + " has a shape whose connection face is not modelled.");
    }
    if (!length || !(length.value() > 0) || !(path.faceHalfStraight >= 0) || !(path.lateralHalfStraight >= 0)) {
        throw CoilException(ErrorCode::COIL_WINDING_ERROR,
            "The turn '" + attachTurn.get_name() + "' of " + terminal_route_label(route) + " has no positive length or face straight.");
    }
    path.length = length.value() + rideLength;
    path.rideLength = rideLength;
    path.rideOwnFace = rides.ride_at(turnX, route.side);
    path.rideBackFace = rides.ride_at(turnX, 1 - route.side);
    return path;
}

double Coil::lead_sleeve_outer_diameter(const ConnectionSleeve& sleeve) {
    return sleeve.get_inner_diameter() + 2 * sleeve.get_wall_thickness() * double(sleeve.get_number_layers().value_or(1));
}

double Coil::terminal_lead_bend_radius(const std::string& windingName, End end, int64_t parallel) {
    // The wire this lead IS (ABT #1496): the stated lead conductor, else the winding's own.
    std::optional<Wire> leadConductor = resolve_lead_conductor(windingName, end, parallel);
    const Wire leadWire = leadConductor ? leadConductor.value() : resolve_wire(get_winding_index_by_name(windingName));
    const auto sleeve = resolve_lead_sleeve(windingName, end, parallel);
    double sweptDiameter = std::max(leadWire.get_maximum_outer_width(), leadWire.get_maximum_outer_height());
    if (sleeve) {
        sweptDiameter = std::max(sweptDiameter, lead_sleeve_outer_diameter(sleeve.value()));
    }
    if (!(sweptDiameter > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The " + std::string(end == End::START ? "start" : "finish") + " lead of winding '" + windingName + "' parallel " +
            std::to_string(parallel) + " has no outer dimension, so the radius it bends at cannot be planned.");
    }
    return lead_bend_radius(leadWire, sleeve, sweptDiameter / 2, is_real_winding_blocking_applied());
}

double Coil::lead_bend_radius(const Wire& wire, const std::optional<ConnectionSleeve>& sleeve, double sweptRadius,
                              bool realWinding) {
    double radius = lead_isotropic_bend_radius(sleeve, sweptRadius, realWinding);
    if (!realWinding) {
        return radius;   // an ideal winding keeps the geometry it always had (Alf, 2026-09-21)
    }

    // Physics of the wire: its insulation must survive the bend (IEC 60317-0-1 Table 6 /
    // IEC 60317-0-2 Table 6). A lead leaving for a pin turns in 3D, so a rectangular wire bends
    // on both of its axes somewhere along the run. Its FLATWISE limit rules the lead. Its EDGEWISE
    // one does not (ABT #2012, Alf 2026-10-23, as for the turns in ABT #2004): the edgewise mandrel
    // of Table 6 qualifies the coating, it is not a forming limit, so an edgewise lead bend below it
    // is drawn and analyze_results reports it (Coil::judge_leads_against_wire_flexibility).
    std::vector<BendAxis> axes = wire.get_type() == WireType::RECTANGULAR
                                     ? std::vector<BendAxis>{BendAxis::FLATWISE}
                                     : std::vector<BendAxis>{BendAxis::ROUND};
    for (auto axis : axes) {
        if (auto wireMinimum = WireBend::get_flexibility_bend_radius_if_standardised(wire, axis)) {
            radius = std::max(radius, wireMinimum.value());
        }
    }
    return radius;
}

double Coil::lead_isotropic_bend_radius(const std::optional<ConnectionSleeve>& sleeve, double sweptRadius,
                                        bool realWinding) {
    // Buildability: what the consumer drawing the corner declared, on the radius it sweeps.
    // Throws on a non-positive swept radius.
    double radius = Settings::resolve_lead_bend_radius(sweptRadius);
    if (!realWinding) {
        return radius;
    }

    // Physics of the sleeve: the tubing's own rated minimum, when its material sources one.
    if (sleeve) {
        const auto& materialUnion = sleeve->get_material();
        InsulationMaterial material = std::holds_alternative<MAS::InsulationMaterial>(materialUnion)
            ? InsulationMaterial(std::get<MAS::InsulationMaterial>(materialUnion))
            : find_insulation_material_by_name(std::get<std::string>(materialUnion));
        if (material.get_minimum_bend_radius()) {
            // The tubing actually fitted is the smallest RATED size that holds this sleeve, so
            // that is the size whose rating applies -- the standards' own "next larger size" rule.
            const double innerDiameter = sleeve->get_inner_diameter();
            const double wall = sleeve->get_wall_thickness();
            std::optional<MAS::MinimumBendRadiusElement> rated;
            const auto ratedPoints = material.get_minimum_bend_radius().value();  // returned by value
            for (const auto& point : ratedPoints) {
                if (point.get_inner_diameter() < innerDiameter - 1e-12 || point.get_wall_thickness() < wall - 1e-12) {
                    continue;
                }
                if (!rated || point.get_inner_diameter() < rated->get_inner_diameter() - 1e-12 ||
                    (std::abs(point.get_inner_diameter() - rated->get_inner_diameter()) <= 1e-12 &&
                     point.get_wall_thickness() < rated->get_wall_thickness())) {
                    rated = point;
                }
            }
            if (!rated) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Sleeve material '" + material.get_name() + "' rates its minimum bend radius, but for no tubing with an inner diameter of at least " +
                    std::to_string(innerDiameter * 1e3) + " mm and a wall of at least " + std::to_string(wall * 1e3) +
                    " mm, which is the sleeve this lead needs");
            }
            radius = std::max(radius, rated->get_value());
        }
    }
    return radius;
}

double Coil::lead_ramp_governing_radius(const Wire& leadWire, double isotropicRadius, double surfaceCurvature) {
    const double inSurfaceCurvature2 = 1.0 / (isotropicRadius * isotropicRadius) - surfaceCurvature * surfaceCurvature;
    if (!(inSurfaceCurvature2 > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "A ramp that must bend no tighter than " + std::to_string(isotropicRadius * 1e3) +
            " mm cannot climb a surface of curvature " + std::to_string(surfaceCurvature) + " 1/m");
    }
    double inSurfaceRadius = 1.0 / std::sqrt(inSurfaceCurvature2);
    // A rectangular wire's limits are per axis (IEC 60317-0-2 Table 6 tests flatwise and
    // edgewise separately). The column's curvature bends the turn about the AXIAL axis:
    // that is the turn's own bend, already wound, and the ramp does not add to it (on a
    // cylinder the normal curvature of any direction is at most 1/rho). The S bends it
    // about the RADIAL axis, with the axial dimension -- the wire's height in this frame
    // (the ramp is only planned for non-contiguous layers) -- in the bend plane. So the S
    // answers to that one axis's minimum alone, not to the tighter of the two (ABT #1481).
    // ABT #2012: an EDGEWISE S answers to no Table 6 minimum -- that mandrel qualifies the
    // coating, it is not a forming limit -- and analyze_results reports the lead instead.
    if (leadWire.get_type() == WireType::RECTANGULAR &&
        WireBend::axis_from_bend_plane_dimension(leadWire, false) != BendAxis::EDGEWISE) {
        const auto inSurfaceAxis = WireBend::axis_from_bend_plane_dimension(leadWire, false);
        if (auto axisMinimum = WireBend::get_flexibility_bend_radius_if_standardised(leadWire, inSurfaceAxis)) {
            inSurfaceRadius = std::max(inSurfaceRadius, axisMinimum.value());
        }
    }
    return inSurfaceRadius;
}

double Coil::lead_ramp_length(double stubHeight, double governingRadius) {
    return std::numbers::pi * std::sqrt(stubHeight * governingRadius / 2);
}

// The same segment-to-segment distance, for the closed-core lead exit
// (CoilLeadExit.cpp), which checks its runs against each other exactly as the pin runs are.
double lead_polyline_distance(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b) {
    return polyline_distance(a, b);
}

int64_t Coil::pin_wrap_turns() {
    static const int64_t wrapTurns = [] {
        auto fs = cmrc::dfmData::get_filesystem();
        auto data = fs.open("src/data/dfm_rules.json");
        auto all = json::parse(std::string(data.begin(), data.end()));
        if (!all.contains("R4") || !all.at("R4").contains("wrapTurns")) {
            throw std::runtime_error("src/data/dfm_rules.json has no R4.wrapTurns (turns of wire wrapped around a pin per wire end)");
        }
        return all.at("R4").at("wrapTurns").get<int64_t>();
    }();
    return wrapTurns;
}

std::vector<int> Coil::terminal_lane_halves(const std::vector<ConnectionRoute>& routes,
                                            const std::vector<std::pair<std::string, IsolationSide>>& windings, bool realWinding,
                                            const std::vector<int>& pinHalves) {
    if (!pinHalves.empty() && pinHalves.size() != routes.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "terminal_lane_halves needs one pin half per route: " + std::to_string(routes.size()) + " routes, " +
            std::to_string(pinHalves.size()) + " pin halves.");
    }
    if (!realWinding) {
        for (int half : pinHalves) {
            if (half != 0) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "terminal_lane_halves was given the halves of pinned leads for an ideal winding, whose lanes all pack on x >= 0.");
            }
        }
    }
    // Spec section 12 (P-half): on a pinned face a lead takes the half its pin lies on, over the S0/S1 default.
    auto with_pin_halves = [&](std::vector<int> halves) {
        for (size_t index = 0; index < pinHalves.size(); ++index) {
            if (pinHalves[index] != 0 && is_terminal_route(routes[index])) {
                halves[index] = pinHalves[index];
            }
        }
        return halves;
    };
    std::vector<int> halves(routes.size(), 0);
    bool anyTerminal = false;
    for (const auto& route : routes) {
        anyTerminal = anyTerminal || is_terminal_route(route);
    }
    if (!anyTerminal) {
        return halves;
    }
    if (!realWinding) {
        for (size_t index = 0; index < routes.size(); ++index) {
            if (is_terminal_route(routes[index])) {
                halves[index] = 1;
            }
        }
        return halves;
    }
    if (windings.empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The terminal lanes under real winding need the magnetic's windings and their isolation sides, to put each winding's "
            "terminals on its half of the connection face; none were given.");
    }
    auto side_of = [&](const ConnectionRoute& route) {
        for (const auto& [name, side] : windings) {
            if (name == route.winding) {
                return side;
            }
        }
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The terminal of " + terminal_route_label(route) + " belongs to no winding given with its isolation side, so the half of "
            "the connection face it packs on is unknown.");
    };
    if (windings.size() == 1) {
        // Spec S1, an inductor: entrances on x <= 0, exits on x >= 0.
        for (size_t index = 0; index < routes.size(); ++index) {
            if (is_terminal_route(routes[index])) {
                side_of(routes[index]);
                halves[index] = routes[index].kind == ConnectionKind::TERMINAL_ENTRANCE ? -1 : 1;
            }
        }
        return with_pin_halves(halves);
    }
    // Spec S0: the windings are isolated by half -- the primary's isolation group on x <= 0, every
    // other isolation group on x >= 0.
    bool primaryGiven = false;
    for (const auto& winding : windings) {
        primaryGiven = primaryGiven || winding.second == IsolationSide::PRIMARY;
    }
    if (!primaryGiven) {
        std::string names;
        for (const auto& [name, side] : windings) {
            names += (names.empty() ? "'" : ", '") + name + "' (" + std::string(magic_enum::enum_name(side)) + ")";
        }
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "A multi-winding magnetic puts the primary's isolation group on one half of the connection face and every other group on "
            "the other, but no winding is on the PRIMARY isolation side: " + names + ".");
    }
    // Spec section 24 (Alf, 2026-10-28): once section 23 has put the isolation groups on separate faces, a
    // face that carries the leads of only ONE isolation group lays them by the one-winding rule S1 --
    // entrances on x <= 0, exits on x >= 0 of that face's frame -- instead of S0's group halves. A magnetic
    // whose windings are all in one isolation group keeps S0 (section 24 is about groups moved apart).
    std::set<IsolationSide> magneticGroups;
    for (const auto& winding : windings) {
        magneticGroups.insert(winding.second);
    }
    std::map<int, std::set<IsolationSide>> groupsOnFace;
    for (const auto& route : routes) {
        if (is_terminal_route(route)) {
            groupsOnFace[route.side].insert(side_of(route));
        }
    }
    for (size_t index = 0; index < routes.size(); ++index) {
        if (is_terminal_route(routes[index])) {
            if (magneticGroups.size() > 1 && groupsOnFace.at(routes[index].side).size() == 1) {
                halves[index] = routes[index].kind == ConnectionKind::TERMINAL_ENTRANCE ? -1 : 1;
            }
            else {
                halves[index] = side_of(routes[index]) == IsolationSide::PRIMARY ? -1 : 1;
            }
        }
    }
    return with_pin_halves(halves);
}

double Coil::terminal_lead_length(const ConnectionRoute& route, double advance, const TerminalTurnPath& path, double lane) {
    if (!is_terminal_route(route) || route.waypoints.size() < 2) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "terminal_lead_length needs a terminal route with legs, not " + terminal_route_label(route) + ".");
    }
    const bool entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
    const bool tangential = std::abs(lane) >= path.laneBound - 1e-12;
    const double along = terminal_lead_along(route, path, lane);
    const double turnHeight = terminal_turn_height_along(route, advance, path, along);
    // The lead in the section from its turn outwards: a stub starts where the turn really is, a route
    // that leaves radially at once is joined by that axial step.
    std::vector<std::vector<double>> lead;
    for (size_t k = 0; k < route.waypoints.size(); ++k) {
        const auto& waypoint = route.waypoints[entrance ? route.waypoints.size() - 1 - k : k];
        lead.push_back({waypoint.at(0), waypoint.at(1)});
    }
    // Spec section 21: an off-face exit's lift (liftLength) leaves its turn where the turn really is, too.
    if (std::abs(lead[1][0] - lead[0][0]) <= 1e-12 || route.liftLength) {
        lead[0][1] = turnHeight;
    }
    else {
        lead.insert(lead.begin(), {lead[0][0], turnHeight});
    }
    double length = 0;
    for (size_t k = 0; k + 1 < lead.size(); ++k) {
        length += std::hypot(lead[k + 1][0] - lead[k][0], lead[k + 1][1] - lead[k][1]);
    }
    if (tangential) {
        const auto plan = plan_terminal_tangential_lead(route, path, lane);
        length += plan.straightLength + route.plannedBendRadius * plan.cornerAngle;
    }
    // The turn's own copper between the crossing and where the lead leaves it: an exit's last turn
    // arrives from +along, so it stops short of the crossing at along > 0 and runs on past it at
    // along < 0; an entrance's first turn the other way round.
    length += entrance ? along : -along;
    return length;
}

std::vector<std::optional<double>> Coil::terminal_exit_slots(const std::vector<ConnectionRoute>& routes,
                                                             const std::vector<double>& diameters,
                                                             const std::vector<double>& attachAxial,
                                                             bool realWinding,
                                                             const std::vector<double>& attachAdvance,
                                                             const std::vector<std::optional<TerminalTurnPath>>& turnPaths,
                                                             const std::vector<std::pair<std::string, IsolationSide>>& windings) {
    return terminal_lane_plan(routes, diameters, attachAxial, realWinding, attachAdvance, turnPaths, windings).slots;
}

TerminalLanePlan Coil::terminal_lane_plan(const std::vector<ConnectionRoute>& routes, const std::vector<double>& diameters,
                                          const std::vector<double>& attachAxial, bool realWinding,
                                          const std::vector<double>& attachAdvance,
                                          const std::vector<std::optional<TerminalTurnPath>>& turnPaths,
                                          const std::vector<std::pair<std::string, IsolationSide>>& windings,
                                          const std::vector<int>& pinHalves) {
    if (diameters.size() != routes.size() || attachAxial.size() != routes.size() || attachAdvance.size() != routes.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "terminal_exit_slots needs one diameter, one attach coordinate and one attach advance per route: " +
            std::to_string(routes.size()) + " routes, " + std::to_string(diameters.size()) + " diameters, " +
            std::to_string(attachAxial.size()) + " attach coordinates, " + std::to_string(attachAdvance.size()) + " attach advances.");
    }
    if (realWinding && turnPaths.size() != routes.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "terminal_exit_slots under real winding needs one turn path per route, to know where each turn stands at the lanes: " +
            std::to_string(routes.size()) + " routes, " + std::to_string(turnPaths.size()) + " turn paths.");
    }
    auto is_terminal = is_terminal_route;
    auto label_of = terminal_route_label;
    auto polyline_of = route_polyline;
    for (size_t index = 0; index < routes.size(); ++index) {
        if (!is_terminal(routes[index])) {
            continue;
        }
        if (!(diameters[index] > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
                "The terminal lead of " + label_of(routes[index]) + " has no outer diameter, so its exit slot pitch is unknown.");
        }
        if (routes[index].waypoints.size() < 2) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "The terminal route of " + label_of(routes[index]) + " has no legs.");
        }
        if (!std::isfinite(attachAxial[index])) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The terminal route of " + label_of(routes[index]) + " has no attach turn coordinate.");
        }
        if (realWinding && !turnPaths[index]) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The terminal route of " + label_of(routes[index]) + " has no turn path, so where its turn stands at the lanes "
                "is unknown.");
        }
        if (realWinding && !std::isfinite(attachAdvance[index])) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The terminal route of " + label_of(routes[index]) + " has no attach turn advance, so the side of its "
                "lane its turn passes on is unknown.");
        }
        // ABT #1705: under real winding the lanes are planned for the lead's bends, so the bend must
        // be known. Coil::lead_bend_radius never returns less than the coated radius it sweeps.
        if (realWinding && !(routes[index].plannedBendRadius >= diameters[index] / 2 - 1e-12)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The terminal route of " + label_of(routes[index]) + " carries a planned bend radius of " +
                std::to_string(routes[index].plannedBendRadius * 1e3) + " mm, below its coated radius " +
                std::to_string(diameters[index] / 2 * 1e3) + " mm, so its lane cannot be planned for its bends.");
        }
    }
    TerminalLanePlan plan;
    plan.slots.assign(routes.size(), std::nullopt);
    plan.halves = terminal_lane_halves(routes, windings, realWinding, pinHalves);
    plan.bundleRanks.assign(routes.size(), std::nullopt);
    plan.findings.assign(routes.size(), {});
    plan.overflowed.assign(routes.size(), false);
    plan.turnSpans.assign(routes.size(), std::nullopt);
    plan.tangential.assign(routes.size(), std::nullopt);
    auto& slots = plan.slots;
    // Spec section 22 (Alf, 2026-10-27, "Overflow to spare half"): under real winding, in a multi-winding
    // magnetic (S0: one half per isolation group), a lead that does not fit its group's half continues on
    // the OTHER half, past the last lane of the group packed there. An inductor (S1) never overflows, and
    // spec section 12 rules over section 22: a winding whose pins decide its half keeps it.
    const bool overflowFace = realWinding && windings.size() > 1;
    auto winding_pinned = [&](const std::string& winding) {
        for (size_t index = 0; index < pinHalves.size(); ++index) {
            if (pinHalves[index] != 0 && is_terminal(routes[index]) && routes[index].winding == winding) {
                return true;
            }
        }
        return false;
    };
    const double tolerance = 1e-9;
    std::set<int> sides;
    for (const auto& route : routes) {
        if (is_terminal(route)) {
            sides.insert(route.side);
        }
    }
    for (int side : sides) {
        double pitch = 0;
        std::vector<std::string> windingOrder;
        for (size_t index = 0; index < routes.size(); ++index) {
            const auto& route = routes[index];
            if (route.side != side || !is_terminal(route)) {
                continue;
            }
            pitch = std::max(pitch, diameters[index]);
            if (std::find(windingOrder.begin(), windingOrder.end(), route.winding) == windingOrder.end()) {
                windingOrder.push_back(route.winding);
            }
        }
        // What a copper piece occupies ALONG the face: see lane_footprint.
        auto ramp_span = [&](size_t index) {
            const auto footprint = lane_footprint(routes[index], 0.0);
            return footprint.second - footprint.first;
        };
        // The stretch of its TURN a ramped lead climbs along (its lane standoff, section 18, is face, not turn).
        auto turn_climb = [&](size_t index) {
            return routes[index].rampLength ? routes[index].rampLength.value() + routes[index].plannedBendRadius : 0.0;
        };
        auto footprint_at = [&](size_t index, double x) { return lane_footprint(routes[index], x); };
        // Spec RW-size (Alf, 2026-10-25): under real winding a lead occupies its REAL extent along the
        // lane, which depends on whether it leaves on the face or tangentially (lead_lane_extent); the
        // ideal winding keeps the coated (or sleeve) diameter.
        // Found once per lead and tangential flag (the placement asks it for every order and anchor).
        std::vector<std::optional<double>> laneExtents(2 * routes.size());
        auto extent_of = [&](size_t index, bool tangential) {
            if (!realWinding) {
                return diameters[index];
            }
            auto& extent = laneExtents[2 * index + (tangential ? 1 : 0)];
            if (!extent) {
                extent = lead_lane_extent(routes[index], turnPaths[index].value(), tangential);
            }
            return extent.value();
        };
        auto beyond_half_extent = [&](size_t index, double lane) {
            return realWinding && std::abs(lane) >= turnPaths[index]->laneBound - 1e-12;
        };
        // ABT #1705 (Alf, 2026-10-24, "lanes only where leads can meet"): under real winding whether two
        // leads need disjoint lanes is decided by their envelopes (terminal_leads_can_meet, with each
        // piece's real thickness, spec RW-size). Each lead's attach turn runs past the other leads near
        // the crossing, so its envelope takes the heights that turn climbs (terminal_turn_height_along,
        // pitch-true) within `turnReach` of the crossing along it: as far along the face as another
        // lead, its bend and its ramp lie from its lane, plus the side's largest lane extent.
        std::map<size_t, std::pair<double, double>> turnSpan;
        if (realWinding) {
            double laneReach = 0;
            for (size_t index = 0; index < routes.size(); ++index) {
                if (routes[index].side == side && is_terminal(routes[index])) {
                    laneReach = std::max({laneReach, extent_of(index, false), extent_of(index, true)});
                }
            }
            double turnReach = laneReach;
            for (size_t index = 0; index < routes.size(); ++index) {
                if (routes[index].side == side && is_terminal(routes[index])) {
                    turnReach = std::max(turnReach, laneReach + routes[index].plannedBendRadius + turn_climb(index));
                }
            }
            for (size_t index = 0; index < routes.size(); ++index) {
                if (routes[index].side != side || !is_terminal(routes[index])) {
                    continue;
                }
                const auto& path = turnPaths[index].value();
                const double low = terminal_turn_height_along(routes[index], attachAdvance[index], path, -turnReach);
                const double high = terminal_turn_height_along(routes[index], attachAdvance[index], path, turnReach);
                turnSpan[index] = {std::min(low, high), std::max(low, high)};
                plan.turnSpans[index] = turnSpan[index];
            }
        }
        std::map<std::pair<size_t, bool>, ThickEnvelope> envelopes;
        auto envelope_of = [&](size_t index, bool tangential) -> const ThickEnvelope& {
            const auto key = std::make_pair(index, tangential);
            auto found = envelopes.find(key);
            if (found == envelopes.end()) {
                found = envelopes.emplace(key, thick_lead_envelope(routes[index], turnPaths[index].value(), tangential, turnSpan.at(index))).first;
            }
            return found->second;
        };
        // Whether two leads meet depends only on their two envelopes, so it is decided once per pair
        // (the placement asks it again for every candidate anchor and every round).
        // Indexed by (a, tangentialA, b, tangentialB) in the order asked: -1 not decided yet, else 0 / 1
        // (allocated at the first question: most faces never ask).
        std::vector<signed char> meetings;
        auto leads_can_meet = [&](size_t a, bool tangentialA, size_t b, bool tangentialB) {
            if (meetings.empty()) {
                meetings.assign(4 * routes.size() * routes.size(), -1);
            }
            auto& meet = meetings[(2 * a + (tangentialA ? 1 : 0)) * 2 * routes.size() + 2 * b + (tangentialB ? 1 : 0)];
            if (meet < 0) {
                meet = envelopes_meet(envelope_of(a, tangentialA), envelope_of(b, tangentialB)) ? 1 : 0;
            }
            return meet == 1;
        };
        struct Occupant {
            size_t route;
            std::pair<double, double> footprint;
            double extent;      // along the lane (real winding), else the diameter
            bool tangential;    // whether its lead leaves tangentially at its lane (real winding)
        };
        std::vector<Occupant> occupants;
        for (size_t index = 0; index < routes.size(); ++index) {
            const auto kind = routes[index].kind;
            if (routes[index].side == side &&
                (kind == ConnectionKind::Z_DRAGBACK || kind == ConnectionKind::EDGE_CONTINUATION ||
                 kind == ConnectionKind::U_ADJACENT || kind == ConnectionKind::U_TANGENTIAL)) {
                occupants.push_back({index, footprint_at(index, 0.0), diameters[index], false});
            }
        }
        // The distance between two routes' polylines depends only on the two routes, so it is measured
        // once per pair (the placement asks it again for every candidate anchor and every round).
        std::map<std::pair<size_t, size_t>, double> polylineDistances;
        auto polyline_distance_of = [&](size_t a, size_t b) {
            const auto key = std::make_pair(a, b);   // in the order asked: the measure is not reordered
            auto found = polylineDistances.find(key);
            if (found == polylineDistances.end()) {
                found = polylineDistances.emplace(key, polyline_distance(polyline_of(routes[a]), polyline_of(routes[b]))).first;
            }
            return found->second;
        };
        // The gap along the face that keeps a lead clear of an occupant, measured exactly (no
        // tolerance), or nothing when the two can never meet: a winding's own link, or sections
        // farther apart than the two can reach.
        auto required_gap = [&](const Occupant& occupant, size_t lead, bool leadTangential) -> std::optional<double> {
            const auto& other = routes[occupant.route];
            const bool link = other.kind == ConnectionKind::U_ADJACENT || other.kind == ConnectionKind::U_TANGENTIAL;
            if (link && other.winding == routes[lead].winding) {
                return std::nullopt;   // a winding's leads never yield to its own links
            }
            if (realWinding && is_terminal(other)) {
                // Rule 4 (ABT #1705; spec C3 and RW-size, Alf 2026-10-25): another bundle's lead that can
                // meet this one stands their two half extents along the lane off, and no more.
                // ABT #1336 (coordinator decision 2026-10-27, after the tie-lane change put 14_dab's
                // Secondary p3 exit on the Primary p3 exit's lane): a placed lead of ANOTHER winding always
                // counts, on either half and across the centre line; "lanes only where leads can meet"
                // spares only a winding's own leads.
                if (other.winding == routes[lead].winding && !leads_can_meet(lead, leadTangential, occupant.route, occupant.tangential)) {
                    return std::nullopt;
                }
                return (extent_of(lead, leadTangential) + occupant.extent) / 2;
            }
            const double requiredGap = (diameters[lead] + diameters[occupant.route]) / 2;   // along the face
            if (polyline_distance_of(lead, occupant.route) >= requiredGap - tolerance) {
                return std::nullopt;
            }
            return requiredGap;
        };
        // Two pieces clash when their sections come within reach AND their footprints along the face
        // stand closer than the required gap. Without ramps a footprint is a point and pitch >= the sum
        // of the radii, so this is exactly the old same-lane test.
        auto blocks_lead = [&](const Occupant& occupant, size_t lead, bool leadTangential, const std::pair<double, double>& footprint) {
            const auto requiredGap = required_gap(occupant, lead, leadTangential);
            return requiredGap && footprint_gap(footprint, occupant.footprint) < requiredGap.value() - tolerance;
        };
        struct BundlePlacement {
            std::vector<size_t> order;
            std::vector<double> lanes;
            std::vector<bool> tangential;
            std::optional<std::string> refusal;   // a tie no order clears (thrown when committed)
        };
        auto place_bundle = [&](std::vector<size_t> members, const std::vector<Occupant>& placed) -> BundlePlacement {
            const auto kind = routes[members.front()].kind;
            const std::string& winding = routes[members.front()].winding;
            std::stable_sort(members.begin(), members.end(), [&](size_t a, size_t b) { return routes[a].parallel < routes[b].parallel; });
            auto attach_radius = [&](size_t index) {
                double radius = std::numeric_limits<double>::max();
                for (const auto& waypoint : routes[index].waypoints) {
                    radius = std::min(radius, waypoint.at(0));
                }
                return radius;
            };
            auto axial_span = [&](size_t index) {
                double low = std::numeric_limits<double>::max();
                double high = std::numeric_limits<double>::lowest();
                for (const auto& waypoint : routes[index].waypoints) {
                    low = std::min(low, waypoint.at(1));
                    high = std::max(high, waypoint.at(1));
                }
                return high - low;
            };
            bool sameRadius = true;
            for (size_t member : members) {
                if (std::abs(attach_radius(member) - attach_radius(members.front())) > diameters[members.front()] / 2) {
                    sameRadius = false;
                }
            }
            const bool entrance = kind == ConnectionKind::TERMINAL_ENTRANCE;
            // The runs of members whose stubs TIE, as [begin, end) into `members` once sorted.
            std::vector<std::pair<size_t, size_t>> ties;
            if (sameRadius && members.size() > 1) {
                // Members stand in the order of their stub spans, and members whose stubs tie in
                // the order of their attach heights, ascending. Ideal winding keeps exactly that.
                std::stable_sort(members.begin(), members.end(), [&](size_t a, size_t b) {
                    const double spanA = axial_span(a);
                    const double spanB = axial_span(b);
                    if (std::abs(spanA - spanB) > 1e-12) {
                        return entrance ? spanA < spanB : spanA > spanB;
                    }
                    return attachAxial[a] < attachAxial[b];
                });
                if (realWinding) {
                    // Spec section 25 (Alf, 2026-10-28, "No first turn crosses"): under real winding the order
                    // along the face comes from each member's wrap direction (ConnectionRoute::wrapDirection,
                    // section 21, in the exitSlot frame). An entrance's first turn covers the face from its lane
                    // on in its wrap direction; an exit's last turn covers it from the far end up to its lane,
                    // i.e. the side against its wrap direction. A sibling whose stub spans farther than the
                    // member's (it climbs past the member's turn) must stand on the side that turn does not
                    // cover. In u = half * x (outwards from the crossing) that puts the larger-span sibling
                    // OUTWARDS of the member when the member's turn covers the inward side, else inwards.
                    // Members whose spans tie carry no constraint (their order is section 20's).
                    const int half = plan.halves[members.front()];
                    auto larger_goes_outwards = [&](size_t index) {
                        if (!routes[index].wrapDirection) {
                            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                                "The terminal lead of " + label_of(routes[index]) + " has no wrap direction, so the side of the "
                                "face its turn covers, and with it its bundle's order (spec section 25), is unknown.");
                        }
                        const int direction = routes[index].wrapDirection.value() * half;   // +1: the turn runs outwards
                        return entrance ? direction < 0 : direction > 0;
                    };
                    // must_precede(a, b): a stands nearer the crossing than b.
                    auto must_precede = [&](size_t a, size_t b) {
                        const double spanA = axial_span(a);
                        const double spanB = axial_span(b);
                        if (std::abs(spanA - spanB) <= 1e-12) {
                            return false;
                        }
                        const size_t smaller = spanA < spanB ? a : b;
                        const bool outwards = larger_goes_outwards(smaller);
                        return smaller == a ? outwards : !outwards;
                    };
                    // The constraints form a partial order when every member wraps the same way (a sort by
                    // span); a mixed bundle is ordered by them one placement at a time, keeping the order
                    // above among the free members, and throws when they cycle.
                    std::vector<size_t> remaining = members;
                    std::vector<size_t> ordered;
                    while (!remaining.empty()) {
                        auto free = std::find_if(remaining.begin(), remaining.end(), [&](size_t candidate) {
                            return std::none_of(remaining.begin(), remaining.end(),
                                                [&](size_t other) { return other != candidate && must_precede(other, candidate); });
                        });
                        if (free == remaining.end()) {
                            std::string names;
                            for (size_t index : remaining) {
                                names += (names.empty() ? "" : ", ") + label_of(routes[index]);
                            }
                            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                                "No order of the " + std::string(entrance ? "entrance" : "exit") + " leads " + names +
                                " keeps every first or last turn off its siblings' stubs (spec section 25): their wrap "
                                "directions ask for a cycle.");
                        }
                        ordered.push_back(*free);
                        remaining.erase(free);
                    }
                    members = ordered;
                }
                if (realWinding) {
                    for (size_t begin = 0; begin < members.size();) {
                        size_t end = begin + 1;
                        while (end < members.size() &&
                               std::abs(axial_span(members[end]) - axial_span(members[begin])) <= 1e-12) {
                            ++end;
                        }
                        if (end - begin > 1) {
                            ties.push_back({begin, end});
                        }
                        begin = end;
                    }
                }
            }
            // The lanes of the group's members in `order`: side by side from the anchor -- plus, between
            // two of them, whatever ramp lies in that gap. THE ANCHOR IS EXACT (ABT #1705, Alf: the lanes
            // are not snapped to the pitch grid). Each (member, occupant) pair that can meet blocks the
            // open interval of anchors at which the member's footprint stands closer than the required
            // gap to the occupant's; above it the pair is clear from the anchor at which the gap is exactly
            // met. The free anchors are therefore a finite union of intervals whose lowest point is 0 or
            // one of those edges: try them in order and take the first that clears every pair.
            // Spec S0/S1 (Alf, 2026-10-25): under real winding each bundle packs outwards from the crossing
            // on ITS HALF of the face (terminal_lane_halves): an inductor's entrances on x <= 0 and exits on
            // x >= 0, a multi-winding magnetic's primary isolation group on x <= 0 and every other group on
            // x >= 0. The packing works in u = mirror * x, outwards from the crossing. Ideal winding keeps
            // every lane at x >= 0.
            double mirror = 1.0;
            if (realWinding) {
                mirror = plan.halves[members.front()];
                for (size_t member : members) {
                    if (plan.halves[member] != plan.halves[members.front()]) {
                        // Spec section 12: only pins can split a bundle (S0/S1 give one half per winding and end).
                        auto pin_text = [&](size_t index) {
                            return label_of(routes[index]) + (routes[index].pinName.empty() ? std::string() : " (pin '" + routes[index].pinName + "')") +
                                   " on the " + (plan.halves[index] < 0 ? "x <= 0" : "x >= 0") + " half";
                        };
                        throw InvalidInputException(ErrorCode::INVALID_INPUT,
                            "The leads of one terminal bundle leave on both halves of face " + std::to_string(side) + ": " +
                            pin_text(members.front()) + ", " + pin_text(member) + ". A bundle's strands leave side by side on "
                            "one half, so their pins must lie on one side of the crossing.");
                    }
                }
            }
            auto to_u = [&](const std::pair<double, double>& footprint) {
                return mirror > 0 ? footprint : std::pair<double, double>{-footprint.second, -footprint.first};
            };
            // Spec section 22: the spare half is the other half of this face, worked in u' = -mirror * x,
            // outwards from the crossing. A member overflows only when its lane on its own half would
            // leave the face (beyond its turn's half extent); it then stands past EVERY lead already on
            // the spare half -- the other group's last lane and this bundle's earlier overflowed members --
            // by their two on-face extents along the lane (spec RW-size), clear of every placed piece
            // (blocks_lead), and is taken only if that lane stays on its turn's face.
            const bool overflowAllowed = overflowFace && !winding_pinned(winding);
            const double spareMirror = -mirror;
            auto to_spare_u = [&](const std::pair<double, double>& footprint) {
                return spareMirror > 0 ? footprint : std::pair<double, double>{-footprint.second, -footprint.first};
            };
            auto lane_of = [&](size_t index, const std::pair<double, double>& footprint) {
                return routes[index].kind == ConnectionKind::TERMINAL_ENTRANCE ? footprint.second : footprint.first;
            };
            // Returns the lane and whether it stays on the turn's face there.
            auto spare_half_lane = [&](const std::vector<size_t>& order, const std::vector<double>& lanes,
                                       const std::vector<bool>& flags, size_t position) -> std::pair<double, bool> {
                const size_t lead = order[position];
                const double leadExtent = extent_of(lead, false);
                const auto atCrossing = to_spare_u(footprint_at(lead, 0.0));   // the footprint at u' is atCrossing + u'
                double u = 0.0;
                auto stand_past = [&](const std::pair<double, double>& footprint, double extent) {
                    u = std::max(u, to_spare_u(footprint).second + (leadExtent + extent) / 2 - atCrossing.first);
                };
                for (const auto& occupant : placed) {
                    if (is_terminal(routes[occupant.route]) && spareMirror * lane_of(occupant.route, occupant.footprint) >= -tolerance) {
                        stand_past(occupant.footprint, occupant.extent);
                    }
                }
                // An earlier member on the spare side counts as a placed occupant does, the crossing (lane 0)
                // included: on a face whose spare half is empty the first overflowing member stands there.
                for (size_t earlier = 0; earlier < position; ++earlier) {
                    if (spareMirror * lanes[earlier] >= -tolerance) {
                        stand_past(footprint_at(order[earlier], lanes[earlier]), extent_of(order[earlier], flags[earlier]));
                    }
                }
                // Clear of every placed piece the lead can meet; each blocking edge moves it outwards, so this ends.
                while (true) {
                    const auto footprint = footprint_at(lead, spareMirror * u);
                    std::optional<double> edge;
                    for (const auto& occupant : placed) {
                        if (blocks_lead(occupant, lead, false, footprint)) {
                            const double past = to_spare_u(occupant.footprint).second + required_gap(occupant, lead, false).value() -
                                                atCrossing.first;
                            edge = std::max(edge.value_or(past), past);
                        }
                    }
                    if (!edge) {
                        break;
                    }
                    if (!(edge.value() > u)) {
                        throw std::logic_error("terminal_exit_slots: the spare-half lane of " + terminal_route_label(routes[lead]) +
                                               " is blocked although it already stands past the blocking edge.");
                    }
                    u = edge.value();
                }
                const double lane = spareMirror * u;
                return {lane, !beyond_half_extent(lead, lane)};
            };
            // A footprint's ramp lies OUTWARDS (towards larger u) for an entrance on x <= 0 and an exit on
            // x >= 0 (lane_footprint: an entrance's ramp runs towards lower x, an exit's towards higher x):
            // then the ramp between two members is the inner (earlier) one's, else the outer one's.
            const bool rampOutwards = entrance ? mirror < 0 : mirror > 0;
            // ABT #1705 (Alf, 2026-10-24, "lanes only where leads can meet"): under real winding a member
            // stands beyond every EARLIER member it can meet in 3D, by their two half extents along the
            // lane (spec RW-size) plus the ramp between, and shares the anchor's lane with members it
            // cannot meet. Ideal winding always stands them one pitch apart. Under real winding a lead's
            // extent and envelope depend on whether its lane leaves the turn's face (tangential), so the
            // placement is repeated until no lead newly lands beyond its turn's half extent (a lead once
            // found tangential stays so: the iteration is monotone and ends).
            // Spec section 20 (Alf, 2026-10-27: "Lanes from placed leads"): a member's lane depends only
            // on the leads placed before it -- the occupants of earlier bundles and the members earlier in
            // `order` -- never on members still to come. Member by member: it stands beyond every earlier
            // member it can meet by their spacing plus the ramp between (the earlier lanes are final), then
            // takes the lowest u from there at which it clears every placed occupant, and is found
            // tangential when that lane lies beyond its turn's half extent (then its extent and envelope
            // change, so its own lane is found again; a lead once found tangential stays so, so this ends
            // within two rounds). Its lane and flag are therefore final as soon as it is placed.
            // The lowest clear shift past an anchor depends only on the lead, its tangential flag, the
            // anchor `offset` and the occupants placed before this bundle (fixed for this call), so it is
            // found once per (lead, flag, exact offset) -- the tie-order search asks it again for every
            // order that shares those three.
            struct ShiftKeyHash {
                size_t operator()(const std::tuple<size_t, bool, uint64_t>& key) const {
                    uint64_t h = std::get<2>(key) * 0x9E3779B97F4A7C15ULL;
                    h ^= (static_cast<uint64_t>(std::get<0>(key)) << 1 | (std::get<1>(key) ? 1 : 0)) + 0x632BE59BD9B4E019ULL + (h << 6) + (h >> 2);
                    return static_cast<size_t>(h);
                }
            };
            std::unordered_map<std::tuple<size_t, bool, uint64_t>, double, ShiftKeyHash> clearShifts;
            auto place_member = [&](const std::vector<size_t>& order, std::vector<double>& lanes, std::vector<bool>& flags,
                                    size_t position) {
                const size_t lead = order[position];
                // Spec section 22: once decided (ConnectionRoute::laneOverflow, set by the first pass of the
                // lane settle) a lead that overflowed stays on the spare half -- leaving tangentially there if
                // its lane grows past its turn's face, drawn and reported -- and one that did not never does.
                const bool overflowDecided = routes[lead].laneOverflow.has_value();
                if (overflowAllowed && routes[lead].laneOverflow == true) {
                    const auto [spare, onFace] = spare_half_lane(order, lanes, flags, position);
                    lanes[position] = spare;
                    flags[position] = !onFace;
                    return;
                }
                bool flag = false;
                while (true) {
                    double offset = 0.0;   // in u, outwards from the crossing
                    for (size_t earlier = 0; earlier < position; ++earlier) {
                        if (mirror * lanes[earlier] < -tolerance) {
                            continue;   // spec section 22: it overflowed onto the spare half
                        }
                        if (realWinding && !leads_can_meet(lead, flag, order[earlier], flags[earlier])) {
                            continue;
                        }
                        const double spacing = realWinding ? (extent_of(lead, flag) + extent_of(order[earlier], flags[earlier])) / 2
                                                           : pitch;
                        const double between = rampOutwards ? ramp_span(order[earlier]) : ramp_span(lead);
                        offset = std::max(offset, mirror * lanes[earlier] + spacing + between);
                    }
                    const auto shiftKey = std::make_tuple(lead, flag, std::bit_cast<uint64_t>(offset));
                    auto knownShift = clearShifts.find(shiftKey);
                    if (knownShift == clearShifts.end()) {
                        // The free u at and beyond `offset` are a finite union of intervals whose lowest point is
                        // `offset` or one of the occupants' clearing edges: try them in order (ABT #1705, exact).
                        const auto footprintAtOffset = to_u(footprint_at(lead, mirror * offset));
                        std::vector<double> candidates{0.0};
                        for (const auto& occupant : placed) {
                            if (const auto requiredGap = required_gap(occupant, lead, flag)) {
                                const double edge = to_u(occupant.footprint).second + requiredGap.value() - footprintAtOffset.first;
                                if (edge > 0.0) {
                                    candidates.push_back(edge);
                                }
                            }
                        }
                        std::sort(candidates.begin(), candidates.end());
                        std::optional<double> shift;
                        for (double candidate : candidates) {
                            const auto footprint = footprint_at(lead, mirror * (offset + candidate));
                            bool free = true;
                            for (const auto& occupant : placed) {
                                if (blocks_lead(occupant, lead, flag, footprint)) {
                                    free = false;
                                    break;
                                }
                            }
                            if (free) {
                                shift = candidate;
                                break;
                            }
                        }
                        if (!shift) {
                            // The highest edge clears every occupant by construction, so this is a defect here.
                            throw std::logic_error("terminal_exit_slots: no lane clears " + terminal_route_label(routes[lead]) +
                                                   " of its bundle, although the highest blocking edge should.");
                        }
                        knownShift = clearShifts.emplace(shiftKey, shift.value()).first;
                    }
                    lanes[position] = mirror * (offset + knownShift->second);
                    if (!flag && beyond_half_extent(lead, lanes[position])) {
                        if (overflowAllowed && !overflowDecided) {
                            const auto [spare, onFace] = spare_half_lane(order, lanes, flags, position);
                            if (onFace) {
                                lanes[position] = spare;
                                flags[position] = false;
                                return;
                            }
                        }
                        flag = true;
                        continue;
                    }
                    flags[position] = flag;
                    return;
                }
            };
            auto place = [&](const std::vector<size_t>& order) {
                std::vector<double> lanes(order.size(), 0.0);
                std::vector<bool> flags(order.size(), false);
                for (size_t position = 0; position < order.size(); ++position) {
                    place_member(order, lanes, flags, position);
                }
                return std::make_pair(lanes, flags);
            };
            // ABT #1336 (Alf, 2026-10-19: "MKF plans pitch-true"): WHICH ORDER a tie takes is decided
            // by where each member's own turn really stands. An exit's last turn arrives at its lane
            // from upstream, so it passes the leads that leave upstream of its own; an entrance's first
            // turn runs on downstream, past those that leave downstream of it. At each sibling lead it
            // passes, the member's turn stands at its pitch-true height there (terminal_turn_height_along),
            // and it must stand clear of the sibling's lead: its stub from ITS turn's pitch-true height
            // there, and its run at its row. The clearance is the two sections' real reach towards each
            // other (spec RW-size): the member's turn section and, where the closest point is the
            // sibling's own attach point, its turn section, else its lead's section across its run --
            // the wires' own sizes, not a sleeve's (a sleeve starts beyond the stub). Returns the worst
            // pair of an order.
            struct TiePair {
                double margin;      // found - required: negative when the turn enters the sibling's envelope
                size_t member;      // the route whose turn passes
                size_t sibling;     // the route whose lane it passes
                double found;
                double required;
                double turnHeight;  // the member's turn at the sibling's lane
            };
            auto attach_waypoint_index = [&](size_t index) {
                return routes[index].kind == ConnectionKind::TERMINAL_ENTRANCE ? routes[index].waypoints.size() - 1 : size_t{0};
            };
            // Distance from a point to a polyline, the unit direction from the closest point to the point,
            // and whether the closest point is the polyline's first point.
            struct Closest {
                double distance;
                double ur;
                double ua;
                bool atStart;
            };
            auto point_to_polyline = [](const std::array<double, 2>& point, const std::vector<std::array<double, 2>>& polyline) {
                Closest best{std::numeric_limits<double>::max(), 1.0, 0.0, false};
                for (size_t k = 0; k + 1 < polyline.size(); ++k) {
                    const double ux = polyline[k + 1][0] - polyline[k][0];
                    const double uy = polyline[k + 1][1] - polyline[k][1];
                    const double length2 = ux * ux + uy * uy;
                    double t = 0.0;
                    if (length2 > 0.0) {
                        t = ((point[0] - polyline[k][0]) * ux + (point[1] - polyline[k][1]) * uy) / length2;
                        t = std::clamp(t, 0.0, 1.0);
                    }
                    const double dx = point[0] - polyline[k][0] - t * ux;
                    const double dy = point[1] - polyline[k][1] - t * uy;
                    const double distance = std::hypot(dx, dy);
                    if (distance < best.distance) {
                        best = {distance, distance > 0 ? dx / distance : 1.0, distance > 0 ? dy / distance : 0.0, k == 0 && t <= 0.0};
                    }
                }
                return best;
            };
            // One tie pair's score: the member's turn where it passes the sibling's lead against that lead,
            // from the two routes and where each leaves its turn (nothing when the member's turn does not
            // pass where the sibling leaves).
            std::vector<std::array<double, 2>> siblingLeadPoints;   // score_tie_pair's buffer, reused
            auto score_tie_pair = [&](size_t member, size_t sibling, double memberAlong, double siblingAlong) -> std::optional<TiePair> {
                const double memberRadius = routes[member].waypoints.at(attach_waypoint_index(member)).at(0);
                const auto& memberPath = turnPaths[member].value();
                const auto& siblingPath = turnPaths[sibling].value();
                // An exit's last turn arrives from upstream (+along), so it passes the leads
                // that leave further upstream; an entrance's first turn runs on downstream
                // (-along), past the leads that leave further downstream.
                if (!(entrance ? siblingAlong < memberAlong : siblingAlong > memberAlong)) {
                    return std::nullopt;   // the member's turn does not pass where this sibling leaves
                }
                const double turnHeight = terminal_turn_height_along(routes[member], attachAdvance[member], memberPath, siblingAlong);
                // The sibling's lead as drawn at its lane: it leaves ITS turn where that
                // turn really is there and runs along the turn axis to its route -- a
                // stub of its route then starts from that height, a route that leaves
                // radially at once is joined by that axial step.
                const auto& siblingWaypoints = routes[sibling].waypoints;
                if (siblingWaypoints.size() < 2) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "The terminal route of " + label_of(routes[sibling]) + " has fewer than two waypoints.");
                }
                auto& siblingLead = siblingLeadPoints;   // from its turn outwards
                siblingLead.clear();
                for (size_t k = 0; k < siblingWaypoints.size(); ++k) {
                    const auto& waypoint = siblingWaypoints[entrance ? siblingWaypoints.size() - 1 - k : k];
                    siblingLead.push_back({waypoint.at(0), waypoint.at(1)});
                }
                const double siblingTurnHeight =
                    terminal_turn_height_along(routes[sibling], attachAdvance[sibling], siblingPath, siblingAlong);
                // Spec section 21: an off-face exit's lift (liftLength) leaves its turn
                // there too, so it starts where its turn is.
                if (std::abs(siblingLead[1][0] - siblingLead[0][0]) <= 1e-12 || routes[sibling].liftLength) {
                    siblingLead[0][1] = siblingTurnHeight;   // its stub starts where its turn is
                }
                else {
                    const std::array<double, 2> stubStart{siblingLead[0][0], siblingTurnHeight};
                    siblingLead.insert(siblingLead.begin(), stubStart);
                }
                const auto closest = point_to_polyline({memberRadius, turnHeight}, siblingLead);
                if (!(memberPath.turnWireDiameter > 0) || !(siblingPath.stubDiameter > 0) || !(siblingPath.turnWireDiameter > 0)) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "The turn path of " + label_of(routes[member]) + " or " + label_of(routes[sibling]) +
                        " carries no wire size (TerminalTurnPath::turnWireDiameter / stubDiameter).");
                }
                // The member's TURN against the sibling's copper there (spec RW-size): each
                // section's reach along the line between them.
                const double memberReach = section_support(memberPath.turnRectangular, memberPath.turnRadialExtent,
                                                           memberPath.turnAxialExtent, memberPath.turnWireDiameter, closest.ur, closest.ua);
                double siblingReach = 0;
                if (closest.atStart) {
                    siblingReach = section_support(siblingPath.turnRectangular, siblingPath.turnRadialExtent, siblingPath.turnAxialExtent,
                                                   siblingPath.turnWireDiameter, closest.ur, closest.ua);
                }
                else if (siblingPath.stubRectangular) {
                    const auto [dr, da] = first_leg_direction(routes[sibling]);
                    siblingReach = (std::abs(da) * siblingPath.stubRadialExtent + std::abs(dr) * siblingPath.stubAxialExtent) / 2;
                }
                else {
                    siblingReach = siblingPath.stubDiameter / 2;
                }
                const double required = memberReach + siblingReach;
                return TiePair{closest.distance - required, member, sibling, closest.distance, required, turnHeight};
            };
            // Only the pairs whose two members are both among the first `placedCount` positions count.
            auto worst_tie_pair = [&](const std::vector<size_t>& order, const std::vector<double>& lanes, size_t placedCount) {
                std::optional<TiePair> worst;
                for (const auto& [begin, end] : ties) {
                    for (size_t p = begin; p < std::min(end, placedCount); ++p) {
                        const size_t member = order[p];
                        // ABT #1705: where along its turn each lead leaves -- at its lane on the face,
                        // or at its tangential lead's tangent point off it.
                        const double memberAlong = terminal_lead_along(routes[member], turnPaths[member].value(), lanes[p]);
                        for (size_t q = begin; q < std::min(end, placedCount); ++q) {
                            if (q == p) {
                                continue;
                            }
                            const size_t sibling = order[q];
                            const double siblingAlong = terminal_lead_along(routes[sibling], turnPaths[sibling].value(), lanes[q]);
                            const auto pair = score_tie_pair(member, sibling, memberAlong, siblingAlong);
                            if (pair && (!worst || pair->margin < worst->margin)) {
                                worst = pair;
                            }
                        }
                    }
                }
                return worst;
            };
            BundlePlacement result;
            result.order = members;
            std::tie(result.lanes, result.tangential) = place(result.order);
            if (!ties.empty()) {
                // WHICH order each tie takes (spec section 20, Alf 2026-10-27): the answer is the one the
                // odometer of permutations (one wheel per tie, wheel 0 turning fastest, every wheel starting
                // ascending, the ideal-winding order) would reach first -- its first order that clears, else
                // its first order of the largest worst margin. It is found by branch and bound over the
                // positions instead of by visiting all k! orders of every tie.
                // THE BOUND IS EXACT because a member's lane and tangential flag depend only on the leads
                // placed before it (place_member), and a tie pair's margin only on the two members' lanes
                // and their own routes and turns (worst_tie_pair). So once both members of a pair are placed,
                // no later choice changes that pair's margin, and the worst margin over the placed pairs is
                // an upper bound of the worst margin of every completion: a prefix with a pair that meets
                // has no completion that clears, and a prefix whose worst is already below the best
                // complete order's has no completion that beats it.
                std::vector<size_t> tieOf(members.size(), ties.size());
                for (size_t t = 0; t < ties.size(); ++t) {
                    for (size_t position = ties[t].first; position < ties[t].second; ++position) {
                        tieOf[position] = t;
                    }
                }
                // The odometer's rank of an order: its last wheel is the most significant.
                auto rank_key = [&](const std::vector<size_t>& chosen) {
                    std::vector<size_t> key;
                    for (size_t t = ties.size(); t-- > 0;) {
                        key.insert(key.end(), chosen.begin() + static_cast<std::ptrdiff_t>(ties[t].first),
                                   chosen.begin() + static_cast<std::ptrdiff_t>(ties[t].second));
                    }
                    return key;
                };
                struct Found {
                    std::vector<size_t> order;
                    std::vector<double> lanes;
                    std::vector<bool> tangential;
                    std::vector<size_t> key;
                    std::optional<TiePair> worst;
                };
                std::vector<size_t> candidate = members;
                std::vector<size_t> chosen(members.size(), 0);
                std::vector<bool> taken(members.size(), false);
                std::vector<double> lanes(members.size(), 0.0);
                std::vector<bool> flags(members.size(), false);
                std::optional<Found> best;
                bool clearing = true;   // the first pass looks for an order that clears, the second for the closest
                // `prefixWorst` is the worst margin of the pairs among the positions before `position` (none
                // when no pair has been found yet): a placed pair's margin never changes (see above), so each
                // position adds only the pairs it forms with the earlier members of its tie, from where each
                // placed member leaves its turn (kept per position: its lane is final while it stays placed).
                std::vector<double> placedAlongs(members.size(), 0.0);
                auto new_pairs_worst = [&](size_t position) {
                    const size_t begin = ties[tieOf[position]].first;
                    const size_t lead = candidate[position];
                    placedAlongs[position] = terminal_lead_along(routes[lead], turnPaths[lead].value(), lanes[position]);
                    std::optional<double> worst;
                    auto count = [&](const std::optional<TiePair>& pair) {
                        if (pair && (!worst || pair->margin < worst.value())) {
                            worst = pair->margin;
                        }
                    };
                    for (size_t earlier = begin; earlier < position; ++earlier) {
                        count(score_tie_pair(candidate[earlier], lead, placedAlongs[earlier], placedAlongs[position]));
                    }
                    for (size_t earlier = begin; earlier < position; ++earlier) {
                        count(score_tie_pair(lead, candidate[earlier], placedAlongs[position], placedAlongs[earlier]));
                    }
                    return worst;
                };
                std::function<void(size_t, std::optional<double>)> extend = [&](size_t position, std::optional<double> prefixWorst) {
                    if (position == members.size()) {
                        const auto worst = worst_tie_pair(candidate, lanes, members.size());
                        auto key = rank_key(chosen);
                        if (clearing) {
                            if (!best || key < best->key) {
                                best = Found{candidate, lanes, flags, std::move(key), worst};
                            }
                        }
                        else if (!best || worst->margin > best->worst->margin ||
                                 (worst->margin == best->worst->margin && key < best->key)) {
                            best = Found{candidate, lanes, flags, std::move(key), worst};
                        }
                        return;
                    }
                    const size_t t = tieOf[position];
                    if (t == ties.size()) {
                        candidate[position] = members[position];
                        place_member(candidate, lanes, flags, position);
                        extend(position + 1, prefixWorst);
                        return;
                    }
                    for (size_t k = 0; k < ties[t].second - ties[t].first; ++k) {
                        const size_t source = ties[t].first + k;
                        if (taken[source]) {
                            continue;
                        }
                        candidate[position] = members[source];
                        chosen[position] = k;
                        place_member(candidate, lanes, flags, position);
                        std::optional<double> worst = prefixWorst;
                        if (const auto added = new_pairs_worst(position)) {
                            if (!worst || added.value() < worst.value()) {
                                worst = added;
                            }
                        }
                        const bool pruned = clearing ? (worst && worst.value() < -tolerance)
                                                     : (worst && best && worst.value() < best->worst->margin);
                        if (!pruned) {
                            taken[source] = true;
                            extend(position + 1, worst);
                            taken[source] = false;
                        }
                    }
                };
                extend(0, std::nullopt);
                if (!best) {
                    clearing = false;
                    extend(0, std::nullopt);
                }
                result.order = best->order;
                result.lanes = best->lanes;
                result.tangential = best->tangential;
                if (!clearing) {
                    const auto& closest = best->worst.value();
                    std::ostringstream out;
                    out << std::fixed << std::setprecision(4);
                    out << "No order of the tied " << (entrance ? "entrance" : "exit") << " leads of winding '" << winding
                        << "' keeps every member's turn clear of its siblings' leads. Closest order:";
                    for (size_t index : best->order) {
                        out << " parallel " << routes[index].parallel;
                    }
                    out << "; there " << label_of(routes[closest.member]) << "'s turn passes the lane of "
                        << label_of(routes[closest.sibling]) << " at height " << closest.turnHeight * 1e3 << " mm, "
                        << closest.found * 1e3 << " mm from that lead against the " << closest.required * 1e3
                        << " mm their sections need: " << (closest.required - closest.found) * 1e3 << " mm short.";
                    // Spec section 17: the closest order is drawn and the tie reported (TIED_LEADS_COLLIDE).
                    result.refusal = out.str();
                }
            }
            return result;
        };
        std::map<int, int64_t> nextRank;   // per half of this face
        auto commit = [&](const BundlePlacement& placement) {
            if (placement.refusal) {
                for (size_t index : placement.order) {
                    LeadFinding finding;
                    finding.kind = LeadFindingKind::TIED_LEADS_COLLIDE;
                    finding.witness = placement.refusal.value();
                    plan.findings[index].push_back(finding);
                }
            }
            const int half = realWinding ? plan.halves[placement.order.front()] : 1;
            const int64_t rank = nextRank[half]++;
            for (size_t position = 0; position < placement.order.size(); ++position) {
                const size_t index = placement.order[position];
                slots[index] = placement.lanes[position];
                if (realWinding) {
                    plan.bundleRanks[index] = rank;
                    plan.tangential[index] = placement.tangential[position];
                    plan.overflowed[index] = half * placement.lanes[position] < -tolerance;
                    // Spec section 22: a lead allowed to overflow that still leaves tangentially found no
                    // lane on either half: drawn on its own half and reported (section 17).
                    if (overflowFace && !winding_pinned(routes[index].winding) && placement.tangential[position]) {
                        std::ostringstream out;
                        out << std::fixed << std::setprecision(4) << "Real winding: the terminal lead of " << label_of(routes[index])
                            << " fits neither its group's half of face " << side << " nor the unused lanes of the other half (its "
                            << "turn's half extent is " << turnPaths[index]->laneBound * 1e3 << " mm); it is drawn at lane x = "
                            << placement.lanes[position] * 1e3 << " mm, off the face.";
                        LeadFinding finding;
                        finding.kind = LeadFindingKind::LANES_DO_NOT_FIT;
                        finding.witness = out.str();
                        plan.findings[index].push_back(finding);
                    }
                }
                occupants.push_back({index, footprint_at(index, placement.lanes[position]), extent_of(index, placement.tangential[position]),
                                     placement.tangential[position]});
            }
        };
        auto place_windings = [&](const std::vector<std::string>& windingSequence) {
        for (const auto& winding : windingSequence) {
            auto members_of = [&](ConnectionKind kind) {
                std::vector<size_t> members;
                for (size_t index = 0; index < routes.size(); ++index) {
                    if (routes[index].side == side && routes[index].winding == winding && routes[index].kind == kind) {
                        members.push_back(index);
                    }
                }
                return members;
            };
            const auto entrances = members_of(ConnectionKind::TERMINAL_ENTRANCE);
            const auto exits = members_of(ConnectionKind::TERMINAL_EXIT);
            // Spec S0 (Alf, 2026-10-25, "MKF picks the shorter lead"): when a winding's entrance and exit
            // bundles share a half, they are placed in the order whose leads are SHORTER in total
            // (terminal_lead_length: lead copper plus the turn copper it adds or saves between the
            // crossing and where it leaves), among the orders MKF can plan and, before that, those whose
            // ties clear. An exact tie keeps the entrance first. When neither order can be planned, the
            // entrance-first one is placed and its refusal thrown.
            const bool decide = realWinding && !entrances.empty() && !exits.empty() &&
                                plan.halves[entrances.front()] == plan.halves[exits.front()];
            if (!decide) {
                if (!entrances.empty()) {
                    commit(place_bundle(entrances, occupants));
                }
                if (!exits.empty()) {
                    commit(place_bundle(exits, occupants));
                }
                continue;
            }
            struct Trial {
                std::optional<BundlePlacement> first;
                std::optional<BundlePlacement> second;
                double length = 0;
                // Spec section 17: a tangential lead of this order would cross the far side of its turn. It is
                // still a plan (drawn and reported), but ranks like a refused one: a clean order wins.
                bool crossesFarSide = false;
            };
            auto trial = [&](const std::vector<size_t>& firstMembers, const std::vector<size_t>& secondMembers) -> std::optional<Trial> {
                try {
                    Trial result;
                    auto placed = occupants;
                    result.first = place_bundle(firstMembers, placed);
                    for (size_t position = 0; position < result.first->order.size(); ++position) {
                        const size_t index = result.first->order[position];
                        placed.push_back({index, footprint_at(index, result.first->lanes[position]),
                                          extent_of(index, result.first->tangential[position]), result.first->tangential[position]});
                    }
                    result.second = place_bundle(secondMembers, placed);
                    for (const auto* placement : {&result.first.value(), &result.second.value()}) {
                        for (size_t position = 0; position < placement->order.size(); ++position) {
                            const size_t index = placement->order[position];
                            result.length += terminal_lead_length(routes[index], attachAdvance[index], turnPaths[index].value(),
                                                                  placement->lanes[position]);
                            if (placement->tangential[position] &&
                                plan_terminal_tangential_lead(routes[index], turnPaths[index].value(), placement->lanes[position]).crossesFarSide) {
                                result.crossesFarSide = true;
                            }
                        }
                    }
                    return result;
                }
                catch (const OpenMagneticsException&) {
                    return std::nullopt;   // this order cannot be planned; the other may
                }
            };
            const auto entranceFirst = trial(entrances, exits);
            const auto exitFirst = trial(exits, entrances);
            if (!entranceFirst && !exitFirst) {
                commit(place_bundle(entrances, occupants));   // throws what made both trials fail
                commit(place_bundle(exits, occupants));
                continue;
            }
            auto refused = [](const Trial& t) {
                return t.first->refusal.has_value() || t.second->refusal.has_value() || t.crossesFarSide;
            };
            bool takeExitFirst = false;
            if (!entranceFirst || !exitFirst) {
                takeExitFirst = !entranceFirst;
            }
            else if (refused(entranceFirst.value()) != refused(exitFirst.value())) {
                takeExitFirst = refused(entranceFirst.value());
            }
            else {
                takeExitFirst = exitFirst->length < entranceFirst->length - 1e-12;
            }
            const auto& chosen = takeExitFirst ? exitFirst.value() : entranceFirst.value();
            commit(chosen.first.value());
            commit(chosen.second.value());
        }
        };
        // Spec section 22: an overflowing group stands past the last lane of the group on the spare half,
        // so that group must be placed first. Place in the given order; when a winding that overflowed was
        // placed before one that did not, place again with every overflowing winding after the others.
        const auto occupantsBefore = occupants;
        const auto slotsBefore = slots;
        const auto ranksBefore = plan.bundleRanks;
        const auto findingsBefore = plan.findings;
        const auto overflowedBefore = plan.overflowed;
        place_windings(windingOrder);
        if (overflowFace) {
            std::vector<std::string> settled;
            std::vector<std::string> overflowing;
            for (const auto& winding : windingOrder) {
                bool overflowed = false;
                for (size_t index = 0; index < routes.size(); ++index) {
                    overflowed = overflowed || (routes[index].side == side && routes[index].winding == winding && plan.overflowed[index]);
                }
                (overflowed ? overflowing : settled).push_back(winding);
            }
            auto sequence = settled;
            sequence.insert(sequence.end(), overflowing.begin(), overflowing.end());
            if (sequence != windingOrder) {
                occupants = occupantsBefore;
                slots = slotsBefore;
                plan.bundleRanks = ranksBefore;
                plan.findings = findingsBefore;
                plan.overflowed = overflowedBefore;
                nextRank.clear();
                place_windings(sequence);
            }
        }
    }
    const auto laneFindings = check_terminal_lanes(routes, diameters, slots, realWinding, turnPaths);
    for (size_t index = 0; index < routes.size(); ++index) {
        plan.findings[index].insert(plan.findings[index].end(), laneFindings[index].begin(), laneFindings[index].end());
    }
    return plan;
}

std::vector<std::vector<LeadFinding>> Coil::check_terminal_lanes(const std::vector<ConnectionRoute>& routes, const std::vector<double>& diameters,
                                const std::vector<std::optional<double>>& slots, bool realWinding,
                                const std::vector<std::optional<TerminalTurnPath>>& turnPaths) {
    if (diameters.size() != routes.size() || slots.size() != routes.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "check_terminal_lanes needs one diameter and one slot per route: " + std::to_string(routes.size()) + " routes, " +
            std::to_string(diameters.size()) + " diameters, " + std::to_string(slots.size()) + " slots.");
    }
    if (realWinding && turnPaths.size() != routes.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "check_terminal_lanes under real winding needs one turn path per route, for each lead's real sections: " +
            std::to_string(routes.size()) + " routes, " + std::to_string(turnPaths.size()) + " turn paths.");
    }
    // A clearance is met when the found distance is at least the required one less this: the lanes
    // are sums of lengths in metres, so a pair placed at exactly the required distance must pass.
    const double tolerance = 1e-9;
    auto millimetres = [](double metres) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(4) << metres * 1e3 << " mm";
        return out.str();
    };
    // Spec section 17: under real winding a collision is drawn and reported on both leads; under ideal
    // winding the planner's own lanes colliding is a defect, thrown as before.
    std::vector<std::vector<LeadFinding>> findings(routes.size());
    size_t i = 0, j = 0;
    auto refuse = [&](const std::string& refusal) {
        if (!realWinding) {
            throw std::logic_error(refusal);
        }
        LeadFinding finding;
        finding.kind = LeadFindingKind::LANES_COLLIDE;
        finding.witness = refusal;
        findings[i].push_back(finding);
        findings[j].push_back(finding);
    };
    for (i = 0; i < routes.size(); ++i) {
        if (!is_terminal_route(routes[i]) || !slots[i]) {
            continue;
        }
        if (realWinding && !turnPaths[i]) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "check_terminal_lanes: the terminal lead of " + terminal_route_label(routes[i]) + " has a lane but no turn path.");
        }
        for (j = i + 1; j < routes.size(); ++j) {
            if (!is_terminal_route(routes[j]) || !slots[j] || routes[j].side != routes[i].side) {
                continue;
            }
            if (realWinding && !turnPaths[j]) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "check_terminal_lanes: the terminal lead of " + terminal_route_label(routes[j]) + " has a lane but no turn path.");
            }
            const double gap = footprint_gap(lane_footprint(routes[i], slots[i].value()), lane_footprint(routes[j], slots[j].value()));
            auto witness = [&]() {
                return terminal_route_label(routes[i]) + " at x " + millimetres(slots[i].value()) + " and " +
                       terminal_route_label(routes[j]) + " at x " + millimetres(slots[j].value()) + ": ";
            };
            if (!realWinding) {
                const double radii = (diameters[i] + diameters[j]) / 2;
                const double section = polyline_distance(route_polyline(routes[i]), route_polyline(routes[j]));
                if (gap < radii - tolerance && section < radii - tolerance) {
                    refuse("Terminal lanes collide: " + witness() + "their footprints along the connection face are " +
                           millimetres(gap) + " apart and their routes in the winding window " + millimetres(section) +
                           ", both under the " + millimetres(radii) + " their coated radii need.");
                }
                continue;
            }
            // Real winding: two leads that can meet -- their drawn paths at their stations, bends and
            // each piece's real section included (no lane-dependent turn heights are known here) --
            // stand their two half extents along the lane apart along the face (spec C3: no bend-radius
            // gap; spec RW-size: each lead's real extent). Siblings and different bundles alike.
            const auto& pathI = turnPaths[i].value();
            const auto& pathJ = turnPaths[j].value();
            const bool tangentialI = std::abs(slots[i].value()) >= pathI.laneBound - 1e-12;
            const bool tangentialJ = std::abs(slots[j].value()) >= pathJ.laneBound - 1e-12;
            const auto& attachI = routes[i].kind == ConnectionKind::TERMINAL_ENTRANCE ? routes[i].waypoints.back() : routes[i].waypoints.front();
            const auto& attachJ = routes[j].kind == ConnectionKind::TERMINAL_ENTRANCE ? routes[j].waypoints.back() : routes[j].waypoints.front();
            const double required = (lead_lane_extent(routes[i], pathI, tangentialI) + lead_lane_extent(routes[j], pathJ, tangentialJ)) / 2;
            if (!(gap < required - tolerance)) {
                continue;
            }
            const bool meet = envelopes_meet(thick_lead_envelope(routes[i], pathI, tangentialI, {attachI.at(1), attachI.at(1)}),
                                             thick_lead_envelope(routes[j], pathJ, tangentialJ, {attachJ.at(1), attachJ.at(1)}));
            // As terminal_lane_plan: another winding's lead always counts; only a winding's own leads that
            // cannot meet may share the face.
            if (!meet && routes[i].winding == routes[j].winding) {
                continue;
            }
            const double section = polyline_distance(route_polyline(routes[i]), route_polyline(routes[j]));
            if (same_bundle(routes[i], routes[j])) {
                refuse("Terminal lanes collide: " + witness() + "the two siblings can meet in 3D (their drawn paths, bends and real " +
                       "sections included, touch in the winding window), yet their footprints along the connection face stand only " +
                       millimetres(gap) + " apart, under the " + millimetres(required) + " their half extents along the lane need.");
            }
            else if (!meet) {
                refuse("Terminal lanes collide: " + witness() + "the two leads belong to different windings, which never share a lane " +
                       "(their routes come " + millimetres(section) + " apart in the winding window), yet their footprints along the " +
                       "connection face stand only " + millimetres(gap) + " apart, under the " + millimetres(required) +
                       " their half extents along the lane need.");
            }
            else {
                refuse("Terminal lanes collide: " + witness() + "the two leads belong to different bundles and can meet in 3D (their " +
                       "drawn paths, bends and real sections included, touch in the winding window; their routes come " +
                       millimetres(section) + " apart), yet their footprints along the connection face stand only " + millimetres(gap) +
                       " apart, under the " + millimetres(required) + " their half extents along the lane need.");
            }
        }
    }
    return findings;
}

std::vector<PinLeadRoute> Coil::route_leads_to_pins(const std::vector<MAS::Pin>& bobbinPins,
                                                    const std::vector<PinLeadRequest>& leads,
                                                    double frontFaceOffset, int64_t wrapTurns,
                                                    const std::vector<Bobbin::PinRailBlock>& pinRails) {
    std::vector<PinLeadRoute> result(leads.size());
    if (leads.empty()) {
        return result;
    }
    if (wrapTurns < 1) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Routing leads to pins needs at least one wrap turn per lead.");
    }
    const auto placedPins = place_pins(bobbinPins);
    std::map<std::string, size_t> placedIndexByName;
    for (size_t index = 0; index < placedPins.size(); ++index) {
        placedIndexByName[placedPins[index].name] = index;
    }
    // Each row's centre along x: a lead approaches its pin from the side away from it first.
    std::map<size_t, std::pair<double, size_t>> rowSum;
    for (const auto& placedPin : placedPins) {
        rowSum[placedPin.row].first += placedPin.centre[0];
        rowSum[placedPin.row].second += 1;
    }
    auto row_centre = [&](size_t row) { return rowSum.at(row).first / double(rowSum.at(row).second); };

    const double tolerance = 1e-9;
    double pitch = 0;
    std::vector<size_t> pinOfLead(leads.size());
    for (size_t index = 0; index < leads.size(); ++index) {
        const auto& lead = leads[index];
        if (lead.windowExit.size() < 2) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Lead " + lead.label + " has no window exit {radial, axial}.");
        }
        if (!(lead.diameter > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
                "Lead " + lead.label + " has no outer diameter, so its clearance to other leads and pins is unknown.");
        }
        if (lead.lift < 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Lead " + lead.label + " has a negative ride-over lift.");
        }
        if (!lead.exitX || !std::isfinite(lead.exitX.value())) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Lead " + lead.label + " has no exit slot x (Coil::terminal_exit_slots): where its in-window run leaves the window is unknown.");
        }
        if (lead.side != 0 && lead.side != 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Lead " + lead.label + " leaves on connection face " + std::to_string(lead.side) + "; the faces are 0 (-Z) and 1 (+Z).");
        }
        if (!lead.pin.get_name()) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Lead " + lead.label + " is assigned to a pin with no name.");
        }
        auto found = placedIndexByName.find(lead.pin.get_name().value());
        if (found == placedIndexByName.end()) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "Lead " + lead.label + " is assigned to pin '" + lead.pin.get_name().value() + "', which is not among the bobbin's pins.");
        }
        pinOfLead[index] = found->second;
        pitch = std::max(pitch, lead.diameter);
    }

    // Planning order. A lead planned later takes a lane further out, so its leg out of the window
    // crosses the lanes of every lead planned before it: that leg must not lie between an earlier
    // lead's exit and its pin. So leads are planned by how far their exit slot sits TOWARD their pin
    // (exit x for a lead heading -x, minus exit x for one heading +x), and a later lead leaves the
    // window behind every earlier one (ABT #1237: with the exit slots fixed, the boost's parallel 1
    // start at x = 0.943 mm planned first sent its lane over parallel 0's exit at x = 0). Ties: inner
    // pins first, outer leads then pass outside the inner ones (Wuerth DFM 2024-04-29 08:42).
    std::vector<size_t> order(leads.size());
    for (size_t index = 0; index < leads.size(); ++index) {
        order[index] = index;
    }
    // Spec section 12: a lead leaves on its own face (PinLeadRequest::side). Its exit slot is in that
    // face's frame; face 1 is the bobbin frame turned 180 degrees about the column, so there the run
    // starts at bobbin x = -exitX and steps out along +Z. Everything below is in the bobbin frame.
    auto bobbin_exit_x = [&](size_t index) { return leads[index].side == 0 ? leads[index].exitX.value() : -leads[index].exitX.value(); };
    auto outward_z = [&](size_t index) { return leads[index].side == 0 ? -1.0 : 1.0; };
    auto toward_pin = [&](size_t index) {
        const double exitX = bobbin_exit_x(index);
        return placedPins[pinOfLead[index]].centre[0] < exitX ? exitX : -exitX;
    };
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (std::abs(toward_pin(a) - toward_pin(b)) > tolerance) {
            return toward_pin(a) < toward_pin(b);
        }
        const auto& pinA = placedPins[pinOfLead[a]];
        const auto& pinB = placedPins[pinOfLead[b]];
        if (pinA.row != pinB.row) {
            return pinA.row < pinB.row;
        }
        const double innerA = std::abs(pinA.centre[0] - row_centre(pinA.row));
        const double innerB = std::abs(pinB.centre[0] - row_centre(pinB.row));
        if (std::abs(innerA - innerB) > tolerance) {
            return innerA < innerB;
        }
        return false;
    });

    const size_t numberLeads = leads.size();
    std::vector<RoutedPinLead> routed;

    // Every lead's in-window run where it leaves the window: at its exit slot x and exit row, from the
    // window border out to its exit depth (the ride-over lift). The in-window lead is drawn there
    // (MKF's slot), so no other lead's run may pass through it, whichever lead is planned first.
    std::vector<std::vector<std::vector<double>>> exitRuns(numberLeads);
    for (size_t index = 0; index < numberLeads; ++index) {
        const auto& lead = leads[index];
        const double x = bobbin_exit_x(index);
        const double zSign = outward_z(index);
        exitRuns[index] = {{x, lead.windowExit[1], zSign * (lead.windowExit[0] + frontFaceOffset)},
                           {x, lead.windowExit[1], zSign * (lead.windowExit[0] + frontFaceOffset + lead.lift)}};
    }
    for (size_t i = 0; i < numberLeads; ++i) {
        for (size_t j = i + 1; j < numberLeads; ++j) {
            const double found = polyline_distance(exitRuns[i], exitRuns[j]);
            const double required = (leads[i].diameter + leads[j].diameter) / 2;
            if (found < required - 1e-9) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The exit slots of " + leads[i].label + " and " + leads[j].label + " share copper where they leave the window: centreline distance " +
                    format_millimetres(found) + " < " + format_millimetres(required) + ". Fix the exit slots (Coil::terminal_exit_slots).");
            }
        }
    }
    // ABT #1336 (2026-10-25): two leads whose exit slots stand closer along the flange than their two
    // radii SHARE A LANE. terminal_exit_slots lets them only because their leads cannot meet
    // (terminal_leads_can_meet: an exact touch counts as meeting), and the pin run is the same lead
    // carried on past the window border. So the runs of a lane-sharing pair keep that premise: they
    // may not even touch, each other or the other's in-window run. Leads on lanes of their own keep
    // the copper rule (an exact touch is two neighbouring lanes, as in the window).
    std::vector<std::vector<bool>> sharesLane(numberLeads, std::vector<bool>(numberLeads, false));
    for (size_t i = 0; i < numberLeads; ++i) {
        for (size_t j = 0; j < numberLeads; ++j) {
            sharesLane[i][j] = i != j && leads[i].side == leads[j].side &&
                               std::abs(leads[i].exitX.value() - leads[j].exitX.value()) < (leads[i].diameter + leads[j].diameter) / 2 - tolerance;
        }
    }
    auto too_close = [&](size_t a, size_t b, double found, double required) {
        return sharesLane[a][b] ? found <= required + tolerance : found < required - tolerance;
    };
    auto shortfall_text = [&](size_t a, size_t b) {
        return std::string(sharesLane[a][b] ? " (they share a lane, so their leads may not even touch)" : "");
    };

    auto axis_point = [](const PlacedPin& pin, size_t axisIndex, double coordinate) {
        auto point = pin.centre;
        point[axisIndex] = coordinate;
        return point;
    };
    auto millimetres = [](double metres) { return format_millimetres(metres); };

    for (size_t leadIndex : order) {
        const auto& lead = leads[leadIndex];
        const auto& pin = placedPins[pinOfLead[leadIndex]];
        const size_t axisIndex = pin.hangsAlongZ ? 2 : 1;   // the coordinate the pin leaves along (decreasing)
        const size_t rowIndex = pin.hangsAlongZ ? 1 : 2;    // the coordinate that locates its row
        const double radius = lead.diameter / 2;
        // ABT #1172: every corner of a pin run turns 90 degrees around an obstacle edge (the rail,
        // the flange, the core), so one clearance serves them all. Equals `radius` when no bend is
        // declared, which is the historical geometry.
        if (!(lead.bendRadius >= radius)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The lead " + lead.label + " carries a planned bend radius of " + format_millimetres(lead.bendRadius) +
                ", below the radius it sweeps (" + format_millimetres(radius) + "): set it from Coil::lead_bend_radius");
        }
        const double legClearance = Settings::lead_leg_clearance(radius, lead.bendRadius, std::numbers::pi / 2);
        const double plannedBendRadius = lead.bendRadius;
        const double pinRadius = pin.diameter / 2;
        const double wrapRadius = pinRadius + radius;
        const double approach = wrapRadius + lead.diameter;   // clear of a neighbour wrap of this wire
        const double tip = 2 * pin.centre[axisIndex] - pin.base[axisIndex];
        const double zSign = outward_z(leadIndex);   // the direction away from the window, on this lead's face
        const double exitDepth = zSign * (lead.windowExit[0] + frontFaceOffset + lead.lift);
        auto beyond_exit = [&](double depth) { return zSign * depth > zSign * exitDepth + tolerance; };
        const double exitAxial = lead.windowExit[1];
        const double outward = pin.centre[0] - row_centre(pin.row);
        const double preferredSide = outward > tolerance ? 1.0 : -1.0;

        // Lanes: stepping out from the exit, and stepping out from each face of every pin rail
        // beyond the exit (a drop has to pass outside a rail it cannot pass through).
        std::vector<double> laneDepths;
        for (size_t lane = 0; lane <= numberLeads; ++lane) {
            laneDepths.push_back(exitDepth + zSign * double(lane) * pitch);
        }
        for (const auto& rail : pinRails) {
            if (rail.centre.size() < 3 || rail.halfExtents.size() < 3) {
                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Pin rail '" + rail.name + "' has no 3D centre and half extents.");
            }
            // ABT #1172: step off the face by the BEND clearance, not the bare wire radius. The
            // drop turns 90 degrees around the rail's edge, and a corner drawn with a centreline
            // radius R cuts the inside of that turn; d >= R - (R - r)/sqrt(2) keeps the copper out
            // of the rail. With no bend declared this is exactly r and the geometry is unchanged.
            for (double face : {rail.centre[2] - rail.halfExtents[2] - legClearance, rail.centre[2] + rail.halfExtents[2] + legClearance}) {
                for (size_t lane = 0; lane <= numberLeads; ++lane) {
                    if (beyond_exit(face + zSign * double(lane) * pitch)) {
                        laneDepths.push_back(face + zSign * double(lane) * pitch);
                    }
                }
            }
        }
        // ... and stepping out from beyond every other lead's in-window run (a lead lifted further out
        // than this one leaves the window across this lead's lanes).
        for (size_t other = 0; other < numberLeads; ++other) {
            if (other == leadIndex) {
                continue;
            }
            const double beyond = exitRuns[other].back()[2] + zSign * (radius + leads[other].diameter / 2);
            for (size_t lane = 0; lane <= numberLeads; ++lane) {
                if (beyond_exit(beyond + zSign * double(lane) * pitch)) {
                    laneDepths.push_back(beyond + zSign * double(lane) * pitch);
                }
            }
        }
        // Nearest the exit first.
        std::sort(laneDepths.begin(), laneDepths.end(), [&](double a, double b) { return zSign * a < zSign * b; });
        laneDepths.erase(std::unique(laneDepths.begin(), laneDepths.end(), [](double a, double b) { return std::abs(a - b) < 1e-12; }),
                         laneDepths.end());
        std::optional<std::string> directConflict;
        std::optional<RoutedPinLead> directCandidate;   // the most direct run, drawn when none clears (spec section 17)
        std::optional<RoutedPinLead> chosen;
        size_t candidates = 0;
        const size_t maximumLevel = (size_t(wrapTurns) + 1) * numberLeads;
        for (size_t level = 0; level <= maximumLevel && !chosen; ++level) {
            // ABT #1172: same clearance on the leg that runs under the rail — it is the other
            // half of the same corner.
            const double wrapTop = pin.base[axisIndex] - legClearance - double(level) * pitch;
            const double wrapBottom = wrapTop - double(wrapTurns) * lead.diameter;
            if (wrapBottom < tip - tolerance) {
                break;   // the wrap would run off the pin's tip
            }
            for (double side : {preferredSide, -preferredSide}) {
                if (chosen) {
                    break;
                }
                for (size_t besideRank = 0; besideRank <= numberLeads && !chosen; ++besideRank) {
                const double besidePin = pin.centre[0] + side * (approach + double(besideRank) * pitch);
                for (size_t lane = 0; lane < laneDepths.size() && !chosen; ++lane) {
                    const double laneDepth = laneDepths[lane];
                    {
                        const double exitX = bobbin_exit_x(leadIndex);
                        std::vector<std::vector<double>> raw;
                        raw.push_back({exitX, exitAxial, exitDepth});
                        raw.push_back({exitX, exitAxial, laneDepth});
                        raw.push_back({besidePin, exitAxial, laneDepth});
                        auto railPoint = raw.back();
                        railPoint[axisIndex] = wrapTop;
                        raw.push_back(railPoint);
                        // ABT #1640: the last leg arrives TANGENT to the wrap circle, a wrap radius
                        // to the side of the pin's row it comes from, and ends where the wrap
                        // begins. A leg into the pin's axis met the wrap at 90 degrees, a corner no
                        // bend can round without cutting into the pin (an arc from a radial leg onto
                        // the circle curves inside it); a tangent arrival has none.
                        const double rowSide = railPoint[rowIndex] >= pin.centre[rowIndex] ? 1.0 : -1.0;
                        auto rowPoint = railPoint;
                        rowPoint[rowIndex] = pin.centre[rowIndex] + rowSide * wrapRadius;
                        raw.push_back(rowPoint);
                        auto wrapStart = axis_point(pin, axisIndex, wrapTop);
                        wrapStart[rowIndex] += rowSide * wrapRadius;
                        raw.push_back(wrapStart);
                        std::vector<std::vector<double>> points;
                        for (auto& point : raw) {
                            if (points.empty() || distance(points.back(), point) > 1e-12) {
                                points.push_back(point);
                            }
                        }
                        ++candidates;
                        std::optional<std::string> conflict;
                        // A leg that doubles back over the previous one is not a route.
                        for (size_t k = 0; k + 2 < points.size() && !conflict; ++k) {
                            double dotProduct = 0;
                            for (int c = 0; c < 3; ++c) {
                                dotProduct += (points[k + 1][c] - points[k][c]) * (points[k + 2][c] - points[k + 1][c]);
                            }
                            const double lengths = distance(points[k], points[k + 1]) * distance(points[k + 1], points[k + 2]);
                            if (dotProduct < -0.99 * lengths) {
                                conflict = "its own leg " + std::to_string(k) + " (it folds back on itself)";
                            }
                        }
                        RoutedPinLead candidate;
                        candidate.request = leadIndex;
                        candidate.pinIndex = pinOfLead[leadIndex];
                        candidate.radius = radius;
                        candidate.wrapRadius = wrapRadius;
                        candidate.points = points;
                        candidate.wrapAxis = {axis_point(pin, axisIndex, wrapTop), axis_point(pin, axisIndex, wrapBottom)};
                        candidate.wrapTop = wrapTop;
                        candidate.wrapBottom = wrapBottom;
                        candidate.dropX = besidePin;
                        // Pin rails: the wire goes round and under them, never through.
                        for (const auto& rail : pinRails) {
                            if (conflict) {
                                break;
                            }
                            for (size_t k = 0; k + 1 < points.size(); ++k) {
                                const double found = segment_box_distance(points[k], points[k + 1], rail);
                                if (found < radius - tolerance) {
                                    conflict = "pin rail '" + rail.name + "': leg " + std::to_string(k) + " comes " + millimetres(found) +
                                               " from it, a wire radius " + millimetres(radius) + " needed";
                                    break;
                                }
                            }
                        }
                        // Pins: every pin but its own keeps the approach clearance from the whole run;
                        // its own pin, from everything but the last leg, which touches it: that leg
                        // runs tangent to the wrap circle, a pin radius plus a wire radius from the
                        // axis (ABT #1640), and may come no closer.
                        for (size_t pinIndex = 0; pinIndex < placedPins.size() && !conflict; ++pinIndex) {
                            const auto& other = placedPins[pinIndex];
                            const size_t otherAxis = other.hangsAlongZ ? 2 : 1;
                            std::vector<std::vector<double>> pinAxis = {other.base, axis_point(other, otherAxis, 2 * other.centre[otherAxis] - other.base[otherAxis])};
                            const bool own = pinIndex == candidate.pinIndex;
                            const double required = other.diameter / 2 + radius + lead.diameter;
                            const double found = polyline_distance(points, pinAxis, own);
                            if (found < required - tolerance) {
                                conflict = "pin '" + other.name + "': centreline " + millimetres(found) + " from its axis, " +
                                           millimetres(required) + " needed (pin radius + wire radius + wire diameter)";
                            }
                            if (own && !conflict) {
                                const auto& last = points.size() >= 2 ? points[points.size() - 2] : points.back();
                                const double touching = segment_distance(last, points.back(), pinAxis[0], pinAxis[1]);
                                if (touching < wrapRadius - tolerance) {
                                    conflict = "its own pin '" + other.name + "': the arriving leg comes " + millimetres(touching) +
                                               " from its axis, inside the wrap radius " + millimetres(wrapRadius);
                                }
                            }
                        }
                        // Other leads' in-window runs at their exit slots, planned or not.
                        for (size_t other = 0; other < numberLeads && !conflict; ++other) {
                            if (other == leadIndex) {
                                continue;
                            }
                            const double required = radius + leads[other].diameter / 2;
                            const double found = polyline_distance(points, exitRuns[other]);
                            if (too_close(leadIndex, other, found, required)) {
                                conflict = "the in-window run of " + leads[other].label + " at its exit slot x = " +
                                           millimetres(leads[other].exitX.value()) + ": centreline distance " + millimetres(found) +
                                           (sharesLane[leadIndex][other] ? " <= " : " < ") + millimetres(required) +
                                           shortfall_text(leadIndex, other);
                            }
                        }
                        for (const auto& placed : routed) {
                            if (conflict) {
                                break;
                            }
                            const auto& otherLabel = leads[placed.request].label;
                            const double required = radius + placed.radius;
                            // Every lead drops to the rail in its own plane along the flange: two drops
                            // at one x would only be kept apart by their lanes (ABT #1237).
                            if (std::abs(besidePin - placed.dropX) < required - tolerance) {
                                conflict = "the drop of " + otherLabel + " at x = " + millimetres(placed.dropX) +
                                           ": " + millimetres(std::abs(besidePin - placed.dropX)) + " apart along the flange, " +
                                           millimetres(required) + " needed";
                                break;
                            }
                            double found = polyline_distance(points, placed.points);
                            if (too_close(leadIndex, placed.request, found, required)) {
                                conflict = "the run of " + otherLabel + ": centreline distance " + millimetres(found) +
                                           (sharesLane[leadIndex][placed.request] ? " <= " : " < ") + millimetres(required) +
                                           shortfall_text(leadIndex, placed.request);
                                break;
                            }
                            const bool samePin = placed.pinIndex == candidate.pinIndex;
                            // This run against that wrap, and that run against this wrap.
                            found = polyline_distance(points, placed.wrapAxis, samePin) - placed.wrapRadius;
                            if (found < required - tolerance) {
                                conflict = "the wrap of " + otherLabel + " on pin '" + placedPins[placed.pinIndex].name +
                                           "': centreline distance " + millimetres(found) + " < " + millimetres(required);
                                break;
                            }
                            found = polyline_distance(placed.points, candidate.wrapAxis, samePin) - wrapRadius;
                            if (found < required - tolerance) {
                                conflict = "the run of " + otherLabel + " against this lead's wrap on pin '" + pin.name +
                                           "': centreline distance " + millimetres(found) + " < " + millimetres(required);
                                break;
                            }
                            if (samePin) {
                                // Two wraps on one pin (strands sharing it) are separated along the axis;
                                // each arrival leg lies at its wrap's top, so that gap covers them too.
                                const double gap = std::max(placed.wrapBottom - wrapTop, wrapBottom - placed.wrapTop);
                                if (gap < required - tolerance) {
                                    conflict = "the wrap of " + otherLabel + " on the same pin '" + pin.name + "': " +
                                               millimetres(std::max(gap, 0.0)) + " apart along the pin, " + millimetres(required) + " needed";
                                    break;
                                }
                            }
                            else {
                                found = segment_distance(candidate.wrapAxis[0], candidate.wrapAxis[1], placed.wrapAxis[0], placed.wrapAxis[1]) -
                                        wrapRadius - placed.wrapRadius;
                                if (found < required - tolerance) {
                                    conflict = "the wrap of " + otherLabel + " on pin '" + placedPins[placed.pinIndex].name +
                                               "' against this lead's wrap: centreline distance " + millimetres(found) + " < " + millimetres(required);
                                    break;
                                }
                            }
                        }
                        if (!conflict) {
                            chosen = candidate;
                        }
                        else if (!directConflict) {
                            directConflict = conflict;
                            directCandidate = candidate;
                        }
                    }
                }
                }
            }
        }
        if (!chosen) {
            const std::string refusal =
                "Cannot route " + lead.label + " to pin '" + pin.name + "' clear of the other leads and pins: all " +
                std::to_string(candidates) + " candidate runs from its exit slot x = " + millimetres(lead.exitX.value()) +
                " (wrap levels on the pin, approach sides, drop offsets, lanes) share "
                "copper, and the most direct one collides with " + directConflict.value_or("nothing recorded") +
                ". Assign that end to another pin, use a thinner wire, or a former with a wider pin pitch or longer pins.";
            // Spec section 17 (Alf, 2026-10-26): the most direct run is drawn, collision and all, and the
            // collision reported on the run (PIN_RUN_COLLIDES). Only when no candidate run could be built
            // at all is there nothing to draw: that still throws.
            if (!directCandidate) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, refusal);
            }
            chosen = directCandidate;
            chosen->collision = refusal;
        }
        routed.push_back(chosen.value());
    }

    for (const auto& lead : routed) {
        PinLeadRoute route;
        route.pinName = placedPins[lead.pinIndex].name;
        route.waypoints = lead.points;
        // ABT #1172/#1296: the bend this run's corners were planned for, the lead's own.
        route.plannedBendRadius = leads[lead.request].bendRadius;
        if (lead.collision) {
            LeadFinding finding;
            finding.kind = LeadFindingKind::PIN_RUN_COLLIDES;
            finding.witness = lead.collision.value();
            route.findings.push_back(finding);
        }
        for (size_t index = 0; index + 1 < route.waypoints.size(); ++index) {
            route.length += distance(route.waypoints[index], route.waypoints[index + 1]);
        }
        result[lead.request] = route;
    }
    return result;
}

double Coil::lead_front_face_offset(Bobbin bobbin) {
    if (!bobbin.get_processed_description()) {
        throw CoilNotProcessedException("the lead's front-face offset needs a processed bobbin");
    }
    auto processed = bobbin.get_processed_description().value();
    auto shape = processed.get_column_shape();
    if (shape == ColumnShape::ROUND) {
        return 0.0;
    }
    if (shape == ColumnShape::RECTANGULAR || shape == ColumnShape::OBLONG) {
        if (!processed.get_column_width()) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "The bobbin's column is rectangular or oblong but states no columnWidth.");
        }
        // columnDepth and columnWidth are HALF dimensions. A rectangular column's front face is
        // their difference deeper than its side face; an oblong one with no straight is round.
        const double offset = processed.get_column_depth() - processed.get_column_width().value();
        if (shape == ColumnShape::OBLONG && offset <= 0) {
            return 0.0;
        }
        return offset;
    }
    throw NotImplementedException("Routing a lead to a pin around an irregular column");
}

bool Coil::has_complete_pin_connections() const {
    for (const auto& winding : get_functional_description()) {
        if (!winding.get_connections()) {
            return false;
        }
        auto connections = winding.get_connections().value();
        if (connections.size() < 2) {
            return false;
        }
        bool hasEnds = false;
        bool hasStart = false;
        bool hasFinish = false;
        for (const auto& connection : connections) {
            if (connection.get_end()) {
                hasEnds = true;
                if (connection.get_end().value() == End::START) {
                    hasStart = true;
                }
                if (connection.get_end().value() == End::FINISH) {
                    hasFinish = true;
                }
            }
            if (!connection.get_pin_name() && is_pin_connection_type(connection.get_type())) {
                return false;
            }
        }
        if (hasEnds && !(hasStart && hasFinish)) {
            return false;
        }
    }
    return true;
}

namespace {

// ABT #1173 (WP4): the toroid rule of RFC 0013 ("Toroids on a base"). A toroid base's pins stand in a
// rectangle, two rows, one per side of the base. Each of two windings takes the side its half of the ring
// faces (design order only when neither side faces one winding more than the other), and walks that row
// from the end its first turn lies towards. For windings with no such preference every row is walked
// from -X, so on the counter-clockwise numbering of Bobbin::expand_pinout a common-mode choke on a
// 4-pin base reads winding 1 = pins 1 (start) and 2 (finish), winding 2 = pins 4 (start) and 3
// (finish), the WE-CMB footprint "1 2 / 4 3". Two windings wound on opposite halves of the ring thus
// land on opposite sides of the base. Inside a row the pins are allocated as on a former's row: the
// start half towards -X, the finish half towards +X, inner pins first (adjacent start and finish),
// strands per rule 4 (already decided in the plan), never a corner pin for insulated or margin-wound
// leads (rule 7). A single winding whose ends do not fit one side puts its starts (and taps) on side 0
// and its finishes on side 1. More windings than sides, a winding with user pins on both sides, or two
// windings claiming one side throw. The creepage between isolation sides is checked by the caller,
// exactly as on a former.
void assign_pins_on_toroid_base(std::vector<WindingPinPlan>& plans, const std::vector<PlacedPin>& placedPins,
                                const std::vector<std::vector<size_t>>& rowWalk,
                                const std::map<std::string, size_t>& placedIndexByName,
                                std::vector<std::optional<std::pair<size_t, std::string>>>& pinOwner,
                                const std::string& bobbinName) {
    const size_t numberRows = rowWalk.size();
    if (numberRows != 2) {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
            "Toroid base '" + bobbinName + "' has " + std::to_string(numberRows) + " pin rows. The toroid rule gives each "
            "winding one of the base's two sides (RFC 0013, 'Toroids on a base'), so the pins must stand in two rows.");
    }
    auto known_pins_text = [&]() {
        std::string text;
        for (const auto& placedPin : placedPins) {
            text += (text.empty() ? "" : ", ") + placedPin.name;
        }
        return text;
    };

    // --- user pins: exist, one owner each, and the side they put their winding on ------------------
    std::map<size_t, std::set<size_t>> userRowsByPlan;
    for (size_t planIndex = 0; planIndex < plans.size(); ++planIndex) {
        for (auto& slot : plans[planIndex].slots) {
            if (!slot.userPin) {
                continue;
            }
            auto found = placedIndexByName.find(slot.pinName.value());
            if (found == placedIndexByName.end()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Winding '" + plans[planIndex].name + "' terminates on pin '" + slot.pinName.value() +
                    "', which the toroid base '" + bobbinName + "' does not have. It has: " + known_pins_text() + ".");
            }
            const size_t pinIndex = found->second;
            if (pinOwner[pinIndex] && pinOwner[pinIndex]->first != planIndex) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Pin '" + slot.pinName.value() + "' of toroid base '" + bobbinName + "' is given to winding '" +
                    pinOwner[pinIndex]->second + "' and to winding '" + plans[planIndex].name + "'.");
            }
            pinOwner[pinIndex] = std::make_pair(planIndex, plans[planIndex].name);
            userRowsByPlan[planIndex].insert(placedPins[pinIndex].row);
        }
    }

    // --- windings in the order the design lists them ------------------------------------------------
    std::vector<size_t> planned;
    for (size_t planIndex = 0; planIndex < plans.size(); ++planIndex) {
        if (!plans[planIndex].legacy) {
            planned.push_back(planIndex);
        }
    }
    std::stable_sort(planned.begin(), planned.end(),
                     [&](size_t a, size_t b) { return plans[a].windingIndex < plans[b].windingIndex; });
    if (planned.size() > numberRows) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "The coil has " + std::to_string(planned.size()) + " windings to terminate on toroid base '" + bobbinName +
            "', which has two sides; the toroid rule gives each winding a side of its own (RFC 0013, 'Toroids on a "
            "base'). Use a base per winding pair or terminate the extra windings as flying leads.");
    }

    // --- one winding's slots onto pins --------------------------------------------------------------
    auto free_pins = [&](size_t row) {
        std::vector<size_t> freeWalk;
        for (auto pinIndex : rowWalk[row]) {
            if (!pinOwner[pinIndex]) {
                freeWalk.push_back(pinIndex);
            }
        }
        return freeWalk;
    };
    auto slots_of = [&](WindingPinPlan& plan, bool startHalf) {
        std::vector<PinSlot*> slots;
        for (auto& slot : plan.slots) {
            if (!slot.userPin && !slot.offPin && slot.startHalf == startHalf) {
                slots.push_back(&slot);
            }
        }
        return slots;
    };
    auto give = [&](size_t planIndex, PinSlot* slot, size_t pinIndex) {
        slot->pinName = placedPins[pinIndex].name;
        pinOwner[pinIndex] = std::make_pair(planIndex, plans[planIndex].name);
    };
    auto usable = [&](size_t planIndex, const std::vector<size_t>& pins) {
        std::vector<size_t> result;
        for (auto pinIndex : pins) {
            if (!(plans[planIndex].avoidCorners && placedPins[pinIndex].corner)) {
                result.push_back(pinIndex);
            }
        }
        return result;
    };
    auto row_too_small = [&](size_t planIndex, size_t row, size_t demand, size_t available) {
        return InvalidInputException(ErrorCode::INVALID_INPUT,
            "Side " + std::to_string(row) + " of toroid base '" + bobbinName + "' has " + std::to_string(available) +
            " free pins for winding '" + plans[planIndex].name + "', which needs " + std::to_string(demand) +
            (plans[planIndex].avoidCorners ? " (its insulated or margin-wound leads may not take a corner pin, RFC 0013 rule 7)" : std::string("")) +
            ". Use a base with more pins, or fewer separate strands and taps.");
    };
    // Starts and finishes on one side: the start half towards the end of the row the winding's first turn
    // lies towards, the finish half towards the other, inner first. A winding whose first turn is further
    // along +X than its last walks the row from +X, so neither lead crosses the other on its way to the row
    // (ABT #1173; a winding that starts and finishes at the same x keeps the -X walk).
    const double sameEnd = 1e-6;
    auto place_on_one_side = [&](size_t planIndex, size_t row) {
        auto& plan = plans[planIndex];
        auto startSlots = slots_of(plan, true);
        auto finishSlots = slots_of(plan, false);
        auto freeWalk = free_pins(row);
        if (plan.firstTurnX > plan.lastTurnX + sameEnd) {
            std::reverse(freeWalk.begin(), freeWalk.end());
        }
        const size_t demand = startSlots.size() + finishSlots.size();
        if (demand > freeWalk.size()) {
            throw row_too_small(planIndex, row, demand, freeWalk.size());
        }
        const size_t startHalfSize = startSlots.size() + (freeWalk.size() - demand) / 2;
        std::vector<size_t> startHalf(freeWalk.begin(), freeWalk.begin() + startHalfSize);
        std::vector<size_t> finishHalf(freeWalk.begin() + startHalfSize, freeWalk.end());
        std::reverse(startHalf.begin(), startHalf.end());  // inner first
        startHalf = usable(planIndex, startHalf);
        finishHalf = usable(planIndex, finishHalf);
        if (startHalf.size() < startSlots.size() || finishHalf.size() < finishSlots.size()) {
            throw row_too_small(planIndex, row, demand, std::min(startHalf.size(), startSlots.size()) + std::min(finishHalf.size(), finishSlots.size()));
        }
        // Inner to outer on the start half: the last tap first, so the row reads start, tap 1, ... from -X.
        std::reverse(startSlots.begin(), startSlots.end());
        for (size_t index = 0; index < startSlots.size(); ++index) {
            give(planIndex, startSlots[index], startHalf[index]);
        }
        for (size_t index = 0; index < finishSlots.size(); ++index) {
            give(planIndex, finishSlots[index], finishHalf[index]);
        }
    };
    // All of one half on a side of its own, centred along the row, in slot order from -X.
    auto place_centred = [&](size_t planIndex, size_t row, std::vector<PinSlot*> slots) {
        auto candidates = usable(planIndex, free_pins(row));
        if (slots.size() > candidates.size()) {
            throw row_too_small(planIndex, row, slots.size(), candidates.size());
        }
        const size_t offset = (candidates.size() - slots.size()) / 2;
        for (size_t index = 0; index < slots.size(); ++index) {
            give(planIndex, slots[index], candidates[offset + index]);
        }
    };

    if (planned.size() == 1) {
        const size_t planIndex = planned.front();
        auto& plan = plans[planIndex];
        std::optional<size_t> userRow;
        if (userRowsByPlan.count(planIndex) && userRowsByPlan.at(planIndex).size() == 1) {
            userRow = *userRowsByPlan.at(planIndex).begin();
        }
        const size_t demand = slots_of(plan, true).size() + slots_of(plan, false).size();
        const bool spansBothSides = userRowsByPlan.count(planIndex) && userRowsByPlan.at(planIndex).size() > 1;
        const size_t row = userRow.value_or(0);
        if (!spansBothSides && demand <= usable(planIndex, free_pins(row)).size()) {
            place_on_one_side(planIndex, row);
        }
        else {
            place_centred(planIndex, 0, slots_of(plan, true));
            place_centred(planIndex, 1, slots_of(plan, false));
        }
        return;
    }

    // Two windings: winding k on side k, unless its user pins already put it on the other side.
    std::map<size_t, size_t> rowOfPlan;
    for (auto planIndex : planned) {
        if (!userRowsByPlan.count(planIndex)) {
            continue;
        }
        const auto& rows = userRowsByPlan.at(planIndex);
        if (rows.size() > 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Winding '" + plans[planIndex].name + "' has user pins on both sides of toroid base '" + bobbinName +
                "'; with two windings on the base, each winding takes the two ends of one side (RFC 0013, 'Toroids on a base').");
        }
        rowOfPlan[planIndex] = *rows.begin();
    }
    if (rowOfPlan.size() == 2 && rowOfPlan.at(planned[0]) == rowOfPlan.at(planned[1])) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Windings '" + plans[planned[0]].name + "' and '" + plans[planned[1]].name + "' both have user pins on side " +
            std::to_string(rowOfPlan.at(planned[0])) + " of toroid base '" + bobbinName + "'; each winding takes a side of its own.");
    }
    // Without user pins each winding takes the side its half of the ring faces (ABT #1173): a common-mode
    // choke wound on opposite halves otherwise sends both windings' leads across the ring to the far row.
    // A row's side is the sign of its pins' z; a winding's half is the sign of its turns' y centroid (MKF's
    // y is the MVB++ z). Only when the two windings do not sit on opposite halves (stacked in one sector, or
    // straddling the rows' axis) does neither side face one more than the other, and design order decides.
    if (rowOfPlan.empty()) {
        const double rowZ0 = placedPins[rowWalk[0].front()].centre[2];
        const double rowZ1 = placedPins[rowWalk[1].front()].centre[2];
        const double y0 = plans[planned[0]].meanTurnY;
        const double y1 = plans[planned[1]].meanTurnY;
        if ((rowZ0 < 0) != (rowZ1 < 0) && std::abs(y0) > sameEnd && std::abs(y1) > sameEnd && (y0 < 0) != (y1 < 0)) {
            const size_t rowForNegative = rowZ0 < 0 ? 0 : 1;
            rowOfPlan[planned[0]] = y0 < 0 ? rowForNegative : 1 - rowForNegative;
            rowOfPlan[planned[1]] = 1 - rowOfPlan[planned[0]];
        }
    }
    for (size_t k = 0; k < planned.size(); ++k) {
        if (rowOfPlan.count(planned[k])) {
            continue;
        }
        const size_t other = planned[1 - k];
        rowOfPlan[planned[k]] = rowOfPlan.count(other) ? 1 - rowOfPlan.at(other) : k;
    }
    for (auto planIndex : planned) {
        place_on_one_side(planIndex, rowOfPlan.at(planIndex));
    }
}

}  // namespace

PinAssignmentResult Coil::assign_pins(const Bobbin& bobbin, const Core& core) {
    return assign_pins_impl(bobbin, &core);
}

PinAssignmentResult Coil::assign_pins_impl(const Bobbin& bobbin, const Core* core) {
    PinAssignmentResult result;
    Bobbin bobbinCopy = bobbin;
    if (!bobbinCopy.get_processed_description()) {
        throw CoilNotProcessedException("assign_pins needs a processed bobbin");
    }
    auto processedBobbin = bobbinCopy.get_processed_description().value();
    if (!processedBobbin.get_pins() || processedBobbin.get_pins()->empty()) {
        result.skipped = true;
        result.skippedReason = "the bobbin has no processedDescription.pins[], so winding ends stay at the window border";
        return result;
    }
    // ABT #1173 (WP4): a round window is a toroid on a family t base; its pins go by the toroid rule below.
    const bool toroidBase = bobbinCopy.get_winding_window_shape() == WindingWindowShape::ROUND;
    if (!get_sections_description() || !get_turns_description()) {
        throw CoilNotProcessedException("assign_pins runs after winding: the coil has no sections or turns yet");
    }
    auto placedPins = place_pins(processedBobbin.get_pins().value());
    std::map<std::string, size_t> placedIndexByName;
    for (size_t index = 0; index < placedPins.size(); ++index) {
        if (placedIndexByName.count(placedPins[index].name)) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Bobbin has two pins named '" + placedPins[index].name + "'.");
        }
        placedIndexByName[placedPins[index].name] = index;
    }
    size_t numberRows = 0;
    for (const auto& placedPin : placedPins) {
        numberRows = std::max(numberRows, placedPin.row + 1);
    }
    const bool horizontal = placedPins.front().hangsAlongZ;
    // Each row as pin indexes walked from its START end to its FINISH end.
    std::vector<std::vector<size_t>> rowWalk(numberRows);
    for (size_t index = 0; index < placedPins.size(); ++index) {
        rowWalk[placedPins[index].row].push_back(index);
    }
    for (size_t rowIndex = 0; rowIndex < numberRows; ++rowIndex) {
        // Horizontal: every row starts at -X. Vertical: row 0 at -X, row 1 at +X (diagonal).
        // Toroid base: every row from -X, so all windings start at the same angular side.
        if (!horizontal && !toroidBase && rowIndex % 2 == 1) {
            std::reverse(rowWalk[rowIndex].begin(), rowWalk[rowIndex].end());
        }
    }

    // --- windings in build order -----------------------------------------------------------
    auto sections = get_sections_description().value();
    auto turns = get_turns_description().value();
    auto& windings = get_mutable_functional_description();
    std::set<std::string> shieldWindings;
    for (const auto& section : sections) {
        if (section.get_type() == ElectricalType::SHIELDING) {
            for (const auto& partialWinding : section.get_partial_windings()) {
                shieldWindings.insert(partialWinding.get_winding());
            }
        }
    }
    if (get_layers_description()) {
        const auto coilLayers = get_layers_description().value();  // by value: the getter returns the optional by value
        for (const auto& layer : coilLayers) {
            if (layer.get_type() == ElectricalType::SHIELDING) {
                for (const auto& partialWinding : layer.get_partial_windings()) {
                    shieldWindings.insert(partialWinding.get_winding());
                }
            }
        }
    }

    const double thickRoundDiameter = pin_assignment_number("thickRoundConductingDiameter");
    const auto maximumWiresPerPinThick = static_cast<int64_t>(pin_assignment_number("maximumWiresPerPinThickRound"));
    const auto maximumWiresPerPinFine = static_cast<int64_t>(pin_assignment_number("maximumWiresPerPinFineRound"));

    std::vector<WindingPinPlan> plans;
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        WindingPinPlan plan;
        plan.windingIndex = windingIndex;
        plan.name = windings[windingIndex].get_name();
        plan.side = windings[windingIndex].get_isolation_side();
        plan.numberParallels = windings[windingIndex].get_number_parallels();
        plan.isShield = shieldWindings.count(plan.name) > 0;
        bool found = false;
        bool marginWound = false;
        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            const auto& section = sections[sectionIndex];
            if (section.get_type() == ElectricalType::INSULATION) {
                continue;
            }
            bool contains = false;
            for (const auto& partialWinding : section.get_partial_windings()) {
                if (partialWinding.get_winding() == plan.name) {
                    contains = true;
                }
            }
            if (!contains) {
                continue;
            }
            if (!found) {
                plan.buildRank = sectionIndex;
                found = true;
            }
            ++plan.numberSections;
            auto margin = resolve_margin(section);
            for (double side : margin) {
                if (side > 0) {
                    marginWound = true;
                }
            }
        }
        if (!found) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR,
                "Winding '" + plan.name + "' is in no section, so its build order and ends are unknown to pin assignment.");
        }
        bool firstTurnFound = false;
        for (const auto& turn : turns) {
            if (turn.get_winding() == plan.name) {
                plan.firstTurnY = turn.get_coordinates()[1];
                firstTurnFound = true;
                break;
            }
        }
        if (!firstTurnFound) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR, "Winding '" + plan.name + "' has no turns; pin assignment runs after winding.");
        }
        if (toroidBase) {
            double sumY = 0;
            size_t count = 0;
            bool firstOfParallel0 = true;
            for (const auto& turn : turns) {
                if (turn.get_winding() != plan.name) {
                    continue;
                }
                if (turn.get_coordinate_system() && turn.get_coordinate_system().value() != CoordinateSystem::CARTESIAN) {
                    throw CoilException(ErrorCode::COIL_WINDING_ERROR,
                        "Turn '" + turn.get_name() + "' of toroidal winding '" + plan.name + "' is not in cartesian coordinates; "
                        "the toroid base rule reads which half of the ring each winding occupies from them.");
                }
                sumY += turn.get_coordinates()[1];
                ++count;
                if (turn.get_parallel() == 0) {
                    if (firstOfParallel0) {
                        plan.firstTurnX = turn.get_coordinates()[0];
                        firstOfParallel0 = false;
                    }
                    plan.lastTurnX = turn.get_coordinates()[0];
                }
            }
            plan.meanTurnY = sumY / static_cast<double>(count);
        }
        auto wire = resolve_wire(windingIndex);
        plan.avoidCorners = is_insulated_wire(wire) || marginWound;

        // Strands on one pin while the wraps fit (rule 4); non-round wire has no wrap count in
        // RFC 0013, so its strands never share.
        bool strandsFit = plan.numberParallels == 1;
        if (!strandsFit && wire.get_type() == WireType::ROUND) {
            if (!wire.get_conducting_diameter()) {
                throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
                    "Winding '" + plan.name + "' is round wire without a conducting diameter, so the wires per pin it allows are unknown.");
            }
            double conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
            int64_t maximumWires = conductingDiameter >= thickRoundDiameter ? maximumWiresPerPinThick : maximumWiresPerPinFine;
            strandsFit = plan.numberParallels <= maximumWires;
        }
        plan.separateStrands = !strandsFit;

        // The design's own connections.
        if (windings[windingIndex].get_connections()) {
            plan.originalConnections = windings[windingIndex].get_connections().value();
        }
        bool anyEnd = false;
        bool anyParallel = false;
        bool allPinned = !plan.originalConnections.empty();
        for (const auto& connection : plan.originalConnections) {
            anyEnd = anyEnd || connection.get_end().has_value();
            anyParallel = anyParallel || connection.get_parallel().has_value();
            allPinned = allPinned && connection.get_pin_name().has_value();
        }
        if (!plan.originalConnections.empty() && !anyEnd) {
            if (!allPinned) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Winding '" + plan.name + "' has connections without `end` and without `pinName`; pin assignment cannot tell "
                    "which end of the winding each one terminates. Give every connection an `end` (start, finish or tap).");
            }
            plan.legacy = true;
        }
        if (!plan.originalConnections.empty() && anyEnd) {
            // The design's data decides the strand split: `parallel` absent means all strands together.
            plan.separateStrands = anyParallel;
            for (const auto& connection : plan.originalConnections) {
                if (!connection.get_end()) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "Winding '" + plan.name + "' mixes connections with and without `end`; every connection needs one.");
                }
                if (anyParallel && !connection.get_parallel()) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "Winding '" + plan.name + "' mixes connections with and without `parallel`.");
                }
            }
        }

        if (!plan.legacy) {
            int64_t groups = plan.separateStrands ? plan.numberParallels : 1;
            auto make_slot = [&](End end, int64_t group, size_t tapIndex, bool startHalf) {
                PinSlot slot;
                slot.end = end;
                if (plan.separateStrands) {
                    slot.parallel = group;
                }
                slot.tapIndex = tapIndex;
                slot.startHalf = startHalf;
                return slot;
            };
            for (int64_t group = 0; group < groups; ++group) {
                plan.slots.push_back(make_slot(End::START, group, 0, true));
            }
            for (size_t tapIndex = 1; tapIndex < plan.numberSections; ++tapIndex) {
                for (int64_t group = 0; group < groups; ++group) {
                    plan.slots.push_back(make_slot(End::TAP, group, tapIndex, true));
                }
            }
            for (int64_t group = 0; group < groups; ++group) {
                plan.slots.push_back(make_slot(End::FINISH, group, 0, false));
            }
            // Map the design's connections onto the slots.
            std::map<std::pair<int64_t, int>, size_t> tapsSeen;
            for (size_t connectionIndex = 0; connectionIndex < plan.originalConnections.size(); ++connectionIndex) {
                const auto& connection = plan.originalConnections[connectionIndex];
                End end = connection.get_end().value();
                int64_t parallel = connection.get_parallel().value_or(-1);
                size_t tapIndex = 0;
                if (end == End::TAP) {
                    tapIndex = ++tapsSeen[{parallel, 0}];
                }
                bool mapped = false;
                for (auto& slot : plan.slots) {
                    if (slot.end == end && slot.parallel.value_or(-1) == parallel && slot.tapIndex == tapIndex) {
                        if (slot.connectionIndex) {
                            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                                "Winding '" + plan.name + "' has two connections for its " + end_name(end) +
                                (parallel >= 0 ? " of parallel " + std::to_string(parallel) : std::string("")) + ".");
                        }
                        slot.connectionIndex = connectionIndex;
                        if (connection.get_pin_name()) {
                            slot.pinName = connection.get_pin_name().value();
                            slot.userPin = true;
                        }
                        else if (!is_pin_connection_type(connection.get_type())) {
                            slot.offPin = true;
                        }
                        mapped = true;
                        break;
                    }
                }
                if (!mapped) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "Winding '" + plan.name + "' has a connection for its " + end_name(end) +
                        (parallel >= 0 ? " of parallel " + std::to_string(parallel) : std::string("")) +
                        (end == End::TAP ? " number " + std::to_string(tapIndex) : std::string("")) +
                        " that the winding does not have (" + std::to_string(plan.numberSections) + " series sections, " +
                        std::to_string(plan.numberParallels) + " parallels).");
                }
            }
        }
        else {
            for (const auto& connection : plan.originalConnections) {
                PinSlot slot;
                slot.pinName = connection.get_pin_name().value();
                slot.userPin = true;
                plan.slots.push_back(slot);
            }
        }
        plans.push_back(plan);
    }
    std::stable_sort(plans.begin(), plans.end(), [](const WindingPinPlan& a, const WindingPinPlan& b) {
        return a.buildRank < b.buildRank;
    });

    // --- owners, sides, required creepage (shared by the former and the toroid base) ---------
    std::vector<std::optional<std::pair<size_t, std::string>>> pinOwner(placedPins.size());  // plan index, what
    std::vector<IsolationSide> sidesInBuildOrder;
    for (const auto& plan : plans) {
        if (std::find(sidesInBuildOrder.begin(), sidesInBuildOrder.end(), plan.side) == sidesInBuildOrder.end()) {
            sidesInBuildOrder.push_back(plan.side);
        }
    }
    const size_t numberSides = sidesInBuildOrder.size();

    // Required creepage between isolation sides.
    if (core == nullptr) {
        result.creepageNotCheckedReason = "planned without the core (Coil::terminal_pin_sides): the creepage check decides no pin";
    }
    else if (numberSides < 2) {
        result.creepageNotCheckedReason = "a single isolation side needs no creepage between pins";
    }
    else if (!_inputs) {
        result.creepageNotCheckedReason = "the coil carries no inputs, so the required creepage is unknown";
    }
    else if (!_inputs->get_design_requirements().get_insulation()) {
        result.creepageNotCheckedReason = "the design requirements declare no insulation, so no creepage is required";
    }
    else if (!_inputs->has_insulation_coordination_requirements()) {
        // ABT #1329: an insulation block that names no standard has no creepage table to read.
        result.creepageNotCheckedReason = "the insulation requirements name no standards, so no creepage can be derived";
    }
    else {
        auto inputs = _inputs.value();
        result.requiredCreepage = InsulationCoordinator().calculate_creepage_distance(inputs, false);
    }
    if (core != nullptr && !core->get_processed_description()) {
        throw CoreNotProcessedException("needed for the creepage path through the core when placing pins");
    }
    // A toroid's pins are in the MVB++ toroid frame (ring axis Y, ring in XZ), while the processed core
    // is MKF's ring frame (ring in XY, width = height = outer diameter, depth = ring height): swap.
    // Without a core (terminal_pin_sides) no creepage path is measured, so these are never read.
    const double halfWidth = core != nullptr ? core->get_width() / 2 : std::numeric_limits<double>::quiet_NaN();
    const double halfHeight = core != nullptr ? (toroidBase ? core->get_depth() / 2 : core->get_height() / 2)
                                              : std::numeric_limits<double>::quiet_NaN();
    const double halfDepth = core != nullptr ? (toroidBase ? core->get_height() / 2 : core->get_depth() / 2)
                                             : std::numeric_limits<double>::quiet_NaN();
    auto creepage_between = [&](size_t pinA, size_t pinB) {
        PinCreepagePath path;
        path.pinA = placedPins[pinA].name;
        path.pinB = placedPins[pinB].name;
        path.surfaceDistance = distance(placedPins[pinA].base, placedPins[pinB].base);
        path.throughCoreDistance = distance_to_core_box(placedPins[pinA].base, halfWidth, halfHeight, halfDepth) +
                                   distance_to_core_box(placedPins[pinB].base, halfWidth, halfHeight, halfDepth);
        return path;
    };

    if (toroidBase) {
        assign_pins_on_toroid_base(plans, placedPins, rowWalk, placedIndexByName, pinOwner,
                                   bobbinCopy.get_name().value_or("<unnamed>"));
    }
    else {
        // --- user pins: exist, and one owner each -----------------------------------------------
        std::vector<std::string> knownPins;
        for (const auto& placedPin : placedPins) {
            knownPins.push_back(placedPin.name);
        }
        auto known_pins_text = [&]() {
            std::string text;
            for (const auto& name : knownPins) {
                text += (text.empty() ? "" : ", ") + name;
            }
            return text;
        };
        std::map<IsolationSide, std::set<size_t>> userRowsBySide;
        std::map<IsolationSide, std::pair<std::string, std::string>> firstUserPinBySide;  // pin, winding
        for (size_t planIndex = 0; planIndex < plans.size(); ++planIndex) {
            auto& plan = plans[planIndex];
            for (auto& slot : plan.slots) {
                if (!slot.userPin) {
                    continue;
                }
                auto found = placedIndexByName.find(slot.pinName.value());
                if (found == placedIndexByName.end()) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "Winding '" + plan.name + "' terminates on pin '" + slot.pinName.value() + "', which the bobbin '" +
                        bobbinCopy.get_name().value_or("<unnamed>") + "' does not have. It has: " + known_pins_text() + ".");
                }
                size_t pinIndex = found->second;
                if (pinOwner[pinIndex] && plans[pinOwner[pinIndex]->first].side != plan.side) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "Pin '" + slot.pinName.value() + "' is given to winding '" + plans[pinOwner[pinIndex]->first].name +
                        "' and to winding '" + plan.name + "', which are on different isolation sides.");
                }
                pinOwner[pinIndex] = std::make_pair(planIndex, plan.name);
                userRowsBySide[plan.side].insert(placedPins[pinIndex].row);
                if (!firstUserPinBySide.count(plan.side)) {
                    firstUserPinBySide[plan.side] = {slot.pinName.value(), plan.name};
                }
            }
        }

        // --- rule 1: rows are isolation groups ----------------------------------------------------
        std::map<IsolationSide, size_t> rowOfSide;
        std::vector<std::vector<IsolationSide>> sidesOfRow(numberRows);
        auto pin_on_row_text = [&](IsolationSide side) {
            auto [pinName, windingName] = firstUserPinBySide.at(side);
            return "pin '" + pinName + "' (winding '" + windingName + "', side " + side_name(side) + ")";
        };
        for (auto side : sidesInBuildOrder) {
            if (!userRowsBySide.count(side)) {
                continue;
            }
            const auto& rows = userRowsBySide.at(side);
            if (rows.size() > 1 && numberSides <= numberRows) {
                // A side spread over two rows while every side could have its own row.
                std::string rowList;
                for (auto row : rows) {
                    rowList += (rowList.empty() ? "" : " and ") + std::to_string(row);
                }
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Isolation side " + side_name(side) + " has user pins on rows " + rowList + " of bobbin '" +
                    bobbinCopy.get_name().value_or("<unnamed>") + "' (" + pin_on_row_text(side) + " among them); rows are "
                    "isolation groups (RFC 0013 rule 1), so one side's pins belong on one row.");
            }
            size_t row = *rows.begin();
            if (numberSides <= numberRows && !sidesOfRow[row].empty()) {
                auto other = sidesOfRow[row].front();
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Wrong row: " + pin_on_row_text(side) + " is on row " + std::to_string(row) + ", which already holds " +
                    pin_on_row_text(other) + ". The bobbin has " + std::to_string(numberRows) + " rows for " +
                    std::to_string(numberSides) + " isolation sides, so each side gets a row of its own (RFC 0013 rule 1).");
            }
            rowOfSide[side] = row;
            sidesOfRow[row].push_back(side);
        }
        for (auto side : sidesInBuildOrder) {
            if (rowOfSide.count(side)) {
                continue;
            }
            // First-wound winding of this side decides which flange-side row is nearer on a
            // horizontal former (rows sit one per end flange, at y = -+rowDistance/2).
            double firstTurnY = 0;
            for (const auto& plan : plans) {
                if (plan.side == side) {
                    firstTurnY = plan.firstTurnY;
                    break;
                }
            }
            std::optional<size_t> chosen;
            for (size_t row = 0; row < numberRows; ++row) {
                if (numberSides <= numberRows && !sidesOfRow[row].empty()) {
                    continue;
                }
                if (!chosen) {
                    chosen = row;
                    continue;
                }
                // More sides than rows: fewest sides first. Otherwise, on a horizontal former, the
                // row nearer the flange the winding starts at; ties keep the lower row.
                if (numberSides > numberRows) {
                    if (sidesOfRow[row].size() < sidesOfRow[chosen.value()].size()) {
                        chosen = row;
                    }
                }
                else if (horizontal) {
                    double rowY = placedPins[rowWalk[row].front()].centre[1];
                    double chosenY = placedPins[rowWalk[chosen.value()].front()].centre[1];
                    if (std::abs(rowY - firstTurnY) < std::abs(chosenY - firstTurnY) - pinCoordinateTolerance) {
                        chosen = row;
                    }
                }
            }
            if (!chosen) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "No row is left for isolation side " + side_name(side) + ".");
            }
            rowOfSide[side] = chosen.value();
            sidesOfRow[chosen.value()].push_back(side);
        }
        // Keep each row's sides in build order.
        for (auto& rowSides : sidesOfRow) {
            std::stable_sort(rowSides.begin(), rowSides.end(), [&](IsolationSide a, IsolationSide b) {
                return std::find(sidesInBuildOrder.begin(), sidesInBuildOrder.end(), a) <
                       std::find(sidesInBuildOrder.begin(), sidesInBuildOrder.end(), b);
            });
        }

        // --- shield finishes: on a spare pin when one exists, else buried --------------------------
        auto pins_needed_on_row = [&](size_t row, bool withShieldFinishes) {
            size_t needed = 0;
            for (const auto& plan : plans) {
                if (plan.legacy || rowOfSide.at(plan.side) != row) {
                    continue;
                }
                for (const auto& slot : plan.slots) {
                    if (slot.userPin || slot.offPin) {
                        continue;
                    }
                    if (plan.isShield && slot.end == End::FINISH && !withShieldFinishes) {
                        continue;
                    }
                    ++needed;
                }
            }
            return needed;
        };
        auto free_pins_on_row = [&](size_t row) {
            size_t freePins = 0;
            for (auto pinIndex : rowWalk[row]) {
                if (!pinOwner[pinIndex]) {
                    ++freePins;
                }
            }
            return freePins;
        };
        for (size_t row = 0; row < numberRows; ++row) {
            if (pins_needed_on_row(row, true) <= free_pins_on_row(row)) {
                continue;
            }
            for (auto& plan : plans) {
                if (!plan.isShield || plan.legacy || rowOfSide.at(plan.side) != row) {
                    continue;
                }
                for (auto& slot : plan.slots) {
                    if (slot.end == End::FINISH && !slot.userPin && !slot.offPin) {
                        slot.offPin = true;
                        result.notes.push_back("Shield winding '" + plan.name + "' has no spare pin on row " + std::to_string(row) +
                                               ", so its finish is buried as a flying lead (WP8 R14 prefers both ends on pins).");
                    }
                }
            }
        }

        // --- allocate, row by row, block by block ------------------------------------------------
        for (size_t row = 0; row < numberRows; ++row) {
            const auto& rowSides = sidesOfRow[row];
            if (rowSides.empty()) {
                continue;
            }
            const auto& walk = rowWalk[row];
            // Free pins of the row in walk order, and each side's demand.
            std::vector<size_t> demand;
            for (auto side : rowSides) {
                size_t sideDemand = 0;
                for (const auto& plan : plans) {
                    if (plan.legacy || plan.side != side) {
                        continue;
                    }
                    for (const auto& slot : plan.slots) {
                        if (!slot.userPin && !slot.offPin) {
                            ++sideDemand;
                        }
                    }
                }
                demand.push_back(sideDemand);
            }
            std::vector<size_t> freeWalk;
            for (auto pinIndex : walk) {
                if (!pinOwner[pinIndex]) {
                    freeWalk.push_back(pinIndex);
                }
            }
            size_t totalDemand = 0;
            for (auto sideDemand : demand) {
                totalDemand += sideDemand;
            }
            if (totalDemand > freeWalk.size()) {
                std::string sideList;
                for (size_t sideIndex = 0; sideIndex < rowSides.size(); ++sideIndex) {
                    sideList += (sideList.empty() ? "" : ", ") + side_name(rowSides[sideIndex]) + " needs " + std::to_string(demand[sideIndex]);
                }
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Row " + std::to_string(row) + " of bobbin '" + bobbinCopy.get_name().value_or("<unnamed>") + "' has " +
                    std::to_string(freeWalk.size()) + " free pins but " + sideList + " (" + std::to_string(totalDemand) +
                    " in all). Use a former with more pins or fewer separate strands and taps.");
            }
            // Blocks along the row, sides in build order, the spare pins spread over the gaps so the
            // sides sit as far apart as the row allows.
            size_t spare = freeWalk.size() - totalDemand;
            size_t numberGaps = rowSides.size() - 1;
            std::vector<std::pair<size_t, size_t>> blocks;  // [first, last+1) into freeWalk
            size_t cursor = 0;
            for (size_t sideIndex = 0; sideIndex < rowSides.size(); ++sideIndex) {
                size_t blockSize = demand[sideIndex];
                if (rowSides.size() == 1) {
                    blockSize = freeWalk.size();  // a row of its own: the whole row is the block
                }
                blocks.push_back({cursor, cursor + blockSize});
                cursor += blockSize;
                if (sideIndex < numberGaps) {
                    // Pin COUNTS: the spare pins split over the gaps, the remainder to the first gaps.
                    size_t gap = spare / numberGaps + (sideIndex < spare % numberGaps ? 1 : 0);
                    cursor += gap;
                }
            }
            for (size_t sideIndex = 0; sideIndex < rowSides.size(); ++sideIndex) {
                auto side = rowSides[sideIndex];
                std::vector<size_t> block(freeWalk.begin() + blocks[sideIndex].first, freeWalk.begin() + blocks[sideIndex].second);
                // Demand per half, over this side's windings.
                size_t startDemand = 0;
                size_t finishDemand = 0;
                for (const auto& plan : plans) {
                    if (plan.legacy || plan.side != side) {
                        continue;
                    }
                    for (const auto& slot : plan.slots) {
                        if (slot.userPin || slot.offPin) {
                            continue;
                        }
                        (slot.startHalf ? startDemand : finishDemand) += 1;
                    }
                }
                // Split the block so the used pins sit centred in it: the spare pins go half to
                // each outer end (pin counts, the odd one to the finish end).
                size_t blockSpare = block.size() - startDemand - finishDemand;
                size_t startHalfSize = startDemand + blockSpare / 2;
                std::vector<size_t> startHalf(block.begin(), block.begin() + startHalfSize);   // walk order
                std::vector<size_t> finishHalf(block.begin() + startHalfSize, block.end());     // walk order
                std::reverse(startHalf.begin(), startHalf.end());  // inner (centre) first
                size_t startCursor = 0;
                size_t finishCursor = 0;
                for (auto& plan : plans) {
                    if (plan.legacy || plan.side != side) {
                        continue;
                    }
                    // Start half, inner to outer: the last tap first, the first start last, so that
                    // walking the row from its start end reads start, tap 1, ..., tap K-1.
                    std::vector<PinSlot*> startOrder;
                    std::vector<PinSlot*> finishOrder;
                    for (auto& slot : plan.slots) {
                        if (slot.userPin || slot.offPin) {
                            continue;
                        }
                        (slot.startHalf ? startOrder : finishOrder).push_back(&slot);
                    }
                    std::reverse(startOrder.begin(), startOrder.end());
                    auto take = [&](std::vector<size_t>& half, size_t& halfCursor, PinSlot* slot, const std::string& halfName) {
                        while (halfCursor < half.size()) {
                            size_t pinIndex = half[halfCursor];
                            ++halfCursor;
                            if (plan.avoidCorners && placedPins[pinIndex].corner) {
                                continue;  // rule 7: left free, never given to a stripped/margin lead
                            }
                            slot->pinName = placedPins[pinIndex].name;
                            pinOwner[pinIndex] = std::make_pair(&plan - plans.data(), plan.name);
                            return;
                        }
                        throw InvalidInputException(ErrorCode::INVALID_INPUT,
                            "Row " + std::to_string(row) + " of bobbin '" + bobbinCopy.get_name().value_or("<unnamed>") +
                            "' has no free pin left on its " + halfName + " half for the " + end_name(slot->end) +
                            " of winding '" + plan.name + "'" +
                            (plan.avoidCorners ? " (its insulated or margin-wound leads may not take a corner pin, RFC 0013 rule 7)" : std::string("")) +
                            ". Nesting windings from the inner pins outward keeps leads from crossing (rule 6), so the pins cannot be borrowed from the other half.");
                    };
                    for (auto* slot : startOrder) {
                        take(startHalf, startCursor, slot, "start");
                    }
                    for (auto* slot : finishOrder) {
                        take(finishHalf, finishCursor, slot, "finish");
                    }
                }
            }
        }

    }

    // --- creepage between every pair of sides ----------------------------------------------
    for (size_t pinA = 0; core != nullptr && pinA < placedPins.size(); ++pinA) {
        if (!pinOwner[pinA]) {
            continue;
        }
        for (size_t pinB = pinA + 1; pinB < placedPins.size(); ++pinB) {
            if (!pinOwner[pinB]) {
                continue;
            }
            auto sideA = plans[pinOwner[pinA]->first].side;
            auto sideB = plans[pinOwner[pinB]->first].side;
            if (sideA == sideB) {
                continue;
            }
            auto path = creepage_between(pinA, pinB);
            path.sideA = side_name(sideA);
            path.sideB = side_name(sideB);
            if (!result.worstPath || path.get_path() < result.worstPath->get_path()) {
                result.worstPath = path;
            }
        }
    }
    if (result.requiredCreepage && result.worstPath && result.worstPath->get_path() < result.requiredCreepage.value() - 1e-9) {
        const auto& worst = result.worstPath.value();
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Creepage between isolation sides " + worst.sideA + " and " + worst.sideB + " is short at pins '" + worst.pinA +
            "' and '" + worst.pinB + "' of bobbin '" + bobbinCopy.get_name().value_or("<unnamed>") + "': the shortest path is " +
            format_millimetres(worst.get_path()) + " (surface " + format_millimetres(worst.surfaceDistance) + ", through the core " +
            format_millimetres(worst.throughCoreDistance) + ") against " + format_millimetres(result.requiredCreepage.value()) +
            " required, a shortfall of " + format_millimetres(result.requiredCreepage.value() - worst.get_path()) + ".");
    }

    // --- write the connections -------------------------------------------------------------
    for (auto& plan : plans) {
        if (plan.legacy) {
            continue;
        }
        std::vector<ConnectionElement> connections;
        for (auto& slot : plan.slots) {
            ConnectionElement connection;
            if (slot.connectionIndex) {
                connection = plan.originalConnections[slot.connectionIndex.value()];
            }
            connection.set_end(slot.end);
            if (slot.parallel) {
                connection.set_parallel(slot.parallel.value());
            }
            if (slot.offPin) {
                if (!connection.get_type()) {
                    connection.set_type(ConnectionType::FLYING_LEAD);
                }
            }
            else {
                const auto& placedPin = placedPins[placedIndexByName.at(slot.pinName.value())];
                connection.set_pin_name(slot.pinName.value());
                if (!connection.get_type()) {
                    connection.set_type(placedPin.type == MAS::PinDescriptionType::SMD ? ConnectionType::SMT : ConnectionType::THT);
                }
                if (!connection.get_diameter()) {
                    connection.set_diameter(placedPin.diameter);
                }
            }
            connections.push_back(connection);
        }
        if (plan.separateStrands && plan.numberParallels > 1) {
            result.notes.push_back("Winding '" + plan.name + "': " + std::to_string(plan.numberParallels) +
                                   " strands terminate on separate adjacent pins, shorted on the PCB.");
        }
        windings[plan.windingIndex].set_connections(connections);
    }
    return result;
}

Coil::TerminalPinSides Coil::terminal_pin_sides() {
    TerminalPinSides sides;
    auto bobbin = resolve_bobbin();
    if (!bobbin.get_processed_description() || !bobbin.get_processed_description()->get_pins() ||
        bobbin.get_processed_description()->get_pins()->empty()) {
        return sides;   // (iii) no pins: the S0/S1 default
    }
    const auto bobbinPins = bobbin.get_processed_description()->get_pins().value();
    // (ii) ends with no pin named: the pins assign_pins gives them, planned on a copy. The copy is
    // never written back: with coil_connect_leads_to_pins off the design keeps its connections and
    // its leads end in free air, on the face and half those pins put them.
    std::optional<Coil> assigned;
    const Coil* source = this;
    if (!has_complete_pin_connections()) {
        assigned = *this;
        assigned->assign_pins_impl(bobbin, nullptr);
        source = &assigned.value();
    }
    const auto placedPins = place_pins(bobbinPins);
    std::map<std::string, const PlacedPin*> placedByName;
    for (const auto& placed : placedPins) {
        placedByName[placed.name] = &placed;
    }
    auto millimetres = [](double metres) { return format_millimetres(metres); };
    // The first pinned end seen of each winding, to name it when another end disagrees on the face.
    std::map<std::string, std::string> firstEndByWinding;
    const auto& windings = source->get_functional_description();
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        const auto& winding = windings[windingIndex];
        if (!winding.get_connections()) {
            continue;
        }
        const auto connections = winding.get_connections().value();
        for (const auto& connection : connections) {
            // The same ends Coil::get_connection_reserved_spaces routes to pins (terminalPinByConductor):
            // a pinName with its end; taps are junctions inside the coil.
            if (!connection.get_pin_name() || !connection.get_end() || connection.get_end().value() == End::TAP) {
                continue;
            }
            const std::string pinName = connection.get_pin_name().value();
            auto found = placedByName.find(pinName);
            if (found == placedByName.end()) {
                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                    "Winding '" + winding.get_name() + "' ends on pin '" + pinName + "', which is not among the bobbin's pins.");
            }
            const auto& pin = *found->second;
            const double z = pin.centre[2];
            if (std::abs(z) < 1e-9) {
                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                    "Pin '" + pinName + "' of winding '" + winding.get_name() + "' stands at z = " + millimetres(z) +
                    ", on neither side of the column, so the connection face its lead leaves on is unknown.");
            }
            const int face = z < 0 ? 0 : 1;
            const bool entrance = connection.get_end().value() == End::START;
            const std::string endText = std::string(entrance ? "start" : "finish") + " on pin '" + pinName + "' (x = " +
                                        millimetres(pin.centre[0]) + ", z = " + millimetres(z) + ", row " + std::to_string(pin.row) + ")";
            auto known = sides.faceByWinding.find(winding.get_name());
            if (known != sides.faceByWinding.end() && known->second != face) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Winding '" + winding.get_name() + "' has its ends on pins on both sides of the column: its " +
                    firstEndByWinding.at(winding.get_name()) + " leaves on face " + std::to_string(known->second) + ", its " + endText +
                    " on face " + std::to_string(face) + ". One winding's leads leave on one connection face (spec section 12 "
                    "covers no winding pinned to both rails).");
            }
            sides.faceByWinding[winding.get_name()] = face;
            if (!firstEndByWinding.count(winding.get_name())) {
                firstEndByWinding[winding.get_name()] = endText;
            }
            // The half in the face's own frame: face 1 is the bobbin frame turned 180 degrees about the column.
            const double faceX = face == 0 ? pin.centre[0] : -pin.centre[0];
            const int half = std::abs(faceX) < 1e-9 ? 0 : (faceX < 0 ? -1 : 1);
            std::vector<int64_t> parallels;
            if (connection.get_parallel()) {
                parallels.push_back(connection.get_parallel().value());
            }
            else {
                for (int64_t parallel = 0; parallel < int64_t(source->get_number_parallels(windingIndex)); ++parallel) {
                    parallels.push_back(parallel);
                }
            }
            for (int64_t parallel : parallels) {
                sides.halfByEnd[{winding.get_name(), parallel, entrance}] = half;
                sides.pinByEnd[{winding.get_name(), parallel, entrance}] = pinName;
            }
        }
    }
    return sides;
}

}  // namespace OpenMagnetics
