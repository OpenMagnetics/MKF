// The terminal-lead exit of a CLOSED core -- a core whose ferrite closes all the way
// around the winding, so no azimuth is free for a lead to leave through. Kept out of Coil.cpp so the
// winding code and this plan can be read (and reverted) separately.
//
// The exit follows what manufacturers build (Alf 2026-09-24; checked online):
//
//   DRUM + RING, ring SEATED on a square base plate (the common shielded SMD
//   power-inductor construction). The drum's bottom flange is a square plate of side M (corner chamfer
//   leg N) extending past the round ring at its corners; the ring stands on the plate and has a
//   NOTCH at its foot per wire end; the wire leaves the winding through the notch and runs over the
//   plate's corner area to the terminal metal. Planned here, per winding end:
//     in-window run  the terminal lead MKF already lays on the plate-side edge row (Coil.cpp forces
//                    that edge), trimmed to end at r_arc = the outermost layer's outer radius plus
//                    MKF's ride-over there (compute_ride_levels: a return laid under the outer
//                    face lifts it) plus half the lead's coated diameter (the lead lies against the
//                    winding as wound);
//     arc            round the column axis at r_arc, at the lead's own height, from the lead's exit
//                    slot on the connection plane (-Z) to the plate-corner diagonal on its side:
//                    entrance -> the (-X, -Z) corner, exit -> the (+X, -Z) corner. The leads of one
//                    end (N-filar) take one lane each, r_arc + k x pitch, leading lead innermost;
//     drop           an outermost-layer end leaves the winding at its own turn's height, not on an
//                    edge row: it descends there, in the band between winding and ring, to lie on the
//                    plate (an edge-row lead keeps its row);
//     notch run      radially along the diagonal, through the notch, to J/2 + its coated radius;
//     corner run     straight over the plate's corner area to the terminal plate.
//   Notch: one per winding end, on that diagonal, through the ring wall (K/2 .. J/2); WIDTH = the
//   end's parallels x the largest coated OD (the Coil::terminal_exit_slots pitch rule: side by side,
//   one pitch apart); HEIGHT = the tallest of its leads above the plate (for an edge-row lead that is
//   its row's edgeDepth, the row allocator's own number).
//   Terminal plate (ELECTRODE): the L-shaped winding/mounting terminal -- a mounting
//   pad padWidth x padDepth on the plate's bottom face (connections[].padWidth / padDepth, type smt)
//   and a leg up the plate's side face that the wire end is soldered onto. Stock thickness and
//   material from src/data/dfm_rules.json "leadExitStock" (sourced; throws when missing). The
//   entrance's terminal is on the -X side, the exit's on the +X side, each at the connection-plane
//   CORNER END of its side: on a seated-ring plate the only plate a lead can reach from a notch is the
//   corner area outside the ring disc (the ring covers the side's middle when M = J).
//
//   MOLDED body. The coil ends
//   are resistance-welded to a metal LEAD FRAME pressed into the body with the coil; the frame
//   carries the current to terminal pads on the side and bottom faces. Planned here, per winding end:
//     in-window run  trimmed to r_tab = outermost layer's outer radius + its ride-over + half the coated diameter
//                    (in the powder when that is outside the cavity E: the body is pressed round the
//                    coil and its leads, so E bounds the winding, not the leads);
//     arc            at r_tab on the lead's own height, from its exit slot to its TAB: each tab on
//                    its own side, one pitch clear of the other end's slot (its descent passes that
//                    end's run), and the two frames' inner ends at least one stock thickness apart
//                    (a sheet-metal slot is at least one stock thickness wide);
//     drop           to lie on the cavity floor at the tab.
//   LEAD_FRAME: the strip in the frame plane directly under the cavity floor, y in [-D/2 - t, -D/2],
//   from the tab (width = the coated OD) to the side face x = -A/2 (entrance) or +A/2 (exit)
//   (width = padWidth, centred on that face), then an
//   ELECTRODE side leg down that face and an ELECTRODE bottom pad padWidth x padDepth (the board
//   face), both of the frame's own stock. WELD: the nugget at the tab, one coated diameter long and
//   wide, from the frame to the top of the wire lying on it (stated rule: the weld nugget is one wire
//   diameter long).
//
// Every dimension comes from the wire (coated OD, planned bend radius), the core's letters, the
// connection's pad data or the sourced dfm entry; the only numbers written here are the geometric
// tolerance the exit slots already use (1e-9 m) and the arc sampling that tolerance implies. Missing
// data and constructions the plan cannot lay out at all THROW naming the lead and the numbers; copper
// it lays where it collides (with the ring, the winding, the plate's edge, the body's side or another
// lead) is laid exactly there and RECORDED (LeadExitClash), so the design is drawn as built and
// analyze_results reports it -- never refused (main's draw-and-report rule, Alf 2026-10-24/25).
// Nothing is shrunk, moved or skipped to fit.
//
// NOT modelled (stated, not guessed): the frame's and terminal plate's resistance (MKF has no
// conductivity for "phosphor bronze"; routedLength carries the WIRE runs only), several windings
// sharing a closed core (two plate corners / two frames carry ONE winding's two ends), N-filar ends
// on a molded lead frame (one frame per end is one wire), a column that is not round.

#include "constructive_models/Coil.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <cmrc/cmrc.hpp>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>

CMRC_DECLARE(dfmData);

namespace OpenMagnetics {

namespace {

// The clearance tolerance Coil::terminal_exit_slots and the pin runs check copper against.
constexpr double clearanceTolerance = 1e-9;

std::string mm(double metres) {
    std::ostringstream stream;
    stream.setf(std::ios::fixed);
    stream.precision(3);
    stream << metres * 1000.0 << " mm";
    return stream.str();
}

// src/data/dfm_rules.json "leadExitStock": the sourced stock the terminal metal is made of.
const json& lead_exit_stock() {
    static const json block = [] {
        auto fs = cmrc::dfmData::get_filesystem();
        auto data = fs.open("src/data/dfm_rules.json");
        auto all = json::parse(std::string(data.begin(), data.end()));
        if (!all.contains("leadExitStock")) {
            throw std::runtime_error("src/data/dfm_rules.json has no leadExitStock block: the lead frame / terminal plate "
                                     "stock of a closed core's lead exit has no sourced value");
        }
        return all.at("leadExitStock");
    }();
    return block;
}

double lead_exit_stock_thickness(const std::string& key) {
    const auto& block = lead_exit_stock();
    if (!block.contains(key) || !block.at(key).is_number()) {
        throw std::runtime_error("src/data/dfm_rules.json leadExitStock has no sourced number '" + key +
                                 "': the stock thickness is not a default MKF may invent");
    }
    const double value = block.at(key).get<double>();
    if (!(value > 0)) {
        throw std::runtime_error("src/data/dfm_rules.json leadExitStock '" + key + "' must be a positive thickness");
    }
    return value;
}

std::string lead_exit_stock_material() {
    const auto& block = lead_exit_stock();
    if (!block.contains("material") || !block.at("material").is_string()) {
        throw std::runtime_error("src/data/dfm_rules.json leadExitStock has no material name");
    }
    return block.at("material").get<std::string>();
}

// Point (R sin a, y, -R cos a): azimuth a measured from -Z towards +X (the pinWaypoints frame).
std::vector<double> at_azimuth(double radius, double azimuth, double axial) {
    return {radius * std::sin(azimuth), axial, -radius * std::cos(azimuth)};
}

// The arc sampled finely enough that no chord strays more than the clearance tolerance from it:
// sagitta R (1 - cos(d/2)) <= tol  <=>  d <= 2 acos(1 - tol / R). For CHECKS only; the route
// carries the exact arc.
std::vector<std::vector<double>> sample_arc(double radius, double from, double to, double axial) {
    const double step = 2 * std::acos(1 - clearanceTolerance / radius);
    const size_t segments = std::max<size_t>(1, size_t(std::ceil(std::abs(to - from) / step)));
    std::vector<std::vector<double>> points;
    for (size_t k = 0; k <= segments; ++k) {
        points.push_back(at_azimuth(radius, from + (to - from) * double(k) / double(segments), axial));
    }
    return points;
}

double distance(const std::vector<double>& a, const std::vector<double>& b) {
    return std::sqrt(std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) + std::pow(a[2] - b[2], 2));
}

// Smallest distance from the column axis (x = z = 0) to the XZ projection of segment a -> b.
double axis_distance_xz(const std::vector<double>& a, const std::vector<double>& b) {
    const double dx = b[0] - a[0];
    const double dz = b[2] - a[2];
    const double length2 = dx * dx + dz * dz;
    double t = length2 > 0 ? -(a[0] * dx + a[2] * dz) / length2 : 0.0;
    t = std::clamp(t, 0.0, 1.0);
    return std::hypot(a[0] + t * dx, a[2] + t * dz);
}

ExtrudedPolygon box_xz(double x0, double x1, double z0, double z1, double y0, double y1) {
    return ExtrudedPolygon{{{std::min(x0, x1), std::min(z0, z1)}, {std::max(x0, x1), std::min(z0, z1)},
                            {std::max(x0, x1), std::max(z0, z1)}, {std::min(x0, x1), std::max(z0, z1)}},
                           std::min(y0, y1), std::max(y0, y1)};
}

struct ExitLead {
    size_t route;
    size_t space;
    std::string label;
    std::string winding;
    int64_t parallel;
    bool entrance;
    double diameter;        // coated (or sleeve) outer diameter
    double bendRadius;      // ConnectionRoute::plannedBendRadius
    double slot;            // ConnectionRoute::exitSlot
    double axial;           // the axial level of the in-window run where it ends
    bool edgeRow;           // laid on an edge row (edgeDepth > 0), not an outermost-layer radial exit
    size_t borderIndex;     // index of the in-window run's outer end in route.waypoints
};

}  // namespace

std::vector<LeadExitSolid> Coil::plan_closed_core_lead_exits(std::vector<ConnectionRoute>& routes,
                                                             std::vector<ConnectionReservedSpace>& spaces,
                                                             const std::map<size_t, size_t>& terminalSpaceIndexByRoute,
                                                             const std::vector<ConnectionRideLevel>& rideLevels,
                                                             std::vector<LeadExitClash>& clashes) {
    std::vector<LeadExitSolid> solids;
    // Copper the plan lays where it collides (with the core, the winding, the plate's edge or
    // another lead) is still laid exactly where the plan puts it, and recorded: the design is drawn
    // as built and the clash reported (analyze_results), never refused.
    auto clash = [&](LeadExitClash::Kind kind, const std::string& winding, const std::string& message) {
        clashes.push_back(LeadExitClash{kind, winding, message});
    };
    if (!_coreLeadExit) {
        return solids;
    }
    const auto& exit = _coreLeadExit.value();
    auto bobbin = resolve_bobbin();
    if (bobbin.get_winding_window_shape() != WindingWindowShape::RECTANGULAR) {
        throw NotImplementedException("A closed core's lead exit is planned for a concentric (rectangular) winding window only");
    }
    const auto processedBobbin = bobbin.get_processed_description().value();
    if (processedBobbin.get_column_shape() != ColumnShape::ROUND) {
        throw NotImplementedException("A closed core's lead exit is planned around a ROUND column only: its arcs run round "
                                      "the column axis at a constant radius");
    }
    if (get_functional_description().size() != 1) {
        throw NotImplementedException("A closed core's lead exit carries ONE winding's two ends (two plate corners, two lead "
                                      "frames); this coil has " + std::to_string(get_functional_description().size()) + " windings");
    }
    auto layers = get_layers_description();
    auto turns = get_turns_description();
    if (!layers || !turns) {
        throw CoilNotProcessedException("A closed core's lead exit needs the wound layers and turns");
    }

    // The winding's outer radius: the outermost layer's outer face (every layer: an insulation
    // wrap outside the copper is still in the way).
    double windingOuterRadius = std::numeric_limits<double>::lowest();
    for (const auto& layer : layers.value()) {
        if (layer.get_orientation() != WindingOrientation::OVERLAPPING) {
            throw NotImplementedException("A closed core's lead exit is planned for OVERLAPPING layers only (layer '" +
                                          layer.get_name() + "' is contiguous)");
        }
        windingOuterRadius = std::max(windingOuterRadius, layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2);
    }

    // The winding's terminal leads, as MKF laid them.
    std::vector<ExitLead> leads;
    auto wires = get_wires();
    for (size_t index = 0; index < routes.size(); ++index) {
        auto& route = routes[index];
        if (route.kind != ConnectionKind::TERMINAL_ENTRANCE && route.kind != ConnectionKind::TERMINAL_EXIT) {
            continue;
        }
        ExitLead lead;
        lead.route = index;
        lead.entrance = route.kind == ConnectionKind::TERMINAL_ENTRANCE;
        lead.winding = route.winding;
        lead.parallel = route.parallel;
        lead.label = "winding '" + route.winding + "' parallel " + std::to_string(route.parallel) + (lead.entrance ? " entrance" : " exit");
        if (!route.pinName.empty()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The terminal lead of " + lead.label + " is assigned to bobbin pin '" + route.pinName + "', but the core is closed: "
                "its leads end on the core's own terminals (lead exit), not on pins");
        }
        auto spaceIt = terminalSpaceIndexByRoute.find(index);
        if (spaceIt == terminalSpaceIndexByRoute.end()) {
            throw std::logic_error("The terminal route of " + lead.label + " has no reservation to carry its exit run");
        }
        lead.space = spaceIt->second;
        auto wire = wires[get_winding_index_by_name(route.winding)];
        lead.diameter = std::max({wire.get_maximum_outer_width(), wire.get_maximum_outer_height(), route.sleeveOuterDiameter.value_or(0.0)});
        lead.bendRadius = route.plannedBendRadius;
        if (!(lead.bendRadius >= lead.diameter / 2 - clearanceTolerance)) {
            throw std::logic_error("The terminal route of " + lead.label + " carries no planned bend radius");
        }
        if (route.tangentDeparture) {
            // The tangential departure laid this lead off the lateral face (no slot on the connection plane, no
            // in-window run to end), and the exit run starts where that run ends.
            throw NotImplementedException("The terminal lead of " + lead.label + " leaves tangentially off the lateral face "
                                          "(ConnectionRoute::tangentDeparture); a closed core's lead exit is planned from a "
                                          "connection-plane exit slot only");
        }
        if (!route.exitSlot) {
            throw std::logic_error("The terminal route of " + lead.label + " has no exit slot (Coil::terminal_exit_slots)");
        }
        lead.slot = route.exitSlot.value();
        // The in-window run's OUTER end: waypoints are in electrical order, so an entrance starts at
        // the terminal end and an exit ends there.
        lead.borderIndex = lead.entrance ? 0 : route.waypoints.size() - 1;
        lead.axial = route.waypoints[lead.borderIndex][1];
        lead.edgeRow = spaces[lead.space].edgeDepth > 0;
        leads.push_back(lead);
    }
    if (leads.empty()) {
        return solids;
    }

    // The window floor = the plate's top face / the cavity floor, in the coil's frame.
    const double floor = exit.kind == CoreLeadExitGeometry::Kind::RING_NOTCH ? exit.ringNotch->plateTop : exit.leadFrame->cavityFloor;
    const auto window = processedBobbin.get_winding_windows()[0];
    const double windowCentre = window.get_coordinates().value()[1];

    for (const auto& lead : leads) {
        const double radius = lead.diameter / 2;
        if (lead.axial - radius < floor - clearanceTolerance) {
            clash(LeadExitClash::Kind::BELOW_FLOOR, lead.winding,
                "the terminal lead of " + lead.label + " runs at " + mm(lead.axial) + ", its copper below the plate / cavity floor at " + mm(floor));
        }
        if (lead.edgeRow) {
            // Coil.cpp forces a closed core's terminal rows onto the plate-side edge. A row that is
            // not there means that forcing did not reach it: refuse rather than route it anyway.
            if (!(lead.axial < windowCentre)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The terminal lead of " + lead.label + " is laid on the edge row at " + mm(lead.axial) +
                    ", on the top edge: a closed core's leads leave over its base, so its row must be at the plate edge");
            }
            // Its stub from the turn to that row runs inside the turn's own layer, at the turn's x:
            // refuse when another turn of the SAME conductor in that layer lies on it (its centre
            // within one coated diameter of the stub) -- which is what forcing the plate-side edge
            // on an end whose own turns stand between it and the plate would produce. The other
            // parallels of an N-filar bundle share the row and are separated along the flange by
            // their exit slots (Coil::terminal_exit_slots); their stubs past each other's turns are
            // the bundle's ordinary layout, not this guard's business.
            const auto& route = routes[lead.route];
            const std::string turnName = lead.entrance ? route.toTurn : route.fromTurn;
            std::optional<Turn> attached;
            for (const auto& turn : turns.value()) {
                if (turn.get_name() == turnName) {
                    attached = turn;
                    break;
                }
            }
            if (!attached || !attached->get_layer()) {
                throw std::logic_error("The terminal route of " + lead.label + " attaches to turn '" + turnName + "', which has no layer");
            }
            const std::string layerName = attached->get_layer().value();
            const double stubLow = std::min(lead.axial, attached->get_coordinates()[1]);
            const double stubHigh = std::max(lead.axial, attached->get_coordinates()[1]);
            for (const auto& turn : turns.value()) {
                if (turn.get_name() == turnName || !turn.get_layer() || turn.get_layer().value() != layerName ||
                    turn.get_winding() != lead.winding || turn.get_parallel() != attached->get_parallel()) {
                    continue;
                }
                const double axial = turn.get_coordinates()[1];
                const double gap = std::abs(axial - std::clamp(axial, stubLow, stubHigh));
                if (gap < lead.diameter - clearanceTolerance && stubHigh - stubLow > clearanceTolerance) {
                    clash(LeadExitClash::Kind::STUB_CROSSES_OWN_TURN, lead.winding,
                        "the terminal lead of " + lead.label + " attaches to turn '" + turnName + "' at " +
                        mm(attached->get_coordinates()[1]) + " and drops to its plate-side row at " + mm(lead.axial) +
                        " through turn '" + turn.get_name() + "' of its own layer '" + layerName + "' (at " + mm(axial) +
                        "): the stub crosses it (a closed core's leads leave over its base)");
                }
            }
        }
    }

    // One pitch per winding: the largest coated OD among its leads (terminal_exit_slots' rule).
    double pitch = 0;
    for (const auto& lead : leads) {
        pitch = std::max(pitch, lead.diameter);
    }
    // The lead lies against the winding: its centre one coated radius outside the winding's face --
    // the face as WOUND, i.e. lifted by MKF's own ride-over (compute_ride_levels) where a return lays
    // a lane under it: a turn at the outer face riding over such a lane stands one ride further out,
    // and the lead leaving past it must clear it there, not at the unlifted face. The ride is read at
    // the winding's outer face on the face (side) of each lead's route, and the leads of one END take
    // the largest of theirs: they share one set of lanes, one pitch apart (below).
    ConnectionLayout rides;
    rides.rideLevels = rideLevels;
    std::map<bool, double> rideByEnd;   // keyed by `entrance`
    for (const auto& lead : leads) {
        const int side = routes[lead.route].side;
        if (side != 0 && side != 1) {
            throw std::logic_error("The terminal route of " + lead.label + " routes on face " + std::to_string(side) +
                                   ", which has no ride levels (0 = the primary isolation side's face, 1 = every other): the ride "
                                   "its lane must clear cannot be determined");
        }
        const double ride = rides.ride_at(windingOuterRadius, side);
        auto [it, inserted] = rideByEnd.emplace(lead.entrance, ride);
        if (!inserted) {
            it->second = std::max(it->second, ride);
        }
    }
    auto ride_of = [&](bool entrance) {
        auto it = rideByEnd.find(entrance);
        if (it == rideByEnd.end()) {
            throw std::logic_error(std::string("The ") + (entrance ? "entrance" : "exit") + " end has no lead, so it has no ride");
        }
        return it->second;
    };
    auto lane_radius_of = [&](bool entrance) { return windingOuterRadius + ride_of(entrance) + pitch / 2; };
    // Names the lane's terms in an error message.
    auto lane_terms = [&](bool entrance) {
        return "winding outer radius " + mm(windingOuterRadius) + " + ride-over " + mm(ride_of(entrance)) + " + half its " + mm(pitch) + " pitch";
    };

    // Trim each in-window run to end ON the lane circle, at its slot x: (slot, axial, -radial) with
    // slot^2 + radial^2 = laneRadius^2. The run beyond it (to the old window border) is replaced by
    // the exit run; the difference is taken off routedLength.
    auto start_radial = [&](const ExitLead& lead, double lane) {
        if (!(std::abs(lead.slot) < lane)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The exit slot of " + lead.label + " (x = " + mm(lead.slot) + ") lies outside the lane round the winding (radius " +
                mm(lane) + ")");
        }
        return std::sqrt(lane * lane - lead.slot * lead.slot);
    };
    auto trim_in_window_run = [&](const ExitLead& lead, double lane) {
        auto& route = routes[lead.route];
        auto& border = route.waypoints[lead.borderIndex];
        const double newRadial = start_radial(lead, lane);
        const size_t neighbour = lead.entrance ? 1 : route.waypoints.size() - 2;
        if (!(newRadial > route.waypoints[neighbour][0])) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The in-window run of " + lead.label + " would end at radius " + mm(newRadial) + ", inside its own last corner at " +
                mm(route.waypoints[neighbour][0]));
        }
        const double change = newRadial - border[0];
        border[0] = newRadial;
        route.routedLength = roundFloat(route.routedLength + change, 9);
        auto& space = spaces[lead.space];
        if (!space.routedLength) {
            throw std::logic_error("The reservation of " + lead.label + " carries no routedLength");
        }
        space.routedLength = roundFloat(space.routedLength.value() + change, 9);
    };
    auto add_exit_run = [&](const ExitLead& lead, std::vector<std::vector<double>> run, ConnectionRoute::ExitRunArc arc,
                            ConnectionRoute::ExitRunEnd end) {
        auto& route = routes[lead.route];
        double length = arc.radius * std::abs(arc.toAzimuth - arc.fromAzimuth);
        for (size_t k = 1; k + 1 < run.size(); ++k) {
            length += distance(run[k], run[k + 1]);
        }
        route.exitRun = std::move(run);
        route.exitRunArc = arc;
        route.exitRunEnd = end;
        route.exitRunRide = ride_of(lead.entrance);
        route.routedLength = roundFloat(route.routedLength + length, 9);
        auto& space = spaces[lead.space];
        space.routedLength = roundFloat(space.routedLength.value() + length, 9);
    };
    // The reserved-space pad data of a winding end (decision 2: connections[].padWidth / padDepth).
    auto pad_of = [&](const ExitLead& lead) -> std::pair<double, double> {
        // By VALUE: get_functional_description() and get_connections() return copies.
        const auto winding = get_functional_description()[get_winding_index_by_name(lead.winding)];
        const End wanted = lead.entrance ? End::START : End::FINISH;
        const auto connections = winding.get_connections();
        if (connections) {
            for (const auto& connection : connections.value()) {
                if (!connection.get_end() || connection.get_end().value() != wanted) {
                    continue;
                }
                if (connection.get_parallel() && connection.get_parallel().value() != lead.parallel) {
                    continue;
                }
                if (!connection.get_pad_width() || !connection.get_pad_depth()) {
                    break;
                }
                if (connection.get_type() && connection.get_type().value() != ConnectionType::SMT) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "The " + std::string(lead.entrance ? "start" : "finish") + " connection of winding '" + lead.winding +
                        "' states a pad but is not an smt terminal; a closed core's terminal is an SMD pad");
                }
                return {connection.get_pad_width().value(), connection.get_pad_depth().value()};
            }
        }
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Winding '" + lead.winding + "' has no connection for its " + (lead.entrance ? "start" : "finish") +
            " end stating padWidth and padDepth: the size of a closed core's terminal pad is the connection's data "
            "(connections[].padWidth / padDepth, type smt), and MKF does not invent one");
    };
    const std::string stockMaterial = lead_exit_stock_material();

    // Full 3D path of a lead for the pairwise clearance check: the in-window run at its slot, the
    // exit run with its arc sampled.
    std::vector<std::vector<std::vector<double>>> paths(leads.size());
    auto record_path = [&](size_t which) {
        const auto& lead = leads[which];
        const auto& route = routes[lead.route];
        std::vector<std::vector<double>> inWindow;
        for (const auto& waypoint : route.waypoints) {
            inWindow.push_back({lead.slot, waypoint[1], -waypoint[0]});
        }
        if (!lead.entrance) {
            // exit: turn -> border; keep that order, the exit run continues from the border.
        }
        else {
            std::reverse(inWindow.begin(), inWindow.end());   // entrance waypoints run terminal -> turn
        }
        auto path = inWindow;
        const auto& arc = route.exitRunArc.value();
        auto arcPoints = sample_arc(arc.radius, arc.fromAzimuth, arc.toAzimuth, arc.axial);
        path.insert(path.end(), arcPoints.begin() + 1, arcPoints.end());
        path.insert(path.end(), route.exitRun.begin() + 2, route.exitRun.end());
        paths[which] = path;
    };

    if (exit.kind == CoreLeadExitGeometry::Kind::RING_NOTCH) {
        const auto& ring = exit.ringNotch.value();
        const double ringInnerRadius = ring.ringInnerDiameter / 2;
        const double ringOuterRadius = ring.ringOuterDiameter / 2;
        const double half = ring.plateSide / 2;
        const double plateBottom = ring.plateTop - ring.plateThickness;
        const double terminalThickness = lead_exit_stock_thickness("terminalPlateThickness");
        // The band between the winding and the ring bore must hold the lanes AND their bends. The
        // leads of one end take one lane each, one pitch apart (below), so that end's outermost lane
        // is its lane radius (ride included) + (its leads - 1) x pitch.
        for (const auto& lead : leads) {
            size_t endCount = 0;
            for (const auto& other : leads) {
                endCount += other.entrance == lead.entrance;
            }
            const double outermostLane = lane_radius_of(lead.entrance) + double(endCount - 1) * pitch;
            const double needed = outermostLane + std::max(lead.diameter / 2, lead.bendRadius);
            if (needed > ringInnerRadius + clearanceTolerance) {
                clash(LeadExitClash::Kind::LANE_PAST_RING_BORE, lead.winding,
                    "the lead of " + lead.label + " does not fit round the winding inside the ring: its outermost lane at " + mm(outermostLane) +
                    " (" + lane_terms(lead.entrance) + ", one pitch per lead of an end) plus its " +
                    mm(std::max(lead.diameter / 2, lead.bendRadius)) + " bend reaches " + mm(needed) +
                    ", past the ring bore K/2 = " + mm(ringInnerRadius));
            }
        }
        for (bool entrance : {true, false}) {
            std::vector<size_t> members;
            for (size_t which = 0; which < leads.size(); ++which) {
                if (leads[which].entrance == entrance) {
                    members.push_back(which);
                }
            }
            if (members.empty()) {
                continue;
            }
            const double side = entrance ? -1.0 : 1.0;
            // Several leads of one end (an N-filar bundle) share an edge row, side by side along the
            // flange at their exit slots, so their arcs cannot share one circle. Each takes its own
            // LANE, one pitch further out: the LEADING lead (the one already furthest towards its
            // corner) the innermost lane and the FAR line through the notch, each following lead the
            // next lane out and the next line nearer -- so a following lead's arc passes outside
            // every lead ahead of it and ends before their radial runs through the notch begin.
            std::stable_sort(members.begin(), members.end(), [&](size_t a, size_t b) {
                if (std::abs(leads[a].slot - leads[b].slot) > clearanceTolerance) {
                    return side * leads[a].slot > side * leads[b].slot;
                }
                return leads[a].parallel < leads[b].parallel;
            });
            const double diagonal = side * std::numbers::pi / 4;   // towards the (side X, -Z) corner
            const std::array<double, 2> along = {std::sin(diagonal), -std::cos(diagonal)};   // (x, z), outward
            const std::array<double, 2> across = {std::cos(diagonal), std::sin(diagonal)};   // perpendicular
            const double notchWidth = double(members.size()) * pitch;
            if (!(notchWidth / 2 < ringInnerRadius)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The notch for winding '" + leads[members[0]].winding + "'s " + (entrance ? "entrance" : "exit") + " would be " +
                    mm(notchWidth) + " wide, not narrower than the ring bore");
            }
            // The terminal plate on this side: on the side face x = side * M/2, from the connection-
            // plane corner (past its chamfer) along the side by the pad's width.
            auto [padWidth, padDepth] = pad_of(leads[members[0]]);
            for (size_t k = 1; k < members.size(); ++k) {
                auto other = pad_of(leads[members[k]]);
                if (std::abs(other.first - padWidth) > clearanceTolerance || std::abs(other.second - padDepth) > clearanceTolerance) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "The parallels of winding '" + leads[members[0]].winding + "' state different pads for one terminal");
                }
            }
            const double sideStart = -(half - ring.plateChamferLeg);
            const double sideEnd = sideStart + padWidth;
            if (sideEnd > half - ring.plateChamferLeg + clearanceTolerance) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The pad of winding '" + leads[members[0]].winding + "' (" + mm(padWidth) + " wide) is longer than the plate's side between its chamfers");
            }
            if (padDepth > half) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The pad of winding '" + leads[members[0]].winding + "' (" + mm(padDepth) + " deep) reaches past the plate's centre");
            }
            if (padWidth < double(members.size()) * pitch) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The terminal of winding '" + leads[members[0]].winding + "' is " + mm(padWidth) + " wide, less than its " +
                    std::to_string(members.size()) + " lead(s) side by side (" + mm(double(members.size()) * pitch) + ")");
            }
            // Positive `across` points away from the connection plane for the exit and towards it
            // for the entrance, so side * offset > 0 is the FAR line.
            auto offset_of = [&](size_t position) {
                return side * (double(members.size() - 1) / 2 - double(position)) * pitch;
            };
            // Where each lead crosses onto the plate: just outside the ring, on its notch line, at the
            // height it passes the notch (an edge-row lead keeps its row; an outermost-layer end has
            // descended to lie on the plate).
            const double legCentreX = side * (half + terminalThickness / 2);
            std::vector<std::vector<double>> crossing(members.size());
            for (size_t position = 0; position < members.size(); ++position) {
                const auto& lead = leads[members[position]];
                const double radius = lead.diameter / 2;
                const double offset = offset_of(position);
                const double notchAxial = lead.edgeRow ? lead.axial : ring.plateTop + radius;
                const double lineOutside = std::sqrt(std::pow(ringOuterRadius + radius, 2) - offset * offset);
                crossing[position] = {offset * across[0] + lineOutside * along[0], notchAxial, offset * across[1] + lineOutside * along[1]};
            }
            // Where each lands on the terminal's leg (soldered onto its top, lying on it), and how it
            // gets there: a STAIRCASE. The crossing highest along the leg (z) runs straight across the
            // plate to the leg (along x) at its own z; each lower one first continues along its notch
            // line -- parallel to the one above, so still one pitch from it -- until it is one pitch
            // below that one's landing, then runs across to the leg. Every run is then one pitch from
            // its neighbours along its whole length (the pitch rule), with no search and no number
            // of our own. A landing that falls off the leg, or a turn that reaches the leg first,
            // throws.
            std::vector<size_t> alongLeg(members.size());
            for (size_t position = 0; position < members.size(); ++position) {
                alongLeg[position] = position;
            }
            std::stable_sort(alongLeg.begin(), alongLeg.end(), [&](size_t a, size_t b) { return crossing[a][2] > crossing[b][2]; });
            std::vector<double> landingZ(members.size());
            std::vector<std::optional<std::vector<double>>> turnOntoLeg(members.size());
            auto target_at = [&](size_t position, double z) {
                return std::vector<double>{legCentreX, ring.plateTop + leads[members[position]].diameter / 2, z};
            };
            std::optional<double> landingAbove;
            for (size_t position : alongLeg) {
                const auto& lead = leads[members[position]];
                const double radius = lead.diameter / 2;
                double landing = crossing[position][2];
                if (landingAbove) {
                    landing = std::min(landing, landingAbove.value() - pitch);
                }
                landing = std::min(landing, sideEnd - radius);
                if (landing < sideStart + radius - clearanceTolerance) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "The lead of " + lead.label + " cannot land on its terminal's leg (" + mm(sideStart) + " .. " + mm(sideEnd) +
                        ") one pitch beside the lead(s) of its end above it: it would land at z = " + mm(landing));
                }
                if (landing < crossing[position][2] - clearanceTolerance) {
                    // Continue along the notch line (z falls by along[1] per unit) down to the landing's z.
                    const double run = (landing - crossing[position][2]) / along[1];
                    std::vector<double> turn = {crossing[position][0] + run * along[0], crossing[position][1], landing};
                    if (side * turn[0] > half - radius + clearanceTolerance) {
                        throw InvalidInputException(ErrorCode::INVALID_INPUT,
                            "The lead of " + lead.label + " reaches the plate's side edge before it is one pitch beside the "
                            "lead(s) of its end above it");
                    }
                    turnOntoLeg[position] = turn;
                }
                landingZ[position] = landing;
                landingAbove = landing;
            }
            double notchHeight = 0;
            for (size_t position = 0; position < members.size(); ++position) {
                const auto& lead = leads[members[position]];
                const double radius = lead.diameter / 2;
                const double offset = offset_of(position);
                const double lane = lane_radius_of(entrance) + double(position) * pitch;
                if (!(std::abs(offset) < lane)) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "The notch lane of " + lead.label + " is offset " + mm(offset) + ", beyond the lane round the winding");
                }
                trim_in_window_run(lead, lane);
                const double startAzimuth = std::atan2(lead.slot, start_radial(lead, lane));
                // Where this member's line through the notch meets its lane circle.
                const double lineAtLane = std::sqrt(lane * lane - offset * offset);
                const std::vector<double> arcEndXz = {offset * across[0] + lineAtLane * along[0], offset * across[1] + lineAtLane * along[1]};
                const double endAzimuth = std::atan2(arcEndXz[0], -arcEndXz[1]);
                if (side * (endAzimuth - startAzimuth) < 0) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "The exit slot of " + lead.label + " lies beyond its plate-corner diagonal: its arc would run backwards");
                }
                // An edge-row lead keeps its row; an outermost-layer end descends to lie on the plate.
                const double notchAxial = crossing[position][1];
                notchHeight = std::max(notchHeight, notchAxial + radius - ring.plateTop);
                const std::vector<double> outside = crossing[position];
                const std::vector<double> target = target_at(position, landingZ[position]);
                // The corner run must stay on the free plate: outside the ring disc by a radius, and
                // inside the chamfer and the plate's far edge.
                const double ringClearance = std::min(axis_distance_xz(outside, turnOntoLeg[position].value_or(outside)),
                                                      axis_distance_xz(turnOntoLeg[position].value_or(outside), target));
                if (ringClearance < ringOuterRadius + radius - clearanceTolerance) {
                    clash(LeadExitClash::Kind::CORNER_RUN_OVER_RING, lead.winding,
                        "the corner run of " + lead.label + " from its notch to its terminal passes " + mm(ringClearance) +
                        " from the axis, through the ring (J/2 + its radius = " + mm(ringOuterRadius + radius) + ")");
                }
                // Where it leaves the plate: over the plate top inside the square by its radius and
                // inside the chamfer line side*x - z = M - N by its radius (a line at 45 degrees:
                // distance (M - N - (side*x - z)) / sqrt(2)); then across the side edge x = side*M/2
                // onto the terminal's leg, on the flat side past the chamfer.
                for (const auto& point : {outside, turnOntoLeg[position].value_or(outside)}) {
                    if (std::abs(point[0]) > half - radius + clearanceTolerance || std::abs(point[2]) > half - radius + clearanceTolerance ||
                        side * point[0] - point[2] > ring.plateSide - ring.plateChamferLeg - radius * std::numbers::sqrt2 + clearanceTolerance) {
                        clash(LeadExitClash::Kind::CORNER_RUN_OFF_PLATE, lead.winding,
                            "the corner run of " + lead.label + " at (" + mm(point[0]) + ", " + mm(point[2]) +
                            ") is not on the plate's corner area (plate side M = " + mm(ring.plateSide) + ", chamfer leg N = " +
                            mm(ring.plateChamferLeg) + ")");
                    }
                }
                const std::vector<double> lastCorner = turnOntoLeg[position].value_or(outside);
                if (std::abs(target[0] - lastCorner[0]) > 0) {
                    const double t = (side * half - lastCorner[0]) / (target[0] - lastCorner[0]);
                    const double crossingZ = lastCorner[2] + t * (target[2] - lastCorner[2]);
                    if (t < 0 || t > 1 || crossingZ < sideStart + radius - clearanceTolerance || crossingZ > sideEnd - radius + clearanceTolerance) {
                        clash(LeadExitClash::Kind::CORNER_RUN_OFF_TERMINAL, lead.winding,
                            "the corner run of " + lead.label + " crosses the plate's side edge at z = " + mm(crossingZ) +
                            ", not over its terminal (z " + mm(sideStart) + " .. " + mm(sideEnd) + ")");
                    }
                }
                std::vector<std::vector<double>> run;
                run.push_back(at_azimuth(lane, startAzimuth, lead.axial));
                run.push_back({arcEndXz[0], lead.axial, arcEndXz[1]});
                if (std::abs(notchAxial - lead.axial) > clearanceTolerance) {
                    run.push_back({arcEndXz[0], notchAxial, arcEndXz[1]});
                }
                run.push_back(outside);
                if (turnOntoLeg[position]) {
                    run.push_back(turnOntoLeg[position].value());
                }
                run.push_back(target);
                add_exit_run(lead, run, ConnectionRoute::ExitRunArc{lane, startAzimuth, endAzimuth, lead.axial},
                             ConnectionRoute::ExitRunEnd::ELECTRODE);
                // The joint between the wire and its terminal. A wound part's coil end is
                // joined to its terminal electrode by THERMOCOMPRESSION bonding: the insulating film of
                // the end portion is removed, and the end is pressed while heated, plastically deformed
                // and bonded to the electrode. The wire was
                // laid on the leg's top face at its COATED radius, so its bare conductor stops one
                // coating thickness short of the metal; the bond closes that gap with the wire's own
                // conductor: from the leg's top face up to the wire's centreline (where the bare
                // conductor is its full width), the bare conductor's width across the run, over the
                // length of the run that lies on the leg (from the plate's side edge x = side*M/2 to
                // the run's end). Every number is the wire's or the terminal's; none is chosen here.
                {
                    const auto& wire = wires[get_winding_index_by_name(lead.winding)];
                    if (wire.get_type() != WireType::ROUND && wire.get_type() != WireType::LITZ) {
                        throw NotImplementedException("The bond of " + lead.label + " to its terminal is modelled for a round "
                                                      "conductor (round or litz wire), whose bare width across the run is its "
                                                      "conducting diameter; this wire is not round");
                    }
                    const double bareDiameter = wire.get_maximum_conducting_width();
                    if (!(bareDiameter > 0) || bareDiameter > lead.diameter + clearanceTolerance) {
                        throw InvalidInputException(ErrorCode::INVALID_INPUT,
                            "The wire of " + lead.label + " states a conducting diameter of " + mm(bareDiameter) +
                            ", which cannot lie inside its " + mm(lead.diameter) + " coated diameter");
                    }
                    const double dx = target[0] - lastCorner[0];
                    const double dz = target[2] - lastCorner[2];
                    if (!(side * dx > 0)) {
                        throw std::logic_error("The run of " + lead.label + " onto its terminal does not travel towards it");
                    }
                    const double t = (side * half - lastCorner[0]) / dx;
                    const std::array<double, 2> onLeg = {side * half, lastCorner[2] + t * dz};
                    const double length = std::hypot(dx, dz);
                    const std::array<double, 2> normal = {-dz / length * bareDiameter / 2, dx / length * bareDiameter / 2};
                    ExtrudedPolygon bond;
                    bond.xz = {{onLeg[0] - normal[0], onLeg[1] - normal[1]},
                               {target[0] - normal[0], target[2] - normal[1]},
                               {target[0] + normal[0], target[2] + normal[1]},
                               {onLeg[0] + normal[0], onLeg[1] + normal[1]}};
                    bond.y0 = ring.plateTop;
                    bond.y1 = target[1];
                    const std::string end = entrance ? "entrance" : "exit";
                    solids.push_back(LeadExitSolid{LeadExitSolid::Kind::BOND,
                                                   lead.winding + " " + end + " parallel " + std::to_string(lead.parallel) + " bond",
                                                   lead.winding, lead.parallel, entrance, bond,
                                                   resolve_wire(get_winding_index_by_name(lead.winding)).resolve_material().get_name(),
                                                   std::nullopt});
                }
            }
            const auto& first = leads[members[0]];
            const std::optional<int64_t> servesParallel = members.size() == 1 ? std::optional<int64_t>(first.parallel) : std::nullopt;
            const std::string end = entrance ? "entrance" : "exit";
            // The notch: a slot of the notch's width along the diagonal, through the ring wall. Its
            // inner end reaches the bore across the whole width (the bore meets the slot's edges at
            // sqrt((K/2)^2 - (w/2)^2)), its outer end is the ring's outer face.
            {
                const double inner = std::sqrt(ringInnerRadius * ringInnerRadius - notchWidth * notchWidth / 4);
                const double outer = ringOuterRadius;
                ExtrudedPolygon body;
                for (auto [s, o] : std::vector<std::pair<double, double>>{{inner, -notchWidth / 2}, {outer, -notchWidth / 2},
                                                                        {outer, notchWidth / 2}, {inner, notchWidth / 2}}) {
                    body.xz.push_back({o * across[0] + s * along[0], o * across[1] + s * along[1]});
                }
                body.y0 = ring.plateTop;
                body.y1 = ring.plateTop + notchHeight;
                solids.push_back(LeadExitSolid{LeadExitSolid::Kind::CORE_CUT, first.winding + " " + end + " ring notch", first.winding,
                                               servesParallel, entrance, body, "", std::nullopt});
            }
            // The terminal plate: its side leg (x = side*M/2 .. + t, from under the plate to its top)
            // and its mounting pad (padDepth in from the side, t thick under the plate: the board face).
            {
                const double faceX = side * half;
                const double outerX = side * (half + terminalThickness);
                ExtrudedPolygon leg = box_xz(faceX, outerX, sideStart, sideEnd, plateBottom - terminalThickness, ring.plateTop);
                solids.push_back(LeadExitSolid{LeadExitSolid::Kind::ELECTRODE, first.winding + " " + end + " terminal leg", first.winding,
                                               servesParallel, entrance, leg, stockMaterial, std::nullopt});
                ExtrudedPolygon pad = box_xz(side * (half - padDepth), outerX, sideStart, sideEnd, plateBottom - terminalThickness, plateBottom);
                ExtrudedPolygon board = pad;
                board.y1 = board.y0;
                solids.push_back(LeadExitSolid{LeadExitSolid::Kind::ELECTRODE, first.winding + " " + end + " terminal pad", first.winding,
                                               servesParallel, entrance, pad, stockMaterial, board});
            }
        }
    }
    else {
        const auto& frame = exit.leadFrame.value();
        const double frameThickness = lead_exit_stock_thickness("moldedLeadFrameThickness");
        if (frameThickness > frame.basePlateThickness) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "The lead frame (" + mm(frameThickness) + " stock, src/data/dfm_rules.json leadExitStock) is thicker than the "
                "molded body's base under the cavity ((B - D)/2 = " + mm(frame.basePlateThickness) + ")");
        }
        // The molded body is pressed AROUND the coil and its leads (powder and resin pressed
        // around coil and lead frame): outside the winding cavity E is powder, not a wall,
        // and the leads run in it. What bounds them is the body itself: the lane and its bends
        // must stay inside the body's footprint.
        const double footprintHalf = std::min(frame.bodyWidth, frame.bodyDepth) / 2;
        // (The bends -- radial run onto the arc, arc down onto the tab -- round their corners on the
        // INSIDE, so the lane's own coated radius is the outermost copper.)
        for (const auto& lead : leads) {
            const double laneRadius = lane_radius_of(lead.entrance);
            const double needed = laneRadius + lead.diameter / 2;
            if (needed > footprintHalf + clearanceTolerance) {
                clash(LeadExitClash::Kind::LANE_PAST_BODY_SIDE, lead.winding,
                    "the lead of " + lead.label + " runs to its lead-frame tab outside the molded body: its lane at " + mm(laneRadius) +
                    " (" + lane_terms(lead.entrance) + ") reaches " + mm(needed) + " with its coated radius, past the body's side at " +
                    mm(footprintHalf));
            }
        }
        size_t entranceIndex = leads.size();
        size_t exitIndex = leads.size();
        for (size_t which = 0; which < leads.size(); ++which) {
            size_t& slotIndex = leads[which].entrance ? entranceIndex : exitIndex;
            if (slotIndex != leads.size()) {
                throw NotImplementedException("A molded lead frame carries one wire per winding end; winding '" + leads[which].winding +
                                              "' has several parallels at its " + (leads[which].entrance ? "start" : "finish"));
            }
            slotIndex = which;
        }
        if (entranceIndex == leads.size() || exitIndex == leads.size()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "A molded lead frame needs both winding ends' leads");
        }
        const auto& entranceLead = leads[entranceIndex];
        const auto& exitLead = leads[exitIndex];
        // Where the tabs go. Each end's tab is on its own side (entrance -X, exit +X) and
        //   (a) one pitch clear of the OTHER end's slot: a lead descending to its tab passes the
        //       other end's in-window run where it leaves the winding, and two leads stand at least
        //       one pitch apart side by side (the Coil::terminal_exit_slots rule);
        //   (b) the two frames' inner ends -- each tab reaching one coated radius past its weld
        //       centre -- stand at least one stock thickness apart (a sheet-metal slot is at least
        //       one stock thickness wide); when they do not, both tabs move out symmetrically.
        double entranceTab = std::min(entranceLead.slot, exitLead.slot - pitch);
        double exitTab = std::max(exitLead.slot, entranceLead.slot + pitch);
        const double separation = pitch + frameThickness;
        if (exitTab - entranceTab < separation) {
            const double shortfall = separation - (exitTab - entranceTab);
            entranceTab -= shortfall / 2;
            exitTab += shortfall / 2;
        }
        for (auto [leadIndex, tab] : std::vector<std::pair<size_t, double>>{{entranceIndex, entranceTab}, {exitIndex, exitTab}}) {
            const auto& lead = leads[leadIndex];
            const double radius = lead.diameter / 2;
            const double side = lead.entrance ? -1.0 : 1.0;
            const double laneRadius = lane_radius_of(lead.entrance);
            if (!(std::abs(tab) < laneRadius)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The lead-frame tab of " + lead.label + " (x = " + mm(tab) + ") lies outside the lane round the winding");
            }
            auto [padWidth, padDepth] = pad_of(lead);
            if (padWidth < lead.diameter) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The pad of " + lead.label + " (" + mm(padWidth) + " wide) is narrower than the lead frame's tab (one coated OD, " +
                    mm(lead.diameter) + ")");
            }
            trim_in_window_run(lead, laneRadius);
            const double startAzimuth = std::atan2(lead.slot, start_radial(lead, laneRadius));
            const double tabZ = -std::sqrt(laneRadius * laneRadius - tab * tab);
            const double endAzimuth = std::atan2(tab, -tabZ);
            const double weldAxial = frame.cavityFloor + radius;
            std::vector<std::vector<double>> run;
            run.push_back(at_azimuth(laneRadius, startAzimuth, lead.axial));
            run.push_back({tab, lead.axial, tabZ});
            if (std::abs(lead.axial - weldAxial) > clearanceTolerance) {
                run.push_back({tab, weldAxial, tabZ});
            }
            add_exit_run(lead, run, ConnectionRoute::ExitRunArc{laneRadius, startAzimuth, endAzimuth, lead.axial},
                         ConnectionRoute::ExitRunEnd::LEAD_FRAME_WELD);
            // The terminal is centred on its side face (one terminal per opposite side,
            // centred), padWidth along it: the body's depth must hold it.
            if (padWidth / 2 > frame.bodyDepth / 2 + clearanceTolerance) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The pad of " + lead.label + " (" + mm(padWidth) + " wide) is wider than the body's side (C = " + mm(frame.bodyDepth) + ")");
            }
            if (padDepth > frame.bodyWidth / 2) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "The pad of " + lead.label + " (" + mm(padDepth) + " deep) reaches past the body's centre");
            }
            const double faceX = side * frame.bodyWidth / 2;
            const double tabInnerX = tab - side * radius;   // the tab's end towards the other frame
            const std::string end = lead.entrance ? "entrance" : "exit";
            const double frameTop = frame.cavityFloor;
            const double frameBottom = frame.cavityFloor - frameThickness;
            const double bodyBottom = frame.cavityFloor - frame.basePlateThickness;
            // The strip: one coated OD wide at the tab, padWidth wide where it reaches the side face.
            ExtrudedPolygon strip{{{tabInnerX, tabZ - radius}, {tabInnerX, tabZ + radius}, {faceX, padWidth / 2}, {faceX, -padWidth / 2}},
                                  frameBottom, frameTop};
            solids.push_back(LeadExitSolid{LeadExitSolid::Kind::LEAD_FRAME, lead.winding + " " + end + " lead frame", lead.winding,
                                           lead.parallel, lead.entrance, strip, stockMaterial, std::nullopt});
            ExtrudedPolygon weld = box_xz(tab - radius, tab + radius, tabZ - radius, tabZ + radius, frameTop, frameTop + lead.diameter);
            solids.push_back(LeadExitSolid{LeadExitSolid::Kind::WELD, lead.winding + " " + end + " weld", lead.winding, lead.parallel,
                                           lead.entrance, weld, resolve_wire(get_winding_index_by_name(lead.winding)).resolve_material().get_name(),
                                           std::nullopt});
            const double outerX = side * (frame.bodyWidth / 2 + frameThickness);
            ExtrudedPolygon leg = box_xz(faceX, outerX, -padWidth / 2, padWidth / 2, bodyBottom - frameThickness, frameTop);
            solids.push_back(LeadExitSolid{LeadExitSolid::Kind::ELECTRODE, lead.winding + " " + end + " terminal leg", lead.winding,
                                           lead.parallel, lead.entrance, leg, stockMaterial, std::nullopt});
            ExtrudedPolygon pad = box_xz(side * (frame.bodyWidth / 2 - padDepth), outerX, -padWidth / 2, padWidth / 2,
                                         bodyBottom - frameThickness, bodyBottom);
            ExtrudedPolygon board = pad;
            board.y1 = board.y0;
            solids.push_back(LeadExitSolid{LeadExitSolid::Kind::ELECTRODE, lead.winding + " " + end + " terminal pad", lead.winding,
                                           lead.parallel, lead.entrance, pad, stockMaterial, board});
        }
        // The frames' inner ends, measured.
        const double gap = (exitTab - pitch / 2) - (entranceTab + pitch / 2);
        if (gap < frameThickness - clearanceTolerance) {
            throw std::logic_error("The two lead frames stand " + mm(gap) + " apart, less than one stock thickness");
        }
    }

    // Every pair of leads keeps its coated clearance along the whole path, the exit runs included.
    for (size_t which = 0; which < leads.size(); ++which) {
        record_path(which);
    }
    for (size_t i = 0; i < leads.size(); ++i) {
        for (size_t j = i + 1; j < leads.size(); ++j) {
            const double found = lead_polyline_distance(paths[i], paths[j]);
            const double required = (leads[i].diameter + leads[j].diameter) / 2;
            if (found < required - clearanceTolerance) {
                clash(LeadExitClash::Kind::LEADS_SHARE_COPPER, leads[i].winding,
                    "the exit runs of " + leads[i].label + " and " + leads[j].label + " share copper: centreline distance " + mm(found) +
                    " < " + mm(required));
            }
        }
    }
    return solids;
}

}  // namespace OpenMagnetics
