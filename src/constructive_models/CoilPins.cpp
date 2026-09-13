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
//   * Round winding windows (toroids on a family `t` base) are WP4's and throw here.
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
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <cmrc/cmrc.hpp>
#include <limits>
#include <map>
#include <set>
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

PinLeadRoute Coil::route_lead_to_pin(const MAS::Pin& pin, const std::vector<double>& windowExit, double frontFaceOffset) {
    if (windowExit.size() < 2) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "route_lead_to_pin needs the window exit as {x, y}.");
    }
    auto placed = place_pins({pin}).front();
    // The lead leaves the coil on the core's open front face (-Z, where a two-piece core has no
    // lateral leg and where MVB++ runs every terminal lead parallel to Z): MKF's radial exit
    // coordinate becomes the depth -(radial + frontFaceOffset) there, its axial coordinate stays y.
    const double ex = windowExit[0] + frontFaceOffset;
    const double ey = windowExit[1];
    std::vector<std::vector<double>> points;
    if (placed.hangsAlongZ) {
        const double zBase = placed.base[2];
        points = {{0.0, ey, -ex}, {0.0, ey, zBase}, {0.0, placed.centre[1], zBase}, {placed.centre[0], placed.centre[1], zBase}};
    }
    else {
        const double yBase = placed.base[1];
        points = {{0.0, ey, -ex}, {0.0, yBase, -ex}, {0.0, yBase, placed.centre[2]}, {placed.centre[0], yBase, placed.centre[2]}};
    }
    PinLeadRoute route;
    route.pinName = placed.name;
    for (auto& point : points) {
        if (route.waypoints.empty() || distance(route.waypoints.back(), point) > 1e-12) {
            route.waypoints.push_back(point);
        }
    }
    for (size_t index = 0; index + 1 < route.waypoints.size(); ++index) {
        route.length += distance(route.waypoints[index], route.waypoints[index + 1]);
    }
    return route;
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

PinAssignmentResult Coil::assign_pins(const Bobbin& bobbin, const Core& core) {
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
    if (bobbinCopy.get_winding_window_shape() == WindingWindowShape::ROUND) {
        throw NotImplementedException("Pin assignment on a round winding window (a toroid on a family t base, RFC 0013 "
                                      "'Toroids on a base', WP4)");
    }
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
        if (!horizontal && rowIndex % 2 == 1) {
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

    // --- user pins: exist, and one owner each -----------------------------------------------
    std::vector<std::optional<std::pair<size_t, std::string>>> pinOwner(placedPins.size());  // plan index, what
    std::vector<IsolationSide> sidesInBuildOrder;
    for (const auto& plan : plans) {
        if (std::find(sidesInBuildOrder.begin(), sidesInBuildOrder.end(), plan.side) == sidesInBuildOrder.end()) {
            sidesInBuildOrder.push_back(plan.side);
        }
    }
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
    const size_t numberSides = sidesInBuildOrder.size();
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

    // Required creepage between isolation sides.
    if (numberSides < 2) {
        result.creepageNotCheckedReason = "a single isolation side needs no creepage between pins";
    }
    else if (!_inputs) {
        result.creepageNotCheckedReason = "the coil carries no inputs, so the required creepage is unknown";
    }
    else if (!_inputs->get_design_requirements().get_insulation()) {
        result.creepageNotCheckedReason = "the design requirements declare no insulation, so no creepage is required";
    }
    else {
        auto inputs = _inputs.value();
        result.requiredCreepage = InsulationCoordinator().calculate_creepage_distance(inputs, false);
    }
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("assign_pins needs the processed core for the through-core creepage path");
    }
    const double halfWidth = core.get_width() / 2;
    const double halfHeight = core.get_height() / 2;
    const double halfDepth = core.get_depth() / 2;
    auto creepage_between = [&](size_t pinA, size_t pinB) {
        PinCreepagePath path;
        path.pinA = placedPins[pinA].name;
        path.pinB = placedPins[pinB].name;
        path.surfaceDistance = distance(placedPins[pinA].base, placedPins[pinB].base);
        path.throughCoreDistance = distance_to_core_box(placedPins[pinA].base, halfWidth, halfHeight, halfDepth) +
                                   distance_to_core_box(placedPins[pinB].base, halfWidth, halfHeight, halfDepth);
        return path;
    };

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

    // --- creepage between every pair of sides ----------------------------------------------
    for (size_t pinA = 0; pinA < placedPins.size(); ++pinA) {
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

}  // namespace OpenMagnetics
