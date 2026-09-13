// ABT #1172 (WP3): the design-for-manufacturing rules that read the pin assignment - R2 (start
// and finish rail convention, inner pins first, no crossings), R3 (gauge spread on one rail),
// R4 (wire wrap against the standoff) and R14 (shield ends on pins). They live apart from
// Manufacturability.cpp because they depend on Coil::assign_pins and its pin placement.
// ABT #1173 (WP4): R16, the toroid rule (coating, barrier, TIW/FIW for reinforced, base limits),
// which reads the toroid base the pins stand on.

#include "advisers/Manufacturability.h"

#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Core.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <sstream>

namespace OpenMagnetics {

namespace {

std::string join_texts(const std::vector<std::string>& items, const std::string& separator) {
    std::string result;
    for (size_t index = 0; index < items.size(); ++index) {
        result += (index > 0 ? separator : std::string()) + items[index];
    }
    return result;
}

std::string format_number(double value, int precision) {
    std::ostringstream stream;
    stream.setf(std::ios::fixed);
    stream.precision(precision);
    stream << value;
    return stream.str();
}

// American Wire Gauge of a conductor of this diameter, by the gauge's definition (ASTM B258):
// d(n) = 0.127 mm * 92^((36 - n) / 39). Continuous, so a spread can be measured in gauges.
double awg_from_diameter(double diameter) {
    return 36.0 - 39.0 * std::log(diameter / 0.000127) / std::log(92.0);
}

// The processed pins of the coil's bobbin, or nullopt when it carries none.
std::optional<std::vector<MAS::Pin>> bobbin_pins(Coil& coil) {
    auto bobbin = coil.resolve_bobbin();
    if (!bobbin.get_processed_description() || !bobbin.get_processed_description()->get_pins() ||
        bobbin.get_processed_description()->get_pins()->empty()) {
        return std::nullopt;
    }
    return bobbin.get_processed_description()->get_pins().value();
}

// Windings in build order: the index of the first section each appears in.
std::vector<size_t> windings_in_build_order(Coil& coil) {
    auto windings = coil.get_functional_description();
    std::vector<std::pair<size_t, size_t>> rankAndIndex;
    auto sections = coil.get_sections_description();
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        size_t rank = std::numeric_limits<size_t>::max();
        if (sections) {
            for (size_t sectionIndex = 0; sectionIndex < sections->size(); ++sectionIndex) {
                for (const auto& partialWinding : sections.value()[sectionIndex].get_partial_windings()) {
                    if (partialWinding.get_winding() == windings[windingIndex].get_name()) {
                        rank = std::min(rank, sectionIndex);
                    }
                }
            }
        }
        rankAndIndex.push_back({rank, windingIndex});
    }
    std::stable_sort(rankAndIndex.begin(), rankAndIndex.end());
    std::vector<size_t> order;
    for (auto& [rank, index] : rankAndIndex) {
        order.push_back(index);
    }
    return order;
}

}  // namespace

// ---------------------------------------------------------------------------------------
// R2 - starts at one rail end, finishes at the other, inner pins first, no crossings
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r2_pin_order(Magnetic& magnetic) {
    auto finding = make_finding("R2");
    finding.set_unit("violations");
    finding.set_threshold_value(0);
    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so there are no winding ends");
        finding.set_message("Pin order not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& coil = magnetic.get_mutable_coil();
    auto pins = bobbin_pins(coil);
    if (!pins) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the bobbin has no pins[] (no pin geometry), so there is no rail to order ends along");
        finding.set_message("Pin order not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto placed = Coil::place_pins(pins.value());
    std::map<std::string, const PlacedPin*> pinByName;
    for (const auto& placedPin : placed) {
        pinByName[placedPin.name] = &placedPin;
    }
    auto windings = coil.get_functional_description();

    struct Ends {
        std::string winding;
        const PlacedPin* start = nullptr;
        const PlacedPin* finish = nullptr;
        double minimumX = std::numeric_limits<double>::max();
        double maximumX = std::numeric_limits<double>::lowest();
        std::set<size_t> rows;
    };
    std::vector<Ends> ends;
    for (auto windingIndex : windings_in_build_order(coil)) {
        const auto& winding = windings[windingIndex];
        if (!winding.get_connections()) {
            continue;
        }
        Ends windingEnds;
        windingEnds.winding = winding.get_name();
        const auto windingConnections = winding.get_connections().value();  // by value: the getter returns the optional by value
        for (const auto& connection : windingConnections) {
            if (!connection.get_pin_name() || !connection.get_end()) {
                continue;
            }
            auto found = pinByName.find(connection.get_pin_name().value());
            if (found == pinByName.end()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Winding '" + winding.get_name() + "' terminates on pin '" + connection.get_pin_name().value() +
                    "', which the bobbin does not have.");
            }
            const PlacedPin* pin = found->second;
            // The strand-0 (or shared) start and finish define the winding's direction.
            bool firstStrand = connection.get_parallel().value_or(0) == 0;
            if (connection.get_end().value() == End::START && firstStrand) {
                windingEnds.start = pin;
            }
            if (connection.get_end().value() == End::FINISH && firstStrand) {
                windingEnds.finish = pin;
            }
            windingEnds.minimumX = std::min(windingEnds.minimumX, pin->centre[0]);
            windingEnds.maximumX = std::max(windingEnds.maximumX, pin->centre[0]);
            windingEnds.rows.insert(pin->row);
        }
        if (windingEnds.start && windingEnds.finish) {
            ends.push_back(windingEnds);
        }
    }
    if (ends.empty()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("no winding has both its start and its finish assigned to a pin (Coil::assign_pins)");
        finding.set_message("Pin order not evaluated: " + finding.get_reason().value());
        return finding;
    }

    std::vector<std::string> violations;
    const bool horizontal = placed.front().hangsAlongZ;
    // The direction start -> finish along the row, normalised so that a vertical former's diagonal
    // convention (row 1 runs the other way) reads the same as a horizontal former's.
    std::optional<double> referenceDirection;
    std::string referenceWinding;
    for (const auto& windingEnds : ends) {
        if (windingEnds.rows.size() > 1) {
            violations.push_back("'" + windingEnds.winding + "' has ends on more than one row");
            continue;
        }
        double run = windingEnds.finish->centre[0] - windingEnds.start->centre[0];
        if (std::abs(run) < 1e-9) {
            continue;  // start and finish on one pin: no direction to compare
        }
        double direction = (run > 0 ? 1.0 : -1.0) * ((!horizontal && windingEnds.start->row % 2 == 1) ? -1.0 : 1.0);
        if (!referenceDirection) {
            referenceDirection = direction;
            referenceWinding = windingEnds.winding;
        }
        else if (direction != referenceDirection.value()) {
            violations.push_back("'" + windingEnds.winding + "' starts at the end of its rail where '" + referenceWinding +
                                 "' finishes (" + (horizontal ? std::string("horizontal formers put every start at one rail end")
                                                              : std::string("vertical formers put the starts in one diagonal corner")) + ")");
        }
    }
    // Inner pins first, no crossings: on one row, a winding wound later must enclose or avoid
    // every earlier winding's pins; a partial overlap is a crossing, and a later winding inside
    // an earlier one took the inner pins out of order.
    for (size_t first = 0; first < ends.size(); ++first) {
        for (size_t later = first + 1; later < ends.size(); ++later) {
            const auto& a = ends[first];
            const auto& b = ends[later];
            if (a.rows != b.rows) {
                continue;
            }
            bool disjoint = b.minimumX > a.maximumX + 1e-9 || b.maximumX < a.minimumX - 1e-9;
            bool bEnclosesA = b.minimumX < a.minimumX - 1e-9 && b.maximumX > a.maximumX + 1e-9;
            bool aEnclosesB = a.minimumX < b.minimumX - 1e-9 && a.maximumX > b.maximumX + 1e-9;
            if (disjoint || bEnclosesA) {
                continue;
            }
            if (aEnclosesB) {
                violations.push_back("'" + b.winding + "' is wound after '" + a.winding + "' but takes pins inside it (inner pins go to the first-wound winding)");
            }
            else {
                violations.push_back("the leads of '" + a.winding + "' and '" + b.winding + "' cross: their pins interleave on one rail");
            }
        }
    }
    finding.set_measured_value(static_cast<double>(violations.size()));
    if (violations.empty()) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("Every start sits at one rail end and every finish at the other, the first-wound winding takes the inner pins and no leads cross.");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::WARNING);
    finding.set_message("Pin order against Wuerth's practice: " + join_texts(violations, "; ") + ".");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R3 - gauge spread on one rail
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r3_rail_gauge_spread(Magnetic& magnetic) {
    auto finding = make_finding("R3");
    finding.set_unit("AWG");
    const double maximumSpread = rule_number("R3", "maximumGaugeSpreadAwg");
    finding.set_threshold_value(maximumSpread);
    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so there is no wire on any rail");
        finding.set_message("Rail gauge spread not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& coil = magnetic.get_mutable_coil();
    auto pins = bobbin_pins(coil);
    if (!pins) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the bobbin has no pins[] (no pin geometry), so no winding is on a rail");
        finding.set_message("Rail gauge spread not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto placed = Coil::place_pins(pins.value());
    std::map<std::string, size_t> rowByPin;
    for (const auto& placedPin : placed) {
        rowByPin[placedPin.name] = placedPin.row;
    }
    struct RailWire {
        std::string winding;
        double awg;
        bool litz;
    };
    std::map<size_t, std::vector<RailWire>> wiresByRow;
    auto windings = coil.get_functional_description();
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        const auto& winding = windings[windingIndex];
        if (!winding.get_connections()) {
            continue;
        }
        std::set<size_t> rows;
        const auto windingConnections = winding.get_connections().value();  // by value: the getter returns the optional by value
        for (const auto& connection : windingConnections) {
            if (connection.get_pin_name() && rowByPin.count(connection.get_pin_name().value())) {
                rows.insert(rowByPin.at(connection.get_pin_name().value()));
            }
        }
        if (rows.empty()) {
            continue;
        }
        auto wire = coil.resolve_wire(windingIndex);
        double diameter;
        if (wire.get_type() == WireType::ROUND) {
            if (!wire.get_conducting_diameter()) {
                throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
                    "Winding '" + winding.get_name() + "' is round wire without a conducting diameter; its gauge is unknown.");
            }
            diameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
        }
        else {
            // Litz, rectangular and foil: the gauge of the round wire of the same copper area.
            diameter = std::sqrt(4.0 * wire.calculate_conducting_area() / std::numbers::pi);
        }
        for (auto row : rows) {
            wiresByRow[row].push_back({winding.get_name(), awg_from_diameter(diameter), wire.get_type() == WireType::LITZ});
        }
    }
    if (wiresByRow.empty()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("no winding end is assigned to a pin (Coil::assign_pins)");
        finding.set_message("Rail gauge spread not evaluated: " + finding.get_reason().value());
        return finding;
    }
    double worstSpread = 0;
    std::vector<std::string> problems;
    for (auto& [row, railWires] : wiresByRow) {
        double finest = std::numeric_limits<double>::lowest();
        double heaviest = std::numeric_limits<double>::max();
        std::string finestName;
        std::string heaviestName;
        for (const auto& railWire : railWires) {
            if (railWire.awg > finest) {
                finest = railWire.awg;
                finestName = railWire.winding;
            }
            if (railWire.awg < heaviest) {
                heaviest = railWire.awg;
                heaviestName = railWire.winding;
            }
        }
        double spread = finest - heaviest;
        worstSpread = std::max(worstSpread, spread);
        if (spread > maximumSpread + 1e-9) {
            problems.push_back("rail " + std::to_string(row) + " carries '" + heaviestName + "' (" + format_number(heaviest, 1) +
                               " AWG) and '" + finestName + "' (" + format_number(finest, 1) + " AWG), " + format_number(spread, 1) +
                               " gauges apart");
        }
        if (railWires.size() > 1) {
            for (const auto& railWire : railWires) {
                if (railWire.litz) {
                    problems.push_back("rail " + std::to_string(row) + " carries litz '" + railWire.winding + "' with other windings");
                }
            }
        }
    }
    finding.set_measured_value(worstSpread);
    if (problems.empty()) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("Every rail carries wires within " + format_number(maximumSpread, 0) + " gauges of each other (worst " +
                            format_number(worstSpread, 1) + ").");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::WARNING);
    finding.set_message("Heavy and fine wire share a rail, which needs two soldering operations: " + join_texts(problems, "; ") + ".");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R4 - wire wrap against the standoff
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r4_wrap_height(Magnetic& magnetic) {
    auto finding = make_finding("R4");
    finding.set_unit("m");
    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so there is no wire wrap");
        finding.set_message("Wrap height not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& coil = magnetic.get_mutable_coil();
    auto pins = bobbin_pins(coil);
    auto bobbin = coil.resolve_bobbin();
    std::optional<double> standoff;
    if (bobbin.get_functional_description() && bobbin.get_functional_description()->get_base()) {
        standoff = resolve_dimensional_values(bobbin.get_functional_description()->get_base()->get_standoff());
    }
    if (!pins || !standoff) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("bobbin has no pin length/standoff");
        finding.set_message("Wrap height not evaluated: " + finding.get_reason().value() +
                            (pins ? " (the bobbin has pins but states no standoff)" : " (the bobbin has no pins[])") + ".");
        return finding;
    }
    finding.set_threshold_value(standoff.value());
    const double wrapTurns = rule_number("R4", "wrapTurns");

    std::map<std::string, double> wrapHeightByPin;
    std::map<std::string, std::vector<std::string>> windingsByPin;
    auto windings = coil.get_functional_description();
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        const auto& winding = windings[windingIndex];
        if (!winding.get_connections()) {
            continue;
        }
        auto wire = coil.resolve_wire(windingIndex);
        double outerSize;
        if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
            if (!wire.get_outer_diameter()) {
                throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
                    "Winding '" + winding.get_name() + "' has no wire outer diameter, so its wrap height is unknown.");
            }
            outerSize = resolve_dimensional_values(wire.get_outer_diameter().value());
        }
        else {
            outerSize = std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height());
        }
        const auto windingConnections = winding.get_connections().value();  // by value: the getter returns the optional by value
        for (const auto& connection : windingConnections) {
            if (!connection.get_pin_name()) {
                continue;
            }
            // Wires wrapped on this pin by this connection: every strand when `parallel` is absent,
            // and a tap carries two ends (the finish of one section and the start of the next).
            double wires = connection.get_parallel() ? 1.0 : static_cast<double>(winding.get_number_parallels());
            if (connection.get_end() && connection.get_end().value() == End::TAP) {
                wires *= 2;
            }
            wrapHeightByPin[connection.get_pin_name().value()] += wires * wrapTurns * outerSize;
            windingsByPin[connection.get_pin_name().value()].push_back(winding.get_name());
        }
    }
    if (wrapHeightByPin.empty()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("no winding end is assigned to a pin (Coil::assign_pins)");
        finding.set_message("Wrap height not evaluated: " + finding.get_reason().value());
        return finding;
    }
    double worstHeight = 0;
    std::string worstPin;
    for (auto& [pinName, height] : wrapHeightByPin) {
        if (height > worstHeight) {
            worstHeight = height;
            worstPin = pinName;
        }
    }
    finding.set_measured_value(worstHeight);
    finding.set_scope("pin " + worstPin);
    if (worstHeight <= standoff.value() + 1e-12) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("The tallest wrap (pin " + worstPin + ", " + format_number(worstHeight * 1000, 3) + " mm) stays within the " +
                            format_number(standoff.value() * 1000, 3) + " mm standoff.");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::FAIL);
    finding.set_message("The wrap on pin " + worstPin + " (" + join_texts(windingsByPin.at(worstPin), ", ") + ") is " +
                        format_number(worstHeight * 1000, 3) + " mm tall against a " + format_number(standoff.value() * 1000, 3) +
                        " mm standoff: the part would sit on its terminations instead of its standoffs. Use thinner wire, split the strands over more pins, or a bobbin with more pins.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R14 - shield windings terminated on pins at both ends
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r14_shield_terminations(Magnetic& magnetic) {
    auto finding = make_finding("R14");
    finding.set_unit("buried shield ends");
    finding.set_threshold_value(0);
    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so there are no shields");
        finding.set_message("Shield terminations not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& coil = magnetic.get_mutable_coil();
    std::set<std::string> shieldWindings;
    if (coil.get_sections_description()) {
        const auto coilSections = coil.get_sections_description().value();  // by value: the getter returns the optional by value
        for (const auto& section : coilSections) {
            if (section.get_type() == ElectricalType::SHIELDING) {
                for (const auto& partialWinding : section.get_partial_windings()) {
                    shieldWindings.insert(partialWinding.get_winding());
                }
            }
        }
    }
    if (coil.get_layers_description()) {
        const auto coilLayers = coil.get_layers_description().value();  // by value: the getter returns the optional by value
        for (const auto& layer : coilLayers) {
            if (layer.get_type() == ElectricalType::SHIELDING) {
                for (const auto& partialWinding : layer.get_partial_windings()) {
                    shieldWindings.insert(partialWinding.get_winding());
                }
            }
        }
    }
    if (shieldWindings.empty()) {
        finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
        finding.set_message("The design has no shield winding (no section or layer of type shielding).");
        return finding;
    }
    std::vector<std::string> buried;
    std::vector<std::string> unassigned;
    for (const auto& winding : coil.get_functional_description()) {
        if (!shieldWindings.count(winding.get_name())) {
            continue;
        }
        bool startOnPin = false;
        bool finishOnPin = false;
        bool anyEnd = false;
        if (winding.get_connections()) {
            const auto windingConnections = winding.get_connections().value();  // by value: the getter returns the optional by value
            for (const auto& connection : windingConnections) {
                if (!connection.get_end()) {
                    continue;
                }
                anyEnd = true;
                if (connection.get_end().value() == End::START && connection.get_pin_name()) {
                    startOnPin = true;
                }
                if (connection.get_end().value() == End::FINISH && connection.get_pin_name()) {
                    finishOnPin = true;
                }
            }
        }
        if (!anyEnd) {
            unassigned.push_back(winding.get_name());
            continue;
        }
        if (!startOnPin) {
            buried.push_back("'" + winding.get_name() + "' start");
        }
        if (!finishOnPin) {
            buried.push_back("'" + winding.get_name() + "' finish");
        }
    }
    if (!unassigned.empty()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("shield winding(s) " + join_texts(unassigned, ", ") + " carry no connections with an end (Coil::assign_pins)");
        finding.set_message("Shield terminations not evaluated: " + finding.get_reason().value());
        return finding;
    }
    finding.set_measured_value(static_cast<double>(buried.size()));
    if (buried.empty()) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("Every shield winding terminates on pins at both ends, which the winding machine can automate.");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::WARNING);
    finding.set_message("Shield ends not on a pin: " + join_texts(buried, ", ") +
                        ". Terminating both ends on pins is automatable; a buried end is manual work.");
    return finding;
}


// ---------------------------------------------------------------------------------------
// R16 - toroids: coated core for basic, barrier between sides, TIW/FIW for double/reinforced,
//       the base's outer-diameter and height limits
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r16_toroid_insulation_and_header(Magnetic& magnetic, Inputs& inputs) {
    auto finding = make_finding("R16");
    if (!magnetic.has_core()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no core, so it is unknown whether it is a toroid");
        finding.set_message("Toroid rule not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& core = magnetic.get_mutable_core();
    if (core.get_shape_family() != CoreShapeFamily::T) {
        finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
        finding.set_message("The core is not a toroid.");
        return finding;
    }
    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so there is no winding to insulate or base to seat");
        finding.set_message("Toroid rule not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& coil = magnetic.get_mutable_coil();
    std::vector<std::string> failures;
    std::vector<std::string> passes;
    std::vector<std::string> notEvaluated;

    // --- insulation: coating, barrier, TIW/FIW ------------------------------------------------
    const auto insulation = inputs.get_design_requirements().get_insulation();
    if (!insulation || !insulation->get_insulation_type()) {
        notEvaluated.push_back("coating, barrier and wire grade: the design requirements state no insulation type");
    }
    else {
        const auto insulationType = insulation->get_insulation_type().value();
        if (insulationType == IsolationClass::FUNCTIONAL) {
            passes.push_back("functional insulation needs no coating, barrier or insulated wire");
        }
        else {
            // The coating counts only when the core DECLARES it: Core::get_coating_thickness models an
            // undeclared toroid as coated, which is a winding-geometry assumption, not insulation.
            if (!core.get_functional_description().get_coating()) {
                failures.push_back("the toroid core declares no coating, and the core coating is what provides the basic insulation between core and wire");
            }
            else {
                passes.push_back("the core declares its coating (basic insulation between core and wire)");
            }

            auto windings = coil.get_functional_description();
            std::set<IsolationSide> sides;
            for (const auto& winding : windings) {
                sides.insert(winding.get_isolation_side());
            }
            if (sides.size() > 1) {
                if (!coil.get_sections_description()) {
                    notEvaluated.push_back("barrier between isolation sides: the coil has not been wound, so it has no sections");
                }
                else {
                    // A mechanical separation between two sides: an insulation section between every pair of
                    // conduction sections of different sides, in the coil's section order.
                    auto sections = coil.get_sections_description().value();
                    std::map<std::string, IsolationSide> sideOfWinding;
                    for (const auto& winding : windings) {
                        sideOfWinding[winding.get_name()] = winding.get_isolation_side();
                    }
                    std::vector<std::string> unseparated;
                    std::optional<std::pair<IsolationSide, std::string>> previous;
                    bool barrierSincePrevious = false;
                    for (const auto& section : sections) {
                        if (section.get_type() == ElectricalType::INSULATION) {
                            barrierSincePrevious = true;
                            continue;
                        }
                        if (section.get_partial_windings().empty()) {
                            continue;
                        }
                        const auto side = sideOfWinding.at(section.get_partial_windings().front().get_winding());
                        if (previous && previous->first != side && !barrierSincePrevious) {
                            unseparated.push_back("'" + previous->second + "' and '" + section.get_name() + "'");
                        }
                        previous = std::make_pair(side, section.get_name());
                        barrierSincePrevious = false;
                    }
                    if (unseparated.empty()) {
                        passes.push_back("an insulation section separates every pair of sections of different isolation sides");
                    }
                    else {
                        failures.push_back("no mechanical separation (insulation section) between sections " + join_texts(unseparated, ", "));
                    }
                }
                if (insulationType == IsolationClass::DOUBLE || insulationType == IsolationClass::REINFORCED) {
                    // Double/reinforced on a toroid: one isolation side wound entirely in TIW or FIW.
                    std::map<IsolationSide, bool> sideAllInsulatedWire;
                    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
                        auto wire = coil.resolve_wire(windingIndex);
                        const bool insulatedWire = is_triple_insulated_wire(wire) || is_fully_insulated_wire(wire);
                        const auto side = windings[windingIndex].get_isolation_side();
                        sideAllInsulatedWire[side] = (sideAllInsulatedWire.count(side) ? sideAllInsulatedWire.at(side) : true) && insulatedWire;
                    }
                    bool anySide = false;
                    for (const auto& [side, allInsulated] : sideAllInsulatedWire) {
                        anySide = anySide || allInsulated;
                    }
                    if (anySide) {
                        passes.push_back("one isolation side is wound entirely in triple-insulated or fully insulated wire");
                    }
                    else {
                        failures.push_back("double or reinforced insulation on a toroid needs one isolation side wound in triple-insulated (TIW) or fully insulated (FIW) wire, and no side is");
                    }
                }
            }
        }
    }

    // --- header: the base's limits against the coated ring ----------------------------------------
    auto bobbin = coil.resolve_bobbin();
    if (!bobbin.has_base()) {
        notEvaluated.push_back("header limits: the coil's bobbin has no toroid base");
    }
    else {
        const auto base = bobbin.get_base();
        const auto shapeDimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
        const double coating = core.get_coating_thickness();
        const double ringOuterDiameter = shapeDimensions.at("A") + 2 * coating;
        const double ringHeight = shapeDimensions.at("C") * static_cast<double>(core.get_number_stacks()) + 2 * coating;
        finding.set_unit("m");
        finding.set_measured_value(ringOuterDiameter);
        if (base.get_maximum_core_outer_diameter()) {
            finding.set_threshold_value(base.get_maximum_core_outer_diameter().value());
            if (ringOuterDiameter > base.get_maximum_core_outer_diameter().value()) {
                failures.push_back("the coated ring's outer diameter " + format_number(ringOuterDiameter * 1000, 2) + " mm exceeds the base's maximum " +
                                   format_number(base.get_maximum_core_outer_diameter().value() * 1000, 2) + " mm");
            }
            else {
                passes.push_back("the coated ring's outer diameter fits the base");
            }
        }
        else {
            notEvaluated.push_back("header outer-diameter limit: the base states no maximumCoreOuterDiameter");
        }
        if (base.get_maximum_core_height()) {
            if (ringHeight > base.get_maximum_core_height().value()) {
                failures.push_back("the coated ring's height " + format_number(ringHeight * 1000, 2) + " mm exceeds the base's maximum " +
                                   format_number(base.get_maximum_core_height().value() * 1000, 2) + " mm");
            }
            else {
                passes.push_back("the coated ring's height fits the base");
            }
        }
        else {
            notEvaluated.push_back("header height limit: the base states no maximumCoreHeight");
        }
        notEvaluated.push_back("header creepage and cavity depth: MAS bobbin.json base has no field for them");
    }

    std::string tail = notEvaluated.empty() ? std::string() : " Not evaluated: " + join_texts(notEvaluated, "; ") + ".";
    if (!failures.empty()) {
        finding.set_status(ManufacturabilityStatus::FAIL);
        finding.set_message("Toroid: " + join_texts(failures, "; ") + "." + tail);
        return finding;
    }
    if (passes.empty()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason(join_texts(notEvaluated, "; "));
        finding.set_message("Toroid rule not evaluated: " + finding.get_reason().value() + ".");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::PASS);
    finding.set_message("Toroid: " + join_texts(passes, "; ") + "." + tail);
    return finding;
}

}  // namespace OpenMagnetics
