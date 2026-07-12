#include <cmath>
#include "constructive_models/MasMigration.h"
#include <algorithm>
#include <map>
#include <set>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <cfloat>
#include <cmath>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Utils.h"
#include "constructive_models/Coil.h"
#include "json.hpp"
#include "constructive_models/InsulationMaterial.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Exceptions.h"
#include "support/Logger.h"

using json = nlohmann::json;

namespace OpenMagnetics {

// FIX H-COIL-2: Helper to safely compute equal proportion per winding
static inline std::vector<double> make_equal_proportion_per_winding(size_t numWindings) {
    if (numWindings == 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Cannot compute proportion for zero windings");
    }
    return std::vector<double>(numWindings, 1.0 / numWindings);
}

// FIX L-COIL-3: Integer factorial to replace fragile tgamma() usage
static inline size_t factorial(size_t n) {
    size_t result = 1;
    for (size_t i = 2; i <= n; ++i) result *= i;
    return result;
}



std::vector<double> Coil::cartesian_to_polar(std::vector<double> value) {
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        throw InvalidInputException("Not supposed to convert for these cores");
    }
    else {
        std::vector<double> convertedValue;
        double angle = atan2(value[1], value[0]) * 180 / std::numbers::pi;
        if (angle < 0) {
            angle += 360;
        }
        double radius = hypot(value[0], value[1]);
        double radialHeight = windingWindows[0].get_radial_height().value() - radius;
        return {radialHeight, angle};
    }
}

std::vector<double> Coil::cartesian_to_polar(std::vector<double> value, double radialHeight) {
    std::vector<double> convertedValue;
    double angle = atan2(value[1], value[0]) * 180 / std::numbers::pi;
    if (angle < 0) {
        angle += 360;
    }
    double radius = hypot(value[0], value[1]);
    double turnRadialHeight = radialHeight - radius;
    return {turnRadialHeight, angle};
}

std::vector<double> Coil::polar_to_cartesian(std::vector<double> value) {
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        throw InvalidInputException("Not supposed to convert for these cores");
    }
    else {
        std::vector<double> convertedValue;
        double radius = windingWindows[0].get_radial_height().value() - value[0];
        double angleRadians = value[1] / 180 * std::numbers::pi;
        double x = radius * cos(angleRadians);
        double y = radius * sin(angleRadians);
        return {x, y};
    }
}

std::vector<double> Coil::polar_to_cartesian(std::vector<double> value, double radialHeight) {
    std::vector<double> convertedValue;
    double radius = radialHeight - value[0];
    double angleRadians = value[1] / 180 * std::numbers::pi;
    double x = radius * cos(angleRadians);
    double y = radius * sin(angleRadians);
    return {x, y};
}

void Coil::convert_turns_to_cartesian_coordinates() {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    if (!get_turns_description()) {
        throw CoilNotProcessedException("Missing turns");
    }

    auto turns = get_turns_description().value();
    if (turns[0].get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        return;
    }

    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto cartesianCoordinates = polar_to_cartesian(turns[turnIndex].get_coordinates(), windingWindowRadialHeight);
        turns[turnIndex].set_coordinates(cartesianCoordinates);
        turns[turnIndex].set_coordinate_system(CoordinateSystem::CARTESIAN);
        if (turns[turnIndex].get_additional_coordinates()) {
            auto additionalCoordinates = turns[turnIndex].get_additional_coordinates().value();
            for (size_t additionalTurnIndex = 0; additionalTurnIndex < additionalCoordinates.size(); ++additionalTurnIndex) {
                additionalCoordinates[additionalTurnIndex] = polar_to_cartesian(additionalCoordinates[additionalTurnIndex], windingWindowRadialHeight);
            }
            turns[turnIndex].set_additional_coordinates(additionalCoordinates);
        }
    }

    set_turns_description(turns);
}

void Coil::convert_turns_to_polar_coordinates() {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    if (!get_turns_description()) {
        throw CoilNotProcessedException("Missing turns");
    }

    auto turns = get_turns_description().value();
    if (turns[0].get_coordinate_system().value() == CoordinateSystem::POLAR) {
        return;
    }

    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto cartesianCoordinates = cartesian_to_polar(turns[turnIndex].get_coordinates(), windingWindowRadialHeight);
        turns[turnIndex].set_coordinates(cartesianCoordinates);
        turns[turnIndex].set_coordinate_system(CoordinateSystem::POLAR);
        if (turns[turnIndex].get_additional_coordinates()) {
            auto additionalCoordinates = turns[turnIndex].get_additional_coordinates().value();
            for (size_t additionalTurnIndex = 0; additionalTurnIndex < additionalCoordinates.size(); ++additionalTurnIndex) {
                additionalCoordinates[additionalTurnIndex] = cartesian_to_polar(additionalCoordinates[additionalTurnIndex], windingWindowRadialHeight);
            }
            turns[turnIndex].set_additional_coordinates(additionalCoordinates);
        }
    }

    set_turns_description(turns);
}

bool Coil::is_planar() {
    if (resolve_wire(0).get_type() == WireType::PLANAR) {
        return true;
    }
    return false;
}

Coil::Coil(json j, size_t interleavingLevel,
                               WindingOrientation windingOrientation,
                               WindingOrientation layersOrientation,
                               CoilAlignment turnsAlignment,
                               CoilAlignment sectionAlignment) {
        OpenMagnetics::compat::migrate_pre_1_0(j);
    _interleavingLevel = interleavingLevel;
    _windingOrientation = windingOrientation;
    _layersOrientation = layersOrientation;
    _turnsAlignment = turnsAlignment;
    _sectionAlignment = sectionAlignment;
    _sectionAlignmentExplicit = true;
    from_json(j, *this);

    if (!is_planar()) {
        wind();
    }
}

Coil::Coil(json j, bool windInConstructor) {
        OpenMagnetics::compat::migrate_pre_1_0(j);
    from_json(j, *this);

    if (windInConstructor) {
        if (!is_planar()) {
            wind();
        }
    }
}

Coil::Coil(const MAS::Coil coil) {
    bool hasSectionsData = false;
    bool hasLayersData = false;
    bool hasTurnsData = false;

    set_functional_description({});
    for (auto winding : coil.get_functional_description()) {
        get_mutable_functional_description().push_back(winding);
    }

    auto bobbinVariant = coil.get_bobbin();
    if (std::holds_alternative<std::string>(bobbinVariant)) {
        auto bobbinName = std::get<std::string>(bobbinVariant);
        set_bobbin(bobbinName);
    }
    else {
        auto bobbin = OpenMagnetics::Bobbin(bobbinVariant);
        set_bobbin(bobbin);
    }

    if (coil.get_sections_description()) {
        hasSectionsData = true;
        set_sections_description(coil.get_sections_description());
    }
    if (coil.get_layers_description()) {
        hasLayersData = true;
        set_layers_description(coil.get_layers_description());
    }
    if (coil.get_turns_description()) {
        hasTurnsData = true;
        set_turns_description(coil.get_turns_description());
    }
    auto delimitAndCompact = settings.get_coil_delimit_and_compact();

    if (!hasSectionsData || !hasLayersData || (!hasTurnsData && are_sections_and_layers_fitting())) {
        if (wind() && delimitAndCompact) {
            delimit_and_compact();
        }
    }

}

void Coil::log(std::string entry) {
    coilLog += entry + "\n";
}

std::string Coil::read_log() {
    return coilLog;
}

void Coil::set_strict(bool value) {
    _strict = value;
}

void Coil::set_inputs(Inputs inputs) {
    _inputs = inputs;
}

void Coil::set_interleaving_level(uint8_t interleavingLevel) {
    _interleavingLevel = interleavingLevel;
    _marginsPerSection = std::vector<std::vector<double>>(interleavingLevel, {0, 0});
}

void Coil::reset_margins_per_section() {
    _marginsPerSection.clear();
}

void Coil::reset_insulation() {
    _marginsPerSection.clear();
    _insulationSections.clear();
}

size_t Coil::get_interleaving_level() const {
    return _currentRepetitions;
}

void Coil::set_winding_orientation(WindingOrientation windingOrientation) {
    _windingOrientation = windingOrientation;
    auto bobbin = resolve_bobbin();
    if (bobbin.get_processed_description()) {
        bobbin.set_winding_orientation(windingOrientation);
        set_bobbin(bobbin);
    }
}

void Coil::set_layers_orientation(WindingOrientation layersOrientation, std::optional<std::string> sectionName) {
    if (sectionName) {
        _layersOrientationPerSection[sectionName.value()] = layersOrientation;
    }
    else {
        _layersOrientation = layersOrientation;
    }
}

void Coil::set_turns_alignment(CoilAlignment turnsAlignment, std::optional<std::string> sectionName) {
    if (sectionName) {
        _turnsAlignmentPerSection[sectionName.value()] = turnsAlignment;
    }
    else {
        _turnsAlignment = turnsAlignment;
    }
}

void Coil::set_section_alignment(CoilAlignment sectionAlignment) {
    _sectionAlignment = sectionAlignment;
    _sectionAlignmentExplicit = true;
}

WindingOrientation Coil::get_winding_orientation() {
    auto bobbin = resolve_bobbin();
    auto windingOrientationFromBobbin = bobbin.get_winding_orientation();

    if (!windingOrientationFromBobbin) {
        return _windingOrientation;
    }
    else {
        return windingOrientationFromBobbin.value();
    }
}

WindingOrientation Coil::get_layers_orientation() const {
    return _layersOrientation;
}

CoilAlignment Coil::get_turns_alignment(std::optional<std::string> sectionName) const {
    if (sectionName) {
        if (_turnsAlignmentPerSection.find(sectionName.value()) != _turnsAlignmentPerSection.end()) {
            return _turnsAlignmentPerSection.at(sectionName.value());
        }
        else {
            return _turnsAlignment;
        }
    }
    else {
        return _turnsAlignment;
    }
}

WindingOrder Coil::get_winding_order(const std::string& sectionName) const {
    // Per-section override wins.
    if (get_sections_description()) {
        auto sections = get_sections_description().value();
        for (const auto& section : sections) {
            if (section.get_name() == sectionName) {
                if (section.get_winding_order()) {
                    return section.get_winding_order().value();
                }
                break;
            }
        }
    }

    // Else the bobbin winding window's default.
    if (std::holds_alternative<Bobbin>(get_bobbin())) {
        auto bobbin = std::get<Bobbin>(get_bobbin());
        if (bobbin.get_processed_description()) {
            auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
            if (windingWindows.size() > 0 && windingWindows[0].get_winding_order()) {
                return windingWindows[0].get_winding_order().value();
            }
        }
    }

    // Default preserves the historical behaviour (every layer wound the same direction).
    return WindingOrder::Z;
}

// Connection leads for a TOROIDAL winding. Toroidal turns are stored in cartesian coordinates on
// circles around the core centre; layers are concentric polar rings. Each parallel is its own
// conductor: its entrance/exit terminal leads run radially out to the winding-window border, and its
// inter-layer continuations run straight (cartesian) from the last turn of one ring to the first turn
// of the next. Drawn (and fed to the connection loss) but no turn-displacement blocking yet — toroidal
// blocking is angular and is a separate follow-up.
static std::vector<ConnectionReservedSpace> toroidal_connection_reserved_spaces(Coil& coil) {
    std::vector<ConnectionReservedSpace> spaces;
    if (!coil.get_turns_description() || !coil.get_layers_description()) {
        return spaces;
    }
    auto turns = coil.get_turns_description().value();
    auto wires = coil.get_wires();
    auto allLayersDescription = coil.get_layers_description().value();

    std::vector<Layer> allLayers;
    for (const auto& layer : allLayersDescription) {
        if (layer.get_type() == ElectricalType::CONDUCTION && !layer.get_partial_windings().empty()) {
            allLayers.push_back(layer);
        }
    }

    // Electrical order of the rings = order their turns are first wound.
    std::map<std::string, size_t> layerElectricalOrder;
    size_t order = 0;
    for (const auto& turn : turns) {
        if (turn.get_layer() && layerElectricalOrder.find(turn.get_layer().value()) == layerElectricalOrder.end()) {
            layerElectricalOrder[turn.get_layer().value()] = order++;
        }
    }

    // Entrance/exit turn per (winding, parallel); first/last turn per (layer, parallel).
    std::map<std::pair<std::string, int64_t>, Turn> entranceTurn, exitTurn, firstTurnInLayer, lastTurnInLayer;
    double maxTurnRadius = 0;
    for (const auto& turn : turns) {
        auto windingKey = std::make_pair(turn.get_winding(), turn.get_parallel());
        if (!entranceTurn.count(windingKey)) {
            entranceTurn[windingKey] = turn;
        }
        exitTurn[windingKey] = turn;
        if (turn.get_layer()) {
            auto layerKey = std::make_pair(turn.get_layer().value(), turn.get_parallel());
            if (!firstTurnInLayer.count(layerKey)) {
                firstTurnInLayer[layerKey] = turn;
            }
            lastTurnInLayer[layerKey] = turn;
        }
        auto c = turn.get_coordinates();
        maxTurnRadius = std::max(maxTurnRadius, std::hypot(c[0], c[1]));
    }

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto windingName = coil.get_functional_description()[windingIndex].get_name();
        double wireOuterWidth = wires[windingIndex].get_maximum_outer_width();
        double wireOuterHeight = wires[windingIndex].get_maximum_outer_height();
        int64_t numberParallels = int64_t(coil.get_number_parallels(windingIndex));
        // Terminal leads route radially out past the outermost turn to the window border.
        double radialBorder = maxTurnRadius + 1.5 * wireOuterWidth;

        auto addTerminalLead = [&](const Turn& connectingTurn, int64_t parallel) {
            auto c = connectingTurn.get_coordinates();
            double radius = std::hypot(c[0], c[1]);
            if (radius <= 1e-9 || radialBorder <= radius) {
                return;
            }
            double angle = std::atan2(c[1], c[0]);
            double radiusMid = (radius + radialBorder) / 2;
            ConnectionReservedSpace lead;
            lead.isTerminal = true;
            lead.winding = windingName;
            lead.parallel = parallel;
            lead.section = connectingTurn.get_section().value_or("");
            lead.layer = "";
            lead.coordinates = {roundFloat(radiusMid * std::cos(angle), 9), roundFloat(radiusMid * std::sin(angle), 9)};
            lead.dimensions = {roundFloat(radialBorder - radius, 9), wireOuterHeight};
            lead.rotation = roundFloat(angle * 180.0 / std::numbers::pi, 6);
            spaces.push_back(lead);
        };
        for (int64_t parallel = 0; parallel < numberParallels; ++parallel) {
            auto key = std::make_pair(windingName, parallel);
            if (entranceTurn.count(key)) {
                addTerminalLead(entranceTurn.at(key), parallel);
            }
            if (exitTurn.count(key)) {
                addTerminalLead(exitTurn.at(key), parallel);
            }
        }

        std::vector<Layer> windingLayers;
        for (const auto& layer : allLayers) {
            if (layer.get_partial_windings()[0].get_winding() == windingName) {
                windingLayers.push_back(layer);
            }
        }
        if (windingLayers.size() < 2) {
            continue;
        }
        std::sort(windingLayers.begin(), windingLayers.end(), [&](const Layer& a, const Layer& b) {
            size_t orderA = layerElectricalOrder.count(a.get_name()) ? layerElectricalOrder.at(a.get_name()) : 0;
            size_t orderB = layerElectricalOrder.count(b.get_name()) ? layerElectricalOrder.at(b.get_name()) : 0;
            return orderA < orderB;
        });

        for (int64_t parallel = 0; parallel < numberParallels; ++parallel) {
            for (size_t i = 0; i + 1 < windingLayers.size(); ++i) {
                auto exitKey = std::make_pair(windingLayers[i].get_name(), parallel);
                auto entryKey = std::make_pair(windingLayers[i + 1].get_name(), parallel);
                if (!lastTurnInLayer.count(exitKey) || !firstTurnInLayer.count(entryKey)) {
                    continue;
                }
                auto a = lastTurnInLayer.at(exitKey).get_coordinates();
                auto b = firstTurnInLayer.at(entryKey).get_coordinates();
                double deltaX = b[0] - a[0];
                double deltaY = b[1] - a[1];
                double length = std::hypot(deltaX, deltaY);
                if (length <= 1e-9) {
                    continue;
                }
                ConnectionReservedSpace link;
                link.winding = windingName;
                link.parallel = parallel;
                link.section = windingLayers[i].get_section().value_or("");
                link.layer = "";
                link.coordinates = {roundFloat((a[0] + b[0]) / 2, 9), roundFloat((a[1] + b[1]) / 2, 9)};
                link.dimensions = {roundFloat(length, 9), wireOuterHeight};
                link.rotation = roundFloat(std::atan2(deltaY, deltaX) * 180.0 / std::numbers::pi, 6);
                spaces.push_back(link);
            }
        }
    }
    return spaces;
}

std::vector<ConnectionReservedSpace> Coil::get_connection_reserved_spaces() {
    // Model (first approximation, to be validated): each winding's wire routes through its
    // conduction layers in electrical (wound) order. The lead between two electrically-consecutive
    // layers of the same winding reserves one wire-thick rectangle on every conduction layer it
    // passes radially. For an interleaved winding (its halves separated by another winding's layer)
    // that continuation lead crosses, and squeezes, the intervening layer; for adjacent layers it
    // reserves at the single boundary it steps over. Rectangular (concentric), overlapping layers
    // only; contiguous and toroidal windows are a TODO.
    std::vector<ConnectionReservedSpace> spaces;
    if (!get_layers_description() || !get_sections_description() || !get_turns_description()) {
        return spaces;
    }
    auto bobbin = resolve_bobbin();
    if (bobbin.get_winding_window_shape() == WindingWindowShape::ROUND) {
        return toroidal_connection_reserved_spaces(*this);
    }
    if (bobbin.get_winding_window_shape() != WindingWindowShape::RECTANGULAR) {
        return spaces;
    }
    auto wires = get_wires();
    auto allLayersDescription = get_layers_description().value();
    auto turns = get_turns_description().value();

    // Rectangular winding window. OVERLAPPING layers stack radially (along x, the "layer axis") with
    // turns stacking axially (along y); CONTIGUOUS layers are the 90° transpose — they stack along y
    // with turns along x. We support both by running the model in a virtual frame where the layer axis
    // is always x: for the contiguous case we transpose (swap x↔y of) the turns, layers, window and
    // wire on the way in and transpose the produced rectangles back on the way out, reusing the tested
    // overlapping logic verbatim.
    std::vector<Layer> allLayers;
    for (const auto& layer : allLayersDescription) {
        if (layer.get_type() == ElectricalType::CONDUCTION
            && (layer.get_orientation() == WindingOrientation::OVERLAPPING
                || layer.get_orientation() == WindingOrientation::CONTIGUOUS)) {
            allLayers.push_back(layer);
        }
    }
    // Single-layer windings still have entrance/exit TERMINAL leads (drawn, and feeding
    // the connection loss); only the inter-layer links need >= 2 layers, and that loop
    // no-ops naturally. Bail out only when there is no conduction layer at all.
    if (allLayers.empty()) {
        return spaces;
    }
    bool layersAreContiguous = (allLayers[0].get_orientation() == WindingOrientation::CONTIGUOUS);
    if (layersAreContiguous) {
        for (auto& turn : turns) {
            auto coordinates = turn.get_coordinates();
            std::swap(coordinates[0], coordinates[1]);
            turn.set_coordinates(coordinates);
        }
        for (auto& layer : allLayers) {
            auto coordinates = layer.get_coordinates();
            std::swap(coordinates[0], coordinates[1]);
            layer.set_coordinates(coordinates);
        }
    }
    std::sort(allLayers.begin(), allLayers.end(), [](const Layer& a, const Layer& b) {
        return a.get_coordinates()[0] < b.get_coordinates()[0];
    });

    // Electrical order of layers = the order in which their turns are first wound.
    std::map<std::string, size_t> layerElectricalOrder;
    size_t order = 0;
    for (const auto& turn : turns) {
        if (turn.get_layer() && layerElectricalOrder.find(turn.get_layer().value()) == layerElectricalOrder.end()) {
            layerElectricalOrder[turn.get_layer().value()] = order++;
        }
    }

    // Outer (radial) border of the winding window: where terminal leads route to their terminals. In
    // the virtual frame the layer axis is x and the turn axis is y, so for contiguous layers we swap
    // the window centre and width/height to match the transposed turns/layers above.
    auto windingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
    double windowCenterLayerAxis = windingWindow.get_coordinates().value()[0];
    double windowCenterTurnAxis = windingWindow.get_coordinates().value()[1];
    double windowSizeLayerAxis = windingWindow.get_width().value();
    double windowSizeTurnAxis = windingWindow.get_height().value();
    if (layersAreContiguous) {
        std::swap(windowCenterLayerAxis, windowCenterTurnAxis);
        std::swap(windowSizeLayerAxis, windowSizeTurnAxis);
    }
    double windowOuterX = windowCenterLayerAxis + windowSizeLayerAxis / 2;
    // Axial extent of the window: terminal leads run along its top edge (entrance) or bottom edge
    // (exit), i.e. above/below all the (blocking-shrunk) conduction layers.
    double windowCenterY = windowCenterTurnAxis;
    double windowTopY = windowCenterY + windowSizeTurnAxis / 2;
    double windowBottomY = windowCenterY - windowSizeTurnAxis / 2;

    // First-wound (entrance) and last-wound (exit) turn of each (winding, parallel): each parallel of
    // a bifilar/N-filar group is its own physical conductor and gets its own entrance/exit terminal
    // leads. Likewise the first/last turn of each (layer, parallel), used to tell a U turnaround (the
    // connecting turns sit at the same axial end) from a Z dragback (opposite ends) — per parallel.
    std::map<std::pair<std::string, int64_t>, Turn> entranceTurnByWindingParallel;
    std::map<std::pair<std::string, int64_t>, Turn> exitTurnByWindingParallel;
    std::map<std::pair<std::string, int64_t>, Turn> firstTurnByLayerParallel;
    std::map<std::pair<std::string, int64_t>, Turn> lastTurnByLayerParallel;
    for (const auto& turn : turns) {
        auto windingKey = std::make_pair(turn.get_winding(), turn.get_parallel());
        if (entranceTurnByWindingParallel.find(windingKey) == entranceTurnByWindingParallel.end()) {
            entranceTurnByWindingParallel[windingKey] = turn;
        }
        exitTurnByWindingParallel[windingKey] = turn;
        if (turn.get_layer()) {
            auto layerKey = std::make_pair(turn.get_layer().value(), turn.get_parallel());
            if (firstTurnByLayerParallel.find(layerKey) == firstTurnByLayerParallel.end()) {
                firstTurnByLayerParallel[layerKey] = turn;
            }
            lastTurnByLayerParallel[layerKey] = turn;
        }
    }

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        auto windingName = get_functional_description()[windingIndex].get_name();
        // Wire footprint in the virtual frame: width is along the layer axis, height along the turn
        // axis, so swap them for contiguous layers (where the wire runs along x within a layer).
        double wireOuterWidth = wires[windingIndex].get_maximum_outer_width();
        double wireOuterHeight = wires[windingIndex].get_maximum_outer_height();
        if (layersAreContiguous) {
            std::swap(wireOuterWidth, wireOuterHeight);
        }
        int64_t numberParallels = int64_t(get_number_parallels(windingIndex));

        // Entrance and exit terminal leads route a winding end radially out to the outer window border.
        // Each PARALLEL of a bifilar/N-filar group is its own conductor, so each emits its own entrance
        // and exit lead. They run ALONG the nearest window edge — so each lead occupies the extreme
        // (top-most / bottom-most) turn slot of every layer it crosses, which is exactly the slot
        // turn-blocking frees. The lead therefore: (1) squeezes each crossed layer at that layer's
        // extreme edge, (2) is drawn as an L — a short vertical stub from the connecting turn up/down to
        // its own layer's edge, then a horizontal run along the edge to the border. Routing at the
        // connecting turn's interior level instead would clip the wound turns of the layers it crosses.
        auto addTerminalLead = [&](const Turn& connectingTurn, int64_t parallel) {
            double turnX = connectingTurn.get_coordinates()[0];
            double turnY = connectingTurn.get_coordinates()[1];
            if (windowOuterX <= turnX) {
                return;
            }
            // The layers the lead routes OVER (outward of the connecting turn) to reach the border.
            std::vector<const Layer*> crossedLayers;
            for (const auto& crossed : allLayers) {
                double crossedX = crossed.get_coordinates()[0];
                if (crossedX > turnX + 1e-9 && crossedX < windowOuterX) {
                    crossedLayers.push_back(&crossed);
                }
            }

            if (crossedLayers.empty()) {
                // Outermost end: the lead crosses nothing, so it just leaves radially at its OWN axial
                // level straight to the border — no edge routing, no vertical stub (which would cross its
                // own column when the end sits mid-window, e.g. a half-full outer interleaved section).
                ConnectionReservedSpace lead;
                lead.isTerminal = true;
                lead.winding = windingName;
                lead.parallel = parallel;
                lead.section = connectingTurn.get_section().value_or("");
                lead.layer = "";
                lead.coordinates = {roundFloat((turnX + windowOuterX) / 2, 9), roundFloat(turnY, 9)};
                lead.dimensions = {roundFloat(windowOuterX - turnX + wireOuterWidth, 9), wireOuterHeight};
                spaces.push_back(lead);
                return;
            }

            // Crosses outer layers: route along the window edge NEAREST the end so the lead sits in the
            // extreme slots that turn-blocking frees on the crossed layers. A short stub bridges the end
            // turn up/down to that edge within its own column.
            bool turnAtTop = (turnY >= windowCenterY);
            double edgeY = turnAtTop ? roundFloat(windowTopY - wireOuterHeight / 2, 9)
                                     : roundFloat(windowBottomY + wireOuterHeight / 2, 9);
            for (const Layer* crossed : crossedLayers) {
                ConnectionReservedSpace space;
                space.isTerminal = true;
                space.winding = windingName;
                space.parallel = parallel;
                space.section = crossed->get_section().value_or("");
                space.layer = crossed->get_name();
                space.coordinates = {crossed->get_coordinates()[0], edgeY};
                space.dimensions = {wireOuterWidth, wireOuterHeight};
                spaces.push_back(space);
            }
            if (std::abs(edgeY - turnY) > wireOuterHeight / 2) {
                double stubDirection = (edgeY >= turnY) ? 1.0 : -1.0;
                double stubFarEnd = edgeY + stubDirection * wireOuterHeight / 2;
                ConnectionReservedSpace stub;
                stub.isTerminal = true;
                stub.winding = windingName;
                stub.parallel = parallel;
                stub.section = connectingTurn.get_section().value_or("");
                stub.layer = "";
                stub.coordinates = {turnX, roundFloat((turnY + stubFarEnd) / 2, 9)};
                stub.dimensions = {wireOuterWidth, roundFloat(std::abs(stubFarEnd - turnY), 9)};
                spaces.push_back(stub);
            }
            ConnectionReservedSpace lead;
            lead.isTerminal = true;
            lead.winding = windingName;
            lead.parallel = parallel;
            lead.section = connectingTurn.get_section().value_or("");
            lead.layer = "";
            lead.coordinates = {roundFloat((turnX + windowOuterX) / 2, 9), edgeY};
            lead.dimensions = {roundFloat(windowOuterX - turnX + wireOuterWidth, 9), wireOuterHeight};
            spaces.push_back(lead);
        };
        for (int64_t parallel = 0; parallel < numberParallels; ++parallel) {
            auto key = std::make_pair(windingName, parallel);
            if (entranceTurnByWindingParallel.count(key)) {
                addTerminalLead(entranceTurnByWindingParallel.at(key), parallel);
            }
            if (exitTurnByWindingParallel.count(key)) {
                addTerminalLead(exitTurnByWindingParallel.at(key), parallel);
            }
        }

        std::vector<Layer> windingLayers;
        for (const auto& layer : allLayers) {
            if (layer.get_partial_windings()[0].get_winding() == windingName) {
                windingLayers.push_back(layer);
            }
        }
        if (windingLayers.size() < 2) {
            continue;
        }
        std::sort(windingLayers.begin(), windingLayers.end(), [&](const Layer& a, const Layer& b) {
            size_t orderA = layerElectricalOrder.count(a.get_name()) ? layerElectricalOrder.at(a.get_name()) : 0;
            size_t orderB = layerElectricalOrder.count(b.get_name()) ? layerElectricalOrder.at(b.get_name()) : 0;
            return orderA < orderB;
        });

        for (size_t i = 0; i + 1 < windingLayers.size(); ++i) {
            double radialA = windingLayers[i].get_coordinates()[0];
            double radialB = windingLayers[i + 1].get_coordinates()[0];
            double radialLow = std::min(radialA, radialB);
            double radialHigh = std::max(radialA, radialB);

            // The intervening conduction layers a continuation crosses are fixed by the layer geometry
            // (same for every parallel); each parallel routes its own wire across them, so the squeezes
            // are emitted once PER PARALLEL below and stack under the height-based blocking.
            std::vector<const Layer*> interveningLayers;
            for (const auto& crossed : allLayers) {
                double radial = crossed.get_coordinates()[0];
                if (radial > radialLow + 1e-12 && radial < radialHigh - 1e-12) {
                    interveningLayers.push_back(&crossed);
                }
            }
            bool crossesIntervening = !interveningLayers.empty();

            // Convention: layer-to-layer connections route along the TOP window edge. Windings start
            // from the bottom and wind up, so a layer finishes at the top and the connection to the
            // next layer naturally sits there; squeeze the crossed layers at the top edge to match.
            // (Squeeze at the WINDOW edge, not the crossed layer's own edge — once the layer is
            // centred/shrunk its edge is where its end turn sits, which would put the slot on a turn.)
            double routeEdgeY = roundFloat(windowTopY - wireOuterHeight / 2, 9);
            WindingOrder windingOrder = get_winding_order(windingLayers[i].get_section().value());

            // Each parallel is its own conductor: it has its own last-turn-of-layer-i and
            // first-turn-of-layer-(i+1), its own crossing squeezes, and its own drawn link.
            for (int64_t parallel = 0; parallel < numberParallels; ++parallel) {
                auto exitKey = std::make_pair(windingLayers[i].get_name(), parallel);
                auto entryKey = std::make_pair(windingLayers[i + 1].get_name(), parallel);
                if (!lastTurnByLayerParallel.count(exitKey) || !firstTurnByLayerParallel.count(entryKey)) {
                    continue;
                }
                const auto& exitTurn = lastTurnByLayerParallel.at(exitKey);
                const auto& entryTurn = firstTurnByLayerParallel.at(entryKey);

                // Per-layer squeeze: this parallel's interleaved continuation crosses (and squeezes)
                // each intervening layer. These entries (layer set) drive the filling factor and are
                // NOT drawn — the link itself is drawn below.
                for (const Layer* crossed : interveningLayers) {
                    ConnectionReservedSpace squeeze;
                    squeeze.section = crossed->get_section().value();
                    squeeze.layer = crossed->get_name();
                    squeeze.winding = windingName;
                    squeeze.parallel = parallel;
                    squeeze.coordinates = {crossed->get_coordinates()[0], routeEdgeY};
                    squeeze.dimensions = {wireOuterWidth, wireOuterHeight};
                    spaces.push_back(squeeze);
                }
                // The continuation does NOT reserve a slot on its endpoint layers (the source's end turn
                // and the destination's start turn): the link is those turns' own wire continuing, so
                // they sit at the edge and the link is drawn between them. Only the intervening layers
                // it routes OVER lose a slot. This lets the endpoint layers fill to the edge.

                // Draw this parallel's connection from the last turn of this layer to the first turn of
                // the next: a U winding turns around (orthogonal L), a Z winding runs straight back to
                // the next layer's start (single diagonal — adjacent = dragback, interleaved =
                // continuation). Every consecutive layer pair of every parallel is drawn.
                double x1 = exitTurn.get_coordinates()[0];
                double y1 = exitTurn.get_coordinates()[1];
                double x2 = entryTurn.get_coordinates()[0];
                double y2 = entryTurn.get_coordinates()[1];

                if (windingOrder == WindingOrder::Z) {
                    // Z: the wire runs straight from one turn to the next, so draw a single diagonal
                    // link (a rotated rectangle from centre to centre).
                    double deltaX = x2 - x1;
                    double deltaY = y2 - y1;
                    double length = std::sqrt(deltaX * deltaX + deltaY * deltaY);
                    ConnectionReservedSpace diagonal;
                    diagonal.winding = windingName;
                    diagonal.parallel = parallel;
                    diagonal.section = windingLayers[i].get_section().value_or("");
                    diagonal.layer = "";  // drawn link; squeeze handled by the per-layer entries above
                    diagonal.coordinates = {roundFloat((x1 + x2) / 2, 9), roundFloat((y1 + y2) / 2, 9)};
                    diagonal.dimensions = {roundFloat(length, 9), wireOuterHeight};
                    diagonal.rotation = roundFloat(std::atan2(deltaY, deltaX) * 180.0 / std::numbers::pi, 6);
                    spaces.push_back(diagonal);
                }
                else if (crossesIntervening) {
                    // U interleaved continuation: route along the window edge it reserves — a vertical
                    // stub from each turn up/down to the edge, and a horizontal run along the edge
                    // across the intervening layer(s) — so the wire never cuts through the crossed
                    // layer's turns (centre-to-centre would clip them).
                    auto pushLink = [&](double cx, double cy, double w, double h) {
                        ConnectionReservedSpace seg;
                        seg.winding = windingName;
                        seg.parallel = parallel;
                        seg.section = windingLayers[i].get_section().value_or("");
                        seg.layer = "";
                        seg.coordinates = {roundFloat(cx, 9), roundFloat(cy, 9)};
                        seg.dimensions = {roundFloat(w, 9), roundFloat(h, 9)};
                        spaces.push_back(seg);
                    };
                    // Verticals start at the turn CENTRE (y1 / y2), not covering the whole turn, and
                    // overlap the horizontal by half a wire at the corner.
                    if (std::abs(routeEdgeY - y1) > 0.5 * wireOuterHeight) {
                        double far1 = routeEdgeY + ((routeEdgeY >= y1) ? 1.0 : -1.0) * wireOuterHeight / 2;
                        pushLink(x1, (y1 + far1) / 2, wireOuterWidth, std::abs(far1 - y1));
                    }
                    pushLink((x1 + x2) / 2, routeEdgeY, std::abs(x2 - x1) + wireOuterWidth, wireOuterHeight);
                    if (std::abs(y2 - routeEdgeY) > 0.5 * wireOuterHeight) {
                        double far2 = routeEdgeY + ((routeEdgeY >= y2) ? 1.0 : -1.0) * wireOuterHeight / 2;
                        pushLink(x2, (y2 + far2) / 2, wireOuterWidth, std::abs(far2 - y2));
                    }
                }
                else {
                    // U adjacent layers: route the wire orthogonally — a horizontal stretch from the
                    // source turn out to the destination layer's radial position, then a vertical
                    // stretch down/up to the destination turn when the two are at different heights.
                    // The horizontal runs half a wire past the corner and the vertical is pulled back
                    // half a wire so the bend reads as one continuous wire.
                    bool needVertical = std::abs(y2 - y1) > 0.5 * wireOuterHeight;
                    ConnectionReservedSpace horizontal;
                    horizontal.winding = windingName;
                    horizontal.parallel = parallel;
                    horizontal.section = windingLayers[i].get_section().value_or("");
                    horizontal.layer = "";
                    if (needVertical) {
                        double horizontalDirection = (x2 >= x1) ? 1.0 : -1.0;
                        horizontal.coordinates = {roundFloat((x1 + x2) / 2 + horizontalDirection * wireOuterWidth / 4, 9), roundFloat(y1, 9)};
                        horizontal.dimensions = {roundFloat(std::abs(x2 - x1) + wireOuterWidth / 2, 9), wireOuterHeight};
                    }
                    else {
                        horizontal.coordinates = {roundFloat((x1 + x2) / 2, 9), roundFloat(y1, 9)};
                        horizontal.dimensions = {roundFloat(std::abs(x2 - x1), 9), wireOuterHeight};
                    }
                    spaces.push_back(horizontal);

                    if (needVertical) {
                        double verticalDirection = (y2 >= y1) ? 1.0 : -1.0;
                        ConnectionReservedSpace vertical;
                        vertical.winding = windingName;
                        vertical.parallel = parallel;
                        vertical.section = windingLayers[i].get_section().value_or("");
                        vertical.layer = "";
                        vertical.coordinates = {roundFloat(x2, 9), roundFloat((y1 + y2) / 2 + verticalDirection * wireOuterHeight / 4, 9)};
                        vertical.dimensions = {wireOuterWidth, roundFloat(std::abs(y2 - y1) - wireOuterHeight / 2, 9)};
                        spaces.push_back(vertical);
                    }
                }
            }
        }
    }
    // Transpose the produced rectangles back to real coordinates for the contiguous case. A reflection
    // across y=x swaps the centre's x↔y and maps a rectangle drawn at angle θ (width×height) to one at
    // angle (90°−θ) with the same width×height — so swap the coordinates and set rotation = 90 − θ.
    if (layersAreContiguous) {
        for (auto& space : spaces) {
            if (space.coordinates.size() >= 2) {
                std::swap(space.coordinates[0], space.coordinates[1]);
            }
            space.rotation = roundFloat(90.0 - space.rotation, 6);
        }
    }
    return spaces;
}

std::map<std::string, std::pair<uint64_t, uint64_t>> Coil::compute_connection_blocked_slots_per_layer() {
    std::map<std::string, std::pair<uint64_t, uint64_t>> blockedSlotsPerLayer;  // layer name -> {top, bottom}
    if (!get_layers_description()) {
        return blockedSlotsPerLayer;
    }
    auto layers = get_layers_description().value();
    auto wires = get_wires();
    std::map<std::string, double> layerCenterHeight;
    std::map<std::string, double> layerWireHeight;  // crossed layer's own wire height (turn pitch)
    for (const auto& layer : layers) {
        // Turn blocking is modelled along the axial (y) direction, i.e. for OVERLAPPING layers only;
        // contiguous-layer leads are drawn but do not (yet) displace turns, so skip those layers here.
        if (layer.get_orientation() != WindingOrientation::OVERLAPPING) {
            continue;
        }
        layerCenterHeight[layer.get_name()] = layer.get_coordinates()[1];
        // Insulation layers have no partial windings; only conduction layers have a wire/turn pitch.
        if (layer.get_type() == ElectricalType::CONDUCTION && !layer.get_partial_windings().empty()) {
            size_t windingIndex = get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
            layerWireHeight[layer.get_name()] = wires[windingIndex].get_maximum_outer_height();
        }
    }
    // At a given edge of a layer, the leads of ONE parallel pile together along the window margin, so
    // that parallel's height cost is the TALLEST of its connections there (the max). The leads of
    // DIFFERENT parallels of a bifilar/N-filar group are separate physical conductors that stack on
    // top of each other, so their heights SUM. The stacked height is then converted to turn slots of
    // the CROSSED layer's own wire: a thick stack over a thin layer displaces several thin turns
    // (ceil), a thin stack over a thick layer still costs at least one thick turn. Z dragbacks route
    // diagonally and do not displace turns; only terminals block in Z. (For a single parallel this
    // reduces to the previous max-then-ceil rule.)
    std::map<std::string, std::map<int64_t, std::pair<double, double>>> perParallelMaxHeight;  // layer -> parallel -> {top, bottom}
    for (const auto& space : get_connection_reserved_spaces()) {
        if (space.layer.empty()) {
            continue;
        }
        if (!space.isTerminal && get_winding_order(space.section) != WindingOrder::U) {
            continue;
        }
        auto found = layerCenterHeight.find(space.layer);
        if (found == layerCenterHeight.end()) {
            continue;
        }
        auto& edges = perParallelMaxHeight[space.layer][space.parallel];
        if (space.coordinates[1] >= found->second) {
            edges.first = std::max(edges.first, space.dimensions[1]);
        }
        else {
            edges.second = std::max(edges.second, space.dimensions[1]);
        }
    }
    for (const auto& [layerName, perParallel] : perParallelMaxHeight) {
        double crossedWireHeight = layerWireHeight.count(layerName) ? layerWireHeight.at(layerName) : 0.0;
        if (crossedWireHeight <= 0) {
            continue;
        }
        double topHeight = 0, bottomHeight = 0;
        for (const auto& [parallel, edges] : perParallel) {
            topHeight += edges.first;
            bottomHeight += edges.second;
        }
        auto slots = [&](double connectionHeight) -> uint64_t {
            return connectionHeight > 1e-12 ? uint64_t(std::ceil(connectionHeight / crossedWireHeight - 1e-9)) : 0u;
        };
        blockedSlotsPerLayer[layerName] = {slots(topHeight), slots(bottomHeight)};
    }
    return blockedSlotsPerLayer;
}

void Coil::redistribute_section_turns_for_blocking() {
    if (!get_sections_description()) {
        return;
    }
    auto sections = get_sections_description().value();
    auto wirePerWinding = get_wires();
    size_t numberWindings = get_functional_description().size();

    auto blockedFor = [&](const std::string& sectionName, size_t layer) -> uint64_t {
        auto it = _connectionBlockedSlotsPerLayer.find(sectionName + " layer " + std::to_string(layer));
        if (it == _connectionBlockedSlotsPerLayer.end()) {
            return 0u;
        }
        return it->second.first + it->second.second;
    };

    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        double wireHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
        if (wireHeight <= 0) {
            continue;
        }
        // Each parallel of a bifilar/N-filar group is wound side by side, so a layer holds an equal
        // number of turns of every parallel: its capacity in PER-PARALLEL turns ("rows") is the
        // physical capacity divided by the parallel count. We redistribute in per-parallel turns
        // (matching get_number_turns, which is per parallel) so interior sections end on whole blocked
        // layers and the remainder is pushed to the outermost section. (K=1 keeps the previous result.)
        int64_t numberParallels = int64_t(get_number_parallels(windingIndex));

        // This winding's conduction sections, in wound (radial) order.
        std::vector<size_t> windingSections;
        for (size_t s = 0; s < sections.size(); ++s) {
            if (sections[s].get_type() == ElectricalType::CONDUCTION
                && sections[s].get_layers_orientation() == WindingOrientation::OVERLAPPING
                && get_winding_index_by_name(sections[s].get_partial_windings()[0].get_winding()) == windingIndex) {
                windingSections.push_back(s);
            }
        }
        if (windingSections.size() < 2) {
            continue;  // nothing to redistribute
        }

        uint64_t totalTurns = get_number_turns(windingIndex);  // per parallel
        uint64_t remaining = totalTurns;

        // Per-parallel capacity ("rows") of a blocked layer: physical capacity (maxTpl − blocked) split
        // evenly across the parallels wound side by side.
        auto rowsCapacity = [&](const std::string& sectionName, size_t layer, uint64_t maximumTurnsPerLayer) -> uint64_t {
            uint64_t blocked = std::min<uint64_t>(blockedFor(sectionName, layer), maximumTurnsPerLayer - 1);
            return uint64_t((maximumTurnsPerLayer - blocked) / numberParallels);
        };

        for (size_t k = 0; k < windingSections.size(); ++k) {
            auto& section = sections[windingSections[k]];
            // Physical turns per layer; at least one row of every parallel.
            uint64_t maximumTurnsPerLayer = std::max<uint64_t>(numberParallels, uint64_t(std::floor(section.get_dimensions()[1] / wireHeight)));
            size_t sectionsRemaining = windingSections.size() - k;

            uint64_t sectionTurns;  // per parallel
            if (sectionsRemaining == 1) {
                // Outermost section of the winding absorbs whatever is left (a partial outer layer
                // is acceptable; only interior orphans are the problem).
                sectionTurns = remaining;
            }
            else {
                // Fill complete blocked layers up to this section's fair share, so it ends on a layer
                // boundary (no interior orphan). Always leave at least one turn for each later section.
                uint64_t fairShare = uint64_t(std::round(double(remaining) / double(sectionsRemaining)));
                uint64_t turns = 0;
                size_t layer = 0;
                while (true) {
                    uint64_t capacity = rowsCapacity(section.get_name(), layer, maximumTurnsPerLayer);
                    if (capacity == 0) {
                        break;  // layer too thin to hold one row of every parallel
                    }
                    if (turns + capacity <= fairShare && (remaining - (turns + capacity)) >= (sectionsRemaining - 1)) {
                        turns += capacity;
                        layer++;
                    }
                    else {
                        break;
                    }
                }
                if (turns == 0) {
                    // Even a single full layer exceeds the fair share: take one layer anyway so this
                    // interior section is not left with a fractional layer.
                    uint64_t capacity0 = std::max<uint64_t>(rowsCapacity(section.get_name(), 0, maximumTurnsPerLayer), 1);
                    turns = std::min<uint64_t>(capacity0, remaining - (sectionsRemaining - 1));
                    if (turns == 0) {
                        turns = 1;
                    }
                }
                sectionTurns = turns;
            }
            remaining -= sectionTurns;

            // Side by side: every parallel gets the same per-parallel share of this section.
            std::vector<double> proportion(numberParallels, double(sectionTurns) / double(totalTurns));
            auto partialWindings = section.get_partial_windings();
            partialWindings[0].set_parallels_proportion(proportion);
            section.set_partial_windings(partialWindings);
        }
    }
    set_sections_description(sections);
}

void Coil::align_blocked_layer_turns() {
    if (_connectionBlockedSlotsPerLayer.empty() || !get_layers_description() || !get_turns_description()) {
        return;
    }
    auto bobbin = resolve_bobbin();
    if (bobbin.get_winding_window_shape() != WindingWindowShape::RECTANGULAR) {
        return;
    }
    auto windingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
    double windowCenterY = windingWindow.get_coordinates().value()[1];
    double windowHalfHeight = windingWindow.get_height().value() / 2;
    double windowTopY = windowCenterY + windowHalfHeight;
    double windowBottomY = windowCenterY - windowHalfHeight;
    auto wires = get_wires();
    auto layers = get_layers_description().value();
    auto turns = get_turns_description().value();

    for (auto& layer : layers) {
        if (layer.get_type() != ElectricalType::CONDUCTION || layer.get_orientation() != WindingOrientation::OVERLAPPING) {
            continue;
        }
        // Reposition this layer's turns to leave exactly the blocked slots free at each edge: spread the
        // turns evenly across the UNBLOCKED band [windowBottom + blockedBottom slots, windowTop −
        // blockedTop slots]. Even spacing packs a full layer (step == wire height) and spreads a partial
        // one — in both cases keeping every turn clear of the slots the connection leads route through.
        // delimit re-centres each layer over the full window height, which would otherwise drop the edge
        // turns back under the leads, so this runs after delimit and overrides its vertical placement.
        auto found = _connectionBlockedSlotsPerLayer.find(layer.get_name());
        if (found == _connectionBlockedSlotsPerLayer.end()) {
            continue;  // unblocked layer: leave delimit's centring (fills the full height)
        }
        uint64_t blockedTop = found->second.first;
        uint64_t blockedBottom = found->second.second;
        if (blockedTop == 0 && blockedBottom == 0) {
            continue;
        }
        size_t windingIndex = get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
        double wireHeight = wires[windingIndex].get_maximum_outer_height();

        // This layer's turns, in current (wound) vertical order.
        std::vector<size_t> layerTurns;
        for (size_t t = 0; t < turns.size(); ++t) {
            if (turns[t].get_layer() && turns[t].get_layer().value() == layer.get_name()) {
                layerTurns.push_back(t);
            }
        }
        if (layerTurns.empty()) {
            continue;
        }
        std::sort(layerTurns.begin(), layerTurns.end(), [&](size_t a, size_t b) {
            return turns[a].get_coordinates()[1] < turns[b].get_coordinates()[1];
        });

        // Centres of the first and last usable slots once the blocked slots are reserved at each edge.
        double bandBottom = roundFloat(windowBottomY + double(blockedBottom) * wireHeight + wireHeight / 2, 9);
        double bandTop = roundFloat(windowTopY - double(blockedTop) * wireHeight - wireHeight / 2, 9);
        size_t numberTurnsInLayer = layerTurns.size();
        double step = (numberTurnsInLayer > 1) ? (bandTop - bandBottom) / double(numberTurnsInLayer - 1) : 0.0;
        if (numberTurnsInLayer > 1 && step < wireHeight) {
            step = wireHeight;  // capacity should make this unreachable; never overlap turns
        }
        for (size_t k = 0; k < numberTurnsInLayer; ++k) {
            double newY = (numberTurnsInLayer == 1)
                ? roundFloat((bandBottom + bandTop) / 2, 9)
                : roundFloat(bandBottom + double(k) * step, 9);
            auto coords = turns[layerTurns[k]].get_coordinates();
            coords[1] = newY;
            turns[layerTurns[k]].set_coordinates(coords);
        }
        layer.set_coordinates(std::vector<double>{layer.get_coordinates()[0], roundFloat((bandBottom + bandTop) / 2, 9), 0});
    }
    set_layers_description(layers);
    set_turns_description(turns);
}

void Coil::apply_connection_reserved_space() {
    if (!get_sections_description() || !get_layers_description()) {
        return;
    }
    auto spaces = get_connection_reserved_spaces();
    if (spaces.empty()) {
        return;
    }

    auto layers = get_layers_description().value();

    // Axial height occupied by connection leads on each conduction layer. Every reserved space that
    // names a layer (an inter-layer transition, or a terminal lead passing over that layer) occupies
    // one wire diameter of that layer's available height. Free-space terminal segments (no layer)
    // are drawn but reserve no layer space.
    std::map<std::string, double> heightReservedPerLayer;
    for (const auto& space : spaces) {
        if (space.layer.empty()) {
            continue;
        }
        heightReservedPerLayer[space.layer] += space.dimensions[1];
    }

    // Reduce each affected layer's available height: the filling factor scales inversely with the
    // remaining height. A value above 1 means the leads no longer fit alongside the turns (the layer
    // is over-subscribed and the build needs more space).
    std::map<std::string, double> reservedAreaPerSection;
    for (auto& layer : layers) {
        auto it = heightReservedPerLayer.find(layer.get_name());
        if (it == heightReservedPerLayer.end()) {
            continue;
        }
        double layerWidth = layer.get_dimensions()[0];
        double layerHeight = layer.get_dimensions()[1];
        if (layerHeight <= 0) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR, "Non-positive layer height while applying connection reserved space to layer " + layer.get_name());
        }
        if (!layer.get_filling_factor()) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR, "Layer filling factor not set before applying connection reserved space to layer " + layer.get_name());
        }
        double reservedHeight = it->second;
        double availableHeight = std::max(layerHeight - reservedHeight, layerHeight * 0.01);
        layer.set_filling_factor(layer.get_filling_factor().value() * layerHeight / availableHeight);
        if (layer.get_section()) {
            reservedAreaPerSection[layer.get_section().value()] += reservedHeight * layerWidth;
        }
    }
    set_layers_description(layers);

    auto sections = get_sections_description().value();
    for (auto& section : sections) {
        auto it = reservedAreaPerSection.find(section.get_name());
        if (it == reservedAreaPerSection.end()) {
            continue;
        }
        double sectionArea = section.get_dimensions()[0] * section.get_dimensions()[1];
        if (sectionArea <= 0) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR, "Non-positive section area while applying connection reserved space to section " + section.get_name());
        }
        if (!section.get_filling_factor()) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR, "Section filling factor not set before applying connection reserved space to section " + section.get_name());
        }
        section.set_filling_factor(section.get_filling_factor().value() + it->second / sectionArea);
    }
    set_sections_description(sections);
}

CoilAlignment Coil::get_section_alignment() {
    auto bobbin = resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        return _sectionAlignment;
    }
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    // Multi-window: returns alignment from window 0 (single alignment for the
    // whole coil; per-window alignment is a v2 refinement).
    if (windingWindows.size() > 0) {
        if (windingWindows[0].get_sections_alignment()) {
            auto alignment = windingWindows[0].get_sections_alignment().value();
            return alignment;
        }
        if (bobbin.get_winding_window_shape() == WindingWindowShape::ROUND) {
            return defaults.defaultRoundWindowSectionsAlignment;
        }
    }
    return _sectionAlignment;
}

bool Coil::fast_wind() {
    _strict = false;

    wind_by_sections();
    if (!get_sections_description()) {
        return false;
    }
    wind_by_layers();
    if (!get_layers_description()) {
        return false;
    }
    {
        // RAII (ABT #113 sweep): exception-safe replacement for the manual
        // save/set/restore — wind_by_turns can throw.
        SettingsGuard<bool> includeAdditionalCoordinatesGuard(settings, &Settings::get_coil_include_additional_coordinates, &Settings::set_coil_include_additional_coordinates, false);
        wind_by_turns();
    }

    if (!get_turns_description()) {
        return false;
    }
    // Multi-column winding: mirror negative-x-window groups into place (no-op for
    // single-window coils).
    apply_group_window_sides();
    return true;
}

bool Coil::unwind() {
    set_sections_description(std::nullopt);
    set_layers_description(std::nullopt);
    set_turns_description(std::nullopt);
    return true;
}

bool Coil::wind() {
    std::vector<double> proportionPerWinding;

    proportionPerWinding = make_equal_proportion_per_winding(get_functional_description().size());
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    return wind(proportionPerWinding, pattern, _interleavingLevel);
}

bool Coil::wind(size_t repetitions){
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    auto proportionPerWinding = make_equal_proportion_per_winding(get_functional_description().size());
    return wind(proportionPerWinding, pattern, repetitions);
}

bool Coil::wind(std::vector<size_t> pattern, size_t repetitions){
    auto proportionPerWinding = make_equal_proportion_per_winding(get_functional_description().size());
    return wind(proportionPerWinding, pattern, repetitions);
}

std::vector<size_t> Coil::extract_stack_up(std::vector<Section> sections) {
    std::vector<size_t> stackUp;
    for (auto section : sections) {
        size_t windingIndex = get_winding_index_by_name(section.get_partial_windings()[0].get_winding());
        stackUp.push_back(windingIndex);
    }
    return stackUp;
}

bool Coil::wind(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    bool windEvenIfNotFit = settings.get_coil_wind_even_if_not_fit();
    bool delimitAndCompact = settings.get_coil_delimit_and_compact();
    bool tryRewind = settings.get_coil_try_rewind();

    std::string bobbinName = "";
    if (std::holds_alternative<std::string>(get_bobbin())) {
        bobbinName = std::get<std::string>(get_bobbin());
        if (bobbinName != "Dummy") {
            auto bobbinData = find_bobbin_by_name(std::get<std::string>(get_bobbin()));
            set_bobbin(bobbinData);
        }
    }
    _currentProportionPerWinding = proportionPerWinding;
    _currentPattern = pattern;
    _currentRepetitions = repetitions;

    if (bobbinName != "Dummy") {
        bool wind = true;                
        for (auto& winding : get_mutable_functional_description()) {
            if (std::holds_alternative<std::string>(winding.get_wire())) {
                std::string wireName = std::get<std::string>(winding.get_wire());
                if (wireName == "Dummy" || wireName.empty()) {
                    wind = false;
                    break;
                }
                auto wire = find_wire_by_name(wireName);
                winding.set_wire(wire);
            }
        }

        if (wind) {
            set_sections_description(std::nullopt);
            set_layers_description(std::nullopt);
            set_turns_description(std::nullopt);
            // Start every wind from the ideal geometry: real-winding turn blocking, if any, is
            // re-derived and re-applied below only when the real-geometry setting is on.
            _applyConnectionBlocking = false;
            _connectionBlockedSlotsPerLayer.clear();

            // Special case: toroid with one physical turn whose wire OD
            // exceeds the inner-hole radius. The wire cannot be wound on
            // the inner wall, so place it at the geometric centre of the
            // hole. Skips the sections/layers fit pipeline entirely.
            if (can_build_centered_single_turn_toroidal()) {
                logEntry("Building centered single-turn toroidal", "Coil", 2);
                return build_centered_single_turn_toroidal();
            }

            if (_insulationSections.size() == 0) {

                if (_inputs) {
                    if (_inputs->get_design_requirements().get_insulation()) {
                        logEntry("Calculating Required Insulation", "Coil", 2);
                        calculate_insulation();
                    }
                    else {
                        logEntry("Calculating Mechanical Insulation", "Coil", 2);
                        calculate_mechanical_insulation();
                    }
                }
                else {
                    logEntry("Calculating Mechanical Insulation", "Coil", 2);
                    calculate_mechanical_insulation();
                }
            }
            logEntry("Winding by sections", "Coil", 2);
            wind_by_sections(proportionPerWinding, pattern, repetitions);
            logEntry("Winding by layers", "Coil", 2);
            wind_by_layers();

            if (!get_layers_description()) {
                return false;
            }

            auto sections = get_sections_description().value();

            if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
                logEntry("Winding by turns", "Coil", 2);
                wind_by_turns();
                if (delimitAndCompact) {
                    logEntry("Delimiting and compacting", "Coil", 2);
                    delimit_and_compact();
                }
            }

            if (tryRewind && (!are_sections_and_layers_fitting() || !get_turns_description())) {
                logEntry("Trying to rewind", "Coil", 2);
                try_rewind();
            }
        }
    }
    // The wind success is decided on the ideal geometry. When real winding geometry is enabled, the
    // space reserved by connection leads is layered on afterwards (filling factors, Painter, losses)
    // so it never changes whether the ideal winding fit.
    bool result = are_sections_and_layers_fitting() && bool(get_turns_description());
    if (result && settings.get_coil_use_real_winding_geometry()) {
        // Turn blocking is GLOBAL to the winding window: a connection lead routes through the whole
        // window and removes a turn slot from every conduction layer it crosses, regardless of which
        // section/winding owns that layer. Derive that per-layer incidence from the wound geometry,
        // then re-wind so each crossed layer frees its blocked top/bottom slots and the affected
        // sections grow by extra layers to fit. Added layers shift radial positions (changing which
        // leads cross what), so iterate to a fixpoint. Entirely gated behind the real-winding flag —
        // ideal winding never enters this loop, so its geometry is unchanged.
        logEntry("Applying real winding geometry (global turn blocking)", "Coil", 2);
        const size_t maximumBlockingIterations = 8;
        for (size_t blockingIteration = 0; blockingIteration < maximumBlockingIterations; ++blockingIteration) {
            auto freshBlocked = compute_connection_blocked_slots_per_layer();
            // Accumulate the blocked slots MONOTONICALLY (element-wise max) instead of replacing them.
            // Freeing the outermost layer's exit slot spills a turn into a new outer layer, which then
            // makes the previous layer's top look unblocked — so a plain replace flip-flops between an
            // N-layer and an (N+1)-layer state and never converges. Taking the running max makes the
            // reserved slots non-decreasing, so the build is stable and converges in a few iterations
            // (at worst it reserves a slot that ends unused — conservative, never a collision).
            bool changed = false;
            for (const auto& [layerName, edges] : freshBlocked) {
                auto& accumulated = _connectionBlockedSlotsPerLayer[layerName];
                uint64_t newTop = std::max(accumulated.first, edges.first);
                uint64_t newBottom = std::max(accumulated.second, edges.second);
                if (newTop != accumulated.first || newBottom != accumulated.second) {
                    accumulated = {newTop, newBottom};
                    changed = true;
                }
            }
            if (!changed) {
                break;
            }
            _applyConnectionBlocking = true;
            wind_by_sections(proportionPerWinding, pattern, repetitions);
            // Re-split each winding's turns across its sections so interior layers fill completely
            // under the blocking (no orphan spillover turns); the outermost section absorbs the rest.
            redistribute_section_turns_for_blocking();
            wind_by_layers();
            if (!get_layers_description()) {
                break;
            }
            if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
                wind_by_turns();
                if (delimitAndCompact) {
                    delimit_and_compact();
                }
            }
        }
        result = are_sections_and_layers_fitting() && bool(get_turns_description());
        // For Z-wound layers, pack the turns against the unblocked edge so the slots freed by blocking
        // line up with the terminal leads at the window edges (delimit otherwise centres them, leaving
        // a turn-sized gap on the unblocked side and clipping the terminal slot on the blocked side).
        align_blocked_layer_turns();
        logEntry("Applying real winding geometry (connection reserved space)", "Coil", 2);
        apply_connection_reserved_space();
    }
    if (result) {
        // Multi-column winding: groups were wound in the +x window-local frame;
        // mirror the ones whose winding window sits on the negative-x side into
        // their real position. No-op for single-window coils.
        apply_group_window_sides();
    }
    return result;
}

bool Coil::wind_planar(std::vector<size_t> stackUp, std::optional<double> borderToWireDistance, std::map<size_t, double> wireToWireDistance, std::map<std::pair<size_t, size_t>, double> insulationThickness, double coreToLayerDistance) {
    bool windEvenIfNotFit = settings.get_coil_wind_even_if_not_fit();
    bool delimitAndCompact = settings.get_coil_delimit_and_compact();
    std::string bobbinName = "";
    if (std::holds_alternative<std::string>(get_bobbin())) {
        bobbinName = std::get<std::string>(get_bobbin());
        if (bobbinName != "Dummy") {
            auto bobbinData = find_bobbin_by_name(std::get<std::string>(get_bobbin()));
            set_bobbin(bobbinData);
        }
    }

    if (bobbinName != "Dummy") {
        bool wind = true;                
        for (auto& winding : get_mutable_functional_description()) {
            if (std::holds_alternative<std::string>(winding.get_wire())) {
                std::string wireName = std::get<std::string>(winding.get_wire());
                if (wireName == "Dummy" || wireName.empty()) {
                    wind = false;
                    break;
                }
                auto wire = find_wire_by_name(wireName);
                winding.set_wire(wire);
            }
        }

        if (wind) {
            set_groups_description(std::nullopt);
            set_sections_description(std::nullopt);
            set_layers_description(std::nullopt);
            set_turns_description(std::nullopt);

            // Note: Insulation clearance from InsulationCoordinator is for conductor-to-conductor
            // isolation (primary-secondary), not for conductor-to-core or turn-to-turn spacing.
            // Both borderToWireDistance and wireToWireDistance should use manufacturing defaults,
            // not the insulation clearance which applies to inter-winding spacing.

            if (!borderToWireDistance) {
                borderToWireDistance = defaults.minimumBorderToWireDistance;
            }
            for (size_t i = 0; i < get_functional_description().size(); ++i) {
                if (!wireToWireDistance.count(i)) {
                    wireToWireDistance[i] = defaults.minimumWireToWireDistance;
                } 
            }

            logEntry("Winding by sections", "Coil", 2);
            wind_by_planar_sections(stackUp, insulationThickness, coreToLayerDistance);
            logEntry("Winding by layers", "Coil", 2);
            wind_by_planar_layers();

            if (!get_layers_description()) {
                return false;
            }

            if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
                logEntry("Winding by turns", "Coil", 2);
                wind_by_planar_turns(borderToWireDistance.value(), wireToWireDistance);
                if (delimitAndCompact) {
                    logEntry("Delimiting and compacting", "Coil", 2);
                    delimit_and_compact();
                }
            }
        }
    }
    return are_sections_and_layers_fitting() && get_turns_description();
}

/**
 * @brief Determine optimal winding style for each winding based on geometry and AC performance.
 *
 * Winding styles:
 * - WIND_BY_CONSECUTIVE_TURNS: Groups turns of same parallel together (P0T0, P0T1, ... P1T0, P1T1, ...)
 * - WIND_BY_CONSECUTIVE_PARALLELS: Groups parallels of same turn together (P0T0, P1T0, ... P0T1, P1T1, ...)
 *
 * For AC performance, WIND_BY_CONSECUTIVE_PARALLELS is generally preferred when numberParallels > 1:
 * 1. Current sharing: Parallels of same turn see similar flux linkage → better current balance
 * 2. Proximity effect: Interleaved parallels reduce effective layer count → lower AC losses
 *
 * However, geometric constraints (divisibility) may override this preference.
 */
std::vector<WindingStyle> Coil::wind_by_consecutive_turns(std::vector<uint64_t> numberTurns, std::vector<uint64_t> numberParallels, std::vector<size_t> numberSlots) {
    std::vector<WindingStyle> windByConsecutiveTurns;
    for (size_t i = 0; i < numberTurns.size(); ++i) {
        if (numberSlots[i] <= 0) {
            throw InvalidInputException("Number of slots cannot be less than 1, please verify your isolation sides requirement");
        }
        
        // When turns < slots, we MUST use CONSECUTIVE_TURNS to distribute physical turns (parallels)
        // across slots. CONSECUTIVE_PARALLELS would put all turns in the first slots, leaving rest empty.
        if (numberTurns[i] < numberSlots[i] && numberParallels[i] > 1) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            log("Winding " + std::to_string(i) + ": CONSECUTIVE_TURNS (turns < slots, must distribute parallels across slots).");
            continue;
        }

        // Case 1: Perfect fit - one turn per slot
        if (numberTurns[i] == numberSlots[i]) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS);
            log("Winding " + std::to_string(i) + ": CONSECUTIVE_PARALLELS (turns == slots, groups parallels of same turn).");
            continue;
        }
        
        // Case 2: Perfect fit - one parallel per slot
        if (numberParallels[i] == numberSlots[i]) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            log("Winding " + std::to_string(i) + ": CONSECUTIVE_TURNS (parallels == slots, groups turns of same parallel).");
            continue;
        }
        
        // Case 3: Multiple parallels - prefer CONSECUTIVE_PARALLELS for better current sharing
        // unless geometric constraints prevent it
        if (numberParallels[i] > 1) {
            // Check if CONSECUTIVE_PARALLELS is geometrically feasible
            // It works well when turns divide evenly into slots
            if (numberTurns[i] % numberSlots[i] == 0) {
                windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS);
                log("Winding " + std::to_string(i) + ": CONSECUTIVE_PARALLELS (multiple parallels, turns divisible by slots - better current sharing).");
                continue;
            }
            // Also prefer CONSECUTIVE_PARALLELS when parallels can be evenly distributed per turn group
            if (numberSlots[i] > 1 && (numberTurns[i] * numberParallels[i]) % numberSlots[i] == 0) {
                windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS);
                log("Winding " + std::to_string(i) + ": CONSECUTIVE_PARALLELS (multiple parallels, physical turns divisible - better proximity effect).");
                continue;
            }
        }
        
        // Case 4: Check divisibility for geometric fit
        if (numberParallels[i] % numberSlots[i] == 0) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            log("Winding " + std::to_string(i) + ": CONSECUTIVE_TURNS (parallels divisible by slots).");
            continue;
        }
        if (numberTurns[i] % numberSlots[i] == 0) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS);
            log("Winding " + std::to_string(i) + ": CONSECUTIVE_PARALLELS (turns divisible by slots).");
            continue;
        }
        
        // Case 5: Default - use CONSECUTIVE_PARALLELS if multiple parallels (for current sharing),
        // otherwise CONSECUTIVE_TURNS
        if (numberParallels[i] > 1) {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS);
            log("Winding " + std::to_string(i) + ": CONSECUTIVE_PARALLELS (default for multiple parallels - prioritizes current sharing).");
        } else {
            windByConsecutiveTurns.push_back(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            log("Winding " + std::to_string(i) + ": CONSECUTIVE_TURNS (single parallel, default style).");
        }
    }
    return windByConsecutiveTurns;
}

/**
 * @brief Determine optimal winding style for a single winding/layer.
 *
 * See multi-winding version for detailed explanation of style tradeoffs.
 * For AC performance with multiple parallels, CONSECUTIVE_PARALLELS is preferred
 * for better current sharing and reduced proximity effect.
 */
WindingStyle Coil::wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, size_t numberSlots) {
    // When turns < slots, we MUST use CONSECUTIVE_TURNS to distribute physical turns (parallels)
    // across slots. CONSECUTIVE_PARALLELS would put all turns in the first slots, leaving rest empty.
    if (numberTurns < numberSlots && numberParallels > 1) {
        log("Layer: CONSECUTIVE_TURNS (turns < slots, must distribute parallels across slots).");
        return WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
    }

    // Perfect fit cases
    if (numberTurns == numberSlots) {
        log("Layer: CONSECUTIVE_PARALLELS (turns == slots).");
        return WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS;
    }
    if (numberParallels == numberSlots) {
        log("Layer: CONSECUTIVE_TURNS (parallels == slots).");
        return WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
    }
    
    // Multiple parallels - prefer CONSECUTIVE_PARALLELS for current sharing
    if (numberParallels > 1) {
        if (numberTurns % numberSlots == 0) {
            log("Layer: CONSECUTIVE_PARALLELS (multiple parallels, turns divisible - better current sharing).");
            return WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS;
        }
        if ((numberTurns * numberParallels) % numberSlots == 0) {
            log("Layer: CONSECUTIVE_PARALLELS (multiple parallels, physical turns divisible - better proximity).");
            return WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS;
        }
    }
    
    // Divisibility-based selection
    if (numberParallels % numberSlots == 0) {
        log("Layer: CONSECUTIVE_TURNS (parallels divisible by slots).");
        return WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
    }
    if (numberTurns % numberSlots == 0) {
        log("Layer: CONSECUTIVE_PARALLELS (turns divisible by slots).");
        return WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS;
    }
    
    // Default: prefer CONSECUTIVE_PARALLELS for multiple parallels, otherwise CONSECUTIVE_TURNS
    if (numberParallels > 1) {
        log("Layer: CONSECUTIVE_PARALLELS (default for multiple parallels - current sharing priority).");
        return WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS;
    }
    log("Layer: CONSECUTIVE_TURNS (single parallel, default).");
    return WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
}


uint64_t Coil::get_number_turns(size_t windingIndex) const {
    return get_functional_description()[windingIndex].get_number_turns();
}

uint64_t Coil::get_number_parallels(size_t windingIndex) const {
    return get_functional_description()[windingIndex].get_number_parallels();
}

uint64_t Coil::get_number_turns(Section section) {
    uint64_t physicalTurnsInSection = 0;
    auto partialWinding = section.get_partial_windings()[0];  // TODO: Support multiwinding in layers
    auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
        physicalTurnsInSection += round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
    }
    return physicalTurnsInSection;
}

uint64_t Coil::get_number_turns(Layer layer) {
    uint64_t physicalTurnsInLayer = 0;
    auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
    auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
        physicalTurnsInLayer += round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
    }
    return physicalTurnsInLayer;
}

std::string Coil::get_name(size_t windingIndex) const {
    return get_functional_description()[windingIndex].get_name();
}

std::vector<uint64_t> Coil::get_number_turns() const {
    std::vector<uint64_t> numberTurns;
    for (auto & winding : get_functional_description()) {
        numberTurns.push_back(winding.get_number_turns());
    }
    return numberTurns;
}

void Coil::set_number_turns(std::vector<uint64_t> numberTurns) {
    for (size_t i=0; i< get_functional_description().size(); ++i) {
        get_mutable_functional_description()[i].set_number_turns(numberTurns[i]);
    }
}

std::vector<IsolationSide> Coil::get_isolation_sides() const {
    std::vector<IsolationSide> isolationSides;
    for (auto & winding : get_functional_description()) {
        isolationSides.push_back(winding.get_isolation_side());
    }
    return isolationSides;
}

void Coil::set_isolation_sides(std::vector<IsolationSide> isolationSides) {
    for (size_t i=0; i< get_functional_description().size(); ++i) {
        get_mutable_functional_description()[i].set_isolation_side(isolationSides[i]);
    }
}

std::vector<Layer> Coil::get_layers_by_section(std::string sectionName) const {
    auto layers = get_layers_description().value();
    std::vector<Layer> foundLayers;
    for (auto & layer : layers) {
        auto layerSectionName = layer.get_section().value();
        if (layerSectionName == sectionName) {
            foundLayers.push_back(layer);
        }
    }
    return foundLayers;
}

std::vector<Turn> Coil::get_turns_by_layer(std::string layerName) const {
    auto turns = get_turns_description().value();
    std::vector<Turn> foundTurns;
    for (auto & turn : turns) {
        auto turnLayerName = turn.get_layer().value();
        if (turnLayerName == layerName) {
            foundTurns.push_back(turn);
        }
    }
    return foundTurns;
}

std::vector<Turn> Coil::get_turns_by_winding(std::string windingName) const {
    auto turns = get_turns_description().value();
    std::vector<Turn> foundTurns;
    for (auto & turn : turns) {
        auto turnSectionName = turn.get_winding();
        if (turnSectionName == windingName) {
            foundTurns.push_back(turn);
        }
    }
    return foundTurns;
}

std::vector<Turn> Coil::get_turns_by_section(std::string sectionName) const {
    auto turns = get_turns_description().value();
    std::vector<Turn> foundTurns;
    for (auto & turn : turns) {
        auto turnSectionName = turn.get_section().value();
        if (turnSectionName == sectionName) {
            foundTurns.push_back(turn);
        }
    }
    return foundTurns;
}

std::vector<std::string> Coil::get_layers_names_by_winding(std::string windingName) const {
    auto layers = get_layers_description().value();
    std::vector<std::string> foundLayers;
    for (auto & layer : layers) {
        auto layerWindings = layer.get_partial_windings();
        for (auto & winding : layerWindings) {
            if (winding.get_winding() == windingName) {
                foundLayers.push_back(layer.get_name());
                break;
            }
        }
    }
    return foundLayers;
}

std::vector<std::string> Coil::get_layers_names_by_section(std::string sectionName) const {
    auto layers = get_layers_description().value();
    std::vector<std::string> foundLayers;
    for (auto & layer : layers) {
        auto layerSectionName = layer.get_section().value();
        if (layerSectionName == sectionName) {
            foundLayers.push_back(layer.get_name());
        }
    }
    return foundLayers;
}

std::vector<std::string> Coil::get_turns_names_by_layer(std::string layerName) const {
    auto turns = get_turns_description().value();
    std::vector<std::string> foundTurns;
    for (auto & turn : turns) {
        auto turnLayerName = turn.get_layer().value();
        if (turnLayerName == layerName) {
            foundTurns.push_back(turn.get_name());
        }
    }
    return foundTurns;
}

std::vector<std::string> Coil::get_turns_names_by_winding(std::string windingName) const {
    auto turns = get_turns_description().value();
    std::vector<std::string> foundTurns;
    for (auto & turn : turns) {
        auto turnWindingName = turn.get_winding();
        if (turnWindingName == windingName) {
            foundTurns.push_back(turn.get_name());
        }
    }
    return foundTurns;
}

std::vector<std::string> Coil::get_turns_names_by_section(std::string sectionName) const {
    auto turns = get_turns_description().value();
    std::vector<std::string> foundTurns;
    for (auto & turn : turns) {
        auto turnSectionName = turn.get_section().value();
        if (turnSectionName == sectionName) {
            foundTurns.push_back(turn.get_name());
        }
    }
    return foundTurns;
}
    
std::vector<size_t> Coil::get_turns_indexes_by_layer(std::string layerName) const {
    auto turns = get_turns_description().value();
    std::vector<size_t> foundTurns;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto turnLayerName = turns[turnIndex].get_layer().value();
        if (turnLayerName == layerName) {
            foundTurns.push_back(turnIndex);
        }
    }
    return foundTurns;
}

std::vector<size_t> Coil::get_turns_indexes_by_winding(std::string windingName) const {
    auto turns = get_turns_description().value();
    std::vector<size_t> foundTurns;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto turnWindingName = turns[turnIndex].get_winding();
        if (turnWindingName == windingName) {
            foundTurns.push_back(turnIndex);
        }
    }
    return foundTurns;
}

std::vector<size_t> Coil::get_turns_indexes_by_section(std::string sectionName) const {
    auto turns = get_turns_description().value();
    std::vector<size_t> foundTurns;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto turnSectionName = turns[turnIndex].get_section().value();
        if (turnSectionName == sectionName) {
            foundTurns.push_back(turnIndex);
        }
    }
    return foundTurns;
}

std::vector<Section> Coil::get_sections_by_group(std::string groupName) const {
    auto sections = get_sections_description().value();
    std::vector<Section> foundSections;
    for (auto & section : sections) {
        if (section.get_group()) {
            auto sectionSectionGroup = section.get_group().value();
            if (sectionSectionGroup == groupName) {
                foundSections.push_back(section);
            }
        }
    }
    return foundSections;
}

const std::vector<Section> Coil::get_sections_by_type(ElectricalType electricalType) const {
    auto sections = get_sections_description().value();
    std::vector<Section> foundSections;
    for (auto & section : sections) {
        auto sectionSectionType = section.get_type();
        if (sectionSectionType == electricalType) {
            foundSections.push_back(section);
        }
    }
    return foundSections;
}

const std::vector<Section> Coil::get_sections_by_winding(std::string windingName) const {
    auto sections = get_sections_description().value();
    std::vector<Section> foundSections;
    for (auto & section : sections) {
        for (auto & winding : section.get_partial_windings()) {
            if (winding.get_winding() == windingName) {
                foundSections.push_back(section);
            }
        }
    }
    return foundSections;
}

const Section Coil::get_section_by_name(std::string name) const {
    auto sections = get_sections_description().value();
    for (auto & section : sections) {
        if (section.get_name() == name) {
            return section;
        }
    }
    throw CoilException(ErrorCode::COIL_WINDING_ERROR, "Not found section with name:" + name);
}

const Layer Coil::get_layer_by_name(std::string name) const {
    if (!get_layers_description()) {
        throw CoilNotProcessedException("Coil is missing layers description");
    }

    auto layers = get_layers_description().value();
    for (auto & layer : layers) {
        if (layer.get_name() == name) {
            return layer;
        }
    }
    throw CoilException(ErrorCode::COIL_WINDING_ERROR, "Not found layer with name:" + name);
}


Turn Coil::get_turn_by_name(std::string name){
    if (_turnByName.count(name) == 0) {

        if (!get_turns_description()) {
            throw CoilNotProcessedException("Turns description not set, did you forget to wind?");
        }
        auto turns = get_turns_description().value();
        bool found = false;
        for (auto turn : turns) {
            if (turn.get_name() == name) {
                _turnByName[name] = turn;
                found = true;
                break;
            }
        }
        if (!found) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR, "No such a turn name: " + name);
        }
    }
    return _turnByName[name];
}

const std::vector<Layer> Coil::get_layers_by_type(ElectricalType electricalType) const {
    auto layers = get_layers_description().value();
    std::vector<Layer> foundLayers;
    for (auto & layer : layers) {
        auto layerSectionType = layer.get_type();
        if (layerSectionType == electricalType) {
            foundLayers.push_back(layer);
        }
    }
    return foundLayers;
}

std::vector<Layer> Coil::get_layers_by_winding_index(size_t windingIndex) {
    auto layers = get_layers_by_type(ElectricalType::CONDUCTION);
    std::vector<Layer> foundLayers;
    for (auto & layer : layers) {
        auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
        auto winding = get_winding_by_name(partialWinding.get_winding());
        auto layerWindingIndex = get_winding_index_by_name(partialWinding.get_winding());
        if (layerWindingIndex == windingIndex) {
            foundLayers.push_back(layer);
        }
    }
    return foundLayers;
}

std::vector<uint64_t> Coil::get_number_parallels() const {
    std::vector<uint64_t> numberParallels;
    for (auto & winding : get_functional_description()) {
        numberParallels.push_back(winding.get_number_parallels());
    }
    return numberParallels;
}

void Coil::set_number_parallels(std::vector<uint64_t> numberParallels){
    for (size_t i=0; i< get_functional_description().size(); ++i) {
        get_mutable_functional_description()[i].set_number_parallels(numberParallels[i]);
    }
}

Winding Coil::get_winding_by_name(std::string name) const {
    for (auto& Winding : get_functional_description()) {
        if (Winding.get_name() == name) {
            return Winding;
        }
    }
    throw CoilException(ErrorCode::COIL_WINDING_ERROR, "No such a winding name: " + name);
}

size_t Coil::get_winding_index_by_name(const std::string& name) {
    auto it = _windingIndexByName.find(name);
    if (it != _windingIndexByName.end()) {
        return it->second;
    }
    size_t index = get_winding_index_by_name(get_functional_description(), name);
    _windingIndexByName[name] = index;
    return index;
}

size_t Coil::get_winding_index_by_name(const std::vector<Winding>& functionalDescription, const std::string& name) {
    for (size_t i=0; i<functionalDescription.size(); ++i) {
        if (functionalDescription[i].get_name() == name) {
            return i;
        }
    }
    throw CoilException(ErrorCode::COIL_WINDING_ERROR, "No such a winding name in functionalDescription: " + name);
}

size_t Coil::get_turn_index_by_name(std::string name) {
    if (!get_turns_description()) {
        throw CoilNotProcessedException("Turns description not set, did you forget to wind?");
    }
    // Note: get_turns_description() returns an optional by value; .value() returns
    // a reference into that temporary. Copy out to avoid dangling reference.
    auto turns = get_turns_description().value();

    // Validate cache: the turns vector may have been replaced since the cache was
    // populated (e.g. re-winding). A stale index would cause out-of-bounds writes
    // in downstream code that does turns[turnIndex].set_*(...).
    auto it = _turnIndexByName.find(name);
    if (it != _turnIndexByName.end()) {
        if (it->second < turns.size() && turns[it->second].get_name() == name) {
            return it->second;
        }
        // Stale cache — clear it entirely since other entries may also be stale
        _turnIndexByName.clear();
    }

    for (size_t i = 0; i < turns.size(); ++i) {
        if (turns[i].get_name() == name) {
            _turnIndexByName[name] = i;
            return i;
        }
    }
    throw CoilException(ErrorCode::COIL_WINDING_ERROR, "No such a turn name: " + name);
}

size_t Coil::get_layer_index_by_name(std::string name) const {
    if (!get_layers_description()) {
        throw CoilNotProcessedException("Layers description not set, did you forget to wind?");
    }
    auto layers = get_layers_description().value();
    for (size_t i=0; i<layers.size(); ++i) {
        if (layers[i].get_name() == name) {
            return i;
        }
    }
    throw CoilException(ErrorCode::COIL_WINDING_ERROR, "No such a layer name: " + name);
}

size_t Coil::get_section_index_by_name(std::string name) const {
    if (!get_sections_description()) {
        throw CoilNotProcessedException("Sections description not set, did you forget to wind?");
    }
    auto sections = get_sections_description().value();
    for (size_t i=0; i<sections.size(); ++i) {
        if (sections[i].get_name() == name) {
            return i;
        }
    }
    throw CoilException(ErrorCode::COIL_WINDING_ERROR, "No such a section name: " + name);
}

bool Coil::are_sections_and_layers_fitting() {
    bool windTurns = true;
    if (!get_sections_description()) {
        return false;
    }
    if (!get_layers_description()) {
        return false;
    }
    auto sections = get_sections_description().value();
    auto layers = get_layers_description().value();

    for (auto& section: sections) {
        if (roundFloat(section.get_filling_factor().value(), 6) > 1 || roundFloat(overlapping_filling_factor(section), 6) > 1 || roundFloat(contiguous_filling_factor(section), 6) > 1) {
            windTurns = false;
        }
    }
    for (auto& layer: layers) {
        if (roundFloat(layer.get_filling_factor().value(), 6) > 1) {
            windTurns = false;
        }
    }

    return windTurns;
}

double Coil::overlapping_filling_factor(Section section) {
    auto bobbin = resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto layers = get_layers_by_section(section.get_name());

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        double sectionWidth = section.get_dimensions()[0];
        double layersWidth = 0;
        for (auto& layer : layers) {
            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                layersWidth += layer.get_dimensions()[0];
            }
            else {
                layersWidth = std::max(layersWidth, layer.get_dimensions()[0]);
            }
        }
        return layersWidth / sectionWidth;
    }
    else {
        double sectionRadialHeight = section.get_dimensions()[0];
        double layersRadialHeight = 0;
        for (auto& layer : layers) {
            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                layersRadialHeight += layer.get_dimensions()[0];
            }
            else {
                layersRadialHeight = std::max(layersRadialHeight, layer.get_dimensions()[0]);
            }
        }
        return layersRadialHeight / sectionRadialHeight;
    }
}

double Coil::contiguous_filling_factor(Section section) {
    auto bobbin = resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto layers = get_layers_by_section(section.get_name());

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        double sectionHeight = section.get_dimensions()[1];
        double layersHeight = 0;
        for (auto& layer : layers) {
            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                layersHeight = std::max(layersHeight, layer.get_dimensions()[1]);
            }
            else {
                layersHeight += layer.get_dimensions()[1];
            }
        }
        return layersHeight / sectionHeight;
    }
    else {
        double sectionAngle = section.get_dimensions()[1];
        double layersAngle = 0;
        for (auto& layer : layers) {
            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                layersAngle = std::max(layersAngle, layer.get_dimensions()[1]);
            }
            else {
                layersAngle += layer.get_dimensions()[1];
            }
        }
        return layersAngle / sectionAngle;

    }
}

std::pair<double, std::pair<double, double>> Coil::calculate_filling_factor(size_t groupIndex) {
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    auto windingOrientation = get_winding_orientation();

    if (!get_layers_description()) {
        throw CoilNotProcessedException("Missing layers to calculate the filling factor.");
    }
    if (!get_turns_description()) {
        throw CoilNotProcessedException("Missing turns to calculate the filling factor.");
    }

    auto layers = get_layers_description().value();
    auto sections = get_sections_description().value();

    double area = 0;
    double availableArea = windingWindows[0].get_area().value();
    double availableContiguousDimension;
    double availableOverlappingDimension;
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        availableContiguousDimension = windingWindows[0].get_height().value();
        availableOverlappingDimension = windingWindows[0].get_width().value();
    }
    else {
        availableOverlappingDimension = windingWindows[0].get_radial_height().value();
        availableContiguousDimension = windingWindows[0].get_angle().value();
    }
    double maximumLayerFillingFactor = 0;
    double contiguousDimension = 0;
    double overlappingDimension = 0;

    for (auto section : sections) {
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            if (section.get_type() == ElectricalType::CONDUCTION) {
                contiguousDimension = std::max(contiguousDimension, section.get_dimensions()[1]);
            }
            overlappingDimension += section.get_dimensions()[0];
        }
        else {
            if (section.get_type() == ElectricalType::CONDUCTION) {
                overlappingDimension = std::max(overlappingDimension, section.get_dimensions()[0]);
            }
            contiguousDimension += section.get_dimensions()[1];

        }
    }

    for (auto section : sections) {
        if (section.get_margin()) {
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                area += (resolve_margin(section)[0] + resolve_margin(section)[1]) * section.get_dimensions()[0];
            }
            else {
                area += (resolve_margin(section)[0] + resolve_margin(section)[1]) * section.get_dimensions()[1];
            }
        }
    }

    for (auto layer : layers) {
        if (layer.get_filling_factor()) {
            if (layer.get_filling_factor().value() > 1) {
                maximumLayerFillingFactor = std::max(maximumLayerFillingFactor, layer.get_filling_factor().value());
            }
        }
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            auto turns = get_turns_by_layer(layer.get_name());
            for (auto turn : turns) {
                area += turn.get_dimensions().value()[0] * turn.get_dimensions().value()[1];
            }
        }
        else {
            if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                area += layer.get_dimensions()[0] * layer.get_dimensions()[1];
            }
            else {
                auto arc = angle_to_wound_distance(layer.get_dimensions()[1], (availableOverlappingDimension - layer.get_coordinates()[0]));
                area += layer.get_dimensions()[0] * arc;
            }
        }
    }

    double areaFillingFactor = area / availableArea;
    areaFillingFactor = std::max(maximumLayerFillingFactor, areaFillingFactor);
    double contiguousFillingFactor = contiguousDimension / availableContiguousDimension;
    double overlappingFillingFactor = overlappingDimension / availableOverlappingDimension;
    return {areaFillingFactor, {overlappingFillingFactor, contiguousFillingFactor}};
}

std::pair<uint64_t, std::vector<double>> get_parallels_proportions(size_t slotIndex, size_t slots, uint64_t numberTurns, uint64_t numberParallels, 
                                                                   std::vector<double> remainingParallelsProportion, WindingStyle windByConsecutiveTurns,
                                                                   std::vector<double> totalParallelsProportion, double slotRelativeProportion=1,
                                                                   std::optional<double> slotAbsolutePhysicalTurns = std::nullopt) {
    uint64_t physicalTurnsThisSlot = 0;
    std::vector<double> slotParallelsProportion(numberParallels, 0);
    if (windByConsecutiveTurns == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
        size_t remainingPhysicalTurns = 0;
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            remainingPhysicalTurns += round(remainingParallelsProportion[parallelIndex] * numberTurns);
        }
        if (slotAbsolutePhysicalTurns)
            physicalTurnsThisSlot = slotAbsolutePhysicalTurns.value();
        else
            physicalTurnsThisSlot = std::min(uint64_t(remainingPhysicalTurns), uint64_t(ceil(double(remainingPhysicalTurns) / (slots - slotIndex) * slotRelativeProportion)));
        uint64_t remainingPhysicalTurnsThisSection = physicalTurnsThisSlot;

        size_t currentParallel = 0;
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            if (remainingParallelsProportion[parallelIndex] > 0) {
                currentParallel = parallelIndex;
                break;
            }
        }

        while (remainingPhysicalTurnsThisSection > 0) {
            uint64_t numberTurnsToFitInCurrentParallel = round(remainingParallelsProportion[currentParallel] * numberTurns);
            if (remainingPhysicalTurnsThisSection >= numberTurnsToFitInCurrentParallel) {
                remainingPhysicalTurnsThisSection -= numberTurnsToFitInCurrentParallel;
                slotParallelsProportion[currentParallel] = double(numberTurnsToFitInCurrentParallel) / numberTurns;
                currentParallel++;
            }
            else {
                double proportionParallelsThisSection = double(remainingPhysicalTurnsThisSection) / numberTurns;
                slotParallelsProportion[currentParallel] += proportionParallelsThisSection;
                remainingPhysicalTurnsThisSection = 0;
            }
        }
    }
    else if (slotAbsolutePhysicalTurns) {
        // Real-winding blocking forces an exact physical-turn count for this slot (layer). Distribute
        // those turns across the parallels round-robin — exactly how WIND_BY_CONSECUTIVE_PARALLELS lays
        // them side by side (P0,P1,...,P0,P1) — so each parallel keeps an equal share and a blocked
        // layer's forced split survives (the default branch below ignores the forced count and reverts
        // to an even per-slot division). The caller passes a multiple of the active-parallel count, so
        // the share comes out equal; round-robin still degrades gracefully if it does not.
        uint64_t turnsToPlace = uint64_t(slotAbsolutePhysicalTurns.value());
        std::vector<uint64_t> remainingPerParallel(numberParallels, 0);
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            remainingPerParallel[parallelIndex] = uint64_t(std::round(remainingParallelsProportion[parallelIndex] * numberTurns));
        }
        std::vector<uint64_t> placePerParallel(numberParallels, 0);
        uint64_t placed = 0;
        bool progress = true;
        while (placed < turnsToPlace && progress) {
            progress = false;
            for (size_t parallelIndex = 0; parallelIndex < numberParallels && placed < turnsToPlace; ++parallelIndex) {
                if (placePerParallel[parallelIndex] < remainingPerParallel[parallelIndex]) {
                    placePerParallel[parallelIndex]++;
                    placed++;
                    progress = true;
                }
            }
        }
        physicalTurnsThisSlot = placed;
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            slotParallelsProportion[parallelIndex] = double(placePerParallel[parallelIndex]) / numberTurns;
        }
    }
    else {
        for (size_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            double remainingSlots = slots - slotIndex;
            double remainingTurnsBeforeThisParallel = numberTurns * remainingParallelsProportion[parallelIndex];
            double numberTurnsToAddToCurrentParallel = ceil(roundFloat(remainingTurnsBeforeThisParallel / remainingSlots * slotRelativeProportion, 10));
            // double numberTurnsToAddToCurrentParallel = ceil(numberTurns * totalParallelsProportion[parallelIndex] / slots);
            double remainingTurnsAfterThisParallel = remainingTurnsBeforeThisParallel - numberTurnsToAddToCurrentParallel;
            double remainingSlotsAfterThisOne = remainingSlots - 1;
            if (remainingTurnsAfterThisParallel < remainingSlotsAfterThisOne) {
                numberTurnsToAddToCurrentParallel = ceil(roundFloat(remainingTurnsBeforeThisParallel / remainingSlots, 10));
            }
            double proportionParallelsThisSection = std::min(remainingParallelsProportion[parallelIndex], numberTurnsToAddToCurrentParallel / numberTurns);
            physicalTurnsThisSlot += numberTurnsToAddToCurrentParallel;
            slotParallelsProportion[parallelIndex] = proportionParallelsThisSection;
        }
    }

    return {physicalTurnsThisSlot, slotParallelsProportion};
}

double get_area_used_in_wires(OpenMagnetics::Wire wire, uint64_t physicalTurns) {
    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        double wireDiameter = wire.get_maximum_outer_width();
        return physicalTurns * pow(wireDiameter, 2);
    }
    else {
        double wireWidth = wire.get_maximum_outer_width();
        double wireHeight = wire.get_maximum_outer_height();
        return physicalTurns * wireWidth * wireHeight;
    }
}

void Coil::set_insulation_layers(std::map<std::pair<size_t, size_t>, std::vector<Layer>> insulationLayers) {
    _insulationInterSectionsLayers = insulationLayers;
}

bool Coil::calculate_custom_thickness_insulation(double thickness) {
    // Insulation layers just for mechanical reasons, one layer between sections at least
    auto wirePerWinding = get_wires();

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    auto layersOrientation = _layersOrientation;

    // TODO: Properly think about insulation layers with weird windings
    auto windingOrientation = get_winding_orientation();

    if (windingOrientation == WindingOrientation::CONTIGUOUS && _layersOrientation == WindingOrientation::OVERLAPPING) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::CONTIGUOUS;
        }
    }
    if (windingOrientation == WindingOrientation::OVERLAPPING && _layersOrientation == WindingOrientation::CONTIGUOUS) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::OVERLAPPING;
        }
    }

    for (size_t leftTopWindingIndex = 0; leftTopWindingIndex < get_functional_description().size(); ++leftTopWindingIndex) {
        for (size_t rightBottomWindingIndex = 0; rightBottomWindingIndex < get_functional_description().size(); ++rightBottomWindingIndex) {
            // if (leftTopWindingIndex == rightBottomWindingIndex) {
            //     continue;
            // }
            auto wireLeftTopWinding = wirePerWinding[leftTopWindingIndex];
            auto wireRightBottomWinding = wirePerWinding[rightBottomWindingIndex];
            auto windingsMapKey = std::pair<size_t, size_t>{leftTopWindingIndex, rightBottomWindingIndex};

            CoilSectionInterface coilSectionInterface;
            coilSectionInterface.set_number_layers_insulation(1);
            InsulationMaterial defaultInsulationMaterial = find_insulation_material_by_name(defaults.defaultInsulationMaterial);
            coilSectionInterface.set_solid_insulation_thickness(defaultInsulationMaterial.get_thinner_tape_thickness());
            coilSectionInterface.set_total_margin_tape_distance(0);
            coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);

            _insulationInterSectionsLayers[windingsMapKey] = std::vector<Layer>();
            _coilSectionInterfaces[windingsMapKey] = coilSectionInterface;

            Layer layer;
            layer.set_partial_windings(std::vector<PartialWinding>{});
            // layer.set_section(section.get_name());
            layer.set_type(ElectricalType::INSULATION);
            layer.set_name("temp");
            layer.set_orientation(layersOrientation);
            layer.set_turns_alignment(CoilAlignment::SPREAD); // HARDCODED, maybe in the future configure for shields made of turns?

            if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                double windingWindowHeight = windingWindows[0].get_height().value();
                double windingWindowWidth = windingWindows[0].get_width().value();
                if (layersOrientation == WindingOrientation::OVERLAPPING) {
                    layer.set_dimensions(std::vector<double>{thickness, windingWindowHeight});
                }
                else if (layersOrientation == WindingOrientation::CONTIGUOUS) {
                    layer.set_dimensions(std::vector<double>{windingWindowWidth, thickness});
                }
            }
            else {
                layer.set_coordinate_system(CoordinateSystem::POLAR);
                double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                double windingWindowAngle = windingWindows[0].get_angle().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    layer.set_dimensions(std::vector<double>{thickness, windingWindowAngle});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    double tapeThicknessInAngle = wound_distance_to_angle(thickness, windingWindowRadialHeight);
                    layer.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                }
            }
            // layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
            layer.set_filling_factor(1);
            // Custom-thickness mechanical-only insulation: no material was
            // chosen above (this path only knows the requested thickness),
            // so fall back to the default layer insulation material so
            // downstream consumers (Temperature, StrayCapacitance) can
            // read thermal_conductivity / permittivity from the MAS.
            layer.set_insulation_material(defaults.defaultLayerInsulationMaterial);
            _insulationInterSectionsLayers[windingsMapKey].push_back(layer);

            Section section;
            section.set_name("temp");
            section.set_partial_windings(std::vector<PartialWinding>{});
            section.set_layers_orientation(layersOrientation);
            section.set_type(ElectricalType::INSULATION);

            if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                section.set_coordinate_system(CoordinateSystem::CARTESIAN);
                double windingWindowHeight = windingWindows[0].get_height().value();
                double windingWindowWidth = windingWindows[0].get_width().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{thickness, windingWindowHeight});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    section.set_dimensions(std::vector<double>{windingWindowWidth, thickness});
                }
            }
            else {
                section.set_coordinate_system(CoordinateSystem::POLAR);
                double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                double windingWindowAngle = windingWindows[0].get_angle().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{thickness, windingWindowAngle});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    double tapeThicknessInAngle = wound_distance_to_angle(thickness, windingWindowRadialHeight);
                    section.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                }
            }
            // section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(1);
            _insulationSections[windingsMapKey] = section;
        }
    }
    return true;
}

bool Coil::calculate_mechanical_insulation() {
    // Insulation layers just for mechanical reasons, one layer between sections at least
    auto wirePerWinding = get_wires();

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    auto layersOrientation = _layersOrientation;

    // TODO: Properly think about insulation layers with weird windings
    auto windingOrientation = get_winding_orientation();

    if (windingOrientation == WindingOrientation::CONTIGUOUS && _layersOrientation == WindingOrientation::OVERLAPPING) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::CONTIGUOUS;
        }
    }
    if (windingOrientation == WindingOrientation::OVERLAPPING && _layersOrientation == WindingOrientation::CONTIGUOUS) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::OVERLAPPING;
        }
    }

    for (size_t leftTopWindingIndex = 0; leftTopWindingIndex < get_functional_description().size(); ++leftTopWindingIndex) {
        for (size_t rightBottomWindingIndex = 0; rightBottomWindingIndex < get_functional_description().size(); ++rightBottomWindingIndex) {
            if (leftTopWindingIndex == rightBottomWindingIndex) {
                continue;
            }
            auto wireLeftTopWinding = wirePerWinding[leftTopWindingIndex];
            auto wireRightBottomWinding = wirePerWinding[rightBottomWindingIndex];
            auto windingsMapKey = std::pair<size_t, size_t>{leftTopWindingIndex, rightBottomWindingIndex};

            CoilSectionInterface coilSectionInterface;
            coilSectionInterface.set_number_layers_insulation(1);
            InsulationMaterial defaultInsulationMaterial = find_insulation_material_by_name(defaults.defaultInsulationMaterial);
            coilSectionInterface.set_solid_insulation_thickness(defaultInsulationMaterial.get_thinner_tape_thickness());
            coilSectionInterface.set_total_margin_tape_distance(0);
            coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);

            _insulationInterSectionsLayers[windingsMapKey] = std::vector<Layer>();
            _coilSectionInterfaces[windingsMapKey] = coilSectionInterface;

            for (size_t layerIndex = 0; layerIndex < coilSectionInterface.get_number_layers_insulation(); ++layerIndex) {
                Layer layer;
                layer.set_partial_windings(std::vector<PartialWinding>{});
                // layer.set_section(section.get_name());
                layer.set_type(ElectricalType::INSULATION);
                layer.set_name("temp");
                layer.set_orientation(layersOrientation);
                layer.set_turns_alignment(CoilAlignment::SPREAD); // HARDCODED, maybe in the future configure for shields made of turns?

                if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                    layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                    double windingWindowHeight = windingWindows[0].get_height().value();
                    double windingWindowWidth = windingWindows[0].get_width().value();
                    if (layersOrientation == WindingOrientation::OVERLAPPING) {
                        layer.set_dimensions(std::vector<double>{defaultInsulationMaterial.get_thinner_tape_thickness(), windingWindowHeight});
                    }
                    else if (layersOrientation == WindingOrientation::CONTIGUOUS) {
                        layer.set_dimensions(std::vector<double>{windingWindowWidth, defaultInsulationMaterial.get_thinner_tape_thickness()});
                    }
                }
                else {
                    layer.set_coordinate_system(CoordinateSystem::POLAR);
                    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                    double windingWindowAngle = windingWindows[0].get_angle().value();
                    if (windingOrientation == WindingOrientation::OVERLAPPING) {
                        layer.set_dimensions(std::vector<double>{defaultInsulationMaterial.get_thinner_tape_thickness(), windingWindowAngle});
                    }
                    else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                        double tapeThicknessInAngle = wound_distance_to_angle(defaultInsulationMaterial.get_thinner_tape_thickness(), windingWindowRadialHeight);
                        layer.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                    }
                }
                // layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                layer.set_filling_factor(1);
                // Default-material mechanical insulation: propagate the
                // chosen default so Temperature/StrayCapacitance can read
                // thermal_conductivity / permittivity from the MAS.
                layer.set_insulation_material(static_cast<MAS::InsulationMaterial>(defaultInsulationMaterial));
                _insulationInterSectionsLayers[windingsMapKey].push_back(layer);
            }
            // _insulationInterSectionsLayersLog[windingsMapKey] = "Adding " + std::to_string(coilSectionInterface.get_number_layers_insulation()) + " insulation layers, as we need a thickness of " + std::to_string(smallestInsulationThicknessCoveringRemaining * 1000) + " mm to achieve " + neededInsulationTypeString + " insulation";

            Section section;
            section.set_name("temp");
            section.set_partial_windings(std::vector<PartialWinding>{});
            section.set_layers_orientation(layersOrientation);
            section.set_type(ElectricalType::INSULATION);

            if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                section.set_coordinate_system(CoordinateSystem::CARTESIAN);
                double windingWindowHeight = windingWindows[0].get_height().value();
                double windingWindowWidth = windingWindows[0].get_width().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowHeight});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    section.set_dimensions(std::vector<double>{windingWindowWidth, coilSectionInterface.get_solid_insulation_thickness()});
                }
            }
            else {
                section.set_coordinate_system(CoordinateSystem::POLAR);
                double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                double windingWindowAngle = windingWindows[0].get_angle().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowAngle});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    double tapeThicknessInAngle = wound_distance_to_angle(coilSectionInterface.get_solid_insulation_thickness(), windingWindowRadialHeight);
                    section.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                }
            }
            // section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(1);
            _insulationSections[windingsMapKey] = section;
        }
    }
    return true;
}

bool Coil::calculate_insulation(bool simpleMode) {
    auto inputs = _inputs.value();

    if (!inputs.get_design_requirements().get_insulation()) {
        return false;
    }

    auto wirePerWinding = get_wires();

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto layersOrientation = _layersOrientation;
    auto windingOrientation = get_winding_orientation();

    for (size_t leftTopWindingIndex = 0; leftTopWindingIndex < get_functional_description().size(); ++leftTopWindingIndex) {
        for (size_t rightBottomWindingIndex = 0; rightBottomWindingIndex < get_functional_description().size(); ++rightBottomWindingIndex) {
            if (leftTopWindingIndex == rightBottomWindingIndex) {
                continue;
            }
            auto wireLeftTopWinding = wirePerWinding[leftTopWindingIndex];
            auto wireRightBottomWinding = wirePerWinding[rightBottomWindingIndex];
            auto windingsMapKey = std::pair<size_t, size_t>{leftTopWindingIndex, rightBottomWindingIndex};

            CoilSectionInterface coilSectionInterface;
            coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::INSULATING);
            InsulationMaterial chosenInsulationMaterial;

            if (simpleMode) {
                InsulationMaterial defaultInsulationMaterial = find_insulation_material_by_name(defaults.defaultInsulationMaterial);
                chosenInsulationMaterial = defaultInsulationMaterial;
                coilSectionInterface.set_solid_insulation_thickness(defaultInsulationMaterial.get_thinner_tape_thickness());
                if (settings.get_coil_allow_margin_tape()) {
                    coilSectionInterface.set_number_layers_insulation(1);
                    coilSectionInterface.set_total_margin_tape_distance(_standardCoordinator.calculate_creepage_distance(inputs, true));
                }
                else {
                    coilSectionInterface.set_number_layers_insulation(3);
                    coilSectionInterface.set_total_margin_tape_distance(0);
                }
            }
            else {
                coilSectionInterface.set_solid_insulation_thickness(DBL_MAX);
                coilSectionInterface.set_number_layers_insulation(ULONG_MAX);

                if (insulationMaterialDatabase.empty()) {
                    load_insulation_materials();
                }

                for (auto& insulationMaterial : insulationMaterialDatabase) {
                    auto auxCoilSectionInterface = _standardCoordinator.calculate_coil_section_interface_layers(inputs, wireLeftTopWinding, wireRightBottomWinding, insulationMaterial.second);
                    if (auxCoilSectionInterface) {
                        if (auxCoilSectionInterface.value().get_solid_insulation_thickness() < coilSectionInterface.get_solid_insulation_thickness()) {
                            coilSectionInterface = auxCoilSectionInterface.value();
                            chosenInsulationMaterial = insulationMaterial.second;
                        }
                    }
                }

                if (coilSectionInterface.get_number_layers_insulation() == ULONG_MAX) {
                    throw InvalidInputException("No insulation found with current requirements");
                }
            }

            _insulationInterSectionsLayers[windingsMapKey] = std::vector<Layer>();
            _coilSectionInterfaces[windingsMapKey] = coilSectionInterface;

            for (size_t layerIndex = 0; layerIndex < coilSectionInterface.get_number_layers_insulation(); ++layerIndex) {
                Layer layer;
                layer.set_partial_windings(std::vector<PartialWinding>{});
                // layer.set_section(section.get_name());
                layer.set_type(ElectricalType::INSULATION);
                layer.set_name("temp");
                layer.set_orientation(_layersOrientation);
                layer.set_turns_alignment(CoilAlignment::SPREAD); // HARDCODED, maybe in the future configure for shields made of turns?

                if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                    layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                    double windingWindowHeight = windingWindows[0].get_height().value();
                    double windingWindowWidth = windingWindows[0].get_width().value();
                    if (layersOrientation == WindingOrientation::OVERLAPPING) {
                        layer.set_dimensions(std::vector<double>{chosenInsulationMaterial.get_thinner_tape_thickness(), windingWindowHeight});
                    }
                    else if (layersOrientation == WindingOrientation::CONTIGUOUS) {
                        layer.set_dimensions(std::vector<double>{windingWindowWidth, chosenInsulationMaterial.get_thinner_tape_thickness()});
                    }
                }
                else {
                    layer.set_coordinate_system(CoordinateSystem::POLAR);
                    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                    double windingWindowAngle = windingWindows[0].get_angle().value();
                    if (windingOrientation == WindingOrientation::OVERLAPPING) {
                        layer.set_dimensions(std::vector<double>{chosenInsulationMaterial.get_thinner_tape_thickness(), windingWindowAngle});
                    }
                    else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                        double tapeThicknessInAngle = wound_distance_to_angle(chosenInsulationMaterial.get_thinner_tape_thickness(), windingWindowRadialHeight);
                        layer.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                    }
                }

                layer.set_filling_factor(1);
                // Propagate the insulation material chosen above so downstream
                // consumers (Temperature::getInsulationLayerThermalResistance,
                // StrayCapacitance) have the dielectric/thermal properties
                // they need. Without this the Temperature plot throws once
                // it tries to read the layer's thermal_conductivity.
                layer.set_insulation_material(static_cast<MAS::InsulationMaterial>(chosenInsulationMaterial));
                _insulationInterSectionsLayers[windingsMapKey].push_back(layer);
            }
            // _insulationInterSectionsLayersLog[windingsMapKey] = "Adding " + std::to_string(coilSectionInterface.get_number_layers_insulation()) + " insulation layers, as we need a thickness of " + std::to_string(smallestInsulationThicknessCoveringRemaining * 1000) + " mm to achieve " + neededInsulationTypeString + " insulation";

            Section section;
            section.set_name("temp");
            section.set_partial_windings(std::vector<PartialWinding>{});
            section.set_layers_orientation(_layersOrientation);
            section.set_type(ElectricalType::INSULATION);

            if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
                section.set_coordinate_system(CoordinateSystem::CARTESIAN);
                double windingWindowHeight = windingWindows[0].get_height().value();
                double windingWindowWidth = windingWindows[0].get_width().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowHeight});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    section.set_dimensions(std::vector<double>{windingWindowWidth, coilSectionInterface.get_solid_insulation_thickness()});
                }
            }
            else {
                section.set_coordinate_system(CoordinateSystem::POLAR);
                double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
                double windingWindowAngle = windingWindows[0].get_angle().value();
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{coilSectionInterface.get_solid_insulation_thickness(), windingWindowAngle});
                }
                else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                    double tapeThicknessInAngle = wound_distance_to_angle(coilSectionInterface.get_solid_insulation_thickness(), windingWindowRadialHeight);
                    section.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
                }
            }
            // section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            section.set_filling_factor(1);
            _insulationSections[windingsMapKey] = section;
        }
    }
    return true;
}

std::vector<std::pair<size_t, double>> Coil::get_ordered_sections(double spaceForSections, std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    std::vector<std::pair<size_t, double>> orderedSections;
    double numberWindings = get_functional_description().size();
    auto numberSectionsPerWinding = std::vector<size_t>(numberWindings, 0);
    for (auto windingIndex : pattern) {
        if (windingIndex >= numberWindings) {
            throw std::invalid_argument("Winding index does not exist in winding");
        }
        numberSectionsPerWinding[windingIndex] += repetitions;
    }

    for (size_t repetitionIndex = 0; repetitionIndex < repetitions; ++repetitionIndex) {
        for (auto windingIndex : pattern) {
            if (roundFloat(proportionPerWinding[windingIndex], 6) > 1) {
                throw std::invalid_argument("proportionPerWinding[windingIndex] cannot be greater than 1: " + std::to_string(proportionPerWinding[windingIndex]));
            }
            double spaceForSection = roundFloat(spaceForSections * proportionPerWinding[windingIndex] / numberSectionsPerWinding[windingIndex], 9);
            orderedSections.push_back({windingIndex, spaceForSection});
        }
    }

    return orderedSections;
}

std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> Coil::add_insulation_to_sections(std::vector<std::pair<size_t, double>> orderedSections){
    std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation;
    auto windingOrientation = get_winding_orientation();
    for (size_t sectionIndex = 1; sectionIndex < orderedSections.size(); ++sectionIndex) {
        auto leftWindingIndex = orderedSections[sectionIndex - 1].first;
        auto rightWindingIndex = orderedSections[sectionIndex].first;
        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex}; 
        if (!_insulationSections.contains(windingsMapKey)) {
            continue;
        }
        auto currentSpaceForLeftSection = orderedSections[sectionIndex - 1].second;
        auto currentSpaceForRightSection = orderedSections[sectionIndex].second;

        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            std::pair<size_t, double> leftSectionInfo = {leftWindingIndex, currentSpaceForLeftSection - _insulationSections[windingsMapKey].get_dimensions()[0] / 2};
            orderedSections[sectionIndex - 1] = leftSectionInfo;
            std::pair<size_t, double> rightSectionInfo = {rightWindingIndex, currentSpaceForRightSection - _insulationSections[windingsMapKey].get_dimensions()[0] / 2};
            orderedSections[sectionIndex] = rightSectionInfo;
        }
        else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            std::pair<size_t, double> leftSectionInfo = {leftWindingIndex, currentSpaceForLeftSection - _insulationSections[windingsMapKey].get_dimensions()[1] / 2};
            orderedSections[sectionIndex - 1] = leftSectionInfo;
            std::pair<size_t, double> rightSectionInfo = {rightWindingIndex, currentSpaceForRightSection - _insulationSections[windingsMapKey].get_dimensions()[1] / 2};
            orderedSections[sectionIndex] = rightSectionInfo;
        }
    }

    orderedSectionsWithInsulation.push_back({ElectricalType::CONDUCTION, orderedSections[0]});
    for (size_t sectionIndex = 1; sectionIndex < orderedSections.size(); ++sectionIndex) {
        auto leftWindingIndex = orderedSections[sectionIndex - 1].first;
        auto rightWindingIndex = orderedSections[sectionIndex].first;
        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex}; 
        if (_insulationSections.contains(windingsMapKey)) {
            std::pair<size_t, double> insulationSectionInfo;
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[0]};
            }
            else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[1]};
            }

            orderedSectionsWithInsulation.push_back({ElectricalType::INSULATION, insulationSectionInfo});
        }
        orderedSectionsWithInsulation.push_back({ElectricalType::CONDUCTION, orderedSections[sectionIndex]});
    }

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();


    // last insulation layer we compare between last and first
    if (windingOrientation != WindingOrientation::CONTIGUOUS || bobbinWindingWindowShape != WindingWindowShape::RECTANGULAR) {
        // We don't add one in the sections are contiguous, as they end in the bobbin
        auto leftWindingIndex = orderedSections.back().first;
        auto rightWindingIndex = orderedSections[0].first;
        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex}; 

        if (_insulationSections.contains(windingsMapKey)) {
            std::pair<size_t, double> insulationSectionInfo;
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[0]};
            }
            else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                insulationSectionInfo = {SIZE_MAX, _insulationSections[windingsMapKey].get_dimensions()[1]};
            }

            orderedSectionsWithInsulation.push_back({ElectricalType::INSULATION, insulationSectionInfo});
        }
    }

    return orderedSectionsWithInsulation;
}

std::vector<double> Coil::get_proportion_per_winding_based_on_wires() {
    std::vector<double> physicalTurnsAreaPerWinding;
    double totalPhysicalTurnsArea = 0;
    auto wirePerWinding = get_wires();
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex){
        double wireOuterRectangularArea = wirePerWinding[windingIndex].get_maximum_outer_width() * wirePerWinding[windingIndex].get_maximum_outer_height();
        double totalAreaThisWinding = wireOuterRectangularArea * get_functional_description()[windingIndex].get_number_turns() * get_functional_description()[windingIndex].get_number_parallels();
        physicalTurnsAreaPerWinding.push_back(totalAreaThisWinding);
        totalPhysicalTurnsArea += totalAreaThisWinding;
    }
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex){
        physicalTurnsAreaPerWinding[windingIndex] /= totalPhysicalTurnsArea;
    }

    return physicalTurnsAreaPerWinding;
}

std::pair<double, double> get_section_round_dimensions(std::pair<ElectricalType, std::pair<size_t, double>> sectionWithInsulationScaledWithArea,
                                                 WindingOrientation windingOrientation, double windingWindowRadialHeight, double windingWindowAngle) {

    auto sectionInfo = sectionWithInsulationScaledWithArea.second;
    auto spaceForSection = sectionInfo.second;

    double currentSectionRadialHeight = 0;
    double currentSectionAngle = 0;
    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        currentSectionRadialHeight = spaceForSection;
        currentSectionAngle = windingWindowAngle;
    }
    else {
        currentSectionRadialHeight = windingWindowRadialHeight;
        currentSectionAngle = spaceForSection;
    }

    return {currentSectionRadialHeight, currentSectionAngle};
}

std::vector<double> get_physical_turns_proportions(std::vector<int64_t> physicalTurns) {
    std::vector<double> physicalTurnsProportions;
    double average = 0;
    for (size_t index = 0; index < physicalTurns.size(); ++index) {
        average += double(physicalTurns[index]);
    }
    average /= physicalTurns.size();

    for (size_t index = 0; index < physicalTurns.size(); ++index) {
        if (index + 1 < physicalTurns.size())
            physicalTurnsProportions.push_back(double(physicalTurns[index]) / average);
        else
            physicalTurnsProportions.push_back(1 + double(physicalTurns[index]) / average);
    }

    return physicalTurnsProportions;
}

std::vector<double> get_length_proportions(std::vector<double> lengths, std::vector<size_t> windingIndexes) {
    // Map each distinct RAW winding index to a compact 0-based position. The
    // averages/numberSectionsPerWinding arrays are sized by the number of
    // distinct windings, so they must be indexed by the compact position, NOT by
    // the raw winding index — otherwise a non-contiguous pattern (e.g. {0,2})
    // writes out of bounds (heap corruption) on toroidal/round windings.
    std::vector<size_t> uniqueIndexes;
    std::map<size_t, size_t> compactPosition;
    for (size_t windingIndex = 0; windingIndex < windingIndexes.size(); ++windingIndex) {
        if (compactPosition.find(windingIndexes[windingIndex]) == compactPosition.end()) {
            compactPosition[windingIndexes[windingIndex]] = uniqueIndexes.size();
            uniqueIndexes.push_back(windingIndexes[windingIndex]);
        }
    }

    std::vector<double> lengthProportions;
    std::vector<double> averages(uniqueIndexes.size(), 0);
    std::vector<double> numberSectionsPerWinding(uniqueIndexes.size(), 0);

    for (size_t index = 0; index < lengths.size(); ++index) {
        size_t pos = compactPosition[windingIndexes[index]];
        averages[pos] += lengths[index];
        numberSectionsPerWinding[pos]++;
    }

    for (size_t windingIndex = 0; windingIndex < averages.size(); ++windingIndex) {
        averages[windingIndex] /= numberSectionsPerWinding[windingIndex];
    }

    for (size_t index = 0; index < lengths.size(); ++index) {
        size_t pos = compactPosition[windingIndexes[index]];
        if (index + 1 < lengths.size())
            lengthProportions.push_back(lengths[index] / averages[pos]);
        else
            lengthProportions.push_back(1 + lengths[index] / averages[pos]);
    }

    return lengthProportions;
}

std::vector<double> get_section_lengths(std::vector<double> currentSectionRadialHeights, std::vector<double> currentSectionAngles, double windingWindowRadialHeight) {
    std::vector<double> sectionLengths;
    double radialHeightIncrease = windingWindowRadialHeight / currentSectionRadialHeights.size();
    for (size_t sectionIndex = 0; sectionIndex < currentSectionRadialHeights.size(); ++sectionIndex) {
        double radius = windingWindowRadialHeight - radialHeightIncrease * sectionIndex - radialHeightIncrease;
        sectionLengths.push_back(2 * std::numbers::pi * radius * currentSectionAngles[sectionIndex] / 360);
    }
    return sectionLengths;
}

std::vector<double> get_section_areas(std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulationScaledWithArea, std::vector<double> currentSectionAngles, double windingWindowRadialHeight) {
    std::vector<double> sectionAreas;
    double currentRadius = windingWindowRadialHeight;
    size_t currentConductionSectionIndex = 0;
    for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulationScaledWithArea.size(); ++sectionIndex) {
        if (orderedSectionsWithInsulationScaledWithArea[sectionIndex].first == ElectricalType::CONDUCTION) {
            auto sectionInfo = orderedSectionsWithInsulationScaledWithArea[sectionIndex].second;
            auto spaceForSection = sectionInfo.second;
            double outerRadius = currentRadius;
            double innerRadius = currentRadius - spaceForSection;
            currentRadius -= spaceForSection;
            sectionAreas.push_back(std::numbers::pi * (pow(outerRadius, 2) - pow(innerRadius, 2)) * currentSectionAngles[currentConductionSectionIndex] / 360);
            currentConductionSectionIndex++;
        }

    }
    return sectionAreas;
}

std::pair<size_t, std::vector<int64_t>> get_number_layers_needed_and_number_physical_turns(double radialHeight, double angle, Wire wire, int64_t physicalTurnsInSection, double windingWindowRadius) {
    int64_t reaminingPhysicalTurnsInSection = physicalTurnsInSection;
    double wireWidth = resolve_dimensional_values(wire.get_maximum_outer_width());
    double wireHeight = resolve_dimensional_values(wire.get_maximum_outer_height());
    double currentRadialHeight = radialHeight;
    double currentRadius;
    if (wire.get_type() == WireType::FOIL){
        throw NotImplementedException("Foil is not supported in toroids");
    }
    if (wire.get_type() == WireType::PLANAR){
        throw NotImplementedException("Planar is not supported in toroids");
    }
    if (wire.get_type() == WireType::RECTANGULAR){
        currentRadius = windingWindowRadius - wireWidth - currentRadialHeight;
    }
    else {
        currentRadius = windingWindowRadius - wireWidth / 2 - currentRadialHeight;
    }
    double sectionAvailableAngle = angle;
    std::vector<int64_t> layerPhysicalTurns;
    size_t numberLayers = 0;
    while (reaminingPhysicalTurnsInSection > 0) {
        double wireAngle = wound_distance_to_angle(wireHeight, std::max(wireWidth, currentRadius));
        int64_t numberTurnsFittingThisLayer = std::max(1.0, floor(sectionAvailableAngle / wireAngle));
        reaminingPhysicalTurnsInSection -= numberTurnsFittingThisLayer;

        layerPhysicalTurns.push_back(numberTurnsFittingThisLayer);
        numberLayers++;
        if (currentRadius > wireWidth) {
            currentRadius -= wireWidth;
        }
    }

    int64_t numberTurnsToCorrect = -reaminingPhysicalTurnsInSection;
    size_t currentIndex = numberLayers - 1;
    while (numberTurnsToCorrect > 0) {
        layerPhysicalTurns[currentIndex]--;
        numberTurnsToCorrect--;
        if (currentIndex == 0)
            currentIndex = numberLayers - 1;
        else
            currentIndex--;
    }

    return {numberLayers, layerPhysicalTurns};
}

std::pair<size_t, std::vector<int64_t>> get_number_layers_needed_and_number_physical_turns(Section section, Wire wire, int64_t physicalTurnsInSection, double windingWindowRadius) {
    return get_number_layers_needed_and_number_physical_turns(section.get_coordinates()[0] - section.get_dimensions()[0] / 2, section.get_dimensions()[1], wire, physicalTurnsInSection, windingWindowRadius);
}

void Coil::apply_margin_tape(std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation, size_t sectionIndexOffset) {
    if (_marginsPerSection.size() < sectionIndexOffset + orderedSectionsWithInsulation.size()) {
        // Resize (not replace) so preloaded margins are preserved, matching equalize_margins
        _marginsPerSection.resize(sectionIndexOffset + orderedSectionsWithInsulation.size(), {0, 0});
    }

    for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
        size_t marginIndex = sectionIndexOffset + sectionIndex;
        if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
            if (sectionIndex > 0 && !_coilSectionInterfaces.empty()) {

                if (orderedSectionsWithInsulation[sectionIndex - 1].first != ElectricalType::INSULATION) {
                    // Adjacent conduction sections are allowed when the pair
                    // belongs to the same isolation side (no insulation was
                    // inserted by add_insulation_to_sections) — this happens
                    // for wound_with-grouped center-tap halves, and for
                    // wires within the same winding side generally.
                    auto leftIdx = orderedSectionsWithInsulation[sectionIndex].second.first;
                    auto rightIdxPrev = orderedSectionsWithInsulation[sectionIndex - 1].second.first;
                    auto pairKey = std::pair<size_t, size_t>{rightIdxPrev, leftIdx};
                    if (!_insulationSections.contains(pairKey)) {
                        continue;
                    }
                    throw InvalidInputException("There cannot be two sections without insulation in between");
                }
                auto windingIndex = orderedSectionsWithInsulation[sectionIndex].second.first;
                auto previousWindingIndex = orderedSectionsWithInsulation[sectionIndex - 2].second.first;
                auto windingsMapKey = std::pair<size_t, size_t>{previousWindingIndex, windingIndex};
                auto coilSectionInterface = _coilSectionInterfaces[windingsMapKey];
                _marginsPerSection[marginIndex][0] =  std::max(_marginsPerSection[marginIndex][0], coilSectionInterface.get_total_margin_tape_distance() / 2);
                _marginsPerSection[marginIndex][1] =  std::max(_marginsPerSection[marginIndex][1], coilSectionInterface.get_total_margin_tape_distance() / 2);
                _marginsPerSection[marginIndex - 2][0] =  std::max(_marginsPerSection[marginIndex - 2][0], coilSectionInterface.get_total_margin_tape_distance() / 2);
                _marginsPerSection[marginIndex - 2][1] =  std::max(_marginsPerSection[marginIndex - 2][1], coilSectionInterface.get_total_margin_tape_distance() / 2);
            }
        }
    }
}

void Coil::equalize_margins(std::vector<std::pair<ElectricalType, std::pair<size_t, double>>> orderedSectionsWithInsulation) {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    // Mirror apply_margin_tape's sizing guard: _marginsPerSection is sized
    // lazily by wind_by_*; equalize_margins can be reached on paths where it
    // was never grown to the section count (e.g. PSFB / multi-section bridge
    // topologies). Reading past the end here used to SEGV in CoilAdviser.
    if (_marginsPerSection.size() < orderedSectionsWithInsulation.size()) {
        _marginsPerSection.resize(orderedSectionsWithInsulation.size(), {0, 0});
    }

    for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
        if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
            if (!_coilSectionInterfaces.empty()) {

                size_t indexForMarginLeftSection = sectionIndex;
                size_t indexForMarginRightSection;
                // The "right" section is two ahead (conduction → insulation →
                // conduction). When near the end, wrap to the first section.
                // The original `!= size()-2` test missed the size()-1 case,
                // letting sectionIndex+2 read past the end.
                if (sectionIndex + 2 < orderedSectionsWithInsulation.size()) {
                    indexForMarginRightSection = sectionIndex + 2;
                }
                else {
                    indexForMarginRightSection = 0;
                }

                auto windingIndex = orderedSectionsWithInsulation[indexForMarginLeftSection].second.first;
                auto previousWindingIndex = orderedSectionsWithInsulation[indexForMarginRightSection].second.first;
                auto windingsMapKey = std::pair<size_t, size_t>{previousWindingIndex, windingIndex};
                auto coilSectionInterface = _coilSectionInterfaces[windingsMapKey];
                double totalMargin = _marginsPerSection[indexForMarginLeftSection][1] + _marginsPerSection[indexForMarginRightSection][0];
                double leftAvailableSpace = orderedSectionsWithInsulation[indexForMarginLeftSection].second.second;
                double rightAvailableSpace = orderedSectionsWithInsulation[indexForMarginRightSection].second.second;
                double totalAvailableSpace = leftAvailableSpace + rightAvailableSpace;
                _marginsPerSection[indexForMarginLeftSection][1] = leftAvailableSpace / totalAvailableSpace * totalMargin;
                _marginsPerSection[indexForMarginRightSection][0] = rightAvailableSpace / totalAvailableSpace * totalMargin;
            }
        }
    }
}

bool Coil::wind_by_sections() {
    auto bobbin = resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    std::vector<double> proportionPerWinding;
    auto windingOrientation = get_winding_orientation();
    auto sectionAlignment = get_section_alignment();

    if (bobbinWindingWindowShape == WindingWindowShape::ROUND && windingOrientation == WindingOrientation::CONTIGUOUS && sectionAlignment == CoilAlignment::SPREAD ) {
        proportionPerWinding = make_equal_proportion_per_winding(get_functional_description().size());
    }
    else {
        proportionPerWinding = get_proportion_per_winding_based_on_wires();
    }
    return wind_by_sections(proportionPerWinding);
}

bool Coil::wind_by_sections(size_t repetitions){
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    auto proportionPerWinding = make_equal_proportion_per_winding(get_functional_description().size());
    return wind_by_sections(proportionPerWinding, pattern, repetitions);
}

bool Coil::wind_by_sections(std::vector<size_t> pattern, size_t repetitions) {
    auto proportionPerWinding = make_equal_proportion_per_winding(get_functional_description().size());
    return wind_by_sections(proportionPerWinding, pattern, repetitions);
}

bool Coil::wind_by_sections(std::vector<double> proportionPerWinding) {
    std::vector<size_t> pattern;
    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        pattern.push_back(windingIndex);
    }
    return wind_by_sections(proportionPerWinding, pattern, _interleavingLevel);
}

bool Coil::create_default_group(Bobbin bobbin, WiringTechnology coilType, double coreToLayerDistance) {
    // Single-window helper. Multi-window dispatch lives in
    // create_default_groups() — see the multi-column plan §10.
    Group group;
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    group.set_coordinates({windingWindows[0].get_coordinates().value()[0], windingWindows[0].get_coordinates().value()[1]});
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        group.set_dimensions(std::vector<double>{windingWindows[0].get_width().value() - coreToLayerDistance * 2, windingWindows[0].get_height().value()});
        group.set_coordinate_system(CoordinateSystem::CARTESIAN);
    }
    else {
        group.set_dimensions(std::vector<double>{windingWindows[0].get_radial_height().value() - coreToLayerDistance * 2, windingWindows[0].get_angle().value()});
        group.set_coordinate_system(CoordinateSystem::POLAR);
    }
    group.set_name("Default group");
    std::vector<PartialWinding> partialWindings;

    double numberWindings = get_functional_description().size();
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        PartialWinding partialWinding;
        partialWinding.set_winding(get_name(windingIndex));
        partialWinding.set_parallels_proportion(std::vector<double>(get_number_parallels(windingIndex), 1));
        partialWindings.push_back(partialWinding);
    }
    group.set_partial_windings(partialWindings);
    group.set_sections_orientation(get_winding_orientation());
    group.set_type(coilType);
    set_groups_description(std::vector<Group>{group});

    return true;
}

bool Coil::create_default_groups(Bobbin bobbin, WiringTechnology coilType, double coreToLayerDistance) {
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (windingWindows.size() <= 1) {
        return create_default_group(bobbin, coilType, coreToLayerDistance);
    }

    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    // Distribute the windings by their functional placement: each winding goes to
    // the winding window its windingWindow field names, defaulting to window 0 (the
    // schema-documented default when the field is absent).
    double numberWindings = get_functional_description().size();
    std::vector<std::vector<PartialWinding>> partialWindingsPerWindow(windingWindows.size());
    bool anyExplicitPlacement = false;
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        auto requestedWindow = get_functional_description()[windingIndex].get_winding_window();
        if (requestedWindow) {
            anyExplicitPlacement = true;
            if (requestedWindow.value() < 0 || static_cast<size_t>(requestedWindow.value()) >= windingWindows.size()) {
                throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                    "Winding " + get_name(windingIndex) + " requests winding window " +
                    std::to_string(requestedWindow.value()) + " but the bobbin has " +
                    std::to_string(windingWindows.size()) + " winding windows");
            }
        }
        size_t windowIndex = requestedWindow ? static_cast<size_t>(requestedWindow.value()) : 0;
        PartialWinding partialWinding;
        partialWinding.set_winding(get_name(windingIndex));
        partialWinding.set_parallels_proportion(std::vector<double>(get_number_parallels(windingIndex), 1));
        partialWindingsPerWindow[windowIndex].push_back(partialWinding);
    }

    std::vector<Group> groups;
    for (size_t i = 0; i < windingWindows.size(); ++i) {
        Group g;
        g.set_name("Column " + std::to_string(i));
        g.set_winding_window(static_cast<int64_t>(i));
        g.set_coordinates({windingWindows[i].get_coordinates().value()[0], windingWindows[i].get_coordinates().value()[1]});
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            g.set_dimensions(std::vector<double>{windingWindows[i].get_width().value() - coreToLayerDistance * 2, windingWindows[i].get_height().value()});
            g.set_coordinate_system(CoordinateSystem::CARTESIAN);
        }
        else {
            g.set_dimensions(std::vector<double>{windingWindows[i].get_radial_height().value() - coreToLayerDistance * 2, windingWindows[i].get_angle().value()});
            g.set_coordinate_system(CoordinateSystem::POLAR);
        }
        g.set_sections_orientation(get_winding_orientation());
        g.set_type(coilType);
        g.set_partial_windings(partialWindingsPerWindow[i]);
        groups.push_back(g);
    }
    set_groups_description(groups);

    if (!anyExplicitPlacement) {
        OM_WARNING("Multi-column bobbin detected (" + std::to_string(windingWindows.size()) +
                   " winding windows) and no winding carries a windingWindow placement. All windings "
                   "placed in window 0 by default. Set windingWindow on the windings or call "
                   "assign_windings_to_columns() to distribute.");
    }
    return true;
}

size_t Coil::find_window_index_for_group(const std::string& groupName) const {
    auto groupsOpt = get_groups_description();
    if (!groupsOpt) return 0;
    auto groups = groupsOpt.value();

    Bobbin bobbinResolved = const_cast<Coil*>(this)->resolve_bobbin();
    if (!bobbinResolved.get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed, cannot resolve the winding window of group " + groupName);
    }
    auto windingWindows = bobbinResolved.get_processed_description().value().get_winding_windows();
    if (windingWindows.size() <= 1) return 0;

    for (auto& g : groups) {
        if (g.get_name() == groupName) {
            // Prefer the explicit windingWindow reference (stamped by
            // create_default_groups / assign_windings_to_columns, or supplied
            // in the MAS file).
            if (g.get_winding_window()) {
                auto windowIndex = g.get_winding_window().value();
                if (windowIndex < 0 || static_cast<size_t>(windowIndex) >= windingWindows.size()) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                        "Group " + groupName + " references winding window " + std::to_string(windowIndex) +
                        " but the bobbin has " + std::to_string(windingWindows.size()) + " winding windows");
                }
                return static_cast<size_t>(windowIndex);
            }
            auto gc = g.get_coordinates();
            for (size_t j = 0; j < windingWindows.size(); ++j) {
                auto wwc = windingWindows[j].get_coordinates().value();
                // Match on x and y (group coordinates are 2D, window are 3D).
                if (gc.size() >= 2 && wwc.size() >= 2 &&
                    std::abs(gc[0] - wwc[0]) < 1e-9 &&
                    std::abs(gc[1] - wwc[1]) < 1e-9) {
                    return j;
                }
            }
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                "Group " + groupName + " carries no windingWindow reference and its coordinates match no "
                "winding window of the multi-window bobbin; cannot resolve its placement");
        }
    }
    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
        "Group " + groupName + " does not exist in the groups description");
}

size_t Coil::resolve_section_winding_window_index(const Section& section) const {
    if (section.get_winding_window()) {
        auto windowIndex = section.get_winding_window().value();
        if (windowIndex < 0) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                "Section " + section.get_name() + " references negative winding window " + std::to_string(windowIndex));
        }
        return static_cast<size_t>(windowIndex);
    }
    if (section.get_group()) {
        return find_window_index_for_group(section.get_group().value());
    }
    return 0;
}

WoundColumnFrame Coil::get_wound_column_frame_for_section(const std::string& sectionName) {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();

    auto buildMainColumnFrame = [&]() {
        WoundColumnFrame frame;
        frame.shape = bobbinProcessedDescription.get_column_shape();
        frame.columnDepth = bobbinProcessedDescription.get_column_depth();
        if (bobbinProcessedDescription.get_column_width()) {
            frame.columnWidth = bobbinProcessedDescription.get_column_width().value();
        }
        else {
            auto bobbinWindingWindow = bobbinProcessedDescription.get_winding_windows()[0];
            frame.columnWidth = bobbinWindingWindow.get_coordinates().value()[0] - bobbinWindingWindow.get_width().value() / 2;
        }
        frame.axisX = 0;
        return frame;
    };

    if (!get_sections_description()) {
        throw CoilNotProcessedException("Sections description is missing, cannot resolve the wound column of section " + sectionName);
    }
    auto sections = get_sections_description().value();
    std::optional<Section> section;
    for (auto& candidate : sections) {
        if (candidate.get_name() == sectionName) {
            section = candidate;
            break;
        }
    }
    if (!section) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
            "Section " + sectionName + " does not exist in the sections description");
    }

    size_t windowIndex = resolve_section_winding_window_index(section.value());
    if (windowIndex == 0) {
        return buildMainColumnFrame();
    }

    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (windowIndex >= windingWindows.size()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
            "Section " + sectionName + " references winding window " + std::to_string(windowIndex) +
            " but the bobbin has " + std::to_string(windingWindows.size()) + " winding windows");
    }
    auto columnEdge = windingWindows[windowIndex].get_column();
    if (!columnEdge) {
        // Schema default: a window without a column edge wraps the main column
        // (e.g. the stacked chambers of a split bobbin).
        return buildMainColumnFrame();
    }
    if (!_coreColumns) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
            "Section " + sectionName + " is placed in winding window " + std::to_string(windowIndex) +
            " wound around core column " + std::to_string(columnEdge.value()) +
            ", but the core columns were not provided; call set_core_columns before winding");
    }
    auto columns = _coreColumns.value();
    if (columnEdge.value() < 0 || static_cast<size_t>(columnEdge.value()) >= columns.size()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
            "Winding window " + std::to_string(windowIndex) + " references core column " +
            std::to_string(columnEdge.value()) + " but the core has " + std::to_string(columns.size()) + " columns");
    }
    auto column = columns[static_cast<size_t>(columnEdge.value())];
    if (column.get_coordinates()[0] == 0) {
        // The window's column IS the main column (shared-region entries pointing back).
        return buildMainColumnFrame();
    }
    WoundColumnFrame frame;
    frame.shape = column.get_shape();
    // Non-main columns carry no bobbin: wrap the bare column.
    frame.columnWidth = column.get_width() / 2;
    frame.columnDepth = column.get_depth() / 2;
    // The winding frame is the +x side; mirrored (negative-x) windows are wound
    // against the mirrored column and flipped into place afterwards.
    frame.axisX = std::abs(column.get_coordinates()[0]);
    return frame;
}

std::optional<double> Coil::get_turn_length_in_frame(const WoundColumnFrame& frame, double turnX) {
    double radius = frame.axisX == 0 ? turnX : std::abs(turnX - frame.axisX);
    double length;
    if (frame.shape == ColumnShape::ROUND) {
        length = 2 * std::numbers::pi * radius;
    }
    else if (frame.shape == ColumnShape::OBLONG) {
        length = 2 * std::numbers::pi * radius + 4 * (frame.columnDepth - frame.columnWidth);
    }
    else if (frame.shape == ColumnShape::RECTANGULAR || frame.shape == ColumnShape::IRREGULAR) {
        length = 4 * frame.columnDepth + 4 * frame.columnWidth + 2 * std::numbers::pi * (radius - frame.columnWidth);
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "only round or rectangular columns supported for bobbins");
    }
    if (length < 0) {
        return std::nullopt;
    }
    return length;
}

void Coil::apply_group_window_sides() {
    if (!get_groups_description() || !get_sections_description()) {
        return;
    }
    Bobbin bobbinResolved = resolve_bobbin();
    if (!bobbinResolved.get_processed_description()) {
        return;
    }
    auto windingWindows = bobbinResolved.get_processed_description().value().get_winding_windows();
    if (windingWindows.size() <= 1) {
        return;
    }

    // Collect the sections whose winding window sits on the negative-x side.
    std::set<std::string> mirroredSections;
    auto sections = get_sections_description().value();
    for (auto& section : sections) {
        size_t windowIndex = resolve_section_winding_window_index(section);
        if (windowIndex >= windingWindows.size()) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                "Section " + section.get_name() + " references winding window " + std::to_string(windowIndex) +
                " but the bobbin has " + std::to_string(windingWindows.size()) + " winding windows");
        }
        if (windingWindows[windowIndex].get_coordinates().value()[0] < 0) {
            mirroredSections.insert(section.get_name());
        }
    }
    if (mirroredSections.empty()) {
        return;
    }

    auto mirrorX = [](std::vector<double> coordinates) {
        coordinates[0] = -coordinates[0];
        return coordinates;
    };

    for (auto& section : sections) {
        if (mirroredSections.contains(section.get_name())) {
            section.set_coordinates(mirrorX(section.get_coordinates()));
        }
    }
    set_sections_description(sections);

    if (get_layers_description()) {
        auto layers = get_layers_description().value();
        for (auto& layer : layers) {
            if (layer.get_section() && mirroredSections.contains(layer.get_section().value())) {
                layer.set_coordinates(mirrorX(layer.get_coordinates()));
                if (layer.get_additional_coordinates()) {
                    auto additionalCoordinates = layer.get_additional_coordinates().value();
                    for (auto& coordinates : additionalCoordinates) {
                        coordinates = mirrorX(coordinates);
                    }
                    layer.set_additional_coordinates(additionalCoordinates);
                }
            }
        }
        set_layers_description(layers);
    }

    if (get_turns_description()) {
        auto turns = get_turns_description().value();
        for (auto& turn : turns) {
            if (turn.get_section() && mirroredSections.contains(turn.get_section().value())) {
                turn.set_coordinates(mirrorX(turn.get_coordinates()));
                if (turn.get_additional_coordinates()) {
                    auto additionalCoordinates = turn.get_additional_coordinates().value();
                    for (auto& coordinates : additionalCoordinates) {
                        coordinates = mirrorX(coordinates);
                    }
                    turn.set_additional_coordinates(additionalCoordinates);
                }
                // A mirrored turn is wound the opposite way around its column.
                if (turn.get_orientation()) {
                    turn.set_orientation(turn.get_orientation().value() == TurnOrientation::CLOCKWISE
                                             ? TurnOrientation::COUNTER_CLOCKWISE
                                             : TurnOrientation::CLOCKWISE);
                }
            }
        }
        set_turns_description(turns);
    }
}

void Coil::assign_windings_to_columns(const std::vector<std::vector<size_t>>& windingIndicesPerColumn) {
    auto bobbin = resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    if (windingIndicesPerColumn.size() != windingWindows.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "windingIndicesPerColumn size (" + std::to_string(windingIndicesPerColumn.size()) +
            ") must match number of winding windows (" + std::to_string(windingWindows.size()) + ")");
    }

    auto& functionalDesc = get_functional_description();
    std::vector<Group> groups;

    for (size_t col = 0; col < windingWindows.size(); ++col) {
        Group g;
        g.set_name("Column " + std::to_string(col));
        g.set_winding_window(static_cast<int64_t>(col));
        g.set_coordinates({windingWindows[col].get_coordinates().value()[0], windingWindows[col].get_coordinates().value()[1]});
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            g.set_dimensions(std::vector<double>{windingWindows[col].get_width().value(), windingWindows[col].get_height().value()});
            g.set_coordinate_system(CoordinateSystem::CARTESIAN);
        }
        else {
            g.set_dimensions(std::vector<double>{windingWindows[col].get_radial_height().value(), windingWindows[col].get_angle().value()});
            g.set_coordinate_system(CoordinateSystem::POLAR);
        }
        g.set_sections_orientation(get_winding_orientation());
        g.set_type(WiringTechnology::WOUND);

        std::vector<PartialWinding> partialWindings;
        for (auto windingIndex : windingIndicesPerColumn[col]) {
            if (windingIndex >= functionalDesc.size()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "winding index " + std::to_string(windingIndex) + " out of range");
            }
            PartialWinding pw;
            pw.set_winding(get_name(windingIndex));
            pw.set_parallels_proportion(std::vector<double>(get_number_parallels(windingIndex), 1));
            partialWindings.push_back(pw);
        }
        g.set_partial_windings(partialWindings);
        groups.push_back(g);
    }
    set_groups_description(groups);
}

bool Coil::wind_by_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    _currentProportionPerWinding = proportionPerWinding;
    _currentPattern = pattern;
    _currentRepetitions = repetitions;

    if (repetitions <= 0) {
        throw InvalidInputException("Interleaving levels must be greater than 0");
    }

    auto bobbin = resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (!bobbin.get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    // Multi-window supported: initialise sections_orientation and
    // sections_alignment for ALL winding windows.
    for (size_t i = 0; i < windingWindows.size(); ++i) {
        if (!windingWindows[i].get_sections_orientation()) {
            windingWindows[i].set_sections_orientation(_windingOrientation);
        }
        if (!windingWindows[i].get_sections_alignment()) {
            windingWindows[i].set_sections_alignment(_sectionAlignmentExplicit ? _sectionAlignment : get_section_alignment());
        }
    }
    bobbinProcessedDescription.set_winding_windows(windingWindows);
    bobbin.set_processed_description(bobbinProcessedDescription);
    set_bobbin(bobbin);

    if (!get_groups_description()) {
        create_default_groups(bobbin);
    }

    set_sections_description(std::nullopt);
    set_layers_description(std::nullopt);
    set_turns_description(std::nullopt);

    std::vector<size_t> maybeVirtualizedPattern = pattern;
    std::vector<double> maybeVirtualizedProportionPerWinding = proportionPerWinding;
    auto functionalDescription = get_functional_description();
    auto needsVirtualization = needs_virtualization();

    if (needsVirtualization) {
        create_virtualization_map();
        auto virtualFunctionalDescription = virtualize_functional_description();
        maybeVirtualizedPattern = virtualize_pattern(pattern);
        maybeVirtualizedProportionPerWinding = virtualize_proportion_per_winding(proportionPerWinding);
        set_functional_description(virtualFunctionalDescription);
        _windingIndexByName.clear();
        _turnIndexByName.clear();
    }

    bool result;

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        result = wind_by_rectangular_sections(maybeVirtualizedProportionPerWinding, maybeVirtualizedPattern, repetitions);
    }
    else {
        result = wind_by_round_sections(maybeVirtualizedProportionPerWinding, maybeVirtualizedPattern, repetitions);
    }

    if (needsVirtualization) {
        _windingIndexByName.clear();
        _turnIndexByName.clear();
        set_functional_description(functionalDescription);
        // wind_by_(rectangular|round)_sections returns false when wires don't
        // fit the available section space, leaving sections_description in
        // its initial nullopt state. Skip devirtualize in that case — there
        // is nothing to devirtualize, and our caller (fast_wind) treats the
        // missing sections_description as the "couldn't wind" signal.
        if (result && get_sections_description()) {
            devirtualize_sections_description();
        }
    }


    return result;
}

bool Coil::needs_virtualization() {
    for (auto winding : get_functional_description()) {
        if (winding.get_wound_with()) {
            if (winding.get_wound_with()->size() > 0) {
                return true;
            }
        }
    }
    return false;
}

std::vector<Winding> Coil::virtualize_functional_description() {
    std::vector<Winding> newFunctionalDescription;
    for (auto [virtualWindingIndex, windingIndexes] : _virtualizationMap) {
        std::string name = "";
        int64_t numberTurns = 0;
        int64_t numberParallels = 0;
        std::optional<IsolationSide> isolationSide = std::nullopt;
        std::optional<Wire> wire = std::nullopt;
        std::vector<ConnectionElement> connections;
        for (auto windingIndex : windingIndexes) {
            auto winding = get_functional_description()[windingIndex];
            numberTurns += winding.get_number_turns();

            if (numberParallels == 0) {
                numberParallels = winding.get_number_parallels();
            }
            else {
                if (numberParallels != winding.get_number_parallels()) {
                    throw InvalidInputException("Windings wound together must have the same number of parallels");
                }
            }

            if (!isolationSide) {
                isolationSide = winding.get_isolation_side();
            }
            else {
                if (isolationSide.value() != winding.get_isolation_side()) {
                    throw InvalidInputException("Windings wound together must have the same isolation side");
                }
            }

            if (!wire) {
                wire = winding.resolve_wire();
            }
            else {
                if (wire.value() != winding.resolve_wire()) {
                    throw InvalidInputException("Windings wound together must have the same wire");
                }
            }

            if (winding.get_connections()) {
                auto windingConnections = winding.get_connections().value();
                for (auto connection : windingConnections) {
                    connections.push_back(connection);
                }
            }

            if (name == "") {
                size_t minimumWindingGroupIndex = get_winding_group_minimum_index(windingIndex);
                name = get_functional_description()[minimumWindingGroupIndex].get_name();
            }
        }

        Winding newWinding;
        newWinding.set_connections(connections);
        newWinding.set_isolation_side(isolationSide.value());
        newWinding.set_name(name);
        newWinding.set_number_parallels(numberParallels);
        newWinding.set_number_turns(numberTurns);
        newWinding.set_wire(wire.value());
        newFunctionalDescription.push_back(newWinding);
    }
    return newFunctionalDescription;
}

Section Coil::devirtualize_section(Section section) {
    if (section.get_type() == ElectricalType::INSULATION) {
        return section;
    }
    std::vector<PartialWinding> newPartialWindings;
    for (auto partialWinding : section.get_partial_windings()) {
        auto virtualWindingIndex = SIZE_MAX;
        for (auto [auxVirtualWindingIndex, virtualWindingName]  : _virtualWindingNames) {
            if (partialWinding.get_winding() == virtualWindingName) {
                virtualWindingIndex = auxVirtualWindingIndex;
            }
        }

        if (virtualWindingIndex == SIZE_MAX) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened looking for virtual indexes");
        }

        for (auto windingIndex : _virtualizationMap[virtualWindingIndex]) {
            auto newPartialWinding = partialWinding;
            auto name = _windingNames[windingIndex];
            newPartialWinding.set_winding(name);
            newPartialWindings.push_back(newPartialWinding);
        }
    }
    section.set_partial_windings(newPartialWindings);
    return section;
}

Layer Coil::devirtualize_layer(Layer layer) {
    if (layer.get_type() == ElectricalType::INSULATION) {
        return layer;
    }
    std::vector<PartialWinding> newPartialWindings;
    for (auto partialWinding : layer.get_partial_windings()) {
        auto virtualWindingIndex = SIZE_MAX;
        for (auto [auxVirtualWindingIndex, virtualWindingName]  : _virtualWindingNames) {
            if (partialWinding.get_winding() == virtualWindingName) {
                virtualWindingIndex = auxVirtualWindingIndex;
            }
        }

        if (virtualWindingIndex == SIZE_MAX) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened looking for virtual indexes");
        }

        for (auto windingIndex : _virtualizationMap[virtualWindingIndex]) {
            auto newPartialWinding = partialWinding;
            auto name = _windingNames[windingIndex];
            newPartialWinding.set_winding(name);
            newPartialWindings.push_back(newPartialWinding);
        }
    }
    layer.set_partial_windings(newPartialWindings);
    return layer;
}

Turn Coil::devirtualize_turn(Turn turn, std::string virtualWindingName, std::string windingName, size_t newParallelIndex) {
    auto name = turn.get_name();
    name = std::regex_replace(name, std::regex(virtualWindingName), windingName);
    turn.set_name(name);
    turn.set_winding(windingName);
    turn.set_parallel(newParallelIndex);
    return turn;
}

Section Coil::virtualize_section(Section section) {
    if (section.get_type() == ElectricalType::INSULATION) {
        return section;
    }
    auto partialWindings = section.get_partial_windings();
    auto firstWindingName = partialWindings[0].get_winding();
    auto firstWindingIndex = get_winding_index_by_name(firstWindingName);
    auto virtualWindingIndex = 0;
    bool found = false;
    for (auto [virtualIndex, indexes] : _virtualizationMap) {
        for (auto index : indexes) {
            if (firstWindingIndex == index) {
                found = true;
                virtualWindingIndex = virtualIndex;
                break;
            }
        }
    }
    if (!found) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened virtualizing section");
    }
    found = false;
    for (auto partialWinding : partialWindings) {
        if (partialWinding.get_winding() == _virtualWindingNames[virtualWindingIndex]) {
            found = true;
            section.set_partial_windings({partialWinding});
            break;
        }
    }
    if (!found) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened virtualizing section 2");
    }

    return section;
}

Layer Coil::virtualize_layer(Layer layer) {
    if (layer.get_type() == ElectricalType::INSULATION) {
        return layer;
    }
    auto partialWindings = layer.get_partial_windings();
    auto firstWindingName = partialWindings[0].get_winding();
    auto firstWindingIndex = get_winding_index_by_name(firstWindingName);
    auto virtualWindingIndex = 0;
    bool found = false;
    for (auto [virtualIndex, indexes] : _virtualizationMap) {
        for (auto index : indexes) {
            if (firstWindingIndex == index) {
                found = true;
                virtualWindingIndex = virtualIndex;
                break;
            }
        }
    }
    if (!found) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened virtualizing layer");
    }
    found = false;
    for (auto partialWinding : partialWindings) {
        if (partialWinding.get_winding() == _virtualWindingNames[virtualWindingIndex]) {
            found = true;
            layer.set_partial_windings({partialWinding});
            break;
        }
    }
    if (!found) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened virtualizing layer 2");
    }
    return layer;
}

Turn Coil::virtualize_turn(Turn turn, std::string virtualWindingName, std::string windingName) {
    auto name = turn.get_name();
    name = std::regex_replace(name, std::regex(windingName), virtualWindingName);
    turn.set_name(name);
    turn.set_winding(virtualWindingName);
    return turn;
}

void Coil::devirtualize_sections_description() {
    if (!get_sections_description()) {
        // Caller is responsible for not invoking devirtualize when wind_by_*
        // returned false. Be defensive: keep nullopt and return.
        return;
    }
    std::vector<Section> newSectionsDescription;
    auto sections = get_sections_description().value();
    for (auto section : sections) {
        auto newSection = devirtualize_section(section);
        newSectionsDescription.push_back(newSection);
    }
    set_sections_description(newSectionsDescription);
}

void Coil::devirtualize_layers_description() {
    if (!get_layers_description()) {
        // wind_by_round_layers()/wind_by_rectangular_layers() return false
        // (leaving layers_description unset) when the winding does not fit
        // the window — but wind_by_layers() still runs this devirtualization
        // cleanup unconditionally for virtualized (wound_with) coils. Mirror
        // the guard in devirtualize_sections_description(): a failed wind has
        // nothing to devirtualize, so keep nullopt and let the caller treat
        // the candidate as not-wound instead of throwing bad_optional_access.
        return;
    }
    std::vector<Layer> newLayersDescription;
    auto layers = get_layers_description().value();
    for (auto layer : layers) {
        auto newLayer = devirtualize_layer(layer);
        newLayersDescription.push_back(newLayer);
    }
    set_layers_description(newLayersDescription);
}

void Coil::devirtualize_turns_description() {
    std::vector<Turn> newTurnsDescription;
    auto turns = get_turns_description().value();

    for (auto [virtualWindingIndex, windingIndexes] : _virtualizationMap) {
        auto turnsInVirtualWinding = get_turns_by_winding(_virtualWindingNames[virtualWindingIndex]);
        
        std::map<size_t, int64_t> remainingNumberTurnsPerWoundTogetherWinding;
        std::map<size_t, int64_t> assignedTurnsPerWinding;  // Track how many turns assigned to each winding
        int64_t minimumNumberTurns = std::numeric_limits<int64_t>::max();
        int64_t totalNumberTurns = 0;
        for (auto windingIndex : windingIndexes) {
            auto numberTurns = get_functional_description()[windingIndex].get_number_turns() * get_functional_description()[windingIndex].get_number_parallels();
            minimumNumberTurns = std::min(minimumNumberTurns, numberTurns);
            totalNumberTurns += numberTurns;
            remainingNumberTurnsPerWoundTogetherWinding[windingIndex] = numberTurns;
            assignedTurnsPerWinding[windingIndex] = 0;
        }

        if (size_t(totalNumberTurns) != turnsInVirtualWinding.size()) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened devirtualizing turns 1");
        }

        std::vector<size_t> devirtualizingPattern = {};
        for (auto windingIndex : windingIndexes) {
            double numberTurns = get_functional_description()[windingIndex].get_number_turns() * get_functional_description()[windingIndex].get_number_parallels();
            size_t numberTurnsPerPattern = round(numberTurns / minimumNumberTurns);
            if (numberTurnsPerPattern == 0) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened devirtualizing turns 2");
            }
            for (size_t index = 0; index < numberTurnsPerPattern; ++index) {
                devirtualizingPattern.push_back(windingIndex);
            }
        }

        size_t devirtualizingPatternIndex = 0;
        for (auto turn : turnsInVirtualWinding) {
            size_t timeout = devirtualizingPattern.size() + 1;
            while (timeout > 0) {
                if (remainingNumberTurnsPerWoundTogetherWinding[devirtualizingPattern[devirtualizingPatternIndex]] > 0) {
                    auto targetWindingIndex = devirtualizingPattern[devirtualizingPatternIndex];
                    auto newWindingName = get_functional_description()[targetWindingIndex].get_name();
                    auto oldWindingName = _virtualWindingNames[virtualWindingIndex];
                    
                    // Calculate the correct parallel index for this winding
                    // Each winding has numberTurns turns per parallel
                    auto numberTurnsPerWinding = get_functional_description()[targetWindingIndex].get_number_turns();
                    auto numberParallels = get_functional_description()[targetWindingIndex].get_number_parallels();
                    size_t newParallelIndex = assignedTurnsPerWinding[targetWindingIndex] / numberTurnsPerWinding;
                    // Clamp to valid range in case of rounding issues
                    if (newParallelIndex >= size_t(numberParallels)) {
                        newParallelIndex = numberParallels - 1;
                    }
                    
                    remainingNumberTurnsPerWoundTogetherWinding[targetWindingIndex]--;
                    assignedTurnsPerWinding[targetWindingIndex]++;
                    devirtualizingPatternIndex = (devirtualizingPatternIndex + 1) % devirtualizingPattern.size();
                    auto newTurn = devirtualize_turn(turn, oldWindingName, newWindingName, newParallelIndex);
                    newTurnsDescription.push_back(newTurn);
                    break;
                }
                else {
                    devirtualizingPatternIndex = (devirtualizingPatternIndex + 1) % devirtualizingPattern.size();
                }
                timeout--;
            }
            if (timeout == 0) {
                throw CalculationException(ErrorCode::CALCULATION_DIVERGED, "Something wrong happened devirtualizing turns 3");
            }
        }
    }

    set_turns_description(newTurnsDescription);
}

std::vector<Section> Coil::virtualize_sections_description() {
    std::vector<Section> newSectionsDescription;
    auto sections = get_sections_description().value();
    for (auto section : sections) {
        auto newSection = virtualize_section(section);
        newSectionsDescription.push_back(newSection);
    }
    return newSectionsDescription;
}

std::vector<Layer> Coil::virtualize_layers_description() {
    std::vector<Layer> newLayersDescription;
    auto layers = get_layers_description().value();
    for (auto layer : layers) {
        auto newLayer = virtualize_layer(layer);
        newLayersDescription.push_back(newLayer);
    }
    return newLayersDescription;
}

std::map<std::pair<size_t, size_t>, std::vector<Layer>> Coil::virtualize_insulation_intersections_layers() {
    std::map<std::pair<size_t, size_t>, std::vector<Layer>> newInsulationInterSectionsLayers;
    for (auto [key, layers] : _insulationInterSectionsLayers) {
        size_t newFirstIndex = 0;
        size_t newSecondIndex = 0;
        bool newFirstIndexFound = false;
        bool newSecondIndexFound = false;
        for (auto [virtualIndex , indexes] : _virtualizationMap) {
            for (auto index : indexes) {
                if (key.first == index) {
                    newFirstIndex = virtualIndex;
                    newFirstIndexFound = true;
                }
                if (key.second == index) {
                    newSecondIndex = virtualIndex;
                    newSecondIndexFound = true;
                }
            }
        }

        if (!newFirstIndexFound || !newSecondIndexFound) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened looking for virtual indexes for insulation intersections layers");
        }
        std::pair<size_t, size_t> virtualKey = {newFirstIndex, newSecondIndex};
        newInsulationInterSectionsLayers[virtualKey] = layers;
    }
    return newInsulationInterSectionsLayers;
}

std::map<size_t, std::vector<size_t>> Coil::create_virtualization_map() {
    std::map<size_t, size_t> inversedVirtualizationMap;

    _windingIndexByName.clear();
    _turnIndexByName.clear();
    _virtualizationMap.clear();
    size_t currentVirtualIndex = 0;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        _windingNames[windingIndex] = get_functional_description()[windingIndex].get_name();
        size_t minimumWindingGroupIndex = get_winding_group_minimum_index(windingIndex);
        if (inversedVirtualizationMap.count(minimumWindingGroupIndex)){
            inversedVirtualizationMap[windingIndex] = inversedVirtualizationMap[minimumWindingGroupIndex];
        }
        else {
            inversedVirtualizationMap[windingIndex] = currentVirtualIndex;
            currentVirtualIndex++;
        }
    }

    for (auto [windingIndex, virtualWindingIndex] : inversedVirtualizationMap) {
        _virtualWindingNames[virtualWindingIndex] = get_functional_description()[get_winding_group_minimum_index(windingIndex)].get_name();
        _virtualizationMap[virtualWindingIndex].push_back(windingIndex);

    }
    return _virtualizationMap;
}

size_t Coil::get_winding_group_minimum_index(size_t currentWindingIndex) {
    size_t minimumWindingIndex = currentWindingIndex;
    if (get_functional_description()[currentWindingIndex].get_wound_with()) {
        auto windingsWoundWith = get_functional_description()[currentWindingIndex].get_wound_with().value();
        for (auto windingWoundWith : windingsWoundWith) {
            minimumWindingIndex = std::min(get_winding_index_by_name(windingWoundWith), minimumWindingIndex);
        }
    }
    return minimumWindingIndex;
}

std::vector<size_t> Coil::virtualize_pattern(std::vector<size_t> pattern) {
    std::vector<size_t> newPattern;
    if (_virtualizationMap.size() == 0) {
        throw CoilNotProcessedException("No virtualization loaded. Did you forget to call create_virtualization_map()?");
    }
    for (auto windingIndex : pattern) {
        auto winding = get_functional_description()[windingIndex];

        size_t windingIndexForSearch;
        if (!winding.get_wound_with()) {
            windingIndexForSearch = windingIndex;
        }
        else {
            windingIndexForSearch = get_winding_group_minimum_index(windingIndex);
        }

        size_t virtualWindingIndex = SIZE_MAX;
        for (auto [auxVirtualWindingIndex, windingIndexes] : _virtualizationMap) {
            if(std::find(windingIndexes.begin(), windingIndexes.end(), windingIndexForSearch) != windingIndexes.end()) {
                virtualWindingIndex = auxVirtualWindingIndex;
            }
        }

        if (virtualWindingIndex == SIZE_MAX) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happened looking for virtual indexes");
        }

        newPattern.push_back(virtualWindingIndex);
    }
    newPattern = compress_pattern(newPattern);
    return newPattern;
}

std::vector<size_t> Coil::compress_pattern(std::vector<size_t> pattern) {
    std::vector<size_t> newPattern;
    for (auto windingIndex : pattern) {
        if (newPattern.size() == 0) {
            newPattern.push_back(windingIndex);
        }
        else if (newPattern.back() != windingIndex) {
            newPattern.push_back(windingIndex);
        }
    }
    return newPattern;
}

std::vector<double> Coil::virtualize_proportion_per_winding(std::vector<double> proportionPerWinding) {
    std::vector<double> newProportionPerWinding;
    if (_virtualizationMap.size() == 0) {
        throw CoilNotProcessedException("No virtualization loaded. Did you forget to call create_virtualization_map()?");
    }

    for (auto [virtualWindingIndex, windingIndexes] : _virtualizationMap) {
        double newProportion = 0;
        for (auto windingIndex : windingIndexes) {
            newProportion += proportionPerWinding[windingIndex];
        }
        newProportionPerWinding.push_back(newProportion);
    }
    return newProportionPerWinding;
}

bool Coil::wind_by_rectangular_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    set_sections_description(std::nullopt);
    std::vector<Section> sectionsDescription;

    if (!get_groups_description()) {
        throw CoilNotProcessedException("At least default group must be defined at this point.");
    }

    auto groups = get_groups_description().value();
    std::vector<std::vector<double>> remainingParallelsProportion;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        remainingParallelsProportion.push_back(std::vector<double>(get_number_parallels(windingIndex), 1));
    }
    // Multi-column winding: each group is wound independently with the subset of the
    // pattern that belongs to its windings, in a window-local frame on the +x side of
    // the main column (apply_group_window_sides mirrors negative-x windows at the end
    // of wind()). Single-group coils keep the exact historical inputs.
    bool multiGroup = groups.size() > 1;
    size_t marginSectionOffset = 0;

    for (auto group : groups) {
        std::vector<size_t> groupPattern = pattern;
        std::vector<double> groupProportionPerWinding = proportionPerWinding;
        std::set<size_t> groupWindingIndexes;
        std::optional<size_t> groupWindowIndex;
        if (multiGroup) {
            for (auto& groupPartialWinding : group.get_partial_windings()) {
                groupWindingIndexes.insert(get_winding_index_by_name(groupPartialWinding.get_winding()));
            }
            if (groupWindingIndexes.empty()) {
                continue;
            }
            groupPattern.clear();
            for (auto windingIndex : pattern) {
                if (groupWindingIndexes.contains(windingIndex)) {
                    groupPattern.push_back(windingIndex);
                }
            }
            if (groupPattern.empty()) {
                throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                    "Group " + group.get_name() + " has windings but none of them appear in the winding pattern");
            }
            // The group's windings share this group's full window space.
            double groupProportionSum = 0;
            for (auto windingIndex : groupWindingIndexes) {
                groupProportionSum += proportionPerWinding[windingIndex];
            }
            for (auto windingIndex : groupWindingIndexes) {
                groupProportionPerWinding[windingIndex] = proportionPerWinding[windingIndex] / groupProportionSum;
            }
            groupWindowIndex = find_window_index_for_group(group.get_name());
            // Wind in the +x window-local frame; mirrored back afterwards.
            if (group.get_coordinates()[0] < 0) {
                auto groupCoordinates = group.get_coordinates();
                groupCoordinates[0] = -groupCoordinates[0];
                group.set_coordinates(groupCoordinates);
            }
        }

        double availableWidth = group.get_dimensions()[0];
        double availableHeight = group.get_dimensions()[1];
        double spaceForSections = 0;
        auto windingOrientation = group.get_sections_orientation();

        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            spaceForSections = availableWidth;
        }
        else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            spaceForSections = availableHeight;
        }

        auto orderedSections = get_ordered_sections(spaceForSections, groupProportionPerWinding, groupPattern, repetitions);

        auto orderedSectionsWithInsulation = add_insulation_to_sections(orderedSections);

        double numberWindings = get_functional_description().size();
        auto numberSectionsPerWinding = std::vector<size_t>(numberWindings, 0);
        auto currentSectionPerWinding = std::vector<size_t>(numberWindings, 0);
        for (auto orderedSection : orderedSectionsWithInsulation) {
            if (orderedSection.first == ElectricalType::CONDUCTION) {
                auto windingIndex = orderedSection.second.first;
                numberSectionsPerWinding[windingIndex]++;
            }
        }

        // wound_with grouping (e.g. AHB / Push-Pull / forward-derived center-tap
        // halves "Sec a" + "Sec b"): the partner winding shares the
        // representative's section, so the pattern legitimately omits its
        // index — leaving its slot count at zero, which then trips
        // wind_by_consecutive_turns's "Number of slots cannot be less than 1"
        // guard. Inherit the representative's slot count for grouped windings
        // so the wire-layout planner sees the same section count its sibling
        // has. Without this, every Path-B-Load → re-wind cycle on a
        // center-tapped magnetic throws an exception.
        for (size_t wIdx = 0; wIdx < numberWindings; ++wIdx) {
            if (numberSectionsPerWinding[wIdx] != 0) continue;
            if (multiGroup && !groupWindingIndexes.contains(wIdx)) continue;
            const auto& wwOpt = get_functional_description()[wIdx].get_wound_with();
            if (!wwOpt || wwOpt->empty()) continue;
            for (const auto& partnerName : wwOpt.value()) {
                bool found = false;
                for (size_t pIdx = 0; pIdx < numberWindings; ++pIdx) {
                    if (get_functional_description()[pIdx].get_name() == partnerName &&
                        numberSectionsPerWinding[pIdx] > 0) {
                        numberSectionsPerWinding[wIdx] = numberSectionsPerWinding[pIdx];
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
        }

        // Fallback inheritance via _virtualizationMap. When a magnetic comes
        // in via the MagneticAdviser auto-wind path (AHB/Push-Pull CT),
        // `wound_with` is sometimes empty even though the coil's virtual map
        // already groups the center-tap halves under the same virtual index.
        // In that case, inherit the slot count from any sibling within the
        // same virtual group that has nonzero sections.
        if (!_virtualizationMap.empty()) {
            for (size_t wIdx = 0; wIdx < numberWindings; ++wIdx) {
                if (numberSectionsPerWinding[wIdx] != 0) continue;
                if (multiGroup && !groupWindingIndexes.contains(wIdx)) continue;
                for (const auto& [virtualIdx, members] : _virtualizationMap) {
                    if (std::find(members.begin(), members.end(), wIdx) == members.end()) continue;
                    for (auto pIdx : members) {
                        if (pIdx < numberWindings && numberSectionsPerWinding[pIdx] > 0) {
                            numberSectionsPerWinding[wIdx] = numberSectionsPerWinding[pIdx];
                            break;
                        }
                    }
                    if (numberSectionsPerWinding[wIdx] != 0) break;
                }
            }
        }

        if (multiGroup) {
            // Windings belonging to other groups have zero sections here; give them a
            // placeholder slot count so the style chooser's zero-slot guard doesn't
            // fire. Their style entry is never read: only this group's windings appear
            // in this group's ordered sections.
            for (size_t wIdx = 0; wIdx < numberWindings; ++wIdx) {
                if (!groupWindingIndexes.contains(wIdx) && numberSectionsPerWinding[wIdx] == 0) {
                    numberSectionsPerWinding[wIdx] = 1;
                }
            }
        }
        auto windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(), get_number_parallels(), numberSectionsPerWinding);

        auto wirePerWinding = get_wires();
        double currentSectionCenterWidth = DBL_MAX;
        double currentSectionCenterHeight = DBL_MAX;

        apply_margin_tape(orderedSectionsWithInsulation, marginSectionOffset);

        for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
            // Margins are indexed flat across ALL groups in winding order.
            size_t marginIndex = marginSectionOffset + sectionIndex;
            if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
                auto sectionInfo = orderedSectionsWithInsulation[sectionIndex].second;
                auto windingIndex = sectionInfo.first;
                auto spaceForSection = sectionInfo.second;

                double currentSectionHeight = 0;
                double currentSectionWidth = 0;
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionWidth = spaceForSection;
                    // currentSectionHeight = availableHeight - _marginsPerSection[sectionIndex][0] - _marginsPerSection[sectionIndex][1];
                    currentSectionHeight = availableHeight;

                    if (currentSectionCenterWidth == DBL_MAX) {
                        currentSectionCenterWidth = group.get_coordinates()[0] - availableWidth / 2;
                    }
                    if (currentSectionCenterHeight == DBL_MAX) {
                        currentSectionCenterHeight = group.get_coordinates()[1];
                    }
                }
                else {
                    currentSectionWidth = availableWidth;
                    currentSectionHeight = spaceForSection;
                    if (currentSectionCenterWidth == DBL_MAX) {
                        currentSectionCenterWidth = group.get_coordinates()[0];
                    }
                    if (currentSectionCenterHeight == DBL_MAX) {
                        currentSectionCenterHeight = group.get_coordinates()[1] + availableHeight / 2;
                    }
                }

                PartialWinding partialWinding;
                Section section;
                partialWinding.set_winding(get_name(windingIndex));

                auto parallelsProportions = get_parallels_proportions(currentSectionPerWinding[windingIndex],
                                                                       numberSectionsPerWinding[windingIndex],
                                                                       get_number_turns(windingIndex),
                                                                       get_number_parallels(windingIndex),
                                                                       remainingParallelsProportion[windingIndex],
                                                                       windByConsecutiveTurns[windingIndex],
                                                                       std::vector<double>(get_number_parallels(windingIndex), 1));

                std::vector<double> sectionParallelsProportion = parallelsProportions.second;
                uint64_t physicalTurnsThisSection = parallelsProportions.first;

                partialWinding.set_parallels_proportion(sectionParallelsProportion);
                if (get_functional_description()[windingIndex].get_wound_with()) {
                    if (get_functional_description()[windingIndex].get_wound_with()->size() > 0) {

                    }
                }

                section.set_name(get_name(windingIndex) +  " section " + std::to_string(currentSectionPerWinding[windingIndex]));
                section.set_partial_windings(std::vector<PartialWinding>{partialWinding});  // TODO: support more than one winding per section?
                section.set_group(group.get_name());
                if (groupWindowIndex) {
                    section.set_winding_window(static_cast<int64_t>(groupWindowIndex.value()));
                }
                section.set_type(ElectricalType::CONDUCTION);
                section.set_margin(_marginsPerSection[marginIndex]);
                section.set_layers_orientation(_layersOrientation);
                section.set_coordinate_system(CoordinateSystem::CARTESIAN);
                
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{currentSectionWidth, currentSectionHeight - _marginsPerSection[marginIndex][0] - _marginsPerSection[marginIndex][1]});
                }
                else {
                    section.set_dimensions(std::vector<double>{currentSectionWidth - _marginsPerSection[marginIndex][0] - _marginsPerSection[marginIndex][1], currentSectionHeight});
                }

                if (wirePerWinding[windingIndex].get_type() == WireType::FOIL && !wirePerWinding[windingIndex].get_conducting_height()) {
                    wirePerWinding[windingIndex].cut_foil_wire_to_section(section);
                    get_mutable_functional_description()[windingIndex].set_wire(wirePerWinding[windingIndex]);
                }

                if (wirePerWinding[windingIndex].get_type() == WireType::PLANAR && !wirePerWinding[windingIndex].get_conducting_width()) {
                    wirePerWinding[windingIndex].cut_planar_wire_to_section(section);
                    get_mutable_functional_description()[windingIndex].set_wire(wirePerWinding[windingIndex]);
                }

                if (windingOrientation == WindingOrientation::OVERLAPPING) {

                    if ((resolve_margin(section)[0] + resolve_margin(section)[1] + resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_height())) > currentSectionHeight) {
                        std::string wireType = to_string(wirePerWinding[windingIndex].get_type());
                        return false;
                        // throw std::runtime_error("Margin plus a turn cannot larger than winding window" + 
                        //                          std::string{", margin:"} + std::to_string(resolve_margin(section)[0] + resolve_margin(section)[1]) + 
                        //                          ", wire type: " + wireType + 
                        //                          ", wire height: " + std::to_string(resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_height())) + 
                        //                          ", section height: " + std::to_string(currentSectionHeight)
                        //                          );
                    }
                }
                else {
                    if ((resolve_margin(section)[0] + resolve_margin(section)[1] + resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_width())) > currentSectionWidth) {
                        return false;
                        // throw std::runtime_error("Margin plus a turn cannot larger than winding window" + 
                        //                          std::string{", margin:"} + std::to_string(resolve_margin(section)[0] + resolve_margin(section)[1]) + 
                        //                          ", wire width:" + std::to_string(resolve_dimensional_values(wirePerWinding[windingIndex].get_maximum_outer_width())) + 
                        //                          ", section width:" + std::to_string(currentSectionWidth)
                        //                          );
                    }
                }

                if (section.get_dimensions()[0] < 0) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Something wrong happened in section dimensions 0: " + std::to_string(section.get_dimensions()[0]) +
                                             " availableWidth: " + std::to_string(availableWidth) +
                                             " currentSectionWidth: " + std::to_string(currentSectionWidth) +
                                             " currentSectionHeight: " + std::to_string(currentSectionHeight) + 
                                             " _marginsPerSection[marginIndex][0]: " + std::to_string(_marginsPerSection[marginIndex][0])
                                             );
                }
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_coordinates(std::vector<double>{currentSectionCenterWidth + currentSectionWidth / 2, currentSectionCenterHeight, 0});
                }
                else {
                    section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight - currentSectionHeight / 2, 0});
                }

                if (section.get_coordinates()[0] < -1) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Something wrong happened in section coordiantes 0: " + std::to_string(section.get_coordinates()[0]) +
                                             " currentSectionCenterWidth: " + std::to_string(currentSectionCenterWidth) +
                                             " group.get_coordinates()[0]: " + std::to_string(group.get_coordinates()[0]) +
                                             " group.get_dimensions()[0]: " + std::to_string(group.get_dimensions()[0]) +
                                             " availableWidth: " + std::to_string(availableWidth) +
                                             " currentSectionWidth: " + std::to_string(currentSectionWidth) +
                                             " currentSectionCenterHeight: " + std::to_string(currentSectionCenterHeight)
                                             );
                }

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / (currentSectionWidth * (currentSectionHeight - _marginsPerSection[marginIndex][0] - _marginsPerSection[marginIndex][1])));
                }
                else {
                    section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / ((currentSectionWidth - _marginsPerSection[marginIndex][0] - _marginsPerSection[marginIndex][1]) * currentSectionHeight));
                }
                section.set_winding_style(windByConsecutiveTurns[windingIndex]);
                sectionsDescription.push_back(section);

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingParallelsProportion[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
                }

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterWidth += currentSectionWidth;
                }
                else {
                    currentSectionCenterHeight -= currentSectionHeight;
                }

                currentSectionPerWinding[windingIndex]++;
            }
            else {
                if (sectionIndex == 0) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Insulation layers cannot be the first one (for now)");
                }

                auto previousWindingIndex = orderedSectionsWithInsulation[sectionIndex - 1].second.first;
                size_t nextWindingIndex;
                if (sectionIndex + 1 != orderedSectionsWithInsulation.size()) {
                    nextWindingIndex = orderedSectionsWithInsulation[sectionIndex + 1].second.first;
                }
                else {
                    nextWindingIndex = orderedSectionsWithInsulation[0].second.first;
                }

                auto windingsMapKey = std::pair<size_t, size_t>{previousWindingIndex, nextWindingIndex};
                if (!_insulationSections.contains(windingsMapKey)) {
                    log(_insulationSectionsLog[windingsMapKey]);
                    continue;
                }

                auto insulationSection = _insulationSections[windingsMapKey];

                insulationSection.set_group(group.get_name());
                if (groupWindowIndex) {
                    insulationSection.set_winding_window(static_cast<int64_t>(groupWindowIndex.value()));
                }
                insulationSection.set_name("Insulation between " + get_name(previousWindingIndex) + " and " + get_name(nextWindingIndex) + " section " + std::to_string(sectionIndex));
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth + insulationSection.get_dimensions()[0] / 2,
                                                                          currentSectionCenterHeight,
                                                                          0});
                }
                else {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth,
                                                                          currentSectionCenterHeight - insulationSection.get_dimensions()[1] / 2,
                                                                          0});
                }

                sectionsDescription.push_back(insulationSection);
                log(_insulationSectionsLog[windingsMapKey]);

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterWidth += insulationSection.get_dimensions()[0];
                }
                else {
                    currentSectionCenterHeight -= insulationSection.get_dimensions()[1];
                }

            }
        }
        marginSectionOffset += orderedSectionsWithInsulation.size();
    }

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            if (roundFloat(remainingParallelsProportion[windingIndex][parallelIndex], 9) != 0) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "There are unassigned parallel proportion in a rectangular section, something went wrong");
            }
        }
    }

    set_sections_description(sectionsDescription);
    return true;
}

void Coil::remove_insulation_if_margin_is_enough(std::vector<std::pair<size_t, double>> orderedSections) {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    size_t multiplier;
    if (_marginsPerSection.size() > orderedSections.size()) {
        multiplier = 2;
    }
    else {
        multiplier = 1;
    }

    for (size_t sectionIndex = 0; sectionIndex < orderedSections.size(); ++sectionIndex) {
        size_t indexForMarginLeftSection = sectionIndex * multiplier;
        size_t indexForMarginRightSection;
        if (sectionIndex + 1 != orderedSections.size()) {
            indexForMarginRightSection = (sectionIndex + 1) * multiplier;
        }
        else {
            indexForMarginRightSection = 0;
        }
        while (indexForMarginLeftSection >= _marginsPerSection.size() || indexForMarginRightSection >= _marginsPerSection.size()) {
            _marginsPerSection.push_back({0, 0});
        }
    }    

    for (size_t sectionIndex = 0; sectionIndex < orderedSections.size(); ++sectionIndex) {
        size_t indexForMarginLeftSection = sectionIndex * multiplier;
        // size_t indexForMarginLeftSection = sectionIndex;
        size_t indexForMarginRightSection;
        size_t leftWindingIndex = orderedSections[sectionIndex].first;
        size_t rightWindingIndex;
        if (sectionIndex + 1 != orderedSections.size()) {
            indexForMarginRightSection = (sectionIndex + 1) * multiplier;
            // indexForMarginRightSection = sectionIndex + 1;
            rightWindingIndex = orderedSections[sectionIndex + 1].first;
        }
        else {
            indexForMarginRightSection = 0;
            rightWindingIndex = orderedSections[0].first;
        }

        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex};
        double totalMargin = 0;
        if (_insulationSections.contains(windingsMapKey)) {
            auto coilSectionInterface = _coilSectionInterfaces[windingsMapKey];
            totalMargin = coilSectionInterface.get_total_margin_tape_distance();
        }

        if (_marginsPerSection.size() != 0) {
            double leftMargin = _marginsPerSection[indexForMarginLeftSection][1];
            double rightMargin = _marginsPerSection[indexForMarginRightSection][0];
            totalMargin = std::max(totalMargin, leftMargin + rightMargin);
        }

        double totalMarginAngle = wound_distance_to_angle(totalMargin, windingWindowRadialHeight);

        double totalInsulationDimension = 0;
        if (_insulationSections.contains(windingsMapKey)) {
            totalInsulationDimension = _insulationSections[windingsMapKey].get_dimensions()[1];

            if (totalMarginAngle > totalInsulationDimension) {
                _insulationSections[windingsMapKey].set_dimensions({_insulationSections[windingsMapKey].get_dimensions()[0], 0});
            }
        }
    }
}

bool Coil::wind_by_round_sections(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    set_sections_description(std::nullopt);
    std::vector<Section> sectionsDescription;

    if (!get_groups_description()) {
        throw CoilNotProcessedException("At least default group must be defined at this point.");
    }

    auto groups = get_groups_description().value();
    std::vector<std::vector<double>> remainingParallelsProportion;

    for (auto group : groups) {

        auto bobbin = resolve_bobbin();
        auto bobbinProcessedDescription = bobbin.get_processed_description().value();
        auto windingWindows = bobbinProcessedDescription.get_winding_windows();
        double availableRadialHeight = group.get_dimensions()[0];
        double availableAngle = group.get_dimensions()[1];

        double spaceForSections = 0;
        auto windingOrientation = get_winding_orientation();

        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            spaceForSections = availableRadialHeight;
        }
        else {
            spaceForSections = availableAngle;
        }

        auto orderedSections = get_ordered_sections(spaceForSections, proportionPerWinding, pattern, repetitions);

        if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            remove_insulation_if_margin_is_enough(orderedSections);
        }
        auto orderedSectionsWithInsulation = add_insulation_to_sections(orderedSections);

        double numberWindings = get_functional_description().size();
        auto numberSectionsPerWinding = std::vector<size_t>(numberWindings, 0);
        auto currentSectionPerWinding = std::vector<size_t>(numberWindings, 0);
        for (auto orderedSection : orderedSectionsWithInsulation) {
            if (orderedSection.first == ElectricalType::CONDUCTION) {
                auto windingIndex = orderedSection.second.first;
                numberSectionsPerWinding[windingIndex]++;
            }
        }
        // wound_with grouping (see rectangular sibling for full rationale):
        // center-tap partner winding shares its representative's section, so
        // the pattern omits its index — inherit the representative's count to
        // avoid the "Number of slots cannot be less than 1" throw on Path-B
        // rewind of toroidal center-tapped coils.
        for (size_t wIdx = 0; wIdx < numberWindings; ++wIdx) {
            if (numberSectionsPerWinding[wIdx] != 0) continue;
            const auto& wwOpt = get_functional_description()[wIdx].get_wound_with();
            if (!wwOpt || wwOpt->empty()) continue;
            for (const auto& partnerName : wwOpt.value()) {
                bool found = false;
                for (size_t pIdx = 0; pIdx < numberWindings; ++pIdx) {
                    if (get_functional_description()[pIdx].get_name() == partnerName &&
                        numberSectionsPerWinding[pIdx] > 0) {
                        numberSectionsPerWinding[wIdx] = numberSectionsPerWinding[pIdx];
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
        }
        // Fallback inheritance via _virtualizationMap (see rectangular sibling).
        if (!_virtualizationMap.empty()) {
            for (size_t wIdx = 0; wIdx < numberWindings; ++wIdx) {
                if (numberSectionsPerWinding[wIdx] != 0) continue;
                for (const auto& [virtualIdx, members] : _virtualizationMap) {
                    if (std::find(members.begin(), members.end(), wIdx) == members.end()) continue;
                    for (auto pIdx : members) {
                        if (pIdx < numberWindings && numberSectionsPerWinding[pIdx] > 0) {
                            numberSectionsPerWinding[wIdx] = numberSectionsPerWinding[pIdx];
                            break;
                        }
                    }
                    if (numberSectionsPerWinding[wIdx] != 0) break;
                }
            }
        }
        auto windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(), get_number_parallels(), numberSectionsPerWinding);
       
        auto wirePerWinding = get_wires();
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            remainingParallelsProportion.push_back(std::vector<double>(get_number_parallels(windingIndex), 1));
        }
        double currentSectionCenterAngle = DBL_MAX;
        double currentSectionCenterRadialHeight = DBL_MAX;

        apply_margin_tape(orderedSectionsWithInsulation);
        if (settings.get_coil_equalize_margins()) {
            equalize_margins(orderedSectionsWithInsulation);
        }

        std::vector<double> currentSectionRadialHeights;
        std::vector<double> currentSectionAngles;
        std::vector<size_t> windingIndexes;

        for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
            if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
                auto sectionInfo = orderedSectionsWithInsulation[sectionIndex].second;
                auto windingIndex = sectionInfo.first;
                auto aux = get_section_round_dimensions(orderedSectionsWithInsulation[sectionIndex], windingOrientation, availableRadialHeight, availableAngle);
                currentSectionRadialHeights.push_back(aux.first);
                currentSectionAngles.push_back(aux.second);
                windingIndexes.push_back(windingIndex);
            }
        }
        
        // std::vector<double> sectionLengths = get_section_lengths(currentSectionRadialHeights, currentSectionAngles, availableRadialHeight);

        std::vector<double> sectionPhysicalTurnsProportions;
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            std::vector<double> sectionAreas = get_section_areas(orderedSectionsWithInsulation, currentSectionAngles, availableRadialHeight);
            sectionPhysicalTurnsProportions = get_length_proportions(sectionAreas, windingIndexes);
        }
        else {
            sectionPhysicalTurnsProportions = std::vector<double>(orderedSectionsWithInsulation.size(), 1);
        }

        size_t conductingSectionIndex = 0;
        for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
            if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
                auto sectionInfo = orderedSectionsWithInsulation[sectionIndex].second;
                auto windingIndex = sectionInfo.first;

                double currentSectionRadialHeight = currentSectionRadialHeights[conductingSectionIndex];
                double currentSectionAngle = currentSectionAngles[conductingSectionIndex];

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    if (currentSectionCenterRadialHeight == DBL_MAX) {
                        currentSectionCenterRadialHeight = 0;
                    }
                    if (currentSectionCenterAngle == DBL_MAX) {
                        currentSectionCenterAngle = 180;
                    }
                }
                else {
                    if (currentSectionCenterRadialHeight == DBL_MAX) {
                        currentSectionCenterRadialHeight = 0;
                    }
                    if (currentSectionCenterAngle == DBL_MAX) {
                        currentSectionCenterAngle = 0;
                    }
                }
                
                PartialWinding partialWinding;
                Section section;
                partialWinding.set_winding(get_name(windingIndex));

                auto parallelsProportions = get_parallels_proportions(currentSectionPerWinding[windingIndex],
                                                                       numberSectionsPerWinding[windingIndex],
                                                                       get_number_turns(windingIndex),
                                                                       get_number_parallels(windingIndex),
                                                                       remainingParallelsProportion[windingIndex],
                                                                       windByConsecutiveTurns[windingIndex],
                                                                       std::vector<double>(get_number_parallels(windingIndex), 1),
                                                                       sectionPhysicalTurnsProportions[conductingSectionIndex]);

                std::vector<double> sectionParallelsProportion = parallelsProportions.second;
                uint64_t physicalTurnsThisSection = parallelsProportions.first;

                double marginAngle0 = 0;
                double marginAngle1 = 0;
                size_t numberLayers = ULONG_MAX;
                size_t prevNumberLayers = 0;

                // We correct the radial height to exactly what we need, so afterwards we can calculate exactly how many turns we need
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    auto aux = get_number_layers_needed_and_number_physical_turns(currentSectionCenterRadialHeight + _marginsPerSection[sectionIndex][0], currentSectionAngle, wirePerWinding[windingIndex], physicalTurnsThisSection, availableRadialHeight);
                    numberLayers = aux.first;
                    currentSectionRadialHeight = numberLayers * wirePerWinding[windingIndex].get_maximum_outer_width();

                    if (_insulationInterLayers.contains(windingIndex)) {
                        auto insulationLayer = _insulationInterLayers[windingIndex];
                        currentSectionRadialHeight += (numberLayers - 1) * insulationLayer.get_dimensions()[0];
                    }
                }
                else {
                    while (numberLayers != prevNumberLayers) {
                        prevNumberLayers = numberLayers;
                        double currentSectionAngleMinusMargin = currentSectionAngle - marginAngle0 - marginAngle1;
                        auto aux = get_number_layers_needed_and_number_physical_turns(currentSectionCenterRadialHeight, currentSectionAngleMinusMargin, wirePerWinding[windingIndex], physicalTurnsThisSection, availableRadialHeight);
                        numberLayers = aux.first;
                        if (_strict) {
                            currentSectionRadialHeight = numberLayers * wirePerWinding[windingIndex].get_maximum_outer_width();
                        }
                        double lastLayerMaximumRadius = availableRadialHeight - (currentSectionCenterRadialHeight + numberLayers * wirePerWinding[windingIndex].get_maximum_outer_width());
                        if (lastLayerMaximumRadius < 0) {
                            break;
                        }
                        marginAngle0 = wound_distance_to_angle(_marginsPerSection[sectionIndex][0], lastLayerMaximumRadius);
                        marginAngle1 = wound_distance_to_angle(_marginsPerSection[sectionIndex][1], lastLayerMaximumRadius);
                    }                
                    currentSectionAngle -= marginAngle0 + marginAngle1;
                }


                // if (currentSectionRadialHeight > availableRadialHeight) {
                //     return false;
                // }
                if (currentSectionAngle < 0) {
                    return false;
                }

                partialWinding.set_parallels_proportion(sectionParallelsProportion);
                section.set_name(get_name(windingIndex) +  " section " + std::to_string(currentSectionPerWinding[windingIndex]));
                section.set_partial_windings(std::vector<PartialWinding>{partialWinding});  // TODO: support more than one winding per section?
                section.set_type(ElectricalType::CONDUCTION);
                section.set_group(group.get_name());
                section.set_margin(_marginsPerSection[sectionIndex]);
                section.set_layers_orientation(_layersOrientation);
                section.set_coordinate_system(CoordinateSystem::POLAR);
                
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{currentSectionRadialHeight, currentSectionAngle});
                    section.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + currentSectionRadialHeight / 2 + _marginsPerSection[sectionIndex][0], currentSectionCenterAngle, 0});
                }
                else {
                    section.set_dimensions(std::vector<double>{currentSectionRadialHeight, currentSectionAngle});
                    section.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + currentSectionRadialHeight / 2, currentSectionCenterAngle  + currentSectionAngle / 2 + marginAngle0, 0});
                }

                if (section.get_dimensions()[0] < 0) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Something wrong happened in section dimensions 0: " + std::to_string(section.get_dimensions()[0]) +
                                             " currentSectionRadialHeight: " + std::to_string(currentSectionRadialHeight) +
                                             " currentSectionAngle: " + std::to_string(currentSectionAngle)
                                             );
                }


                if (section.get_dimensions()[1] < 0) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Something wrong happened in section dimensions 1: " + std::to_string(section.get_dimensions()[1]) +
                                             " currentSectionRadialHeight: " + std::to_string(currentSectionRadialHeight) +
                                             " currentSectionAngle: " + std::to_string(currentSectionAngle)
                                             );
                }
                
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    double ringArea = std::numbers::pi * pow(availableRadialHeight - currentSectionCenterRadialHeight, 2) - std::numbers::pi * pow(availableRadialHeight - (currentSectionCenterRadialHeight + currentSectionRadialHeight), 2);

                    section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / ringArea);
                } 
                else {
                    double ringArea = std::numbers::pi * pow(availableRadialHeight - currentSectionCenterRadialHeight, 2) * currentSectionAngle / 360;
                    section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / ringArea);
                }
                section.set_winding_style(windByConsecutiveTurns[windingIndex]);
                sectionsDescription.push_back(section);

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingParallelsProportion[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
                }

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterRadialHeight += currentSectionRadialHeight + _marginsPerSection[sectionIndex][0] + _marginsPerSection[sectionIndex][1];
                }
                else {
                    currentSectionCenterAngle += currentSectionAngle + marginAngle0 + marginAngle1;
                }

                currentSectionPerWinding[windingIndex]++;
                conductingSectionIndex++;
            }
            else {
                if (sectionIndex == 0) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Insulation layers cannot be the first one (for now)");
                }
     
                auto previousWindingIndex = orderedSectionsWithInsulation[sectionIndex - 1].second.first;
                size_t nextWindingIndex;
                if (sectionIndex + 1 != orderedSectionsWithInsulation.size()) {
                    nextWindingIndex = orderedSectionsWithInsulation[sectionIndex + 1].second.first;
                }
                else {
                    nextWindingIndex = orderedSectionsWithInsulation[0].second.first;
                }

                auto windingsMapKey = std::pair<size_t, size_t>{previousWindingIndex, nextWindingIndex};
                if (!_insulationSections.contains(windingsMapKey)) {
                    log(_insulationSectionsLog[windingsMapKey]);
                    continue;
                }

                auto insulationSection = _insulationSections[windingsMapKey];

                insulationSection.set_group(group.get_name());
                insulationSection.set_name("Insulation between " + get_name(previousWindingIndex) + " and " + get_name(nextWindingIndex) + " section " + std::to_string(sectionIndex));
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + insulationSection.get_dimensions()[0] / 2,
                                                                          currentSectionCenterAngle,
                                                                          0});
                }
                else {
                    insulationSection.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + insulationSection.get_dimensions()[0] / 2,
                                                                          currentSectionCenterAngle - insulationSection.get_dimensions()[1] / 2,
                                                                          0});
                }

                sectionsDescription.push_back(insulationSection);
                log(_insulationSectionsLog[windingsMapKey]);

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterRadialHeight += insulationSection.get_dimensions()[0];
                }
                else {
                    currentSectionCenterAngle += insulationSection.get_dimensions()[1];
                }
            }
        }
    }


    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            if (roundFloat(remainingParallelsProportion[windingIndex][parallelIndex], 9) != 0) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "There are unassigned parallel proportion in a round section, something went wrong");
            }
        }
    }

    set_sections_description(sectionsDescription);
    return true;
}

bool Coil::wind_by_planar_sections(std::vector<size_t> stackUpForThisGroup, std::map<std::pair<size_t, size_t>, double> insulationThickness, double coreToLayerDistance) {
    // In planar coils each section will have only one layer
    set_layers_description(std::nullopt);
    std::vector<Section> sections;

    auto bobbin = resolve_bobbin();
    if (!get_groups_description()) {
        create_default_groups(bobbin, WiringTechnology::PRINTED, coreToLayerDistance);
    }

    if (!get_groups_description()) {
        throw CoilNotProcessedException("At least default group must be defined at this point.");
    }

    set_winding_orientation(WindingOrientation::CONTIGUOUS);
    auto sectionAlignment = get_section_alignment();
    set_section_alignment(sectionAlignment);
    set_turns_alignment(CoilAlignment::CENTERED);

    auto groups = get_groups_description().value();
    // Planar coils v1: only the first group is wound here. Multi-column
    // planar layouts are uncommon; defer to v2 if needed.
    auto group = groups[0];

    auto wirePerWinding = get_wires();
    if (wirePerWinding.size() == 0) {
        throw WireNotFoundException("Wires missing");
    }

    std::vector<size_t> numberSectionsPerWinding = std::vector<size_t>(wirePerWinding.size(), 0);
    std::vector<std::vector<double>> totalParallelsProportionPerWinding;
    std::vector<std::vector<double>> remainingParallelsProportionPerWinding;
    std::vector<size_t> currentSectionIndexPerwinding = std::vector<size_t>(wirePerWinding.size(), 0);

    for (auto windingIndex : stackUpForThisGroup) {
        numberSectionsPerWinding[windingIndex]++;
    }

    for (auto winding : group.get_partial_windings()) {
        totalParallelsProportionPerWinding.push_back(winding.get_parallels_proportion());
        remainingParallelsProportionPerWinding.push_back(winding.get_parallels_proportion());
    }

    std::vector<double> sectionWidthPerWinding;
    std::vector<double> sectionHeightPerWinding;
    double totalSectionHeight = 0;
    double totalAvailableHeight = group.get_dimensions()[1];
    for (size_t stackUpIndex = 0; stackUpIndex < stackUpForThisGroup.size(); ++stackUpIndex) {
        if (stackUpIndex + 1 < stackUpForThisGroup.size()) {
            std::pair key = {stackUpForThisGroup[stackUpIndex], stackUpForThisGroup[stackUpIndex + 1]};
            if (insulationThickness.count(key)) {
                totalAvailableHeight -= insulationThickness[key];
            }
            else {
                totalAvailableHeight -= defaults.pcbInsulationThickness;
            }
        }
    }

    for (size_t stackUpIndex = 0; stackUpIndex < stackUpForThisGroup.size(); ++stackUpIndex) {
        sectionWidthPerWinding.push_back(group.get_dimensions()[0]);
        // double sectionHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
        double sectionHeight = totalAvailableHeight / stackUpForThisGroup.size();
        sectionHeightPerWinding.push_back(sectionHeight);
        totalSectionHeight += sectionHeight;

        if (stackUpIndex + 1 < stackUpForThisGroup.size()) {
            std::pair key = {stackUpForThisGroup[stackUpIndex], stackUpForThisGroup[stackUpIndex + 1]};
            if (insulationThickness.count(key)) {
                totalSectionHeight += insulationThickness[key];
            }
            else {
                totalSectionHeight += defaults.pcbInsulationThickness;
            }
        }
    }
    double currentSectionCenterWidth = roundFloat(group.get_coordinates()[0], 9);
    double windowHeight = group.get_dimensions()[1];
    double groupCenterY = group.get_coordinates()[1];
    double currentSectionCenterHeight;
    
    // Calculate starting Y position based on section alignment
    // Sections are placed from top to bottom
    switch (sectionAlignment) {
        case CoilAlignment::OUTER_OR_BOTTOM:
            // Place sections at the bottom of the window
            // Start higher so the stack ends at the bottom
            currentSectionCenterHeight = roundFloat(groupCenterY + windowHeight/2 - totalSectionHeight + sectionHeightPerWinding[0]/2, 9);
            break;
        case CoilAlignment::CENTERED:
            // Center the stack in the window
            currentSectionCenterHeight = roundFloat(groupCenterY - totalSectionHeight/2 + sectionHeightPerWinding[0]/2, 9);
            break;
        case CoilAlignment::SPREAD:
            // Spread sections to fill the window (like centered but using full height)
            currentSectionCenterHeight = roundFloat(groupCenterY + totalSectionHeight/2, 9);
            break;
        case CoilAlignment::INNER_OR_TOP:
        default:
            // Start from top (original behavior)
            currentSectionCenterHeight = roundFloat(groupCenterY + totalSectionHeight/2, 9);
            break;
    }

    for (size_t stackUpIndex = 0; stackUpIndex < stackUpForThisGroup.size(); ++stackUpIndex) {
        Section section;
        auto windingIndex = stackUpForThisGroup[stackUpIndex];
        auto remainingParallelsProportionInWinding = remainingParallelsProportionPerWinding[windingIndex];
        auto totalParallelsProportionInWinding = totalParallelsProportionPerWinding[windingIndex];
        auto numberSections = numberSectionsPerWinding[windingIndex];
        auto winding = get_functional_description()[windingIndex];
        auto sectionIndex = currentSectionIndexPerwinding[windingIndex];
        double sectionWidth = sectionWidthPerWinding[windingIndex];
        double sectionHeight = sectionHeightPerWinding[windingIndex];
        currentSectionCenterHeight -= sectionHeight / 2;

        WindingStyle windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberSections);

        auto parallelsProportions = get_parallels_proportions(sectionIndex,
                                                               numberSections,
                                                               get_number_turns(windingIndex),
                                                               get_number_parallels(windingIndex),
                                                               remainingParallelsProportionInWinding,
                                                               windByConsecutiveTurns,
                                                               totalParallelsProportionInWinding);

        std::vector<double> sectionParallelsProportion = parallelsProportions.second;

        size_t numberParallelsProportionsToZero = 0;
        for (auto parallelProportion : sectionParallelsProportion) {
            if (parallelProportion == 0) {
                numberParallelsProportionsToZero++;
            }
        }

        if (numberParallelsProportionsToZero == sectionParallelsProportion.size()) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Parallel proportion in section cannot be all be 0");
        }

        uint64_t physicalTurnsThisSection = parallelsProportions.first;

        auto partialWinding = group.get_partial_windings()[windingIndex];
        partialWinding.set_parallels_proportion(sectionParallelsProportion);
        section.set_partial_windings(std::vector<PartialWinding>{partialWinding});
        section.set_group(group.get_name());
        section.set_type(ElectricalType::CONDUCTION);
        section.set_name(winding.get_name() + " section " + std::to_string(sectionIndex));
        section.set_layers_orientation(WindingOrientation::CONTIGUOUS);
        section.set_dimensions(std::vector<double>{sectionWidth, sectionHeight});
        section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
        section.set_coordinate_system(CoordinateSystem::CARTESIAN);
        section.set_margin(std::vector<double>{0, 0});


        section.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisSection) / (sectionWidth * sectionHeight));
        section.set_winding_style(windByConsecutiveTurns);
        sections.push_back(section);

        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            remainingParallelsProportionPerWinding[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
        }

        currentSectionCenterHeight = roundFloat(currentSectionCenterHeight -= sectionHeight / 2, 9);
        currentSectionIndexPerwinding[windingIndex]++;

        if (stackUpIndex + 1 < stackUpForThisGroup.size()) {
            std::pair key = {stackUpForThisGroup[stackUpIndex], stackUpForThisGroup[stackUpIndex + 1]};
            double insulationThicknessThisLayer;
            if (insulationThickness.count(key)) {
                insulationThicknessThisLayer = insulationThickness[key];
            }
            else {
                insulationThicknessThisLayer = defaults.pcbInsulationThickness;
            }
            currentSectionCenterHeight -= insulationThicknessThisLayer / 2;

            Section insulationSection;
            insulationSection.set_type(ElectricalType::INSULATION);
            insulationSection.set_name("Insulation section between stack index " + std::to_string(stackUpIndex) + " and " + std::to_string(stackUpIndex + 1));
            insulationSection.set_dimensions(std::vector<double>{sectionWidth, insulationThicknessThisLayer});
            insulationSection.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            insulationSection.set_coordinate_system(CoordinateSystem::CARTESIAN);
            insulationSection.set_layers_orientation(WindingOrientation::CONTIGUOUS);
            insulationSection.set_filling_factor(1);
            insulationSection.set_margin(std::vector<double>{0, 0});
            sections.push_back(insulationSection);
            currentSectionCenterHeight -= insulationThicknessThisLayer / 2;
        }

    }
    set_sections_description(sections);
    return true;
}

bool Coil::wind_by_layers() {
    set_layers_description(std::nullopt);
    if (!get_sections_description()) {
        return false;
    }
    auto bobbin = resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    auto functionalDescription = get_functional_description();
    auto sectionsDescription = get_sections_description().value();
    auto insulationInterSectionsLayers = _insulationInterSectionsLayers;
    auto needsVirtualization = needs_virtualization();

    if (needsVirtualization) {
        create_virtualization_map();
        auto virtualFunctionalDescription = virtualize_functional_description();
        auto virtualSectionsDescription = virtualize_sections_description();
        _insulationInterSectionsLayers = virtualize_insulation_intersections_layers();
        set_functional_description(virtualFunctionalDescription);
        set_sections_description(virtualSectionsDescription);
        _windingIndexByName.clear();
        _turnIndexByName.clear();
    }

    bool result;
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        result = wind_by_rectangular_layers();
    }
    else {
        result = wind_by_round_layers();
    }

    if (needsVirtualization) {
        _windingIndexByName.clear();
        _turnIndexByName.clear();
        set_functional_description(functionalDescription);
        set_sections_description(sectionsDescription);
        _insulationInterSectionsLayers = insulationInterSectionsLayers;
        devirtualize_layers_description();
    }

    return result;
}

bool Coil::wind_by_rectangular_layers() {
    set_layers_description(std::nullopt);
    if (!get_sections_description()) {
        return false;
    }

    auto wirePerWinding = get_wires();

    auto sections = get_sections_description().value();

    std::vector<Layer> layers;
    for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
        if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {

            uint64_t maximumNumberLayersFittingInSection;
            uint64_t maximumNumberPhysicalTurnsPerLayer;
            uint64_t minimumNumberLayerNeeded;
            uint64_t numberLayers;
            uint64_t physicalTurnsInSection = 0;
            double layerWidth = 0;
            double layerHeight = 0;
            auto remainingParallelsProportionInSection = sections[sectionIndex].get_partial_windings()[0].get_parallels_proportion();
            auto totalParallelsProportionInSection = sections[sectionIndex].get_partial_windings()[0].get_parallels_proportion();
            if (sections[sectionIndex].get_partial_windings().size() > 1) {
                throw NotImplementedException("More than one winding per layer not supported yet");
            }
            auto partialWinding = sections[sectionIndex].get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                physicalTurnsInSection += round(remainingParallelsProportionInSection[parallelIndex] * get_number_turns(windingIndex));
            }

            if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                if (!wirePerWinding[windingIndex].get_outer_diameter()) {
                    throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing wire outer diameter");
                }
                double wireDiameter = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_diameter().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[0] / wireDiameter;
                    maximumNumberPhysicalTurnsPerLayer = floor(sections[sectionIndex].get_dimensions()[1] / wireDiameter);
                    layerWidth = wireDiameter;
                    layerHeight = sections[sectionIndex].get_dimensions()[1];
                } else {
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[1] / wireDiameter;
                    maximumNumberPhysicalTurnsPerLayer = floor(sections[sectionIndex].get_dimensions()[0] / wireDiameter);
                    layerWidth = sections[sectionIndex].get_dimensions()[0];
                    layerHeight = wireDiameter;
                }
            }
            else {
                if (!wirePerWinding[windingIndex].get_outer_width()) {
                    throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing wire outer width");
                }
                if (!wirePerWinding[windingIndex].get_outer_height()) {
                    throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing wire outer height");
                }
                double wireWidth = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_width().value());
                double wireHeight = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_height().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[0] / wireWidth;
                    if (wirePerWinding[windingIndex].get_type() == WireType::FOIL) {
                        maximumNumberPhysicalTurnsPerLayer = 1;
                    }
                    else {
                        maximumNumberPhysicalTurnsPerLayer = floor(sections[sectionIndex].get_dimensions()[1] / wireHeight);
                    }
                    layerWidth = wireWidth;
                    layerHeight = sections[sectionIndex].get_dimensions()[1];
                } else {
                    maximumNumberLayersFittingInSection = sections[sectionIndex].get_dimensions()[1] / wireHeight;
                    if (wirePerWinding[windingIndex].get_type() == WireType::RECTANGULAR && settings.get_coil_only_one_turn_per_layer_in_contiguous_rectangular()) {
                        maximumNumberPhysicalTurnsPerLayer = 1; 
                    }
                    else {
                        maximumNumberPhysicalTurnsPerLayer = floor(sections[sectionIndex].get_dimensions()[0] / wireWidth);
                    }
                    layerWidth = sections[sectionIndex].get_dimensions()[0];
                    layerHeight = wireHeight;
                }
            }

            if (sections[sectionIndex].get_number_layers()) {
                numberLayers = sections[sectionIndex].get_number_layers().value();
            }
            else {
                if (maximumNumberLayersFittingInSection == 0) {
                    numberLayers = ceil(double(physicalTurnsInSection) / maximumNumberPhysicalTurnsPerLayer);
                }
                else if (maximumNumberPhysicalTurnsPerLayer == 0) {
                    numberLayers = maximumNumberLayersFittingInSection;
                }
                else {
                    minimumNumberLayerNeeded = ceil(double(physicalTurnsInSection) / maximumNumberPhysicalTurnsPerLayer);
                    numberLayers = std::min(minimumNumberLayerNeeded, maximumNumberLayersFittingInSection);
                }
            }

            // We cannot have more layers than physical turns
            numberLayers = std::min(numberLayers, physicalTurnsInSection);

            // Real winding geometry: turn blocking is GLOBAL to the winding window. wind() computes,
            // per conduction layer, how many connection leads cross its top and bottom
            // (_connectionBlockedSlotsPerLayer, keyed by layer name) and re-winds with
            // _applyConnectionBlocking set. Here we size each layer around its blocked slots —
            // iterating layer by layer with the reduced capacity of each — so the turns are known to
            // fit before placement and the section simply grows by extra layers. Works for any
            // winding (single or interleaved); a lead blocks whatever layer it crosses regardless of
            // section. Overlapping, unconstrained-layer-count only. Works for any parallel count: a
            // bifilar/N-filar winding lays its parallels side by side, so each layer holds an equal
            // number of turns of every parallel (its physical capacity is rounded down to a multiple of
            // the parallel count) and the per-layer split is carried by CONSECUTIVE_PARALLELS below.
            int64_t numberParallels = int64_t(get_number_parallels(windingIndex));
            bool realWindingBlocking = settings.get_coil_use_real_winding_geometry()
                && _applyConnectionBlocking
                && sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING
                && !sections[sectionIndex].get_number_layers()
                && maximumNumberPhysicalTurnsPerLayer > 1
                && maximumNumberPhysicalTurnsPerLayer >= uint64_t(numberParallels);
            double wireAxialSize = wirePerWinding[windingIndex].get_maximum_outer_height();
            auto blockedSlotsForLayer = [&](size_t layerIndexInSection) -> std::pair<uint64_t, uint64_t> {
                auto found = _connectionBlockedSlotsPerLayer.find(
                    sections[sectionIndex].get_name() + " layer " + std::to_string(layerIndexInSection));
                if (found == _connectionBlockedSlotsPerLayer.end()) {
                    return {0u, 0u};
                }
                return found->second;
            };
            std::vector<uint64_t> realPerLayerTurns;
            if (realWindingBlocking) {
                uint64_t remainingTurns = physicalTurnsInSection;
                size_t builtLayers = 0;
                while (remainingTurns > 0) {
                    auto blocked = blockedSlotsForLayer(builtLayers);
                    uint64_t blockedSlots = std::min<uint64_t>(blocked.first + blocked.second, maximumNumberPhysicalTurnsPerLayer - 1);
                    uint64_t capacity = maximumNumberPhysicalTurnsPerLayer - blockedSlots;
                    // Side-by-side parallels: hold an equal number of turns of every parallel, so the
                    // layer capacity is a whole number of parallel rows (a multiple of the parallel
                    // count). K=1 leaves it unchanged.
                    capacity = (capacity / uint64_t(numberParallels)) * uint64_t(numberParallels);
                    if (capacity == 0) {
                        break;  // layer too thin to hold one row of every parallel
                    }
                    uint64_t turnsThisLayer = std::min<uint64_t>(capacity, remainingTurns);
                    realPerLayerTurns.push_back(turnsThisLayer);
                    remainingTurns -= turnsThisLayer;
                    builtLayers++;
                }
                numberLayers = realPerLayerTurns.size();
            }

            auto turnsAlignment = get_turns_alignment(sections[sectionIndex].get_name());
            // Real winding geometry: turns stay tightly packed; the connection leads route along the
            // window edges and align_blocked_layer_turns() (called after delimit) shifts each blocked
            // layer's turns to the unblocked edge so the freed slots sit exactly where the leads are
            // (and no window space is wasted). So no per-section alignment override is needed here.

            double currentLayerCenterWidth;
            double currentLayerCenterHeight;
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerWidth / 2, 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0], 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);

                if (turnsAlignment == CoilAlignment::CENTERED || turnsAlignment == CoilAlignment::SPREAD) {
                    currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + (numberLayers * layerHeight) / 2 - layerHeight / 2, 9);
                }
                else if (turnsAlignment == CoilAlignment::INNER_OR_TOP) {
                    currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);
                }
                else if (turnsAlignment == CoilAlignment::OUTER_OR_BOTTOM) {
                    currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] - sections[sectionIndex].get_dimensions()[1] / 2 + (numberLayers * layerHeight) - layerHeight / 2, 9);
                }
                else {
                    throw std::invalid_argument("Unknown turns alignment");
                }
            }

            WindingStyle windByConsecutiveTurns;
            if (sections[sectionIndex].get_winding_style()) {
                windByConsecutiveTurns = sections[sectionIndex].get_winding_style().value();
            }
            else {
                windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberLayers);
            }

            if (windByConsecutiveTurns == WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS && maximumNumberPhysicalTurnsPerLayer < get_number_parallels(windingIndex)) {
                windByConsecutiveTurns = WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
            }

            // Real-winding blocking forces an exact per-layer physical-turn count (realPerLayerTurns)
            // through get_parallels_proportions' slotAbsolutePhysicalTurns, honored by BOTH the
            // WIND_BY_CONSECUTIVE_TURNS branch and (now) the WIND_BY_CONSECUTIVE_PARALLELS branch. For a
            // single parallel the two styles place identically, so we pin CONSECUTIVE_TURNS (its
            // long-standing path). For a bifilar/N-filar winding we KEEP CONSECUTIVE_PARALLELS so the
            // parallels are laid side by side (better current sharing / proximity) and the forced count
            // is split evenly across them — forcing TURNS here would change the electrical layout.
            if (realWindingBlocking && numberParallels == 1) {
                windByConsecutiveTurns = WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
            }

            for (size_t layerIndex = 0; layerIndex < numberLayers; ++layerIndex) {
                Layer layer;

                auto parallelsProportions = realWindingBlocking
                    ? get_parallels_proportions(layerIndex, numberLayers, get_number_turns(windingIndex), get_number_parallels(windingIndex),
                                                remainingParallelsProportionInSection, windByConsecutiveTurns, totalParallelsProportionInSection,
                                                1.0, double(realPerLayerTurns[layerIndex]))
                    : get_parallels_proportions(layerIndex, numberLayers, get_number_turns(windingIndex), get_number_parallels(windingIndex),
                                                remainingParallelsProportionInSection, windByConsecutiveTurns, totalParallelsProportionInSection);

                std::vector<double> layerParallelsProportion = parallelsProportions.second;

                size_t numberParallelsProportionsToZero = 0;
                for (auto parallelProportion : layerParallelsProportion) {
                    if (parallelProportion == 0) {
                        numberParallelsProportionsToZero++;
                    }
                }

                if (numberParallelsProportionsToZero == layerParallelsProportion.size()) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Parallel proportion in layer cannot be all be 0");
                }

                uint64_t physicalTurnsThisLayer = parallelsProportions.first;

                partialWinding.set_parallels_proportion(layerParallelsProportion);
                layer.set_partial_windings(std::vector<PartialWinding>{partialWinding});
                layer.set_section(sections[sectionIndex].get_name());
                layer.set_type(ElectricalType::CONDUCTION);
                layer.set_name(sections[sectionIndex].get_name() +  " layer " + std::to_string(layerIndex));
                layer.set_orientation(sections[sectionIndex].get_layers_orientation());
                // Real winding geometry: SPREAD every layer's turns over its (blocking-reduced) height.
                // A full layer (turns == capacity) spreads with ~zero gap so it looks filled; a partial
                // (spillover) layer distributes its few turns evenly. The reduced/shifted height already
                // frees the blocked edges, so no separate edge-packing is needed.
                auto layerTurnsAlignment = turnsAlignment;
                if (realWindingBlocking) {
                    layerTurnsAlignment = CoilAlignment::SPREAD;
                }
                layer.set_turns_alignment(layerTurnsAlignment);
                // In real winding geometry, each layer loses the turn slots blocked by connection
                // leads crossing it: shrink the layer height by the blocked slots and shift it away
                // from those slots (top slots push it down, bottom slots push it up) so the turns
                // stop before the connection area, leaving room for the leads.
                double thisLayerHeight = layerHeight;
                double thisLayerCenterHeight = currentLayerCenterHeight;
                if (realWindingBlocking) {
                    auto blocked = blockedSlotsForLayer(layerIndex);
                    uint64_t blockedSlots = std::min<uint64_t>(blocked.first + blocked.second, maximumNumberPhysicalTurnsPerLayer - 1);
                    thisLayerHeight = roundFloat(layerHeight - blockedSlots * wireAxialSize, 9);
                    thisLayerCenterHeight = roundFloat(currentLayerCenterHeight
                        + (double(blocked.second) - double(blocked.first)) * wireAxialSize / 2, 9);
                }
                layer.set_dimensions(std::vector<double>{layerWidth, thisLayerHeight});
                layer.set_coordinates(std::vector<double>{currentLayerCenterWidth, thisLayerCenterHeight, 0});
                layer.set_coordinate_system(CoordinateSystem::CARTESIAN);

                layer.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisLayer) / (layerWidth * thisLayerHeight));
                layer.set_winding_style(windByConsecutiveTurns);
                layers.push_back(layer);

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingParallelsProportionInSection[parallelIndex] -= layerParallelsProportion[parallelIndex];
                }

                if (layerIndex == numberLayers - 1) {
                    break;
                }

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                    currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - layerHeight / 2, 9);
                }
                else {
                    currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + layerWidth / 2, 9);
                }

                if (_insulationInterLayers.contains(windingIndex)) {
     
                    auto insulationLayer = _insulationInterLayers[windingIndex];
                    if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                        currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - insulationLayer.get_dimensions()[1] / 2, 9);
                    }
                    else {
                        currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + insulationLayer.get_dimensions()[0] / 2, 9);
                    }

                    insulationLayer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                    insulationLayer.set_section(sections[sectionIndex].get_name());
                    insulationLayer.set_name(sections[sectionIndex].get_name() +  " insulation layer " + std::to_string(layerIndex));
                    insulationLayer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                    if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                        insulationLayer.set_dimensions(std::vector<double>{layerWidth, insulationLayer.get_dimensions()[1]});
                    }
                    else {
                        insulationLayer.set_dimensions(std::vector<double>{insulationLayer.get_dimensions()[0], layerHeight});
                    }
                    layers.push_back(insulationLayer);

                    if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                        currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - insulationLayer.get_dimensions()[1] / 2, 9);
                    }
                    else {
                        currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + insulationLayer.get_dimensions()[0] / 2, 9);
                    }

                }

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                    currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - layerHeight / 2, 9);
                }
                else {
                    currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + layerWidth / 2, 9);
                }
            }

        }
        else {
            if (sectionIndex == 0) {
                throw NotImplementedException("inner insulation layers not implemented");
            }

            auto partialWinding = sections[sectionIndex - 1].get_partial_windings()[0];
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            Section nextSection;
            if (sectionIndex + 1 != sections.size()) {
                if (sections[sectionIndex - 1].get_type() != ElectricalType::CONDUCTION || sections[sectionIndex + 1].get_type() != ElectricalType::CONDUCTION) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Previous and next sections must be conductive");
                }
                nextSection = sections[sectionIndex + 1];
            }
            else {
                nextSection = sections[0];
            }
            // auto nextSection = sections[sectionIndex + 1];
            auto nextPartialWinding = nextSection.get_partial_windings()[0];
            auto nextWindingIndex = get_winding_index_by_name(nextPartialWinding.get_winding());

            auto windingsMapKey = std::pair<size_t, size_t>{windingIndex, nextWindingIndex};
            if (!_insulationInterSectionsLayers.contains(windingsMapKey)) {
                log(_insulationInterSectionsLayersLog[windingsMapKey]);
                continue;
            }
            auto insulationLayers = _insulationInterSectionsLayers[windingsMapKey];

            if (insulationLayers.size() == 0) {
                continue;
                // throw std::runtime_error("There must be at least one insulation layer between layers");
            }

            double layerWidth = insulationLayers[0].get_dimensions()[0];
            double layerHeight = insulationLayers[0].get_dimensions()[1];

            double currentLayerCenterWidth;
            double currentLayerCenterHeight;
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerWidth / 2, 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                currentLayerCenterWidth = roundFloat(sections[sectionIndex].get_coordinates()[0], 9);
                currentLayerCenterHeight = roundFloat(sections[sectionIndex].get_coordinates()[1] + sections[sectionIndex].get_dimensions()[1] / 2 - layerHeight / 2, 9);
            }

            for (size_t layerIndex = 0; layerIndex < insulationLayers.size(); ++layerIndex) {
                auto insulationLayer = insulationLayers[layerIndex];
                insulationLayer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                insulationLayer.set_section(sections[sectionIndex].get_name());
                insulationLayer.set_name(sections[sectionIndex].get_name() +  " insulation layer " + std::to_string(layerIndex));
                insulationLayer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                layers.push_back(insulationLayer);

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                    currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - layerHeight, 9);
                }
                else {
                    currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + layerWidth, 9);
                }
            }
        }
    }
    set_layers_description(layers);
    return true;
}

bool Coil::wind_by_round_layers() {
    set_layers_description(std::nullopt);
    if (!get_sections_description()) {
        return false;
    }
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    auto wirePerWinding = get_wires();

    auto sections = get_sections_description().value();

    std::vector<Layer> layers;
    for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
        if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
            uint64_t maximumNumberLayersFittingInSection;
            uint64_t maximumNumberPhysicalTurnsPerLayer;
            uint64_t minimumNumberLayerNeeded = 0;
            uint64_t numberLayers;
            std::vector<int64_t> layerPhysicalTurns;
            uint64_t physicalTurnsInSection = 0;
            double layerRadialHeight = 0;
            double layerAngle = 0;
            auto remainingParallelsProportionInSection = sections[sectionIndex].get_partial_windings()[0].get_parallels_proportion();
            auto totalParallelsProportionInSection = sections[sectionIndex].get_partial_windings()[0].get_parallels_proportion();
            if (sections[sectionIndex].get_partial_windings().size() > 1) {
                throw NotImplementedException("More than one winding per layer not supported yet");
            }
            auto partialWinding = sections[sectionIndex].get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                physicalTurnsInSection += round(remainingParallelsProportionInSection[parallelIndex] * get_number_turns(windingIndex));
            }

            if (std::isnan(sections[sectionIndex].get_coordinates()[1]) ||sections[sectionIndex].get_coordinates()[1] < 0) {
                return false;
                // throw std::runtime_error("sections[sectionIndex].get_coordinates()[1] cannot be negative: " + std::to_string(sections[sectionIndex].get_coordinates()[1]));
            }


            if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                double wireDiameter = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_diameter().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    maximumNumberLayersFittingInSection = roundFloat(sections[sectionIndex].get_dimensions()[0] / wireDiameter, 9);
                    double averageLayerPerimeter = 2 * std::numbers::pi * (sections[sectionIndex].get_dimensions()[1] / 360) * (windingWindowRadialHeight - sections[sectionIndex].get_coordinates()[0]);
                    maximumNumberPhysicalTurnsPerLayer = floor(averageLayerPerimeter / wireDiameter);
                    layerRadialHeight = wireDiameter;
                    layerAngle = sections[sectionIndex].get_dimensions()[1];
                } else {
                    throw std::invalid_argument("Only overlapping layers allowed in toroids");
                }
            }
            else {
                double wireWidth = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_width().value());
                double wireHeight = resolve_dimensional_values(wirePerWinding[windingIndex].get_outer_height().value());
                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    maximumNumberLayersFittingInSection = roundFloat(sections[sectionIndex].get_dimensions()[0] / wireWidth, 9);
                    double averageLayerPerimeter = 2 * std::numbers::pi * (sections[sectionIndex].get_dimensions()[1] / 360) * (windingWindowRadialHeight - sections[sectionIndex].get_coordinates()[0]);
                    if (wirePerWinding[windingIndex].get_type() == WireType::FOIL) {
                        throw std::invalid_argument("Cannot have foil in toroids");
                    }
                    else {
                        maximumNumberPhysicalTurnsPerLayer = floor(averageLayerPerimeter / wireHeight);
                    }
                    layerRadialHeight = wireWidth;
                    layerAngle = sections[sectionIndex].get_dimensions()[1];
                } else {
                    throw std::invalid_argument("Only overlapping layers allowed in toroids");
                }
            }

            if (maximumNumberLayersFittingInSection == 0) {
                auto aux = get_number_layers_needed_and_number_physical_turns(sections[sectionIndex], wirePerWinding[windingIndex], physicalTurnsInSection, windingWindowRadialHeight);
                numberLayers = aux.first;
                layerPhysicalTurns = aux.second;
            }
            else if (maximumNumberPhysicalTurnsPerLayer == 0) {
                auto aux = get_number_layers_needed_and_number_physical_turns(sections[sectionIndex], wirePerWinding[windingIndex], physicalTurnsInSection, windingWindowRadialHeight);
                numberLayers = maximumNumberLayersFittingInSection;
                layerPhysicalTurns = aux.second;
            }
            else {
                auto aux = get_number_layers_needed_and_number_physical_turns(sections[sectionIndex], wirePerWinding[windingIndex], physicalTurnsInSection, windingWindowRadialHeight);
                minimumNumberLayerNeeded = aux.first;
                numberLayers = std::min(minimumNumberLayerNeeded, maximumNumberLayersFittingInSection);
                layerPhysicalTurns = aux.second;
            }

            // We cannot have more layers than physical turns
            if (sections[sectionIndex].get_number_layers()) {
                numberLayers = sections[sectionIndex].get_number_layers().value();
            }
            numberLayers = std::min(numberLayers, physicalTurnsInSection);

            if (minimumNumberLayerNeeded > numberLayers) {
                return false;
            }

            double currentLayerCenterRadialHeight;
            double currentLayerCenterAngle;
            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                currentLayerCenterRadialHeight = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerRadialHeight / 2, 9);
                currentLayerCenterAngle = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                throw std::invalid_argument("Only overlapping layers allowed in toroids");
            }

            WindingStyle windByConsecutiveTurns;
            if (sections[sectionIndex].get_winding_style()) {
                windByConsecutiveTurns = sections[sectionIndex].get_winding_style().value();
            }
            else {
                windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberLayers);
            }

            // For toroidal sections, layer capacity varies with radius (outer layers fit more
            // turns than inner ones). WIND_BY_CONSECUTIVE_PARALLELS places one turn per parallel
            // per layer and ignores per-layer capacity, so when any layer cannot fit one turn for
            // each parallel we must fall back to WIND_BY_CONSECUTIVE_TURNS, which honors the
            // per-layer turn count.
            int64_t minimumLayerPhysicalTurns = layerPhysicalTurns.empty() ? 0 : *std::min_element(layerPhysicalTurns.begin(), layerPhysicalTurns.end());
            if (windByConsecutiveTurns == WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS &&
                (maximumNumberPhysicalTurnsPerLayer < get_number_parallels(windingIndex) ||
                 (minimumLayerPhysicalTurns > 0 && uint64_t(minimumLayerPhysicalTurns) < get_number_parallels(windingIndex)))) {
                windByConsecutiveTurns = WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
            }

            for (size_t layerIndex = 0; layerIndex < numberLayers; ++layerIndex) {
                Layer layer;

                // TODO: probably will have to add a factor to have more or less turns per layer than the section average

                auto parallelsProportions = get_parallels_proportions(layerIndex,
                                                                       numberLayers,
                                                                       get_number_turns(windingIndex),
                                                                       get_number_parallels(windingIndex),
                                                                       remainingParallelsProportionInSection,
                                                                       windByConsecutiveTurns,
                                                                       totalParallelsProportionInSection,
                                                                       1,
                                                                       layerPhysicalTurns[layerIndex]);

                std::vector<double> layerParallelsProportion = parallelsProportions.second;

                size_t numberParallelsProportionsToZero = 0;
                for (auto parallelProportion : layerParallelsProportion) {
                    if (parallelProportion == 0) {
                        numberParallelsProportionsToZero++;
                    }
                }

                if (numberParallelsProportionsToZero == layerParallelsProportion.size()) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Parallel proportion in layer cannot be all be 0");
                }

                uint64_t physicalTurnsThisLayer = parallelsProportions.first;
                auto turnsAlignment = get_turns_alignment(sections[sectionIndex].get_name());

                partialWinding.set_parallels_proportion(layerParallelsProportion);
                layer.set_partial_windings(std::vector<PartialWinding>{partialWinding});
                layer.set_section(sections[sectionIndex].get_name());
                layer.set_type(ElectricalType::CONDUCTION);
                layer.set_name(sections[sectionIndex].get_name() +  " layer " + std::to_string(layerIndex));
                layer.set_orientation(sections[sectionIndex].get_layers_orientation());
                layer.set_turns_alignment(turnsAlignment);
                layer.set_dimensions(std::vector<double>{layerRadialHeight, layerAngle});
                layer.set_coordinates(std::vector<double>{currentLayerCenterRadialHeight, currentLayerCenterAngle, 0});
                layer.set_coordinate_system(CoordinateSystem::POLAR);

                double layerPerimeter = 2 * std::numbers::pi * (layerAngle / 360) * (windingWindowRadialHeight - layerRadialHeight / 2);
                layer.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisLayer) / (layerPerimeter * layerRadialHeight));
                layer.set_winding_style(windByConsecutiveTurns);
                layers.push_back(layer);

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingParallelsProportionInSection[parallelIndex] -= layerParallelsProportion[parallelIndex];
                }


                if (layerIndex == numberLayers - 1) {
                    break;
                }

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    currentLayerCenterRadialHeight = roundFloat(currentLayerCenterRadialHeight + layerRadialHeight / 2, 9);
                }
                else {
                    throw std::invalid_argument("Only overlapping layers allowed in toroids");
                }


                if (_insulationInterLayers.contains(windingIndex)) {
     
                    auto insulationLayer = _insulationInterLayers[windingIndex];

                    if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                        currentLayerCenterRadialHeight = roundFloat(currentLayerCenterRadialHeight + insulationLayer.get_dimensions()[0] / 2, 9);
                    }
                    else {
                        throw std::invalid_argument("Only overlapping layers allowed in toroids");
                    }

                    insulationLayer.set_section(sections[sectionIndex].get_name());
                    insulationLayer.set_coordinate_system(CoordinateSystem::POLAR);
                    insulationLayer.set_name(sections[sectionIndex].get_name() +  " insulation layer " + std::to_string(layerIndex));
                    insulationLayer.set_dimensions({insulationLayer.get_dimensions()[0], layerAngle});
                    insulationLayer.set_coordinates(std::vector<double>{currentLayerCenterRadialHeight, currentLayerCenterAngle, 0});
                    layers.push_back(insulationLayer);

                    if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                        currentLayerCenterRadialHeight = roundFloat(currentLayerCenterRadialHeight + insulationLayer.get_dimensions()[0] / 2, 9);
                    }
                    else {
                        throw std::invalid_argument("Only overlapping layers allowed in toroids");
                    }

                }

                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    currentLayerCenterRadialHeight = roundFloat(currentLayerCenterRadialHeight + layerRadialHeight / 2, 9);
                }
                else {
                    throw std::invalid_argument("Only overlapping layers allowed in toroids");
                }
            }
        }
        else {
            if (sectionIndex == 0) {
                throw NotImplementedException("Inner insulation layers not implemented");
            }

            auto partialWinding = sections[sectionIndex - 1].get_partial_windings()[0];
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());

            Section nextSection;
            if (sectionIndex + 1 != sections.size()) {
                if (sections[sectionIndex - 1].get_type() != ElectricalType::CONDUCTION || sections[sectionIndex + 1].get_type() != ElectricalType::CONDUCTION) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Previous and next sections must be conductive");
                }
                nextSection = sections[sectionIndex + 1];
            }
            else {
                nextSection = sections[0];
            }
            // auto nextSection = sections[sectionIndex + 1];
            auto nextPartialWinding = nextSection.get_partial_windings()[0];
            auto nextWindingIndex = get_winding_index_by_name(nextPartialWinding.get_winding());

            // If the angle of the section is 0 it means that the margin is enoguh and we don't need to add insulation layers
            if (sections[sectionIndex].get_dimensions()[1] == 0) {
                continue;
            }

            auto windingsMapKey = std::pair<size_t, size_t>{windingIndex, nextWindingIndex};
            if (!_insulationInterSectionsLayers.contains(windingsMapKey)) {
                log(_insulationInterSectionsLayersLog[windingsMapKey]);
                continue;
            }

            auto insulationLayers = _insulationInterSectionsLayers[windingsMapKey];
            if (insulationLayers.size() == 0) {
                throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "There must be at least one insulation layer between layers");
            }

            double layerRadialHeight = insulationLayers[0].get_dimensions()[0];

            double currentLayerCenterRadialHeight;
            double currentLayerCenterAngle;

            if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                currentLayerCenterRadialHeight = roundFloat(sections[sectionIndex].get_coordinates()[0] - sections[sectionIndex].get_dimensions()[0] / 2 + layerRadialHeight / 2, 9);
                currentLayerCenterAngle = roundFloat(sections[sectionIndex].get_coordinates()[1], 9);
            } else {
                throw std::invalid_argument("Only overlapping layers allowed in toroids");
            }

            for (size_t layerIndex = 0; layerIndex < insulationLayers.size(); ++layerIndex) {
                auto insulationLayer = insulationLayers[layerIndex];
                insulationLayer.set_section(sections[sectionIndex].get_name());
                insulationLayer.set_coordinate_system(CoordinateSystem::POLAR);
                insulationLayer.set_name(sections[sectionIndex].get_name() +  " insulation layer " + std::to_string(layerIndex));
                insulationLayer.set_coordinates(std::vector<double>{currentLayerCenterRadialHeight, currentLayerCenterAngle, 0});
                layers.push_back(insulationLayer);


                if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                    currentLayerCenterRadialHeight = roundFloat(currentLayerCenterRadialHeight + layerRadialHeight, 9);
                }
                else {
                    throw std::invalid_argument("Only overlapping layers allowed in toroids");
                }
            }
        }
    }
    set_layers_description(layers);
    return true;
}

bool Coil::wind_by_planar_layers() {
    set_layers_description(std::nullopt);
    std::vector<Layer> layers;
    if (!get_sections_description()) {
        return false;
    }

    auto sections = get_sections_description().value();

    for (auto section : sections) {
        Layer layer;
        layer.set_partial_windings(section.get_partial_windings());
        layer.set_section(section.get_name());
        layer.set_type(section.get_type());
        layer.set_orientation(section.get_layers_orientation());
        layer.set_dimensions(section.get_dimensions());
        layer.set_coordinates(section.get_coordinates());
        layer.set_coordinate_system(section.get_coordinate_system());
        layer.set_winding_style(section.get_winding_style());
        layer.set_filling_factor(section.get_filling_factor());
        layer.set_name(std::regex_replace(std::string(section.get_name()), std::regex("section"), "layer"));
        if (section.get_type() == ElectricalType::CONDUCTION) {
            layer.set_turns_alignment(CoilAlignment::SPREAD);
        }
        else {
            layer.set_insulation_material(defaults.defaultPcbInsulationMaterial);
        }

        layers.push_back(layer);

    }
    set_layers_description(layers);
    return true;
}

bool Coil::wind_by_turns() {
    set_turns_description(std::nullopt);
    if (!get_layers_description()) {
        return false;
    }
    auto bobbin = resolve_bobbin();

    auto functionalDescription = get_functional_description();
    auto sectionsDescription = get_sections_description().value();
    auto layersDescription = get_layers_description().value();
    auto needsVirtualization = needs_virtualization();

    if (needsVirtualization) {
        create_virtualization_map();
        auto virtualFunctionalDescription = virtualize_functional_description();
        auto virtualSectionsDescription = virtualize_sections_description();
        auto virtualLayersDescription = virtualize_layers_description();
        set_functional_description(virtualFunctionalDescription);
        set_sections_description(virtualSectionsDescription);
        set_layers_description(virtualLayersDescription);
        _windingIndexByName.clear();
        _turnIndexByName.clear();
    }

    bool result;
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        result = wind_by_rectangular_turns();
    }
    else {
        result = wind_by_round_turns();
    }

    if (needsVirtualization) {
        _windingIndexByName.clear();
        _turnIndexByName.clear();
        set_functional_description(functionalDescription);
        set_sections_description(sectionsDescription);
        set_layers_description(layersDescription);
        // Only devirtualize when the underlying wind actually populated
        // a turns description. wind_by_rectangular_turns / wind_by_round_turns
        // can return false (geometry doesn't fit) and leave the turns
        // description std::nullopt; devirtualize_turns_description's
        // unconditional .value() then died with bad_optional_access for
        // push_pull / weinberg with multi-winding virtualized cores.
        if (result) {
            devirtualize_turns_description();
        }
    }
    return result;
}

bool Coil::wind_by_rectangular_turns() {
    set_turns_description(std::nullopt);
    if (!get_layers_description()) {
        return false;
    }
    auto wirePerWinding = get_wires();
    std::vector<std::vector<int64_t>> currentTurnIndex;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        currentTurnIndex.push_back(std::vector<int64_t>(get_number_parallels(windingIndex), 0));
    }
    // Wound-column frame per section (multi-column winding support): the main window
    // resolves to the historical bobbin column scalars; sections placed in other
    // winding windows resolve their window's column edge.
    std::map<std::string, WoundColumnFrame> woundColumnFramePerSection;
    auto getFrameForSection = [&](const std::string& sectionName) -> const WoundColumnFrame& {
        auto frameIterator = woundColumnFramePerSection.find(sectionName);
        if (frameIterator == woundColumnFramePerSection.end()) {
            frameIterator = woundColumnFramePerSection.emplace(sectionName, get_wound_column_frame_for_section(sectionName)).first;
        }
        return frameIterator->second;
    };

    auto layers = get_layers_description().value();

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        if (wirePerWinding[windingIndex].get_type() == WireType::PLANAR) {
            auto conductionLayers = get_layers_by_type(ElectricalType::CONDUCTION);
            if (conductionLayers.size() > settings.get_coil_maximum_layers_planar()) {
                return false;
            }
        }

    }

    // Per-winding ordinal of each conduction layer in wound order. U winding alternates direction
    // layer by layer across the WHOLE winding (continuously, even over interleaved section breaks),
    // so the connection always leaves on the side the previous layer finished.
    std::map<std::string, int64_t> windingLayerOrderCount;
    std::vector<Turn> turns;
    for (auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            double currentTurnCenterWidth = 0;
            double currentTurnCenterHeight = 0;
            double currentTurnWidthIncrement = 0;
            double currentTurnHeightIncrement = 0;
            double totalLayerHeight;
            double totalLayerWidth;
            if (layer.get_partial_windings().size() > 1) {
                throw NotImplementedException("More than one winding per layer not supported yet");
            }
            auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();
            double wireHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
            auto physicalTurnsInLayer = get_number_turns(layer);
            auto alignment = layer.get_turns_alignment().value();

            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                totalLayerWidth = layer.get_dimensions()[0];
                totalLayerHeight = roundFloat(physicalTurnsInLayer * wireHeight, 9);

                currentTurnWidthIncrement = 0;
                currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0], 9);
                switch (alignment) {
                    case CoilAlignment::CENTERED:
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + totalLayerHeight / 2 - wireHeight / 2, 9);
                        currentTurnHeightIncrement = wireHeight;
                        break;

                    case CoilAlignment::INNER_OR_TOP:
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - wireHeight / 2, 9);
                        currentTurnHeightIncrement = wireHeight;
                        break;

                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 + totalLayerHeight - wireHeight / 2, 9);
                        currentTurnHeightIncrement = wireHeight;
                        break;

                    case CoilAlignment::SPREAD:
                        currentTurnHeightIncrement = roundFloat(layer.get_dimensions()[1] / physicalTurnsInLayer, 9);
                        currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - currentTurnHeightIncrement / 2, 9);
                        break;
                }

            } 
            else {
                // Place the turns at the geometry the caller committed to (the preset layer's
                // coordinates/dimensions). Whether that geometry actually fits the winding window is a
                // SEPARATE decision, owned by are_sections_and_layers_fitting() and enforced upstream in
                // wind()/rewind() (which only call wind_by_turns once fit — or windEvenIfNotFit — holds).
                // The old veto here ("contiguous turns wider than the window -> return false") re-made
                // that fit decision inside the turn placer, but ONLY for the contiguous branch: the
                // OVERLAPPING branch above and the round winder both render preset over-full layers
                // without objection. That asymmetry meant a direct wind_by_turns() over a preset,
                // over-full contiguous-rectangular coil (web bug 359 / ABT #160) uniquely produced no
                // turns at all, while every other orientation rendered them. Drop the veto so the
                // rectangular contiguous branch behaves like its siblings; degenerate placements (a turn
                // pushed to a non-positive radius) are still rejected by the turn-length guards below.
                totalLayerWidth = roundFloat(physicalTurnsInLayer * wireWidth, 9);
                totalLayerHeight = layer.get_dimensions()[1];
                currentTurnHeightIncrement = 0;
                currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1], 9);
                switch (alignment) {
                    case CoilAlignment::CENTERED:
                        currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - totalLayerWidth / 2 + wireWidth / 2, 9);
                        currentTurnWidthIncrement = wireWidth;
                        break;

                    case CoilAlignment::INNER_OR_TOP:
                        currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 + wireWidth / 2, 9);
                        currentTurnWidthIncrement = wireWidth;
                        break;

                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 + (layer.get_dimensions()[0] - totalLayerWidth) + wireWidth / 2, 9);
                        currentTurnWidthIncrement = wireWidth;
                        break;

                    case CoilAlignment::SPREAD:
                        currentTurnWidthIncrement = roundFloat(layer.get_dimensions()[0] / physicalTurnsInLayer, 9);
                        currentTurnCenterWidth = roundFloat(layer.get_coordinates()[0] - layer.get_dimensions()[0] / 2 + currentTurnWidthIncrement / 2, 9);
                        break;
                }
            }

            // U/Z winding order. In U order, every other conduction layer of the winding (counted in
            // wound order across all its sections) is wound in the opposite direction, so consecutive
            // electrical turns stay physically adjacent across the layer boundary (back-and-forth) and
            // the connection leaves on the side the previous layer finished. In Z order (default)
            // every layer is wound the same direction, implying a return wire (dragback) between
            // layers. The reversal is generic: start from the last forward position and negate the
            // increment, which holds for any turns alignment.
            std::string windingNameForOrder = layer.get_partial_windings()[0].get_winding();
            int64_t windingLayerOrdinal = windingLayerOrderCount[windingNameForOrder]++;
            // Convention: every winding STARTS FROM THE BOTTOM. The alignment above lays the first turn
            // at the top, so reverse it to start at the bottom and wind up — for every layer EXCEPT the
            // odd ordinals of a U winding, whose boustrophedon turnaround means they start at the top
            // (so consecutive layers stay adjacent and the layer-to-layer connection stays on top).
            if (!(get_winding_order(layer.get_section().value()) == WindingOrder::U && (windingLayerOrdinal % 2 == 1))) {
                currentTurnCenterWidth = roundFloat(currentTurnCenterWidth + (int64_t(physicalTurnsInLayer) - 1) * currentTurnWidthIncrement, 9);
                currentTurnCenterHeight = roundFloat(currentTurnCenterHeight - (int64_t(physicalTurnsInLayer) - 1) * currentTurnHeightIncrement, 9);
                currentTurnWidthIncrement = -currentTurnWidthIncrement;
                currentTurnHeightIncrement = -currentTurnHeightIncrement;
            }

            if (!layer.get_winding_style()) {
                layer.set_winding_style(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            }


            if (layer.get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    int64_t numberTurns = round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
                    for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                        Turn turn;
                        turn.set_coordinates(std::vector<double>{currentTurnCenterWidth, currentTurnCenterHeight});
                        turn.set_layer(layer.get_name());
                        {
                            auto turnLength = get_turn_length_in_frame(getFrameForSection(layer.get_section().value()), currentTurnCenterWidth);
                            if (!turnLength) {
                                return false;
                            }
                            turn.set_length(turnLength.value());
                        }
                        turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                        turn.set_orientation(TurnOrientation::CLOCKWISE);
                        turn.set_parallel(parallelIndex);
                        turn.set_section(layer.get_section().value());
                        turn.set_winding(partialWinding.get_winding());
                        turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                        turn.set_rotation(0);
                        if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                            turn.set_cross_sectional_shape(TurnCrossSectionalShape::ROUND);
                        }
                        else {
                            turn.set_cross_sectional_shape(TurnCrossSectionalShape::RECTANGULAR);
                        }
                        turn.set_coordinate_system(CoordinateSystem::CARTESIAN);
                        turns.push_back(turn);
                        currentTurnCenterWidth += currentTurnWidthIncrement;
                        currentTurnCenterHeight -= currentTurnHeightIncrement;
                        currentTurnIndex[windingIndex][parallelIndex]++; 
                    }
                }
            }
            else {
                int64_t firstParallelIndex = 0;
                while (roundFloat(partialWinding.get_parallels_proportion()[firstParallelIndex], 10) == 0) {
                    firstParallelIndex++;
                }
                int64_t numberTurns = round(partialWinding.get_parallels_proportion()[firstParallelIndex] * get_number_turns(windingIndex));
                for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                        if (roundFloat(partialWinding.get_parallels_proportion()[parallelIndex], 10) > 0) {
                            Turn turn;
                            turn.set_coordinates(std::vector<double>{currentTurnCenterWidth, currentTurnCenterHeight});
                            turn.set_layer(layer.get_name());
                            {
                                auto turnLength = get_turn_length_in_frame(getFrameForSection(layer.get_section().value()), currentTurnCenterWidth);
                                if (!turnLength) {
                                    return false;
                                }
                                turn.set_length(turnLength.value());
                            }
                            turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                            turn.set_orientation(TurnOrientation::CLOCKWISE);
                            turn.set_parallel(parallelIndex);
                            turn.set_section(layer.get_section().value());
                            turn.set_winding(partialWinding.get_winding());
                            turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                            turn.set_rotation(0);
                            if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                                turn.set_cross_sectional_shape(TurnCrossSectionalShape::ROUND);
                            }
                            else {
                                turn.set_cross_sectional_shape(TurnCrossSectionalShape::RECTANGULAR);
                            }
                            turn.set_coordinate_system(CoordinateSystem::CARTESIAN);

                            turns.push_back(turn);
                            currentTurnCenterWidth += currentTurnWidthIncrement;
                            currentTurnCenterHeight -= currentTurnHeightIncrement;
                            currentTurnIndex[windingIndex][parallelIndex]++; 
                        }
                    }
                }
            }
        }
    }

    set_turns_description(turns);
    return true;
}

bool Coil::can_build_centered_single_turn_toroidal() {
    // Bobbin must be set up with a round winding window (toroid).
    if (!std::holds_alternative<Bobbin>(get_bobbin())) {
        return false;
    }
    Bobbin bobbin = std::get<Bobbin>(get_bobbin());
    if (!bobbin.get_processed_description()) {
        return false;
    }
    auto pd = bobbin.get_processed_description().value();
    auto windingWindows = pd.get_winding_windows();
    if (windingWindows.empty() || !windingWindows[0].get_radial_height()) {
        return false;
    }
    auto shape = bobbin.get_winding_window_shape();
    if (shape != WindingWindowShape::ROUND) {
        return false;
    }
    // Exactly one functional winding with 1 turn × 1 parallel.
    const auto& fd = get_functional_description();
    if (fd.size() != 1) return false;
    if (get_number_turns(0) != 1 || get_number_parallels(0) != 1) return false;

    // Wire OD > inner radius (won't fit against the inner wall) and ≤ B
    // (still fits inside the hole). When _strict is false we tolerate
    // oversize-beyond-B wires for visualisation purposes.
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
    auto wires = get_wires();
    double wireOuterDiameter = wires[0].get_maximum_outer_width();
    if (wireOuterDiameter <= windingWindowRadialHeight) {
        return false;
    }
    if (_strict && wireOuterDiameter > 2 * windingWindowRadialHeight) {
        return false;
    }
    return true;
}

bool Coil::build_centered_single_turn_toroidal() {
    auto bobbin = resolve_bobbin();
    auto bobbinPd = bobbin.get_processed_description().value();
    auto bobbinColumnShape = bobbinPd.get_column_shape();
    auto bobbinColumnDepth = bobbinPd.get_column_depth();
    if (!bobbinPd.get_column_width()) {
        throw CoilNotProcessedException("Toroids must have their bobbin column set");
    }
    double bobbinColumnWidth = bobbinPd.get_column_width().value();
    auto windingWindows = bobbinPd.get_winding_windows();
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    auto wires = get_wires();
    auto& wire = wires[0];
    double wireOuterDiameter = wire.get_maximum_outer_width();
    double wireWidth = wireOuterDiameter;
    double wireHeight = wire.get_maximum_outer_height();
    auto windingName = get_name(0);

    PartialWinding partialWinding;
    partialWinding.set_winding(windingName);
    partialWinding.set_parallels_proportion({1.0});

    Section section;
    section.set_name(windingName + " section 0");
    section.set_partial_windings({partialWinding});
    section.set_type(ElectricalType::CONDUCTION);
    section.set_layers_orientation(WindingOrientation::OVERLAPPING);
    section.set_coordinate_system(CoordinateSystem::POLAR);
    // Section covers the whole inner-hole disk: full radial extent and 360°.
    section.set_dimensions(std::vector<double>{windingWindowRadialHeight, 360.0});
    // Section centre in polar: midway between bobbin inner surface and toroid axis.
    section.set_coordinates(std::vector<double>{windingWindowRadialHeight / 2.0, 0.0, 0.0});
    section.set_margin(std::vector<double>{0.0, 0.0});
    set_sections_description(std::vector<Section>{section});

    Layer layer;
    layer.set_name(windingName + " section 0 layer 0");
    layer.set_partial_windings({partialWinding});
    layer.set_type(ElectricalType::CONDUCTION);
    layer.set_section(section.get_name());
    layer.set_orientation(WindingOrientation::OVERLAPPING);
    layer.set_coordinate_system(CoordinateSystem::POLAR);
    layer.set_dimensions(std::vector<double>{wireWidth, 360.0});
    // Layer at polar radialHeight = B/2 → cartesian (0,0) via polar_to_cartesian.
    layer.set_coordinates(std::vector<double>{windingWindowRadialHeight, 0.0});
    layer.set_winding_style(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
    layer.set_turns_alignment(CoilAlignment::CENTERED);
    set_layers_description(std::vector<Layer>{layer});

    Turn turn;
    turn.set_coordinates(std::vector<double>{windingWindowRadialHeight, 0.0});
    turn.set_layer(layer.get_name());
    if (bobbinColumnShape == ColumnShape::ROUND) {
        turn.set_length(2 * std::numbers::pi * (windingWindowRadialHeight + bobbinColumnWidth));
    }
    else if (bobbinColumnShape == ColumnShape::OBLONG) {
        turn.set_length(2 * std::numbers::pi * (windingWindowRadialHeight + bobbinColumnWidth) + 4 * (bobbinColumnDepth - bobbinColumnWidth));
    }
    else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
        turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * windingWindowRadialHeight);
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "only round or rectangular columns supported for bobbins");
    }
    if (turn.get_length() < 0) {
        return false;
    }
    turn.set_name(windingName + " parallel 0 turn 0");
    turn.set_orientation(TurnOrientation::CLOCKWISE);
    turn.set_parallel(0);
    turn.set_section(section.get_name());
    turn.set_winding(windingName);
    turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
    turn.set_rotation(0.0);
    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        turn.set_cross_sectional_shape(TurnCrossSectionalShape::ROUND);
    }
    else {
        turn.set_cross_sectional_shape(TurnCrossSectionalShape::RECTANGULAR);
    }
    turn.set_coordinate_system(CoordinateSystem::POLAR);

    // Outer XY-plane crossing, same polar-mirror convention as wind_toroidal_additional_turns
    // ({-2*columnWidth - radialDepth, same angle}): the centered turn crosses the hole plane at
    // the toroid centre AND at the mirrored point outside the ring. Every other toroidal wind
    // path emits this via wind_toroidal_additional_turns; without it downstream consumers
    // (Painter, 3D builders) cannot know where the wire wraps the ring.
    if (settings.get_coil_include_additional_coordinates()) {
        turn.set_additional_coordinates(std::vector<std::vector<double>>{
            {-2 * bobbinColumnWidth - windingWindowRadialHeight, 0.0}});
    }

    set_turns_description(std::vector<Turn>{turn});
    convert_turns_to_cartesian_coordinates();
    return true;
}

bool Coil::wind_by_round_turns() {
    set_turns_description(std::nullopt);
    if (!get_layers_description()) {
        return false;
    }
    auto wirePerWinding = get_wires();
    std::vector<std::vector<int64_t>> currentTurnIndex;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        currentTurnIndex.push_back(std::vector<int64_t>(get_number_parallels(windingIndex), 0));
    }
    auto bobbinColumnShape = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_shape();
    auto bobbinColumnDepth = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_depth();
    double bobbinColumnWidth;
    if (std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_width()) {
        bobbinColumnWidth = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_width().value();
    }
    else {
        throw CoilNotProcessedException("Toroids must have their bobbin column set");
    }

    auto layers = get_layers_description().value();

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        if (wirePerWinding[windingIndex].get_type() == WireType::RECTANGULAR) {
            auto layersInWinding = get_layers_by_winding_index(windingIndex);
            if (layersInWinding.size() > 1) {
                return false;
            }
        }
    }

    // Per-winding ordinal of each conduction layer in wound order, so U winding alternates direction
    // continuously across the whole winding (see the rectangular winder for the rationale).
    std::map<std::string, int64_t> windingLayerOrderCount;
    std::vector<Turn> turns;
    for (auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            double currentTurnCenterRadialHeight = 0;
            double currentTurnCenterAngle = 0;
            double currentTurnRadialHeightIncrement = 0;
            double currentTurnAngleIncrement = 0;
            double totalLayerAngle;
            if (layer.get_partial_windings().size() > 1) {
                throw NotImplementedException("More than one winding per layer not supported yet");
            }
            auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();
            double wireHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
            auto physicalTurnsInLayer = get_number_turns(layer);
            auto alignment = layer.get_turns_alignment().value();

            auto bobbin = resolve_bobbin();
            auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
            double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

            double wireRadius;
            if (wirePerWinding[windingIndex].get_type() == WireType::RECTANGULAR) {
                wireRadius = windingWindowRadialHeight - layer.get_coordinates()[0] - wireWidth / 2;
            }
            else {
                wireRadius = windingWindowRadialHeight - layer.get_coordinates()[0];
            }
            double wireAngle = wound_distance_to_angle(wireHeight, wireRadius);
            // Physically-impossible placement is never a valid "loose fit": a
            // non-positive winding radius means the layer sits at or past the
            // toroid centre, and wireAngle > 180 (or NaN) means a single turn's
            // chord exceeds the diameter at this radius. Either way the wire
            // cannot be laid here, so reject regardless of _strict — otherwise
            // the non-strict adviser path emits turns at negative radii with a
            // 360 deg fallback angle that pile onto adjacent windings and only
            // get caught downstream by the collision check.
            if (wireRadius <= 0 || wireAngle > 180 || std::isnan(wireAngle)) {
                // Turns won't fit
                return false;
            }

            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                totalLayerAngle = physicalTurnsInLayer * wireAngle;

                currentTurnRadialHeightIncrement = 0;
                currentTurnCenterRadialHeight = roundFloat(layer.get_coordinates()[0], 9);
                switch (alignment) {
                    case CoilAlignment::CENTERED:
                        currentTurnCenterAngle = roundFloat(layer.get_coordinates()[1] - totalLayerAngle / 2 + wireAngle / 2, 9);
                        currentTurnAngleIncrement = wireAngle;
                        break;

                    case CoilAlignment::INNER_OR_TOP:
                        currentTurnCenterAngle = roundFloat(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 + wireAngle / 2, 9);
                        currentTurnAngleIncrement = wireAngle;
                        break;

                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentTurnCenterAngle = roundFloat(layer.get_coordinates()[1] + layer.get_dimensions()[1] / 2 - totalLayerAngle + wireAngle / 2, 9);
                        currentTurnAngleIncrement = wireAngle;
                        break;

                    case CoilAlignment::SPREAD:
                        currentTurnAngleIncrement = roundFloat(layer.get_dimensions()[1] / physicalTurnsInLayer, 9);
                        currentTurnCenterAngle = roundFloat(layer.get_coordinates()[1] - layer.get_dimensions()[1] / 2 + currentTurnAngleIncrement / 2, 9);
                        break;
                }

                // Real winding geometry distributes the turns around the core rather than packing them
                // into a tight arc. Only when the sections span the WHOLE winding window (overlapping
                // winding orientation) is the layer's full angle available — for contiguous sections
                // (angular sectors) the layer's angle is still the full window at this point (delimit
                // compacts it to the sector afterwards), so spreading over it would overrun the sector
                // into the neighbour; those keep the centred packing, which already fills the sector.
                // Spread ONLY when the even spacing is at least one wire (a toroid's wire-angle grows at
                // smaller radii); when a ring is at capacity, even spacing would be tighter than the wire
                // and overrun the circle, so keep the centred packing the winder already fit.
                if (settings.get_coil_use_real_winding_geometry() && get_winding_orientation() == WindingOrientation::OVERLAPPING) {
                    // Spread over the angle MINUS one wire, leaving a one-wire seam where the winding
                    // starts/ends (physically real, and it keeps a full ring from closing onto its own
                    // first turn — the consecutive-parallels placement can round to one extra turn).
                    double spreadIncrement = roundFloat((layer.get_dimensions()[1] - wireAngle) / physicalTurnsInLayer, 9);
                    if (spreadIncrement >= wireAngle) {
                        currentTurnAngleIncrement = spreadIncrement;
                        currentTurnCenterAngle = roundFloat(layer.get_coordinates()[1] - physicalTurnsInLayer * spreadIncrement / 2 + spreadIncrement / 2, 9);
                    }
                }

            }
            else {
                throw std::invalid_argument("Only overlapping layers allowed in toroids");
            }

            // U/Z winding order (toroidal). Same generic reversal as the rectangular winder: for U
            // order, every other conduction layer in a section starts from the last forward
            // position and winds back. Radial height and angle both advance with += here, so a
            // single formula reverses both.
            std::string windingNameForOrder = layer.get_partial_windings()[0].get_winding();
            int64_t windingLayerOrdinal = windingLayerOrderCount[windingNameForOrder]++;
            if (get_winding_order(layer.get_section().value()) == WindingOrder::U && (windingLayerOrdinal % 2 == 1)) {
                currentTurnCenterRadialHeight = roundFloat(currentTurnCenterRadialHeight + (int64_t(physicalTurnsInLayer) - 1) * currentTurnRadialHeightIncrement, 9);
                currentTurnCenterAngle = roundFloat(currentTurnCenterAngle + (int64_t(physicalTurnsInLayer) - 1) * currentTurnAngleIncrement, 9);
                currentTurnRadialHeightIncrement = -currentTurnRadialHeightIncrement;
                currentTurnAngleIncrement = -currentTurnAngleIncrement;
            }

            if (!layer.get_winding_style()) {
                layer.set_winding_style(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            }

            if (layer.get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    int64_t numberTurns = round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
                    for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                        Turn turn;
                        turn.set_coordinates(std::vector<double>{currentTurnCenterRadialHeight, currentTurnCenterAngle});
                        turn.set_layer(layer.get_name());
                        if (bobbinColumnShape == ColumnShape::ROUND) {
                            turn.set_length(2 * std::numbers::pi * (currentTurnCenterRadialHeight + bobbinColumnWidth));
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::OBLONG) {
                            turn.set_length(2 * std::numbers::pi * (currentTurnCenterRadialHeight + bobbinColumnWidth) + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                            double currentTurnCornerRadius = turn.get_coordinates()[0];
                            turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else {
                            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "only round or rectangular columns supported for bobbins");
                        }
                        turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                        turn.set_orientation(TurnOrientation::CLOCKWISE);
                        turn.set_parallel(parallelIndex);
                        turn.set_section(layer.get_section().value());
                        turn.set_winding(partialWinding.get_winding());
                        turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                        turn.set_rotation(currentTurnCenterAngle);
                        if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                            turn.set_cross_sectional_shape(TurnCrossSectionalShape::ROUND);
                        }
                        else {
                            turn.set_cross_sectional_shape(TurnCrossSectionalShape::RECTANGULAR);
                        }
                        turn.set_coordinate_system(CoordinateSystem::POLAR);

                        turns.push_back(turn);
                        currentTurnCenterRadialHeight += currentTurnRadialHeightIncrement;
                        currentTurnCenterAngle += currentTurnAngleIncrement;
                        currentTurnIndex[windingIndex][parallelIndex]++; 
                    }
                }
            }
            else {
                int64_t firstParallelIndex = 0;
                while (roundFloat(partialWinding.get_parallels_proportion()[firstParallelIndex], 10) == 0) {
                    firstParallelIndex++;
                }
                int64_t numberTurns = round(partialWinding.get_parallels_proportion()[firstParallelIndex] * get_number_turns(windingIndex));
                for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                        if (roundFloat(partialWinding.get_parallels_proportion()[parallelIndex], 10) > 0) {
                            Turn turn;
                            turn.set_coordinates(std::vector<double>{currentTurnCenterRadialHeight, currentTurnCenterAngle});
                            turn.set_layer(layer.get_name());
                            if (bobbinColumnShape == ColumnShape::ROUND) {
                                turn.set_length(2 * std::numbers::pi * (currentTurnCenterRadialHeight + bobbinColumnWidth));
                                    if (turn.get_length() < 0) {
                                        return false;
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::OBLONG) {
                                turn.set_length(2 * std::numbers::pi * (currentTurnCenterRadialHeight + bobbinColumnWidth) + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                                    if (turn.get_length() < 0) {
                                        return false;
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                double currentTurnCornerRadius = currentTurnCenterRadialHeight;
                                turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                                if (turn.get_length() < 0) {
                                    return false;
                                }
                            }
                            else {
                                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "only round or rectangular columns supported for bobbins");
                            }
                            turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                            turn.set_orientation(TurnOrientation::CLOCKWISE);
                            turn.set_parallel(parallelIndex);
                            turn.set_section(layer.get_section().value());
                            turn.set_winding(partialWinding.get_winding());
                            turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                            turn.set_rotation(currentTurnCenterAngle);
                            if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                                turn.set_cross_sectional_shape(TurnCrossSectionalShape::ROUND);
                            }
                            else {
                                turn.set_cross_sectional_shape(TurnCrossSectionalShape::RECTANGULAR);
                            }
                            turn.set_coordinate_system(CoordinateSystem::POLAR);


                            turns.push_back(turn);
                            currentTurnCenterRadialHeight += currentTurnRadialHeightIncrement;
                            currentTurnCenterAngle += currentTurnAngleIncrement;
                            currentTurnIndex[windingIndex][parallelIndex]++; 
                        }
                    }
                }
            }
        }
    }

    set_turns_description(turns);

    convert_turns_to_cartesian_coordinates();
    return true;
}

bool Coil::wind_by_planar_turns(double borderToWireDistance, std::map<size_t, double> wireToWireDistance) {
    set_turns_description(std::nullopt);
    if (!get_layers_description()) {
        return false;
    }
    auto wirePerWinding = get_wires();

    std::vector<std::vector<int64_t>> currentTurnIndex;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        currentTurnIndex.push_back(std::vector<int64_t>(get_number_parallels(windingIndex), 0));
    }
    auto bobbin = resolve_bobbin();
    auto bobbinColumnShape = bobbin.get_processed_description().value().get_column_shape();
    auto bobbinColumnDepth = bobbin.get_processed_description().value().get_column_depth();
    double bobbinColumnWidth;
    if (bobbin.get_processed_description().value().get_column_width()) {
        bobbinColumnWidth = bobbin.get_processed_description().value().get_column_width().value();
    }
    else {
        auto bobbinWindingWindow = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows()[0];
        double bobbinWindingWindowWidth = bobbinWindingWindow.get_width().value();
        double bobbinWindingWindowCenterWidth = bobbinWindingWindow.get_coordinates().value()[0];
        bobbinColumnWidth = bobbinWindingWindowCenterWidth - bobbinWindingWindowWidth / 2;
    }

    auto layers = get_layers_description().value();

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        if (wirePerWinding[windingIndex].get_type() == WireType::PLANAR) {
            auto conductionLayers = get_layers_by_type(ElectricalType::CONDUCTION);
            if (conductionLayers.size() > settings.get_coil_maximum_layers_planar()) {
                return false;
            }
        }
    }

    std::vector<Turn> turns;
    for (auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            if (layer.get_partial_windings().size() > 1) {
                throw NotImplementedException("More than one winding per layer not supported yet");
            }
            auto partialWinding = layer.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();
            double wireHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
            double layerTurnsClearance;

            if (wireToWireDistance.count(windingIndex)) {
                layerTurnsClearance = wireToWireDistance[windingIndex];
            }
            else {
                layerTurnsClearance = defaults.minimumWireToWireDistance;
            }
            double currentTurnWidthIncrement = wireWidth + layerTurnsClearance;
            double currentTurnHeightIncrement = 0;

            // Center the planar turns within the layer instead of left-justifying them.
            // The layer is already inset from the core walls by coreToLayerDistance, so
            // splitting the remaining width evenly keeps borderToWireDistance as a per-edge
            // minimum while giving symmetric core-to-copper clearance on BOTH sides. Without
            // this the copper hugs the left layer edge and drifts toward (and, as
            // coreToLayerDistance grows, past) the opposite core wall.
            int64_t physicalTurnsThisLayer = 0;
            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                physicalTurnsThisLayer += round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
            }
            double turnsBlockWidth = 0;
            if (physicalTurnsThisLayer > 0) {
                turnsBlockWidth = physicalTurnsThisLayer * wireWidth + (physicalTurnsThisLayer - 1) * layerTurnsClearance;
            }
            double layerWidth = layer.get_dimensions()[0];
            double layerLeftEdge = layer.get_coordinates()[0] - layerWidth / 2;
            double turnsBlockMargin = (layerWidth - turnsBlockWidth) / 2;
            if (turnsBlockMargin < borderToWireDistance) {
                turnsBlockMargin = borderToWireDistance;
            }
            double currentTurnCenterWidth = roundFloat(layerLeftEdge + turnsBlockMargin + wireWidth / 2, 9);
            double currentTurnCenterHeight = roundFloat(layer.get_coordinates()[1], 9);

            if (!layer.get_winding_style()) {
                layer.set_winding_style(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            }

            if (layer.get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    int64_t numberTurns = round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
                    double totalWidthNeeded = borderToWireDistance * 2 + numberTurns * wireWidth + (numberTurns - 1) * layerTurnsClearance;
                    if(_strict && totalWidthNeeded > layer.get_dimensions()[0]) {
                        return false;
                    }

                    for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                        Turn turn;
                        turn.set_coordinates(std::vector<double>{currentTurnCenterWidth, currentTurnCenterHeight});
                        turn.set_layer(layer.get_name());
                        if (bobbinColumnShape == ColumnShape::ROUND) {
                            turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth);
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::OBLONG) {
                            turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                            double currentTurnCornerRadius = currentTurnCenterWidth - bobbinColumnWidth;
                            turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);

                            if (turn.get_length() < 0) {
                                return false;
                            }
                        }
                        else {
                            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "only round or rectangular columns supported for bobbins");
                        }
                        turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                        turn.set_orientation(TurnOrientation::CLOCKWISE);
                        turn.set_parallel(parallelIndex);
                        turn.set_section(layer.get_section().value());
                        turn.set_winding(partialWinding.get_winding());
                        turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                        turn.set_rotation(0);
                        turn.set_cross_sectional_shape(TurnCrossSectionalShape::RECTANGULAR);
                        turn.set_coordinate_system(CoordinateSystem::CARTESIAN);

                        turns.push_back(turn);
                        currentTurnCenterWidth += currentTurnWidthIncrement;
                        currentTurnCenterHeight -= currentTurnHeightIncrement;
                        currentTurnIndex[windingIndex][parallelIndex]++; 
                    }
                }
            }
            else {
                int64_t firstParallelIndex = 0;
                while (roundFloat(partialWinding.get_parallels_proportion()[firstParallelIndex], 10) == 0) {
                    firstParallelIndex++;
                }
                int64_t numberTurns = round(partialWinding.get_parallels_proportion()[firstParallelIndex] * get_number_turns(windingIndex));
                for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                    for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                        if (roundFloat(partialWinding.get_parallels_proportion()[parallelIndex], 10) > 0) {
                            Turn turn;
                            turn.set_coordinates(std::vector<double>{currentTurnCenterWidth, currentTurnCenterHeight});
                            turn.set_layer(layer.get_name());
                            if (bobbinColumnShape == ColumnShape::ROUND) {
                                turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth);
                                    if (turn.get_length() < 0) {
                                        return false;
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::OBLONG) {
                                turn.set_length(2 * std::numbers::pi * currentTurnCenterWidth + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                                    if (turn.get_length() < 0) {
                                        return false;
                                    }
                            }
                            else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                double currentTurnCornerRadius = currentTurnCenterWidth - bobbinColumnWidth;
                                turn.set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                                if (turn.get_length() < 0) {
                                    return false;
                                }
                            }
                            else {
                                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "only round or rectangular columns supported for bobbins");
                            }
                            turn.set_name(partialWinding.get_winding() + " parallel " + std::to_string(parallelIndex) + " turn " + std::to_string(currentTurnIndex[windingIndex][parallelIndex]));
                            turn.set_orientation(TurnOrientation::CLOCKWISE);
                            turn.set_parallel(parallelIndex);
                            turn.set_section(layer.get_section().value());
                            turn.set_winding(partialWinding.get_winding());
                            turn.set_dimensions(std::vector<double>{wireWidth, wireHeight});
                            turn.set_rotation(0);
                            turn.set_cross_sectional_shape(TurnCrossSectionalShape::RECTANGULAR);
                            turn.set_coordinate_system(CoordinateSystem::CARTESIAN);

                            turns.push_back(turn);
                            currentTurnCenterWidth += currentTurnWidthIncrement;
                            currentTurnCenterHeight -= currentTurnHeightIncrement;
                            currentTurnIndex[windingIndex][parallelIndex]++; 
                        }
                    }
                }
            }
        }
    }

    set_turns_description(turns);
    return true;
}

std::vector<std::pair<double, std::vector<double>>> Coil::get_collision_distances(std::vector<double> turnCoordinates, std::vector<std::vector<double>> placedTurnsCoordinates, double wireHeight) {
    std::vector<std::pair<double, std::vector<double>>> collisions;
    auto turnCartesianCoordinates = polar_to_cartesian(turnCoordinates);
    for (auto placedTurnCoordinates : placedTurnsCoordinates) {
        auto placedTurnCartesianCoordinates = polar_to_cartesian(placedTurnCoordinates);
        double distance = sqrt(pow(turnCartesianCoordinates[0] - placedTurnCartesianCoordinates[0], 2) + pow(turnCartesianCoordinates[1] - placedTurnCartesianCoordinates[1], 2));
        // Use a small tolerance to account for floating point precision
        if (distance < wireHeight - 1e-9) {
            double collisionDistance = wireHeight - distance;
            auto placedCoordinates = placedTurnCoordinates;
            collisions.push_back({collisionDistance, placedCoordinates});
        }

        if (collisions.size() == 2) {
            break;
        }
    }

    return collisions;
}

bool Coil::wind_toroidal_additional_turns() {
    if (!get_layers_description()) {
        return false;
    }
    if (!get_turns_description()) {
        return false;
    }
    auto wirePerWinding = get_wires();
    std::vector<std::vector<int64_t>> currentTurnIndex;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        currentTurnIndex.push_back(std::vector<int64_t>(get_number_parallels(windingIndex), 0));
    }
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    double bobbinColumnWidth;
    if (bobbin.get_processed_description().value().get_column_width()) {
        bobbinColumnWidth = bobbin.get_processed_description().value().get_column_width().value();
    }
    else {
        throw CoilNotProcessedException("Toroids must have their bobbin column set");
    }
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
    auto bobbinColumnShape = bobbin.get_processed_description().value().get_column_shape();
    auto bobbinColumnDepth = bobbin.get_processed_description().value().get_column_depth();

    auto sections = get_sections_description().value();
    auto layers = get_layers_description().value();
    auto turns = get_turns_description().value();
    for (auto& t : turns) {
        t.set_additional_coordinates(std::nullopt);
    }
    set_turns_description(turns);
    double currentBaseRadialHeight = -bobbinColumnWidth * 2;
    std::vector<std::pair<Layer, double>> maximumAdditionalRadialHeightPerInsulationLayerByIndex;
    auto windingOrientation = get_winding_orientation();

    for (auto section : sections) { 
        if (section.get_type() == ElectricalType::CONDUCTION) {
            std::vector<std::vector<double>> placedTurnsCoordinates;
            auto turnsInSection = get_turns_by_section(section.get_name());
            auto partialWinding = section.get_partial_windings()[0];  // TODO: Support multiwinding in layers
            auto winding = get_winding_by_name(partialWinding.get_winding());
            auto windingIndex = get_winding_index_by_name(partialWinding.get_winding());
            // double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();
            double wireHeight = wirePerWinding[windingIndex].get_maximum_outer_height();
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentBaseRadialHeight -= turnsInSection[0].get_dimensions().value()[0] / 2;
            }
            else {
                currentBaseRadialHeight = -bobbinColumnWidth * 2 - turnsInSection[0].get_dimensions().value()[0] / 2;
            }
            double currentSectionMaximumAdditionalRadialHeight = 0;
            double currentBaseRadialHeightForLayers = currentBaseRadialHeight + turnsInSection[0].get_dimensions().value()[0] / 2;

            auto layersThisSection = get_layers_by_section(section.get_name());
            // Check if there are at least 2 conduction layers with NO real insulation between them
            // Real insulation = thickness > MIN_INSULATION_THICKNESS
            const double MIN_INSULATION_THICKNESS = 1e-9;
            bool areLayersTaped = true;
            
            // Find all conduction layer indices
            std::vector<size_t> conductionLayerIndices;
            for (size_t i = 0; i < layersThisSection.size(); ++i) {
                if (layersThisSection[i].get_type() == ElectricalType::CONDUCTION) {
                    conductionLayerIndices.push_back(i);
                }
            }
            
            // If we have 2+ conduction layers, check if there's real insulation between any pair
            if (conductionLayerIndices.size() >= 2) {
                areLayersTaped = false;  // Assume no tape until we find real insulation
                
                for (size_t i = 0; i < conductionLayerIndices.size() - 1; ++i) {
                    size_t firstCondIdx = conductionLayerIndices[i];
                    size_t secondCondIdx = conductionLayerIndices[i + 1];
                    
                    // Check all layers between these two conduction layers
                    bool hasRealInsulation = false;
                    for (size_t j = firstCondIdx + 1; j < secondCondIdx; ++j) {
                        if (layersThisSection[j].get_type() == ElectricalType::INSULATION) {
                            double insulationThickness = layersThisSection[j].get_dimensions()[0];
                            if (insulationThickness > MIN_INSULATION_THICKNESS) {
                                hasRealInsulation = true;
                                break;
                            }
                        }
                    }
                    
                    if (hasRealInsulation) {
                        areLayersTaped = true;
                    }
                }
            }
            
            size_t conductionLayerCount = 0;
            for (auto layer : layersThisSection) {
                if (layer.get_type() == ElectricalType::CONDUCTION) {
                    auto turnsThisLayer = get_turns_by_layer(layer.get_name());
                    bool isFirstConductionLayer = (conductionLayerCount == 0);
                    conductionLayerCount++;
                    for (auto turn : turnsThisLayer) {
                        auto turnIndex = get_turn_index_by_name(turn.get_name());
                        std::vector<double> additionalCoordinates = {-bobbinColumnWidth * 2 - turn.get_coordinates()[0], turn.get_coordinates()[1]};

                        if (!areLayersTaped) {

                            if (!isFirstConductionLayer) {
                                {
                                    auto collisions = get_collision_distances(additionalCoordinates, placedTurnsCoordinates, wireHeight);
                                    if (collisions.size() > 0) {
                                        double currentRadius = windingWindowRadialHeight - additionalCoordinates[0];
                                        double currentWireAngle = ceilFloat(wound_distance_to_angle(wireHeight, currentRadius), 3);
                                        double sectionMinAngle = section.get_coordinates()[1] - section.get_dimensions()[1] / 2;
                                        double sectionMaxAngle = section.get_coordinates()[1] + section.get_dimensions()[1] / 2;
                                        double scanStep = currentWireAngle / 2;
                                        double defaultAngle = additionalCoordinates[1];
                                        bool foundSlot = false;
                                        for (double offset = scanStep; !foundSlot; offset += scanStep) {
                                            for (int sign : {-1, +1}) {
                                                double testAngle = defaultAngle + sign * offset;
                                                if (testAngle < sectionMinAngle + currentWireAngle / 2 ||
                                                    testAngle > sectionMaxAngle - currentWireAngle / 2) continue;
                                                std::vector<double> testCoords = {additionalCoordinates[0], testAngle};
                                                auto testCollisions = get_collision_distances(testCoords, placedTurnsCoordinates, wireHeight);
                                                if (testCollisions.size() == 0) {
                                                    additionalCoordinates[1] = testAngle;
                                                    foundSlot = true;
                                                    break;
                                                }
                                            }
                                            if (offset > section.get_dimensions()[1]) break;
                                        }
                                    }
                                }

                                std::vector<double> newCoordinates = {additionalCoordinates[0], additionalCoordinates[1]};
                                newCoordinates[0] = currentBaseRadialHeight;
                                auto collisions = get_collision_distances(newCoordinates, placedTurnsCoordinates, wireHeight);

                                if (collisions.size() > 0) {
                                    bool tryAngularMove = collisions.size() > 0;
                                    bool tryReversedAngularMove = collisions.size() > 0;
                                    bool previouslyAdditionAngularMovement = false;
                                    bool try0Degrees = true;
                                    bool tryMinus0Degrees = true;
                                    bool try30Degrees = true;
                                    bool tryMinus30Degrees = true;
                                    bool try45Degrees = true;
                                    bool tryMinus45Degrees = true;
                                    bool try60Degrees = true;
                                    bool tryMinus60Degrees = true;
                                    bool tryAvoidingCollisionDistance = true;
                                    double previousCollisionDistance = 0;
                                    std::vector<double> originalCollidedCoordinate;
                                    double restoredHeightAfter60Degrees = 0;

                                    double collisionDistance = collisions[0].first;
                                    auto collidedCoordinate = collisions[0].second;

                                    uint64_t timeout = 1000;
                                    while (newCoordinates[0] > additionalCoordinates[0]) {
                                        timeout--;
                                        if (timeout == 0) {
                                            throw CalculationException(ErrorCode::CALCULATION_TIMEOUT, "timeout in wind_toroidal_additional_turns");
                                        }
                                        if (tryAvoidingCollisionDistance && collisionDistance < 1e-6) {
                                            tryAvoidingCollisionDistance = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            double collisionAngle = ceilFloat(wound_distance_to_angle(collisionDistance, currentRadius), 3);
                                            if (collidedCoordinate[1] > newCoordinates[1]) {
                                                newCoordinates[1] -= collisionAngle;
                                            }
                                            else {
                                                newCoordinates[1] += collisionAngle;
                                            }
                                        }
                                        else if (tryAngularMove) {
                                            tryAngularMove = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            double increment = ceilFloat(wound_distance_to_angle(collisionDistance, currentRadius), 3);
                                            if (collidedCoordinate[1] > newCoordinates[1]) {
                                                previouslyAdditionAngularMovement = false;
                                                if (newCoordinates[1] - increment > (section.get_coordinates()[1] - section.get_dimensions()[1] / 2))
                                                    newCoordinates[1] -= increment;
                                            }
                                            else {
                                                previouslyAdditionAngularMovement = true;
                                                if (newCoordinates[1] + increment < (section.get_coordinates()[1] + section.get_dimensions()[1] / 2))
                                                    newCoordinates[1] += increment;
                                            }
                                        }
                                        else if (tryReversedAngularMove) {
                                            tryReversedAngularMove = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            double currentAngleCollision = ceilFloat(wound_distance_to_angle(previousCollisionDistance, currentRadius), 3);
                                            double currentWireAngle = ceilFloat(wound_distance_to_angle(wireHeight, currentRadius), 3);
                                            double currentAngleMovement = currentWireAngle + (currentWireAngle - currentAngleCollision);

                                            if (previouslyAdditionAngularMovement) {
                                                if (newCoordinates[1] - currentAngleMovement > (section.get_coordinates()[1] - section.get_dimensions()[1] / 2))
                                                    newCoordinates[1] -= currentAngleMovement;
                                            }
                                            else {
                                                if (newCoordinates[1] + currentAngleMovement < (section.get_coordinates()[1] + section.get_dimensions()[1] / 2))
                                                    newCoordinates[1] += currentAngleMovement;
                                            }
                                        }
                                        else if (try0Degrees) {
                                            try0Degrees = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            restoredHeightAfter60Degrees = newCoordinates[0];
                                            newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(0);
                                            newCoordinates[1] = originalCollidedCoordinate[1] + ceilFloat(wound_distance_to_angle(wireHeight * cos(0), currentRadius), 3);
                                        }
                                        else if (tryMinus0Degrees) {
                                            tryMinus0Degrees = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(0);
                                            newCoordinates[1] = originalCollidedCoordinate[1] - ceilFloat(wound_distance_to_angle(wireHeight * cos(0), currentRadius), 3);
                                        }
                                        else if (try30Degrees) {
                                            try30Degrees = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 6);
                                            newCoordinates[1] = originalCollidedCoordinate[1] + ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 6), currentRadius), 3);
                                        }
                                        else if (tryMinus30Degrees) {
                                            tryMinus30Degrees = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 6);
                                            newCoordinates[1] = originalCollidedCoordinate[1] - ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 6), currentRadius), 3);
                                        }
                                        else if (try45Degrees) {
                                            try45Degrees = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 4);
                                            newCoordinates[1] = originalCollidedCoordinate[1] + ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 4), currentRadius), 3);
                                        }
                                        else if (tryMinus45Degrees) {
                                            tryMinus45Degrees = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 4);
                                            newCoordinates[1] = originalCollidedCoordinate[1] - ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 4), currentRadius), 3);
                                        }
                                        else if (try60Degrees) {
                                            try60Degrees = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 3);
                                            newCoordinates[1] = originalCollidedCoordinate[1] + ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 3), currentRadius), 3);
                                        }
                                        else if (tryMinus60Degrees) {
                                            tryMinus60Degrees = false;
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            newCoordinates[0] = originalCollidedCoordinate[0] - wireHeight * sin(std::numbers::pi / 3);
                                            newCoordinates[1] = originalCollidedCoordinate[1] - ceilFloat(wound_distance_to_angle(wireHeight * cos(std::numbers::pi / 3), currentRadius), 3);
                                        }
                                        else {
                                            // Before falling back to a new radial layer, try to find an empty slot at the first layer
                                            // by scanning through all angular positions
                                            double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                            double currentWireAngle = ceilFloat(wound_distance_to_angle(wireHeight, currentRadius), 3);
                                            double sectionMinAngle = section.get_coordinates()[1] - section.get_dimensions()[1] / 2;
                                            double sectionMaxAngle = section.get_coordinates()[1] + section.get_dimensions()[1] / 2;
                                            
                                            bool foundSlot = false;
                                            // Only try slot scanning if we're still at the first additional layer
                                            if (std::abs(newCoordinates[0] - currentBaseRadialHeight) < turn.get_dimensions().value()[0] / 4) {
                                                // Scan from minimum angle to maximum angle looking for an empty slot
                                                // Use a smaller step for denser packing - half the wire angle for better slot finding
                                                double scanStep = currentWireAngle / 2;
                                                for (double testAngle = sectionMinAngle + currentWireAngle / 2; testAngle <= sectionMaxAngle - currentWireAngle / 2; testAngle += scanStep) {
                                                    std::vector<double> testCoords = {currentBaseRadialHeight, testAngle};
                                                    auto testCollisions = get_collision_distances(testCoords, placedTurnsCoordinates, wireHeight);
                                                    if (testCollisions.size() == 0) {
                                                        newCoordinates = testCoords;
                                                        foundSlot = true;
                                                        break;
                                                    }
                                                }
                                            }
                                            
                                            if (!foundSlot) {
                                                // Fall back to a new radial layer
                                                try0Degrees = true;
                                                tryMinus0Degrees = true;
                                                try30Degrees = true;
                                                tryMinus30Degrees = true;
                                                try45Degrees = true;
                                                tryMinus45Degrees = true;
                                                try60Degrees = true;
                                                tryMinus60Degrees = true;
                                                tryAngularMove = true;
                                                tryAvoidingCollisionDistance = true;
                                                previousCollisionDistance = 0;
                                                if (restoredHeightAfter60Degrees != 0) {
                                                    newCoordinates[0] = restoredHeightAfter60Degrees;
                                                    restoredHeightAfter60Degrees = 0;
                                                }
                                                newCoordinates[0] -= turn.get_dimensions().value()[0] / 2;
                                                newCoordinates[1] = additionalCoordinates[1];
                                            }
                                        }
                                        double currentRadius = windingWindowRadialHeight - currentBaseRadialHeight;
                                        double currentWireAngle = ceilFloat(wound_distance_to_angle(wireHeight, currentRadius), 3);

                                        // Normalize angles for comparison to handle wrap-around cases
                                        double sectionMinAngle = section.get_coordinates()[1] - section.get_dimensions()[1] / 2;
                                        double sectionMaxAngle = section.get_coordinates()[1] + section.get_dimensions()[1] / 2;
                                        
                                        // Normalize all angles to [0, 360) range for consistent comparison
                                        auto normalizeAngle = [](double angle) -> double {
                                            while (angle < 0) angle += 360;
                                            while (angle >= 360) angle -= 360;
                                            return angle;
                                        };
                                        
                                        double normNewAngle = normalizeAngle(newCoordinates[1]);
                                        double normSectionMin = normalizeAngle(sectionMinAngle);
                                        double normSectionMax = normalizeAngle(sectionMaxAngle);
                                        double normAdditionalAngle = normalizeAngle(additionalCoordinates[1]);
                                        
                                        // Check if angle is outside section bounds
                                        // For full 360 degree sections, skip the check
                                        bool isFullCircle = (section.get_dimensions()[1] >= 360 - 1e-6);
                                        
                                        if (!isFullCircle) {
                                            bool outsideBounds = false;
                                            if (normSectionMin < normSectionMax) {
                                                // Normal case: min < max
                                                outsideBounds = (normNewAngle < normSectionMin + currentWireAngle / 2) || 
                                                               (normNewAngle > normSectionMax - currentWireAngle / 2);
                                            } else {
                                                // Wrap-around case: min > max (section crosses 0/360 boundary)
                                                outsideBounds = (normNewAngle < normSectionMin + currentWireAngle / 2) && 
                                                               (normNewAngle > normSectionMax - currentWireAngle / 2);
                                            }
                                            
                                            if (outsideBounds) {
                                                newCoordinates[1] = normAdditionalAngle;
                                            }
                                        }

                                        collisions = get_collision_distances(newCoordinates, placedTurnsCoordinates, wireHeight);
                                        if (collisions.size() == 0) {
                                            break;
                                        }
                                        collidedCoordinate = collisions[0].second;
                                        if (previousCollisionDistance == 0) {
                                            originalCollidedCoordinate = collidedCoordinate;
                                        }
                                        previousCollisionDistance = collisionDistance;
                                        collisionDistance = collisions[0].first;
                                    }
                                }
                                additionalCoordinates = newCoordinates;
                            }
                        }
                        currentSectionMaximumAdditionalRadialHeight = std::min(currentSectionMaximumAdditionalRadialHeight, additionalCoordinates[0]);
                        turn.set_additional_coordinates(std::vector<std::vector<double>>{additionalCoordinates});

                        if (bobbinColumnShape == ColumnShape::ROUND) {
                            double b = (turn.get_coordinates()[0] - turn.get_additional_coordinates().value()[0][0]) / 2;
                            double a = turn.get_coordinates()[0];
                            // Ramanujan  approximation for ellipse perimeter
                            double perimeter = std::numbers::pi * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b)));
                            turns[turnIndex].set_length(perimeter);
                            if (turns[turnIndex].get_length() < 0) {
                                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Something wrong happened in turn length 1: " + std::to_string(turns[turnIndex].get_length()) + " turns[turnIndex].get_coordinates()[0]: " + std::to_string(turns[turnIndex].get_coordinates()[0]));
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::OBLONG) {
                            double b = (turn.get_coordinates()[0] - turn.get_additional_coordinates().value()[0][0]) / 2;
                            double a = turn.get_coordinates()[0];
                            // Ramanujan  approximation for ellipse perimeter
                            double perimeter = std::numbers::pi * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b))) + 4 * (bobbinColumnDepth - bobbinColumnWidth);
                            turns[turnIndex].set_length(perimeter);
                            if (turns[turnIndex].get_length() < 0) {
                                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Something wrong happened in turn length 1: " + std::to_string(turns[turnIndex].get_length()) + " turns[turnIndex].get_coordinates()[0]: " + std::to_string(turns[turnIndex].get_coordinates()[0]));
                            }
                        }
                        else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                            double currentInternalTurnCornerRadius = turns[turnIndex].get_coordinates()[0];
                            double currentExternalTurnCornerRadius = -turn.get_additional_coordinates().value()[0][0] - 2 * bobbinColumnWidth;
                            double maximumVerticalDistance = currentInternalTurnCornerRadius * 2 + 2 * bobbinColumnDepth;
                            double externalVerticalStraightDistance = maximumVerticalDistance - 2 * currentExternalTurnCornerRadius;
                            turns[turnIndex].set_length(2 * bobbinColumnDepth + 4 * bobbinColumnWidth + externalVerticalStraightDistance + std::numbers::pi * currentInternalTurnCornerRadius + std::numbers::pi * currentExternalTurnCornerRadius);

                            if (turns[turnIndex].get_length() < 0) {
                                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Something wrong happened in turn length 1: " + std::to_string(turns[turnIndex].get_length()) + " bobbinColumnDepth: " + std::to_string(bobbinColumnDepth)  + " bobbinColumnWidth: " + std::to_string(bobbinColumnWidth)  + " currentExternalTurnCornerRadius: " + std::to_string(currentExternalTurnCornerRadius));
                            }
                        }
                        else {
                            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "only round or rectangular columns supported for bobbins");
                        }

                        turns[turnIndex] = turn;
                        placedTurnsCoordinates.push_back(additionalCoordinates);
                    }
                    currentBaseRadialHeightForLayers -= turnsInSection[0].get_dimensions().value()[0];
                }
                else {
                    maximumAdditionalRadialHeightPerInsulationLayerByIndex.push_back({layer, currentBaseRadialHeightForLayers});
                    currentBaseRadialHeightForLayers -= layer.get_dimensions()[0];
                }
            }

            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentSectionMaximumAdditionalRadialHeight -= turnsInSection[0].get_dimensions().value()[0] / 2;
                currentBaseRadialHeight = currentSectionMaximumAdditionalRadialHeight;
            }
        }
        else {
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                auto layersThisSection = get_layers_by_section(section.get_name());
                for (auto layer : layersThisSection) {
                    maximumAdditionalRadialHeightPerInsulationLayerByIndex.push_back({layer, currentBaseRadialHeight});
                    currentBaseRadialHeight -= layer.get_dimensions()[0];
                }
            }
        }

    }
    set_turns_description(turns);

    for (auto [layer, currentRadialHeight] : maximumAdditionalRadialHeightPerInsulationLayerByIndex) {
        if (layer.get_type() == ElectricalType::INSULATION) {
            auto layerIndex = get_layer_index_by_name(layer.get_name());
            currentRadialHeight -= layer.get_dimensions()[0] / 2;
            std::vector<double> additionalCoordinates = {currentRadialHeight, layer.get_coordinates()[1]};
            layers[layerIndex].set_additional_coordinates(std::vector<std::vector<double>>{additionalCoordinates});
            currentRadialHeight -= layer.get_dimensions()[0] / 2;
        }
    }
    set_layers_description(layers);

    return true;
}

std::vector<double> Coil::get_aligned_section_dimensions_rectangular_window(size_t sectionIndex) {
    auto sections = get_sections_description().value();
    if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
        sections[sectionIndex].set_margin(std::vector<double>{0, 0});
    }

    auto allWindingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    // Multi-column: redirect "[0]" to the section's group's window so the
    // alignment math uses that column's geometry. Single-column behaviour is
    // unchanged because find_window_index_for_group returns 0 in that case.
    size_t windowIndex = 0;
    auto sectionGroupOpt = sections[sectionIndex].get_group();
    if (sectionGroupOpt) {
        windowIndex = find_window_index_for_group(sectionGroupOpt.value());
    }
    if (windowIndex >= allWindingWindows.size()) windowIndex = 0;
    std::vector<WindingWindowElement> windingWindows = {allWindingWindows[windowIndex]};
    // Winding frame: groups are wound on the +x side of the main column and mirrored
    // into place at the end of wind(); align against the window's +x image.
    if (windingWindows[0].get_coordinates() && windingWindows[0].get_coordinates().value()[0] < 0) {
        auto windowCoordinates = windingWindows[0].get_coordinates().value();
        windowCoordinates[0] = -windowCoordinates[0];
        windingWindows[0].set_coordinates(windowCoordinates);
    }
    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();
    auto windingOrientation = get_winding_orientation();

    if (sections.size() == 0) {
        throw CoilNotProcessedException("No sections in coil");
    }
    // With more than one group, alignment totals only cover the sections sharing this
    // section's group (each group fills its own winding window).
    bool scopeToGroup = get_groups_description() && get_groups_description()->size() > 1;
    double totalSectionsWidth = 0;
    double totalSectionsHeight = 0;
    for (size_t auxSectionIndex = 0; auxSectionIndex < sections.size(); ++auxSectionIndex) {
        if (scopeToGroup && sections[auxSectionIndex].get_group() != sections[sectionIndex].get_group()) {
            continue;
        }
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            totalSectionsWidth += sections[auxSectionIndex].get_dimensions()[0];
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsHeight = std::max(totalSectionsHeight, sections[auxSectionIndex].get_dimensions()[1]);
            }
        }
        else {
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsWidth = std::max(totalSectionsWidth, sections[auxSectionIndex].get_dimensions()[0]);
            }
            totalSectionsHeight += sections[auxSectionIndex].get_dimensions()[1];
        }
    }

    double currentCoilWidth;
    double currentCoilHeight;
    double paddingAmongSectionWidth = 0;
    double paddingAmongSectionHeight = 0;
    auto turnsAlignment = get_turns_alignment(sections[sectionIndex].get_name());

    auto sectionAlignment = get_section_alignment();
    switch (sectionAlignment) {
        case CoilAlignment::INNER_OR_TOP:

            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilHeight = 0;
                            double currentCoilHeightTop = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                            double currentCoilHeightBottom = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                            currentCoilHeight = std::min(currentCoilHeight, currentCoilHeightTop);
                            currentCoilHeight = std::max(currentCoilHeight, currentCoilHeightBottom);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilHeight = -resolve_margin(sections[sectionIndex])[0] / 2 + resolve_margin(sections[sectionIndex])[1] / 2;
                        break;
                    default:
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
                }
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - sections[sectionIndex].get_dimensions()[0] / 2;
                            double currentCoilWidthLeft = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                            double currentCoilWidthRight = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
                            currentCoilWidth = std::max(currentCoilWidth, currentCoilWidthLeft);
                            currentCoilWidth = std::min(currentCoilWidth, currentCoilWidthRight);
                            break;
                        }
                    case CoilAlignment::SPREAD:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                        break;
                    default:
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
                }
            }
            break;
        case CoilAlignment::OUTER_OR_BOTTOM:
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - totalSectionsWidth;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilHeight = 0;
                            double currentCoilHeightTop = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                            double currentCoilHeightBottom = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                            currentCoilHeight = std::min(currentCoilHeight, currentCoilHeightTop);
                            currentCoilHeight = std::max(currentCoilHeight, currentCoilHeightBottom);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilHeight = -resolve_margin(sections[sectionIndex])[0] / 2 + resolve_margin(sections[sectionIndex])[1] / 2;
                        break;
                    default:
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
                }
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + totalSectionsHeight;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - sections[sectionIndex].get_dimensions()[0] / 2;
                            double currentCoilWidthLeft = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                            double currentCoilWidthRight = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
                            currentCoilWidth = std::max(currentCoilWidth, currentCoilWidthLeft);
                            currentCoilWidth = std::min(currentCoilWidth, currentCoilWidthRight);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                        break;
                    default:
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
                }
            }
            break;
        case CoilAlignment::SPREAD:
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilHeight = 0;
                            double currentCoilHeightTop = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                            double currentCoilHeightBottom = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                            currentCoilHeight = std::min(currentCoilHeight, currentCoilHeightTop);
                            currentCoilHeight = std::max(currentCoilHeight, currentCoilHeightBottom);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilHeight = -resolve_margin(sections[sectionIndex])[0] / 2 + resolve_margin(sections[sectionIndex])[1] / 2;
                        break;
                    default:
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
                }
                paddingAmongSectionWidth = windingWindows[0].get_width().value() - totalSectionsWidth;
                if (sections.size() > 1) {
                    paddingAmongSectionWidth /= sections.size() - 1;
                }
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2;
                paddingAmongSectionHeight = windingWindows[0].get_height().value() - totalSectionsHeight;
                if (sections.size() > 1) {
                    paddingAmongSectionHeight /= sections.size() - 1;
                }
                else {
                    currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + totalSectionsHeight / 2;
                }

                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - sections[sectionIndex].get_dimensions()[0] / 2;
                            double currentCoilWidthLeft = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                            double currentCoilWidthRight = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
                            currentCoilWidth = std::max(currentCoilWidth, currentCoilWidthLeft);
                            currentCoilWidth = std::min(currentCoilWidth, currentCoilWidthRight);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                        break;
                    default:
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
                }
            }
            break;
        case CoilAlignment::CENTERED:
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilHeight = 0;
                            double currentCoilHeightTop = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                            double currentCoilHeightBottom = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                            currentCoilHeight = std::min(currentCoilHeight, currentCoilHeightTop);
                            currentCoilHeight = std::max(currentCoilHeight, currentCoilHeightBottom);
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilHeight = -resolve_margin(sections[sectionIndex])[0] / 2 + resolve_margin(sections[sectionIndex])[1] / 2;
                        break;
                    default:
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
                }
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + totalSectionsHeight / 2;
                switch (turnsAlignment) {
                    case CoilAlignment::INNER_OR_TOP:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                        break;
                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
                        break;
                    case CoilAlignment::CENTERED:
                        {
                            currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - sections[sectionIndex].get_dimensions()[0] / 2;
                            double currentCoilWidthLeft = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                            double currentCoilWidthRight = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
                            if (currentCoilWidthLeft < 0) {
                                throw std::invalid_argument("currentCoilWidthLeft cannot be less than 0: " + std::to_string(currentCoilWidthLeft));
                            }
                            if (currentCoilWidthRight < 0) {
                                throw std::invalid_argument("currentCoilWidthRight cannot be less than 0: " + std::to_string(currentCoilWidthRight));
                            }
                            currentCoilWidth = std::max(currentCoilWidth, currentCoilWidthLeft);
                            if (currentCoilWidthRight >= 0) {
                                currentCoilWidth = std::min(currentCoilWidth, currentCoilWidthRight);
                            }
                            break;
                        }
                        break;
                    case CoilAlignment::SPREAD:
                        currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                        break;
                    default:
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
                }
            }
            break;
        default:
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");

    }

    return {currentCoilWidth, currentCoilHeight, paddingAmongSectionWidth, paddingAmongSectionHeight};
}

std::vector<double> Coil::get_aligned_section_dimensions_round_window(size_t sectionIndex) {
    auto sections = get_sections_description().value();
    if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
        sections[sectionIndex].set_margin(std::vector<double>{0, 0});
    }

    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    double windingWindowAngle = windingWindows[0].get_angle().value();
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
    auto windingOrientation = get_winding_orientation();

    if (sections.size() == 0) {
        throw CoilNotProcessedException("No sections in coil");
    }
    double totalSectionsRadialHeight = 0;
    double totalSectionsAngle = 0;
    for (size_t auxSectionIndex = 0; auxSectionIndex < sections.size(); ++auxSectionIndex) {
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            totalSectionsRadialHeight += sections[auxSectionIndex].get_dimensions()[0];
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsAngle = std::max(totalSectionsAngle, sections[auxSectionIndex].get_dimensions()[1]);
            }
        }
        else {
            double marginAngle0 = 0;
            double marginAngle1 = 0;
            if (sections[auxSectionIndex].get_type() == ElectricalType::CONDUCTION) {
                totalSectionsRadialHeight = std::max(totalSectionsRadialHeight, sections[auxSectionIndex].get_dimensions()[0]);
                double lastLayerMaximumRadius = windingWindowRadialHeight - (sections[auxSectionIndex].get_coordinates()[0] + sections[auxSectionIndex].get_dimensions()[0] / 2);
                marginAngle0 = wound_distance_to_angle(resolve_margin(sections[auxSectionIndex])[0], lastLayerMaximumRadius);
                marginAngle1 = wound_distance_to_angle(resolve_margin(sections[auxSectionIndex])[1], lastLayerMaximumRadius);
            }
            totalSectionsAngle += sections[auxSectionIndex].get_dimensions()[1] + marginAngle0 + marginAngle1;
        }
    }

    double currentCoilRadialHeight = 0;
    double currentCoilAngle;
    double paddingAmongSectionRadialHeight = 0;
    double paddingAmongSectionAngle = 0;
    double marginAngle0 = 0;

    if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
        double lastLayerMaximumRadius = windingWindowRadialHeight - (sections[sectionIndex].get_coordinates()[0] + sections[sectionIndex].get_dimensions()[0] / 2);
        marginAngle0 = wound_distance_to_angle(resolve_margin(sections[sectionIndex])[0], lastLayerMaximumRadius);
    }
    auto turnsAlignment = get_turns_alignment(sections[sectionIndex].get_name());

    if (windingOrientation == WindingOrientation::OVERLAPPING) {
        currentCoilRadialHeight = 0;
        switch (turnsAlignment) {
            case CoilAlignment::INNER_OR_TOP:
                currentCoilAngle = sections[sectionIndex].get_dimensions()[1] / 2;
                break;
            case CoilAlignment::OUTER_OR_BOTTOM:
                currentCoilAngle = windingWindowAngle - sections[sectionIndex].get_dimensions()[1] / 2;
                break;
            case CoilAlignment::CENTERED:
                {
                    currentCoilAngle = 180;
                    break;
                }
                break;
            case CoilAlignment::SPREAD:
                currentCoilAngle = sections[sectionIndex].get_dimensions()[1] / 2;
                break;
            default:
                throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
        }
    }
    else {
        auto sectionAlignment = get_section_alignment();
        switch (sectionAlignment) {
            case CoilAlignment::INNER_OR_TOP:
                currentCoilAngle = sections[sectionIndex].get_coordinates()[1] - sections[sectionIndex].get_dimensions()[1] / 2 - marginAngle0;
                break;
            case CoilAlignment::OUTER_OR_BOTTOM:
                currentCoilAngle = windingWindowAngle - totalSectionsAngle;
                break;
            case CoilAlignment::SPREAD:
                currentCoilAngle = sections[sectionIndex].get_coordinates()[1];
                paddingAmongSectionAngle = windingWindows[0].get_angle().value() - totalSectionsAngle;
                if (sections.size() > 1) {
                    paddingAmongSectionAngle /= sections.size() - 1;
                }
                else {
                    currentCoilAngle = windingWindowAngle / 2 + totalSectionsAngle / 2;
                }
                break;
            case CoilAlignment::CENTERED:
                currentCoilAngle = windingWindowAngle / 2 - totalSectionsAngle / 2;
                break;
            default:
                throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
        }
    }

    return {currentCoilRadialHeight, currentCoilAngle, paddingAmongSectionRadialHeight, paddingAmongSectionAngle};
}

bool Coil::delimit_and_compact() {

    auto bobbin = resolve_bobbin();

    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        return delimit_and_compact_rectangular_window();
    }
    else {
        return delimit_and_compact_round_window();
    }
}

WiringTechnology Coil::get_coil_type(size_t groupIndex) const {
    if (!get_groups_description()) {
        return WiringTechnology::WOUND;
    }
    auto groups = get_groups_description().value();
    if (groupIndex >= groups.size()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Non existing group index");
    }
    auto group = get_groups_description().value()[groupIndex];
    return group.get_type();
}

bool Coil::delimit_and_compact_rectangular_window() {
    // Delimit
    auto groupType = get_coil_type();

    if (!get_sections_description()) {
        throw CoilNotProcessedException("No sections to delimit");
    }

    if (get_layers_description()) {
        auto layers = get_layers_description().value();
        if (get_turns_description()) {
            for (size_t i = 0; i < layers.size(); ++i) {
                if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                    auto turnsInLayer = get_turns_by_layer(layers[i].get_name());
                    auto layerCoordinates = layers[i].get_coordinates();
                    double currentLayerMaximumWidth = (turnsInLayer[0].get_coordinates()[0] - layerCoordinates[0]) + turnsInLayer[0].get_dimensions().value()[0] / 2;
                    double currentLayerMinimumWidth = (turnsInLayer[0].get_coordinates()[0] - layerCoordinates[0]) - turnsInLayer[0].get_dimensions().value()[0] / 2;
                    double currentLayerMaximumHeight = (turnsInLayer[0].get_coordinates()[1] - layerCoordinates[1]) + turnsInLayer[0].get_dimensions().value()[1] / 2;
                    double currentLayerMinimumHeight = (turnsInLayer[0].get_coordinates()[1] - layerCoordinates[1]) - turnsInLayer[0].get_dimensions().value()[1] / 2;
                    for (auto& turn : turnsInLayer) {
                        currentLayerMaximumWidth = std::max(currentLayerMaximumWidth, (turn.get_coordinates()[0] - layerCoordinates[0]) + turn.get_dimensions().value()[0] / 2);
                        currentLayerMinimumWidth = std::min(currentLayerMinimumWidth, (turn.get_coordinates()[0] - layerCoordinates[0]) - turn.get_dimensions().value()[0] / 2);
                        currentLayerMaximumHeight = std::max(currentLayerMaximumHeight, (turn.get_coordinates()[1] - layerCoordinates[1]) + turn.get_dimensions().value()[1] / 2);
                        currentLayerMinimumHeight = std::min(currentLayerMinimumHeight, (turn.get_coordinates()[1] - layerCoordinates[1]) - turn.get_dimensions().value()[1] / 2);
                    }
                    if (groupType == WiringTechnology::PRINTED) {
                        layers[i].set_coordinates(std::vector<double>({layerCoordinates[0],
                                                                   layerCoordinates[1] + (currentLayerMaximumHeight + currentLayerMinimumHeight) / 2}));
                        layers[i].set_dimensions(std::vector<double>({layers[i].get_dimensions()[0],
                                                                   currentLayerMaximumHeight - currentLayerMinimumHeight}));
                    }
                    else {
                        layers[i].set_coordinates(std::vector<double>({layerCoordinates[0] + (currentLayerMaximumWidth + currentLayerMinimumWidth) / 2,
                                                                   layerCoordinates[1] + (currentLayerMaximumHeight + currentLayerMinimumHeight) / 2}));
                        layers[i].set_dimensions(std::vector<double>({currentLayerMaximumWidth - currentLayerMinimumWidth,
                                                                   currentLayerMaximumHeight - currentLayerMinimumHeight}));
                    }
                    if (i + 1 < layers.size()) {
                        layerCoordinates = layers[i + 1].get_coordinates();
                        if (layers[i + 1].get_type() == ElectricalType::INSULATION && layers[i + 1].get_section() == layers[i].get_section() && layers[i + 1].get_orientation() == WindingOrientation::CONTIGUOUS) {
                            layers[i + 1].set_coordinates(std::vector<double>({layerCoordinates[0] + (currentLayerMaximumWidth + currentLayerMinimumWidth) / 2,
                                                                       layerCoordinates[1]}));
                            layers[i + 1].set_dimensions(std::vector<double>({currentLayerMaximumWidth - currentLayerMinimumWidth, layers[i + 1].get_dimensions()[1]}));
                        }
                    }
                }
                set_layers_description(layers);
            }
        }

        auto sections = get_sections_description().value();
        for (size_t i = 0; i < sections.size(); ++i) {
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                auto layersInSection = get_layers_by_section(sections[i].get_name());
                if (layersInSection.size() == 0) {
                    throw CoilNotProcessedException("No layers in section: " + sections[i].get_name());
                }
                auto sectionCoordinates = sections[i].get_coordinates();
                double currentSectionMaximumWidth = (layersInSection[0].get_coordinates()[0] - sectionCoordinates[0]) + layersInSection[0].get_dimensions()[0] / 2;
                double currentSectionMinimumWidth = (layersInSection[0].get_coordinates()[0] - sectionCoordinates[0]) - layersInSection[0].get_dimensions()[0] / 2;
                double currentSectionMaximumHeight = (layersInSection[0].get_coordinates()[1] - sectionCoordinates[1]) + layersInSection[0].get_dimensions()[1] / 2;
                double currentSectionMinimumHeight = (layersInSection[0].get_coordinates()[1] - sectionCoordinates[1]) - layersInSection[0].get_dimensions()[1] / 2;

                for (auto& layer : layersInSection) {
                    if (layer.get_type() == ElectricalType::CONDUCTION) {
                        currentSectionMaximumWidth = std::max(currentSectionMaximumWidth, (layer.get_coordinates()[0] - sectionCoordinates[0]) + layer.get_dimensions()[0] / 2);
                        currentSectionMinimumWidth = std::min(currentSectionMinimumWidth, (layer.get_coordinates()[0] - sectionCoordinates[0]) - layer.get_dimensions()[0] / 2);
                        currentSectionMaximumHeight = std::max(currentSectionMaximumHeight, (layer.get_coordinates()[1] - sectionCoordinates[1]) + layer.get_dimensions()[1] / 2);
                        currentSectionMinimumHeight = std::min(currentSectionMinimumHeight, (layer.get_coordinates()[1] - sectionCoordinates[1]) - layer.get_dimensions()[1] / 2);
                    }
                }

                if (groupType == WiringTechnology::PRINTED) {
                    sections[i].set_coordinates(std::vector<double>({sectionCoordinates[0],
                                                               sectionCoordinates[1] + (currentSectionMaximumHeight + currentSectionMinimumHeight) / 2}));
                    sections[i].set_dimensions(std::vector<double>({sections[i].get_dimensions()[0],
                                                                    currentSectionMaximumHeight - currentSectionMinimumHeight}));
                }
                else {
                    sections[i].set_coordinates(std::vector<double>({sectionCoordinates[0] + (currentSectionMaximumWidth + currentSectionMinimumWidth) / 2,
                                                               sectionCoordinates[1] + (currentSectionMaximumHeight + currentSectionMinimumHeight) / 2}));
                    sections[i].set_dimensions(std::vector<double>({currentSectionMaximumWidth - currentSectionMinimumWidth,
                                                                    currentSectionMaximumHeight - currentSectionMinimumHeight}));
                }
            }
        }
        set_sections_description(sections);
    }

     // Compact
    if (get_sections_description()) {
        auto sections = get_sections_description().value();

        std::vector<std::vector<double>> alignedSectionDimensionsPerSection;

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            alignedSectionDimensionsPerSection.push_back(get_aligned_section_dimensions_rectangular_window(sectionIndex));
        }

        double currentCoilWidth = alignedSectionDimensionsPerSection[0][0];
        double currentCoilHeight = alignedSectionDimensionsPerSection[0][1];
        double paddingAmongSectionWidth = alignedSectionDimensionsPerSection[0][2];
        double paddingAmongSectionHeight = alignedSectionDimensionsPerSection[0][3];

        std::vector<Turn> turns;
        if (get_turns_description()) {
            turns = get_turns_description().value();
        }

        std::vector<Layer> layers;
        if (get_layers_description()) {
            layers = get_layers_description().value();
        }

        auto bobbinColumnShape = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_column_shape();

        // Wound-column frame per section for the turn-length recompute after the
        // compacting shift (multi-column winding support).
        std::map<std::string, WoundColumnFrame> woundColumnFramePerSection;
        auto getFrameForSection = [&](const std::string& sectionName) -> const WoundColumnFrame& {
            auto frameIterator = woundColumnFramePerSection.find(sectionName);
            if (frameIterator == woundColumnFramePerSection.end()) {
                frameIterator = woundColumnFramePerSection.emplace(sectionName, get_wound_column_frame_for_section(sectionName)).first;
            }
            return frameIterator->second;
        };

        auto windingOrientation = get_winding_orientation();
        bool multiWindowCompaction = get_groups_description() && get_groups_description()->size() > 1;

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            // A new group means a new winding window: restart the compaction cursor
            // from this section's aligned position (multi-column winding support).
            // Guarded on both sections carrying a group: sections without one (e.g.
            // planar insulation sections) must not break their neighbours' run.
            if (multiWindowCompaction && sectionIndex > 0 &&
                sections[sectionIndex].get_group() && sections[sectionIndex - 1].get_group() &&
                sections[sectionIndex].get_group().value() != sections[sectionIndex - 1].get_group().value()) {
                currentCoilWidth = alignedSectionDimensionsPerSection[sectionIndex][0];
                currentCoilHeight = alignedSectionDimensionsPerSection[sectionIndex][1];
                paddingAmongSectionWidth = alignedSectionDimensionsPerSection[sectionIndex][2];
                paddingAmongSectionHeight = alignedSectionDimensionsPerSection[sectionIndex][3];
            }
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilHeight = alignedSectionDimensionsPerSection[sectionIndex][1];
                currentCoilWidth += sections[sectionIndex].get_dimensions()[0] / 2;
            }
            else {
                currentCoilHeight -= sections[sectionIndex].get_dimensions()[1] / 2;
                currentCoilWidth = alignedSectionDimensionsPerSection[sectionIndex][0];
            }

            double compactingShiftWidth = sections[sectionIndex].get_coordinates()[0] - currentCoilWidth;
            double compactingShiftHeight = sections[sectionIndex].get_coordinates()[1] - currentCoilHeight;

            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                // compactingShiftWidth += sections[sectionIndex].get_dimensions()[0] / 2;
                if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
                    compactingShiftHeight = 0;
                }

            }
            else {
                compactingShiftWidth -= sections[sectionIndex].get_dimensions()[0] / 2;
                if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
                    compactingShiftWidth = 0;
                }
            }

            if (compactingShiftWidth != 0 || compactingShiftHeight != 0) {
                if (groupType == WiringTechnology::PRINTED) {
                    sections[sectionIndex].set_coordinates(std::vector<double>({
                        sections[sectionIndex].get_coordinates()[0],
                        sections[sectionIndex].get_coordinates()[1] - compactingShiftHeight
                    }));
                }
                else {
                    sections[sectionIndex].set_coordinates(std::vector<double>({
                        sections[sectionIndex].get_coordinates()[0] - compactingShiftWidth,
                        sections[sectionIndex].get_coordinates()[1] - compactingShiftHeight
                    }));
                }

                for (size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
                    if (layers[layerIndex].get_section().value() == sections[sectionIndex].get_name()){
                        if (groupType == WiringTechnology::PRINTED) {
                            layers[layerIndex].set_coordinates(std::vector<double>({
                                layers[layerIndex].get_coordinates()[0],
                                layers[layerIndex].get_coordinates()[1] - compactingShiftHeight
                            }));
                        }
                        else {
                            layers[layerIndex].set_coordinates(std::vector<double>({
                                layers[layerIndex].get_coordinates()[0] - compactingShiftWidth,
                                layers[layerIndex].get_coordinates()[1] - compactingShiftHeight
                            }));
                        }
                        for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {

                            if (turns[turnIndex].get_layer().value() == layers[layerIndex].get_name()){


                                turns[turnIndex].set_coordinate_system(CoordinateSystem::CARTESIAN);
        
                                if (bobbinColumnShape == ColumnShape::ROUND || bobbinColumnShape == ColumnShape::OBLONG || bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                    if (turns[turnIndex].get_coordinates()[0] < compactingShiftWidth) {
                                        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Something wrong happened with compactingShiftWidth: " + std::to_string(compactingShiftWidth) +
                                                                 "\nsections[sectionIndex].get_coordinates()[0]: " + std::to_string(sections[sectionIndex].get_coordinates()[0]) +
                                                                 "\ncurrentCoilWidth: " + std::to_string(currentCoilWidth) +
                                                                 "\nturns[turnIndex].get_coordinates()[0]: " + std::to_string(turns[turnIndex].get_coordinates()[0])
                                                                 );
                                    }
                                }
                                else {
                                    throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "only round or rectangular columns supported for bobbins");
                                }

                                if (groupType == WiringTechnology::PRINTED) {
                                    turns[turnIndex].set_coordinates(std::vector<double>({
                                        turns[turnIndex].get_coordinates()[0],
                                        turns[turnIndex].get_coordinates()[1] - compactingShiftHeight
                                    }));
                                }
                                else {
                                    turns[turnIndex].set_coordinates(std::vector<double>({
                                        turns[turnIndex].get_coordinates()[0] - compactingShiftWidth,
                                        turns[turnIndex].get_coordinates()[1] - compactingShiftHeight
                                    }));
                                }

                                {
                                    auto turnLength = get_turn_length_in_frame(getFrameForSection(sections[sectionIndex].get_name()), turns[turnIndex].get_coordinates()[0]);
                                    if (!turnLength) {
                                        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Something wrong happened in turn length 1: negative length for turn " + turns[turnIndex].get_name() + " at x: " + std::to_string(turns[turnIndex].get_coordinates()[0]));
                                    }
                                    turns[turnIndex].set_length(turnLength.value());
                                }
                            }
                        }
                    }
                }
            }
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth += sections[sectionIndex].get_dimensions()[0] / 2 + paddingAmongSectionWidth;
            }
            else {
                currentCoilHeight -= sections[sectionIndex].get_dimensions()[1] / 2 + paddingAmongSectionHeight;
            }
        }
        if (get_turns_description()) {
            set_turns_description(turns);
        }
        if (get_layers_description()) {
            set_layers_description(layers);
        }
        set_sections_description(sections);
    }

    // Add extra margin for support if required
    bool fillCoilSectionsWithMarginTape = settings.get_coil_fill_sections_with_margin_tape();

    // Compact groups in planar case

    if (get_layers_description() && groupType == WiringTechnology::PRINTED) {

        if (get_groups_description()) {
            auto groups = get_groups_description().value();
            for (size_t i = 0; i < groups.size(); ++i) {
                auto sectionsInGroup = get_sections_by_group(groups[i].get_name());
                if (sectionsInGroup.size() == 0) {
                    throw CoilNotProcessedException("No sections in group: " + groups[i].get_name());
                }
                auto groupCoordinates = groups[i].get_coordinates();
                double currentGroupMaximumHeight = (sectionsInGroup[0].get_coordinates()[1] - groupCoordinates[1]) + sectionsInGroup[0].get_dimensions()[1] / 2;
                double currentGroupMinimumHeight = (sectionsInGroup[0].get_coordinates()[1] - groupCoordinates[1]) - sectionsInGroup[0].get_dimensions()[1] / 2;

                for (auto& section : sectionsInGroup) {
                    currentGroupMaximumHeight = std::max(currentGroupMaximumHeight, (section.get_coordinates()[1] - groupCoordinates[1]) + section.get_dimensions()[1] / 2);
                    currentGroupMinimumHeight = std::min(currentGroupMinimumHeight, (section.get_coordinates()[1] - groupCoordinates[1]) - section.get_dimensions()[1] / 2);
                }
                groups[i].set_coordinates(std::vector<double>({groupCoordinates[0], groupCoordinates[1] + (currentGroupMaximumHeight + currentGroupMinimumHeight) / 2}));
                groups[i].set_dimensions(std::vector<double>({groups[i].get_dimensions()[0], currentGroupMaximumHeight - currentGroupMinimumHeight}));
            }
            set_groups_description(groups);
        }
    }

    if (fillCoilSectionsWithMarginTape) {
        auto bobbin = resolve_bobbin();
        auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
        auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
        double windingWindowHeight = windingWindowDimensions[1];
        double windingWindowWidth = windingWindowDimensions[0];
        auto sections = get_sections_description().value();
        for (size_t i = 0; i < sections.size(); ++i) {
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                auto sectionOrientation = bobbin.get_winding_window_sections_orientation(0);
                if (sectionOrientation == WindingOrientation::OVERLAPPING) {
                    auto topSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[1] + windingWindowHeight / 2) - (sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2));
                    auto bottomSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[1] - windingWindowHeight / 2) - (sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2));
                    sections[i].set_margin(std::vector<double>{topSpaceBetweenSectionAndBobbin, bottomSpaceBetweenSectionAndBobbin});
                }
                else if (sectionOrientation == WindingOrientation::CONTIGUOUS) {
                    auto innerSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[0] - windingWindowWidth / 2) - (sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2));
                    auto outerSpaceBetweenSectionAndBobbin = fabs((windingWindowCoordinates[0] + windingWindowWidth / 2) - (sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2));
                    sections[i].set_margin(std::vector<double>{innerSpaceBetweenSectionAndBobbin, outerSpaceBetweenSectionAndBobbin});
                }
            }
        }
        set_sections_description(sections);
    }

    return true;
}

bool Coil::delimit_and_compact_round_window() {
    if (get_turns_description()) {
        convert_turns_to_polar_coordinates();
    }

    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto windingWindowsRadius = windingWindows[0].get_radial_height().value();


    // Radial Delimit
    if (get_layers_description()) {
        auto layers = get_layers_description().value();
        if (get_turns_description()) {
            for (size_t i = 0; i < layers.size(); ++i) {
                if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                    auto turnsInLayer = get_turns_by_layer(layers[i].get_name());
                    auto layerCoordinates = layers[i].get_coordinates();
                    auto section = get_section_by_name(layers[i].get_section().value());

                    double currentLayerMaximumRadialHeight = (turnsInLayer[0].get_coordinates()[0] - layerCoordinates[0]) + turnsInLayer[0].get_dimensions().value()[0] / 2;
                    double currentLayerMinimumRadialHeight = (turnsInLayer[0].get_coordinates()[0] - layerCoordinates[0]) - turnsInLayer[0].get_dimensions().value()[0] / 2;

                    for (auto& turn : turnsInLayer) {
                        currentLayerMaximumRadialHeight = std::max(currentLayerMaximumRadialHeight, (turn.get_coordinates()[0] - layerCoordinates[0]) + turn.get_dimensions().value()[0] / 2);
                        currentLayerMinimumRadialHeight = std::min(currentLayerMinimumRadialHeight, (turn.get_coordinates()[0] - layerCoordinates[0]) - turn.get_dimensions().value()[0] / 2);
                    }

                    layers[i].set_coordinates(std::vector<double>({layerCoordinates[0] + (currentLayerMaximumRadialHeight + currentLayerMinimumRadialHeight) / 2, layers[i].get_coordinates()[1]}));
                    layers[i].set_dimensions(std::vector<double>({currentLayerMaximumRadialHeight - currentLayerMinimumRadialHeight, layers[i].get_dimensions()[1]}));
                }
                set_layers_description(layers);
            }
        }

        auto sections = get_sections_description().value();
        for (size_t i = 0; i < sections.size(); ++i) {
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                auto layersInSection = get_layers_by_section(sections[i].get_name());
                if (layersInSection.size() == 0) {
                    throw CoilNotProcessedException("No layers in section: " + sections[i].get_name());
                }
                auto sectionCoordinates = sections[i].get_coordinates();
                double currentSectionMaximumRadialHeight = (layersInSection[0].get_coordinates()[0] - sectionCoordinates[0]) + layersInSection[0].get_dimensions()[0] / 2;
                double currentSectionMinimumRadialHeight = (layersInSection[0].get_coordinates()[0] - sectionCoordinates[0]) - layersInSection[0].get_dimensions()[0] / 2;

                for (auto& layer : layersInSection) {
                    currentSectionMaximumRadialHeight = std::max(currentSectionMaximumRadialHeight, (layer.get_coordinates()[0] - sectionCoordinates[0]) + layer.get_dimensions()[0] / 2);
                    currentSectionMinimumRadialHeight = std::min(currentSectionMinimumRadialHeight, (layer.get_coordinates()[0] - sectionCoordinates[0]) - layer.get_dimensions()[0] / 2);
                }
                sections[i].set_coordinates(std::vector<double>({sectionCoordinates[0] + (currentSectionMaximumRadialHeight + currentSectionMinimumRadialHeight) / 2,
                                                           sections[i].get_coordinates()[1]}));
                sections[i].set_dimensions(std::vector<double>({currentSectionMaximumRadialHeight - currentSectionMinimumRadialHeight,
                                                       sections[i].get_dimensions()[1]}));
            }
        }
        set_sections_description(sections);
    }

    // Angular Delimit
    if (get_layers_description()) {
        auto wirePerWinding = get_wires();
        auto layers = get_layers_description().value();
        if (get_turns_description()) {
            for (size_t i = 0; i < layers.size(); ++i) {
                if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                    auto turnsInLayer = get_turns_by_layer(layers[i].get_name());
                    auto layerCoordinates = layers[i].get_coordinates();
                    auto section = get_section_by_name(layers[i].get_section().value());

                    auto windingIndex = get_winding_index_by_name(turnsInLayer[0].get_winding());
                    double wireWidth = wirePerWinding[windingIndex].get_maximum_outer_width();

                    double wireRadius;
                    if (wirePerWinding[windingIndex].get_type() == WireType::RECTANGULAR) {
                        wireRadius = windingWindows[0].get_radial_height().value() - turnsInLayer[0].get_coordinates()[0] - wireWidth / 2;
                    }
                    else {
                        wireRadius = windingWindows[0].get_radial_height().value() - turnsInLayer[0].get_coordinates()[0];
                    }

                    double turnDimensionAngle = wound_distance_to_angle(turnsInLayer[0].get_dimensions().value()[1], wireRadius);

                    // for (auto& turn : turnsInLayer) {
                        // turnDimensionAngle = wound_distance_to_angle(turn.get_dimensions().value()[1], windingWindows[0].get_radial_height().value() - turn.get_coordinates()[0]);
                    // }
                    double layerAngle = turnDimensionAngle * turnsInLayer.size();
                    double layerCenterAngle = 0;
                    // double marginAngle0 = 0;
                    // double marginAngle1 = 0;

                    // if (section.get_type() == ElectricalType::CONDUCTION) {
                    //     double lastLayerMaximumRadius = windingWindowsRadius - (section.get_coordinates()[0] + section.get_dimensions()[0] / 2);
                    //     marginAngle0 = wound_distance_to_angle(resolve_margin(section)[0], lastLayerMaximumRadius);
                    //     marginAngle1 = wound_distance_to_angle(resolve_margin(section)[1], lastLayerMaximumRadius);
                    // }

                    switch (layers[i].get_turns_alignment().value()) {
                        case CoilAlignment::INNER_OR_TOP:
                            layerCenterAngle = section.get_coordinates()[1] - section.get_dimensions()[1] / 2 + layerAngle / 2;
                            break;
                        case CoilAlignment::OUTER_OR_BOTTOM:
                            layerCenterAngle = section.get_coordinates()[1] + section.get_dimensions()[1] / 2 - layerAngle / 2;
                            break;
                        case CoilAlignment::CENTERED:
                            layerCenterAngle = section.get_coordinates()[1];
                            break;
                        case CoilAlignment::SPREAD:
                            layerCenterAngle = section.get_coordinates()[1];
                            layerAngle = section.get_dimensions()[1];
                            break;
                        default:
                            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
                    }
                    layers[i].set_coordinates(std::vector<double>({layers[i].get_coordinates()[0], layerCenterAngle}));
                    layers[i].set_dimensions(std::vector<double>({layers[i].get_dimensions()[0], layerAngle}));
                }
                set_layers_description(layers);
            }
        }

        auto sections = get_sections_description().value();
        for (size_t i = 0; i < sections.size(); ++i) {
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                auto layersInSection = get_layers_by_section(sections[i].get_name());
                if (layersInSection.size() == 0) {
                    throw CoilNotProcessedException("No layers in section: " + sections[i].get_name());
                }
                auto sectionCoordinates = sections[i].get_coordinates();
                double currentSectionMaximumAngle = (layersInSection[0].get_coordinates()[1] - sectionCoordinates[1]) + layersInSection[0].get_dimensions()[1] / 2;
                double currentSectionMinimumAngle = (layersInSection[0].get_coordinates()[1] - sectionCoordinates[1]) - layersInSection[0].get_dimensions()[1] / 2;

                for (auto& layer : layersInSection) {
                    currentSectionMaximumAngle = std::max(currentSectionMaximumAngle, (layer.get_coordinates()[1] - sectionCoordinates[1]) + layer.get_dimensions()[1] / 2);
                    currentSectionMinimumAngle = std::min(currentSectionMinimumAngle, (layer.get_coordinates()[1] - sectionCoordinates[1]) - layer.get_dimensions()[1] / 2);
                }
                sections[i].set_coordinates(std::vector<double>({sections[i].get_coordinates()[0], sectionCoordinates[1] + (currentSectionMaximumAngle + currentSectionMinimumAngle) / 2}));
                sections[i].set_dimensions(std::vector<double>({sections[i].get_dimensions()[0], currentSectionMaximumAngle - currentSectionMinimumAngle}));
            }
        }
        set_sections_description(sections);
    }

     // Angular Compact
    if (get_sections_description()) {
        auto sections = get_sections_description().value();

        std::vector<std::vector<double>> alignedSectionDimensionsPerSection;

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            alignedSectionDimensionsPerSection.push_back(get_aligned_section_dimensions_round_window(sectionIndex));
        }

        double currentCoilAngle = alignedSectionDimensionsPerSection[0][1];
        double paddingAmongSectionAngle = alignedSectionDimensionsPerSection[0][3];
        std::vector<Turn> turns;
        if (get_turns_description()) {
            turns = get_turns_description().value();
        }

        std::vector<Layer> layers;
        if (get_layers_description()) {
            layers = get_layers_description().value();
        }

        auto bobbinColumnShape = bobbin.get_processed_description().value().get_column_shape();
        auto bobbinColumnDepth = bobbin.get_processed_description().value().get_column_depth();
        double bobbinColumnWidth;
        if (bobbin.get_processed_description().value().get_column_width()) {
            bobbinColumnWidth = bobbin.get_processed_description().value().get_column_width().value();
        }
        else {
            auto bobbinWindingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
            double bobbinWindingWindowWidth = bobbinWindingWindow.get_width().value();
            double bobbinWindingWindowCenterWidth = bobbinWindingWindow.get_coordinates().value()[0];
            bobbinColumnWidth = bobbinWindingWindowCenterWidth - bobbinWindingWindowWidth / 2;
        }

        auto windingOrientation = get_winding_orientation();

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {

            double marginAngle0 = 0;
            double marginAngle1 = 0;

            if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
                double lastLayerMaximumRadius = windingWindowsRadius - (sections[sectionIndex].get_coordinates()[0] + sections[sectionIndex].get_dimensions()[0] / 2);
                marginAngle0 = wound_distance_to_angle(resolve_margin(sections[sectionIndex])[0], lastLayerMaximumRadius);
                marginAngle1 = wound_distance_to_angle(resolve_margin(sections[sectionIndex])[1], lastLayerMaximumRadius);
            }


            auto sectionAlignment = get_section_alignment();
            if (windingOrientation == WindingOrientation::OVERLAPPING || sectionAlignment == CoilAlignment::SPREAD) {
                currentCoilAngle = alignedSectionDimensionsPerSection[sectionIndex][1];
            }
            else {
                currentCoilAngle += sections[sectionIndex].get_dimensions()[1] / 2 + marginAngle0;
            }

            double compactingShiftAngle = sections[sectionIndex].get_coordinates()[1] - currentCoilAngle;

            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                if (sections[sectionIndex].get_type() == ElectricalType::INSULATION) {
                    compactingShiftAngle = 0;
                }
            }

            sections[sectionIndex].set_coordinates(std::vector<double>({
                sections[sectionIndex].get_coordinates()[0],
                sections[sectionIndex].get_coordinates()[1] - compactingShiftAngle
            }));

            for (size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
                if (layers[layerIndex].get_section().value() == sections[sectionIndex].get_name()){
                    layers[layerIndex].set_coordinates(std::vector<double>({
                        layers[layerIndex].get_coordinates()[0],
                        layers[layerIndex].get_coordinates()[1] - compactingShiftAngle
                    }));
                    size_t turnInThisLayerIndex = 0;
                    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
                        if (turns[turnIndex].get_layer().value() == layers[layerIndex].get_name()){

                            turns[turnIndex].set_coordinates(std::vector<double>({
                                turns[turnIndex].get_coordinates()[0],
                                turns[turnIndex].get_coordinates()[1] - compactingShiftAngle
                            }));


                            if (bobbinColumnShape == ColumnShape::ROUND) {
                                turns[turnIndex].set_length(2 * std::numbers::pi * (turns[turnIndex].get_coordinates()[0] + bobbinColumnWidth));
                                if (turns[turnIndex].get_length() < 0) {
                                    return false;
                                }
                            }
                            else if (bobbinColumnShape == ColumnShape::OBLONG) {
                                turns[turnIndex].set_length(2 * std::numbers::pi * (turns[turnIndex].get_coordinates()[0] + bobbinColumnWidth) + 4 * (bobbinColumnDepth - bobbinColumnWidth));
                                if (turns[turnIndex].get_length() < 0) {
                                    return false;
                                }
                            }
                            else if (bobbinColumnShape == ColumnShape::RECTANGULAR || bobbinColumnShape == ColumnShape::IRREGULAR) {
                                double currentTurnCornerRadius = turns[turnIndex].get_coordinates()[0];
                                turns[turnIndex].set_length(4 * bobbinColumnDepth + 4 * bobbinColumnWidth + 2 * std::numbers::pi * currentTurnCornerRadius);
                                if (turns[turnIndex].get_length() < 0) {
                                    return false;
                                }
                            }
                            else {
                                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "only round or rectangular columns supported for bobbins");
                            }

                            turnInThisLayerIndex++;
                        }
                    }
                }
                }
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
            }
            else {
                currentCoilAngle += sections[sectionIndex].get_dimensions()[1] / 2 + paddingAmongSectionAngle + marginAngle1;
            }
        }
        if (get_turns_description()) {
            set_turns_description(turns);
        }
        if (get_layers_description()) {
            set_layers_description(layers);
        }
        set_sections_description(sections);
        
        if (settings.get_coil_include_additional_coordinates()) {
            wind_toroidal_additional_turns();
        }
    }

    if (get_turns_description()) {
        convert_turns_to_cartesian_coordinates();
    }
    return true;
}

std::vector<Wire> Coil::get_wires() {
    std::vector<Wire> wirePerWinding;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        Wire wire = resolve_wire(get_functional_description()[windingIndex]);
        wirePerWinding.push_back(wire);
    }
    return wirePerWinding;
}

void Coil::set_wires(std::vector<Wire> wires) {
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        get_mutable_functional_description()[windingIndex].set_wire(wires[windingIndex]);
    }
}

Wire Coil::resolve_wire(size_t windingIndex) {
    return resolve_wire(get_functional_description()[windingIndex]);
}

Wire Winding::resolve_wire() {
    auto wireOrString = get_wire();
    Wire wire;
    if (std::holds_alternative<std::string>(wireOrString)) {
        try {
            wire = find_wire_by_name(std::get<std::string>(wireOrString));
        }
        catch (const std::exception &exc) {
            // If wire is not found because it is "Dummy", we return a small Round, as it should only happening when get an advised wire
            if (std::get<std::string>(wireOrString) == "Dummy") {
                wire = find_wire_by_name("Round 0.01 - Grade 1");
            }
            else {
                throw WireNotFoundException("wire not found: " + std::get<std::string>(wireOrString));
            }
        }
    }
    else {
        wire = std::get<OpenMagnetics::Wire>(wireOrString);
    }
    return wire;
}

Wire Coil::resolve_wire(Winding winding) {
    return winding.resolve_wire();
}

std::vector<double> Coil::get_wires_length() const {
    std::vector<double> wiresLength;
    if (!get_turns_description()) {
        throw CoilNotProcessedException("Missing turns");
    }
    for (auto winding : get_functional_description()) {
        auto turns = get_turns_by_winding(winding.get_name());
        double wireLength = 0;
        for (auto turn : turns) {
            wireLength += turn.get_length();
        }
        wiresLength.push_back(wireLength);
    }
    return wiresLength;
}

WireType Coil::get_wire_type(Winding winding) {
    return resolve_wire(winding).get_type();
}

WireType Coil::get_wire_type(size_t windingIndex) {
    return get_wire_type(get_functional_description()[windingIndex]);
}

std::string Coil::get_wire_name(Winding winding) {
    auto name = resolve_wire(winding).get_name();
    if (name) {
        return name.value();
    }
    else {
        return "Custom";
    }
}

std::string Coil::get_wire_name(size_t windingIndex) {
    return get_wire_name(get_functional_description()[windingIndex]);
}

Bobbin Coil::resolve_bobbin(Coil coil) { 
    return coil.resolve_bobbin();
}

Bobbin Coil::resolve_bobbin() {
    if (_bobbin_resolved) {
        return _bobbin;
    }

    auto bobbinDataOrNameUnion = get_bobbin();
    if (std::holds_alternative<std::string>(bobbinDataOrNameUnion)) {
        if (std::get<std::string>(bobbinDataOrNameUnion) == "Dummy")
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Bobbin is dummy");

        auto bobbin = find_bobbin_by_name(std::get<std::string>(bobbinDataOrNameUnion));
        _bobbin = bobbin;
        return bobbin;
    }
    else {
        _bobbin = Bobbin(std::get<Bobbin>(bobbinDataOrNameUnion));
        return _bobbin;
    }
}

size_t Coil::convert_conduction_section_index_to_global(size_t conductionSectionIndex) {
    size_t currentConductionSectionIndex = 0;
    if (!get_sections_description()) {
        throw CoilNotProcessedException("In Convert Conduction Sections: Section description empty, wind coil first");
    }
    auto sections = get_sections_description().value();
    for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
        if (sections[sectionIndex].get_type() == ElectricalType::CONDUCTION) {
            if (currentConductionSectionIndex == conductionSectionIndex) {
                return sectionIndex;
            }
            currentConductionSectionIndex++;
        }
    }
    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Index not found");
}

void Coil::clear() {
    set_groups_description(std::nullopt);
    set_sections_description(std::nullopt);
    set_layers_description(std::nullopt);
    set_turns_description(std::nullopt);
}

void Coil::try_rewind() {
    if (!get_sections_description()) {
        return;
    }
    if (!get_layers_description()) {
        return;
    }
 
    if (!get_turns_description()) {
        wind_by_turns();
        delimit_and_compact();
    }
    auto electricalSections = get_sections_by_type(ElectricalType::CONDUCTION);

    if (electricalSections.size() == 1 || get_functional_description().size() == 1) {
        return;
    }

    bool windEvenIfNotFit = settings.get_coil_wind_even_if_not_fit();
    bool delimitAndCompact = settings.get_coil_delimit_and_compact();

    auto sections = get_sections_description().value();
    std::vector<double> extraSpaceNeededPerSection;
    double totalExtraSpaceNeeded = 0;
    auto bobbin = resolve_bobbin();
    auto sectionOrientation = bobbin.get_winding_window_sections_orientation(0);
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions(0);
    double windingWindowRemainingRestrictiveDimension;
    double windingWindowRemainingRestrictiveDimensionAccordingToSections;
    double windingWindowRestrictiveDimension;
    if (sectionOrientation == WindingOrientation::OVERLAPPING) {
        windingWindowRemainingRestrictiveDimensionAccordingToSections = windingWindowDimensions[0];
        windingWindowRemainingRestrictiveDimension = windingWindowDimensions[0];
        windingWindowRestrictiveDimension = windingWindowDimensions[0];
    }
    else {
        windingWindowRemainingRestrictiveDimensionAccordingToSections = windingWindowDimensions[1];
        windingWindowRemainingRestrictiveDimension = windingWindowDimensions[1];
        windingWindowRestrictiveDimension = windingWindowDimensions[1];
    }


    for (auto& section : sections) {
        if (section.get_type() == ElectricalType::INSULATION) {
            if (sectionOrientation == WindingOrientation::OVERLAPPING) {
                windingWindowRestrictiveDimension -= section.get_dimensions()[0];
            }
            else {
                windingWindowRestrictiveDimension -= section.get_dimensions()[1];
            }
        }
    }

    for (auto& section : sections) {
        double sectionRestrictiveDimension;
        double layersRestrictiveDimension = 0;
        double sectionFillingFactor;
        double extraSpaceNeededThisSection = 0;

        auto layers = get_layers_by_section(section.get_name());
        if (sectionOrientation == WindingOrientation::OVERLAPPING) {
            if (section.get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                for (auto layer : layers) {
                    double layerRestrictiveDimension = layer.get_dimensions()[0];
                    double layerFillingFactor = layer.get_filling_factor().value();
                    layersRestrictiveDimension += layerRestrictiveDimension;

                    extraSpaceNeededThisSection += std::max(0.0, (layerFillingFactor - 1) * layerRestrictiveDimension);
                    windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
                }
            }
            if (section.get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                double layerRestrictiveDimension = 0;
                double layerFillingFactor = 0;
                for (auto layer : layers) {
                    layerRestrictiveDimension = std::max(layerRestrictiveDimension, layer.get_dimensions()[0]);
                    layerFillingFactor = std::max(layerFillingFactor, layer.get_filling_factor().value());
                }
                layersRestrictiveDimension = layerRestrictiveDimension;
                extraSpaceNeededThisSection += std::max(0.0, (layerFillingFactor - 1) * layerRestrictiveDimension);
                windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
            }
        }
        else if (sectionOrientation == WindingOrientation::CONTIGUOUS) {
            if (section.get_layers_orientation() == WindingOrientation::OVERLAPPING) {
                double layerRestrictiveDimension = 0;
                double layerFillingFactor = 0;
                for (auto layer : layers) {
                    layerRestrictiveDimension = std::max(layerRestrictiveDimension, layer.get_dimensions()[1]);
                    layerFillingFactor = std::max(layerFillingFactor, layer.get_filling_factor().value());
                }

                layersRestrictiveDimension = layerRestrictiveDimension;
                extraSpaceNeededThisSection += std::max(0.0, (layerFillingFactor - 1) * layerRestrictiveDimension);
                windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
            }
            if (section.get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                for (auto layer : layers) {
                    double layerRestrictiveDimension = layer.get_dimensions()[1];
                    double layerFillingFactor = layer.get_filling_factor().value();
                    layersRestrictiveDimension += layerRestrictiveDimension;

                    extraSpaceNeededThisSection += std::max(0.0, (layerFillingFactor - 1) * layerRestrictiveDimension);
                    windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
                }
            }
        }

        if (sectionOrientation == WindingOrientation::OVERLAPPING) {
            sectionRestrictiveDimension = section.get_dimensions()[0];
            sectionFillingFactor = overlapping_filling_factor(section);
        }
        else {
            sectionRestrictiveDimension = section.get_dimensions()[1];
            sectionFillingFactor = contiguous_filling_factor(section);
        }
        windingWindowRemainingRestrictiveDimensionAccordingToSections -= sectionRestrictiveDimension;


        extraSpaceNeededThisSection = std::max(extraSpaceNeededThisSection, (sectionFillingFactor - 1) * sectionRestrictiveDimension);
        if (extraSpaceNeededThisSection < 0 || std::isnan(extraSpaceNeededThisSection)) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "extraSpaceNeededThisSection cannot be negative or nan: " + std::to_string(extraSpaceNeededThisSection));
        }
        if (section.get_type() == ElectricalType::CONDUCTION) {
        }
        extraSpaceNeededPerSection.push_back(extraSpaceNeededThisSection);
        totalExtraSpaceNeeded += extraSpaceNeededThisSection;
    }

    if (windingWindowRemainingRestrictiveDimensionAccordingToSections <= 0 || totalExtraSpaceNeeded <= 0) {
        return;
    }

    std::vector<double> newProportions;
    double numberWindings = get_functional_description().size();

    if (totalExtraSpaceNeeded < 0 || std::isnan(totalExtraSpaceNeeded)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "totalExtraSpaceNeeded cannot be negative or nan: " + std::to_string(totalExtraSpaceNeeded));
    }

    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        // Windings grouped via wound_with share sections; only the group's
        // representative (minimum index) accumulates the shared space. Non-
        // representatives contribute zero so that virtualize_proportion_per_winding
        // doesn't double-count the shared sections when collapsing the group.
        if (get_winding_group_minimum_index(windingIndex) != windingIndex) {
            newProportions.push_back(0.0);
            continue;
        }
        // double currentProportion = _currentProportionPerWinding[windingIndex];
        double currentSpace = 0;
        double extraSpaceNeededThisWinding = 0;

        for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex) {
            for (auto & winding : sections[sectionIndex].get_partial_windings()) {
                if (winding.get_winding() == get_functional_description()[windingIndex].get_name()) {
                    if (sectionOrientation == WindingOrientation::OVERLAPPING) {
                        currentSpace += sections[sectionIndex].get_dimensions()[0];

                        // We need to add half the insulation space after it, in case there is
                        if (sectionIndex + 1 < sections.size()) {
                            if (sections[sectionIndex + 1].get_type() == ElectricalType::INSULATION) {
                                // throw std::runtime_error("Consecutive layer to CONDUCTION must always be INSULATION");
                                if (sectionIndex == 0) {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[0] / 2;
                                }
                                else if (sectionIndex == sections.size() - 2) {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[0] * 3 / 2;
                                }
                                else {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[0];
                                }
                            }
                        }
                    }
                    else {
                        currentSpace += sections[sectionIndex].get_dimensions()[1];

                        // We need to add half the insulation space after it, in case there is
                        if (sectionIndex + 1 < sections.size()) {
                            if (sections[sectionIndex + 1].get_type() == ElectricalType::INSULATION) {
                                // throw std::runtime_error("Consecutive layer to CONDUCTION must always be INSULATION");
                                if (sectionIndex == 0 || sectionIndex == sections.size() - 2) {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[1] / 2;
                                }
                                else {
                                    currentSpace += sections[sectionIndex + 1].get_dimensions()[1];
                                }
                            }
                        }
                    }

                    extraSpaceNeededThisWinding += extraSpaceNeededPerSection[sectionIndex];
                    continue;
                }
            }
        }
        if (extraSpaceNeededThisWinding < 0 || std::isnan(extraSpaceNeededThisWinding)) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "extraSpaceNeededThisWinding cannot be negative or nan: " + std::to_string(extraSpaceNeededThisWinding));
        }
        // double proportionOfNeededForThisWinding = extraSpaceNeededThisWinding / totalExtraSpaceNeeded;
        double extraSpaceGottenByThisWinding = windingWindowRemainingRestrictiveDimensionAccordingToSections * extraSpaceNeededThisWinding / totalExtraSpaceNeeded;
        double newSpaceGottenByThisWinding = currentSpace + extraSpaceGottenByThisWinding;
        double newProportionGottenByThisWinding = newSpaceGottenByThisWinding / windingWindowRestrictiveDimension;

        if (extraSpaceGottenByThisWinding < 0 || std::isnan(extraSpaceGottenByThisWinding)) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "extraSpaceGottenByThisWinding cannot be negative or nan: " + std::to_string(extraSpaceGottenByThisWinding));
        }
        if (newProportionGottenByThisWinding < 0 || std::isnan(newProportionGottenByThisWinding)) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "newProportionGottenByThisWinding cannot be negative or nan: " + std::to_string(newProportionGottenByThisWinding));
        }
        if (roundFloat(newProportionGottenByThisWinding, 6) > 1 || std::isnan(newProportionGottenByThisWinding)) {
            return;
            // throw std::runtime_error("newProportionGottenByThisWinding cannot be greater than 1 or nan: " + std::to_string(newProportionGottenByThisWinding));
        }

        newProportions.push_back(newProportionGottenByThisWinding);
    }

    wind_by_sections(newProportions, _currentPattern, _currentRepetitions);



    wind_by_layers();

    if (!get_layers_description()) {
        return;
    }
    // set_turns_description(std::nullopt);
    if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
        wind_by_turns();
        if (delimitAndCompact) {
            delimit_and_compact();
        }
    }
}

void Coil::preload_margins(std::vector<std::vector<double>> marginPairs) {
    for (auto margins : marginPairs) {
        _marginsPerSection.push_back(margins);
        // Add an extra one for the insulation layer
        _marginsPerSection.push_back(margins);
    }
}

void Coil::add_margin_to_section_by_index(size_t sectionIndex, std::vector<double> margins) {
    if (!get_sections_description()) {
        throw CoilNotProcessedException("In Add Margin to Section: Section description empty, wind coil first");
    }
    if (margins.size() != 2) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Margin vector must have two elements");
    }
    auto sections = get_sections_description().value();
    auto globalSectionIndex = convert_conduction_section_index_to_global(sectionIndex);
    // _marginsPerSection is only sized through set_interleaving_level() / wind_by_sections().
    // When a Coil is reconstructed from JSON (e.g. via the Coil(json, bool) constructor used by
    // the PyOpenMagnetics bindings and any caller that round-trips sectionsDescription through JSON),
    // _marginsPerSection is empty even though sectionsDescription is populated, which made
    // the indexed assignment below segfault. Grow the vector to match the existing sections
    // and seed any uninitialized entries from the section's own margin, falling back to {0, 0}.
    if (_marginsPerSection.size() < sections.size()) {
        size_t previousSize = _marginsPerSection.size();
        _marginsPerSection.resize(sections.size(), {0, 0});
        for (size_t i = previousSize; i < sections.size(); ++i) {
            auto existingMargin = sections[i].get_margin();
            if (existingMargin) {
                if (std::holds_alternative<std::vector<double>>(existingMargin.value())) {
                    _marginsPerSection[i] = std::get<std::vector<double>>(existingMargin.value());
                }
            }
        }
    }
    _marginsPerSection[globalSectionIndex] = margins;
    sections[globalSectionIndex].set_margin(margins);

    set_sections_description(sections);


    bool windEvenIfNotFit = settings.get_coil_wind_even_if_not_fit();
    bool delimitAndCompact = settings.get_coil_delimit_and_compact();
    bool tryRewind = settings.get_coil_try_rewind();

    wind_by_sections();
    wind_by_layers();
    if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
        wind_by_turns();
        if (delimitAndCompact) {
            delimit_and_compact();
        }
    }
    if (tryRewind && !are_sections_and_layers_fitting()) {
        try_rewind();
    }
}

std::vector<Section> Coil::get_sections_description_conduction() const {
    std::vector<Section> sectionsConduction;
    if (!get_sections_description()) {
        throw CoilNotProcessedException("Not wound by sections");
    }
    std::vector<Section> sections = get_sections_description().value();
    for (auto section : sections) {
        if (section.get_type() == ElectricalType::CONDUCTION) {
            sectionsConduction.push_back(section);
        }
    }

    return sectionsConduction;
}

std::vector<Layer> Coil::get_layers_description_conduction() const {
    std::vector<Layer> layersConduction;
    if (!get_layers_description()) {
        throw CoilNotProcessedException("Not wound by layers");
    }
    std::vector<Layer> layers = get_layers_description().value();
    for (auto layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            layersConduction.push_back(layer);
        }
    }

    return layersConduction;
}

std::vector<Section> Coil::get_sections_description_insulation() const {
    std::vector<Section> sectionsInsulation;
    if (!get_sections_description()) {
        throw CoilNotProcessedException("Not wound by sections");
    }
    std::vector<Section> sections = get_sections_description().value();
    for (auto section : sections) {
        if (section.get_type() == ElectricalType::INSULATION) {
            sectionsInsulation.push_back(section);
        }
    }

    return sectionsInsulation;
}

std::vector<Layer> Coil::get_layers_description_insulation() const {
    std::vector<Layer> layersInsulation;
    if (!get_layers_description()) {
        throw CoilNotProcessedException("Not wound by layers");
    }
    std::vector<Layer> layers = get_layers_description().value();
    for (auto layer : layers) {
        if (layer.get_type() == ElectricalType::INSULATION) {
            layersInsulation.push_back(layer);
        }
    }

    return layersInsulation;
}

double Coil::calculate_external_proportion_for_wires_in_toroidal_cores(Core core, Coil coil) {
    CoreShape shape = std::get<CoreShape>(core.get_functional_description().get_shape());
    auto processedDescription = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});

    double coreWidth = processedDescription.get_width();

    if (!coil.get_turns_description()) {
        return 1;
    }

    auto turns = coil.get_turns_description().value();
    double maximumAdditionalRadialCoordinate = 0;
    for (size_t i = 0; i < turns.size(); ++i){
        if (turns[i].get_additional_coordinates()) {
            auto additionalCoordinates = turns[i].get_additional_coordinates().value();
            for (auto additionalCoordinate : additionalCoordinates){
                maximumAdditionalRadialCoordinate = std::max(maximumAdditionalRadialCoordinate, hypot(additionalCoordinate[0], additionalCoordinate[1]) + turns[i].get_dimensions().value()[0] / 2);
            }
        }
    }
    auto bobbin = coil.resolve_bobbin();

    auto sectionsOrientation = bobbin.get_winding_window_sections_orientation();

    if (maximumAdditionalRadialCoordinate > 0 && sectionsOrientation == WindingOrientation::OVERLAPPING) {
        auto sections = coil.get_sections_by_type(ElectricalType::INSULATION);
        for (auto section : sections){
            maximumAdditionalRadialCoordinate += section.get_dimensions()[0];
        }
    }

    if (maximumAdditionalRadialCoordinate == 0) {
        return 1;
    }

    return (2 * maximumAdditionalRadialCoordinate) / coreWidth;
}


double Coil::get_insulation_section_thickness(std::string sectionName) {
    return get_insulation_section_thickness(*this, sectionName);
}

double Coil::get_insulation_section_thickness(Coil coil, std::string sectionName) {
    if (!coil.get_sections_description()) {
        throw CoilNotProcessedException("Coil is missing sections description");
    }
    if (!coil.get_layers_description()) {
        throw CoilNotProcessedException("Coil is missing layers description");
    }

    auto layers = coil.get_layers_by_section(sectionName);

    double thickness = 0;

    for (auto layer : layers) {
        thickness += coil.get_insulation_layer_thickness(layer);
    }

    return thickness;
}

double Coil::get_insulation_layer_thickness(Coil coil, std::string layerName) {
    return coil.get_insulation_layer_thickness(layerName);
}

double Coil::get_insulation_layer_thickness(std::string layerName) {
    if (!get_layers_description()) {
        throw CoilNotProcessedException("Coil is missing layers description");
    }
    auto layer = get_layer_by_name(layerName);
    return get_insulation_layer_thickness(layer);
}

double Coil::get_insulation_layer_thickness(Layer layer) {
    if (!layer.get_coordinate_system()) {
        layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
    }
    if (layer.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        if (layer.get_orientation() == WindingOrientation::CONTIGUOUS) {
            return layer.get_dimensions()[1];
        }
        else {
            return layer.get_dimensions()[0];
        }
    }
    else {
        if (layer.get_orientation() == WindingOrientation::CONTIGUOUS) {
            auto bobbin = resolve_bobbin();
            auto bobbinProcessedDescription = bobbin.get_processed_description().value();
            auto windingWindows = bobbinProcessedDescription.get_winding_windows();

            double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
            double layerRadialHeight = layer.get_dimensions()[0];
            double radius = windingWindowRadialHeight - layerRadialHeight;
            double layerAngle = layer.get_dimensions()[1];
            double layerPerimeter = std::numbers::pi * (layerAngle / 180) * radius;
            return layerPerimeter;
        }
        else {
            return layer.get_dimensions()[0];
        }
    }
}

InsulationMaterial Coil::resolve_insulation_layer_insulation_material(std::string layerName){
    auto layer = get_layer_by_name(layerName);
    return resolve_insulation_layer_insulation_material(layer);
}

InsulationMaterial Coil::resolve_insulation_layer_insulation_material(Coil coil, std::string layerName) {
    auto layer = coil.get_layer_by_name(layerName);
    return coil.resolve_insulation_layer_insulation_material(layer);
}

InsulationMaterial Coil::resolve_insulation_layer_insulation_material(Layer layer) {
    if (!layer.get_insulation_material()) {
        layer.set_insulation_material(defaults.defaultLayerInsulationMaterial);
        // throw std::runtime_error("Layer is missing material information");
    }

    auto insulationMaterial = layer.get_insulation_material().value();
    // If the material is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(insulationMaterial)) {
        try {
            auto insulationMaterialData = find_insulation_material_by_name(std::get<std::string>(insulationMaterial));
            return insulationMaterialData;
        }
        catch (...) {
            // A NAMED material that fails to resolve is bad data; keep the
            // default so rendering/insulation flows continue, but make it visible
            logEntry("Insulation material '" + std::get<std::string>(insulationMaterial) + "' for layer " + layer.get_name() + " not found, falling back to " + defaults.defaultInsulationMaterial, "Coil", 1);
            return find_insulation_material_by_name(defaults.defaultInsulationMaterial);
        }

    }
    else {
        return InsulationMaterial(std::get<MAS::InsulationMaterial>(insulationMaterial));
    }
}

double Coil::get_insulation_layer_relative_permittivity(std::string layerName) {
    auto layer = get_layer_by_name(layerName);
    return get_insulation_layer_relative_permittivity(layer);
}
double Coil::get_insulation_layer_relative_permittivity(Coil coil, std::string layerName) {
    return coil.get_insulation_layer_relative_permittivity(layerName);
}
double Coil::get_insulation_layer_relative_permittivity(Layer layer) {
    auto coatingInsulationMaterial = resolve_insulation_layer_insulation_material(layer);
    if (!coatingInsulationMaterial.get_relative_permittivity())
        throw InvalidInputException(ErrorCode::INVALID_INSULATION_DATA, "Coating insulation material is missing dielectric constant");
    return coatingInsulationMaterial.get_relative_permittivity().value();
}

double Coil::get_insulation_section_relative_permittivity(std::string sectionName) {
    return get_insulation_section_relative_permittivity(*this, sectionName);
}
double Coil::get_insulation_section_relative_permittivity(Coil coil, std::string sectionName) {
    auto layers = coil.get_layers_by_section(sectionName);
    if (layers.size() == 0)
        throw CoilNotProcessedException("No layers in this section");

    double averagerelativePermittivity = 0;
    for (auto layer : layers) {
        averagerelativePermittivity += coil.get_insulation_layer_relative_permittivity(layer);
    }
    return averagerelativePermittivity / layers.size();
}

std::vector<double> Coil::get_turns_ratios() const {
    std::vector<double>  turnsRatios;
    for (size_t windingIndex = 1; windingIndex < get_functional_description().size(); ++windingIndex) {
        turnsRatios.push_back(double(get_functional_description()[0].get_number_turns()) / get_functional_description()[windingIndex].get_number_turns());
    }

    return turnsRatios;
}

std::vector<double> Coil::get_maximum_dimensions() {
    std::vector<double> bobbinMaximumDimensions = resolve_bobbin().get_maximum_dimensions();

    if (!get_turns_description()) {
        throw CoilNotProcessedException("Missing turns");
    }
    auto turns = get_turns_description().value();

    double width = 0;
    double height = 0;
    double depth = 0;

    for (auto turn : turns) {
        double turnMaxWidthPosition = 0;
        double turnMaxHeightPosition = 0;
        if (turn.get_additional_coordinates()) {
            auto additionalCoordinates = turn.get_additional_coordinates().value();
            turnMaxWidthPosition = fabs(additionalCoordinates[0][0]) + turn.get_dimensions().value()[0] / 2;
            turnMaxHeightPosition = fabs(additionalCoordinates[0][1]) + turn.get_dimensions().value()[1] / 2;
        }
        else {
            turnMaxWidthPosition = fabs(turn.get_coordinates()[0]) + turn.get_dimensions().value()[0] / 2;
            turnMaxHeightPosition = fabs(turn.get_coordinates()[1]) + turn.get_dimensions().value()[1] / 2;
        }

        width = std::max(width, turnMaxWidthPosition);
        height = std::max(height, turnMaxHeightPosition);
    }

    double bobbinExtraDepthDimension = bobbinMaximumDimensions[0] - bobbinMaximumDimensions[2];
    depth = width + bobbinExtraDepthDimension;

    width = std::max(width, bobbinMaximumDimensions[0]);
    height = std::max(height, bobbinMaximumDimensions[1]);
    depth = std::max(depth, bobbinMaximumDimensions[2]);

    return std::vector<double>{width, height, depth};
}


/**
 * @brief Generate winding section patterns for coil construction.
 *
 * For multi-secondary transformers, explores different orderings:
 * - Primary-Secondary ordering variations (P-S1-S2 vs P-S2-S1)
 * - Sandwich configurations for reduced leakage (S1-P-S2)
 *
 * Patterns are limited by defaults.maximumCoilPattern to bound runtime.
 */
std::vector<std::vector<size_t>> Coil::get_patterns(Inputs& inputs, CoreType coreType) {
    auto isolationSidesRequired = inputs.get_isolation_sides_used();

    if (!inputs.get_design_requirements().get_isolation_sides()) {
        throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, "Missing isolation sides requirement");
    }

    auto isolationSidesRequirement = inputs.get_design_requirements().get_isolation_sides().value();
    size_t numWindings = inputs.get_mutable_design_requirements().get_turns_ratios().size() + 1;

    std::vector<std::vector<size_t>> sectionPatterns;
    
    // Generate patterns based on isolation side permutations. Bound = n!/2 but at least 1:
    // the old float form tgamma(n+1)/2 evaluated to 0.5 for a single isolation side and the
    // `i < 0.5` comparison still ran one iteration; integer n!/2 would truncate to 0 and
    // produce NO pattern at all for plain inductors. (FIX L-COIL-3 factorial helper.)
    size_t patternBound = std::max<size_t>(1, factorial(isolationSidesRequired.size()) / 2);
    for(size_t i = 0; i < patternBound; ++i) {
        std::vector<size_t> sectionPattern;
        for (auto isolationSide : isolationSidesRequired) {
            for (size_t windingIndex = 0; windingIndex < numWindings; ++windingIndex) {
                if (isolationSidesRequirement[windingIndex] == isolationSide) {
                    sectionPattern.push_back(windingIndex);
                }
            }
        }
        sectionPatterns.push_back(sectionPattern);
        if (sectionPatterns.size() >= defaults.maximumCoilPattern) {
            break;
        }

        std::next_permutation(isolationSidesRequired.begin(), isolationSidesRequired.end());
    }
    
    // For multi-secondary (3+ windings), add secondary ordering variations
    // This helps find better configurations for leakage inductance and coupling
    if (numWindings >= 3 && sectionPatterns.size() < defaults.maximumCoilPattern) {
        // Find windings on the secondary isolation side
        std::vector<size_t> secondaryWindings;
        for (size_t windingIndex = 1; windingIndex < numWindings; ++windingIndex) {
            if (isolationSidesRequirement[windingIndex] != isolationSidesRequirement[0]) {
                secondaryWindings.push_back(windingIndex);
            }
        }
        
        // If we have multiple secondaries on the same side, try reversed order
        if (secondaryWindings.size() >= 2 && sectionPatterns.size() < defaults.maximumCoilPattern) {
            // Take the first pattern and reverse the secondary ordering
            auto basePattern = sectionPatterns[0];
            std::vector<size_t> reversedPattern;
            std::vector<size_t> secondariesInPattern;
            
            // Extract secondaries positions
            for (auto windingIndex : basePattern) {
                if (std::find(secondaryWindings.begin(), secondaryWindings.end(), windingIndex) != secondaryWindings.end()) {
                    secondariesInPattern.push_back(windingIndex);
                }
            }
            
            // Reverse and rebuild pattern
            if (secondariesInPattern.size() >= 2) {
                std::reverse(secondariesInPattern.begin(), secondariesInPattern.end());
                size_t secIdx = 0;
                for (auto windingIndex : basePattern) {
                    if (std::find(secondaryWindings.begin(), secondaryWindings.end(), windingIndex) != secondaryWindings.end()) {
                        reversedPattern.push_back(secondariesInPattern[secIdx++]);
                    } else {
                        reversedPattern.push_back(windingIndex);
                    }
                }
                
                // Add if it's a new pattern
                bool isNew = true;
                for (const auto& existing : sectionPatterns) {
                    if (existing == reversedPattern) {
                        isNew = false;
                        break;
                    }
                }
                if (isNew) {
                    sectionPatterns.push_back(reversedPattern);
                }
            }
        }
        
        // Add sandwich pattern: S1-P-S2 (primary in middle, reduces leakage)
        // Only for 3 windings with P on one side and S1,S2 on the other
        if (numWindings == 3 && secondaryWindings.size() == 2 && 
            sectionPatterns.size() < defaults.maximumCoilPattern) {
            std::vector<size_t> sandwichPattern = {secondaryWindings[0], 0, secondaryWindings[1]};
            
            bool isNew = true;
            for (const auto& existing : sectionPatterns) {
                if (existing == sandwichPattern) {
                    isNew = false;
                    break;
                }
            }
            if (isNew) {
                sectionPatterns.push_back(sandwichPattern);
            }
        }
    }

    if (coreType == CoreType::TOROIDAL) {
        // We remove the last combination as in toroids they go around
        size_t elementsToKeep = std::max(size_t(1), isolationSidesRequired.size() - 1);
        sectionPatterns = std::vector<std::vector<size_t>>(sectionPatterns.begin(), sectionPatterns.end() - (sectionPatterns.size() - elementsToKeep));
    }

    return sectionPatterns;
}


/**
 * @brief Get valid winding repetition patterns for a given configuration.
 *
 * For Common Mode Chokes (CMC) on toroidal cores, returns {2, 1} to enable
 * bifilar (interleaved) winding which is essential for common-mode rejection.
 * Bifilar winding ensures both windings have identical impedance characteristics.
 *
 * @param inputs Design inputs including sub-application type
 * @param coreType Type of core being wound
 * @return Vector of valid repetition counts to try
 */
std::vector<size_t> Coil::get_repetitions(Inputs& inputs, CoreType coreType) {
    // DMC topology (any winding count, any core type): one section per phase,
    // no bifilar interleave. Single-winding DMCs are plain inductors and
    // also take repetitions={1} (matches turnsRatios.size()==0 fallback
    // below, but make the intent explicit).
    if (inputs.get_design_requirements().get_topology() &&
        inputs.get_design_requirements().get_topology().value() == MAS::Topology::DIFFERENTIAL_MODE_CHOKE) {
        return {1};
    }
    // CMCs on toroids need bifilar winding for common-mode rejection
    if (coreType == CoreType::TOROIDAL) {
        if (inputs.get_design_requirements().get_sub_application() &&
            inputs.get_design_requirements().get_sub_application().value() == SubApplication::COMMON_MODE_NOISE_FILTERING) {
            // Bifilar (interleaved) winding is preferred for CMCs to ensure matched impedance
            return {2, 1};
        }
        // Non-CMC toroids (inductors) don't need interleaving
        if (inputs.get_design_requirements().get_turns_ratios().size() == 0) {
            return {1};
        }
        // Multi-winding toroidal transformers: prefer non-interleaved first
        // even when a leakage_inductance spec is set. Toroids have an
        // intrinsically closed magnetic path so the leakage gain from
        // interleaving is small, but rep=2 splits the inner arc into more
        // sections — on tight inner columns most rep=2 wind() attempts fail
        // geometrically and consume the wind-failure budget for nothing.
        // The {2} pass still runs after, so designs that genuinely need
        // interleaving for leakage are not blocked.
        return {1, 2};
    }
    
    if (inputs.get_design_requirements().get_turns_ratios().size() == 0) {
        return {1};
    }
    if (inputs.get_design_requirements().get_wiring_technology()) {
        if (inputs.get_design_requirements().get_wiring_technology().value() == WiringTechnology::PRINTED) {
            std::vector<size_t> repetitions;
            for (size_t repetition = 1; repetition <= (settings.get_coil_maximum_layers_planar() / (inputs.get_design_requirements().get_turns_ratios().size() + 1)); ++repetition) {
                repetitions.push_back(repetition);
            }
            return repetitions;
        }
    }
    if (inputs.get_design_requirements().get_leakage_inductance()) {
        return {2, 1};
    }
    else{
        return {1, 2};
    }
}

std::pair<std::vector<size_t>, size_t> Coil::check_pattern_and_repetitions_integrity(std::vector<size_t> pattern, size_t repetitions) {
    bool needsMerge = false;
    for (auto winding : get_functional_description()) {
        // TODO expand for more than one winding per layer
        size_t numberPhysicalTurns = winding.get_number_turns() * winding.get_number_parallels();
        if (numberPhysicalTurns < repetitions) {
            needsMerge = true;
        }
    }

    std::vector<size_t> newPattern;
    if (needsMerge) {
        for (size_t repetition = 1; repetition <= repetitions; ++repetition) {
            for (auto windingIndex : pattern) {
                auto winding = get_functional_description()[windingIndex];
                size_t numberPhysicalTurns = winding.get_number_turns() * winding.get_number_parallels();
                if (numberPhysicalTurns >= repetition) {
                    newPattern.push_back(windingIndex);
                }
            }
        }
        return {newPattern, 1};
    }
    return {pattern, repetitions};
}

bool Coil::is_edge_wound_coil() {
    auto wires = get_wires();
    for (auto wire : wires) {
        if (wire.get_type() != WireType::RECTANGULAR) {
            return false;
        }
    }

    return true;
}

void Coil::set_interlayer_insulation(double layerThickness, std::optional<std::string> material, std::optional<std::string> windingName, bool autowind) {
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    auto layersOrientation = _layersOrientation;

    // TODO: Properly think about insulation layers with weird windings
    // auto windingOrientation = get_winding_orientation();

    Layer layer;
    layer.set_partial_windings(std::vector<PartialWinding>{});
    layer.set_type(ElectricalType::INSULATION);
    layer.set_name("custom thickness temp");
    layer.set_orientation(layersOrientation);
    layer.set_turns_alignment(CoilAlignment::SPREAD); // HARDCODED, maybe in the future configure for shields made of turns?

    if (material) {
        layer.set_insulation_material(material.value());
    }
    else {
        layer.set_insulation_material(defaults.defaultLayerInsulationMaterial);
    }
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
        double windingWindowHeight = windingWindows[0].get_height().value();
        double windingWindowWidth = windingWindows[0].get_width().value();
        if (layersOrientation == WindingOrientation::OVERLAPPING) {
            layer.set_dimensions(std::vector<double>{layerThickness, windingWindowHeight});
        }
        else if (layersOrientation == WindingOrientation::CONTIGUOUS) {
            layer.set_dimensions(std::vector<double>{windingWindowWidth, layerThickness});
        }
    }
    else {
        layer.set_coordinate_system(CoordinateSystem::POLAR);
        double windingWindowAngle = windingWindows[0].get_angle().value();
        layer.set_dimensions(std::vector<double>{layerThickness, windingWindowAngle});
    }
    layer.set_filling_factor(1);

    if (windingName) {
        auto windingIndex = get_winding_index_by_name(windingName.value());
        _insulationInterLayers[windingIndex] = layer;
    }
    else {
        for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
            _insulationInterLayers[windingIndex] = layer;
        }
    }

    if (autowind) {
        wind();
    }

}

void Coil::set_intersection_insulation(double layerThickness, size_t numberInsulationLayers, std::optional<std::string> material, std::optional<std::pair<std::string, std::string>> windingNames, bool autowind) {
    auto bobbin = resolve_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();

    auto layersOrientation = _layersOrientation;

    // TODO: Properly think about insulation layers with weird windings
    auto windingOrientation = get_winding_orientation();

    if (windingOrientation == WindingOrientation::CONTIGUOUS && _layersOrientation == WindingOrientation::OVERLAPPING) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::CONTIGUOUS;
        }
    }
    if (windingOrientation == WindingOrientation::OVERLAPPING && _layersOrientation == WindingOrientation::CONTIGUOUS) {
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            layersOrientation = WindingOrientation::OVERLAPPING;
        }
    }

    std::vector<Layer> insulationLayers;
    Layer layer;
    layer.set_partial_windings(std::vector<PartialWinding>{});
    layer.set_type(ElectricalType::INSULATION);
    layer.set_name("custom thickness temp");
    layer.set_orientation(layersOrientation);
    layer.set_turns_alignment(CoilAlignment::SPREAD); // HARDCODED, maybe in the future configure for shields made of turns?

    if (material) {
        layer.set_insulation_material(material.value());
    }
    else {
        layer.set_insulation_material(defaults.defaultLayerInsulationMaterial);
    }
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        layer.set_coordinate_system(CoordinateSystem::CARTESIAN);
        double windingWindowHeight = windingWindows[0].get_height().value();
        double windingWindowWidth = windingWindows[0].get_width().value();
        if (layersOrientation == WindingOrientation::OVERLAPPING) {
            layer.set_dimensions(std::vector<double>{layerThickness, windingWindowHeight});
        }
        else if (layersOrientation == WindingOrientation::CONTIGUOUS) {
            layer.set_dimensions(std::vector<double>{windingWindowWidth, layerThickness});
        }
    }
    else {
        layer.set_coordinate_system(CoordinateSystem::POLAR);
        double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
        double windingWindowAngle = windingWindows[0].get_angle().value();
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            layer.set_dimensions(std::vector<double>{layerThickness, windingWindowAngle});
        }
        else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            double tapeThicknessInAngle = wound_distance_to_angle(layerThickness, windingWindowRadialHeight);
            layer.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
        }
    }
    layer.set_filling_factor(1);
    for (size_t layerIndex = 0; layerIndex < numberInsulationLayers; ++layerIndex) {
        insulationLayers.push_back(layer);
    }

    Section section;
    section.set_name("custom thickness temp");
    section.set_partial_windings(std::vector<PartialWinding>{});
    section.set_layers_orientation(layersOrientation);
    section.set_type(ElectricalType::INSULATION);

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        section.set_coordinate_system(CoordinateSystem::CARTESIAN);
        double windingWindowHeight = windingWindows[0].get_height().value();
        double windingWindowWidth = windingWindows[0].get_width().value();
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            section.set_dimensions(std::vector<double>{layerThickness * numberInsulationLayers, windingWindowHeight});
        }
        else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            section.set_dimensions(std::vector<double>{windingWindowWidth, layerThickness * numberInsulationLayers});
        }
    }
    else {
        section.set_coordinate_system(CoordinateSystem::POLAR);
        double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
        double windingWindowAngle = windingWindows[0].get_angle().value();
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            section.set_dimensions(std::vector<double>{layerThickness * numberInsulationLayers, windingWindowAngle});
        }
        else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            double tapeThicknessInAngle = wound_distance_to_angle(layerThickness * numberInsulationLayers, windingWindowRadialHeight);
            section.set_dimensions(std::vector<double>{windingWindowRadialHeight, tapeThicknessInAngle});
        }
    }
    section.set_filling_factor(1);

    if (windingNames) {
        auto windingIndex = get_winding_index_by_name(windingNames->first);
        auto nextWindingIndex = get_winding_index_by_name(windingNames->second);
        {
            auto windingsMapKey = std::pair<size_t, size_t>{windingIndex, nextWindingIndex};
            _insulationInterSectionsLayers[windingsMapKey] = insulationLayers;
            _insulationSections[windingsMapKey] = section;
        }
        {
            auto windingsMapKey = std::pair<size_t, size_t>{nextWindingIndex, windingIndex};
            _insulationInterSectionsLayers[windingsMapKey] = insulationLayers;
            _insulationSections[windingsMapKey] = section;
        }
    }
    else {
        for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
            for (size_t nextWindingIndex = 0; nextWindingIndex < get_functional_description().size(); ++nextWindingIndex) {
                auto windingsMapKey = std::pair<size_t, size_t>{nextWindingIndex, windingIndex};
                _insulationInterSectionsLayers[windingsMapKey] = insulationLayers;
                _insulationSections[windingsMapKey] = section;
            }
        }
    }

    if (autowind) {
        wind();
    }
}

std::vector<Wire> Coil::guess_round_wire_from_dc_resistance(std::vector<double> dcResistances, double maxError) {
    // RAII (ABT #113 sweep): the loop below winds/simulates repeatedly and can
    // throw; the manual restore at the end leaked wind_even_if_not_fit=true.
    SettingsGuard<bool> windEvenIfNotFitGuard(settings, &Settings::get_coil_wind_even_if_not_fit, &Settings::set_coil_wind_even_if_not_fit, true);

    double maximumError = DBL_MAX;
    size_t timeout = 100;
    while (maximumError > maxError) {
        auto wireLengths = get_wires_length();
        maximumError = 0;
        bool areWiresTheSame = true;
        auto calculatedDcResistances = WindingOhmicLosses::calculate_dc_resistance_per_winding(*this, defaults.ambientTemperature);
        for (size_t index = 0; index < get_functional_description().size(); ++index) {
            double error = fabs(calculatedDcResistances[index] - dcResistances[index]) / dcResistances[index];
            if (error > maxError) {
                auto dcResistancesPerMeter = dcResistances[index] / wireLengths[index];
                auto newWire = Wire::get_wire_for_dc_resistance_per_meter(dcResistancesPerMeter);
                if (newWire.get_name().value() != get_mutable_functional_description()[index].resolve_wire().get_name().value()) {
                    areWiresTheSame = false;
                }
                get_mutable_functional_description()[index].set_wire(newWire);
            }
            maximumError = std::max(maximumError, error);
        }
        if (areWiresTheSame) {
            break;
        }
        set_turns_description(std::nullopt);
        unwind();
        fast_wind();
        timeout--;
        if (timeout == 0) {
            break;
        }
    }

    return get_wires();
}

std::vector<double> Coil::resolve_margin(size_t sectionIndex) {
    if (!get_sections_description()) {
        throw CoilNotProcessedException("Sections not found");
    }
    auto sections = get_sections_description().value();
    return resolve_margin(sections[sectionIndex]);
}

std::vector<double> Coil::resolve_margin(Section section) {
    if (!section.get_margin()) {
        return {0.0, 0.0};
    }
    return resolve_margin(section.get_margin().value());
}

std::vector<double> Coil::resolve_margin(Margin marginVariant) {
    if (std::holds_alternative<std::vector<double>>(marginVariant)) {
        auto margin = std::get<std::vector<double>>(marginVariant);
        return margin;
    }
    else {
        std::vector<double> margin;
        auto marginInfo = std::get<MarginInfo>(marginVariant);
        margin.push_back(marginInfo.get_top_or_left_width());
        margin.push_back(marginInfo.get_bottom_or_right_width());
        return margin;
    }
}

MarginInfo Coil::resolve_margin_info(size_t sectionIndex) {
    if (!get_sections_description()) {
        throw CoilNotProcessedException("Sections not found");
    }
    auto sections = get_sections_description().value();
    return resolve_margin_info(sections[sectionIndex]);
}

MarginInfo Coil::resolve_margin_info(Section section) {
    if (!section.get_margin()) {
        MarginInfo marginInfo;
        marginInfo.set_top_or_left_width(0);
        marginInfo.set_bottom_or_right_width(0);
        marginInfo.set_number_layers(0);
        section.set_margin(marginInfo);
        return marginInfo;
    }
    return resolve_margin_info(section.get_margin().value());
}

MarginInfo Coil::resolve_margin_info(Margin marginVariant) {
    if (std::holds_alternative<std::vector<double>>(marginVariant)) {
        MarginInfo marginInfo;
        auto margin = std::get<std::vector<double>>(marginVariant);
        marginInfo.set_top_or_left_width(margin[0]);
        marginInfo.set_bottom_or_right_width(margin[1]);
        return marginInfo;
    }
    else {
        std::vector<double> margin;
        auto marginInfo = std::get<MarginInfo>(marginVariant);
        return marginInfo;
    }
}

void Winding::set_isolation_side_from_index(size_t windingIndex) {
    set_isolation_side(get_isolation_side_from_index(windingIndex));
}

Coil Coil::create_quick_coil(std::string coreShapeName, std::vector<int64_t> numberTurns, std::vector<int64_t> numberParallels, std::vector<OpenMagnetics::Wire> wires, WindingOrientation windingOrientation, WindingOrientation layersOrientation, CoilAlignment turnsAlignment, CoilAlignment sectionsAlignment, uint8_t interleavingLevel, bool useBobbin, int numberStacks) {
    Coil coil;

    auto core = Core::create_quick_core(coreShapeName, "Dummy");
    OpenMagnetics::Bobbin bobbin;
    if (core.get_shape_family() == CoreShapeFamily::T) {
        bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);
    }
    else {
        bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, useBobbin);
    }
    coil.set_bobbin(bobbin);

    for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex){
        OpenMagnetics::Winding winding;
        winding.set_name("winding " + std::to_string(windingIndex));
        winding.set_number_turns(numberTurns[windingIndex]);
        winding.set_number_parallels(numberParallels[windingIndex]);
        winding.set_isolation_side_from_index(windingIndex);
        if (windingIndex < wires.size()) {
            winding.set_wire(wires[windingIndex]);
        }
        else {
            winding.set_wire("Round 0.475 - Grade 1");
        }
        coil.get_mutable_functional_description().push_back(winding);
    }
    coil.set_interleaving_level(interleavingLevel);
    coil.set_winding_orientation(windingOrientation);
    coil.set_layers_orientation(layersOrientation);
    coil.set_turns_alignment(turnsAlignment);
    coil.set_section_alignment(sectionsAlignment);

    coil.wind();
    return coil;
}

std::vector<Turn> Coil::get_turns_touching_bobbin_column() {
    if (!get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description");
    }
    auto turns = get_turns_description().value();
    return get_turns_touching_bobbin_column(turns);
}

std::vector<Turn> Coil::get_turns_touching_bobbin_column(std::vector<size_t> turnIndexes) {
    if (!get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description");
    }
    auto turns = get_turns_description().value();
    std::vector<Turn> filteredTurns;
    for (auto turnIndex : turnIndexes) {
        filteredTurns.push_back(turns[turnIndex]);
    }
    return get_turns_touching_bobbin_column(filteredTurns);
}

std::vector<Turn> Coil::get_turns_touching_bobbin_column(std::vector<Turn> turns) {
    std::vector<Turn> touchingTurns;
    auto bobbin = resolve_bobbin();
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions();
    auto bobbinLeftCoordinate = windingWindowCoordinates[0] - windingWindowDimensions[0] / 2;
    auto windingWindowShape = bobbin.get_winding_window_shape();
    if (windingWindowShape == WindingWindowShape::RECTANGULAR) {
        for (auto turn : turns) {
            auto leftSideCoordinate = turn.get_coordinates()[0] - turn.get_dimensions().value()[0] / 2;
            if (fabs((leftSideCoordinate - bobbinLeftCoordinate) / bobbinLeftCoordinate) < 0.05) {
                touchingTurns.push_back(turn);
            }
        }
    }
    else {
        throw NotImplementedException("Not implemented yet");

    }
    return touchingTurns;
}

std::vector<Turn> Coil::get_turns_touching_bobbin_walls() {
    if (!get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description");
    }
    auto turns = get_turns_description().value();
    return get_turns_touching_bobbin_walls(turns);
}

std::vector<Turn> Coil::get_turns_touching_bobbin_walls(std::vector<size_t> turnIndexes) {
    if (!get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description");
    }
    auto turns = get_turns_description().value();
    std::vector<Turn> filteredTurns;
    for (auto turnIndex : turnIndexes) {
        filteredTurns.push_back(turns[turnIndex]);
    }
    return get_turns_touching_bobbin_walls(filteredTurns);
}

std::vector<Turn> Coil::get_turns_touching_bobbin_walls(std::vector<Turn> turns) {
    std::vector<Turn> touchingTurns;
    auto bobbin = resolve_bobbin();
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates();
    auto windingWindowDimensions = bobbin.get_winding_window_dimensions();
    auto bobbinTopCoordinate = windingWindowCoordinates[1] + windingWindowDimensions[1] / 2;
    auto bobbinBottomCoordinate = windingWindowCoordinates[1] - windingWindowDimensions[1] / 2;
    auto windingWindowShape = bobbin.get_winding_window_shape();
    if (windingWindowShape != WindingWindowShape::RECTANGULAR) {
        return touchingTurns;
    }
    for (auto turn : turns) {
        auto topSideCoordinate = turn.get_coordinates()[1] + turn.get_dimensions().value()[1] / 2;
        auto bottomSideCoordinate = turn.get_coordinates()[1] - turn.get_dimensions().value()[1] / 2;
        if (fabs((bottomSideCoordinate - bobbinBottomCoordinate) / bobbinBottomCoordinate) < 0.05) {
            touchingTurns.push_back(turn);
        }
        if (fabs((topSideCoordinate - bobbinTopCoordinate) / bobbinTopCoordinate) < 0.05) {
            touchingTurns.push_back(turn);
        }
    }
    return touchingTurns;
}

} // namespace OpenMagnetics
 