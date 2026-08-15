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

// Number of physical turns wound TOGETHER as one N-filar bundle on this layer (ABT #578).
// WIND_BY_CONSECUTIVE_PARALLELS lays the parallels of one electrical turn side by side
// (P0T0, P1T0, ... P0T1, P1T1, ...), so a bundle is one turn of every parallel actually PRESENT
// on the layer — parallels whose proportion is zero are not wound here and must not be counted,
// or the bundle would be padded with turns that never get placed. Every other winding style
// advances one parallel at a time, so each turn is its own bundle.
static int64_t get_layer_bundle_size(const Layer& layer) {
    if (!layer.get_winding_style() ||
        layer.get_winding_style().value() != WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS) {
        return 1;
    }
    int64_t bundleSize = 0;
    auto parallelsProportion = layer.get_partial_windings()[0].get_parallels_proportion();
    for (auto proportion : parallelsProportion) {
        if (roundFloat(proportion, 10) > 0) {
            bundleSize++;
        }
    }
    // A layer with no wound parallel places no turns at all, so the bundle size is moot; 1 keeps
    // the station arithmetic below total instead of dividing by zero.
    return bundleSize > 0 ? bundleSize : 1;
}

// Turn-centre stations for a SPREAD layer along its turn axis, in ASCENDING coordinate order.
//
// SPREAD distributes the layer's slack FENCE-POST (ABT #579): the outermost turns' SURFACES sit
// on the layer edges and the slack lands only in the gaps BETWEEN turn groups — never as a dead
// half-gap margin outside them, which is exactly what centring each turn in an equal
// axisLength/N slot used to leave at both flanges. When turns are wound together N-filar the
// members of a bundle stay TOUCHING and the slack is shared over the gaps BETWEEN bundles only
// (ABT #578), so every terminal connection leaves straight out along its own bundle's row
// instead of threading between parallels that have been pulled apart:
//
//     gap = (axisLength - numberPhysicalTurns * wireSize) / (numberBundles - 1)
//
// A single bundle — or a layer with no slack left to give — falls back to a centred contiguous
// block: there is no gap to distribute, and centring keeps an over-full layer overflowing
// symmetrically instead of hanging off one edge. Both cases are the same expression: with the
// gap above, the block length is exactly axisLength, so centring it lands the first turn's
// surface on the low edge.
static std::vector<double> compute_spread_turn_stations(double axisCenter,
                                                        double axisLength,
                                                        double wireSize,
                                                        int64_t numberPhysicalTurns,
                                                        int64_t bundleSize) {
    std::vector<double> stations;
    if (numberPhysicalTurns <= 0) {
        return stations;
    }
    stations.reserve(numberPhysicalTurns);
    if (bundleSize < 1) {
        bundleSize = 1;
    }
    int64_t numberBundles = (numberPhysicalTurns + bundleSize - 1) / bundleSize;
    double slack = axisLength - double(numberPhysicalTurns) * wireSize;
    double gap = (numberBundles > 1 && slack > 0) ? slack / double(numberBundles - 1) : 0;
    double blockLength = double(numberPhysicalTurns) * wireSize + double(numberBundles - 1) * gap;
    double blockLow = axisCenter - blockLength / 2;
    // Each station is computed from its own index rather than accumulated, so the run carries one
    // rounding instead of N of them. Accumulating drifted the last turn off the far edge by ~1 nm,
    // which is enough to make the margin-tape arithmetic see the turns as overflowing the section.
    // turnIndex / bundleSize is the number of COMPLETED bundles before this turn, i.e. how many
    // inter-bundle gaps precede it; turns inside a bundle share that count and so stay touching.
    for (int64_t turnIndex = 0; turnIndex < numberPhysicalTurns; ++turnIndex) {
        double position = blockLow + wireSize / 2 + double(turnIndex) * wireSize + double(turnIndex / bundleSize) * gap;
        stations.push_back(roundFloat(position, 9));
    }
    return stations;
}

// Move the running turn centre onto the next SPREAD station. Only the spread axis is driven from
// the stations; the other axis keeps the fixed value the alignment switch set. A no-op for every
// non-SPREAD alignment (no stations), whose uniform increment still carries placement.
// Running out of stations means the placer laid more physical turns than the layer reported, which
// would silently stack turns on the last station's coordinates — so say so instead.
static void take_next_spread_station(const std::vector<double>& stations,
                                     size_t& stationIndex,
                                     size_t stationAxis,
                                     const std::string& layerName,
                                     double& currentTurnCenterWidth,
                                     double& currentTurnCenterHeight) {
    if (stations.empty()) {
        return;
    }
    if (stationIndex >= stations.size()) {
        throw std::runtime_error("SPREAD turn stations exhausted on layer " + layerName + ": the placer laid more physical turns than the layer reports (" +
                                 std::to_string(stations.size()) + ")");
    }
    if (stationAxis == 1) {
        currentTurnCenterHeight = stations[stationIndex];
    }
    else {
        currentTurnCenterWidth = stations[stationIndex];
    }
    stationIndex++;
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
    if (hasSectionsData) {
        // Externally-provided descriptions always sit at their FINAL multi-window
        // positions (the +x winding frame is transient inside wind()); without
        // this a later delimit_and_compact would re-compact mirrored-window
        // sections as if they were frame-local.
        set_group_window_sides_applied(true);
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
    // Clear, don't seed: a level change invalidates any per-section margins (the winders
    // resize lazily), but seeding {0,0} entries here left the vector NON-empty, which
    // (a) defeated the ABT #676 margin-recovery gate `_marginsPerSection.empty()` in
    // wind() — set_interleaving_level-then-rewind silently dropped persisted margins —
    // and (b) made a subsequent preload_margins() append AFTER the seeded entries,
    // landing every preloaded pair at the wrong index.
    _marginsPerSection.clear();
}

void Coil::reset_margins_per_section() {
    _marginsPerSection.clear();
    // See _marginsExplicitlyCleared in Coil.h (ABT #724): without this, wind()'s ABT #676
    // recovery resurrected the persisted section margins and an explicit reset could never
    // actually clear them.
    _marginsExplicitlyCleared = true;
}

void Coil::reset_insulation() {
    _marginsPerSection.clear();
    _insulationSections.clear();
}

size_t Coil::get_interleaving_level() const {
    return _interleavingLevel;
}

size_t Coil::get_current_repetitions() const {
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

WindingOrientation Coil::get_layers_orientation(std::optional<std::string> sectionName) const {
    // Mirror get_turns_alignment: honour the per-section override map. Until 2026-08 the
    // per-section set_layers_orientation overload wrote _layersOrientationPerSection but
    // nothing ever read it, so the override (exposed in PyOM and the web winding studio)
    // was a silent no-op.
    if (sectionName) {
        auto it = _layersOrientationPerSection.find(sectionName.value());
        if (it != _layersOrientationPerSection.end()) {
            return it->second;
        }
    }
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

    // ABT #187: mean radius of each conduction ring (turns are cartesian here), used to find which
    // rings a radial terminal lead passes over so those rings get angular squeeze markers.
    std::map<std::string, double> ringRadiusByLayer;
    {
        std::map<std::string, std::pair<double, size_t>> radiusAccumulator;
        for (const auto& turn : turns) {
            if (!turn.get_layer()) {
                continue;
            }
            auto& acc = radiusAccumulator[turn.get_layer().value()];
            acc.first += std::hypot(turn.get_coordinates()[0], turn.get_coordinates()[1]);
            acc.second++;
        }
        for (const auto& [layerName, acc] : radiusAccumulator) {
            if (acc.second > 0) {
                ringRadiusByLayer[layerName] = acc.first / double(acc.second);
            }
        }
    }
    // Ring -> owning winding, for the entrance-corridor emission (ABT #723): the input
    // connection's corridor applies to its OWN winding's rings only.
    std::map<std::string, std::string> windingByRingName;
    for (const auto& layer : allLayers) {
        windingByRingName[layer.get_name()] = layer.get_partial_windings()[0].get_winding();
    }
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto windingName = coil.get_functional_description()[windingIndex].get_name();
        double wireOuterWidth = wires[windingIndex].get_maximum_outer_width();
        double wireOuterHeight = wires[windingIndex].get_maximum_outer_height();
        int64_t numberParallels = int64_t(coil.get_number_parallels(windingIndex));
        // Terminal leads route radially out past the outermost turn to the window border.
        //
        // ABT #230: the border is CAPPED at the bore wall. maxTurnRadius + 1.5 * wireOuterWidth
        // overshoots it whenever the outermost ring is wall-adjacent (maxTurnRadius ~ bore - wr),
        // which put the rect's far edge up to ~2 wire ODs INSIDE the core annulus — e.g. 12.96 mm
        // against a 12.0 mm bore on T 40/24/16 with 0.959 mm OD. MVB++ replays these rects verbatim
        // as 3D lead routes, so the overrun would place copper inside the core.
        double uncappedBorder = maxTurnRadius + 1.5 * wireOuterWidth;
        double boreRadius = uncappedBorder;  // no cap available -> keep the historical border
        {
            auto windingWindows = coil.resolve_bobbin().get_processed_description().value().get_winding_windows();
            if (!windingWindows.empty() && windingWindows[0].get_radial_height()) {
                boreRadius = windingWindows[0].get_radial_height().value();
            }
        }
        double radialBorder = std::min(uncappedBorder, boreRadius);

        auto addTerminalLead = [&](const Turn& connectingTurn, int64_t parallel, bool isEntrance) {
            auto c = connectingTurn.get_coordinates();
            double radius = std::hypot(c[0], c[1]);
            double angle = std::atan2(c[1], c[0]);

            // ABT #723 (owner ruling, Alf 2026-08-14): the INPUT internal terminal connection
            // owns its angular corridor through the WHOLE winding depth — its below-core
            // vertical rises at this azimuth past every deeper ring, and those rings'
            // below-core returns must clear it, so no turn may sit in that angle on ANY ring
            // of this winding. The outward loop below only marks rings the drawn radial run
            // crosses (and the wall-adjacent entrance crosses none); emit the corridor marker
            // on each ring INWARD of the entrance too. align_blocked_ring_turns then
            // displaces those rings' stations out of the corridor (or reports capacity
            // deficits for the blocking re-wind), and the final strict sweep finds clear
            // return azimuths. Emitted before the radial-run guard: the entrance lead often
            // has no drawable radial run at all.
            if (isEntrance && radius > 1e-9) {
                for (const auto& [ringName, ringRadius] : ringRadiusByLayer) {
                    if (connectingTurn.get_layer() && connectingTurn.get_layer().value() == ringName) {
                        continue;   // its own ring: the entrance turn IS the connection here
                    }
                    auto ringWindingIt = windingByRingName.find(ringName);
                    if (ringWindingIt == windingByRingName.end() || ringWindingIt->second != windingName) {
                        continue;   // the corridor lives in this winding's own sector
                    }
                    if (ringRadius >= radius - wireOuterWidth / 2) {
                        continue;   // not inward of the entrance turn
                    }
                    ConnectionReservedSpace corridor;
                    corridor.coordinateSystem = CoordinateSystem::POLAR;
                    corridor.isTerminal = true;
                    corridor.winding = windingName;
                    corridor.parallel = parallel;
                    corridor.section = connectingTurn.get_section().value_or("");
                    corridor.layer = ringName;
                    corridor.coordinates = {roundFloat(ringRadius * std::cos(angle), 9), roundFloat(ringRadius * std::sin(angle), 9)};
                    corridor.dimensions = {roundFloat(wireOuterWidth, 9), roundFloat(wireOuterHeight, 9)};
                    corridor.routedLength = 0;  // space-only: the vertical's copper is billed by the terminal lead itself
                    corridor.rotation = roundFloat(angle * 180.0 / std::numbers::pi, 6);
                    spaces.push_back(corridor);
                }
            }

            if (radius <= 1e-9 || radialBorder <= radius) {
                return;
            }
            // ABT #230: the NEAR edge is the wire envelope, not the crossing centreline. The
            // concentric path already spans [turnX - w/2, borderX + w/2]; taking the centreline here
            // reserved only half the wire at the connecting turn, so the two conventions disagreed.
            double radiusNear = radius - wireOuterWidth / 2;
            double radiusMid = (radiusNear + radialBorder) / 2;
            ConnectionReservedSpace lead;
            lead.coordinateSystem = CoordinateSystem::POLAR;
            lead.isTerminal = true;
            lead.winding = windingName;
            lead.parallel = parallel;
            lead.section = connectingTurn.get_section().value_or("");
            lead.layer = "";
            lead.coordinates = {roundFloat(radiusMid * std::cos(angle), 9), roundFloat(radiusMid * std::sin(angle), 9)};
            lead.dimensions = {roundFloat(radialBorder - radiusNear, 9), wireOuterHeight};
            lead.routedLength = roundFloat(radialBorder - radiusNear, 9);  // the radial run itself
            lead.rotation = roundFloat(angle * 180.0 / std::numbers::pi, 6);
            spaces.push_back(lead);

            // ABT #187: the radial run passes over every ring OUTWARD of the connecting turn. Emit an
            // angular squeeze marker on each crossed ring: centred at the crossing point, radial
            // extent = the ring band it crosses, azimuthal extent = the lead's wire thickness
            // (dimensions[1], same convention as the drawn lead). align_blocked_ring_turns turns
            // these into blocked angular corridors that the ring's turns must clear; 3D consumers
            // can read the same markers as MKF's routing decision. Rings at the crossed radius whose
            // sector does not reach the lead azimuth simply have no turns near the corridor, so the
            // marker is inert for them (align_blocked_ring_turns confines displacement to each
            // ring's own occupied arc).
            for (const auto& [ringName, ringRadius] : ringRadiusByLayer) {
                if (connectingTurn.get_layer() && connectingTurn.get_layer().value() == ringName) {
                    continue;  // its own ring: the lead starts here, no crossing
                }
                if (ringRadius <= radius + wireOuterWidth / 2 || ringRadius >= radialBorder) {
                    continue;  // not outward of the connecting turn
                }

                ConnectionReservedSpace crossing;
                crossing.coordinateSystem = CoordinateSystem::POLAR;
                crossing.isTerminal = true;
                crossing.winding = windingName;
                crossing.parallel = parallel;
                crossing.section = connectingTurn.get_section().value_or("");
                crossing.layer = ringName;
                crossing.coordinates = {roundFloat(ringRadius * std::cos(angle), 9), roundFloat(ringRadius * std::sin(angle), 9)};
                crossing.dimensions = {roundFloat(wireOuterWidth, 9), roundFloat(wireOuterHeight, 9)};
                crossing.routedLength = 0;  // space-only squeeze: the lead's copper is the drawn radial run
                crossing.rotation = roundFloat(angle * 180.0 / std::numbers::pi, 6);
                spaces.push_back(crossing);
            }
        };
        for (int64_t parallel = 0; parallel < numberParallels; ++parallel) {
            auto key = std::make_pair(windingName, parallel);
            if (entranceTurn.count(key)) {
                addTerminalLead(entranceTurn.at(key), parallel, true);
            }
            if (exitTurn.count(key)) {
                addTerminalLead(exitTurn.at(key), parallel, false);
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
                link.coordinateSystem = CoordinateSystem::POLAR;
                link.winding = windingName;
                link.parallel = parallel;
                link.section = windingLayers[i].get_section().value_or("");
                link.layer = "";
                link.coordinates = {roundFloat((a[0] + b[0]) / 2, 9), roundFloat((a[1] + b[1]) / 2, 9)};
                link.dimensions = {roundFloat(length, 9), wireOuterHeight};
                link.routedLength = roundFloat(length, 9);  // centre-to-centre inter-ring hop
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
    // A winding that does not FIT still has to be drawn honestly: when MKF places layers past
    // the window's outer edge (over-subscribed design), the terminal border and the crossed-layer
    // test must follow the COPPER, not the window. Clipping them at the window let every layer
    // outside it escape the lead's reservation, so an out-of-space design drew as if it fitted --
    // Alf, 2026-08-08: "when the turns go beyond the winding window (even if it's impossible in
    // reality) we have to maintain that restriction to show the customer how much out of space
    // he is".
    double outermostLayerX = windowOuterX;
    for (const auto& layer : allLayers)
        outermostLayerX = std::max(outermostLayerX,
                                   layer.get_coordinates()[0] + layer.get_dimensions()[0] / 2);
    windowOuterX = std::max(windowOuterX, outermostLayerX);
    // Axial extent of the window: terminal leads run along its top edge (entrance) or bottom edge
    // (exit), i.e. above/below all the (blocking-shrunk) conduction layers.
    double windowCenterY = windowCenterTurnAxis;
    double windowTopY = windowCenterY + windowSizeTurnAxis / 2;
    double windowBottomY = windowCenterY - windowSizeTurnAxis / 2;

    // ABT #229: per-(window, edge) row allocator. Every run that routes ALONG a window edge (a
    // terminal lead or a U interleaved continuation) is a separate physical conductor, so each gets
    // its own row, stacked inward from the edge in allocation order — the K parallels of an N-filar
    // group can no longer be drawn on the SAME line (which made 3D consumers' collision gates throw
    // on coincident centrelines). Rows are per window (multi-column groups run in the +x frame here,
    // so different windows must not share a stack) and per edge (top/bottom). The allocator is the
    // single source of truth: the drawn geometry AND the turn blocking both follow its decision.
    std::map<std::string, size_t> windowIndexBySection;
    // get_sections_description() returns the optional BY VALUE; iterating its .value()
    // directly binds the range-for to storage inside a temporary destroyed at the end of
    // the full expression (UB; gcc -Wdangling-reference). Materialize the copy first.
    auto sectionsDescription = get_sections_description();
    if (sectionsDescription) {
        for (const auto& section : sectionsDescription.value()) {
            windowIndexBySection[section.get_name()] = resolve_section_winding_window_index(section);
        }
    }
    auto windowIndexOf = [&](const std::string& sectionName) -> size_t {
        auto found = windowIndexBySection.find(sectionName);
        return found == windowIndexBySection.end() ? 0 : found->second;
    };
    // ABT #684: margin tape is reserved for TAPE. A terminal lead or an edge continuation runs
    // ALONG the window edge, so measuring its row from the window put the copper inside the
    // margin band — visible in the 2D as the magenta lead crossing the yellow tape. A run along
    // an edge crosses the whole window radially, so it must clear the LARGEST margin on that edge
    // in that window. margin[0]/margin[1] are the section's "top or left"/"bottom or right".
    // ABT #726: in the virtual frame "top" is the HIGH turn-axis side. For OVERLAPPING layers the
    // turn axis is the real y, so margin[0] ("top") insets the high side; for CONTIGUOUS layers
    // the frame is the x<->y transpose — the virtual high side is the real RIGHT, which margin[1]
    // ("bottom or right") owns. align_blocked_layer_turns insets with exactly this convention
    // (high -= margin[turnAxis == 1 ? 0 : 1]); mapping margin[0] to the virtual top
    // unconditionally re-created the #684 lead-through-tape defect on the contiguous path: rows
    // inset on the tape-free side, stacked through the tape on the other.
    std::map<size_t, double> topMarginPerWindow;
    std::map<size_t, double> bottomMarginPerWindow;
    if (sectionsDescription) {
        for (const auto& marginSection : sectionsDescription.value()) {
            if (marginSection.get_type() != ElectricalType::CONDUCTION) {
                continue;
            }
            auto sectionMargin = resolve_margin(marginSection);
            size_t windowIndex = windowIndexOf(marginSection.get_name());
            double virtualTopMargin = layersAreContiguous ? sectionMargin[1] : sectionMargin[0];
            double virtualBottomMargin = layersAreContiguous ? sectionMargin[0] : sectionMargin[1];
            topMarginPerWindow[windowIndex] = std::max(topMarginPerWindow[windowIndex], virtualTopMargin);
            bottomMarginPerWindow[windowIndex] = std::max(bottomMarginPerWindow[windowIndex], virtualBottomMargin);
        }
    }
    auto edgeBaseY = [&](size_t windowIndex, bool atTop) -> double {
        const auto& margins = atTop ? topMarginPerWindow : bottomMarginPerWindow;
        auto found = margins.find(windowIndex);
        double margin = found == margins.end() ? 0.0 : found->second;
        return atTop ? windowTopY - margin : windowBottomY + margin;
    };
    // ABT #615: edge rows are SHARED by runs whose RADIAL SPANS don't overlap (Alf, 2026-08-09:
    // the primary's inter-section run covers the secondary's section and vice versa — disjoint
    // spans, ONE row, "which can then be reused by the inter section connection in secondary").
    // Each (window, edge) keeps a stack of rows; a run takes the outermost row of its own wire
    // height whose occupied intervals it does not cross, else opens a new row underneath. Sharing
    // never permutes rows, so the #577 emission-order doctrine is untouched, and blocking stays
    // exact: a crossed layer only sees the runs whose span actually covers it.
    struct EdgeRow {
        double height;                                 // rows shared only between equal wire heights
        double depthBefore;                            // stack depth from the window edge to this row
        std::vector<std::pair<double, double>> spans;  // occupied radial intervals
    };
    std::map<std::pair<size_t, int>, std::vector<EdgeRow>> edgeRows;
    // ABT #615: the one shared inter-section continuation band per (window, edge) — see the
    // continuation allocation below.
    struct ContinuationBand {
        double edgeY;
        double runDepth;
        double height;
        size_t rowIndex;  // into edgeRows[{window, edge}], so reuse can register more spans
    };
    std::map<std::pair<size_t, int>, ContinuationBand> continuationBand;
    // Allocates a row for a run of height `wireHeight` spanning [spanLo, spanHi] radially; returns
    // {edgeY, runDepth} where edgeY is the row's centre and runDepth the distance from the window
    // edge to the row's inner side.
    auto allocateEdgeRow = [&](size_t windowIndex, bool atTop, double wireHeight,
                               double spanLo, double spanHi) -> std::pair<double, double> {
        auto& rows = edgeRows[{windowIndex, atTop ? 0 : 1}];
        EdgeRow* target = nullptr;
        for (auto& row : rows) {
            if (std::abs(row.height - wireHeight) > 1e-12) {
                continue;
            }
            bool overlaps = false;
            for (const auto& [lo, hi] : row.spans) {
                if (spanLo < hi - 1e-12 && lo < spanHi - 1e-12) {
                    overlaps = true;
                    break;
                }
            }
            if (!overlaps) {
                target = &row;
                break;
            }
        }
        if (target == nullptr) {
            double depthBefore = rows.empty() ? 0.0 : rows.back().depthBefore + rows.back().height;
            rows.push_back({wireHeight, depthBefore, {}});
            target = &rows.back();
        }
        target->spans.push_back({spanLo, spanHi});
        // Rows stack inward from the margin's inner face, not from the window edge.
        double rowBaseY = edgeBaseY(windowIndex, atTop);
        double edgeY = atTop ? roundFloat(rowBaseY - target->depthBefore - wireHeight / 2, 9)
                             : roundFloat(rowBaseY + target->depthBefore + wireHeight / 2, 9);
        return {edgeY, target->depthBefore + wireHeight};
    };

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

    // Terminal lead emissions are collected for ALL windings and sorted GLOBALLY before any row
    // is allocated (Alf, 2026-08-09, on 25_psps: sections were displaced because a lead whose span
    // never reaches them forced another lead's row deeper). Primary key: WIDEST-CROSSING FIRST —
    // the lead crossing the most layers takes the shallowest row, so the fewest layers pay the
    // deepest depths and a narrow lead (an exit connecting near the border) stacks below without
    // charging inner sections. Groups (same window/edge/crossed-set) then keep first-seen order,
    // and WITHIN a group leads emit nearest-edge-turn first (the ABT #577 rule, unchanged).
    struct LeadEmission {
        Turn turn;
        int64_t parallel;
        bool atTop;
        std::string crossedSignature;
        size_t crossedCount;
        size_t groupRank;      // first emission index of this lead's group, keeps groups in order
        double edgeDistance;   // ordering WITHIN a group: nearest the edge first
        std::string windingName;
        double wireW;
        double wireH;
    };
    std::vector<LeadEmission> allEmissions;
    auto crossedLayerCountAndSignature = [&](const Turn& connectingTurn) {
        std::string signature;
        size_t count = 0;
        double turnX = connectingTurn.get_coordinates()[0];
        for (const auto& crossed : allLayers) {
            double crossedX = crossed.get_coordinates()[0];
            if (crossedX > turnX + 1e-9 && crossedX < windowOuterX) {
                signature += crossed.get_name() + "|";
                ++count;
            }
        }
        return std::make_pair(count, signature);
    };
    auto addTerminalLead = [&](const std::string& windingName, double wireOuterWidth,
                           double wireOuterHeight, const Turn& connectingTurn, int64_t parallel) {
        double turnX = connectingTurn.get_coordinates()[0];
        double turnY = connectingTurn.get_coordinates()[1];
        if (windowOuterX <= turnX) {
            return;
        }
        // The layers the lead routes OVER (outward of the connecting turn) to reach the border.
        std::vector<const Layer*> crossedLayers;
        for (const auto& crossed : allLayers) {
            double crossedX = crossed.get_coordinates()[0];
            // No window clip: a layer OUTSIDE the window is still crossed and still
            // squeezed (see windowOuterX above).
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
            lead.routedLength = roundFloat(windowOuterX - turnX + wireOuterWidth, 9);  // the radial run
            spaces.push_back(lead);
            return;
        }

        // Crosses outer layers: route along the window edge NEAREST the end so the lead sits in the
        // extreme slots that turn-blocking frees on the crossed layers. A short stub bridges the end
        // turn up/down to that edge within its own column. ABT #229: the row on that edge comes from
        // the per-edge allocator, so each parallel's lead is its own line. Allocating in parallel
        // order matches the order the parallels' end turns stack from the edge, so for uniform wire
        // the k-th parallel's lead lines up with its own turn and the stub degenerates.
        bool turnAtTop = (turnY >= windowCenterY);
        // The lead occupies its row from the connecting turn's stub out to the border.
        auto [edgeY, runDepth] = allocateEdgeRow(windowIndexOf(connectingTurn.get_section().value_or("")),
                                                 turnAtTop, wireOuterHeight,
                                                 turnX - wireOuterWidth / 2, windowOuterX);
        for (const Layer* crossed : crossedLayers) {
            ConnectionReservedSpace space;
            space.isTerminal = true;
            space.winding = windingName;
            space.parallel = parallel;
            space.section = crossed->get_section().value_or("");
            space.layer = crossed->get_name();
            space.coordinates = {crossed->get_coordinates()[0], edgeY};
            space.dimensions = {wireOuterWidth, wireOuterHeight};
            // ABT #240: a lead crossing a layer of ANOTHER winding must clear that winding's
            // turns by the mechanical insulation that separates the two windings — the same
            // insulation the coil already builds between them. Without it the reserved band is
            // exactly one wire deep, so the crossed layer's extreme turn ends up flush against
            // the lead (measured separation 7.6e-13 um on 16_coupled_inductor_e2513_dmr95):
            // legal for same-winding packing, where adjacent turns touch by convention, but two
            // different windings may not touch.
            //
            // The clearance is NOT a margin invented here: it is the summed thickness of the
            // insulation sections the coil placed radially between the connecting turn and the
            // crossed layer, read back through get_insulation_section_thickness.
            double interWindingInsulation = 0;
            if (!crossed->get_partial_windings().empty() &&
                crossed->get_partial_windings()[0].get_winding() != windingName) {
                double crossedX = crossed->get_coordinates()[0];
                for (const auto& insulationSection : get_sections_by_type(ElectricalType::INSULATION)) {
                    // Sections come back in the REAL frame; turnX/crossedX live in the virtual
                    // frame, which is the x<->y transpose of it for contiguous layers.
                    double insulationX = layersAreContiguous ? insulationSection.get_coordinates()[1]
                                                             : insulationSection.get_coordinates()[0];
                    if (insulationX > turnX && insulationX < crossedX) {
                        interWindingInsulation +=
                            get_insulation_section_thickness(insulationSection.get_name());
                    }
                }
            }
            space.edgeDepth = runDepth + interWindingInsulation;
            space.routedLength = 0;  // space-only squeeze: the lead's copper is the drawn stub + edge run
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
            stub.routedLength = roundFloat(std::abs(stubFarEnd - turnY), 9);  // the vertical climb to the edge row
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
        lead.routedLength = roundFloat(windowOuterX - turnX + wireOuterWidth, 9);  // the edge run to the border
        lead.edgeDepth = runDepth;
        spaces.push_back(lead);
    };

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
        // ABT #577: the order leads are EMITTED is the order allocateEdgeRow hands out the
        // stacked rows, so it decides which lead is drawn on which row. In parallel order the
        // rows come out REVERSED against the turns: at the top edge the allocator gives
        // parallel 0 the row nearest the edge, while the fence-post N-filar bundle (#578/#579)
        // puts parallel 0's end turn at the BOTTOM of its bundle -- so every exit lead is drawn
        // on a SIBLING's turn row and its L-stub is driven straight through that sibling's
        // copper (measured on 06/11/14/23/24; worst 0.907 mm of interpenetration on
        // 11_pushpull). Alf, 2026-08-07: "in the output connection you have to invert the order
        // of exit for parallels ... the L segments are not needed anymore".
        //
        // Emit NEAREST-EDGE-TURN FIRST so each lead lands on its OWN turn's row: the fence-post
        // layout puts the outermost turn exactly one wire from the flange and the bundle one
        // wire apart, which is precisely the allocator's own row pitch, so the rows coincide
        // with the turn rows and every stub degenerates (the |edgeY - turnY| guard below).
        //
        // ONLY leads that cross the SAME set of layers may be permuted. Their per-layer
        // squeezes are then the same multiset of depths on the same layers, so the TURN
        // BLOCKING that reads edgeDepth is bit-identical. (A previous attempt sorted every
        // lead of the winding together -- entrance leads start on the innermost layer and
        // exits on the outermost, so they cross DIFFERENT layers; permuting those changed each
        // layer's reserved depth and MOVED THE TURNS: the 24-design sweep fell 18/25 -> 11/25.
        // Reverted in cca5a9e5; this is the surgical version.)
        for (int64_t parallel = 0; parallel < numberParallels; ++parallel) {
            auto key = std::make_pair(windingName, parallel);
            for (bool entrance : {true, false}) {
                const auto& source = entrance ? entranceTurnByWindingParallel
                                              : exitTurnByWindingParallel;
                if (!source.count(key)) {
                    continue;
                }
                const Turn& connectingTurn = source.at(key);
                double turnY = connectingTurn.get_coordinates()[1];
                bool atTop = (turnY >= windowCenterY);
                if (entrance && parallel == 0 && connectingTurn.get_layer()) {
                    // Feed the entrance edge back into the next wind's direction choice
                    // (ABT #616): the winding starts at the edge its own terminal row uses.
                    // OVERLAPPING layers only — a contiguous layer's rows run along X and a
                    // y-half signal is noise.
                    auto layersForEdge = get_layers_description().value();
                    for (const auto& l : layersForEdge) {
                        if (l.get_name() == connectingTurn.get_layer().value()) {
                            if (l.get_orientation() == WindingOrientation::OVERLAPPING) {
                                _terminalEntranceAtTop[windingName] = atTop;
                            }
                            break;
                        }
                    }
                }
                auto [crossedCount, crossedSignature] = crossedLayerCountAndSignature(connectingTurn);
                allEmissions.push_back({connectingTurn, parallel, atTop, crossedSignature,
                                        crossedCount, 0,
                                        atTop ? edgeBaseY(windowIndexOf(connectingTurn.get_section().value_or("")), true) - turnY
                                              : turnY - edgeBaseY(windowIndexOf(connectingTurn.get_section().value_or("")), false),
                                        windingName, wireOuterWidth, wireOuterHeight});
            }
        }
    }
    {
        std::map<std::tuple<size_t, bool, std::string>, size_t> groupRankByKey;
        for (size_t i = 0; i < allEmissions.size(); ++i) {
            auto groupKey = std::make_tuple(windowIndexOf(allEmissions[i].turn.get_section().value_or("")),
                                            allEmissions[i].atTop, allEmissions[i].crossedSignature);
            auto found = groupRankByKey.find(groupKey);
            if (found == groupRankByKey.end()) {
                groupRankByKey[groupKey] = i;
                allEmissions[i].groupRank = i;
            }
            else {
                allEmissions[i].groupRank = found->second;
            }
        }
        std::stable_sort(allEmissions.begin(), allEmissions.end(),
                         [](const LeadEmission& a, const LeadEmission& b) {
                             if (a.crossedCount != b.crossedCount) {
                                 return a.crossedCount > b.crossedCount;  // widest-crossing first
                             }
                             if (a.groupRank != b.groupRank) {
                                 return a.groupRank < b.groupRank;
                             }
                             return a.edgeDistance < b.edgeDistance;
                         });
        for (const auto& emission : allEmissions) {
            addTerminalLead(emission.windingName, emission.wireW, emission.wireH,
                            emission.turn, emission.parallel);
        }
    }

    // SECOND PASS: inter-layer and inter-section links, after every winding's terminal rows are
    // allocated, so the shared continuation bands stack against the final terminal occupancy.
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        auto windingName = get_functional_description()[windingIndex].get_name();
        double wireOuterWidth = wires[windingIndex].get_maximum_outer_width();
        double wireOuterHeight = wires[windingIndex].get_maximum_outer_height();
        if (layersAreContiguous) {
            std::swap(wireOuterWidth, wireOuterHeight);
        }
        int64_t numberParallels = int64_t(get_number_parallels(windingIndex));

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
            WindingOrder windingOrder = get_winding_order(windingLayers[i].get_section().value());
            // ABT #615 (Alf, 2026-08-09): EVERY inter-section continuation crossing intervening
            // sections routes along the window edge — drawn (blue), row-allocated and BLOCKING the
            // crossed sections — regardless of winding order. The ABT #492 FRONT_YZ face dragback
            // is superseded for this case: invisible no-cost returns left the crossed section's end
            // turns inside the return's climb corridor (ABT #612, 0.267 mm vs 0.500 needed on the
            // PSPS pair). Adjacent-layer Z dragbacks (nothing intervening) keep their in-plane
            // diagonal.
            bool routesAlongEdge = crossesIntervening;
            size_t routeWindowIndex = windowIndexOf(windingLayers[i].get_section().value_or(""));

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

                // ABT #615: ALL inter-section continuations on an edge share ONE band (Alf: the
                // corridor blocking part of the crossed section "can then be reused by the inter
                // section connection in secondary"). The band is a HEIGHT reservation; in 3D the
                // connections are tangential chords separated in AZIMUTH, so radial overlap inside
                // the band is fine — unlike terminal leads, which fan on one connection plane and
                // therefore still need one row each. The band claims only its runs' radial span
                // (Alf, 2026-08-09), so it shares the leads' row height wherever they don't meet.
                // ABT #615 stage 2: the connection edge FOLLOWS THE EXIT TURN (Alf's alternation
                // — hop 1 exits top, hop 2 exits bottom, ...). With direction alternation the
                // receiving section starts on the same edge, so both stubs stay short.
                const bool routeAtTop =
                    exitTurn.get_coordinates()[1] >= windowCenterTurnAxis;
                // An edge continuation is copper too: same margin inset as the terminal rows.
                double routeBaseY = edgeBaseY(windowIndexOf(exitTurn.get_section().value_or("")), routeAtTop);
                double routeEdgeY = roundFloat(routeAtTop ? routeBaseY - wireOuterHeight / 2
                                                          : routeBaseY + wireOuterHeight / 2, 9);
                double runDepth = 0;
                if (routesAlongEdge) {
                    // ABT #615, Alf 2026-08-09: the band claims ONLY the radial span its runs
                    // actually cover -- "the terminal wire is not there [inward of its connecting
                    // turn], and we don't need to reserve the space blocker". A terminal row blocks
                    // from its connecting turn OUTWARD, so a band whose runs live inward of the
                    // leads SHARES the leads' row height instead of stacking under it (span-aware
                    // allocateEdgeRow). Reuse registers each new run's span on the band's row, so
                    // later terminal leads still see the true occupancy and stack below only where
                    // they genuinely overlap it.
                    const double spanLo = std::min(exitTurn.get_coordinates()[0],
                                                   entryTurn.get_coordinates()[0]) - wireOuterWidth / 2;
                    const double spanHi = std::max(exitTurn.get_coordinates()[0],
                                                   entryTurn.get_coordinates()[0]) + wireOuterWidth / 2;
                    const int routeEdge = routeAtTop ? 0 : 1;
                    auto bandKey = std::make_pair(routeWindowIndex, routeEdge);
                    auto existingBand = continuationBand.find(bandKey);
                    if (existingBand != continuationBand.end()
                        && existingBand->second.height + 1e-12 >= wireOuterHeight) {
                        routeEdgeY = existingBand->second.edgeY;
                        runDepth = existingBand->second.runDepth;
                        edgeRows[{routeWindowIndex, routeEdge}][existingBand->second.rowIndex]
                            .spans.push_back({spanLo, spanHi});
                    }
                    else {
                        std::tie(routeEdgeY, runDepth) =
                            allocateEdgeRow(routeWindowIndex, routeAtTop, wireOuterHeight, spanLo, spanHi);
                        size_t rowIndex = 0;
                        auto& rows = edgeRows[{routeWindowIndex, routeEdge}];
                        for (size_t r = 0; r < rows.size(); ++r) {
                            if (std::abs((rows[r].depthBefore + rows[r].height) - runDepth) < 1e-12) {
                                rowIndex = r;
                                break;
                            }
                        }
                        continuationBand[bandKey] = {routeEdgeY, runDepth, wireOuterHeight, rowIndex};
                    }
                }

                // Per-layer squeeze: a U parallel's interleaved continuation crosses (and squeezes)
                // each intervening layer in its allocated edge row. These entries (layer set) drive
                // the filling factor and the turn blocking and are NOT drawn — the link itself is
                // drawn below. ABT #615: Z interleaved continuations now squeeze too — their run is
                // in-window like U's, and the blocked corridor is exactly what keeps the crossed
                // section's end turns out of the return's path (ABT #612).
                if (routesAlongEdge) {
                    for (const Layer* crossed : interveningLayers) {
                        ConnectionReservedSpace squeeze;
                        squeeze.section = crossed->get_section().value();
                        squeeze.layer = crossed->get_name();
                        squeeze.winding = windingName;
                        squeeze.parallel = parallel;
                        squeeze.coordinates = {crossed->get_coordinates()[0], routeEdgeY};
                        squeeze.dimensions = {wireOuterWidth, wireOuterHeight};
                        squeeze.routedLength = 0;  // space-only: the copper is the drawn stubs + edge run below
                        squeeze.edgeDepth = runDepth;
                        spaces.push_back(squeeze);
                    }
                }

                // NO ENDPOINT SQUEEZE — Alf, 2026-08-09 (ABT #615): "the blue connections from a
                // section should not block space in its own section, or the receiving section, just
                // the sections in between." The run is the endpoint turns' OWN wire continuing —
                // charging their layers took two top slots from the source section's outer layers
                // (25_psps: P0 layers 1-2 topped at +0.325 instead of ~+0.97), and a first-iteration
                // endpoint marker could seed the monotone depth map with source-section depths that
                // outlived the layer layout. The old #229 corner-overhang rationale dies with the
                // tangential-chord 3D realization: the connection leaves the end turn as a chord at
                // the corridor height, not as a stub-and-corner over the endpoint columns.
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

                // ABT #615: the FRONT_YZ face-dragback emission for Z inter-section returns lived
                // here (ABT #492). Superseded: a return crossing intervening sections now routes
                // in-window along the edge — same drawn run, squeezes and blocking as the U
                // interleaved continuation, emitted by the crossesIntervening branch below.
                if (windingOrder == WindingOrder::Z && !crossesIntervening) {
                    // Z between ADJACENT layers: the classic dragback, drawn as the single in-plane
                    // diagonal from one turn straight to the next (a rotated rectangle from centre
                    // to centre). Nothing intervenes, so it displaces and reserves nothing.
                    double deltaX = x2 - x1;
                    double deltaY = y2 - y1;
                    double length = std::sqrt(deltaX * deltaX + deltaY * deltaY);
                    ConnectionReservedSpace diagonal;
                    diagonal.winding = windingName;
                    diagonal.parallel = parallel;
                    diagonal.section = windingLayers[i].get_section().value_or("");
                    diagonal.layer = "";
                    diagonal.coordinates = {roundFloat((x1 + x2) / 2, 9), roundFloat((y1 + y2) / 2, 9)};
                    diagonal.dimensions = {roundFloat(length, 9), wireOuterHeight};
                    diagonal.routedLength = roundFloat(length, 9);  // centre-to-centre dragback hop
                    diagonal.rotation = roundFloat(std::atan2(deltaY, deltaX) * 180.0 / std::numbers::pi, 6);
                    spaces.push_back(diagonal);
                }
                else if (crossesIntervening) {
                    // U interleaved continuation: route along the window edge it reserves — a vertical
                    // stub from each turn up/down to the edge, and a horizontal run along the edge
                    // across the intervening layer(s) — so the wire never cuts through the crossed
                    // layer's turns (centre-to-centre would clip them).
                    // Each segment's copper length is passed EXPLICITLY (the verticals run along the
                    // TURN axis; inferring their length from an orientation-derived index counted
                    // each stub as one wire width).
                    auto pushLink = [&](double cx, double cy, double w, double h, double copperLength,
                                        double depth = 0) {
                        ConnectionReservedSpace seg;
                        seg.winding = windingName;
                        seg.parallel = parallel;
                        seg.section = windingLayers[i].get_section().value_or("");
                        seg.layer = "";
                        seg.coordinates = {roundFloat(cx, 9), roundFloat(cy, 9)};
                        seg.dimensions = {roundFloat(w, 9), roundFloat(h, 9)};
                        seg.routedLength = roundFloat(copperLength, 9);
                        seg.edgeDepth = depth;
                        spaces.push_back(seg);
                    };
                    // Verticals start at the turn CENTRE (y1 / y2), not covering the whole turn, and
                    // overlap the horizontal by half a wire at the corner.
                    if (std::abs(routeEdgeY - y1) > 0.5 * wireOuterHeight) {
                        double far1 = routeEdgeY + ((routeEdgeY >= y1) ? 1.0 : -1.0) * wireOuterHeight / 2;
                        pushLink(x1, (y1 + far1) / 2, wireOuterWidth, std::abs(far1 - y1), std::abs(far1 - y1));
                    }
                    pushLink((x1 + x2) / 2, routeEdgeY, std::abs(x2 - x1) + wireOuterWidth, wireOuterHeight,
                             std::abs(x2 - x1) + wireOuterWidth, runDepth);
                    if (std::abs(y2 - routeEdgeY) > 0.5 * wireOuterHeight) {
                        double far2 = routeEdgeY + ((routeEdgeY >= y2) ? 1.0 : -1.0) * wireOuterHeight / 2;
                        pushLink(x2, (y2 + far2) / 2, wireOuterWidth, std::abs(far2 - y2), std::abs(far2 - y2));
                    }
                }
                else {
                    // ABT #608 RETRACTED (Alf, 2026-08-08): the U tangential link does NOT cost
                    // the landing layer a turn slot. The "arriving crossing" is not a separate
                    // occupant of the top slot -- it IS the layer's top turn: the turn starts with
                    // the tangential segment at the source turn's height, its revolution is the top
                    // turn, and it ends at its own station at that same height. Reserving a slot
                    // for the crossing AND dropping the first turn a pitch counted the same copper
                    // twice: on the exact cut plane a non-first U layer shows M circles, not M+1
                    // (the tangential segment's only cut-plane copper is its start point, which
                    // coincides with the previous layer's last station), and a real serpentine
                    // loses no capacity per layer -- its first turn nests level with the previous
                    // layer's last. Crossing accounting stays N + 1 per parallel per winding, the
                    // global RealWindingCrossingBump, and nothing else.
                    const bool sameSection =
                        windingLayers[i].get_section() == windingLayers[i + 1].get_section();

                    // U adjacent layers: route the wire orthogonally — a horizontal stretch from the
                    // source turn out to the destination layer's radial position, then a vertical
                    // stretch down/up to the destination turn when the two are at different heights.
                    // The horizontal runs half a wire past the corner and the vertical is pulled back
                    // half a wire so the bend reads as one continuous wire.
                    //
                    // ...EXCEPT for the same-section tangential link, which is horizontal ONLY. It
                    // is manufactured as a tangential run at constant height — Alf, 2026-08-08:
                    // "vertical connections are for dragbacks; the connection to the external next
                    // layer must be tangential". Normally the two stations sit at the same height
                    // and this guard is dormant; when spread layers put them apart, any residual
                    // height difference is the destination turn's own helix, not a drawn stub —
                    // bridging it would draw copper the winder never lays, and would contradict the
                    // 3D, which builds this link as a constant-height tangential segment.
                    bool needVertical = std::abs(y2 - y1) > 0.5 * wireOuterHeight
                                        && !(windingOrder == WindingOrder::U && sameSection);
                    ConnectionReservedSpace horizontal;
                    horizontal.winding = windingName;
                    horizontal.parallel = parallel;
                    horizontal.section = windingLayers[i].get_section().value_or("");
                    horizontal.layer = "";
                    if (needVertical) {
                        double horizontalDirection = (x2 >= x1) ? 1.0 : -1.0;
                        horizontal.coordinates = {roundFloat((x1 + x2) / 2 + horizontalDirection * wireOuterWidth / 4, 9), roundFloat(y1, 9)};
                        horizontal.dimensions = {roundFloat(std::abs(x2 - x1) + wireOuterWidth / 2, 9), wireOuterHeight};
                        horizontal.routedLength = roundFloat(std::abs(x2 - x1) + wireOuterWidth / 2, 9);
                    }
                    else {
                        horizontal.coordinates = {roundFloat((x1 + x2) / 2, 9), roundFloat(y1, 9)};
                        horizontal.dimensions = {roundFloat(std::abs(x2 - x1), 9), wireOuterHeight};
                        horizontal.routedLength = roundFloat(std::abs(x2 - x1), 9);
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
                        // The vertical stretch runs along the TURN axis: its copper is its own
                        // extent, not a wire width.
                        vertical.routedLength = roundFloat(std::abs(y2 - y1) - wireOuterHeight / 2, 9);
                        spaces.push_back(vertical);
                    }
                }
            }
        }
    }
    // Back to the coil's own frame for the contiguous case. Everything else in a rectangular winding
    // window follows ONE convention — dimensions[0] is the X extent and dimensions[1] the Y extent,
    // for sections, layers and turns alike — and these rectangles obey it too. The contiguous pass
    // above ran with the layer axis mapped onto x, so coming back is a reflection across y = x: the
    // centre's x and y trade places, and a rectangle whose long axis sits at angle θ maps to one at
    // 90 − θ.
    //
    // For an axis-aligned rectangle — which is every marker that names a layer, plus the edge runs —
    // that reflection is simply "the X and Y extents trade places", so it is stored that way: swapped
    // extents and NO rotation. Keeping it as an unswapped rectangle rotated 90° would describe the
    // same geometry while quietly introducing a second convention, where dimensions[0] means Y and
    // consumers must know to undo a rotation. Only the Z diagonals are genuinely not axis-aligned, and
    // only they carry an angle.
    if (layersAreContiguous) {
        for (auto& space : spaces) {
            if (space.coordinates.size() >= 2) {
                std::swap(space.coordinates[0], space.coordinates[1]);
            }
            double reflectedRotation = roundFloat(90.0 - space.rotation, 6);
            double angleFromXAxis = std::fmod(std::fmod(reflectedRotation, 180.0) + 180.0, 180.0);
            if (std::abs(angleFromXAxis - 90.0) < 1e-9) {
                // Long axis now lies along Y: express it as an axis-aligned rectangle in (X, Y) order.
                if (space.dimensions.size() >= 2) {
                    std::swap(space.dimensions[0], space.dimensions[1]);
                }
                space.rotation = 0;
            }
            else if (angleFromXAxis < 1e-9) {
                space.rotation = 0;  // already along X, extents already in (X, Y) order
            }
            else {
                space.rotation = reflectedRotation;
            }
        }
    }
    return spaces;
}

std::map<std::string, std::pair<uint64_t, uint64_t>> Coil::compute_connection_blocked_slots_per_layer(
        std::map<std::string, std::pair<double, double>>* freshDepths) {
    std::map<std::string, std::pair<uint64_t, uint64_t>> blockedSlotsPerLayer;  // layer name -> {top, bottom}
    if (!get_layers_description()) {
        return blockedSlotsPerLayer;
    }
    // ABT #187: the top/bottom slot model below is rectangular (axial y). Toroidal blocking is
    // ANGULAR — computed and applied by align_blocked_ring_turns from the ring crossing markers —
    // so a round window contributes nothing here (which also keeps wind()'s rectangular fixpoint
    // loop a no-op for toroids, as before).
    if (resolve_bobbin().get_winding_window_shape() == WindingWindowShape::ROUND) {
        return blockedSlotsPerLayer;
    }
    auto layers = get_layers_description().value();
    auto wires = get_wires();
    // ABT #427: blocking is modelled along each layer's TURN axis — the direction its turns stack, and
    // so the one a crossing lead takes a slot out of. That is Y for an OVERLAPPING layer (turns stack
    // axially) and X for a CONTIGUOUS one, which is one wire tall by construction with its turns
    // running along its width (see wind_by_layers). Coordinates and dimensions are plain X/Y
    // throughout — layers, turns and the connection markers all share the one convention — so the
    // orientation only picks WHICH INDEX to read, never a different frame. The returned pair is
    // {high side, low side} of that axis: top/bottom for overlapping, right/left for contiguous.
    // Indexing by axis rather than mirroring the loop keeps the two orientations from drifting apart.
    std::map<std::string, size_t> layerTurnAxis;
    std::map<std::string, double> layerCenterOnTurnAxis;
    std::map<std::string, double> layerWirePitch;  // crossed layer's own wire, along its turn axis
    std::map<std::string, double> layerExtent;     // crossed layer's own extent, same axis
    for (const auto& layer : layers) {
        size_t turnAxis = (layer.get_orientation() == WindingOrientation::OVERLAPPING) ? 1 : 0;
        layerTurnAxis[layer.get_name()] = turnAxis;
        layerCenterOnTurnAxis[layer.get_name()] = layer.get_coordinates()[turnAxis];
        layerExtent[layer.get_name()] = layer.get_dimensions()[turnAxis];
        // Insulation layers have no partial windings; only conduction layers have a wire/turn pitch.
        if (layer.get_type() == ElectricalType::CONDUCTION && !layer.get_partial_windings().empty()) {
            size_t windingIndex = get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
            layerWirePitch[layer.get_name()] = (turnAxis == 1)
                ? wires[windingIndex].get_maximum_outer_height()
                : wires[windingIndex].get_maximum_outer_width();
        }
    }
    // ABT #229: every edge-routed run (terminal lead or U interleaved continuation) now carries the
    // DEPTH the per-edge row allocator assigned it (distance from the window edge to the run's inner
    // side, see edgeDepth in ConnectionReservedSpace) — the drawn geometry is the single source of
    // truth. A crossed layer must free the band from the window edge down to the DEEPEST run that
    // crosses it: blocked height = max over its markers of edgeDepth, converted to turn slots of the
    // CROSSED layer's own wire (ceil: a deep stack over a thin layer displaces several thin turns, a
    // shallow stack over a thick layer still costs at least one thick turn). Since parallels'
    // (and different windings') runs stack in distinct rows, this replaces the old per-parallel
    // max-then-SUM rule — which coincided all runs on one line and merged different windings that
    // shared a parallel index. ABT #615: EVERY layer-naming marker blocks — Z inter-section
    // continuations route in-window along the edge like U's (the ABT #492 FRONT_YZ model is
    // superseded), so their squeezes carry real corridor depths; only the adjacent-layer Z
    // diagonal (which names no layer) reserves nothing.
    std::map<std::string, std::pair<double, double>> maxRunDepth;  // layer -> {top, bottom}
    for (const auto& space : get_connection_reserved_spaces()) {
        if (space.layer.empty()) {
            continue;
        }
        auto found = layerCenterOnTurnAxis.find(space.layer);
        if (found == layerCenterOnTurnAxis.end()) {
            continue;
        }
        // Defensive: a marker without an allocated depth still costs its own thickness (one row),
        // measured along the crossed layer's turn axis — the axis this lead takes a slot out of.
        double depth = std::max(space.edgeDepth, space.dimensions[layerTurnAxis.at(space.layer)]);
        auto& edges = maxRunDepth[space.layer];
        if (space.coordinates[layerTurnAxis.at(space.layer)] >= found->second) {
            edges.first = std::max(edges.first, depth);
        }
        else {
            edges.second = std::max(edges.second, depth);
        }
    }
    for (const auto& [layerName, edges] : maxRunDepth) {
        double crossedWirePitch = layerWirePitch.count(layerName) ? layerWirePitch.at(layerName) : 0.0;
        if (crossedWirePitch <= 0) {
            continue;
        }
        // The layer loses the CAPACITY the reserved band actually costs it, not the band rounded
        // up to whole slots at each edge. Rounding each edge up charged a THIN lead a whole THICK
        // slot: Alf, 2026-08-08, on a reordered 13_current_sense -- the secondary's 0.1125 mm lead
        // took a full 0.534 mm primary slot at each edge (1.068 mm of an 1.8293 mm window), so the
        // primary kept 0.801 mm, fitting ONE turn, and its two real-winding crossings could not
        // share a layer. Comparing the packing WITH and WITHOUT the bands charges exactly the
        // turns that no longer fit: a fat lead over a fine layer still displaces the several fine
        // turns it really covers (unchanged, 10 slots in the original order), while a thin lead
        // over a thick layer costs a slot only when it genuinely pushes one out.
        // Measured against the layer's UNBLOCKED extent. A layer that has already been shrunk by a
        // previous iteration of the fixpoint carries a smaller extent, while the depths are absolute
        // (from the window edge) — so comparing the two counts the same reservation twice, saturates
        // (`fits(extent - depths)` hits 0), and reports FEWER blocked slots than the first iteration
        // did. The loop only survived that because it accumulates the slot counts monotonically, but
        // the DEPTHS kept growing underneath, and the placement pass hugs the depths: a layer could
        // end up spreading more turns than its own reserved span holds. Adding back the room the
        // layer surrendered (recorded by wind_by_layers, already capped to what it really gave up)
        // restores the ideal extent, so every iteration measures the same geometry the first one did
        // and the two currencies agree.
        double extent = layerExtent.count(layerName) ? layerExtent.at(layerName) : 0.0;
        auto surrendered = _connectionBlockedRoomPerLayer.find(layerName);
        if (surrendered != _connectionBlockedRoomPerLayer.end()) {
            extent += surrendered->second;
        }
        const auto fits = [&](double usable) -> uint64_t {
            return usable > 0.0 ? uint64_t(std::floor(usable / crossedWirePitch + 1e-9)) : 0u;
        };
        const uint64_t freeCapacity = fits(extent);
        const uint64_t blockedCapacity =
            freeCapacity - std::min(freeCapacity, fits(extent - edges.first - edges.second));
        // Attribute the lost slots to the edges that caused them (deepest first), so the placement
        // pass below still knows which side to pack against.
        uint64_t topSlots = 0, bottomSlots = 0;
        if (blockedCapacity > 0) {
            const double totalDepth = edges.first + edges.second;
            if (totalDepth > 1e-12) {
                topSlots = uint64_t(std::llround(blockedCapacity * (edges.first / totalDepth)));
                topSlots = std::min<uint64_t>(topSlots, blockedCapacity);
            }
            bottomSlots = blockedCapacity - topSlots;
            if (topSlots == 0 && edges.first > edges.second) std::swap(topSlots, bottomSlots);
        }
        blockedSlotsPerLayer[layerName] = {topSlots, bottomSlots};
        // MKF_BLOCKING_DIAG: what each layer surrendered and to which edge. The fixpoint is easy to
        // misread from the wound geometry alone (a layer can look one slot short at the wrong edge),
        // so print the currency it actually decided in: unblocked extent, the two reserved depths,
        // and the capacity that comparison costs.
        if (std::getenv("MKF_BLOCKING_DIAG")) {
            std::cerr << "[blocking] " << layerName << " extent=" << extent * 1e3
                      << " pitch=" << crossedWirePitch * 1e3 << " depths={" << edges.first * 1e3
                      << "," << edges.second * 1e3 << "} free=" << freeCapacity
                      << " blocked=" << blockedCapacity << " -> {" << topSlots << "," << bottomSlots << "}\n";
        }
        // Same keys as the slot map, by construction: the continuous depths the ceil above rounded up
        // from. The placement pass hugs THESE, so turns reach the crossing runs exactly instead of
        // stopping at the whole-slot grid.
        if (freshDepths != nullptr) {
            (*freshDepths)[layerName] = edges;
        }
    }
    return blockedSlotsPerLayer;
}

std::map<std::string, std::pair<double, double>> Coil::compute_u_landing_extra_depths() {
    // See the header comment (ABT #608 final form): a non-first U layer's first station sits one
    // wire OD past the tangential arrival — placement only, via the depth map, never a slot.
    std::map<std::string, std::pair<double, double>> extraDepths;
    if (!get_layers_description() || !get_turns_description()) {
        return extraDepths;
    }
    auto bobbin = resolve_bobbin();
    if (bobbin.get_winding_window_shape() != WindingWindowShape::RECTANGULAR) {
        return extraDepths;
    }
    auto windingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
    std::array<double, 2> windowCenterPerAxis = {windingWindow.get_coordinates().value()[0],
                                                 windingWindow.get_coordinates().value()[1]};
    std::array<double, 2> windowHalfSizePerAxis = {windingWindow.get_width().value() / 2,
                                                   windingWindow.get_height().value() / 2};
    auto wires = get_wires();
    auto layers = get_layers_description().value();
    auto turns = get_turns_description().value();
    // ABT #683: the margin the landing layer's own section reserves, held by name.
    std::map<std::string, std::vector<double>> marginBySection;
    if (get_sections_description()) {
        auto sectionsForMargin = get_sections_description().value();
        for (const auto& section : sectionsForMargin) {
            marginBySection[section.get_name()] = resolve_margin(section);
        }
    }

    // Layers in ELECTRICAL order (the order their turns are first wound), with each layer's turns.
    std::map<std::string, size_t> layerElectricalOrder;
    std::map<std::string, std::vector<const Turn*>> turnsByLayer;
    size_t order = 0;
    for (const auto& turn : turns) {
        if (!turn.get_layer()) {
            continue;
        }
        const std::string layerName = turn.get_layer().value();
        if (layerElectricalOrder.find(layerName) == layerElectricalOrder.end()) {
            layerElectricalOrder[layerName] = order++;
        }
        turnsByLayer[layerName].push_back(&turn);
    }

    // Conduction layers grouped by SECTION (the tangential link is a same-section transition),
    // U-order sections only, in electrical order within the section.
    std::map<std::string, std::vector<const Layer*>> sectionLayers;
    for (const auto& layer : layers) {
        if (layer.get_type() != ElectricalType::CONDUCTION || !layer.get_section()
            || !layerElectricalOrder.count(layer.get_name())) {
            continue;
        }
        if (get_winding_order(layer.get_section().value()) != WindingOrder::U) {
            continue;
        }
        sectionLayers[layer.get_section().value()].push_back(&layer);
    }
    for (auto& [sectionName, sorted] : sectionLayers) {
        std::sort(sorted.begin(), sorted.end(), [&](const Layer* a, const Layer* b) {
            return layerElectricalOrder.at(a->get_name()) < layerElectricalOrder.at(b->get_name());
        });
        for (size_t i = 0; i + 1 < sorted.size(); ++i) {
            const Layer& previous = *sorted[i];
            const Layer& landing = *sorted[i + 1];
            size_t windingIndex =
                get_winding_index_by_name(landing.get_partial_windings()[0].get_winding());
            // EXCEPTION (Alf): a landing turn that is the section's LAST stays level — nothing
            // follows it, so it winds in the connection's own height. That is the section's final
            // layer holding one turn per parallel and nothing more.
            bool landingIsSectionEnd =
                (i + 2 == sorted.size())
                && turnsByLayer.at(landing.get_name()).size()
                       <= size_t(get_number_parallels(windingIndex));
            if (landingIsSectionEnd) {
                continue;
            }
            size_t turnAxis = (landing.get_orientation() == WindingOrientation::OVERLAPPING) ? 1 : 0;
            double wireOD = (turnAxis == 1) ? wires[windingIndex].get_maximum_outer_height()
                                            : wires[windingIndex].get_maximum_outer_width();
            // The arrival height is the previous layer's ELECTRICALLY LAST turn (the one the
            // tangential chunk leaves); with parallels wound side by side, the link leaves the
            // turn nearest the landing edge, so take the extreme over the last bundle.
            const auto& previousTurns = turnsByLayer.at(previous.get_name());
            double lastY = previousTurns.back()->get_coordinates()[turnAxis];
            bool landsAtHighSide = lastY >= landing.get_coordinates()[turnAxis];
            double arrivalY = lastY;
            int64_t numberParallels = get_number_parallels(windingIndex);
            for (size_t k = previousTurns.size() - std::min<size_t>(previousTurns.size(),
                                                                    size_t(numberParallels));
                 k < previousTurns.size(); ++k) {
                double y = previousTurns[k]->get_coordinates()[turnAxis];
                arrivalY = landsAtHighSide ? std::max(arrivalY, y) : std::min(arrivalY, y);
            }
            // Place the landing layer's first station one OD past the arrival: the aligned spread
            // puts the end station half an OD inside the span, so the span boundary sits at
            // arrival -/+ OD/2. CAPPED to what the layer's own turns leave free: a FULL layer
            // cannot descend without its far-end turn overflowing into the window edge or the
            // terminal-lead rows (pigeonhole — the drop and the turn count cannot both hold, and
            // capacity is deliberately NOT charged, per the N_layer+1 retraction). Physically
            // that IS the serpentine: where the window has slack the landing descends, and a
            // packed layer lands level with the previous layer's last turn. The cap never goes
            // below the depth already accumulated (lead rows keep their reservation).
            double windowHigh = windowCenterPerAxis[turnAxis] + windowHalfSizePerAxis[turnAxis];
            double windowLow = windowCenterPerAxis[turnAxis] - windowHalfSizePerAxis[turnAxis];
            // ABT #683: measured from the band the turns may actually use, which is the window
            // INSET BY THE SECTION'S MARGIN (ABT #676) — align_blocked_layer_turns subtracts these
            // depths from that same inset band, so taking them from the raw window counted the
            // margin twice and pushed the landing a whole margin further in.
            {
                auto marginIt = marginBySection.find(landing.get_section().value());
                if (marginIt != marginBySection.end()) {
                    windowHigh -= (turnAxis == 1) ? marginIt->second[0] : marginIt->second[1];
                    windowLow += (turnAxis == 1) ? marginIt->second[1] : marginIt->second[0];
                }
            }
            double copperNeeded = double(turnsByLayer.at(landing.get_name()).size()) * wireOD;
            std::pair<double, double> accumulated{0.0, 0.0};
            auto accumulatedIt = _connectionBlockedDepthPerLayer.find(landing.get_name());
            if (accumulatedIt != _connectionBlockedDepthPerLayer.end()) {
                accumulated = accumulatedIt->second;
            }
            auto& edges = extraDepths[landing.get_name()];
            // ABT #683 (Alf): "in U windings the second layer starts just after the first layer,
            // that is the point of U winding, so we just need to connect to it horizontally". The
            // landing layer's END STATION therefore sits AT the arrival, not one wire past it:
            // under the SPREAD alignment real winding forces on every layer that station is half
            // an OD inside the span, so the span boundary is arrival -/+ OD/2 and the depth
            // reaches DOWN to it.
            //
            // ONE CONDUCTOR ONLY. With parallels wound side by side, landing level puts every
            // parallel's link at the station its neighbour departs from and the landing
            // revolutions overlap — which is exactly what the ABT #608 descent was introduced to
            // stop, and what "multi-layer multi-parallel builds collision-free" guards (8t x 2p:
            // the two links converged to 0.822 mm against a 0.9 mm envelope). A single conductor
            // has nothing to overlap, so it lands level; parallels keep the descent until the
            // landing stations can be placed per PARALLEL rather than per layer.
            const bool landsLevel = numberParallels == 1;
            if (landsAtHighSide) {
                double ideal = landsLevel ? windowHigh - arrivalY - wireOD / 2
                                          : windowHigh - arrivalY + wireOD / 2;
                double maxAllowed = (windowHigh - windowLow) - accumulated.second - copperNeeded;
                double depth = std::max(accumulated.first, std::min(ideal, maxAllowed));
                edges.first = std::max(edges.first, roundFloat(depth, 9));
            }
            else {
                // Mirror of the high side.
                double ideal = landsLevel ? arrivalY - windowLow - wireOD / 2
                                          : arrivalY - windowLow + wireOD / 2;
                double maxAllowed = (windowHigh - windowLow) - accumulated.first - copperNeeded;
                double depth = std::max(accumulated.second, std::min(ideal, maxAllowed));
                edges.second = std::max(edges.second, roundFloat(depth, 9));
            }
            if (std::getenv("MKF_BLOCKING_DIAG")) {
                std::cerr << "[u-landing] " << landing.get_name() << " arrival=" << arrivalY * 1e3
                          << (landsAtHighSide ? " (top)" : " (bottom)")
                          << " accumulated={" << accumulated.first * 1e3 << ","
                          << accumulated.second * 1e3 << "} copper=" << copperNeeded * 1e3
                          << " -> depths={" << edges.first * 1e3 << "," << edges.second * 1e3 << "}\n";
            }
        }
    }
    return extraDepths;
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
        // Each parallel of a bifilar/N-filar group is wound side by side, so a layer holds an equal
        // number of turns of every parallel: its capacity in PER-PARALLEL turns ("rows") is the
        // physical capacity divided by the parallel count. We redistribute in per-parallel turns
        // (matching get_number_turns, which is per parallel) so interior sections end on whole blocked
        // layers and the remainder is pushed to the outermost section. (K=1 keeps the previous result.)
        int64_t numberParallels = int64_t(get_number_parallels(windingIndex));

        // This winding's conduction sections, in wound (radial) order. ABT #427: both layer
        // orientations block, so both redistribute; the per-layer capacity below reads the section
        // extent and wire dimension along whichever axis that section's turns run.
        std::vector<size_t> windingSections;
        for (size_t s = 0; s < sections.size(); ++s) {
            if (sections[s].get_type() == ElectricalType::CONDUCTION
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
            // Physical turns per layer; at least one row of every parallel. The turns run along the
            // section's HEIGHT when its layers overlap and along its WIDTH when they are contiguous, so
            // both the section extent and the wire dimension are read on that same axis.
            size_t turnAxis = (section.get_layers_orientation() == WindingOrientation::OVERLAPPING) ? 1 : 0;
            double wirePitch = (turnAxis == 1) ? wirePerWinding[windingIndex].get_maximum_outer_height()
                                               : wirePerWinding[windingIndex].get_maximum_outer_width();
            if (wirePitch <= 0) {
                continue;
            }
            uint64_t maximumTurnsPerLayer = std::max<uint64_t>(numberParallels, uint64_t(std::floor(section.get_dimensions()[turnAxis] / wirePitch)));
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
    // ABT #427: the band the turns are spread along is the layer's TURN axis — the window's HEIGHT for
    // an OVERLAPPING layer (turns stack axially) and its WIDTH for a CONTIGUOUS one (turns run
    // laterally). Indexed by axis rather than mirrored, so the two orientations cannot drift apart.
    std::array<double, 2> windowCenterPerAxis = {windingWindow.get_coordinates().value()[0],
                                                 windingWindow.get_coordinates().value()[1]};
    std::array<double, 2> windowHalfSizePerAxis = {windingWindow.get_width().value() / 2,
                                                   windingWindow.get_height().value() / 2};
    auto wires = get_wires();
    auto layers = get_layers_description().value();
    auto turns = get_turns_description().value();
    // Held by name: the margin each layer must stay clear of belongs to its own section.
    std::vector<Section> sectionsForMargins;
    if (get_sections_description()) {
        sectionsForMargins = get_sections_description().value();
    }

    for (auto& layer : layers) {
        if (layer.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        size_t turnAxis = (layer.get_orientation() == WindingOrientation::OVERLAPPING) ? 1 : 0;
        double windowHighSide = windowCenterPerAxis[turnAxis] + windowHalfSizePerAxis[turnAxis];
        double windowLowSide = windowCenterPerAxis[turnAxis] - windowHalfSizePerAxis[turnAxis];
        // ABT #676: margin tape is not free space. Spreading against the raw window put the
        // re-aligned turns straight into the margin band, and the layout that followed carried the
        // section — and its terminal leads — in with them, while the painter went on drawing a
        // margin the geometry had stopped honouring. The band this layer may use is the window
        // inset by ITS OWN section's margin. margin[0] is "top or left", so it insets the HIGH
        // side on the turn (y) axis and the LOW side on the layer (x) axis; margin[1] mirrors it.
        for (const auto& marginSection : sectionsForMargins) {
            if (!layer.get_section() || marginSection.get_name() != layer.get_section().value()) {
                continue;
            }
            auto sectionMargin = resolve_margin(marginSection);
            windowHighSide -= (turnAxis == 1) ? sectionMargin[0] : sectionMargin[1];
            windowLowSide += (turnAxis == 1) ? sectionMargin[1] : sectionMargin[0];
            break;
        }
        // Reposition this layer's turns to leave exactly the blocked slots free at each edge: spread the
        // turns evenly across the UNBLOCKED band [windowBottom + blockedBottom slots, windowTop −
        // blockedTop slots]. Even spacing packs a full layer (step == wire height) and spreads a partial
        // one — in both cases keeping every turn clear of the slots the connection leads route through.
        // delimit re-centres each layer over the full window height, which would otherwise drop the edge
        // turns back under the leads, so this runs after delimit and overrides its vertical placement.
        // Keyed on the DEPTHS, not the slot counts: a U landing layer (ABT #608 final form) carries
        // a placement-only depth with NO blocked slots — its turns still must move off the arrival
        // band, exactly as lead-crossed layers move off the lead rows. The effective depth per edge
        // is the max of the (monotone) marker depths and the (recomputed) U landing overlay.
        std::pair<double, double> effectiveDepths{0.0, 0.0};
        auto foundDepths = _connectionBlockedDepthPerLayer.find(layer.get_name());
        if (foundDepths != _connectionBlockedDepthPerLayer.end()) {
            effectiveDepths = foundDepths->second;
        }
        auto foundLanding = _uLandingDepthPerLayer.find(layer.get_name());
        if (foundLanding != _uLandingDepthPerLayer.end()) {
            effectiveDepths.first = std::max(effectiveDepths.first, foundLanding->second.first);
            effectiveDepths.second = std::max(effectiveDepths.second, foundLanding->second.second);
        }
        if (effectiveDepths.first <= 1e-12 && effectiveDepths.second <= 1e-12) {
            continue;  // untouched layer: leave delimit's centring (fills the full height)
        }
        size_t windingIndex = get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
        double wirePitch = (turnAxis == 1) ? wires[windingIndex].get_maximum_outer_height()
                                           : wires[windingIndex].get_maximum_outer_width();

        // This layer's turns, in current (wound) order along the turn axis.
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
            return turns[a].get_coordinates()[turnAxis] < turns[b].get_coordinates()[turnAxis];
        });

        // The usable span runs from each window edge to the deepest run that crosses it — the
        // CONTINUOUS depth recorded alongside the slot counts, not the whole-slot quantization.
        // The ceil'd slots are the capacity currency (a partially covered slot still costs a whole
        // turn of capacity), but placing turns on the slot grid parked up to one pitch of dead band
        // against every run row: on 23_llc the primary's 1.093 mm run stack over the 0.679 mm
        // secondary rounded to 2 slots = 1.358 mm, floating the secondary 0.265 mm short of the
        // space it was actually free to reach. The turns may hug the runs exactly — edgeDepth
        // already includes the inter-winding insulation.
        const auto& blockedDepths = effectiveDepths;
        double spanLow = roundFloat(windowLowSide + blockedDepths.second, 9);
        double spanHigh = roundFloat(windowHighSide - blockedDepths.first, 9);
        size_t numberTurnsInLayer = layerTurns.size();
        int64_t bundleSize = get_layer_bundle_size(layer);
        int64_t numberBundles = (int64_t(numberTurnsInLayer) + bundleSize - 1) / bundleSize;
        // ABT #683 (Alf): a U layer that is NOT full must be packed CONTIGUOUSLY from the end
        // the wire arrives at — never fence-post spread across the band. Spreading a partial
        // landing layer parks its bundles at BOTH ends with a hole in the middle, and the
        // conductor then has to fly across that hole: on the 8t x 2p fixture layer 2's two
        // bundles sat at -2.10/-1.14 and 3.06/4.02, a 4.2 mm axial jump that MVB++ reads —
        // correctly — as a Z-style return, lays a dragback for, and makes every layer outside it
        // ride over with a bump. A U winding has no dragbacks and needs no bumps: the turns must
        // sit as close in the turn axis as the layers are in the layer axis. A FULL layer is
        // unaffected (spread with no slack IS packed).
        bool packFromArrival = false;
        bool arrivalAtHighSide = false;
        if (layer.get_section() && get_winding_order(layer.get_section().value()) == WindingOrder::U) {
            auto landingIt = _uLandingDepthPerLayer.find(layer.get_name());
            if (landingIt != _uLandingDepthPerLayer.end()) {
                if (landingIt->second.first > 1e-12) {
                    packFromArrival = true;
                    arrivalAtHighSide = true;
                }
                else if (landingIt->second.second > 1e-12) {
                    packFromArrival = true;
                    arrivalAtHighSide = false;
                }
            }
        }
        std::vector<double> stations;
        if (packFromArrival) {
            double firstCentre = arrivalAtHighSide ? spanHigh - wirePitch / 2
                                                   : spanLow + wirePitch / 2;
            double step = arrivalAtHighSide ? -wirePitch : wirePitch;
            for (size_t k = 0; k < numberTurnsInLayer; ++k) {
                stations.push_back(roundFloat(firstCentre + double(k) * step, 9));
            }
            // The assignment below walks layerTurns sorted along the axis, so hand it ascending
            // stations; packing from the high side produces them in descending order.
            if (arrivalAtHighSide) {
                std::reverse(stations.begin(), stations.end());
            }
        }
        else if (numberBundles == 1) {
            // A lone bundle has no fence-post meaning (no gaps to distribute), and centring it
            // strands it mid-window between the run stacks. Pack it against the LESS-blocked edge
            // instead — directly under the shallower run stack (on 23_llc: right below the primary's
            // two rows, not floating between them and the secondary's four entrance rows). Its own
            // terminal leads then exit at rows adjacent to that shallow stack, away from the deep one.
            bool packAgainstHighSide = blockedDepths.first <= blockedDepths.second;
            double firstCentre = packAgainstHighSide
                ? spanHigh - double(numberTurnsInLayer) * wirePitch + wirePitch / 2
                : spanLow + wirePitch / 2;
            for (size_t k = 0; k < numberTurnsInLayer; ++k) {
                stations.push_back(roundFloat(firstCentre + double(k) * wirePitch, 9));
            }
        }
        else {
            // Re-spread on the SAME rule the winder uses (ABT #578/#579): fence-post over the span,
            // bundles touching. Spacing the turns uniformly here would silently undo the winder's
            // bundle grouping on exactly the layers that carry connection leads — the real-winding
            // layers this whole routine exists for.
            stations = compute_spread_turn_stations((spanLow + spanHigh) / 2,
                                                    spanHigh - spanLow,
                                                    wirePitch,
                                                    int64_t(numberTurnsInLayer),
                                                    bundleSize);
        }
        // ABT #624 (Alf, 26_psps: "why is the turn placed so high when it has space below?").
        // The reserved band is a SOFT constraint — it keeps turns clear of connection rows —
        // while the winding window is HARD. When a layer's copper does not fit the span left
        // by its reservations (26_e3216: the U-landing depth was capped against 34 turns and
        // the layer ended up with 35), the fence-post spread centres the copper on the span
        // and the surplus leaves the window at BOTH ends. Give back reservation instead:
        // slide the whole layer along its axis until its copper is inside the window. If it
        // cannot fit even then, leave it — are_turns_inside_winding_window() will refuse the
        // wind rather than let it ship. Single-window coils only: this frame is window 0's.
        const bool singleWindingWindow =
            bobbin.get_processed_description()->get_winding_windows().size() == 1;
        if (!stations.empty() && singleWindingWindow) {
            const double lowEdge = stations.front() - wirePitch / 2;
            const double highEdge = stations.back() + wirePitch / 2;
            double shift = 0.0;
            // ABT #682: when the copper does not fit the span its reservations leave, the
            // fence-post spread centres it — which eats into BOTH reservations at once. Surrender
            // the DEEPER one and keep the shallower intact instead: the shallow depth is typically
            // a terminal lead ROW, real copper routed along that edge, while the deep one is a
            // landing/placement band with no conductor in it. Centring put a U layer's last turn
            // exactly on the entrance lead's row (E16, 3 mm margin: turn 12 at -4.67 mm, the lead
            // at -4.67 mm), which the 3D conductor builder refused as a collision — U alone and
            // margin alone both built, only the two together failed. Same "pack against the
            // less-blocked edge" rule the lone-bundle branch above already follows.
            if ((highEdge - lowEdge) > (spanHigh - spanLow) + 1e-12) {
                shift = (blockedDepths.second <= blockedDepths.first) ? (spanLow - lowEdge)
                                                                     : (spanHigh - highEdge);
            }
            if (highEdge + shift > windowHighSide + 1e-12) {
                shift = windowHighSide - highEdge;
            }
            if (lowEdge + shift < windowLowSide - 1e-12) {
                shift = windowLowSide - lowEdge;   // clamping down would push it out the bottom
            }
            if (std::abs(shift) > 1e-12 && (highEdge - lowEdge) <= (windowHighSide - windowLowSide) + 1e-12) {
                for (auto& station : stations) {
                    station = roundFloat(station + shift, 9);
                }
            }
        }
        for (size_t k = 0; k < numberTurnsInLayer && k < stations.size(); ++k) {
            auto coords = turns[layerTurns[k]].get_coordinates();
            coords[turnAxis] = stations[k];
            turns[layerTurns[k]].set_coordinates(coords);
        }
        auto layerCoordinates = layer.get_coordinates();
        layerCoordinates[turnAxis] = roundFloat((stations.front() + stations.back()) / 2, 9);
        layer.set_coordinates(std::vector<double>{layerCoordinates[0], layerCoordinates[1], 0});
        // ABT #616, OVERLAPPING layers only: track the re-spread turns' envelope in the layer
        // rect (and rescale the area-ratio filling factor with the extent), so the mid-loop
        // fit gates and the crossing measurement see the geometry the turns actually occupy.
        // CONTIGUOUS layers keep their partition dims untouched: their width difference IS the
        // record of the room surrendered to leads (ABT #449's measurement relies on it).
        if (turnAxis == 1) {
            const double newExtent =
                roundFloat(std::abs(stations.back() - stations.front()) + wirePitch, 9);
            const double oldExtent = layer.get_dimensions()[turnAxis];
            if (newExtent > 0 && std::abs(newExtent - oldExtent) > 1e-12) {
                auto layerDimensions = layer.get_dimensions();
                layerDimensions[turnAxis] = newExtent;
                layer.set_dimensions(layerDimensions);
                if (layer.get_filling_factor()) {
                    layer.set_filling_factor(
                        roundFloat(layer.get_filling_factor().value() * oldExtent / newExtent, 6));
                }
            }
        }
    }
    set_layers_description(layers);
    set_turns_description(turns);

    // A section must contain its layers (ABT #616/#624). Packing a layer against its
    // unblocked edge moves copper past the rect the partition gave the section, and
    // delimit_and_compact cannot repair it — it runs BEFORE this pass and re-centres what
    // this pass deliberately offsets. Grow each section along the turn axis to cover its own
    // conduction layers, so the rects keep describing the geometry that is really there
    // (without it the section filling factors read over 1 on every design whose layers were
    // packed: 06, 11, 19, 23, 24 all reported "does not fit" while their copper fitted).
    if (get_sections_description()) {
        auto sections = get_sections_description().value();
        bool sectionsChanged = false;
        for (auto& section : sections) {
            if (section.get_type() != ElectricalType::CONDUCTION) {
                continue;
            }
            double low = std::numeric_limits<double>::max();
            double high = std::numeric_limits<double>::lowest();
            size_t axis = 1;
            bool any = false;
            for (const auto& layer : layers) {
                if (layer.get_type() != ElectricalType::CONDUCTION || !layer.get_section()
                    || layer.get_section().value() != section.get_name()) {
                    continue;
                }
                axis = (layer.get_orientation() == WindingOrientation::OVERLAPPING) ? 1 : 0;
                low = std::min(low, layer.get_coordinates()[axis] - layer.get_dimensions()[axis] / 2);
                high = std::max(high, layer.get_coordinates()[axis] + layer.get_dimensions()[axis] / 2);
                any = true;
            }
            if (!any) {
                continue;
            }
            const double currentLow = section.get_coordinates()[axis] - section.get_dimensions()[axis] / 2;
            const double currentHigh = section.get_coordinates()[axis] + section.get_dimensions()[axis] / 2;
            const double newLow = std::min(currentLow, low);
            const double newHigh = std::max(currentHigh, high);
            if (newLow < currentLow - 1e-12 || newHigh > currentHigh + 1e-12) {
                auto coordinates = section.get_coordinates();
                auto dimensions = section.get_dimensions();
                const double oldExtent = dimensions[axis];
                coordinates[axis] = roundFloat((newLow + newHigh) / 2, 9);
                dimensions[axis] = roundFloat(newHigh - newLow, 9);
                section.set_coordinates(coordinates);
                section.set_dimensions(dimensions);
                if (section.get_filling_factor() && dimensions[axis] > 0) {
                    section.set_filling_factor(roundFloat(
                        section.get_filling_factor().value() * oldExtent / dimensions[axis], 6));
                }
                sectionsChanged = true;
            }
        }
        if (sectionsChanged) {
            set_sections_description(sections);
        }
    }
}

// ABT #624: is every turn's copper inside the winding window? Reported as part of the wind's
// FINAL verdict only — deliberately NOT inside are_sections_and_layers_fitting(), which the
// fixpoint consults mid-loop: rejecting a transient state there sends the winder down
// try_rewind() and it settles on a layout whose terminal routes cross copper (measured on
// 13_current_sense, which then violates the ABT #577 clearance contract). The trajectory must
// stay exactly as it was; what changes is that a coil whose copper ends up outside its window
// no longer reports success.
bool Coil::are_turns_inside_winding_window() {
    if (!get_turns_description()) {
        return true;
    }
    auto bobbin = resolve_bobbin();
    if (bobbin.get_winding_window_shape() != WindingWindowShape::RECTANGULAR) {
        return true;   // round windows block angularly (ABT #187), not by this envelope
    }
    // MAS getters return BY VALUE: binding a reference through
    // get_processed_description()->get_winding_windows()[0] dangles into a destroyed temporary,
    // which is exactly how the first version of this check silently read garbage and never
    // fired at all. Copy first.
    auto processedDescription = bobbin.get_processed_description();
    if (!processedDescription || processedDescription->get_winding_windows().empty()) {
        return true;
    }
    // MULTI-COLUMN: a coil can have several winding windows and a turn belongs to whichever
    // one its group was wound in (apply_group_window_sides mirrors them into place), so the
    // test is "inside ANY window", never "inside window 0".
    auto windingWindows = processedDescription->get_winding_windows();
    struct WindowBox { double x0, x1, y0, y1; };
    std::vector<WindowBox> boxes;
    for (const auto& windingWindow : windingWindows) {
        if (!windingWindow.get_coordinates() || !windingWindow.get_width() || !windingWindow.get_height()) {
            continue;
        }
        boxes.push_back({(*windingWindow.get_coordinates())[0] - *windingWindow.get_width() / 2,
                         (*windingWindow.get_coordinates())[0] + *windingWindow.get_width() / 2,
                         (*windingWindow.get_coordinates())[1] - *windingWindow.get_height() / 2,
                         (*windingWindow.get_coordinates())[1] + *windingWindow.get_height() / 2});
    }
    if (boxes.empty()) {
        return true;
    }
    const double tolerance = 1e-9;
    auto wires = get_wires();
    auto turnsToCheck = get_turns_description().value();
    for (const auto& turn : turnsToCheck) {
        const size_t windingIndex = get_winding_index_by_name(turn.get_winding());
        const double halfWidth = wires[windingIndex].get_maximum_outer_width() / 2;
        const double halfHeight = wires[windingIndex].get_maximum_outer_height() / 2;
        const auto& coordinates = turn.get_coordinates();
        bool insideAny = false;
        for (const auto& box : boxes) {
            if (coordinates[0] - halfWidth >= box.x0 - tolerance &&
                coordinates[0] + halfWidth <= box.x1 + tolerance &&
                coordinates[1] - halfHeight >= box.y0 - tolerance &&
                coordinates[1] + halfHeight <= box.y1 + tolerance) {
                insideAny = true;
                break;
            }
        }
        if (!insideAny) {
            if (std::getenv("MKF_BLOCKING_DIAG")) {
                std::cerr << "[window] turn " << turn.get_name() << " at (" << coordinates[0]
                          << "," << coordinates[1] << ") lies outside every winding window\n";
            }
            return false;
        }
    }
    return true;
}

std::map<std::string, uint64_t> Coil::align_blocked_ring_turns() {
    std::map<std::string, uint64_t> ringDeficitSlots;
    if (!get_turns_description() || !get_layers_description()) {
        return ringDeficitSlots;
    }
    auto bobbin = resolve_bobbin();
    if (bobbin.get_winding_window_shape() != WindingWindowShape::ROUND) {
        return ringDeficitSlots;
    }
    auto wires = get_wires();

    // Signed smallest angular difference a - b, in degrees, in (-180, 180].
    auto angularDifference = [](double a, double b) {
        return std::fmod(a - b + 540.0, 360.0) - 180.0;
    };
    auto normalizeAngle = [](double a) {
        double n = std::fmod(a, 360.0);
        return n < 0 ? n + 360.0 : n;
    };

    // Displacement-only fixpoint: re-spreading a ring moves its end turns, which moves the leads
    // attached to them (spaces are recomputed from the turns each pass), which moves the corridors
    // on the rings THOSE leads cross. Converges in a couple of passes for realistic windings —
    // but the ABT #723 sector cases chase sub-degree lead movements (a re-spread ring moves its
    // own exit lead's corridor onto a neighbouring turn), so give the chase more budget: extra
    // iterations only run while intrusions persist.
    const size_t maximumIterations = 24;
    // ABT #723: a ring that needs re-spreading over and over is not converging by
    // displacement — its free arc is arithmetically sufficient but no stable discrete
    // arrangement realizes it against the moving lead corridors (the bifilar full-circle
    // limit cycle). After a few re-spreads, escalate it to a capacity deficit of one slot:
    // the blocking re-wind spills a turn to the next ring, which genuinely frees arc and
    // terminates the cycle.
    std::map<std::string, size_t> respreadCountPerRing;
    // Owner ruling (ABT #723, option "spill"): a ring that keeps needing re-spreads is
    // over-full under the no-turn-in-corridor rule; escalate to a capacity deficit so a
    // turn spills inward. The structural consequence (an extra ring adds one inter-ring
    // crossing per parallel) is accepted and the affected pins updated.
    const size_t maximumRespreadsPerRing = 4;
    for (size_t iteration = 0; iteration < maximumIterations; ++iteration) {
        auto turns = get_turns_description().value();
        auto spaces = get_connection_reserved_spaces();

        // Ring membership (turns are cartesian at this point) and mean ring radius.
        std::map<std::string, std::vector<size_t>> ringTurnIndexes;
        for (size_t t = 0; t < turns.size(); ++t) {
            if (turns[t].get_layer()) {
                ringTurnIndexes[turns[t].get_layer().value()].push_back(t);
            }
        }
        // A ring's angular territory is the arc its turns actually OCCUPY (largest-gap analysis —
        // the section's polar envelope can be stale after compaction, cf. ABT #186). A sector ring
        // (contiguous windings) must be re-spread WITHIN its own arc — spilling its turns into a
        // neighbour's sector would collide with that winding. Same rule as the marker filter in
        // toroidal_connection_reserved_spaces so the two stay in agreement.
        std::map<std::string, std::pair<bool, std::pair<double, double>>> ringOccupiedArc;  // ring -> {fullCircle, {start, span}}
        for (const auto& [ringName, turnIdxs] : ringTurnIndexes) {
            std::vector<double> angles;
            for (size_t t : turnIdxs) {
                angles.push_back(normalizeAngle(std::atan2(turns[t].get_coordinates()[1], turns[t].get_coordinates()[0]) * 180.0 / std::numbers::pi));
            }
            std::sort(angles.begin(), angles.end());
            double largestGap = 360.0 - (angles.back() - angles.front());
            double gapEnd = angles.front();
            for (size_t i = 1; i < angles.size(); ++i) {
                double gap = angles[i] - angles[i - 1];
                if (gap > largestGap) {
                    largestGap = gap;
                    gapEnd = angles[i];
                }
            }
            double meanPitch = 360.0 / double(std::max<size_t>(angles.size(), 1));
            bool fullCircle = largestGap <= 2 * meanPitch;
            // Raw arc of turn CENTRES (no margin): margins are added at the use site, clamped so
            // they never invade a neighbouring winding's sector at the same radius.
            double rawSpan = 360.0 - largestGap;
            ringOccupiedArc[ringName] = {fullCircle, {gapEnd, rawSpan}};
        }
        std::map<std::string, double> ringRadius;
        for (const auto& [ringName, turnIdxs] : ringTurnIndexes) {
            double sum = 0;
            for (size_t t : turnIdxs) {
                sum += std::hypot(turns[t].get_coordinates()[0], turns[t].get_coordinates()[1]);
            }
            ringRadius[ringName] = sum / double(turnIdxs.size());
        }

        // Blocked corridors for turn CENTERS per ring: marker azimuth +- (marker angular half-width
        // + the ring's own turn angular half-pitch). dimensions[1] is the lead's azimuthal wire
        // thickness, matching the drawn radial lead's convention.
        std::map<std::string, std::vector<std::pair<double, double>>> corridorsPerRing;  // {center, half}
        for (const auto& space : spaces) {
            if (space.layer.empty() || !ringRadius.count(space.layer)) {
                continue;
            }
            double radius = ringRadius.at(space.layer);
            double markerHalfAngle = wound_distance_to_angle(space.dimensions[1], radius) / 2;
            size_t ringWindingIndex = get_winding_index_by_name(turns[ringTurnIndexes.at(space.layer)[0]].get_winding());
            double turnHalfAngle = wound_distance_to_angle(wires[ringWindingIndex].get_maximum_outer_height(), radius) / 2;
            if (markerHalfAngle >= 180 || turnHalfAngle >= 180) {
                continue;  // wound_distance_to_angle's does-not-fit sentinel — nothing sane to block
            }
            corridorsPerRing[space.layer].push_back({normalizeAngle(space.rotation), markerHalfAngle + turnHalfAngle});
        }

        bool anyIntrusion = false;
        bool anyDisplacement = false;
        for (const auto& [ringName, corridors] : corridorsPerRing) {
            const auto& turnIdxs = ringTurnIndexes.at(ringName);
            double radius = ringRadius.at(ringName);
            size_t ringWindingIndex = get_winding_index_by_name(turns[turnIdxs[0]].get_winding());
            double turnPitchAngle = wound_distance_to_angle(wires[ringWindingIndex].get_maximum_outer_height(), radius);

            auto insideCorridor = [&](double angle) -> const std::pair<double, double>* {
                for (const auto& corridor : corridors) {
                    if (std::abs(angularDifference(angle, corridor.first)) < corridor.second - 1e-9) {
                        return &corridor;
                    }
                }
                return nullptr;
            };

            bool intrusion = false;
            for (size_t t : turnIdxs) {
                double turnAngle = normalizeAngle(std::atan2(turns[t].get_coordinates()[1], turns[t].get_coordinates()[0]) * 180.0 / std::numbers::pi);
                if (insideCorridor(turnAngle)) {
                    intrusion = true;
                    break;
                }
            }
            if (!intrusion) {
                continue;
            }
            anyIntrusion = true;
            if (respreadCountPerRing[ringName] >= maximumRespreadsPerRing) {
                // Escalate the limit cycle to capacity (see the counter's comment above).
                ringDeficitSlots[ringName] = std::max<uint64_t>(ringDeficitSlots[ringName], 1);
                continue;
            }

            // The ring's angular territory. Full-circle rings (overlapping windings) use the
            // whole circle cyclically. Sector rings (contiguous windings) use the space the
            // sector actually OWNS: the occupied arc extended outward until one full pitch
            // short of the nearest FOREIGN turn at the same radius on each side (the
            // neighbouring sector's boundary turns), capped at half the circle per side.
            //
            // ABT #723: the extension used to be capped at half a pitch per side, which made
            // the deficit measure SELF-REFERENTIAL — a blocking re-wind spills turns, the
            // occupied arc shrinks with them, the free space shrinks in proportion, and the
            // same corridor deficit re-fires forever (the observed 2-cycle oscillation
            // ratcheting blocked slots up to the iteration cap). With the territory bounded
            // by the real neighbours instead, spilled turns genuinely free arc and the
            // fixpoint converges. A ring with NO angular neighbour (sole sector at this
            // depth, e.g. a spilled third ring) extends until the neighbour clamp or the
            // half-circle cap, which the corridor sweep then prunes.
            bool fullCircle = true;
            double sectionSpan = 360;
            double territoryStart = 0;
            if (ringOccupiedArc.count(ringName) && !ringOccupiedArc.at(ringName).first) {
                fullCircle = false;
                double rawStart = ringOccupiedArc.at(ringName).second.first;
                double rawSpan = ringOccupiedArc.at(ringName).second.second;
                double rawEnd = normalizeAngle(rawStart + rawSpan);
                double wireRadialWidth = wires[ringWindingIndex].get_maximum_outer_width();
                double maximumExtension = std::max(turnPitchAngle / 2, (360.0 - rawSpan) / 2);
                double marginStart = maximumExtension;
                double marginEnd = maximumExtension;
                for (size_t t = 0; t < turns.size(); ++t) {
                    if (turns[t].get_layer() && turns[t].get_layer().value() == ringName) {
                        continue;
                    }
                    double foreignRadius = std::hypot(turns[t].get_coordinates()[0], turns[t].get_coordinates()[1]);
                    if (std::abs(foreignRadius - radius) > wireRadialWidth / 2) {
                        continue;  // different ring depth: no angular contention
                    }
                    double foreignAngle = normalizeAngle(std::atan2(turns[t].get_coordinates()[1], turns[t].get_coordinates()[0]) * 180.0 / std::numbers::pi);
                    double behindStart = normalizeAngle(rawStart - foreignAngle);
                    double aheadOfEnd = normalizeAngle(foreignAngle - rawEnd);
                    if (behindStart < 180.0) {
                        marginStart = std::min(marginStart, std::max(0.0, behindStart - turnPitchAngle));
                    }
                    if (aheadOfEnd < 180.0) {
                        marginEnd = std::min(marginEnd, std::max(0.0, aheadOfEnd - turnPitchAngle));
                    }
                }
                territoryStart = normalizeAngle(rawStart - marginStart);
                sectionSpan = rawSpan + marginStart + marginEnd;
            }

            // ABT #723: MINIMAL NUDGE first. The full even re-spread below moves every turn of
            // the ring — including the lead-attached END turns, whose corridors then land on
            // OTHER rings' turns, and the displacement ping-pongs (an observed 24-iteration
            // limit cycle on the bifilar overlapping case). Nudging only the INTRUDING turns to
            // the nearest clear corridor edge leaves the leads where they are, so the corridor
            // set stays fixed and the pass converges. Only when a nudge cannot fit (no clear
            // edge with a full pitch to every same-ring neighbour inside the territory) does
            // the ring fall back to the even re-spread / capacity-deficit machinery.
            {
                auto moveTurnToAngle = [&](size_t t, double newAngleDegrees) {
                    double turnRadius = std::hypot(turns[t].get_coordinates()[0], turns[t].get_coordinates()[1]);
                    double newAngleRadians = newAngleDegrees / 180.0 * std::numbers::pi;
                    turns[t].set_coordinates(std::vector<double>{
                        roundFloat(turnRadius * std::cos(newAngleRadians), 9),
                        roundFloat(turnRadius * std::sin(newAngleRadians), 9), 0});
                    if (turns[t].get_additional_coordinates()) {
                        auto additionalCoordinates = turns[t].get_additional_coordinates().value();
                        for (auto& additional : additionalCoordinates) {
                            double additionalRadius = std::hypot(additional[0], additional[1]);
                            additional = {roundFloat(additionalRadius * std::cos(newAngleRadians), 9),
                                          roundFloat(additionalRadius * std::sin(newAngleRadians), 9)};
                        }
                        turns[t].set_additional_coordinates(additionalCoordinates);
                    }
                };
                auto insideTerritory = [&](double angle) {
                    if (fullCircle) {
                        return true;
                    }
                    return normalizeAngle(angle - territoryStart) <= sectionSpan + 1e-9;
                };
                bool nudgedAll = true;
                for (size_t t : turnIdxs) {
                    double turnAngle = normalizeAngle(std::atan2(turns[t].get_coordinates()[1], turns[t].get_coordinates()[0]) * 180.0 / std::numbers::pi);
                    const auto* corridor = insideCorridor(turnAngle);
                    if (corridor == nullptr) {
                        continue;
                    }
                    double lowerEdge = normalizeAngle(corridor->first - corridor->second - 1e-6);
                    double upperEdge = normalizeAngle(corridor->first + corridor->second + 1e-6);
                    bool lowerCloser = std::abs(angularDifference(turnAngle, lowerEdge)) <= std::abs(angularDifference(turnAngle, upperEdge));
                    bool placed = false;
                    for (double candidate : {lowerCloser ? lowerEdge : upperEdge, lowerCloser ? upperEdge : lowerEdge}) {
                        if (!insideTerritory(candidate) || insideCorridor(candidate) != nullptr) {
                            continue;
                        }
                        bool clearOfNeighbours = true;
                        for (size_t other : turnIdxs) {
                            if (other == t) {
                                continue;
                            }
                            double otherAngle = normalizeAngle(std::atan2(turns[other].get_coordinates()[1], turns[other].get_coordinates()[0]) * 180.0 / std::numbers::pi);
                            if (std::abs(angularDifference(candidate, otherAngle)) < turnPitchAngle - 1e-9) {
                                clearOfNeighbours = false;
                                break;
                            }
                        }
                        if (clearOfNeighbours) {
                            moveTurnToAngle(t, candidate);
                            placed = true;
                            break;
                        }
                    }
                    if (!placed) {
                        nudgedAll = false;
                        break;
                    }
                }
                if (nudgedAll) {
                    anyDisplacement = true;
                    continue;   // ring settled by nudges alone; leads untouched
                }
            }

            // Free angular space = territory minus the corridors (union measured by sweeping — the
            // stacked parallels' corridors overlap heavily).
            double blockedUnion = 0;
            {
                const int sweepSteps = 3600;
                int blockedSteps = 0;
                for (int s = 0; s < sweepSteps; ++s) {
                    double sweepAngle = normalizeAngle(territoryStart + s * sectionSpan / sweepSteps);
                    if (insideCorridor(sweepAngle)) {
                        blockedSteps++;
                    }
                }
                blockedUnion = blockedSteps * sectionSpan / sweepSteps;
            }
            double totalFree = sectionSpan - blockedUnion;
            double needed = double(turnIdxs.size()) * turnPitchAngle;
            if (totalFree < needed - 1e-9) {
                // The ring is too full to clear its corridors by displacement alone: report the
                // CAPACITY deficit (in this ring's own turn slots) so wind()'s toroidal blocking
                // loop can re-wind with that many slots reserved (turns spill to the next ring),
                // then displacement succeeds on the re-wound geometry.
                uint64_t deficit = uint64_t(std::ceil((needed - totalFree) / turnPitchAngle - 1e-9));
                ringDeficitSlots[ringName] = std::max(ringDeficitSlots[ringName], deficit);
                logEntry("Toroidal ring '" + ringName + "' needs " + std::to_string(deficit)
                             + " blocked slot(s) to clear its connection corridors ("
                             + std::to_string(totalFree) + " deg free < " + std::to_string(needed) + " deg needed)",
                         "Coil", 2);
                continue;
            }

            // Even spread over the free space: march from the territory start (for a full circle,
            // from the end of the widest corridor), placing each turn (original cyclic order
            // preserved) every totalFree/n of FREE arc, jumping over corridors as they are met.
            double referenceAngle = territoryStart;
            if (fullCircle) {
                const std::pair<double, double>* widest = &corridors[0];
                for (const auto& corridor : corridors) {
                    if (corridor.second > widest->second) {
                        widest = &corridor;
                    }
                }
                referenceAngle = normalizeAngle(widest->first + widest->second);
            }
            double spacing = totalFree / double(turnIdxs.size());

            // Ring turns in cyclic order starting just after the reference angle.
            std::vector<size_t> orderedTurns(turnIdxs.begin(), turnIdxs.end());
            std::sort(orderedTurns.begin(), orderedTurns.end(), [&](size_t a, size_t b) {
                double angleA = normalizeAngle(normalizeAngle(std::atan2(turns[a].get_coordinates()[1], turns[a].get_coordinates()[0]) * 180.0 / std::numbers::pi) - referenceAngle);
                double angleB = normalizeAngle(normalizeAngle(std::atan2(turns[b].get_coordinates()[1], turns[b].get_coordinates()[0]) * 180.0 / std::numbers::pi) - referenceAngle);
                return angleA < angleB;
            });

            double position = referenceAngle;
            for (size_t k = 0; k < orderedTurns.size(); ++k) {
                // Advance by the next chunk of free arc, jumping over any corridor met on the way.
                double toAdvance = (k == 0) ? spacing / 2 : spacing;
                // Safety bound: the free-space guard above should make this unreachable, but a
                // sweep-resolution edge case must degrade to a partial march, not an infinite loop.
                size_t marchGuard = 0;
                while (toAdvance > 1e-12 && marchGuard++ < 1000) {
                    const auto* corridor = insideCorridor(position + 1e-9);
                    if (corridor != nullptr) {
                        position = normalizeAngle(corridor->first + corridor->second);
                        continue;
                    }
                    // Distance to the nearest corridor start ahead of us (cyclically).
                    double nearestAhead = 360.0;
                    for (const auto& c : corridors) {
                        double ahead = normalizeAngle((c.first - c.second) - position);
                        if (ahead > 1e-12) {
                            nearestAhead = std::min(nearestAhead, ahead);
                        }
                    }
                    double step = std::min(toAdvance, nearestAhead);
                    position = normalizeAngle(position + step);
                    toAdvance -= step;
                }
                size_t t = orderedTurns[k];
                double turnRadius = std::hypot(turns[t].get_coordinates()[0], turns[t].get_coordinates()[1]);
                double positionRadians = position / 180.0 * std::numbers::pi;
                turns[t].set_coordinates(std::vector<double>{
                    roundFloat(turnRadius * std::cos(positionRadians), 9),
                    roundFloat(turnRadius * std::sin(positionRadians), 9)});
                // Same invariant as ABT #186: a toroidal turn's rotation is its azimuth.
                turns[t].set_rotation(roundFloat(position, 6));
                if (turns[t].get_additional_coordinates()) {
                    auto additionalCoordinates = turns[t].get_additional_coordinates().value();
                    for (auto& additional : additionalCoordinates) {
                        double additionalRadius = std::hypot(additional[0], additional[1]);
                        additional = {roundFloat(additionalRadius * std::cos(positionRadians), 9),
                                      roundFloat(additionalRadius * std::sin(positionRadians), 9)};
                    }
                    turns[t].set_additional_coordinates(additionalCoordinates);
                }
            }
            anyDisplacement = true;
            respreadCountPerRing[ringName]++;
        }

        if (anyDisplacement) {
            set_turns_description(turns);
        }
        if (std::getenv("MKF_BLOCKING_DIAG")) {
            std::cerr << "[align] iter " << iteration << " rings=" << corridorsPerRing.size()
                      << " intrusion=" << anyIntrusion << " displaced=" << anyDisplacement
                      << " deficits=" << ringDeficitSlots.size() << "\n";
        }
        if (!anyIntrusion || !anyDisplacement) {
            break;  // clean, or stuck on capacity (deficits reported to the caller) — either way stop
        }
        ringDeficitSlots.clear();  // re-measured next pass on the displaced geometry
    }
    return ringDeficitSlots;
}

void Coil::apply_connection_reserved_space() {
    if (!get_sections_description() || !get_layers_description()) {
        return;
    }
    // ABT #187: this function's height/area math is rectangular. Toroidal (round-window) coils
    // historically hit it as a no-op because their spaces carried no layer; now that toroidal
    // crossing markers DO name their ring (angular corridors, consumed by align_blocked_ring_turns),
    // gate explicitly — mixing a ring's polar dimensions with a marker's cartesian height would
    // corrupt the filling factors.
    if (resolve_bobbin().get_winding_window_shape() == WindingWindowShape::ROUND) {
        return;
    }
    auto spaces = get_connection_reserved_spaces();
    if (spaces.empty()) {
        return;
    }

    auto layers = get_layers_description().value();

    // Space occupied by connection leads on each conduction layer, measured along that layer's TURN
    // axis — the direction its turns stack, and so the axis a lead takes a slot out of. Every reserved
    // space that names a layer (an inter-layer transition, or a terminal lead passing over that layer)
    // occupies one wire diameter of it. Free-space terminal segments (no layer) are drawn but reserve
    // no layer space. Marker dimensions are plain {X, Y} like every other rectangle in the coil, so
    // the turn axis is index 1 for an OVERLAPPING layer and index 0 for a CONTIGUOUS one.
    std::map<std::string, size_t> layerTurnAxis;
    for (const auto& layer : layers) {
        layerTurnAxis[layer.get_name()] = (layer.get_orientation() == WindingOrientation::OVERLAPPING) ? 1 : 0;
    }
    std::map<std::string, double> turnAxisReservedPerLayer;
    for (const auto& space : spaces) {
        auto turnAxis = layerTurnAxis.find(space.layer);
        if (space.layer.empty() || turnAxis == layerTurnAxis.end()) {
            continue;
        }
        turnAxisReservedPerLayer[space.layer] += space.dimensions[turnAxis->second];
    }

    // Charge each affected layer for the space its leads take ALONG ITS TURN AXIS. A resulting value
    // above 1 means the leads no longer fit alongside the turns (the layer is over-subscribed and the
    // build needs more space).
    //
    // ABT #424: which of the layer's two dimensions that axis is depends on the layer's orientation,
    // and the two cases are mirror images (see wind_by_layers):
    //     OVERLAPPING: turns stack along the layer's HEIGHT (its width is one wire) -> charge height
    //     CONTIGUOUS:  turns run   along the layer's WIDTH  (its height is one wire) -> charge width
    // Charging the height unconditionally made a CONTIGUOUS layer — which is exactly one wire tall by
    // construction — surrender its whole thickness for a single lead, and the section picked up a whole
    // layer's worth of reserved area instead of one turn slot.
    std::map<std::string, double> reservedAreaPerSection;
    for (auto& layer : layers) {
        auto it = turnAxisReservedPerLayer.find(layer.get_name());
        if (it == turnAxisReservedPerLayer.end()) {
            continue;
        }
        bool turnsStackAlongHeight = (layer.get_orientation() == WindingOrientation::OVERLAPPING);
        double turnAxisExtent = turnsStackAlongHeight ? layer.get_dimensions()[1] : layer.get_dimensions()[0];
        double layerAxisExtent = turnsStackAlongHeight ? layer.get_dimensions()[0] : layer.get_dimensions()[1];
        if (turnAxisExtent <= 0) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR, "Non-positive layer extent along its turn axis while applying connection reserved space to layer " + layer.get_name());
        }
        if (!layer.get_filling_factor()) {
            throw CoilException(ErrorCode::COIL_WINDING_ERROR, "Layer filling factor not set before applying connection reserved space to layer " + layer.get_name());
        }
        double reserved = it->second;
        // ABT #430: charge only the room the layer did NOT already surrender to these leads. When real
        // winding geometry blocks turn slots, wind_by_rectangular_layers shrinks the layer by exactly
        // the leads' room and computes its filling factor against what remains ("shrink the layer
        // height by the blocked slots ... leaving room for the leads") — so the extent below ALREADY
        // excludes them, and charging the full lead extent again counted the same room twice. It stayed
        // invisible while the leads were thin relative to the layer and exploded with fine wire and deep
        // lead stacks: on 13_current_sense_er95_n87 a correctly-packed 0.95-full layer was reported at
        // 2.97, and the coil read as not-fitting when it fits.
        //
        // The remainder is real: the shrink is capped at one turn slot minimum, so a layer asked for
        // more room than it can give keeps the difference as a genuine "the leads do not fit here"
        // signal. Layers with no blocking applied — CONTIGUOUS ones, which have no turn blocking at all
        // (ABT #427) — surrender nothing, so they are charged in full, as they must be.
        double reservedNotYetMadeRoomFor = reserved;
        auto surrendered = _connectionBlockedRoomPerLayer.find(layer.get_name());
        if (surrendered != _connectionBlockedRoomPerLayer.end()) {
            // Never below zero: blocking rounds up to whole turn slots, so the room given up can
            // slightly exceed what the leads need, and negative leftover space is meaningless.
            reservedNotYetMadeRoomFor = std::max(reserved - surrendered->second, 0.0);
        }
        // ABT #616: a layer with CONTINUOUS blocked depths applied was aligned off the rows and
        // compacted to its turn envelope — the leads' room already lies wholly outside its
        // extent, so any further charge counts that room twice (delimit's envelope ff is ~1.0
        // for a full layer, and the double charge pushed it past the fitting threshold).
        if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
            auto blockedDepths = _connectionBlockedDepthPerLayer.find(layer.get_name());
            if (blockedDepths != _connectionBlockedDepthPerLayer.end()
                && (blockedDepths->second.first > 1e-12 || blockedDepths->second.second > 1e-12)) {
                reservedNotYetMadeRoomFor = 0.0;
            }
        }
        // Add the leads' own share of the layer to the turns' share — the same shape as the section's
        // factor below (filling factor + reserved area / section area). Since the filling factor is an
        // area ratio, the leads' share is (reserved * layerAxisExtent) / (turnAxisExtent *
        // layerAxisExtent), i.e. reserved / turnAxisExtent: the layer-axis extent is one wire on both
        // sides and cancels.
        //
        // NOT the leftover space in the denominator. `fill * extent / (extent - reserved)` is singular
        // when the leads take the whole layer and NEGATIVE beyond it, and a negative filling factor
        // reads as FITTING at are_sections_and_layers_fitting's `> 1` test — so the layers whose leads
        // need MORE room than the layer has, the worst ones, would have reported as fitting. That sign
        // flip is what the old `std::max(extent - reserved, extent * 0.01)` clamp was holding shut, and
        // it paid for it by saturating every overflow at exactly 100x the pre-lead factor: a layer 1%
        // short and one 5x over-subscribed reported the same number. It also masked ABT #424 for two
        // months — on a contiguous layer reserved == extent exactly, so without the clamp that bug
        // would have produced an inf on day one instead of a plausible-looking finite number.
        //
        // This form is linear and increasing in the reserved space, always positive, never singular,
        // and crosses 1 at the IDENTICAL point as the old one — A/(W(E-r)) >= 1 <=> A + rW >= WE —
        // so no fitting verdict anywhere changes, only the magnitude reported past the crossing.
        layer.set_filling_factor(layer.get_filling_factor().value() + reservedNotYetMadeRoomFor / turnAxisExtent);
        if (layer.get_section()) {
            // The lead's own footprint on this layer: the slot it takes along the turn axis times the
            // layer's extent along the other axis (one wire either way, whichever axis that is). The
            // SECTION is charged the FULL lead extent even when the layer surrendered room for it —
            // the section's own dimensions never shrank, so that room is still area the leads occupy
            // inside it. Only the layer, whose extent was reduced, must not be charged twice.
            reservedAreaPerSection[layer.get_section().value()] += reserved * layerAxisExtent;
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
    _groupWindowSidesApplied = false;
    set_sections_description(std::nullopt);
    set_layers_description(std::nullopt);
    set_turns_description(std::nullopt);
    return true;
}

bool Coil::rewind_layers_and_turns() {
    if (!get_sections_description()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
            "rewind_layers_and_turns needs a sections description to re-flow into");
    }
    // Sections may arrive at their FINAL multi-window positions (hand-edited or
    // deserialized); the layer/turn placement math lives in the +x winding
    // frame, so unwrap, re-flow, and re-apply — the custom rectangles round-trip
    // through the same transform.
    bool rewrapGroupWindowSides = _groupWindowSidesApplied;
    if (rewrapGroupWindowSides) {
        apply_group_window_sides(true);
    }
    wind_by_layers();
    bool result = false;
    if (get_layers_description()) {
        wind_by_turns();
        result = get_turns_description().has_value();
    }
    if (rewrapGroupWindowSides) {
        apply_group_window_sides(false);
    }
    if (result) {
        generate_toroidal_additional_coordinates();
    }
    return result;
}

void Coil::generate_toroidal_additional_coordinates() {
    // The outer return crossings of toroidal turns (additionalCoordinates) are
    // generated inside delimit_and_compact_round_window; re-flow paths that skip
    // compaction (custom-rect rewinds, compact-off winds) must rebuild them or
    // the external half of every turn silently disappears from the description.
    if (!settings.get_coil_include_additional_coordinates()) {
        return;
    }
    if (!get_turns_description() || !get_layers_description()) {
        return;
    }
    auto bobbin = resolve_bobbin();
    if (bobbin.get_winding_window_shape() != WindingWindowShape::ROUND) {
        return;
    }
    // Same frame dance as delimit_and_compact_round_window: the additional-turn
    // math runs on polar turns, the stored description is cartesian.
    convert_turns_to_polar_coordinates();
    wind_toroidal_additional_turns();
    convert_turns_to_cartesian_coordinates();
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
    for (const auto& section : sections) {
        size_t windingIndex = get_winding_index_by_name(section.get_partial_windings()[0].get_winding());
        stackUp.push_back(windingIndex);
    }
    return stackUp;
}

bool Coil::wind(std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions) {
    // REAL WINDING: a wire that makes N turns crosses the winding-window plane N+1
    // times — the beginning of the first turn occupies its own physical slot in the
    // cross-section (for 3 turns, 4 wire crossings per parallel appear in the 2D
    // projection). Wind one extra turn per winding (one extra slot per parallel, since
    // slots = numberTurns x numberParallels): the first placed turn of each parallel is
    // that beginning crossing; every following turn is the wrap ending at its own
    // crossing. The bump lasts for the whole wind — including the turn-blocking
    // re-winds — and the electrical turn count is restored on exit (RAII), so
    // inductance and turns-ratio semantics are untouched.
    struct RealWindingCrossingBump {
        Coil& coil;
        bool active = false;
        // ABT #728: the destructor's station-zeroing may only touch turns produced by THIS
        // wind. wind() clears turns_description before arming the bump and flips this flag at
        // that clear — so when a wind fails before producing turns, a PREVIOUS wind's intact
        // turns are never mis-zeroed (DC resistance silently shrank by one turn per parallel).
        // Deliberately NOT gated on the wind's fitting result: with windEvenIfNotFit the
        // caller consumes the not-fitting turns and still needs the station carrying no length.
        bool turnsAreFresh = false;
        explicit RealWindingCrossingBump(Coil& c) : coil(c) {}
        void arm() {
            if (active || !settings.get_coil_use_real_winding_geometry()) return;
            for (auto& winding : coil.get_mutable_functional_description()) {
                winding.set_number_turns(winding.get_number_turns() + 1);
            }
            active = true;
        }
        ~RealWindingCrossingBump() {
            if (active) {
                for (auto& winding : coil.get_mutable_functional_description()) {
                    winding.set_number_turns(winding.get_number_turns() - 1);
                }
                // ABT #674: the extra crossing is a STATION, not a turn. Each parallel's first
                // station is where its wire begins; the copper is the WRAP between consecutive
                // stations, so N turns need N+1 stations and every wrap's length belongs to the
                // station it ENDS at. Leaving a full turn's length on the beginning crossing made
                // every consumer that sums turn lengths count one turn too many —
                // WindingOhmicLosses does exactly that, so DC resistance came out ~5% high on a
                // 20-turn winding (1/N) whenever real winding was on. The station stays (MVB++
                // builds one wrap fewer without it: measured 82 wrap primitives against 87); only
                // its length goes, which is what makes the sum the conductor's real length.
                if (turnsAreFresh && coil.get_turns_description()) {
                    auto crossingTurns = coil.get_turns_description().value();
                    // ABT #728: the station is tagged EXPLICITLY — the winders name every turn
                    // "<winding> parallel <p> turn <k>" with k the per-parallel wind-order
                    // counter, so each parallel's beginning crossing is exactly its "turn 0".
                    // First-seen vector order is not a marker: any post-wind pass that reorders
                    // the description would zero a mid-winding wrap instead.
                    for (auto& crossingTurn : crossingTurns) {
                        if (crossingTurn.get_name() == crossingTurn.get_winding() + " parallel "
                                                          + std::to_string(crossingTurn.get_parallel())
                                                          + " turn 0") {
                            crossingTurn.set_length(0);
                        }
                    }
                    coil.set_turns_description(crossingTurns);
                }
            }
        }
    } realWindingCrossingBump(*this);

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

    // ABT #676: margin tape is PERSISTED on the sections but was only ever read back from the
    // transient _marginsPerSection, which a caller arriving with an already-wound coil does not
    // have. Every consumer that re-winds — MVB++'s internal autocomplete behind the 3D view and
    // the STEP export, PyOpenMagnetics, anything round-tripping a MAS — therefore dropped the
    // margins and wound the copper over the tape, which is why the 3D ignored a margin the 2D
    // drew correctly. The wound coil already carries the answer; read it back rather than
    // require it to be handed in again.
    //
    // ABT #724 hardening: an EXPLICIT reset_margins_per_section() means "the next wind is
    // margin-free" (CoilAdviser candidate sweeps rely on it) — the empty vector it leaves
    // behind must not be taken as "nothing was handed in, recover the old ones", or margins
    // could never actually be cleared. And when margin tape is disallowed outright, there is
    // nothing legitimate to recover. (The remaining #724 defect — positional index-matching
    // of recovered margins against a re-layout with a different pattern/repetitions — is
    // solved structurally by the #720 re-key.)
    if (_marginsPerSection.empty() && !_marginsExplicitlyCleared &&
        settings.get_coil_allow_margin_tape() && get_sections_description()) {
        auto sectionsWithMargins = get_sections_description().value();
        std::vector<std::vector<double>> recoveredMargins;
        recoveredMargins.reserve(sectionsWithMargins.size());
        bool anyMargin = false;
        // ABT #720: _marginsPerSection is keyed by CONDUCTION-section ordinal — recover one
        // entry per persisted conduction section, in wound order.
        for (const auto& sectionWithMargin : sectionsWithMargins) {
            if (sectionWithMargin.get_type() != ElectricalType::CONDUCTION) {
                continue;
            }
            auto margin = resolve_margin(sectionWithMargin);
            anyMargin = anyMargin || margin[0] > 0 || margin[1] > 0;
            recoveredMargins.push_back(margin);
        }
        if (anyMargin) {
            _marginsPerSection = recoveredMargins;
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
            // ABT #492 owner ruling: planar wires are PCBs — the real-winding connection model
            // (leads, markers, blocking, YZ-face dragbacks, connection losses) is for WOUND
            // magnetics only, and real winding for planar has not been started. Throw before ANY
            // of that machinery engages (the N+1 crossing bump armed just below is already part of
            // it), so every downstream path — marker emission, blocking, dragbacks, ohmic lead
            // lengths, apply_connection_reserved_space — is covered by one unavoidable gate.
            // Production planar flows are unaffected: the setting defaults to false.
            if (settings.get_coil_use_real_winding_geometry()) {
                for (const auto& wire : get_wires()) {
                    if (wire.get_type() == WireType::PLANAR) {
                        throw std::runtime_error(
                            "Real winding geometry (connection/lead routing) is not implemented for "
                            "planar (PCB) constructions; disable coilUseRealWindingGeometry for "
                            "planar magnetics");
                    }
                }
            }
            set_sections_description(std::nullopt);
            set_layers_description(std::nullopt);
            set_turns_description(std::nullopt);
            // ABT #728: from this point every turn in the description was produced by this
            // wind, so the crossing-bump destructor may zero its stations.
            realWindingCrossingBump.turnsAreFresh = true;
            // Start every wind from the ideal geometry: real-winding turn blocking, if any, is
            // re-derived and re-applied below only when the real-geometry setting is on.
            _applyConnectionBlocking = false;
            _connectionBlockedSlotsPerLayer.clear();
            _connectionBlockedDepthPerLayer.clear();
            _connectionBlockedRoomPerLayer.clear();
            _uLandingDepthPerLayer.clear();
            _terminalEntranceAtTop.clear();

            // Special case: toroid with one physical turn whose wire OD
            // exceeds the inner-hole radius. The wire cannot be wound on
            // the inner wall, so place it at the geometric centre of the
            // hole. Skips the sections/layers fit pipeline entirely.
            if (can_build_centered_single_turn_toroidal()) {
                logEntry("Building centered single-turn toroidal", "Coil", 2);
                return build_centered_single_turn_toroidal();
            }

            // REAL WINDING: arm the N+1-crossing bump (declared at function scope so it
            // also covers the global turn-blocking re-wind below) only after the
            // centered-single-turn special case has had its chance.
            realWindingCrossingBump.arm();

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
    // ABT #650: asking for real winding and silently not getting it is the worst outcome — the
    // caller receives a layout with none of the connection corridors reserved and nothing says so.
    // Real winding is applied only when the IDEAL wind fits (see below); when it does not, say it
    // out loud, and name what failed the fit so the reason is one log line away rather than a day
    // of bisecting.
    if (settings.get_coil_use_real_winding_geometry() && !result) {
        logEntry("Real winding was requested but connection blocking was NOT applied: the ideal "
                 "wind does not fit"
                 + std::string(get_turns_description() ? "" : " (no turns were produced)")
                 + (_lastFitFailure.empty() ? std::string() : " — " + _lastFitFailure),
                 "Coil", 1);
    }
    if (result && settings.get_coil_use_real_winding_geometry()) {
        // Turn blocking is GLOBAL to the winding window: a connection lead routes through the whole
        // window and removes a turn slot from every conduction layer it crosses, regardless of which
        // section/winding owns that layer. Derive that per-layer incidence from the wound geometry,
        // then re-wind so each crossed layer frees its blocked top/bottom slots and the affected
        // sections grow by extra layers to fit. Added layers shift radial positions (changing which
        // leads cross what), so iterate to a fixpoint. Entirely gated behind the real-winding flag —
        // ideal winding never enters this loop, so its geometry is unchanged.
        logEntry("Applying real winding geometry (global turn blocking)", "Coil", 2);
        // Upper bound: each iteration can at most add one blocked slot per layer edge and
        // spill one extra layer per section; window-sized windings converge in a handful,
        // and the post-loop check below turns genuine divergence into a loud error.
        const size_t maximumBlockingIterations = 16;
        size_t directionRegimeResets = 0;
        for (size_t blockingIteration = 0; blockingIteration < maximumBlockingIterations; ++blockingIteration) {
            std::map<std::string, bool> entranceEdgesBefore = _terminalEntranceAtTop;
            std::map<std::string, std::pair<double, double>> freshDepths;
            auto freshBlocked = compute_connection_blocked_slots_per_layer(&freshDepths);
            // ABT #616: the entrance-edge feedback flips a winding's direction, which relays
            // every one of its rows to the other edge — the monotone accumulation must NOT
            // union rows of two different direction regimes (measured: an 18-turn U design
            // accumulated 6-of-6 blocked slots from the superposition). On a regime change,
            // restart the accumulation and re-measure the new-direction geometry.
            if (_terminalEntranceAtTop != entranceEdgesBefore && directionRegimeResets < 2) {
                directionRegimeResets++;
                _connectionBlockedSlotsPerLayer.clear();
                _connectionBlockedDepthPerLayer.clear();
                _uLandingDepthPerLayer.clear();
                _applyConnectionBlocking = true;
                wind_by_sections(proportionPerWinding, pattern, repetitions);
                redistribute_section_turns_for_blocking();
                wind_by_layers();
                if (!get_layers_description()) {
                    break;
                }
                wind_by_turns();
                if (!get_turns_description()) {
                    break;
                }
                if (delimitAndCompact) {
                    delimit_and_compact();
                }
                align_blocked_layer_turns();
                continue;
            }
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
            // The continuous depths accumulate the same way, and drive re-iteration on their own:
            // a deeper run can appear WITHOUT changing the ceil'd slot count (1.05 -> 1.09 pitches
            // both block 2 slots), and the aligned turns hug these depths — so a depth-only growth
            // still moves turns and must be re-measured. Monotone and bounded (row stacks take
            // finitely many values), so convergence is unaffected.
            for (const auto& [layerName, edges] : freshDepths) {
                auto& accumulated = _connectionBlockedDepthPerLayer[layerName];
                if (edges.first > accumulated.first + 1e-9) {
                    accumulated.first = edges.first;
                    changed = true;
                }
                if (edges.second > accumulated.second + 1e-9) {
                    accumulated.second = edges.second;
                    changed = true;
                }
            }
            // ABT #608 (final form): U landing placement — a non-first U layer's span excludes one
            // wire OD past the tangential arrival, so its first station descends from the arrival
            // instead of sitting level with it. Depths ONLY (align spreads against them): no slot
            // counts, no capacity, no markers. REPLACED each iteration, not max-merged: the depth
            // is capped by the layer's current turn count, which the redistribution moves — see
            // _uLandingDepthPerLayer.
            {
                auto freshLanding = compute_u_landing_extra_depths();
                if (freshLanding.size() != _uLandingDepthPerLayer.size()) {
                    changed = true;
                }
                else {
                    for (const auto& [layerName, edges] : freshLanding) {
                        auto previous = _uLandingDepthPerLayer.find(layerName);
                        if (previous == _uLandingDepthPerLayer.end()
                            || std::abs(edges.first - previous->second.first) > 1e-9
                            || std::abs(edges.second - previous->second.second) > 1e-9) {
                            changed = true;
                            break;
                        }
                    }
                }
                _uLandingDepthPerLayer = std::move(freshLanding);
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
            // Wind turns even when the grown layout momentarily stops fitting: the loop MEASURES
            // blocking from the drawn geometry, so skipping the turns here makes the next iteration
            // measure STALE turns against the re-grown layers (phantom crossings) or, once the
            // stale turns are cleared, end the loop early with a turnless coil even though the
            // fixpoint exists (ABT #278: the ER 9.5 current-sense example died mid-loop at a
            // transiently-unfitting state; with the measurement geometry kept alive it converges
            // to a FITTING 15-layer layout). The post-loop `result` check still enforces the real
            // fitting contract; this only keeps the measurement geometry alive.
            wind_by_turns();
            if (!get_turns_description()) {
                break;
            }
            if (delimitAndCompact) {
                delimit_and_compact();
            }
            // Align INSIDE the loop: packing the blocked layers' turns against their
            // unblocked edge moves which slots the leads cross, so the next iteration's
            // blocking must be derived from the ALIGNED geometry — aligning only after
            // the loop leaves silently unblocked slots (turns inside terminal leads).
            align_blocked_layer_turns();
        }
        // RELAXATION (ABT #615, Alf 2026-08-09 on 25_psps: "layer 3 could fit more turns,
        // right?"): the monotone max above converges by keeping the DEEPEST reservation any
        // iteration ever measured, so a band that settled shallower leaves stale slots — S0's
        // layers were charged 3 bottom slots by early iterations whose final need was 2, holding
        // one fewer turn each and floating 0.43 mm off the corridor. After convergence, recompute
        // the blocking FRESH from the final geometry; if the accumulated maps over-reserve
        // anywhere, adopt the fresh values and re-wind once. If the relaxed layout re-introduces
        // residual blocking (the layout genuinely needed the conservative reservation), restore
        // the maps and re-wind back — deterministic, one extra pass, never a collision.
        // Repeated relax ROUNDS (ABT #616): one round settles into SOME self-consistent
        // state, but a state built with stale-fat maps verifies against those same fat maps
        // — wasteful yet 'holding'. Each new round re-measures the settled geometry; while
        // the fresh need is strictly smaller anywhere, adopt it and settle again. Bounded,
        // monotone in the applied maps across rounds, keeps the last state that held.
        for (size_t relaxRound = 0; relaxRound < 3; ++relaxRound) {
            if (std::getenv("MKF_BLOCKING_DIAG")) {
                std::cerr << "[relax] round=" << relaxRound << "\n";
            }
            std::map<std::string, std::pair<double, double>> freshDepths;
            auto freshBlocked = compute_connection_blocked_slots_per_layer(&freshDepths);
            bool overReserved = false;
            for (const auto& [layerName, accumulated] : _connectionBlockedSlotsPerLayer) {
                auto fresh = freshBlocked.find(layerName);
                uint64_t freshTop = fresh == freshBlocked.end() ? 0 : fresh->second.first;
                uint64_t freshBottom = fresh == freshBlocked.end() ? 0 : fresh->second.second;
                if (accumulated.first > freshTop || accumulated.second > freshBottom) {
                    overReserved = true;
                    break;
                }
            }
            if (!overReserved) {
                for (const auto& [layerName, accumulated] : _connectionBlockedDepthPerLayer) {
                    auto fresh = freshDepths.find(layerName);
                    double freshTop = fresh == freshDepths.end() ? 0.0 : fresh->second.first;
                    double freshBottom = fresh == freshDepths.end() ? 0.0 : fresh->second.second;
                    if (accumulated.first > freshTop + 1e-9 || accumulated.second > freshBottom + 1e-9) {
                        overReserved = true;
                        break;
                    }
                }
            }
            if (std::getenv("MKF_BLOCKING_DIAG")) {
                std::cerr << "[relax] overReserved=" << overReserved << "\n";
            }
            if (!overReserved) {
                break;   // nothing left to reclaim: the settled state is tight
            }
            {
                auto backupSlots = _connectionBlockedSlotsPerLayer;
                auto backupDepths = _connectionBlockedDepthPerLayer;
                auto backupLanding = _uLandingDepthPerLayer;
                auto rewindOnce = [&]() {
                    wind_by_sections(proportionPerWinding, pattern, repetitions);
                    redistribute_section_turns_for_blocking();
                    wind_by_layers();
                    if (get_layers_description()) {
                        wind_by_turns();
                        if (get_turns_description() && delimitAndCompact) {
                            delimit_and_compact();
                        }
                        align_blocked_layer_turns();
                    }
                };
                _connectionBlockedSlotsPerLayer = freshBlocked;
                _connectionBlockedDepthPerLayer = freshDepths;
                _uLandingDepthPerLayer = compute_u_landing_extra_depths();
                rewindOnce();
                // ABT #616: a single-shot relax fell back WHOLESALE whenever the relaxed
                // layout re-introduced any residual — 26_psps then shipped the conservative
                // degenerate layout (single-turn layers, 2-turn layers where 3 fit). Instead,
                // re-run the fixpoint FROM THE FRESH SEED: merge each pass's residual into
                // the applied maps (monotone from fresh, not from the old conservative state)
                // and re-wind, a bounded number of times. Fall back only if it never settles.
                bool relaxedHolds = false;
                const size_t maximumRelaxPasses = 4;
                for (size_t relaxPass = 0; relaxPass < maximumRelaxPasses; ++relaxPass) {
                    if (!get_turns_description()) {
                        break;   // relaxed layout failed to wind at all
                    }
                    std::map<std::string, std::pair<double, double>> residDepths;
                    auto relaxedResidual = compute_connection_blocked_slots_per_layer(&residDepths);
                    bool settled = true;
                    for (const auto& [layerName, edges] : relaxedResidual) {
                        auto applied = _connectionBlockedSlotsPerLayer.find(layerName);
                        uint64_t appliedTop = applied == _connectionBlockedSlotsPerLayer.end() ? 0 : applied->second.first;
                        uint64_t appliedBottom = applied == _connectionBlockedSlotsPerLayer.end() ? 0 : applied->second.second;
                        if (edges.first > appliedTop || edges.second > appliedBottom) {
                            settled = false;
                            break;
                        }
                    }
                    if (settled) {
                        relaxedHolds = true;
                        break;
                    }
                    if (relaxPass + 1 == maximumRelaxPasses) {
                        break;   // budget spent still growing: fall back below
                    }
                    for (const auto& [layerName, edges] : relaxedResidual) {
                        auto& applied = _connectionBlockedSlotsPerLayer[layerName];
                        applied.first = std::max(applied.first, edges.first);
                        applied.second = std::max(applied.second, edges.second);
                    }
                    for (const auto& [layerName, depths] : residDepths) {
                        auto& applied = _connectionBlockedDepthPerLayer[layerName];
                        applied.first = std::max(applied.first, depths.first);
                        applied.second = std::max(applied.second, depths.second);
                    }
                    _uLandingDepthPerLayer = compute_u_landing_extra_depths();
                    rewindOnce();
                }
                if (std::getenv("MKF_BLOCKING_DIAG")) {
                    std::cerr << "[relax] relaxedHolds=" << relaxedHolds << "\n";
                }
                if (!relaxedHolds) {
                    _connectionBlockedSlotsPerLayer = backupSlots;
                    _connectionBlockedDepthPerLayer = backupDepths;
                    _uLandingDepthPerLayer = backupLanding;
                    rewindOnce();
                    break;   // this round could not tighten: keep the last state that held
                }
            }
        }
        // Verify the fixpoint actually converged: the last re-wind may have produced NEW
        // blocking that the loop never re-applied (cap exhaustion) — silent residue leaves
        // turns inside reserved lead slots.
        {
            auto residual = compute_connection_blocked_slots_per_layer();
            for (const auto& [layerName, edges] : residual) {
                const auto& accumulated = _connectionBlockedSlotsPerLayer[layerName];
                if (edges.first > accumulated.first || edges.second > accumulated.second) {
                    if (std::getenv("MKF_BLOCKING_DIAG"))
                        std::cerr << "[fixpoint-fail] layer " << layerName << " needs {"
                                  << edges.first << "," << edges.second << "} applied {"
                                  << accumulated.first << "," << accumulated.second << "}\n";
                    throw CoilException(
                        ErrorCode::COIL_WINDING_ERROR,
                        "Real winding turn blocking did not converge for layer '" + layerName +
                        "' (needs {" + std::to_string(edges.first) + "," +
                        std::to_string(edges.second) + "} blocked slots, applied {" +
                        std::to_string(accumulated.first) + "," +
                        std::to_string(accumulated.second) + "})");
                }
            }
        }
        result = are_sections_and_layers_fitting() && bool(get_turns_description());
        // ABT #187: toroidal coils block ANGULAR corridors (radial leads crossing outer rings)
        // instead of top/bottom slots. align_blocked_ring_turns displaces each crossed ring's turns
        // clear of the corridors and reports the rings too full to clear; those get turn slots
        // reserved (capacity spills inward) and the coil re-winds — the toroidal analog of the
        // rectangular fixpoint above (which is a no-op for round windows and vice versa).
        if (resolve_bobbin().get_winding_window_shape() == WindingWindowShape::ROUND) {
            const size_t maximumToroidalBlockingIterations = 16;
            for (size_t blockingIteration = 0; blockingIteration < maximumToroidalBlockingIterations; ++blockingIteration) {
                auto ringDeficits = align_blocked_ring_turns();
                bool changed = false;
                for (const auto& [ringName, deficitSlots] : ringDeficits) {
                    if (deficitSlots == 0) {
                        continue;
                    }
                    // Monotonic accumulation, as in the rectangular loop: the deficit is measured
                    // against the CURRENT (already partially blocked) geometry, so add it on top.
                    _connectionBlockedSlotsPerLayer[ringName].first += deficitSlots;
                    changed = true;
                }
                if (std::getenv("MKF_BLOCKING_DIAG")) {
                    std::cerr << "[torfix] iter " << blockingIteration << " deficits={";
                    for (const auto& [ringName, deficitSlots] : ringDeficits) {
                        std::cerr << ringName << ":" << deficitSlots << ",";
                    }
                    std::cerr << "} changed=" << changed << "\n";
                }
                if (!changed) {
                    break;
                }
                _applyConnectionBlocking = true;
                wind_by_sections(proportionPerWinding, pattern, repetitions);
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
            // Owner ruling (ABT #723): the input terminal's internal connection OWNS its angular
            // corridor — no turn may sit there — and it is the displacement above
            // (align_blocked_ring_turns) that enforces it. The outer-return crossings, however,
            // were placed BEFORE displacement (tolerantly: the pre-blocking pass falls back
            // instead of throwing, because the corridor machinery had not run yet). Regenerate
            // them now on the FINAL displaced stations with the STRICT connection-aware sweep —
            // a failure here is a genuine corridor blockage, not a sequencing artifact.
            //
            // The regen MOVES the outer crossings, and the lead corridors derive from them, so
            // the corridors the displacement validated are stale after it: run one more
            // displacement pass against the regenerated leads, and — if it moved anything —
            // regenerate once more so crossings and corridors leave this function consistent.
            if (get_turns_description()) {
                bool applyConnectionBlockingBackup = _applyConnectionBlocking;
                _applyConnectionBlocking = true;   // strict sweep: no pre-blocking fallback
                generate_toroidal_additional_coordinates();
                auto turnsBeforeSettle = get_turns_description();
                align_blocked_ring_turns();
                bool settleMoved = false;
                if (get_turns_description() && turnsBeforeSettle) {
                    const auto& before = turnsBeforeSettle.value();
                    const auto after = get_turns_description().value();
                    settleMoved = before.size() != after.size();
                    for (size_t t = 0; !settleMoved && t < after.size(); ++t) {
                        settleMoved = std::abs(before[t].get_coordinates()[0] - after[t].get_coordinates()[0]) > 1e-9 ||
                                      std::abs(before[t].get_coordinates()[1] - after[t].get_coordinates()[1]) > 1e-9;
                    }
                }
                if (settleMoved) {
                    generate_toroidal_additional_coordinates();
                }
                _applyConnectionBlocking = applyConnectionBlockingBackup;
            }
            result = are_sections_and_layers_fitting() && bool(get_turns_description());
        }
        logEntry("Applying real winding geometry (connection reserved space)", "Coil", 2);
        apply_connection_reserved_space();
    }
    if (result) {
        // Multi-column winding: groups were wound in the +x window-local frame;
        // mirror the ones whose winding window sits on the negative-x side into
        // their real position. No-op for single-window coils.
        apply_group_window_sides();
        // Hand-drawn section rectangles (winding studio) override the computed
        // placement LAST — after compaction and mirroring — so no later pass
        // can move them; layers+turns are re-flowed inside each drawn rect.
        result = apply_custom_section_rects();
        if (result && !delimitAndCompact) {
            // Compact-off toroid winds skipped delimit_and_compact_round_window,
            // the only place the outer return crossings are normally generated.
            generate_toroidal_additional_coordinates();
        }
    }
    // ABT #624: a wind that leaves copper outside the winding window has not succeeded, however
    // well its filling factors read — nothing else in this function compares the turns against
    // the window at all, which is how 26_psps shipped a turn 0.163 mm past the edge (and, in 3D,
    // inside the bobbin flange).
    if (result && !are_turns_inside_winding_window()) {
        result = false;
    }
    return result;
}

bool Coil::apply_custom_section_rects() {
    if (_customSectionRects.empty() || !get_sections_description()) {
        return true;
    }
    auto sections = get_sections_description().value();
    bool anyApplied = false;
    for (auto& section : sections) {
        auto customRect = _customSectionRects.find(section.get_name());
        if (customRect == _customSectionRects.end()) {
            continue;
        }
        section.set_coordinates(customRect->second.first);
        section.set_dimensions(customRect->second.second);
        // The drawn rectangle invalidates the previously wound layer count —
        // wind_by_layers honors a section's numberLayers when present, which
        // would keep e.g. a single overflowing layer instead of re-packing
        // into the new rect.
        section.set_number_layers(std::nullopt);
        anyApplied = true;
    }
    if (!anyApplied) {
        // Stale names (the pattern changed and the drawn sections no longer
        // exist) are the caller's to clean up; the wound result stays valid.
        return true;
    }
    set_sections_description(sections);
    return rewind_layers_and_turns();
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
std::optional<WindingStyle> Coil::get_winding_style_override(size_t windingIndex) const {
    if (_windingStyleOverridePerWinding.empty() || windingIndex >= get_functional_description().size()) {
        return std::nullopt;
    }
    auto windingName = get_functional_description()[windingIndex].get_name();
    auto it = _windingStyleOverridePerWinding.find(windingName);
    if (it == _windingStyleOverridePerWinding.end()) {
        return std::nullopt;
    }
    return it->second;
}

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

        // User override (winding studio): wins over every heuristic below; the
        // physical must-case above stays dominant.
        if (auto styleOverride = get_winding_style_override(i)) {
            windByConsecutiveTurns.push_back(styleOverride.value());
            log("Winding " + std::to_string(i) + ": user winding-style override.");
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
WindingStyle Coil::wind_by_consecutive_turns(uint64_t numberTurns, uint64_t numberParallels, size_t numberSlots, std::optional<size_t> windingIndex) {
    // When turns < slots, we MUST use CONSECUTIVE_TURNS to distribute physical turns (parallels)
    // across slots. CONSECUTIVE_PARALLELS would put all turns in the first slots, leaving rest empty.
    if (numberTurns < numberSlots && numberParallels > 1) {
        log("Layer: CONSECUTIVE_TURNS (turns < slots, must distribute parallels across slots).");
        return WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
    }

    // User override (winding studio): wins over the heuristics below; the
    // physical must-case above stays dominant.
    if (windingIndex) {
        if (auto styleOverride = get_winding_style_override(windingIndex.value())) {
            log("Layer: user winding-style override.");
            return styleOverride.value();
        }
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
        for (const auto& turn : turns) {
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
            if (std::getenv("MKF_BLOCKING_DIAG"))
                std::cerr << "[fit] section " << section.get_name() << " ff="
                          << section.get_filling_factor().value() << " ovl="
                          << overlapping_filling_factor(section) << " cont="
                          << contiguous_filling_factor(section) << "\n";
            windTurns = false;
        }
    }
    for (auto& layer: layers) {
        if (roundFloat(layer.get_filling_factor().value(), 6) > 1) {
            if (std::getenv("MKF_BLOCKING_DIAG"))
                std::cerr << "[fit] layer " << layer.get_name() << " ff="
                          << layer.get_filling_factor().value() << " dims=("
                          << layer.get_dimensions()[0] * 1e3 << "x"
                          << layer.get_dimensions()[1] * 1e3 << ")\n";
            windTurns = false;
        }
    }

    // ABT #616: nothing above compares against the WINDOW — real-winding blocking can grow a
    // section's layer count until its layers walk radially past the winding window edge, and
    // the coil still reported "fits" (26_psps U: Secondary section 1 reached x=12.13 in a
    // window ending at 10.42 — 1.7 mm of silently overflowing copper). Every conduction
    // layer must lie inside the window envelope.
    //
    // Wound coils only (ABT #675): planar does not support real winding yet, so there is no
    // blocking here to grow anything past the window and this check guards nothing — while it
    // does reject PCB layouts whose outer turns sit past the window edge, which took the planar
    // coil adviser to zero candidates for dense designs the moment the window was read
    // correctly (ABT #650). When planar gains real winding, this needs a planar-aware envelope
    // (the board is not the window), not a straight re-enable.
    if (!is_planar()) {
        auto bobbin = resolve_bobbin();
        if (bobbin.get_winding_window_shape() == WindingWindowShape::RECTANGULAR &&
            bobbin.get_processed_description()) {
            // ABT #650: the generated getter returns the processed description BY VALUE and
            // get_winding_windows() hands out a reference INTO that temporary, so binding the
            // window to `...get_processed_description()->get_winding_windows()[0]` left it
            // dangling the instant the full-expression ended. Natively the freed bytes still
            // read back right; under Emscripten the block is recycled and the window's
            // coordinates came back as (0,0) — the window jumped to the core axis and this
            // check failed on the FIRST turn of every concentric design, silently dropping
            // real-winding blocking in the browser. Keep the description alive by name.
            const auto processedDescription = bobbin.get_processed_description().value();
            if (processedDescription.get_winding_windows().empty()) {
                return windTurns;
            }
            // ABT #730: each turn is tested against ITS OWN winding window's envelope,
            // resolved through its section — a lateral group in window 1+ has a different
            // x-range than window 0, and testing everything against window 0 either refused
            // designs that fit or silently dropped real-winding blocking for them.
            //
            // Frame awareness (the 2026-08-14 attempt broke without it): this check runs
            // BEFORE apply_group_window_sides mirrors lateral groups to negative x, where
            // turns still sit in the +x window-local frame — there the envelope is the
            // window MIRRORED to +x (|x-center|). After the mirror (flag set), turns and
            // windows are both in the real frame. Every caller is a wind-time path, so the
            // flag is accurate for the coordinates being tested.
            const auto& windingWindowsForEnvelope = processedDescription.get_winding_windows();
            auto envelopeForWindow = [&](size_t windowIndex) -> std::optional<std::array<double, 4>> {
                if (windowIndex >= windingWindowsForEnvelope.size()) {
                    // Stale/foreign window index on a section: keep the historical window-0
                    // behaviour rather than rejecting the whole coil on bookkeeping.
                    windowIndex = 0;
                }
                const auto& ww = windingWindowsForEnvelope[windowIndex];
                if (!(ww.get_coordinates() && ww.get_width() && ww.get_height())) {
                    return std::nullopt;
                }
                double xCenter = (*ww.get_coordinates())[0];
                if (!_groupWindowSidesApplied) {
                    xCenter = std::abs(xCenter);
                }
                return std::array<double, 4>{xCenter - *ww.get_width() / 2,
                                             xCenter + *ww.get_width() / 2,
                                             (*ww.get_coordinates())[1] - *ww.get_height() / 2,
                                             (*ww.get_coordinates())[1] + *ww.get_height() / 2};
            };
            std::map<std::string, size_t> windowIndexBySectionName;
            for (const auto& sectionForWindow : sections) {
                windowIndexBySectionName[sectionForWindow.get_name()] =
                    resolve_section_winding_window_index(sectionForWindow);
            }
            // Measured on the TURNS (the actual copper): the layer rects go stale by a
            // few um once align_blocked_layer_turns re-spreads the turns, and a stale rect
            // must not fail a coil whose copper sits exactly at the window edge.
            const double tol = 1e-9;
            if (get_turns_description()) {
                auto wires = get_wires();
                auto turnsForEnvelope = get_turns_description().value();
                for (const auto& turn : turnsForEnvelope) {
                    size_t windowIndex = 0;
                    if (turn.get_section()) {
                        auto foundWindow = windowIndexBySectionName.find(turn.get_section().value());
                        if (foundWindow != windowIndexBySectionName.end()) {
                            windowIndex = foundWindow->second;
                        }
                    }
                    auto envelope = envelopeForWindow(windowIndex);
                    if (!envelope) {
                        continue;
                    }
                    const auto& [x0, x1, y0, y1] = envelope.value();
                    size_t windingIndex = get_winding_index_by_name(turn.get_winding());
                    const double hw = wires[windingIndex].get_maximum_outer_width() / 2;
                    const double hh = wires[windingIndex].get_maximum_outer_height() / 2;
                    const auto& c = turn.get_coordinates();
                    if (c[0] - hw < x0 - tol || c[0] + hw > x1 + tol ||
                        c[1] - hh < y0 - tol || c[1] + hh > y1 + tol) {
                        // Logged, not just cerr'd: getenv/cerr is unreachable from a WASM
                        // consumer, and this is the check that decides whether real-winding
                        // blocking runs at all (ABT #650).
                        _lastFitFailure = turn.get_name() + " at (" + std::to_string(c[0]) + ","
                                        + std::to_string(c[1]) + ") outside window "
                                        + std::to_string(windowIndex) + " x["
                                        + std::to_string(x0) + "," + std::to_string(x1) + "] y["
                                        + std::to_string(y0) + "," + std::to_string(y1) + "]";
                        windTurns = false;
                        break;
                    }
                }
            }
        }
    }

    return windTurns;
}

double Coil::overlapping_filling_factor(const Section& section) {
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

double Coil::contiguous_filling_factor(const Section& section) {
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

Coil::FillingFactorsOutput Coil::calculate_filling_factor(size_t groupIndex) {
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

    for (const auto& section : sections) {
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

    for (const auto& section : sections) {
        if (section.get_margin()) {
            double marginWidth = resolve_margin(section)[0] + resolve_margin(section)[1];
            if (bobbinWindingWindowShape != WindingWindowShape::RECTANGULAR) {
                // ABT #245: on a toroid, section dimensions are [radial height (m),
                // angular span (DEGREES)]. Margins are tape widths along the winding
                // path, so the band they occupy is marginWidth x radial height — using
                // the angular span here multiplied metres by degrees and inflated the
                // area by ~50x per section (the primary of the #306 fixture alone
                // contributed 0.0476 "metre-degrees" against a 0.000452 m2 window,
                // i.e. an area fill of 105). It went unnoticed because the result was
                // then discarded by max(maximumLayerFillingFactor, ...).
                area += marginWidth * section.get_dimensions()[0];
            }
            else if (windingOrientation == WindingOrientation::OVERLAPPING) {
                area += marginWidth * section.get_dimensions()[0];
            }
            else {
                area += marginWidth * section.get_dimensions()[1];
            }
        }
    }

    for (const auto& layer : layers) {
        // Track the true maximum, not just overflows: this value is reported now, and a
        // healthy coil should show its real headroom (e.g. 0.49) rather than a placeholder 0.
        if (layer.get_filling_factor()) {
            maximumLayerFillingFactor = std::max(maximumLayerFillingFactor, layer.get_filling_factor().value());
        }
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            auto turns = get_turns_by_layer(layer.get_name());
            for (const auto& turn : turns) {
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

    FillingFactorsOutput output;
    // ABT #245: the area fraction is reported as the area fraction. It used to be
    // max(maximumLayerFillingFactor, area / availableArea), which smuggled a per-layer
    // OVERFILL RATIO into this slot — a coil whose real areal fill was 2.35% reported
    // 43264% because one degenerate 0.0048-degree section could not hold its turns, and
    // the builder printed that as a percentage. The overfill is still reported, next to
    // it, and windingFits carries the verdict that consumers actually want.
    output.areaFillingFactor = area / availableArea;
    output.maxLayerFillingFactor = maximumLayerFillingFactor;
    output.contiguousFillingFactor = contiguousDimension / availableContiguousDimension;
    output.overlappingFillingFactor = overlappingDimension / availableOverlappingDimension;

    // Only the dimension the sections actually stack along can overflow: sections laid
    // out contiguously grow along the contiguous axis (height, or angle on a toroid),
    // overlapping ones along the overlapping axis (width, or radial height).
    double stackingFillingFactor = windingOrientation == WindingOrientation::CONTIGUOUS
                                       ? output.contiguousFillingFactor
                                       : output.overlappingFillingFactor;
    output.windingFits = output.areaFillingFactor <= 1 &&
                         output.maxLayerFillingFactor <= 1 &&
                         stackingFillingFactor <= 1;
    return output;
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

    // ABT #415: an inter-SECTION insulation strip separates STACKED SECTIONS, so its geometry
    // follows the winding orientation — same normalization its two mechanical siblings
    // (calculate_mechanical_insulation, calculate_custom_thickness_insulation) already apply.
    // This function forgot it, so a coil mixing winding OVERLAPPING with layers CONTIGUOUS built
    // TRANSPOSED insulation layers (a full-window-width horizontal tape inside a thin vertical
    // strip section): overlapping_filling_factor reported 244x and are_sections_and_layers_fitting
    // rejected EVERY isolated multi-winding candidate the CoilAdviser proposed.
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
                layer.set_orientation(layersOrientation);  // ABT #415: the NORMALIZED orientation, like the mechanical siblings
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

            Section section;
            section.set_name("temp");
            section.set_partial_windings(std::vector<PartialWinding>{});
            section.set_layers_orientation(layersOrientation);  // ABT #415: normalized, like the mechanical siblings
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

    // Insulation entries carry SIZE_MAX in the winding-index slot: they belong to no winding.
    // The entries are TYPED (ElectricalType::INSULATION) and every consumer reads the winding
    // index only from CONDUCTION entries, so the placeholder is never interpreted (ABT #720).
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

std::pair<size_t, std::vector<int64_t>> get_number_layers_needed_and_number_physical_turns(double radialHeight, double angle, Wire wire, int64_t physicalTurnsInSection, double windingWindowRadius, const std::vector<int64_t>* blockedSlotsPerLayer = nullptr) {
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
        // ABT #187: connection leads crossing this ring block angular turn slots (real winding
        // geometry) — the ring holds that many fewer turns, spilling them inward to the next ring.
        int64_t blockedSlots = (blockedSlotsPerLayer != nullptr && numberLayers < blockedSlotsPerLayer->size())
            ? (*blockedSlotsPerLayer)[numberLayers] : 0;
        int64_t numberTurnsFittingThisLayer = std::max(1.0, floor(sectionAvailableAngle / wireAngle) - double(blockedSlots));
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

std::pair<size_t, std::vector<int64_t>> get_number_layers_needed_and_number_physical_turns(Section section, Wire wire, int64_t physicalTurnsInSection, double windingWindowRadius, const std::vector<int64_t>* blockedSlotsPerLayer = nullptr) {
    return get_number_layers_needed_and_number_physical_turns(section.get_coordinates()[0] - section.get_dimensions()[0] / 2, section.get_dimensions()[1], wire, physicalTurnsInSection, windingWindowRadius, blockedSlotsPerLayer);
}

void Coil::apply_margin_tape(const std::vector<std::pair<ElectricalType, std::pair<size_t, double>>>& orderedSectionsWithInsulation, size_t conductionSectionOffset) {
    // ABT #720: _marginsPerSection is keyed by CONDUCTION-section ordinal (flat across groups,
    // in wound order). conductionSectionOffset is the number of conduction sections already
    // wound by previous groups.
    size_t conductionCount = 0;
    for (const auto& orderedSection : orderedSectionsWithInsulation) {
        if (orderedSection.first == ElectricalType::CONDUCTION) {
            ++conductionCount;
        }
    }
    if (_marginsPerSection.size() < conductionSectionOffset + conductionCount) {
        // Resize (not replace) so preloaded margins are preserved, matching equalize_margins
        _marginsPerSection.resize(conductionSectionOffset + conductionCount, {0, 0});
    }

    size_t conductionOrdinal = conductionSectionOffset;
    for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
        if (orderedSectionsWithInsulation[sectionIndex].first != ElectricalType::CONDUCTION) {
            continue;
        }
        size_t marginIndex = conductionOrdinal++;
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
            // No interface recorded for this pair means no margin-tape requirement (the
            // custom-insulation path fills _insulationSections without interfaces). Skip
            // instead of operator[]-defaulting, which used to insert a junk entry whose
            // layerPurpose is uninitialized.
            auto coilSectionInterfaceIt = _coilSectionInterfaces.find(windingsMapKey);
            if (coilSectionInterfaceIt == _coilSectionInterfaces.end()) {
                continue;
            }
            double halfMarginTapeDistance = coilSectionInterfaceIt->second.get_total_margin_tape_distance() / 2;
            // This conduction section and the PREVIOUS conduction section (across the
            // insulation between them) both get the tape floor.
            _marginsPerSection[marginIndex][0] =  std::max(_marginsPerSection[marginIndex][0], halfMarginTapeDistance);
            _marginsPerSection[marginIndex][1] =  std::max(_marginsPerSection[marginIndex][1], halfMarginTapeDistance);
            _marginsPerSection[marginIndex - 1][0] =  std::max(_marginsPerSection[marginIndex - 1][0], halfMarginTapeDistance);
            _marginsPerSection[marginIndex - 1][1] =  std::max(_marginsPerSection[marginIndex - 1][1], halfMarginTapeDistance);
        }
    }
}

void Coil::equalize_margins(const std::vector<std::pair<ElectricalType, std::pair<size_t, double>>>& orderedSectionsWithInsulation, size_t conductionSectionOffset) {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    // ABT #720: margins are keyed by CONDUCTION-section ordinal. Map each ordered entry to
    // its ordinal once; the ordered list is still walked for the available-space lookups.
    std::vector<std::optional<size_t>> ordinalByOrderedIndex(orderedSectionsWithInsulation.size(), std::nullopt);
    size_t conductionCount = 0;
    for (size_t index = 0; index < orderedSectionsWithInsulation.size(); ++index) {
        if (orderedSectionsWithInsulation[index].first == ElectricalType::CONDUCTION) {
            ordinalByOrderedIndex[index] = conductionSectionOffset + conductionCount;
            ++conductionCount;
        }
    }
    // Mirror apply_margin_tape's sizing guard: _marginsPerSection is sized
    // lazily by wind_by_*; equalize_margins can be reached on paths where it
    // was never grown to the section count (e.g. PSFB / multi-section bridge
    // topologies). Reading past the end here used to SEGV in CoilAdviser.
    if (_marginsPerSection.size() < conductionSectionOffset + conductionCount) {
        _marginsPerSection.resize(conductionSectionOffset + conductionCount, {0, 0});
    }

    // ABT #721: what "equalizing" means depends on the window's geometry.
    //
    //   ROUND (toroid): adjacent sections are angular SECTORS, so the left sector's END
    //   margin ([1]) and the right sector's START margin ([0]) are collinear along the
    //   angle and face each other directly across the insulation — THEY form the
    //   inter-winding gap, and the pair is redistributed proportionally to each sector's
    //   allotted space. The window closes on itself at 0°, so the last sector's partner
    //   wraps to the first (across the trailing insulation add_insulation_to_sections
    //   appends for exactly this reason).
    //
    //   RECTANGULAR: margins do NOT sit between the sections — they sit on the SAME two
    //   window edges of every section (top/bottom for overlapping, left/right for
    //   contiguous), and the creepage path between two adjacent windings runs over each
    //   edge as (left section's edge margin + insulation + right section's edge margin).
    //   So the redistribution pairs SAME-index margins — left[0] with right[0], left[1]
    //   with right[1] — conserving each edge's path total while splitting it by section
    //   space. And there is NO wrap: the window ends at the bobbin walls, the last and
    //   first sections are not neighbours. (The pre-#721 wrap pairing applied to
    //   rectangular would have redistributed the two window-wall margins between
    //   non-adjacent sections.)
    //
    // Preloaded margins are redistributed like tape-derived ones: the setting's contract
    // is "the coil may re-split inter-winding margins by section size"; callers that
    // need their exact preloaded values disable coilEqualizeMargins.
    bool windowWraps = bobbin.get_winding_window_shape() == WindingWindowShape::ROUND;

    for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
        if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
            if (!_coilSectionInterfaces.empty()) {

                size_t indexForMarginLeftSection = sectionIndex;
                size_t indexForMarginRightSection;
                // The "right" section is two ahead (conduction → insulation →
                // conduction). When near the end, wrap to the first section —
                // but only where the window physically closes on itself.
                // The original `!= size()-2` test missed the size()-1 case,
                // letting sectionIndex+2 read past the end.
                if (sectionIndex + 2 < orderedSectionsWithInsulation.size()) {
                    indexForMarginRightSection = sectionIndex + 2;
                }
                else if (windowWraps) {
                    indexForMarginRightSection = 0;
                }
                else {
                    continue;   // rectangular window: the last section has no next neighbour
                }
                // Rectangular only: redistribute solely across an actual insulation
                // boundary — wound_with center-tap halves sit adjacent WITHOUT
                // insulation (same isolation side) and carry no inter-winding tape to
                // re-split. Toroids keep the historical blind pairing exactly (their
                // pinned geometry predates this check).
                if (!windowWraps && sectionIndex + 1 < orderedSectionsWithInsulation.size() &&
                    orderedSectionsWithInsulation[sectionIndex + 1].first != ElectricalType::INSULATION) {
                    continue;
                }

                double leftAvailableSpace = orderedSectionsWithInsulation[indexForMarginLeftSection].second.second;
                double rightAvailableSpace = orderedSectionsWithInsulation[indexForMarginRightSection].second.second;
                double totalAvailableSpace = leftAvailableSpace + rightAvailableSpace;
                if (totalAvailableSpace <= 0) {
                    continue;
                }
                // ABT #720: the partner entries must both be conduction sections to have
                // margins at all (the +2 step assumed strict alternation; an insulation
                // entry landing there used to read/write a junk margin row).
                if (!ordinalByOrderedIndex[indexForMarginLeftSection] || !ordinalByOrderedIndex[indexForMarginRightSection]) {
                    continue;
                }
                size_t leftOrdinal = ordinalByOrderedIndex[indexForMarginLeftSection].value();
                size_t rightOrdinal = ordinalByOrderedIndex[indexForMarginRightSection].value();
                if (windowWraps) {
                    double totalMargin = _marginsPerSection[leftOrdinal][1] + _marginsPerSection[rightOrdinal][0];
                    _marginsPerSection[leftOrdinal][1] = leftAvailableSpace / totalAvailableSpace * totalMargin;
                    _marginsPerSection[rightOrdinal][0] = rightAvailableSpace / totalAvailableSpace * totalMargin;
                }
                else {
                    for (size_t edge : {size_t(0), size_t(1)}) {
                        double totalMargin = _marginsPerSection[leftOrdinal][edge] + _marginsPerSection[rightOrdinal][edge];
                        _marginsPerSection[leftOrdinal][edge] = leftAvailableSpace / totalAvailableSpace * totalMargin;
                        _marginsPerSection[rightOrdinal][edge] = rightAvailableSpace / totalAvailableSpace * totalMargin;
                    }
                }
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
    split_shared_window_groups(groups, windingWindows);
    set_groups_description(groups);

    if (!anyExplicitPlacement) {
        OM_WARNING("Multi-column bobbin detected (" + std::to_string(windingWindows.size()) +
                   " winding windows) and no winding carries a windingWindow placement. All windings "
                   "placed in window 0 by default. Set windingWindow on the windings or call "
                   "assign_windings_to_columns() to distribute.");
    }
    return true;
}

void Coil::split_shared_window_groups(std::vector<Group>& groups, const std::vector<WindingWindowElement>& windingWindows) {
    // Region sharing: the main-column winding forms an annulus whose two crossings
    // occupy the inner side of BOTH window regions, while a lateral winding hugs its
    // leg on the outer side of its region. When both kinds are wound they share the
    // region width: main-wound groups keep the inner half, lateral-wound groups the
    // outer half. With only one kind wound, every group keeps the full region.
    if (windingWindows.size() <= 1) {
        return;
    }
    if (!windingWindows[0].get_width()) {
        // Radial (toroidal) windows share by angle, not width; nothing to split.
        return;
    }
    auto mainColumnEdge = windingWindows[0].get_column();
    auto isLateralWound = [&](const Group& group) {
        size_t windowIndex = group.get_winding_window() ? static_cast<size_t>(group.get_winding_window().value()) : 0;
        if (windowIndex >= windingWindows.size()) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                "Group " + group.get_name() + " references winding window " + std::to_string(windowIndex) +
                " but the bobbin has " + std::to_string(windingWindows.size()) + " winding windows");
        }
        auto columnEdge = windingWindows[windowIndex].get_column();
        return bool(columnEdge && (!mainColumnEdge || columnEdge.value() != mainColumnEdge.value()));
    };

    bool anyMainWound = false;
    bool anyLateralWound = false;
    for (auto& group : groups) {
        if (group.get_partial_windings().empty()) {
            continue;
        }
        if (isLateralWound(group)) {
            anyLateralWound = true;
        }
        else {
            anyMainWound = true;
        }
    }
    if (!anyLateralWound) {
        return;
    }

    // ABT #228.4: two (or more) lateral-wound groups sharing a region with NO main-wound
    // group present are not handled below -- the split only fires when
    // anyMainWound && anyLateralWound, so two lateral groups whose regions genuinely
    // overlap (same coordinates/dimensions) would otherwise fall through untouched and
    // silently wind into the same physical space. Rather than guess an even N-way split
    // for a configuration no fixture yet exercises, fail loudly: the caller needs an
    // explicit placement (or a real N-way split, once a use case defines what it should
    // look like) instead of overlapping copper it cannot see.
    if (!anyMainWound) {
        std::vector<size_t> woundGroupIndexes;
        for (size_t i = 0; i < groups.size(); ++i) {
            if (!groups[i].get_partial_windings().empty()) {
                woundGroupIndexes.push_back(i);
            }
        }
        for (size_t a = 0; a < woundGroupIndexes.size(); ++a) {
            for (size_t b = a + 1; b < woundGroupIndexes.size(); ++b) {
                auto& groupA = groups[woundGroupIndexes[a]];
                auto& groupB = groups[woundGroupIndexes[b]];
                auto coordinatesA = groupA.get_coordinates();
                auto coordinatesB = groupB.get_coordinates();
                auto dimensionsA = groupA.get_dimensions();
                auto dimensionsB = groupB.get_dimensions();
                bool sameRegion = coordinatesA.size() >= 2 && coordinatesB.size() >= 2 &&
                    std::abs(coordinatesA[0] - coordinatesB[0]) < 1e-9 &&
                    std::abs(coordinatesA[1] - coordinatesB[1]) < 1e-9 &&
                    std::abs(dimensionsA[0] - dimensionsB[0]) < 1e-9 &&
                    std::abs(dimensionsA[1] - dimensionsB[1]) < 1e-9;
                if (sameRegion) {
                    throw NotImplementedException(
                        "Groups " + groupA.get_name() + " and " + groupB.get_name() +
                        " are both lateral-wound and share the same winding-window region "
                        "with no main-wound group to split against: N-way region sharing "
                        "between lateral-only groups is not implemented yet");
                }
            }
        }
        // No two lateral-wound groups actually overlap (the common case: each lateral
        // leg has its own, non-shared window) -- fall through to the per-group loop
        // below, which still trims each lateral group's own column-wall thickness even
        // though anyMainWound is false and the coreWindowMidline branch will not fire.
    }

    // The lateral bobbin's column wall sits between the leg's face and the winding
    // space (mirroring the main bobbin's wall against the main column), so every
    // lateral-wound group gives up one columnThickness at its leg. Under sharing,
    // the region is split at the CORE window's midline (main bobbin wall to leg
    // face), which hands both sides equal winding space: wall + space inward, space
    // + wall outward.
    double columnThickness = resolve_bobbin().get_processed_description().value().get_column_thickness();

    for (auto& group : groups) {
        if (group.get_partial_windings().empty()) {
            continue;
        }
        bool lateralWound = isLateralWound(group);
        auto dimensions = group.get_dimensions();
        auto coordinates = group.get_coordinates();
        double sideSign = coordinates[0] >= 0 ? 1.0 : -1.0;
        double regionInnerEdge = std::abs(coordinates[0]) - dimensions[0] / 2;
        double regionOuterEdge = std::abs(coordinates[0]) + dimensions[0] / 2;
        double newInnerEdge = regionInnerEdge;
        double newOuterEdge = lateralWound ? regionOuterEdge - columnThickness : regionOuterEdge;
        if (anyMainWound && anyLateralWound) {
            // Core-window midline: from the main bobbin's wall face (inner edge minus
            // its wall) to the leg face (region outer edge).
            double coreWindowMidline = (regionInnerEdge - columnThickness + regionOuterEdge) / 2;
            if (lateralWound) {
                newInnerEdge = coreWindowMidline;
            }
            else {
                newOuterEdge = coreWindowMidline;
            }
        }
        coordinates[0] = sideSign * (newInnerEdge + newOuterEdge) / 2;
        dimensions[0] = newOuterEdge - newInnerEdge;
        group.set_coordinates(coordinates);
        group.set_dimensions(dimensions);
    }
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
    // The lateral bobbin wraps its leg like the main one wraps the main column:
    // one column-wall thickness sits between the leg and the winding space.
    frame.columnWidth = column.get_width() / 2 + bobbinProcessedDescription.get_column_thickness();
    frame.columnDepth = column.get_depth() / 2 + bobbinProcessedDescription.get_column_thickness();
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

void Coil::apply_group_window_sides(bool inverse) {
    if (!get_groups_description() || !get_sections_description()) {
        return;
    }
    if (inverse && !_groupWindowSidesApplied) {
        return;
    }
    if (!inverse && _groupWindowSidesApplied) {
        return;
    }
    Bobbin bobbinResolved = resolve_bobbin();
    if (!bobbinResolved.get_processed_description()) {
        return;
    }
    auto windingWindows = bobbinResolved.get_processed_description().value().get_winding_windows();
    // Single-window coils get NO additionalCoordinates from this pass: several
    // physics consumers (StrayCapacitance surrounding-turn search, Temperature
    // outer nodes, CoilMesher phantom conductors) interpret them as REAL
    // conductor positions, so emitting the center-leg mirror here would
    // silently change their results for every classic design. The winding
    // studio synthesizes that mirror as pure display geometry instead.
    // Multi-window (multi-column) coils keep their crossings — those designs
    // are new with this machinery and their consumers are gated (e.g. the
    // LeakageInductance round-window guard). Toroids keep their own
    // outer-return machinery untouched.
    bool wantsBothCrossings = settings.get_coil_include_additional_coordinates() && !inverse
        && windingWindows.size() > 1
        && bobbinResolved.get_winding_window_shape() == WindingWindowShape::RECTANGULAR;
    if (windingWindows.size() <= 1 && !wantsBothCrossings) {
        return;
    }

    // Per-section placement transform out of the +x winding frame:
    // 1. Sections wound around a NON-main column are reflected across their window
    //    center, so the winding hugs the column it actually wraps (the lateral leg's
    //    face is the window's outer edge; the winder laid them against the inner edge
    //    like a main-column winding). Their turn lengths are recomputed for the new
    //    radius around the leg.
    // 2. Sections whose window sits on the negative-x side are mirrored into place.
    // 3. Every turn gets its SECOND cross-section crossing as an additional
    //    coordinate: a turn around any column intersects the drawing plane twice
    //    (main-column turns on the far side of the main column, lateral-column turns
    //    outside the core), mirroring what toroidal turns already carry.
    struct SectionWindowTransform {
        bool reflectAcrossWindowCenter = false;
        double windingFrameWindowCenterX = 0;
        bool mirrorSide = false;
        double finalColumnAxisX = 0;
    };
    std::map<std::string, SectionWindowTransform> transformPerSection;

    auto mainColumnEdge = windingWindows[0].get_column();
    // Reflection pivots on the GROUP's allocated sub-region (which may be half the
    // window when the region is shared with the main winding's annulus), falling
    // back to the window itself for sections without a group.
    std::map<std::string, double> groupCenterXByName;
    auto placementGroups = get_groups_description().value();
    for (auto& group : placementGroups) {
        groupCenterXByName[group.get_name()] = group.get_coordinates()[0];
    }
    auto sections = get_sections_description().value();
    bool anyTransform = false;
    for (auto& section : sections) {
        size_t windowIndex = resolve_section_winding_window_index(section);
        if (windowIndex >= windingWindows.size()) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                "Section " + section.get_name() + " references winding window " + std::to_string(windowIndex) +
                " but the bobbin has " + std::to_string(windingWindows.size()) + " winding windows");
        }
        auto& windingWindow = windingWindows[windowIndex];
        SectionWindowTransform transform;
        transform.mirrorSide = windingWindow.get_coordinates().value()[0] < 0;
        auto columnEdge = windingWindow.get_column();
        // Schema default: a window without a column edge wraps the main column.
        bool lateralWound = columnEdge && (!mainColumnEdge || columnEdge.value() != mainColumnEdge.value());
        transform.reflectAcrossWindowCenter = lateralWound;
        if (section.get_group() && groupCenterXByName.contains(section.get_group().value())) {
            transform.windingFrameWindowCenterX = std::abs(groupCenterXByName[section.get_group().value()]);
        }
        else {
            transform.windingFrameWindowCenterX = std::abs(windingWindow.get_coordinates().value()[0]);
        }
        if (lateralWound) {
            double windingFrameAxisX = get_wound_column_frame_for_section(section.get_name()).axisX;
            transform.finalColumnAxisX = transform.mirrorSide ? -windingFrameAxisX : windingFrameAxisX;
        }
        transformPerSection[section.get_name()] = transform;
        anyTransform = anyTransform || transform.reflectAcrossWindowCenter || transform.mirrorSide;
    }

    bool emitBothCrossings = wantsBothCrossings;
    if (!anyTransform && !emitBothCrossings) {
        return;
    }

    auto transformX = [inverse](double x, const SectionWindowTransform& transform) {
        // Forward: reflect (hug the wound leg) then mirror (negative-x side).
        // Inverse: same involutions in reverse order.
        if (inverse && transform.mirrorSide) {
            x = -x;
        }
        if (transform.reflectAcrossWindowCenter) {
            x = 2 * transform.windingFrameWindowCenterX - x;
        }
        if (!inverse && transform.mirrorSide) {
            x = -x;
        }
        return x;
    };
    auto transformCoordinates = [&](std::vector<double> coordinates, const SectionWindowTransform& transform) {
        coordinates[0] = transformX(coordinates[0], transform);
        return coordinates;
    };

    for (auto& section : sections) {
        section.set_coordinates(transformCoordinates(section.get_coordinates(), transformPerSection[section.get_name()]));
    }
    set_sections_description(sections);

    if (get_layers_description()) {
        auto layers = get_layers_description().value();
        for (auto& layer : layers) {
            if (!layer.get_section() || !transformPerSection.contains(layer.get_section().value())) {
                continue;
            }
            auto& transform = transformPerSection[layer.get_section().value()];
            layer.set_coordinates(transformCoordinates(layer.get_coordinates(), transform));
            if (layer.get_additional_coordinates()) {
                auto additionalCoordinates = layer.get_additional_coordinates().value();
                for (auto& coordinates : additionalCoordinates) {
                    coordinates = transformCoordinates(coordinates, transform);
                }
                layer.set_additional_coordinates(additionalCoordinates);
            }
        }
        set_layers_description(layers);
    }

    if (get_turns_description()) {
        auto turns = get_turns_description().value();
        for (auto& turn : turns) {
            if (!turn.get_section() || !transformPerSection.contains(turn.get_section().value())) {
                continue;
            }
            auto& transform = transformPerSection[turn.get_section().value()];
            turn.set_coordinates(transformCoordinates(turn.get_coordinates(), transform));
            if (transform.reflectAcrossWindowCenter) {
                // The reflection changed the turn's radius around its column: recompute
                // the length in the winding frame (absolute coordinates).
                auto frame = get_wound_column_frame_for_section(turn.get_section().value());
                auto turnLength = get_turn_length_in_frame(frame, std::abs(turn.get_coordinates()[0]));
                if (!turnLength) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                        "Negative turn length after placing turn " + turn.get_name() + " around its column");
                }
                turn.set_length(turnLength.value());
            }
            if (transform.mirrorSide) {
                // A mirrored turn is wound the opposite way around its column.
                if (turn.get_orientation()) {
                    turn.set_orientation(turn.get_orientation().value() == TurnOrientation::CLOCKWISE
                                             ? TurnOrientation::COUNTER_CLOCKWISE
                                             : TurnOrientation::CLOCKWISE);
                }
            }
            if (emitBothCrossings) {
                // Second crossing of the turn with the drawing plane: the reflection of
                // the first crossing across the wound column's axis (axis 0 for the
                // main column: the far side of the center leg; the leg axis for
                // lateral columns: outside the core).
                double secondCrossingX = 2 * transform.finalColumnAxisX - turn.get_coordinates()[0];
                turn.set_additional_coordinates(std::vector<std::vector<double>>{
                    {secondCrossingX, turn.get_coordinates()[1]}});
            }
        }
        set_turns_description(turns);
    }
    _groupWindowSidesApplied = !inverse;
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
    split_shared_window_groups(groups, windingWindows);
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
    // Fresh sections are wound in the +x winding frame.
    _groupWindowSidesApplied = false;

    std::vector<size_t> maybeVirtualizedPattern = pattern;
    std::vector<double> maybeVirtualizedProportionPerWinding = proportionPerWinding;
    auto functionalDescription = get_functional_description();
    auto needsVirtualization = needs_virtualization();

    std::optional<std::vector<Group>> originalGroupsDescription;
    if (needsVirtualization) {
        create_virtualization_map();
        auto virtualFunctionalDescription = virtualize_functional_description();
        maybeVirtualizedPattern = virtualize_pattern(pattern);
        maybeVirtualizedProportionPerWinding = virtualize_proportion_per_winding(proportionPerWinding);
        // Multi-window default groups reference the REAL winding names; remap
        // them to their virtual (group-leader) names BEFORE swapping in the
        // virtual description — the group loop's name lookup throws on merged
        // windings otherwise. Restored together with the description below.
        if (get_groups_description()) {
            originalGroupsDescription = get_groups_description();
            auto virtualGroups = originalGroupsDescription.value();
            for (auto& group : virtualGroups) {
                std::vector<PartialWinding> virtualPartialWindings;
                std::set<std::string> seenVirtualNames;
                for (auto partialWinding : group.get_partial_windings()) {
                    size_t realIndex = get_winding_index_by_name(get_functional_description(), partialWinding.get_winding());
                    auto leaderName = get_functional_description()[get_winding_group_minimum_index(realIndex)].get_name();
                    partialWinding.set_winding(leaderName);
                    if (seenVirtualNames.insert(leaderName).second) {
                        virtualPartialWindings.push_back(partialWinding);
                    }
                }
                group.set_partial_windings(virtualPartialWindings);
            }
            set_groups_description(virtualGroups);
        }
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
        if (originalGroupsDescription) {
            set_groups_description(originalGroupsDescription);
        }
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
    for (const auto& section : sections) {
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
    for (const auto& layer : layers) {
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
    for (const auto& section : sections) {
        auto newSection = virtualize_section(section);
        newSectionsDescription.push_back(newSection);
    }
    return newSectionsDescription;
}

std::vector<Layer> Coil::virtualize_layers_description() {
    std::vector<Layer> newLayersDescription;
    auto layers = get_layers_description().value();
    for (const auto& layer : layers) {
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

Coil::SectionGroupPlan Coil::plan_section_group(Group group, const std::vector<double>& proportionPerWinding,
                                                const std::vector<size_t>& pattern, size_t repetitions,
                                                bool multiGroup, WindingWindowShape windowShape,
                                                size_t conductionSectionOffset) {
    SectionGroupPlan plan;
    std::vector<size_t> groupPattern = pattern;
    std::vector<double> groupProportionPerWinding = proportionPerWinding;
    if (multiGroup) {
        for (auto& groupPartialWinding : group.get_partial_windings()) {
            plan.groupWindingIndexes.insert(get_winding_index_by_name(groupPartialWinding.get_winding()));
        }
        if (plan.groupWindingIndexes.empty()) {
            plan.skip = true;
            plan.group = group;
            return plan;
        }
        groupPattern.clear();
        for (auto windingIndex : pattern) {
            if (plan.groupWindingIndexes.contains(windingIndex)) {
                groupPattern.push_back(windingIndex);
            }
        }
        if (groupPattern.empty()) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                "Group " + group.get_name() + " has windings but none of them appear in the winding pattern");
        }
        // The group's windings share this group's full window space.
        double groupProportionSum = 0;
        for (auto windingIndex : plan.groupWindingIndexes) {
            groupProportionSum += proportionPerWinding[windingIndex];
        }
        for (auto windingIndex : plan.groupWindingIndexes) {
            groupProportionPerWinding[windingIndex] = proportionPerWinding[windingIndex] / groupProportionSum;
        }
        if (windowShape == WindingWindowShape::RECTANGULAR) {
            plan.groupWindowIndex = find_window_index_for_group(group.get_name());
            // Wind in the +x window-local frame; mirrored back afterwards.
            if (group.get_coordinates()[0] < 0) {
                auto groupCoordinates = group.get_coordinates();
                groupCoordinates[0] = -groupCoordinates[0];
                group.set_coordinates(groupCoordinates);
            }
        }
    }
    plan.group = group;

    // The section-stacking axis: rectangular groups carry a per-group orientation,
    // toroidal groups follow the coil-global winding orientation.
    double spaceForSections = 0;
    if (windowShape == WindingWindowShape::RECTANGULAR) {
        auto windingOrientation = group.get_sections_orientation();
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            spaceForSections = group.get_dimensions()[0];
        }
        else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
            spaceForSections = group.get_dimensions()[1];
        }
    }
    else {
        auto windingOrientation = get_winding_orientation();
        if (windingOrientation == WindingOrientation::OVERLAPPING) {
            spaceForSections = group.get_dimensions()[0];
        }
        else {
            spaceForSections = group.get_dimensions()[1];
        }
    }

    auto orderedSections = get_ordered_sections(spaceForSections, groupProportionPerWinding, groupPattern, repetitions);

    if (windowShape != WindingWindowShape::RECTANGULAR && get_winding_orientation() == WindingOrientation::CONTIGUOUS) {
        remove_insulation_if_margin_is_enough(orderedSections, conductionSectionOffset);
    }

    plan.orderedSectionsWithInsulation = add_insulation_to_sections(orderedSections);

    size_t numberWindings = get_functional_description().size();
    plan.numberSectionsPerWinding = std::vector<size_t>(numberWindings, 0);
    for (const auto& orderedSection : plan.orderedSectionsWithInsulation) {
        if (orderedSection.first == ElectricalType::CONDUCTION) {
            plan.numberSectionsPerWinding[orderedSection.second.first]++;
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
        if (plan.numberSectionsPerWinding[wIdx] != 0) continue;
        if (multiGroup && !plan.groupWindingIndexes.contains(wIdx)) continue;
        const auto& wwOpt = get_functional_description()[wIdx].get_wound_with();
        if (!wwOpt || wwOpt->empty()) continue;
        for (const auto& partnerName : wwOpt.value()) {
            bool found = false;
            for (size_t pIdx = 0; pIdx < numberWindings; ++pIdx) {
                if (get_functional_description()[pIdx].get_name() == partnerName &&
                    plan.numberSectionsPerWinding[pIdx] > 0) {
                    plan.numberSectionsPerWinding[wIdx] = plan.numberSectionsPerWinding[pIdx];
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
            if (plan.numberSectionsPerWinding[wIdx] != 0) continue;
            if (multiGroup && !plan.groupWindingIndexes.contains(wIdx)) continue;
            for (const auto& [virtualIdx, members] : _virtualizationMap) {
                if (std::find(members.begin(), members.end(), wIdx) == members.end()) continue;
                for (auto pIdx : members) {
                    if (pIdx < numberWindings && plan.numberSectionsPerWinding[pIdx] > 0) {
                        plan.numberSectionsPerWinding[wIdx] = plan.numberSectionsPerWinding[pIdx];
                        break;
                    }
                }
                if (plan.numberSectionsPerWinding[wIdx] != 0) break;
            }
        }
    }

    if (multiGroup) {
        // Windings belonging to other groups have zero sections here; give them a
        // placeholder slot count so the style chooser's zero-slot guard doesn't
        // fire. Their style entry is never read: only this group's windings appear
        // in this group's ordered sections.
        for (size_t wIdx = 0; wIdx < numberWindings; ++wIdx) {
            if (!plan.groupWindingIndexes.contains(wIdx) && plan.numberSectionsPerWinding[wIdx] == 0) {
                plan.numberSectionsPerWinding[wIdx] = 1;
            }
        }
    }

    plan.windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(), get_number_parallels(), plan.numberSectionsPerWinding);

    apply_margin_tape(plan.orderedSectionsWithInsulation, conductionSectionOffset);
    // ABT #721 (owner ruling 2026-08-14): coilEqualizeMargins applies to rectangular
    // windows too — with rectangular-aware semantics (same-edge pairing, no wrap; see
    // equalize_margins). Enabled by default; tests that pin the old 50/50 geometry
    // disable the setting explicitly.
    if (settings.get_coil_equalize_margins()) {
        equalize_margins(plan.orderedSectionsWithInsulation, conductionSectionOffset);
    }

    return plan;
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
    // ABT #720: margins are keyed by CONDUCTION-section ordinal, flat across groups in
    // wound order. This offset counts the conduction sections of the groups already wound.
    size_t conductionSectionOffset = 0;

    for (auto& groupEntry : groups) {
        // ABT #720: shared orchestration (pattern subsetting, proportion renormalization,
        // ordered sections + insulation, slot inheritance, winding styles, margin
        // application) lives in plan_section_group — one copy for both winders.
        auto plan = plan_section_group(groupEntry, proportionPerWinding, pattern, repetitions, multiGroup,
                                       WindingWindowShape::RECTANGULAR, conductionSectionOffset);
        if (plan.skip) {
            continue;
        }
        auto group = plan.group;
        const auto& orderedSectionsWithInsulation = plan.orderedSectionsWithInsulation;
        const auto& numberSectionsPerWinding = plan.numberSectionsPerWinding;
        const auto& windByConsecutiveTurns = plan.windByConsecutiveTurns;
        std::optional<size_t> groupWindowIndex = plan.groupWindowIndex;
        auto currentSectionPerWinding = std::vector<size_t>(get_functional_description().size(), 0);

        double availableWidth = group.get_dimensions()[0];
        double availableHeight = group.get_dimensions()[1];
        auto windingOrientation = group.get_sections_orientation();

        auto wirePerWinding = get_wires();
        // ABT #720: explicit first-section flag instead of DBL_MAX in-band sentinels.
        bool sectionCentersInitialized = false;
        double currentSectionCenterWidth = 0;
        double currentSectionCenterHeight = 0;

        // Margins are keyed by conduction ordinal, flat across ALL groups in wound order.
        size_t conductionOrdinal = conductionSectionOffset;
        for (size_t sectionIndex = 0; sectionIndex < orderedSectionsWithInsulation.size(); ++sectionIndex) {
            if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
                size_t marginIndex = conductionOrdinal++;
                auto sectionInfo = orderedSectionsWithInsulation[sectionIndex].second;
                auto windingIndex = sectionInfo.first;
                auto spaceForSection = sectionInfo.second;

                double currentSectionHeight = 0;
                double currentSectionWidth = 0;
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionWidth = spaceForSection;
                    currentSectionHeight = availableHeight;
                    if (!sectionCentersInitialized) {
                        currentSectionCenterWidth = group.get_coordinates()[0] - availableWidth / 2;
                        currentSectionCenterHeight = group.get_coordinates()[1];
                        sectionCentersInitialized = true;
                    }
                }
                else {
                    currentSectionWidth = availableWidth;
                    currentSectionHeight = spaceForSection;
                    if (!sectionCentersInitialized) {
                        currentSectionCenterWidth = group.get_coordinates()[0];
                        currentSectionCenterHeight = group.get_coordinates()[1] + availableHeight / 2;
                        sectionCentersInitialized = true;
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

                section.set_name(get_name(windingIndex) +  " section " + std::to_string(currentSectionPerWinding[windingIndex]));
                section.set_partial_windings(std::vector<PartialWinding>{partialWinding});  // TODO: support more than one winding per section?
                section.set_group(group.get_name());
                if (groupWindowIndex) {
                    section.set_winding_window(static_cast<int64_t>(groupWindowIndex.value()));
                }
                section.set_type(ElectricalType::CONDUCTION);
                section.set_margin(_marginsPerSection[marginIndex]);
                section.set_layers_orientation(get_layers_orientation(section.get_name()));
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

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterWidth += insulationSection.get_dimensions()[0];
                }
                else {
                    currentSectionCenterHeight -= insulationSection.get_dimensions()[1];
                }

            }
        }
        conductionSectionOffset = conductionOrdinal;
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

void Coil::remove_insulation_if_margin_is_enough(const std::vector<std::pair<size_t, double>>& orderedSections, size_t conductionSectionOffset) {
    auto bobbin = resolve_bobbin();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

    // ABT #720: _marginsPerSection is keyed by CONDUCTION-section ordinal, and orderedSections
    // here is conduction-only — the index IS the group-local ordinal, offset by the previous
    // groups' conduction count. (The old flat keying forced a runtime "multiplier" inference of
    // this function's own indexing, and the missing offset made every group after the first
    // read the first group's margin rows.)
    if (_marginsPerSection.size() < conductionSectionOffset + orderedSections.size()) {
        _marginsPerSection.resize(conductionSectionOffset + orderedSections.size(), {0, 0});
    }

    for (size_t sectionIndex = 0; sectionIndex < orderedSections.size(); ++sectionIndex) {
        size_t indexForMarginLeftSection = conductionSectionOffset + sectionIndex;
        size_t indexForMarginRightSection;
        size_t leftWindingIndex = orderedSections[sectionIndex].first;
        size_t rightWindingIndex;
        if (sectionIndex + 1 != orderedSections.size()) {
            indexForMarginRightSection = conductionSectionOffset + sectionIndex + 1;
            rightWindingIndex = orderedSections[sectionIndex + 1].first;
        }
        else {
            // The toroidal window closes on itself: the last section's neighbour wraps to
            // this group's first.
            indexForMarginRightSection = conductionSectionOffset;
            rightWindingIndex = orderedSections[0].first;
        }

        auto windingsMapKey = std::pair<size_t, size_t>{leftWindingIndex, rightWindingIndex};
        double totalMargin = 0;
        if (_insulationSections.contains(windingsMapKey)) {
            // find, not operator[]: the custom-insulation path fills _insulationSections
            // without interfaces, and operator[] here used to insert a junk default
            // interface (uninitialized layerPurpose) as a side effect.
            auto coilSectionInterfaceIt = _coilSectionInterfaces.find(windingsMapKey);
            if (coilSectionInterfaceIt != _coilSectionInterfaces.end()) {
                totalMargin = coilSectionInterfaceIt->second.get_total_margin_tape_distance();
            }
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
    // Filled ONCE at function scope, like the rectangular sibling: it is indexed by global
    // winding index, and until 2026-08 each group iteration appended a fresh set of
    // entries, so on a multi-group toroid the second group kept consuming the first
    // group's already-decremented rows while the appended tail was never validated.
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        remainingParallelsProportion.push_back(std::vector<double>(get_number_parallels(windingIndex), 1));
    }
    // Multi-column winding (rectangular sibling's rationale applies here unchanged, ABT
    // #227.6): each group is wound with only the subset of the pattern that belongs to
    // its own windings. Without this, a POLAR winding window with more than one group
    // (create_default_groups emits one per toroidal winding window) fed the FULL global
    // pattern/proportionPerWinding into every group's get_ordered_sections call, so a
    // group ended up trying to place windings that do not belong to it.
    bool multiGroup = groups.size() > 1;
    // ABT #720: margins are keyed by CONDUCTION-section ordinal, flat across ALL groups in
    // wound order (rectangular sibling's convention); without the offset the second group's
    // apply_margin_tape re-indexed from 0 and clobbered the first group's margins.
    size_t conductionSectionOffset = 0;

    for (auto& groupEntry : groups) {
        // ABT #720: shared orchestration — see plan_section_group.
        auto plan = plan_section_group(groupEntry, proportionPerWinding, pattern, repetitions, multiGroup,
                                       WindingWindowShape::ROUND, conductionSectionOffset);
        if (plan.skip) {
            continue;
        }
        auto group = plan.group;
        const auto& orderedSectionsWithInsulation = plan.orderedSectionsWithInsulation;
        const auto& numberSectionsPerWinding = plan.numberSectionsPerWinding;
        const auto& windByConsecutiveTurns = plan.windByConsecutiveTurns;
        auto currentSectionPerWinding = std::vector<size_t>(get_functional_description().size(), 0);

        double availableRadialHeight = group.get_dimensions()[0];
        double availableAngle = group.get_dimensions()[1];
        auto windingOrientation = get_winding_orientation();

        auto wirePerWinding = get_wires();
        // ABT #720: explicit first-section flag instead of DBL_MAX in-band sentinels.
        bool sectionCentersInitialized = false;
        double currentSectionCenterAngle = 0;
        double currentSectionCenterRadialHeight = 0;

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
            // ABT #720: margins are keyed by conduction ordinal, flat across ALL groups.
            size_t marginIndex = conductionSectionOffset + conductingSectionIndex;
            if (orderedSectionsWithInsulation[sectionIndex].first == ElectricalType::CONDUCTION) {
                auto sectionInfo = orderedSectionsWithInsulation[sectionIndex].second;
                auto windingIndex = sectionInfo.first;

                double currentSectionRadialHeight = currentSectionRadialHeights[conductingSectionIndex];
                double currentSectionAngle = currentSectionAngles[conductingSectionIndex];

                if (!sectionCentersInitialized) {
                    currentSectionCenterRadialHeight = 0;
                    currentSectionCenterAngle = windingOrientation == WindingOrientation::OVERLAPPING ? 180 : 0;
                    sectionCentersInitialized = true;
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

                // ABT #187: real-winding angular blocking — size the section for the ring capacities
                // AFTER the blocked slots are removed, so spilled turns get their radial space. The
                // section name is assigned below with this same deterministic convention.
                //
                // ABT #723: include the STATIC input-connection corridor for rings after the first
                // under the same gate as wind_by_round_layers and the wind_by_round_turns placement
                // shift — all three must agree or the section is sized for fewer rings than the
                // layer split later produces and the wind fails with no turns (the escalated-spill
                // U-order re-wind hit exactly that).
                std::vector<int64_t> blockedSlotsPerLayerIndex;
                if (settings.get_coil_use_real_winding_geometry()) {
                    std::string futureSectionName = get_name(windingIndex) + " section " + std::to_string(currentSectionPerWinding[windingIndex]);
                    int64_t inputCorridorSlots = int64_t(get_number_parallels(windingIndex)) + 1;
                    for (size_t k = 0; k < 64; ++k) {
                        int64_t slots = (k >= 1) ? inputCorridorSlots : 0;
                        if (_applyConnectionBlocking) {
                            auto found = _connectionBlockedSlotsPerLayer.find(futureSectionName + " layer " + std::to_string(k));
                            if (found != _connectionBlockedSlotsPerLayer.end()) {
                                slots += int64_t(found->second.first);
                            }
                        }
                        blockedSlotsPerLayerIndex.push_back(slots);
                    }
                }
                const std::vector<int64_t>* blockedSlotsPointer = blockedSlotsPerLayerIndex.empty() ? nullptr : &blockedSlotsPerLayerIndex;

                // We correct the radial height to exactly what we need, so afterwards we can calculate exactly how many turns we need
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    auto aux = get_number_layers_needed_and_number_physical_turns(currentSectionCenterRadialHeight + _marginsPerSection[marginIndex][0], currentSectionAngle, wirePerWinding[windingIndex], physicalTurnsThisSection, availableRadialHeight, blockedSlotsPointer);
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
                        auto aux = get_number_layers_needed_and_number_physical_turns(currentSectionCenterRadialHeight, currentSectionAngleMinusMargin, wirePerWinding[windingIndex], physicalTurnsThisSection, availableRadialHeight, blockedSlotsPointer);
                        numberLayers = aux.first;
                        if (_strict) {
                            currentSectionRadialHeight = numberLayers * wirePerWinding[windingIndex].get_maximum_outer_width();
                        }
                        double lastLayerMaximumRadius = availableRadialHeight - (currentSectionCenterRadialHeight + numberLayers * wirePerWinding[windingIndex].get_maximum_outer_width());
                        if (lastLayerMaximumRadius < 0) {
                            break;
                        }
                        marginAngle0 = wound_distance_to_angle(_marginsPerSection[marginIndex][0], lastLayerMaximumRadius);
                        marginAngle1 = wound_distance_to_angle(_marginsPerSection[marginIndex][1], lastLayerMaximumRadius);
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
                section.set_margin(_marginsPerSection[marginIndex]);
                section.set_layers_orientation(get_layers_orientation(section.get_name()));
                section.set_coordinate_system(CoordinateSystem::POLAR);
                
                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    section.set_dimensions(std::vector<double>{currentSectionRadialHeight, currentSectionAngle});
                    section.set_coordinates(std::vector<double>{currentSectionCenterRadialHeight + currentSectionRadialHeight / 2 + _marginsPerSection[marginIndex][0], currentSectionCenterAngle, 0});
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
                    currentSectionCenterRadialHeight += currentSectionRadialHeight + _marginsPerSection[marginIndex][0] + _marginsPerSection[marginIndex][1];
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

                if (windingOrientation == WindingOrientation::OVERLAPPING) {
                    currentSectionCenterRadialHeight += insulationSection.get_dimensions()[0];
                }
                else {
                    currentSectionCenterAngle += insulationSection.get_dimensions()[1];
                }
            }
        }
        conductionSectionOffset += conductingSectionIndex;
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

        WindingStyle windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberSections, windingIndex);

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
            // section. Unconstrained-layer-count only. Works for any parallel count: a
            // bifilar/N-filar winding lays its parallels side by side, so each layer holds an equal
            // number of turns of every parallel (its physical capacity is rounded down to a multiple of
            // the parallel count) and the per-layer split is carried by CONSECUTIVE_PARALLELS below.
            // ABT #427: both layer orientations block. The turn axis — the direction a lead takes a slot
            // out of — is the layer's HEIGHT when layers overlap (turns stack axially) and its WIDTH
            // when they are contiguous (the layer is one wire tall, turns run along it), so the slot
            // pitch is the wire's outer height or outer width to match.
            int64_t numberParallels = int64_t(get_number_parallels(windingIndex));
            bool layersStackAlongWidth = (sections[sectionIndex].get_layers_orientation() == WindingOrientation::OVERLAPPING);
            bool realWindingBlocking = settings.get_coil_use_real_winding_geometry()
                && _applyConnectionBlocking
                && !sections[sectionIndex].get_number_layers()
                && maximumNumberPhysicalTurnsPerLayer > 1
                && maximumNumberPhysicalTurnsPerLayer >= uint64_t(numberParallels);
            // INPUT-CONNECTION ANGULAR BLOCKER (toroids, Alf's rule -- the round-window analog of
            // the height reserved for connections on concentric windings): the connection enters
            // at the section's start extreme on the FIRST ring, so every SUBSEQUENT ring
            // surrenders the angular slots the connection's parallels occupy on that edge, and
            // -- since a deeper ring's circumference shrinks -- its capacity is measured at ITS
            // OWN radius, not the bore rim's (the uniform rim capacity is what let outer rings
            // overhang the first ring on both edges). The spread placement applies the matching
            // angular shift (wind_by_round_turns).
            auto resolvedBobbinForCorridor = resolve_bobbin();
            bool toroidalConnectionCorridor = settings.get_coil_use_real_winding_geometry()
                && resolvedBobbinForCorridor.get_winding_window_shape() == WindingWindowShape::ROUND
                && !sections[sectionIndex].get_number_layers()
                && maximumNumberPhysicalTurnsPerLayer > 1
                && maximumNumberPhysicalTurnsPerLayer >= uint64_t(numberParallels);
            double wireTurnAxisSize = layersStackAlongWidth
                ? wirePerWinding[windingIndex].get_maximum_outer_height()
                : wirePerWinding[windingIndex].get_maximum_outer_width();
            auto blockedSlotsForLayer = [&](size_t layerIndexInSection) -> std::pair<uint64_t, uint64_t> {
                auto found = _connectionBlockedSlotsPerLayer.find(
                    sections[sectionIndex].get_name() + " layer " + std::to_string(layerIndexInSection));
                if (found == _connectionBlockedSlotsPerLayer.end()) {
                    return {0u, 0u};
                }
                return found->second;
            };
            std::vector<uint64_t> realPerLayerTurns;
            if (realWindingBlocking || toroidalConnectionCorridor) {
                uint64_t remainingTurns = physicalTurnsInSection;
                size_t builtLayers = 0;
                while (remainingTurns > 0) {
                    auto blocked = blockedSlotsForLayer(builtLayers);
                    uint64_t blockedSlots = std::min<uint64_t>(blocked.first + blocked.second, maximumNumberPhysicalTurnsPerLayer - 1);
                    uint64_t capacity = maximumNumberPhysicalTurnsPerLayer - blockedSlots;
                    if (realWindingBlocking && !toroidalConnectionCorridor) {
                        // ABT #616 (Alf, 26_psps: "3 turns should fit in layer 2 3 and 4"): the
                        // slot currency CEILs each edge's run depth to whole turn slots, so a
                        // 25 um insulation overhang on a two-row stack (0.885 mm) costs a whole
                        // third slot and a layer that continuously fits 3 turns is filled with 2.
                        // The placement (align_blocked_layer_turns) already lets turns hug the
                        // runs at their CONTINUOUS depth — size the capacity the same way. The
                        // U-landing overlay stays uncharged (placement-only, per the N_layer+1
                        // retraction).
                        const std::string layerKey = sections[sectionIndex].get_name() +
                                                     " layer " + std::to_string(builtLayers);
                        std::pair<double, double> depths{0.0, 0.0};
                        auto foundDepth = _connectionBlockedDepthPerLayer.find(layerKey);
                        if (foundDepth != _connectionBlockedDepthPerLayer.end()) {
                            depths = foundDepth->second;
                        }
                        const double extent = layersStackAlongWidth
                            ? sections[sectionIndex].get_dimensions()[1]
                            : sections[sectionIndex].get_dimensions()[0];
                        const double freeBand = extent - depths.first - depths.second;
                        if (freeBand > 0.0 && wireTurnAxisSize > 0.0) {
                            uint64_t continuousCapacity =
                                uint64_t(std::floor(freeBand / wireTurnAxisSize + 1e-9));
                            capacity = std::max(capacity, continuousCapacity);
                            capacity = std::min<uint64_t>(capacity, maximumNumberPhysicalTurnsPerLayer);
                        }
                    }
                    if (toroidalConnectionCorridor) {
                        // Deeper rings hold fewer turns (same wire, smaller radius) and, past the
                        // first ring, surrender the connection corridor's slots.
                        auto windingWindowsForCapacity = resolvedBobbinForCorridor.get_processed_description().value().get_winding_windows();
                        double boreRadius = windingWindowsForCapacity[0].get_radial_height().value();
                        double ringRadiusFirst = boreRadius - wireTurnAxisSize / 2;
                        double ringRadiusThis = boreRadius - (double(builtLayers) + 0.5) * wireTurnAxisSize;
                        if (ringRadiusThis <= 0) {
                            break;   // past the toroid centre: nothing more fits
                        }
                        uint64_t radiusCapacity = uint64_t(std::floor(
                            double(maximumNumberPhysicalTurnsPerLayer) * ringRadiusThis / ringRadiusFirst));
                        capacity = std::min(capacity, radiusCapacity);
                        if (builtLayers >= 1) {
                            uint64_t corridorSlots = uint64_t(numberParallels) + 2;
                            capacity = capacity > corridorSlots ? capacity - corridorSlots : 0;
                        }
                    }
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
                if (std::getenv("MKF_BLOCKING_DIAG")) {
                    std::cerr << "[fill] " << sections[sectionIndex].get_name() << " turns="
                              << physicalTurnsInSection << " perLayer={";
                    for (auto t : realPerLayerTurns) std::cerr << t << ",";
                    std::cerr << "} blocked={";
                    for (size_t li = 0; li < realPerLayerTurns.size(); ++li) {
                        auto b = blockedSlotsForLayer(li);
                        std::cerr << b.first << "+" << b.second << ",";
                    }
                    std::cerr << "} maxPerLayer=" << maximumNumberPhysicalTurnsPerLayer << "\n";
                }
            }

            // ABT #229 root-cause fix: an N-filar winding laid CONSECUTIVE_PARALLELS holds its K
            // parallels side by side, so a layer's usable capacity is a whole number of K-turn rows
            // (floor(capacity / K) rows). numberLayers above was derived from the RAW capacity, so
            // whenever the per-parallel ceil split exceeds the rounded capacity the first layer
            // overflows (e.g. 3 parallels x 19 crossings in 20-slot layers: ceil(57/20)=3 layers ->
            // per-parallel 7/6/6 -> 21 > 20, filling factor > 1, and wind() then silently skipped
            // the whole real-winding blocking machinery because the coil "did not fit"). Grow the
            // layer count to what the rows actually need (the compaction pass shrinks the section to
            // its content afterwards, so an extra layer only needs window space, which the fitting
            // check still verifies via maximumNumberLayersFittingInSection).
            if (!realWindingBlocking && numberParallels > 1 && !sections[sectionIndex].get_number_layers()
                && maximumNumberPhysicalTurnsPerLayer > 0 && numberLayers > 0) {
                WindingStyle tentativeStyle = sections[sectionIndex].get_winding_style()
                    ? sections[sectionIndex].get_winding_style().value()
                    : wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberLayers);
                uint64_t rowsPerLayer = maximumNumberPhysicalTurnsPerLayer / uint64_t(numberParallels);
                if (tentativeStyle == WindingStyle::WIND_BY_CONSECUTIVE_PARALLELS && rowsPerLayer > 0) {
                    uint64_t perParallelTurns = uint64_t(std::ceil(double(physicalTurnsInSection) / double(numberParallels) - 1e-9));
                    uint64_t perLayerSplit = uint64_t(std::ceil(double(perParallelTurns) / double(numberLayers) - 1e-9));
                    if (perLayerSplit * uint64_t(numberParallels) > maximumNumberPhysicalTurnsPerLayer) {
                        uint64_t layersNeededForRows = uint64_t(std::ceil(double(perParallelTurns) / double(rowsPerLayer) - 1e-9));
                        numberLayers = std::max(numberLayers, layersNeededForRows);
                        // A zero maximumNumberLayersFittingInSection means "section thinner than one
                        // wire" and the sizing above treats it as unbounded — don't clamp to it.
                        if (maximumNumberLayersFittingInSection > 0) {
                            numberLayers = std::min(numberLayers, maximumNumberLayersFittingInSection);
                        }
                    }
                }
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
                windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberLayers, windingIndex);
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
            if ((realWindingBlocking || toroidalConnectionCorridor) && numberParallels == 1) {
                windByConsecutiveTurns = WindingStyle::WIND_BY_CONSECUTIVE_TURNS;
            }

            for (size_t layerIndex = 0; layerIndex < numberLayers; ++layerIndex) {
                Layer layer;

                auto parallelsProportions = (realWindingBlocking || toroidalConnectionCorridor)
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
                // leads crossing it: shrink the layer along its TURN axis by the blocked slots and
                // shift it away from those slots (high-side slots push it toward the low side and
                // vice versa) so the turns stop before the connection area, leaving room for the
                // leads. ABT #427: that axis is the layer's HEIGHT when layers overlap and its WIDTH
                // when they are contiguous, so this is one shrink driven by an axis flag rather than
                // two mirrored branches that can drift apart. The layer PITCH is untouched either way
                // — only the layer's own extent shrinks, so the stepping below is unaffected.
                double thisLayerHeight = layerHeight;
                double thisLayerWidth = layerWidth;
                double thisLayerCenterHeight = currentLayerCenterHeight;
                double thisLayerCenterWidth = currentLayerCenterWidth;
                if (realWindingBlocking) {
                    auto blocked = blockedSlotsForLayer(layerIndex);
                    uint64_t blockedSlots = std::min<uint64_t>(blocked.first + blocked.second, maximumNumberPhysicalTurnsPerLayer - 1);
                    double surrenderedRoom = blockedSlots * wireTurnAxisSize;
                    // blocked is {high side, low side} of the turn axis, so a high-side block moves the
                    // layer toward the low side (negative) and vice versa.
                    double shiftTowardLowSide = (double(blocked.second) - double(blocked.first)) * wireTurnAxisSize / 2;
                    if (layersStackAlongWidth) {
                        thisLayerHeight = roundFloat(layerHeight - surrenderedRoom, 9);
                        thisLayerCenterHeight = roundFloat(currentLayerCenterHeight + shiftTowardLowSide, 9);
                    }
                    else {
                        thisLayerWidth = roundFloat(layerWidth - surrenderedRoom, 9);
                        thisLayerCenterWidth = roundFloat(currentLayerCenterWidth + shiftTowardLowSide, 9);
                    }
                    // ABT #430: record the room actually surrendered here — AFTER the one-slot-minimum
                    // cap above, so it is what the layer really gave up and not what was asked of it.
                    // apply_connection_reserved_space subtracts this from the leads it charges, because
                    // the shrunken extent (and the filling factor computed from it just below) already
                    // excludes it; charging the full lead extent again counted the same room twice.
                    _connectionBlockedRoomPerLayer[layer.get_name()] = surrenderedRoom;
                }
                layer.set_dimensions(std::vector<double>{thisLayerWidth, thisLayerHeight});
                layer.set_coordinates(std::vector<double>{thisLayerCenterWidth, thisLayerCenterHeight, 0});
                layer.set_coordinate_system(CoordinateSystem::CARTESIAN);

                layer.set_filling_factor(get_area_used_in_wires(wirePerWinding[windingIndex], physicalTurnsThisLayer) / (thisLayerWidth * thisLayerHeight));
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
                    // ABT #725: the template stores its THICKNESS in the dims slot named by its
                    // OWN orientation (stamped from the global at set_interlayer_insulation
                    // time). Reading the slot chosen by this section's — possibly overridden —
                    // orientation read the full window dimension as "thickness" whenever the two
                    // disagreed: a window-sized insulation sheet inserted inside the section.
                    double insulationThickness = insulationLayer.get_orientation() == WindingOrientation::CONTIGUOUS
                                                     ? insulationLayer.get_dimensions()[1]
                                                     : insulationLayer.get_dimensions()[0];
                    if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                        currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - insulationThickness / 2, 9);
                    }
                    else {
                        currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + insulationThickness / 2, 9);
                    }

                    insulationLayer.set_coordinate_system(CoordinateSystem::CARTESIAN);
                    insulationLayer.set_section(sections[sectionIndex].get_name());
                    insulationLayer.set_name(sections[sectionIndex].get_name() +  " insulation layer " + std::to_string(layerIndex));
                    insulationLayer.set_coordinates(std::vector<double>{currentLayerCenterWidth, currentLayerCenterHeight, 0});
                    // The inserted sheet lives in this section: its orientation and dims follow
                    // the section's layers orientation, not the template's global one.
                    insulationLayer.set_orientation(sections[sectionIndex].get_layers_orientation());
                    if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                        insulationLayer.set_dimensions(std::vector<double>{layerWidth, insulationThickness});
                    }
                    else {
                        insulationLayer.set_dimensions(std::vector<double>{insulationThickness, layerHeight});
                    }
                    layers.push_back(insulationLayer);

                    if (sections[sectionIndex].get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                        currentLayerCenterHeight = roundFloat(currentLayerCenterHeight - insulationThickness / 2, 9);
                    }
                    else {
                        currentLayerCenterWidth = roundFloat(currentLayerCenterWidth + insulationThickness / 2, 9);
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

            // ABT #187: real-winding angular blocking — leads crossing a ring reserve turn slots on
            // it, so the ring capacity computation subtracts them (turns spill to deeper rings).
            //
            // ABT #723: ALSO charge, from the very FIRST wind, the input-connection corridor that
            // wind_by_round_turns unconditionally shaves off every ring after the first
            // ((parallels + 1) x wireAngle at the section-start edge). The placement-side shift
            // has always applied whenever real winding is on, but this capacity charge only
            // existed in an unreachable branch of wind_by_rectangular_layers — so deeper rings
            // were handed turn counts sized for the FULL span and packed straight back through
            // the shaved corridor (the section-contiguous 2-ring toroid's turn 29 conflict).
            // Placement and capacity must share one gate or the turns overflow the shifted span.
            std::vector<int64_t> blockedSlotsPerLayerIndex;
            if (settings.get_coil_use_real_winding_geometry()) {
                int64_t inputCorridorSlots = int64_t(get_number_parallels(windingIndex)) + 1;
                for (size_t k = 0; k < 64; ++k) {
                    int64_t slots = (k >= 1) ? inputCorridorSlots : 0;
                    if (_applyConnectionBlocking) {
                        auto found = _connectionBlockedSlotsPerLayer.find(sections[sectionIndex].get_name() + " layer " + std::to_string(k));
                        if (found != _connectionBlockedSlotsPerLayer.end()) {
                            slots += int64_t(found->second.first);
                        }
                    }
                    blockedSlotsPerLayerIndex.push_back(slots);
                }
            }
            const std::vector<int64_t>* blockedSlotsPointer = blockedSlotsPerLayerIndex.empty() ? nullptr : &blockedSlotsPerLayerIndex;

            if (maximumNumberLayersFittingInSection == 0) {
                auto aux = get_number_layers_needed_and_number_physical_turns(sections[sectionIndex], wirePerWinding[windingIndex], physicalTurnsInSection, windingWindowRadialHeight, blockedSlotsPointer);
                numberLayers = aux.first;
                layerPhysicalTurns = aux.second;
            }
            else if (maximumNumberPhysicalTurnsPerLayer == 0) {
                auto aux = get_number_layers_needed_and_number_physical_turns(sections[sectionIndex], wirePerWinding[windingIndex], physicalTurnsInSection, windingWindowRadialHeight, blockedSlotsPointer);
                numberLayers = maximumNumberLayersFittingInSection;
                layerPhysicalTurns = aux.second;
            }
            else {
                auto aux = get_number_layers_needed_and_number_physical_turns(sections[sectionIndex], wirePerWinding[windingIndex], physicalTurnsInSection, windingWindowRadialHeight, blockedSlotsPointer);
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
                windByConsecutiveTurns = wind_by_consecutive_turns(get_number_turns(windingIndex), get_number_parallels(windingIndex), numberLayers, windingIndex);
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

    for (const auto& section : sections) {
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
    // ABT #615 stage 2: each winding's k-th SECTION alternates its winding direction (Alf,
    // 2026-08-09: the connection arrives on top, so the receiving section "should connect to the
    // topmost turn, and from there start the Z winding the other way around: from top to bottom in
    // the layer, and back to top with a dragback" -- and invert again on the next hop).
    std::map<std::string, int64_t> windingSectionOrdinal;
    std::map<std::string, std::string> windingLastSection;
    std::vector<Turn> turns;
    for (auto& layer : layers) {
        if (layer.get_type() == ElectricalType::CONDUCTION) {
            double currentTurnCenterWidth = 0;
            double currentTurnCenterHeight = 0;
            double currentTurnWidthIncrement = 0;
            double currentTurnHeightIncrement = 0;
            double totalLayerHeight;
            double totalLayerWidth;
            // SPREAD lays its turns on explicit fence-post stations (ABT #578/#579) rather than on a
            // uniform increment: bundle members touch while only the gaps BETWEEN bundles carry the
            // slack, so consecutive steps differ and no single increment can express them. Left empty
            // by every other alignment, which keeps its uniform-increment placement untouched (those
            // alignments pack turns at exactly one wire pitch, so their bundles already touch).
            std::vector<double> turnStations;
            size_t turnStationAxis = 0;  // 0 = width, 1 = height; only read when turnStations is set
            size_t turnStationIndex = 0;
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
                        // Fence-post, bundle-aware stations (ABT #578/#579) instead of a uniform
                        // height/N increment: the step is not constant, so it cannot be carried in
                        // currentTurnHeightIncrement. Stations come back ascending; placement here
                        // runs top-down, so reverse them.
                        turnStations = compute_spread_turn_stations(layer.get_coordinates()[1],
                                                                    layer.get_dimensions()[1],
                                                                    wireHeight,
                                                                    physicalTurnsInLayer,
                                                                    get_layer_bundle_size(layer));
                        std::reverse(turnStations.begin(), turnStations.end());
                        turnStationAxis = 1;
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
                        // Mirror of the OVERLAPPING branch above (ABT #578/#579), on the width axis.
                        // Placement here runs left-to-right, which is already the stations' ascending
                        // order, so unlike the height branch these are NOT reversed.
                        turnStations = compute_spread_turn_stations(layer.get_coordinates()[0],
                                                                    layer.get_dimensions()[0],
                                                                    wireWidth,
                                                                    physicalTurnsInLayer,
                                                                    get_layer_bundle_size(layer));
                        turnStationAxis = 0;
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
            if (layer.get_section()
                && windingLastSection[windingNameForOrder] != layer.get_section().value()) {
                if (!windingLastSection[windingNameForOrder].empty()) {
                    windingSectionOrdinal[windingNameForOrder]++;
                }
                windingLastSection[windingNameForOrder] = layer.get_section().value();
            }
            // Convention: every winding STARTS FROM THE BOTTOM. The alignment above lays the first turn
            // at the top, so reverse it to start at the bottom and wind up — for every layer EXCEPT:
            //   - the odd ordinals of a U winding, whose boustrophedon turnaround means they start at
            //     the top (so consecutive layers stay adjacent and the connection stays on top);
            //   - ABT #615 stage 2 (real winding): the odd SECTIONS of a winding's interleave chain.
            //     The inter-section connection arrives at the edge the previous section finished on,
            //     so the receiving section starts at its TOPMOST turn and winds the other way (top to
            //     bottom per layer, dragbacks returning bottom to top), inverting again per hop.
            bool startFromTop = false;
            // ABT #616: under real winding, the alternation chain is BASED at the edge the
            // winding's entrance terminal row occupies (recorded by the previous blocking
            // iteration) — the first section/layer starts adjacent to its own connection,
            // and the U/Z alternation inverts from there.
            bool entranceBase = false;
            if (settings.get_coil_use_real_winding_geometry()
                && layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                auto foundEdge = _terminalEntranceAtTop.find(windingNameForOrder);
                if (foundEdge != _terminalEntranceAtTop.end()) {
                    entranceBase = foundEdge->second;
                }
            }
            if (get_winding_order(layer.get_section().value()) == WindingOrder::U) {
                startFromTop = entranceBase != (windingLayerOrdinal % 2 == 1);
            }
            else if (settings.get_coil_use_real_winding_geometry()) {
                startFromTop = entranceBase != (windingSectionOrdinal[windingNameForOrder] % 2 == 1);
            }
            if (!startFromTop) {
                currentTurnCenterWidth = roundFloat(currentTurnCenterWidth + (int64_t(physicalTurnsInLayer) - 1) * currentTurnWidthIncrement, 9);
                currentTurnCenterHeight = roundFloat(currentTurnCenterHeight - (int64_t(physicalTurnsInLayer) - 1) * currentTurnHeightIncrement, 9);
                currentTurnWidthIncrement = -currentTurnWidthIncrement;
                currentTurnHeightIncrement = -currentTurnHeightIncrement;
                // The increment arithmetic above reverses a uniform run by starting from its last
                // position and negating the step; SPREAD's stations are not uniform, so they reverse
                // by reversing the list. (The arithmetic is a no-op for SPREAD: its increments are 0.)
                std::reverse(turnStations.begin(), turnStations.end());
            }

            if (!layer.get_winding_style()) {
                layer.set_winding_style(WindingStyle::WIND_BY_CONSECUTIVE_TURNS);
            }


            if (layer.get_winding_style().value() == WindingStyle::WIND_BY_CONSECUTIVE_TURNS) {
                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    int64_t numberTurns = round(partialWinding.get_parallels_proportion()[parallelIndex] * get_number_turns(windingIndex));
                    for (int64_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
                        take_next_spread_station(turnStations, turnStationIndex, turnStationAxis, layer.get_name(),
                                                 currentTurnCenterWidth, currentTurnCenterHeight);
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
                            take_next_spread_station(turnStations, turnStationIndex, turnStationAxis, layer.get_name(),
                                                     currentTurnCenterWidth, currentTurnCenterHeight);
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

            // INPUT-CONNECTION ANGULAR BLOCKER (real winding, toroids): rings after the first
            // surrender the connection corridor on the section-start (low-angle) edge -- the
            // entrance runs there on the first ring, and a later ring's turn placed behind it
            // is copper the final 3D cannot avoid crossing. Work in the REDUCED, SHIFTED span
            // for every alignment so no fallback path can reoccupy the corridor. Capacity was
            // already charged in wind_by_layers; here only the angular window moves.
            double layerAngularDimension = layer.get_dimensions()[1];
            double layerAngularCentre = layer.get_coordinates()[1];
            if (settings.get_coil_use_real_winding_geometry() &&
                layer.get_type() == ElectricalType::CONDUCTION) {
                const std::string& layerNameForCorridor = layer.get_name();
                auto marker = layerNameForCorridor.rfind(" layer ");
                if (marker != std::string::npos &&
                    std::stoul(layerNameForCorridor.substr(marker + 7)) >= 1) {
                    double corridorAngle =
                        (double(get_number_parallels(windingIndex)) + 1.0) * wireAngle;
                    corridorAngle = std::min(corridorAngle, layerAngularDimension / 2);
                    layerAngularDimension -= corridorAngle;
                    layerAngularCentre += corridorAngle / 2;
                }
            }

            if (layer.get_orientation() == WindingOrientation::OVERLAPPING) {
                totalLayerAngle = physicalTurnsInLayer * wireAngle;

                currentTurnRadialHeightIncrement = 0;
                currentTurnCenterRadialHeight = roundFloat(layer.get_coordinates()[0], 9);
                switch (alignment) {
                    case CoilAlignment::CENTERED:
                        currentTurnCenterAngle = roundFloat(layerAngularCentre - totalLayerAngle / 2 + wireAngle / 2, 9);
                        currentTurnAngleIncrement = wireAngle;
                        break;

                    case CoilAlignment::INNER_OR_TOP:
                        currentTurnCenterAngle = roundFloat(layerAngularCentre - layerAngularDimension / 2 + wireAngle / 2, 9);
                        currentTurnAngleIncrement = wireAngle;
                        break;

                    case CoilAlignment::OUTER_OR_BOTTOM:
                        currentTurnCenterAngle = roundFloat(layerAngularCentre + layerAngularDimension / 2 - totalLayerAngle + wireAngle / 2, 9);
                        currentTurnAngleIncrement = wireAngle;
                        break;

                    case CoilAlignment::SPREAD:
                        currentTurnAngleIncrement = roundFloat(layerAngularDimension / physicalTurnsInLayer, 9);
                        currentTurnCenterAngle = roundFloat(layerAngularCentre - layerAngularDimension / 2 + currentTurnAngleIncrement / 2, 9);
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
                    double spreadIncrement = roundFloat((layerAngularDimension - wireAngle) / physicalTurnsInLayer, 9);
                    if (spreadIncrement >= wireAngle) {
                        currentTurnAngleIncrement = spreadIncrement;
                        currentTurnCenterAngle = roundFloat(layerAngularCentre - physicalTurnsInLayer * spreadIncrement / 2 + spreadIncrement / 2, 9);
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

std::vector<std::pair<double, std::vector<double>>> Coil::get_collision_distances(const std::vector<double>& turnCoordinates, const std::vector<std::vector<double>>& placedTurnsCoordinates, double wireHeight) {
    std::vector<std::pair<double, std::vector<double>>> collisions;
    auto turnCartesianCoordinates = polar_to_cartesian(turnCoordinates);
    for (const auto& placedTurnCoordinates : placedTurnsCoordinates) {
        auto placedTurnCartesianCoordinates = polar_to_cartesian(placedTurnCoordinates);
        double distance = sqrt(pow(turnCartesianCoordinates[0] - placedTurnCartesianCoordinates[0], 2) + pow(turnCartesianCoordinates[1] - placedTurnCartesianCoordinates[1], 2));
        // Use a small tolerance to account for floating point precision
        if (distance < wireHeight - 1e-9) {
            double collisionDistance = wireHeight - distance;
            collisions.push_back({collisionDistance, placedTurnCoordinates});
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

            // Terminal-connection crossing lines of this section (real winding): each
            // winding-parallel's wire enters at its FIRST station's azimuth and leaves at its
            // LAST station's azimuth, and the connection runs RADIALLY across the faces there.
            // Outer-crossing candidates must keep a lateral wire clearance to these lines --
            // a rested crossing sitting on the connection line is copper through the lead.
            // The lead's own station is exempt (it IS that wire's continuation).
            // The connection VERTICALS in the final 3D: each winding-parallel's wire descends
            // BELOW the core at its FIRST station's exact position (entrance) and ascends ABOVE
            // at its LAST (exit). In the winding plane these are POINTS (the station positions);
            // the implied below-core return runs (outer crossing k -> inner station k+1) must
            // never come within a wire OD of the entrance point, and the top runs (inner k ->
            // outer k) must clear the exit point. Cartesian XY of a polar station:
            // r = windowRadialHeight - radialHeight at the station's angle.
            struct TerminalVertical {
                double x, y;              // station XY (the vertical's position)
                double wireOuterDiameter;
                bool below;               // true: entrance (below core); false: exit (above)
                std::string ownTurnName;
            };
            std::vector<TerminalVertical> terminalVerticals;
            if (settings.get_coil_use_real_winding_geometry()) {
                const double wwRadialHeight = windingWindows[0].get_radial_height().value();
                auto stationXY = [&](const Turn& t) {
                    double r = wwRadialHeight - t.get_coordinates()[0];
                    double aRad = t.get_coordinates()[1] / 180 * std::numbers::pi;
                    return std::make_pair(r * cos(aRad), r * sin(aRad));
                };
                std::map<std::pair<std::string, int64_t>, std::pair<const Turn*, const Turn*>> firstLastPerParallel;
                for (const auto& turn : turnsInSection) {
                    auto key = std::make_pair(turn.get_winding(), turn.get_parallel());
                    auto found = firstLastPerParallel.find(key);
                    if (found == firstLastPerParallel.end()) {
                        firstLastPerParallel[key] = {&turn, &turn};
                    }
                    else {
                        found->second.second = &turn;
                    }
                }
                for (const auto& [key, firstLast] : firstLastPerParallel) {
                    double odTerminal = wirePerWinding[get_winding_index_by_name(key.first)].get_maximum_outer_height();
                    auto [ex, ey] = stationXY(*firstLast.first);
                    terminalVerticals.push_back({ex, ey, odTerminal, true, firstLast.first->get_name()});
                    auto [xx, xy] = stationXY(*firstLast.second);
                    terminalVerticals.push_back({xx, xy, odTerminal, false, firstLast.second->get_name()});
                }
            }
            // Per-turn inner-station XY and the NEXT station's inner XY (same winding-parallel,
            // winding order): the below-core return implied by an outer-crossing candidate runs
            // from the candidate to the next inner station, and the top run from the own inner
            // station to the candidate. Both are swept against the connection verticals.
            std::map<std::string, std::pair<double, double>> ownInnerXYByTurn;
            std::map<std::string, std::pair<double, double>> nextInnerXYByTurn;
            if (settings.get_coil_use_real_winding_geometry()) {
                const double wwRadialHeight = windingWindows[0].get_radial_height().value();
                auto stationXY = [&](const Turn& t) {
                    double r = wwRadialHeight - t.get_coordinates()[0];
                    double aRad = t.get_coordinates()[1] / 180 * std::numbers::pi;
                    return std::make_pair(r * cos(aRad), r * sin(aRad));
                };
                std::map<std::pair<std::string, int64_t>, const Turn*> previousInParallel;
                for (const auto& turn : turnsInSection) {
                    auto key = std::make_pair(turn.get_winding(), turn.get_parallel());
                    ownInnerXYByTurn[turn.get_name()] = stationXY(turn);
                    auto found = previousInParallel.find(key);
                    if (found != previousInParallel.end()) {
                        nextInnerXYByTurn[found->second->get_name()] = stationXY(turn);
                    }
                    previousInParallel[key] = &turn;
                }
            }

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
                    // Winding progression sense of this layer (sign of the inner-azimuth step),
                    // for the outer-crossing monotonicity guard.
                    std::optional<double> previousOuterCrossingAzimuth;
                    double layerWindingDirection = 0.0;
                    if (turnsThisLayer.size() >= 2) {
                        double firstStep = std::remainder(turnsThisLayer[1].get_coordinates()[1] -
                                                          turnsThisLayer[0].get_coordinates()[1], 360.0);
                        layerWindingDirection = firstStep < 0 ? -1.0 : 1.0;
                    }
                    for (auto turn : turnsThisLayer) {
                        auto turnIndex = get_turn_index_by_name(turn.get_name());
                        std::vector<double> additionalCoordinates = {-bobbinColumnWidth * 2 - turn.get_coordinates()[0], turn.get_coordinates()[1]};

                        if (!areLayersTaped) {

                            if (!isFirstConductionLayer) {
                            // ABT #231. Two separate constraints, previously conflated by a single
                            // search that satisfied neither reliably.
                            //
                            // AZIMUTH is not free. A toroidal turn wraps the core at ONE angle, so its
                            // outer crossing must share the inner crossing's azimuth. The old search
                            // scanned angularly for a free slot, which produced outer angles out of
                            // sequence with the inner ones (measured 26.1, 41.1, 56.0, 48.5, 63.5
                            // against monotonic inner 29.9, 38.5, 47.1, 55.7, 64.3) and so crossed
                            // consecutive turns' top chords in 3D. There is no angular search now.
                            //
                            // RADIUS compacts, and only steps outward when it must. The outer face has
                            // a larger circumference than the bore, so a later ring's crossings usually
                            // interleave into the gaps left at the FIRST outer radius; forcing one wire
                            // OD of stacking per ring would be wrong, and is what
                            // Test_Additiona_Turns_Bug guards against. So: start at the innermost outer
                            // ring and step out by one wire only while the fixed-azimuth crossing still
                            // collides with an already-placed one.
                            //
                            // The old radial loop could exit while collisions remained and then accept
                            // the colliding placement, which is where the reported 0.87 OD centre-to-
                            // centre spacing came from. This one never accepts a collision: it either
                            // finds a clear radius or throws.
                            //
                            // COMPACTION fills the previous layers' gaps AS MUCH AS POSSIBLE. The wire's
                            // outer leg is not rigidly locked to the turn's poloidal plane: on a real
                            // part it LEANS azimuthally (up to about one wire OD of arc) into the
                            // nearest gap of the layers below and rests there. So the crossing takes
                            // the DEEPEST rest position within that lean window, where the rest radius
                            // at any azimuth is the packing-surface height: the base radius, pushed out
                            // by tangency against every already-placed crossing (centre distance == one
                            // wire OD: r = rP*cos(dAng) + sqrt(od^2 - rP^2*sin^2(dAng))). Ties resolve
                            // to the smallest lean. MONOTONICITY in winding order is enforced against
                            // the previous same-layer crossing -- the pre-ABT-#231 free angular search
                            // filled gaps too, but out of sequence, which crossed consecutive turns'
                            // top chords in 3D; the lean window plus the order guard keeps the fill
                            // without the crossings.
                            const double windowRadialHeight = windingWindows[0].get_radial_height().value();
                            const double baseRadius = windowRadialHeight - currentBaseRadialHeight;
                            // Lean + connection sweep are REAL-WINDING behaviour: the classic method
                            // keeps the plain fixed-azimuth tangency rest, so a problem in the new
                            // placement can be sidestepped by the flag alone.
                            const bool realWindingPlacement = settings.get_coil_use_real_winding_geometry();
                            const double leanDegrees = realWindingPlacement
                                ? (wireHeight / baseRadius) * 180 / std::numbers::pi : 0.0;
                            const double thetaDeg = additionalCoordinates[1];

                            // Candidate evaluation. A candidate azimuth yields: the tangency REST
                            // radius on the packing surface, and the FINAL-3D runs it implies --
                            // the top run (own inner station -> candidate) and the below-core
                            // return (candidate -> next inner station). A candidate is DISCARDED
                            // when either run comes within a wire clearance of a terminal
                            // connection's vertical (the exit vertical above the core for top
                            // runs, the entrance vertical below it for returns): the turns of the
                            // final 3D must never cross the connection wires. Inner coordinates
                            // are never touched; only this outer crossing moves.
                            auto restAt = [&](double az) {
                                double restRadius = baseRadius;
                                for (auto& placedCoordinates : placedTurnsCoordinates) {
                                    double placedRadius = windowRadialHeight - placedCoordinates[0];
                                    double dAng = std::remainder(placedCoordinates[1] - az, 360.0) / 180 * std::numbers::pi;
                                    double chord = placedRadius * sin(dAng);
                                    double discriminant = wireHeight * wireHeight - chord * chord;
                                    if (discriminant <= 0) {
                                        continue;   // cannot bind at this azimuth
                                    }
                                    double tangentRadius = placedRadius * cos(dAng) + sqrt(discriminant);
                                    restRadius = std::max(restRadius, tangentRadius);
                                }
                                return restRadius;
                            };
                            auto segmentPointDistance = [](double ax, double ay, double bx, double by,
                                                           double px, double py) {
                                double ux = bx - ax, uy = by - ay;
                                double len2 = ux * ux + uy * uy;
                                double t = len2 > 1e-18
                                    ? std::max(0.0, std::min(1.0, ((px - ax) * ux + (py - ay) * uy) / len2))
                                    : 0.0;
                                double cx = ax + ux * t, cy = ay + uy * t;
                                return std::hypot(px - cx, py - cy);
                            };
                            auto candidateAcceptable = [&](double az, double restRadius) {
                                if (previousOuterCrossingAzimuth.has_value() && layerWindingDirection != 0.0) {
                                    double progress = std::remainder(az - previousOuterCrossingAzimuth.value(), 360.0) * layerWindingDirection;
                                    if (progress <= 1e-9) {
                                        return false;   // reordered crossings: chords would cross in 3D
                                    }
                                }
                                double azRad = az / 180 * std::numbers::pi;
                                double cx = restRadius * cos(azRad), cy = restRadius * sin(azRad);
                                auto ownInner = ownInnerXYByTurn.find(turn.get_name());
                                auto nextInner = nextInnerXYByTurn.find(turn.get_name());
                                for (const auto& vertical : terminalVerticals) {
                                    double clearance = (wireHeight + vertical.wireOuterDiameter) / 2 - 1e-9;
                                    if (!vertical.below && vertical.ownTurnName != turn.get_name() &&
                                        ownInner != ownInnerXYByTurn.end()) {
                                        // top run vs an EXIT vertical (the exit's own top run IS
                                        // the connection wire)
                                        if (segmentPointDistance(ownInner->second.first, ownInner->second.second,
                                                                 cx, cy, vertical.x, vertical.y) < clearance) {
                                            return false;
                                        }
                                    }
                                    if (vertical.below && nextInner != nextInnerXYByTurn.end()) {
                                        // below-core return vs an ENTRANCE vertical
                                        if (segmentPointDistance(cx, cy, nextInner->second.first, nextInner->second.second,
                                                                 vertical.x, vertical.y) < clearance) {
                                            return false;
                                        }
                                    }
                                }
                                return true;
                            };

                            double bestRadius = std::numeric_limits<double>::max();
                            double bestAzimuth = thetaDeg;
                            if (!realWindingPlacement) {
                                bestRadius = restAt(thetaDeg);
                            }
                            else {
                                // Phase 1: the physical lean window (deepest nest wins, ties to the
                                // smallest lean).
                                const int leanSteps = 81;
                                for (int leanIndex = 0; leanIndex < leanSteps; ++leanIndex) {
                                    double az = thetaDeg + leanDegrees * (2.0 * leanIndex / (leanSteps - 1) - 1.0);
                                    double restRadius = restAt(az);
                                    if (!candidateAcceptable(az, restRadius)) {
                                        continue;
                                    }
                                    bool deeper = restRadius < bestRadius - 1e-9;
                                    bool tieCloser = std::abs(restRadius - bestRadius) <= 1e-9 &&
                                                     std::abs(az - thetaDeg) < std::abs(bestAzimuth - thetaDeg);
                                    if (deeper || tieCloser) {
                                        bestRadius = restRadius;
                                        bestAzimuth = az;
                                    }
                                }
                                // Phase 2: nothing in the lean window clears the connections --
                                // SWEEP THE WHOLE monotonic-feasible range WITHIN THE SECTION'S
                                // ANGULAR TERRITORY for the NEAREST azimuth whose implied runs are
                                // collision-free (radius is secondary: any rest beats a crossing).
                                // Leaving the section is never allowed: that azimuth belongs to a
                                // neighbouring winding's sector.
                                if (bestRadius == std::numeric_limits<double>::max()) {
                                    const double sectionCentre = section.get_coordinates()[1];
                                    const double sectionHalfSpan = section.get_dimensions()[1] / 2;
                                    const double sweepStep = 0.25;
                                    for (double offset = sweepStep; offset <= 2 * sectionHalfSpan; offset += sweepStep) {
                                        for (double sign : {1.0, -1.0}) {
                                            double az = thetaDeg + sign * offset;
                                            if (std::abs(std::remainder(az - sectionCentre, 360.0)) > sectionHalfSpan) {
                                                continue;
                                            }
                                            double restRadius = restAt(az);
                                            if (!candidateAcceptable(az, restRadius)) {
                                                continue;
                                            }
                                            bestRadius = restRadius;
                                            bestAzimuth = az;
                                            break;
                                        }
                                        if (bestRadius != std::numeric_limits<double>::max()) {
                                            break;
                                        }
                                    }
                                }
                            }
                            if (bestRadius == std::numeric_limits<double>::max() &&
                                realWindingPlacement && !_applyConnectionBlocking) {
                                // ABT #723: this sweep also runs at the end of the IDEAL
                                // (pre-blocking) wind, where the ABT #187 corridor
                                // reservations that would clear these runs have not been
                                // derived yet — they are derived FROM this wind by the
                                // blocking fixpoint that follows. Throwing here killed the
                                // pipeline before its own repair step could run (the
                                // section-contiguous 2-ring toroid: ring 1 fills back into
                                // the connection corridor on the ideal pass). Take the
                                // classic fixed-azimuth tangency rest for this intermediate
                                // pass; the corridor-blocked re-wind redoes the sweep with
                                // the stations moved off the corridor, and a failure THERE
                                // still throws below.
                                bestRadius = restAt(thetaDeg);
                                bestAzimuth = thetaDeg;
                            }
                            if (bestRadius == std::numeric_limits<double>::max()) {
                                std::string stationMap = "; section stations (layer:angle):";
                                for (const auto& sectionTurn : turnsInSection) {
                                    stationMap += " " + (sectionTurn.get_layer() ? sectionTurn.get_layer().value().substr(sectionTurn.get_layer().value().rfind(' ') + 1) : "?") +
                                                  ":" + std::to_string(sectionTurn.get_coordinates()[1]).substr(0, 6);
                                }
                                std::string blockers = "; verticals:";
                                for (const auto& vertical : terminalVerticals) {
                                    blockers += std::string(" ") + (vertical.below ? "in@" : "out@") +
                                                std::to_string(atan2(vertical.y, vertical.x) * 180 / std::numbers::pi).substr(0, 6);
                                }
                                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                                    "wind_toroidal_additional_turns: NO outer-crossing azimuth exists for turn " + turn.get_name() +
                                    " whose implied 3D runs clear the terminal connection verticals while keeping "
                                    "the winding order monotonic -- the inner-station layout leaves the connection "
                                    "corridor blocked (fix the layer spread / turn distribution, not this sweep)" +
                                    stationMap + blockers);
                            }
                            additionalCoordinates[0] = windowRadialHeight - bestRadius;
                            additionalCoordinates[1] = bestAzimuth;
                            if (!get_collision_distances(additionalCoordinates, placedTurnsCoordinates, wireHeight).empty()) {
                                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                                    "wind_toroidal_additional_turns: rested outer crossing still collides for turn " + turn.get_name());
                            }
                            }
                        }
                        currentSectionMaximumAdditionalRadialHeight = std::min(currentSectionMaximumAdditionalRadialHeight, additionalCoordinates[0]);
                        previousOuterCrossingAzimuth = additionalCoordinates[1];
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
    // Multi-window coils align against the group's allocated sub-region (which may be
    // half the window when the region is shared with the main winding's annulus).
    if (allWindingWindows.size() > 1 && sectionGroupOpt && get_groups_description()) {
        auto alignmentGroups = get_groups_description().value();
        for (auto& group : alignmentGroups) {
            if (group.get_name() == sectionGroupOpt.value()) {
                windingWindows[0].set_coordinates(std::vector<double>{std::abs(group.get_coordinates()[0]), group.get_coordinates()[1], 0});
                windingWindows[0].set_width(group.get_dimensions()[0]);
                windingWindows[0].set_height(group.get_dimensions()[1]);
                break;
            }
        }
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
    // ABT #720: the turn-axis coordinate is the SAME inner turnsAlignment switch in every
    // sectionAlignment case (eight hand-mirrored copies with the orientation split), and only
    // ONE copy carried the negative-coordinate guards. One copy each now, guards for all.
    auto turnAxisCoordinateOverlapping = [&]() -> double {
        switch (turnsAlignment) {
            case CoilAlignment::INNER_OR_TOP:
                return windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
            case CoilAlignment::OUTER_OR_BOTTOM:
                return windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
            case CoilAlignment::CENTERED: {
                double currentCoilHeightTop = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - resolve_margin(sections[sectionIndex])[0] - sections[sectionIndex].get_dimensions()[1] / 2;
                double currentCoilHeightBottom = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + resolve_margin(sections[sectionIndex])[1] + sections[sectionIndex].get_dimensions()[1] / 2;
                double currentCoilHeight = std::min(0.0, currentCoilHeightTop);
                return std::max(currentCoilHeight, currentCoilHeightBottom);
            }
            case CoilAlignment::SPREAD:
                return -resolve_margin(sections[sectionIndex])[0] / 2 + resolve_margin(sections[sectionIndex])[1] / 2;
            default:
                throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
        }
    };
    auto turnAxisCoordinateContiguous = [&]() -> double {
        switch (turnsAlignment) {
            case CoilAlignment::INNER_OR_TOP:
                return windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
            case CoilAlignment::OUTER_OR_BOTTOM:
                return windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
            case CoilAlignment::CENTERED: {
                double currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - sections[sectionIndex].get_dimensions()[0] / 2;
                double currentCoilWidthLeft = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
                double currentCoilWidthRight = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - resolve_margin(sections[sectionIndex])[1] - sections[sectionIndex].get_dimensions()[0];
                if (currentCoilWidthLeft < 0) {
                    throw std::invalid_argument("currentCoilWidthLeft cannot be less than 0: " + std::to_string(currentCoilWidthLeft));
                }
                if (currentCoilWidthRight < 0) {
                    throw std::invalid_argument("currentCoilWidthRight cannot be less than 0: " + std::to_string(currentCoilWidthRight));
                }
                currentCoilWidth = std::max(currentCoilWidth, currentCoilWidthLeft);
                return std::min(currentCoilWidth, currentCoilWidthRight);
            }
            case CoilAlignment::SPREAD:
                return windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + resolve_margin(sections[sectionIndex])[0];
            default:
                throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "No such section alignment");
        }
    };
    switch (sectionAlignment) {
        case CoilAlignment::INNER_OR_TOP:

            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                currentCoilHeight = turnAxisCoordinateOverlapping();
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2;
                currentCoilWidth = turnAxisCoordinateContiguous();
            }
            break;
        case CoilAlignment::OUTER_OR_BOTTOM:
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] + windingWindowWidth / 2 - totalSectionsWidth;
                currentCoilHeight = turnAxisCoordinateOverlapping();
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] - windingWindowHeight / 2 + totalSectionsHeight;
                currentCoilWidth = turnAxisCoordinateContiguous();
            }
            break;
        case CoilAlignment::SPREAD:
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                currentCoilHeight = turnAxisCoordinateOverlapping();
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

                currentCoilWidth = turnAxisCoordinateContiguous();
            }
            break;
        case CoilAlignment::CENTERED:
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                currentCoilWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2;
                currentCoilHeight = turnAxisCoordinateOverlapping();
            }
            else {
                currentCoilHeight = windingWindows[0].get_coordinates().value()[1] + totalSectionsHeight / 2;
                currentCoilWidth = turnAxisCoordinateContiguous();
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

    // The compaction cursor math lives in the +x winding frame. Callers may
    // re-compact an already-placed multi-window coil (e.g. CoilAdviser after
    // wind()): unwrap the placement transform, compact, and re-apply it.
    bool rewrapGroupWindowSides = _groupWindowSidesApplied;
    if (rewrapGroupWindowSides) {
        apply_group_window_sides(true);
    }

    bool result;
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        result = delimit_and_compact_rectangular_window();
    }
    else {
        result = delimit_and_compact_round_window();
    }

    if (rewrapGroupWindowSides) {
        apply_group_window_sides(false);
    }
    return result;
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
                        // ABT #616: the dims above are the exact turn envelope, so the stored
                        // filling factor (copper / PARTITION-time extent) is stale the moment the
                        // envelope differs — a 3-turn layer partitioned at 1.08 mm but wound to a
                        // 1.29 mm envelope read ff=1.19 and failed the fit gates while its copper
                        // exactly filled its rect. Refresh the ratio along the layer's turn axis
                        // from the same envelope just measured.
                        if (layers[i].get_type() == ElectricalType::CONDUCTION && !turnsInLayer.empty()) {
                            size_t ffAxis = (layers[i].get_orientation() == WindingOrientation::OVERLAPPING) ? 1 : 0;
                            double copperAlongAxis = 0;
                            for (auto& turn : turnsInLayer) {
                                copperAlongAxis += turn.get_dimensions().value()[ffAxis];
                            }
                            double envelope = (ffAxis == 1)
                                ? currentLayerMaximumHeight - currentLayerMinimumHeight
                                : currentLayerMaximumWidth - currentLayerMinimumWidth;
                            if (envelope > 0) {
                                layers[i].set_filling_factor(roundFloat(copperAlongAxis / envelope, 6));
                            }
                        }
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

                            // ABT #186: a toroidal turn's rotation is the cross-section azimuth and is
                            // set equal to its polar angle at creation. The angular compaction above shifts
                            // the polar angle, so rotation must be re-synced to it — otherwise MagneticField
                            // rotates the induced-field images by a stale angle and painters mis-orient the box.
                            turns[turnIndex].set_rotation(turns[turnIndex].get_coordinates()[1]);


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
            if (windingOrientation != WindingOrientation::OVERLAPPING) {
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
        for (const auto& turn : turns) {
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

Bobbin Coil::merge_per_column_bobbins(const std::vector<BobbinDataOrNameUnion> & perColumnBobbins) {
    auto resolveElement = [](const BobbinDataOrNameUnion & element, size_t columnIndex) -> Bobbin {
        if (std::holds_alternative<std::string>(element)) {
            auto name = std::get<std::string>(element);
            if (name == "Dummy") {
                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "coil.bobbin[" + std::to_string(columnIndex) + "] is Dummy: every element of the per-column bobbin array must be a real bobbin");
            }
            return find_bobbin_by_name(name);
        }
        return Bobbin(std::get<Bobbin>(element));
    };

    Bobbin mergedBobbin = resolveElement(perColumnBobbins[0], 0);
    if (!mergedBobbin.get_processed_description()) {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "coil.bobbin[0] (centre column) has no processedDescription: cannot merge per-column bobbins");
    }
    auto mergedProcessedDescription = mergedBobbin.get_processed_description().value();
    auto mergedWindingWindows = mergedProcessedDescription.get_winding_windows();
    for (size_t columnIndex = 1; columnIndex < perColumnBobbins.size(); ++columnIndex) {
        auto columnBobbin = resolveElement(perColumnBobbins[columnIndex], columnIndex);
        if (!columnBobbin.get_processed_description()) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "coil.bobbin[" + std::to_string(columnIndex) + "] has no processedDescription: cannot merge per-column bobbins");
        }
        auto columnWindingWindows = columnBobbin.get_processed_description().value().get_winding_windows();
        if (columnWindingWindows.empty()) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "coil.bobbin[" + std::to_string(columnIndex) + "] has no winding windows");
        }
        mergedWindingWindows.insert(mergedWindingWindows.end(), columnWindingWindows.begin(), columnWindingWindows.end());
    }
    mergedProcessedDescription.set_winding_windows(mergedWindingWindows);
    mergedBobbin.set_processed_description(mergedProcessedDescription);
    return mergedBobbin;
}

Bobbin Coil::resolve_bobbin() {
    if (_bobbin_resolved) {
        return _bobbin;
    }

    auto bobbinDataOrNameUnion = get_bobbin();
    if (std::holds_alternative<std::string>(bobbinDataOrNameUnion)) {
        if (std::get<std::string>(bobbinDataOrNameUnion) == "Dummy")
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Bobbin is dummy");

        _bobbin = find_bobbin_by_name(std::get<std::string>(bobbinDataOrNameUnion));
    }
    else {
        _bobbin = Bobbin(std::get<Bobbin>(bobbinDataOrNameUnion));
    }
    // The cache flag was never set, so every one of the ~57 call sites (several inside
    // per-section loops) re-ran the name lookup and copied the Bobbin. set_bobbin
    // invalidates; get_mutable_bobbin has no callers, so no mutation path bypasses this.
    _bobbin_resolved = true;
    return _bobbin;
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

    // The reallocation attempt below calls wind_by_sections(), which unconditionally
    // clears sections/layers/turns descriptions before rebuilding them, and its own
    // final re-wind is gated on the new layout actually fitting (unless
    // coilWindEvenIfNotFit is set). If the reallocation does not fix the fit, that
    // gate never fires and turns_description is left unset -- silently discarding
    // the layout just placed above, even though it is a valid (if marginally tight)
    // fallback. Snapshot it so it can be restored when the reallocation is no better.
    auto fallbackSectionsDescription = get_sections_description();
    auto fallbackLayersDescription = get_layers_description();
    auto fallbackTurnsDescription = get_turns_description();

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
                for (const auto& layer : layers) {
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
                for (const auto& layer : layers) {
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
                for (const auto& layer : layers) {
                    layerRestrictiveDimension = std::max(layerRestrictiveDimension, layer.get_dimensions()[1]);
                    layerFillingFactor = std::max(layerFillingFactor, layer.get_filling_factor().value());
                }

                layersRestrictiveDimension = layerRestrictiveDimension;
                extraSpaceNeededThisSection += std::max(0.0, (layerFillingFactor - 1) * layerRestrictiveDimension);
                windingWindowRemainingRestrictiveDimension -= layerRestrictiveDimension;
            }
            if (section.get_layers_orientation() == WindingOrientation::CONTIGUOUS) {
                for (const auto& layer : layers) {
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
        if (fallbackTurnsDescription) {
            set_sections_description(fallbackSectionsDescription);
            set_layers_description(fallbackLayersDescription);
            set_turns_description(fallbackTurnsDescription);
        }
        return;
    }
    // set_turns_description(std::nullopt);
    if (windEvenIfNotFit || are_sections_and_layers_fitting()) {
        wind_by_turns();
        if (delimitAndCompact) {
            delimit_and_compact();
        }
    }
    else if (fallbackTurnsDescription) {
        // The reallocated layout is no better than the fallback -- keep the
        // fallback's turns rather than leaving turns_description unset (ABT #621).
        set_sections_description(fallbackSectionsDescription);
        set_layers_description(fallbackLayersDescription);
        set_turns_description(fallbackTurnsDescription);
    }
}

void Coil::preload_margins(std::vector<std::vector<double>> marginPairs) {
    // Explicit margins re-arm the ABT #676 recovery for later winds (ABT #724).
    _marginsExplicitlyCleared = false;
    // ABT #720: _marginsPerSection is keyed by CONDUCTION-section ordinal — one entry per
    // conduction section, in wound order across all groups. (The old flat interleaved keying
    // forced a duplicated entry here "for the insulation layer", which broke the moment two
    // conduction sections sat adjacent without insulation.)
    for (auto margins : marginPairs) {
        _marginsPerSection.push_back(margins);
    }
}

void Coil::add_margin_to_section_by_index(size_t sectionIndex, std::vector<double> margins) {
    if (!get_sections_description()) {
        throw CoilNotProcessedException("In Add Margin to Section: Section description empty, wind coil first");
    }
    // Explicit margins re-arm the ABT #676 recovery for later winds (ABT #724).
    _marginsExplicitlyCleared = false;
    if (margins.size() != 2) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Margin vector must have two elements");
    }
    auto sections = get_sections_description().value();
    auto globalSectionIndex = convert_conduction_section_index_to_global(sectionIndex);
    // ABT #720: _marginsPerSection is keyed by CONDUCTION-section ordinal, which is exactly the
    // sectionIndex this method takes. It is only sized through the winders; when a Coil is
    // reconstructed from JSON (e.g. via the Coil(json, bool) constructor used by the
    // PyOpenMagnetics bindings and any caller that round-trips sectionsDescription through
    // JSON), it is empty even though sectionsDescription is populated, which made the indexed
    // assignment below segfault. Grow to the conduction-section count and seed any
    // uninitialized entries from the persisted sections' own margins, falling back to {0, 0}.
    std::vector<size_t> conductionGlobalIndexes;
    for (size_t i = 0; i < sections.size(); ++i) {
        if (sections[i].get_type() == ElectricalType::CONDUCTION) {
            conductionGlobalIndexes.push_back(i);
        }
    }
    if (_marginsPerSection.size() < conductionGlobalIndexes.size()) {
        size_t previousSize = _marginsPerSection.size();
        _marginsPerSection.resize(conductionGlobalIndexes.size(), {0, 0});
        for (size_t ordinal = previousSize; ordinal < conductionGlobalIndexes.size(); ++ordinal) {
            auto existingMargin = sections[conductionGlobalIndexes[ordinal]].get_margin();
            if (existingMargin) {
                if (std::holds_alternative<std::vector<double>>(existingMargin.value())) {
                    _marginsPerSection[ordinal] = std::get<std::vector<double>>(existingMargin.value());
                }
            }
        }
    }
    _marginsPerSection[sectionIndex] = margins;
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
    for (const auto& section : sections) {
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
    for (const auto& layer : layers) {
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
    for (const auto& section : sections) {
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
    for (const auto& layer : layers) {
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

    for (const auto& layer : layers) {
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

    // ABT #609: under REAL WINDING geometry, also offer each pattern's REVERSAL, appended AFTER
    // the base set. The n!/2 bound above deliberately drops reversed permutations because an
    // ideal winding is radially symmetric under reversal — but a real winding is not: the leads,
    // margins and blocking load the two orders differently, and 13_current_sense is the proof
    // (order 01 over-subscribes the ER 9.5 window to ~2x while order 10 fits and builds
    // watertight 3D). Appending — never interleaving — keeps the adviser's behaviour and
    // runtime identical whenever the base patterns deliver: the enumeration loop early-terminates
    // on a full candidate pool, so the reversals are only ever REACHED when the base set failed
    // to produce enough fitting candidates (Alf, 2026-08-08: "try more pattern combinations, but
    // make sure that time is not increased"). Gated behind the real-winding setting so every
    // existing flow with the setting off sees the exact historical pattern list.
    if (settings.get_coil_use_real_winding_geometry()) {
        size_t baseCount = sectionPatterns.size();
        for (size_t i = 0; i < baseCount && sectionPatterns.size() < defaults.maximumCoilPattern; ++i) {
            std::vector<size_t> reversed(sectionPatterns[i].rbegin(), sectionPatterns[i].rend());
            bool isNew = true;
            for (const auto& existing : sectionPatterns) {
                if (existing == reversed) {
                    isNew = false;
                    break;
                }
            }
            if (isNew) {
                sectionPatterns.push_back(reversed);
            }
        }
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

std::vector<double> Coil::resolve_margin(const Section& section) {
    if (!section.get_margin()) {
        return {0.0, 0.0};
    }
    return resolve_margin(section.get_margin().value());
}

std::vector<double> Coil::resolve_margin(const Margin& marginVariant) {
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

MarginInfo Coil::resolve_margin_info(const Section& section) {
    if (!section.get_margin()) {
        MarginInfo marginInfo;
        marginInfo.set_top_or_left_width(0);
        marginInfo.set_bottom_or_right_width(0);
        marginInfo.set_number_layers(0);
        // (The set_margin back onto the section that used to sit here only ever
        // mutated a discarded by-value copy.)
        return marginInfo;
    }
    return resolve_margin_info(section.get_margin().value());
}

MarginInfo Coil::resolve_margin_info(const Margin& marginVariant) {
    if (std::holds_alternative<std::vector<double>>(marginVariant)) {
        MarginInfo marginInfo;
        auto margin = std::get<std::vector<double>>(marginVariant);
        marginInfo.set_top_or_left_width(margin[0]);
        marginInfo.set_bottom_or_right_width(margin[1]);
        return marginInfo;
    }
    else {
        return std::get<MarginInfo>(marginVariant);
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
 