#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/Resistivity.h"
#include "Constants.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"
#include "support/Settings.h"

namespace OpenMagnetics {

std::vector<std::vector<double>> WindingOhmicLosses::calculate_connection_length_per_winding_per_parallel(Coil coil) {
    auto windings = coil.get_functional_description();
    std::vector<std::vector<double>> connectionLength;
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        connectionLength.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0.0));
    }

    // Connections only contribute when real winding geometry is enabled (otherwise ideal mode keeps
    // historical results untouched).
    if (!Settings::GetInstance().get_coil_use_real_winding_geometry()) {
        return connectionLength;
    }

    auto wirePerWinding = coil.get_wires();

    // ABT #492 owner ruling: planar wires are PCBs — the real-winding connection model is for WOUND
    // magnetics only. wind() already throws on this combination, but a coil deserialized from JSON
    // reaches this entry without ever winding, so the gate must live here too.
    for (const auto& wire : wirePerWinding) {
        if (wire.get_type() == WireType::PLANAR) {
            throw std::runtime_error(
                "Real winding geometry (connection/lead routing) is not implemented for planar "
                "(PCB) constructions; disable coilUseRealWindingGeometry for planar magnetics");
        }
    }

    // Provided terminal-lead lengths from the design requirements apply to the winding as a whole;
    // split evenly across its parallels (each parallel reaches the same terminal).
    std::vector<double> providedTerminalLength(windings.size(), 0.0);
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        if (windings[windingIndex].get_connections()) {
            auto connections = windings[windingIndex].get_connections().value();
            for (const auto& connection : connections) {
                if (connection.get_length()) {
                    providedTerminalLength[windingIndex] += connection.get_length().value();
                }
            }
        }
    }

    // Geometric lead lengths from the reserved-space model, PER PARALLEL: inter-layer crossings
    // (radial climb) and terminal leads (routed out to the window border). Each parallel is its own
    // conductor, so the lengths are accumulated per (winding, parallel).
    std::vector<std::vector<double>> crossingLength;
    std::vector<std::vector<double>> geometricTerminalLength;
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        crossingLength.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0.0));
        geometricTerminalLength.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0.0));
    }
    // ABT #492 loss-reader audit: every marker carries its centerline copper length EXPLICITLY in
    // routedLength, set by its emitter from its own geometry — nothing is inferred from the
    // rectangle here any more. The old inference was structurally wrong in two ways: an
    // orientation-derived index misread every segment running along the OTHER axis (a U-route
    // vertical stub counted as one wire width instead of its climb), and any shape-based rule
    // (long-axis, per-plane index) misreads a radial climb of a tall RECTANGULAR/FOIL wire whose
    // height exceeds the climb's run. Space-only squeeze markers carry 0 — their copper is the
    // DRAWN segments of the same route. An unset routedLength is an emitter bug (markers are
    // always freshly emitted at runtime) and throws per the no-fallback rule.
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        auto windingIndex = coil.get_winding_index_by_name(space.winding);
        int64_t parallelIndex = space.parallel;
        if (parallelIndex < 0 || parallelIndex >= int64_t(coil.get_number_parallels(windingIndex))) {
            continue;  // a winding-level lead with no parallel; cannot attribute to a branch
        }
        if (!space.routedLength) {
            throw std::runtime_error(
                "Connection marker without routedLength (winding '" + space.winding + "', section '" +
                space.section + "'): every emitter in Coil::get_connection_reserved_spaces must set "
                "the centerline copper length it charges — this is an emitter bug");
        }
        double runLength = space.routedLength.value();
        if (space.isTerminal) {
            geometricTerminalLength[windingIndex][parallelIndex] += runLength;
        }
        else {
            crossingLength[windingIndex][parallelIndex] += runLength;
        }
    }

    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        int64_t numberParallels = int64_t(coil.get_number_parallels(windingIndex));
        for (int64_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            // The terminal-lead length is taken from the design requirements when provided (shared
            // across the parallels), otherwise from the per-parallel geometric routing to the window
            // border (never both, to avoid double counting).
            double terminalLength = providedTerminalLength[windingIndex] > 0
                ? providedTerminalLength[windingIndex] / double(numberParallels)
                : geometricTerminalLength[windingIndex][parallelIndex];
            connectionLength[windingIndex][parallelIndex] = crossingLength[windingIndex][parallelIndex] + terminalLength;
        }
    }

    return connectionLength;
}

std::vector<std::vector<double>> WindingOhmicLosses::calculate_connection_resistance_per_winding_per_parallel(Coil coil, double temperature) {
    // DC resistance of the connection copper. The AC (skin) increment of the same lengths is added
    // by WindingSkinEffectLosses::calculate_skin_effect_losses from the SAME length source, so the
    // DC and AC stages can never disagree about what copper exists.
    auto connectionLength = calculate_connection_length_per_winding_per_parallel(coil);
    auto wirePerWinding = coil.get_wires();
    std::vector<std::vector<double>> connectionResistance;
    for (size_t windingIndex = 0; windingIndex < connectionLength.size(); ++windingIndex) {
        connectionResistance.push_back(std::vector<double>(connectionLength[windingIndex].size(), 0.0));
        bool anyLength = false;
        for (double length : connectionLength[windingIndex]) {
            if (length > 0) {
                anyLength = true;
            }
        }
        if (!anyLength) {
            continue;
        }
        double resistancePerMeter = calculate_dc_resistance_per_meter(wirePerWinding[windingIndex], temperature);
        for (size_t parallelIndex = 0; parallelIndex < connectionLength[windingIndex].size(); ++parallelIndex) {
            connectionResistance[windingIndex][parallelIndex] =
                resistancePerMeter * connectionLength[windingIndex][parallelIndex];
        }
    }

    return connectionResistance;
}

double WindingOhmicLosses::calculate_dc_resistance(Turn turn, const Wire& wire, double temperature) {
    double wireLength = turn.get_length();
    return calculate_dc_resistance(wireLength, wire, temperature);
}

double WindingOhmicLosses::calculate_dc_resistance(double wireLength, const Wire& wire, double temperature) {
    if (std::isnan(wireLength)) {
        throw NaNResultException("NaN found in wireLength value");
    }

    return calculate_dc_resistance_per_meter(wire, temperature) * wireLength;
}

double WindingOhmicLosses::calculate_dc_resistance_per_meter(Wire wire, double temperature) {

    WireMaterial wireMaterial = wire.resolve_material();

    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wireMaterial, temperature);

    double wireConductingArea = wire.calculate_conducting_area();

    double dcResistancePerMeter = resistivity / wireConductingArea;

    if (std::isnan(dcResistancePerMeter)) {
        throw NaNResultException("NaN found in dcResistancePerMeter value");
    }
    if (dcResistancePerMeter <= 0) {
        throw InvalidInputException(ErrorCode::CALCULATION_INVALID_RESULT, "dcResistancePerMeter must be positive");
    }
    return dcResistancePerMeter;
};

double WindingOhmicLosses::calculate_effective_resistance_per_meter(Wire wire, double frequency, double temperature) {
    WireMaterial wireMaterial = wire.resolve_material();

    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wireMaterial, temperature);

    double wireEffectiveConductingArea = wire.calculate_effective_conducting_area(frequency, temperature);

    // Check for invalid conducting area that would cause division issues
    if (std::isnan(wireEffectiveConductingArea) || std::isinf(wireEffectiveConductingArea) || wireEffectiveConductingArea <= 0) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, 
            "Invalid wire effective conducting area: " + std::to_string(wireEffectiveConductingArea));
    }

    double dcResistancePerMeter = resistivity / wireEffectiveConductingArea;
    return dcResistancePerMeter;
};

namespace {

// ABT #246: a parallel branch whose series resistance is zero (or non-finite) has no
// turns of that parallel in turnsDescription — an inconsistent coil, e.g. a stale or
// corrupted wind. The parallel-combination arithmetic then runs 1/0 -> infinite
// conductance -> zero parallel resistance -> 0/0 branch current = NaN, and that NaN
// propagates into the winding-losses TOTAL, which nlohmann serialises as JSON null.
// The emitted MAS then fails MAS's own schema (outputs.windingLosses.windingLosses is
// a required number). Fail loudly, naming the branch, instead of emitting the NaN.
std::map<std::pair<size_t, size_t>, size_t> count_turns_per_winding_per_parallel(Coil& coil, const std::vector<Turn>& turns) {
    std::map<std::pair<size_t, size_t>, size_t> counts;
    for (auto& turn : turns) {
        counts[{coil.get_winding_index_by_name(turn.get_winding()), turn.get_parallel()}]++;
    }
    return counts;
}

void check_parallel_branch(Coil& coil,
                           size_t windingIndex,
                           size_t parallelIndex,
                           double seriesResistance,
                           const std::map<std::pair<size_t, size_t>, size_t>& turnCounts) {
    if (std::isfinite(seriesResistance) && seriesResistance > 0) {
        return;
    }
    auto found = turnCounts.find({windingIndex, parallelIndex});
    size_t turnsThisParallel = found == turnCounts.end() ? 0 : found->second;
    throw CoilException(ErrorCode::COIL_WINDING_ERROR,
        "Winding '" + coil.get_functional_description()[windingIndex].get_name() + "' parallel " +
        std::to_string(parallelIndex) + " has a series resistance of " + std::to_string(seriesResistance) +
        " Ohm: " + std::to_string(turnsThisParallel) + " of its turns are in turnsDescription, while " +
        "functionalDescription declares " + std::to_string(coil.get_number_turns(windingIndex)) +
        " turns over " + std::to_string(coil.get_number_parallels(windingIndex)) +
        " parallels. The coil must be re-wound before its losses can be computed.");
}

}  // namespace

std::vector<double> WindingOhmicLosses::calculate_dc_resistance_per_winding(Coil coil, double temperature) {
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description");
    }
    auto turns = coil.get_turns_description().value();
    std::vector<std::vector<double>> seriesResistancePerWindingPerParallel;
    auto wirePerWinding = coil.get_wires();
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        seriesResistancePerWindingPerParallel.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0));
    }

    std::vector<double> dcResistancePerWinding;
    for (auto& turn : turns) {
        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto parallelIndex = turn.get_parallel();

        double turnResistance = calculate_dc_resistance(turn, wirePerWinding[windingIndex], temperature);
        seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] += turnResistance;
    }

    // Terminal/connection leads add series resistance to each parallel branch (zero in ideal mode):
    // a parallel's leads are in series with its turns, then the branches combine in parallel.
    auto connectionResistancePerWindingPerParallel = calculate_connection_resistance_per_winding_per_parallel(coil, temperature);
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] += connectionResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
    }

    auto turnCounts = count_turns_per_winding_per_parallel(coil, turns);
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        double conductance = 0;
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            check_parallel_branch(coil, windingIndex, parallelIndex,
                                  seriesResistancePerWindingPerParallel[windingIndex][parallelIndex], turnCounts);
            conductance += 1. / seriesResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
        double parallelResistance = 1. / conductance;
        dcResistancePerWinding.push_back(parallelResistance);
    }

    return dcResistancePerWinding;
}

WindingLossesOutput WindingOhmicLosses::calculate_ohmic_losses(Coil coil, OperatingPoint operatingPoint, double temperature) {
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description");
    }
    auto turns = coil.get_turns_description().value();
    std::vector<std::vector<double>> seriesResistancePerWindingPerParallel;
    std::vector<std::vector<double>> dcCurrentPerWindingPerParallel;
    std::vector<double> dcCurrentPerWinding;
    auto wirePerWinding = coil.get_wires();
    
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        seriesResistancePerWindingPerParallel.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0));
        dcCurrentPerWindingPerParallel.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0));
        
        if (windingIndex >= operatingPoint.get_excitations_per_winding().size()) {
            dcCurrentPerWinding.push_back(0);
            continue;
        }
        
        auto& exc = operatingPoint.get_excitations_per_winding()[windingIndex];
        if (!exc.get_current() || !exc.get_current()->get_processed() || !exc.get_current()->get_processed()->get_rms()) {
            // Missing RMS is handled the same way as a missing current spec (the guard
            // used to stop at get_processed() and then .value() the RMS — a raw
            // bad_optional_access instead of the documented absent-data path).
            dcCurrentPerWinding.push_back(0);
            continue;
        }

        double currentRms = exc.get_current()->get_processed()->get_rms().value();
        dcCurrentPerWinding.push_back(currentRms);
    }

    std::vector<double> dcResistancePerTurn;
    std::vector<double> dcResistancePerWinding;
    for (auto& turn : turns) {
        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto parallelIndex = turn.get_parallel();

        double turnResistance = calculate_dc_resistance(turn, wirePerWinding[windingIndex], temperature);
        dcResistancePerTurn.push_back(turnResistance);
        seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] += turnResistance;
    }

    // Terminal/connection leads add series resistance to each parallel branch (zero in ideal mode):
    // each parallel's leads are in series with its own turns, then the branches combine in parallel.
    auto connectionResistancePerWindingPerParallel = calculate_connection_resistance_per_winding_per_parallel(coil, temperature);
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] += connectionResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
    }

    auto turnCounts = count_turns_per_winding_per_parallel(coil, turns);
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        double conductance = 0;
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            check_parallel_branch(coil, windingIndex, parallelIndex,
                                  seriesResistancePerWindingPerParallel[windingIndex][parallelIndex], turnCounts);
            conductance += 1. / seriesResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
        double parallelResistance = 1. / conductance;
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            dcCurrentPerWindingPerParallel[windingIndex][parallelIndex] = dcCurrentPerWinding[windingIndex] * parallelResistance / seriesResistancePerWindingPerParallel[windingIndex][parallelIndex];
        }
        dcResistancePerWinding.push_back(parallelResistance);
    }
    std::vector<WindingLossesPerElement> windingLossesPerTurn;
    std::vector<double> currentDividerPerTurn;
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        Turn turn = turns[turnIndex];
        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto parallelIndex = turn.get_parallel();

        auto currentDividerThisTurn = dcCurrentPerWinding[windingIndex] == 0? 0 : dcCurrentPerWindingPerParallel[windingIndex][parallelIndex] / dcCurrentPerWinding[windingIndex];

        double windingOhmicLossesInTurn = pow(dcCurrentPerWindingPerParallel[windingIndex][parallelIndex], 2) * dcResistancePerTurn[turnIndex];
        OhmicLosses ohmicLosses;
        WindingLossesPerElement windingLossesThisTurn;
        ohmicLosses.set_losses(windingOhmicLossesInTurn);
        ohmicLosses.set_method_used("Ohm");
        ohmicLosses.set_origin(ResultOrigin::SIMULATION);
        windingLossesThisTurn.set_ohmic_losses(ohmicLosses);
        windingLossesThisTurn.set_name(turn.get_name());
        windingLossesPerTurn.push_back(windingLossesThisTurn);

        if (std::isnan(currentDividerThisTurn)) {
            throw NaNResultException("NaN found in currentDividerThisTurn value");
        }
        currentDividerPerTurn.push_back(currentDividerThisTurn);
    }

    double windingOhmicLossesTotal = 0;

    std::vector<WindingLossesPerElement> windingLossesPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        double windingOhmicLossesInWinding = 0;

        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            // seriesResistance already includes this parallel's connection-lead resistance.
            windingOhmicLossesInWinding += seriesResistancePerWindingPerParallel[windingIndex][parallelIndex] * pow(dcCurrentPerWindingPerParallel[windingIndex][parallelIndex], 2);
        }

        OhmicLosses ohmicLosses;
        WindingLossesPerElement windingLossesThisWinding;
        ohmicLosses.set_losses(windingOhmicLossesInWinding);
        ohmicLosses.set_method_used("Ohm");
        ohmicLosses.set_origin(ResultOrigin::SIMULATION);
        windingLossesThisWinding.set_ohmic_losses(ohmicLosses);
        windingLossesThisWinding.set_name(coil.get_functional_description()[windingIndex].get_name());
        windingLossesPerWinding.push_back(windingLossesThisWinding);
        windingOhmicLossesTotal += windingOhmicLossesInWinding;
    }

    WindingLossesOutput result;
    result.set_winding_losses_per_winding(windingLossesPerWinding);
    result.set_winding_losses_per_turn(windingLossesPerTurn);
    result.set_winding_losses(windingOhmicLossesTotal);
    result.set_temperature(temperature);
    result.set_temperature(temperature);
    result.set_origin(ResultOrigin::SIMULATION);
    result.set_dc_resistance_per_turn(dcResistancePerTurn);
    result.set_dc_resistance_per_winding(dcResistancePerWinding);
    result.set_current_per_winding(operatingPoint);
    result.set_current_divider_per_turn(currentDividerPerTurn);

    return result;
}

double WindingOhmicLosses::calculate_ohmic_losses_per_meter(Wire wire, SignalDescriptor current, double temperature) {

    double dcResistancePerMeter = calculate_dc_resistance_per_meter(wire, temperature);
    if (!current.get_processed()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is not processed");
    }
    if (!current.get_processed()->get_rms()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current processed is missing field RMS");
    }
    auto currentRms = current.get_processed()->get_rms().value();

    double windingOhmicLossesPerMeter = pow(currentRms, 2) * dcResistancePerMeter;

    return windingOhmicLossesPerMeter;
};
} // namespace OpenMagnetics
