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

std::vector<std::vector<double>> WindingOhmicLosses::calculate_connection_resistance_per_winding_per_parallel(Coil coil, double temperature) {
    auto windings = coil.get_functional_description();
    std::vector<std::vector<double>> connectionResistance;
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        connectionResistance.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0.0));
    }

    // Connections only contribute when real winding geometry is enabled (otherwise ideal mode keeps
    // historical results untouched).
    if (!Settings::GetInstance().get_coil_use_real_winding_geometry()) {
        return connectionResistance;
    }

    auto wirePerWinding = coil.get_wires();

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
    for (const auto& space : coil.get_connection_reserved_spaces()) {
        auto windingIndex = coil.get_winding_index_by_name(space.winding);
        int64_t parallelIndex = space.parallel;
        if (parallelIndex < 0 || parallelIndex >= int64_t(coil.get_number_parallels(windingIndex))) {
            continue;  // a winding-level lead with no parallel; cannot attribute to a branch
        }
        if (space.isTerminal) {
            geometricTerminalLength[windingIndex][parallelIndex] += space.dimensions[0];
        }
        else {
            crossingLength[windingIndex][parallelIndex] += space.dimensions[0];
        }
    }

    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        int64_t numberParallels = int64_t(coil.get_number_parallels(windingIndex));
        double resistancePerMeter = calculate_dc_resistance_per_meter(wirePerWinding[windingIndex], temperature);
        for (int64_t parallelIndex = 0; parallelIndex < numberParallels; ++parallelIndex) {
            // The terminal-lead length is taken from the design requirements when provided (shared
            // across the parallels), otherwise from the per-parallel geometric routing to the window
            // border (never both, to avoid double counting).
            double terminalLength = providedTerminalLength[windingIndex] > 0
                ? providedTerminalLength[windingIndex] / double(numberParallels)
                : geometricTerminalLength[windingIndex][parallelIndex];
            double leadLength = crossingLength[windingIndex][parallelIndex] + terminalLength;
            if (leadLength > 0) {
                connectionResistance[windingIndex][parallelIndex] += resistancePerMeter * leadLength;
            }
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

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        double conductance = 0;
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
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
        if (!exc.get_current() || !exc.get_current()->get_processed()) {
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

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        double conductance = 0;
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
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
