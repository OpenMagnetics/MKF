#include "processors/CircuitSimulatorInterface.h"
#include "processors/CircuitSimulatorExporterHelpers.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingLosses.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "Defaults.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace OpenMagnetics {

// Forward declarations for free functions defined in CircuitSimulatorInterface.cpp.
std::string to_string(double d, size_t precision);
std::vector<std::string> to_string(std::vector<double> v, size_t precision);

CircuitSimulatorExporterSimbaModel::CircuitSimulatorExporterSimbaModel() {
    // FIXED: Was seeding with random_device then immediately overwriting with time(0)
    // Now use random_device only, which provides proper entropy
    std::random_device rd;
    _gen = std::mt19937(rd());
}
std::string CircuitSimulatorExporterSimbaModel::generate_id() {
    // generator for hex numbers from 0 to F
    std::uniform_int_distribution<int> dis(0, 15);
    
    const char hexChars[] = { '0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F' };
    const std::size_t hexCharsSize = sizeof(hexChars) / sizeof(hexChars[0]);
    std::string id = "";
    for (std::size_t i = 0; i < 8; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    id += "-";
    for (std::size_t i = 0; i < 4; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    id += "-";
    for (std::size_t i = 0; i < 4; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    id += "-";
    for (std::size_t i = 0; i < 4; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    id += "-";
    for (std::size_t i = 0; i < 12; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    return id;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_device(std::string libraryName, std::vector<int> coordinates, int angle, std::string name) {
    ordered_json deviceJson;

    deviceJson["LibraryName"] = libraryName;
    deviceJson["Top"] = coordinates[1];
    deviceJson["Left"] = coordinates[0];
    deviceJson["Angle"] = angle;
    deviceJson["HF"] = false;
    deviceJson["VF"] = false;
    deviceJson["Disabled"] = false;
    deviceJson["Name"] = name;
    deviceJson["ID"] = generate_id();
    deviceJson["Parameters"]["Name"] = name;
    deviceJson["EnabledScopes"] = ordered_json::array();

    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_air_gap(std::vector<int> coordinates, double area, double length, int angle, std::string name) {
    ordered_json deviceJson = create_device("Air Gap", coordinates, angle, name);

    deviceJson["Parameters"]["RelativePermeability"] = "1";
    deviceJson["Parameters"]["Area"] = to_string(area, _precision);
    deviceJson["Parameters"]["Length"] = to_string(length, _precision);
    return deviceJson;
}
ordered_json CircuitSimulatorExporterSimbaModel::create_core(double initialPermeability, std::vector<int> coordinates, double area, double length, int angle, std::string name) {
    ordered_json deviceJson = create_device("Linear Core", coordinates, angle, name);

    deviceJson["Parameters"]["RelativePermeability"] = to_string(initialPermeability, _precision);
    deviceJson["Parameters"]["Area"] = to_string(area, _precision);
    deviceJson["Parameters"]["Length"] = to_string(length, _precision);
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_winding(size_t numberTurns, std::vector<int> coordinates, int angle, std::string name) {
    ordered_json deviceJson = create_device("Winding", coordinates, angle, name);

    deviceJson["Parameters"]["NumberOfTurns"] = to_string(numberTurns, _precision);
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_resistor(double resistance, std::vector<int> coordinates, int angle, std::string name) {

    ordered_json deviceJson = create_device("Resistor", coordinates, angle, name);

    deviceJson["Parameters"]["Value"] = to_string(resistance, _precision);
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_inductor(double inductance, std::vector<int> coordinates, int angle, std::string name) {

    ordered_json deviceJson = create_device("Inductor", coordinates, angle, name);

    deviceJson["Parameters"]["Value"] = to_string(inductance, _precision);
    deviceJson["Parameters"]["Iinit"] = "0";
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_capacitor(double capacitance, std::vector<int> coordinates, int angle, std::string name) {
    ordered_json deviceJson = create_device("Capacitor", coordinates, angle, name);
    deviceJson["Parameters"]["Value"] = to_string(capacitance, _precision);
    deviceJson["Parameters"]["Vinit"] = "0";
    return deviceJson;
}

std::pair<std::vector<ordered_json>, std::vector<ordered_json>> CircuitSimulatorExporterSimbaModel::create_ladder(std::vector<double> ladderCoefficients, std::vector<int> coordinates, std::string name) {
    std::vector<ordered_json> ladderJsons;
    std::vector<ordered_json> ladderConnectorsJsons;
    coordinates[0] -= 6;
    for (size_t ladderIndex = 0; ladderIndex < ladderCoefficients.size(); ladderIndex+=2) {
        std::string ladderIndexs = std::to_string(ladderIndex);
        auto ladderInductorJson = create_inductor(ladderCoefficients[ladderIndex + 1], coordinates, 180, name + " Ladder element " + std::to_string(ladderIndex));
        ladderJsons.push_back(ladderInductorJson);
        coordinates[1] -= 3;
        coordinates[0] += 3;
        auto ladderResistorJson = create_resistor(ladderCoefficients[ladderIndex], coordinates, 90, name + " Ladder element " + std::to_string(ladderIndex));
        ladderJsons.push_back(ladderResistorJson);
        coordinates[1] -= 3;
        coordinates[0] -= 3;
        {
            std::vector<int> finalConnectorTopCoordinates = {coordinates[0], coordinates[1] + 1};
            std::vector<int> finalConnectorBottomCoordinates = {coordinates[0], coordinates[1] + 7};
            auto ladderConnectorJson = create_connector(finalConnectorTopCoordinates, finalConnectorBottomCoordinates, "Connector between winding " + name + " ladder elements" + std::to_string(ladderIndex) + " and " + std::to_string(ladderIndex + 2));
            ladderConnectorsJsons.push_back(ladderConnectorJson);
        }
        if (ladderIndex == ladderCoefficients.size() - 2) {
            std::vector<int> finalConnectorTopCoordinates = {coordinates[0], coordinates[1] + 1};
            std::vector<int> finalConnectorBottomCoordinates = {coordinates[0] + 6, coordinates[1] + 1};
            auto ladderConnectorJson = create_connector(finalConnectorTopCoordinates, finalConnectorBottomCoordinates, "Connector between winding " + name + " ladder elements" + std::to_string(ladderIndex) + " and " + std::to_string(ladderIndex + 2));
            ladderConnectorsJsons.push_back(ladderConnectorJson);
        }
    }

    return {ladderJsons, ladderConnectorsJsons};
}

std::pair<std::vector<ordered_json>, std::vector<ordered_json>> CircuitSimulatorExporterSimbaModel::create_rosano_ladder(
    std::vector<double> coefficients, std::vector<int> coordinates, std::string name) {
    // Rosano topology: Z = R1 + (L1||R2) + (L2||(R3+C1))
    // Three stages in series. Each stage is a two-terminal block.
    // Layout: each stage occupies a horizontal row.
    //   - Stage R: single resistor (horizontal)
    //   - Stage RL: L on left rail, R on right rail, both between same top/bottom nodes
    //   - Stage RLC: L on left rail, R+C on right rail (R then C vertical)
    // coefficients: [R1, R2, L1, R3, L2, C1]
    std::vector<ordered_json> deviceJsons;
    std::vector<ordered_json> connectorJsons;

    if (coefficients.size() < 6) return {deviceJsons, connectorJsons};

    double R1 = coefficients[0];
    double R2 = coefficients[1];
    double L1 = coefficients[2];
    double R3 = coefficients[3];
    double L2 = coefficients[4];
    double C1 = coefficients[5];

    // We use two vertical rails: left rail (x=coordinates[0]-6) and right rail (x=coordinates[0])
    // Components on the left rail are horizontal (angle 180), on the right rail vertical (angle 90)
    // Each stage: left component from nodeA to nodeB, right component from nodeA to nodeB (parallel)
    // Then a connector links nodeB of this stage to nodeA of the next stage

    int x = coordinates[0];
    int y = coordinates[1];

    // Stage 1: R1 series (single horizontal resistor)
    if (R1 > 0) {
        auto rJson = create_resistor(R1, {x - 6, y}, 180, name + " Stage1 R");
        deviceJsons.push_back(rJson);
        y -= 6;
        // Connector from R1 output to Stage 2 input
        connectorJsons.push_back(create_connector({x - 6, y + 1}, {x - 6, y + 7}, name + " S1-S2 connector"));
    }

    // Stage 2: L1 || R2 (two components in parallel between same two nodes)
    if (R2 > 0 && L1 > 0) {
        // L1: horizontal on left rail
        auto l1Json = create_inductor(L1, {x - 6, y}, 180, name + " Stage2 L");
        deviceJsons.push_back(l1Json);
        // R2: horizontal on right rail (offset by +6 in x), same y
        auto r2Json = create_resistor(R2, {x, y}, 180, name + " Stage2 R");
        deviceJsons.push_back(r2Json);
        // Connect top nodes (L1 top to R2 top)
        connectorJsons.push_back(create_connector({x - 6, y + 5}, {x, y + 5}, name + " Stage2 top"));
        // Connect bottom nodes (L1 bottom to R2 bottom)
        connectorJsons.push_back(create_connector({x - 6, y - 1}, {x, y - 1}, name + " Stage2 bot"));
        y -= 6;
        // Connector to Stage 3
        connectorJsons.push_back(create_connector({x - 6, y + 1}, {x - 6, y + 7}, name + " S2-S3 connector"));
    }

    // Stage 3: L2 || (R3 + C1)
    if (R3 > 0 && L2 > 0 && C1 > 0) {
        // L2: horizontal on left rail
        auto l2Json = create_inductor(L2, {x - 6, y}, 180, name + " Stage3 L");
        deviceJsons.push_back(l2Json);
        // R3: vertical on right rail (top of R+C series)
        auto r3Json = create_resistor(R3, {x + 3, y}, 90, name + " Stage3 R");
        deviceJsons.push_back(r3Json);
        // C1: vertical below R3
        auto c1Json = create_capacitor(C1, {x + 3, y - 3}, 90, name + " Stage3 C");
        deviceJsons.push_back(c1Json);
        // Connect L2 top to R3 top
        connectorJsons.push_back(create_connector({x - 6, y + 5}, {x + 3, y + 5}, name + " Stage3 top"));
        // Connect L2 bottom to C1 bottom
        connectorJsons.push_back(create_connector({x - 6, y - 1}, {x + 3, y - 1 - 3}, name + " Stage3 bot"));
        // Final connector to close the circuit
        connectorJsons.push_back(create_connector({x - 6, y - 1}, {x - 6 + 6, y - 1}, name + " Stage3 final"));
    }

    return {deviceJsons, connectorJsons};
}

std::pair<std::vector<ordered_json>, std::vector<ordered_json>> CircuitSimulatorExporterSimbaModel::create_fracpole_ladder(
    const FractionalPoleNetwork& net, std::vector<int> coordinates, const std::string& name) {
    std::vector<ordered_json> ladderJsons;
    std::vector<ordered_json> ladderConnectorsJsons;
    coordinates[0] -= 6;

    for (size_t k = 0; k < net.stages.size(); ++k) {
        std::string ks = std::to_string(k);
        // Resistor (series, horizontal)
        auto resistorJson = create_resistor(net.stages[k].R * net.opts.coef, coordinates, 180,
                                            name + " Fracpole R" + ks);
        ladderJsons.push_back(resistorJson);
        coordinates[1] -= 3;
        coordinates[0] += 3;
        // Capacitor (shunt, vertical) - replaces the inductor in ladder mode
        auto capacitorJson = create_capacitor(net.stages[k].C / net.opts.coef, coordinates, 90,
                                              name + " Fracpole C" + ks);
        ladderJsons.push_back(capacitorJson);
        coordinates[1] -= 3;
        coordinates[0] -= 3;
        // Connector between stages
        {
            std::vector<int> top = {coordinates[0], coordinates[1] + 1};
            std::vector<int> bot = {coordinates[0], coordinates[1] + 7};
            auto connJson = create_connector(top, bot,
                "Fracpole connector " + name + " stage " + ks);
            ladderConnectorsJsons.push_back(connJson);
        }
        if (k == net.stages.size() - 1) {
            std::vector<int> top = {coordinates[0], coordinates[1] + 1};
            std::vector<int> bot = {coordinates[0] + 6, coordinates[1] + 1};
            auto connJson = create_connector(top, bot, "Fracpole final connector " + name);
            ladderConnectorsJsons.push_back(connJson);
        }
    }
    return {ladderJsons, ladderConnectorsJsons};
}

ordered_json CircuitSimulatorExporterSimbaModel::create_pin(std::vector<int> coordinates, int angle, std::string name) {
    return create_device("Electrical Pin", coordinates, angle, name);
}

ordered_json CircuitSimulatorExporterSimbaModel::create_magnetic_ground(std::vector<int> coordinates, int angle, std::string name) {
    return create_device("Magnetic Ground", coordinates, angle, name);
}

ordered_json CircuitSimulatorExporterSimbaModel::create_connector(std::vector<int> startingCoordinates, std::vector<int> endingCoordinates, std::string name) {
    ordered_json connectorJson;

    connectorJson["Name"] = name;
    connectorJson["Segments"] = ordered_json::array();

    if (startingCoordinates[0] == endingCoordinates[0] || startingCoordinates[1] == endingCoordinates[1]) {
        ordered_json segment;
        segment["StartX"] = startingCoordinates[0];
        segment["StartY"] = startingCoordinates[1];
        segment["EndX"] = endingCoordinates[0];
        segment["EndY"] = endingCoordinates[1];
        connectorJson["Segments"].push_back(segment);
    }
    else {
        {
            ordered_json segment;
            segment["StartX"] = startingCoordinates[0];
            segment["StartY"] = startingCoordinates[1];
            segment["EndX"] = endingCoordinates[0];
            segment["EndY"] = startingCoordinates[1];
            connectorJson["Segments"].push_back(segment);
        }
        {
            ordered_json segment;
            segment["StartX"] = endingCoordinates[0];
            segment["StartY"] = startingCoordinates[1];
            segment["EndX"] = endingCoordinates[0];
            segment["EndY"] = endingCoordinates[1];
            connectorJson["Segments"].push_back(segment);
        }
    }

    return connectorJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::merge_connectors(ordered_json connectors) {
    auto mergeConnectors = connectors;
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t firstIndex = 0; firstIndex < mergeConnectors.size(); ++firstIndex) {
            bool sharedPoint = false;
            ordered_json mergedConnector;
            for (size_t secondIndex = firstIndex + 1; secondIndex < mergeConnectors.size(); ++secondIndex) {
                for (auto firstConnectorSegment : mergeConnectors[firstIndex]["Segments"]) {
                    for (auto secondConnectorSegment : mergeConnectors[secondIndex]["Segments"]) {
                        if ((firstConnectorSegment["StartX"] == secondConnectorSegment["StartX"] &&
                            firstConnectorSegment["StartY"] == secondConnectorSegment["StartY"]) ||
                            (firstConnectorSegment["EndX"] == secondConnectorSegment["EndX"] &&
                            firstConnectorSegment["EndY"] == secondConnectorSegment["EndY"])) {

                            sharedPoint = true;
                            break;
                        }
                    }
                    if (sharedPoint) {
                        break;
                    }
                }

                if (sharedPoint) {
                    mergedConnector["Segments"] = ordered_json::array();
                    mergedConnector["Name"] = "Merge of connector: " + std::string{mergeConnectors[firstIndex]["Name"]} + " with " + std::string{mergeConnectors[secondIndex]["Name"]};
                    for (auto firstConnectorSegment : mergeConnectors[firstIndex]["Segments"]) {
                        mergedConnector["Segments"].push_back(firstConnectorSegment);
                    }
                    for (auto secondConnectorSegment : mergeConnectors[secondIndex]["Segments"]) {
                        mergedConnector["Segments"].push_back(secondConnectorSegment);
                    }
                    mergeConnectors.erase(secondIndex);
                    mergeConnectors.erase(firstIndex);
                    mergeConnectors.push_back(mergedConnector);
                    merged = true;
                    break;
                }
            }

            if (sharedPoint) {
                break;
            }
        }
    }
    return mergeConnectors;
}



std::string CircuitSimulatorExporterSimbaModel::export_magnetic_as_subcircuit(Magnetic magnetic, double frequency, double temperature, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    ordered_json simulation;
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    _scale = _modelSize / core.get_width();


    if (filePathOrFile) {
        try {
            if(!std::filesystem::exists(filePathOrFile.value())) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "File not found");
            }
            std::ifstream json_file(filePathOrFile.value());
            if(json_file.is_open()) {
                simulation = ordered_json::parse(json_file);
            }
            else {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "File not found");
            }
        }
        catch(const InvalidInputException& re) {
            std::stringstream json_file(filePathOrFile.value());
            simulation = ordered_json::parse(json_file);
        }
    }

    if (!simulation.contains("Libraries")) {
        simulation["Libraries"] = ordered_json::array();
    }
    if (!simulation.contains("Designs")) {
        simulation["Designs"] = ordered_json::array();
    }
    ordered_json library;
    ordered_json device;
    library["LibraryItemName"] = "OpenMagnetics components";
    // library["Devices"] = ordered_json::array();

    device["LibraryName"] = magnetic.get_reference();
    device["Angle"] = 0;
    device["Disabled"] = false;
    device["Name"] = magnetic.get_reference();
    device["Id"] = generate_id();
    device["Parameters"]["Name"] = magnetic.get_reference();
    device["SubcircuitDefinition"]["Devices"] = ordered_json::array();
    device["SubcircuitDefinition"]["Connectors"] = ordered_json::array();
    device["SubcircuitDefinition"]["Name"] = magnetic.get_reference();
    device["SubcircuitDefinition"]["Id"] = generate_id();
    device["SubcircuitDefinition"]["Variables"] = ordered_json::array();
    device["SubcircuitDefinition"]["VariableFile"] = "";
    device["SubcircuitDefinitionID"] = device["SubcircuitDefinition"]["Id"];

    auto columns = core.get_columns();
    auto coreEffectiveArea = core.get_effective_area();

    auto coreEffectiveLengthMinusColumns = core.get_effective_length();
    if (columns.size() > 1) {
        for (auto column : columns) {
            if (column.get_coordinates()[0] >= 0) {
                coreEffectiveLengthMinusColumns -= column.get_height();
            }
        }
    }

    std::vector<std::vector<int>> columnBottomCoordinates; 
    std::vector<std::vector<int>> columnTopCoordinates;
    
    auto acResistanceCoefficientsPerWinding = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(magnetic, temperature);
    // Resolve AUTO mode from settings (fracpole support)
    auto resolvedMode_sb = CircuitSimulatorExporter::resolve_curve_fitting_mode(mode);
    std::vector<FractionalPoleNetwork> fracpoleNets_sb;
    if (resolvedMode_sb == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
        fracpoleNets_sb = CircuitSimulatorExporter::calculate_fracpole_networks_per_winding(magnetic, temperature);
    }
    // Core losses (calculate coefficients for positioning)
    auto coreLossTopology_sb_pos = static_cast<CoreLossTopology>(Settings::GetInstance().get_circuit_simulator_core_loss_topology());
    auto coreResistanceCoefficients_sb = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature, coreLossTopology_sb_pos);
    double leakageInductance = resolve_dimensional_values(LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0]);
    int numberLadderPairElements = acResistanceCoefficientsPerWinding[0].size() / 2 - 1;
    int numberCoreLadderPairElements = (coreLossTopology_sb_pos == CoreLossTopology::ROSANO)
        ? 5  // R, RL(2), RLC(3) branches ~ 5 visual rows
        : static_cast<int>(coreResistanceCoefficients_sb.size() / 2) + 1;

    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        auto column = columns[columnIndex];
        auto gapsInThisColumn = core.find_gaps_by_column(column);
        ordered_json coreChunkJson;
        std::vector<int> columnCoordinates;
        if (column.get_coordinates()[0] == 0 && column.get_coordinates()[2] != 0) {
            columnCoordinates = {static_cast<int>(column.get_coordinates()[2] * _scale), 0};
        }
        else {
            columnCoordinates = {static_cast<int>(column.get_coordinates()[0] * _scale), 0};
        }

        std::vector<int> columnTopCoordinate = {columnCoordinates[0] + 3, -2 - numberLadderPairElements * 6};  // Don't ask
        std::vector<int> columnBottomCoordinate = {columnCoordinates[0] + 3, 4 + numberCoreLadderPairElements * 6};  // Don't ask
        int currentColumnGapHeight = 6;
        for (size_t gapIndex = 0; gapIndex < gapsInThisColumn.size(); ++gapIndex) {
            currentColumnGapHeight += 6;
            columnBottomCoordinate[1] += 6;
        }
        columnBottomCoordinates.push_back(columnBottomCoordinate);
        columnTopCoordinates.push_back(columnTopCoordinate);
    }


    std::vector<int> windingCoordinates = {columnTopCoordinates[0][0] - 2, columnTopCoordinates[0][1] - 6 + numberLadderPairElements * 6};
    windingCoordinates[1] -= 6 * (coil.get_functional_description().size() - 1);
    int numberLeftWindings = 0;
    int numberRightWindings = 0;
    int maximumWindingCoordinate = windingCoordinates[1];
    int maximumLadderCoordinate = windingCoordinates[1];
    int maximumLeftCoordinate = windingCoordinates[0];
    int maximumRightCoordinate = windingCoordinates[0];

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto dcResistanceThisWinding = WindingLosses::calculate_effective_resistance_of_winding(magnetic, windingIndex, 0.1, temperature);
        auto winding = coil.get_functional_description()[windingIndex];
        std::vector<int> coordinates = windingCoordinates;
        ordered_json windingJson;
        ordered_json topPinJson;
        ordered_json bottomPinJson;
        ordered_json acResistorJson;
        std::vector<ordered_json> ladderJsons;
        std::vector<ordered_json> ladderConnectorsJsons;

        if (winding.get_isolation_side() == IsolationSide::PRIMARY) {
            windingJson = create_winding(winding.get_number_turns(), coordinates, 0, winding.get_name());
            coordinates[0] -= numberLeftWindings * 18;
        }
        else {
            windingJson = create_winding(winding.get_number_turns(), coordinates, 180, winding.get_name());
            coordinates[0] += numberRightWindings * 12;
        }

        if (winding.get_isolation_side() == IsolationSide::PRIMARY) {
            if (windingIndex == 0) {
                coordinates[0] -= 6;
                auto leakageInductanceJson = create_inductor(leakageInductance, coordinates, 0, winding.get_name() + " Leakage inductance");
                device["SubcircuitDefinition"]["Devices"].push_back(leakageInductanceJson);
            }
            coordinates[0] -= 6;
            acResistorJson = create_resistor(dcResistanceThisWinding, coordinates, 180, winding.get_name() + " AC resistance");

            std::pair<std::vector<ordered_json>, std::vector<ordered_json>> aux;
            if (resolvedMode_sb == CircuitSimulatorExporterCurveFittingModes::FRACPOLE &&
                windingIndex < fracpoleNets_sb.size()) {
                aux = create_fracpole_ladder(fracpoleNets_sb[windingIndex], coordinates, winding.get_name());
            } else {
                aux = create_ladder(acResistanceCoefficientsPerWinding[windingIndex], coordinates, winding.get_name());
            }
            if (acResistanceCoefficientsPerWinding[windingIndex].size() > 0) {
                coordinates[0] -= 6;
            }
            ladderJsons = aux.first;
            ladderConnectorsJsons = aux.second;

            coordinates[0] -= 2;
            bottomPinJson = create_pin(coordinates, 0, winding.get_name() + " Input");

            if (numberLeftWindings > 0) {
                auto connectorJson = create_connector({windingCoordinates[0] - numberLeftWindings * 18, windingCoordinates[1] + 1}, {windingCoordinates[0], windingCoordinates[1] + 1}, "Connector between DC resistance in winding " + std::to_string(windingIndex));
                device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
            }
            topPinJson = create_pin({windingCoordinates[0] - 2, windingCoordinates[1] + 4}, 0, winding.get_name() + " Output");
        }
        else {
            coordinates[0] += 4;
            acResistorJson = create_resistor(dcResistanceThisWinding, coordinates, 180, winding.get_name() + " AC resistance");

            std::pair<std::vector<ordered_json>, std::vector<ordered_json>> aux;
            if (resolvedMode_sb == CircuitSimulatorExporterCurveFittingModes::FRACPOLE &&
                windingIndex < fracpoleNets_sb.size()) {
                aux = create_fracpole_ladder(fracpoleNets_sb[windingIndex], {coordinates[0] + 12, coordinates[1]}, winding.get_name());
            } else {
                aux = create_ladder(acResistanceCoefficientsPerWinding[windingIndex], {coordinates[0] + 12, coordinates[1]}, winding.get_name());
            }
            ladderJsons = aux.first;
            ladderConnectorsJsons = aux.second;
            if (acResistanceCoefficientsPerWinding[windingIndex].size() > 0) {
                coordinates[0] += 6;
            }

            coordinates[0] += 6;
            topPinJson = create_pin(coordinates, 180, winding.get_name() + " Input");

            if (numberRightWindings > 0) {
                auto connectorJson = create_connector({windingCoordinates[0] + numberRightWindings * 12 + 4, windingCoordinates[1] + 1}, {windingCoordinates[0] + 4, windingCoordinates[1] + 1}, "Connector between DC resistance in winding " + std::to_string(windingIndex));
                device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
            }
            bottomPinJson = create_pin({windingCoordinates[0] + 4, windingCoordinates[1] + 4}, 180, winding.get_name() + " Output");
        }

        std::vector<int> connectorTopCoordinates = {windingCoordinates[0] + 2, windingCoordinates[1] + 5};
        std::vector<int> connectorBottomCoordinates = {windingCoordinates[0] + 2, windingCoordinates[1] + 7};
        if (windingIndex == coil.get_functional_description().size() - 1) {
            connectorBottomCoordinates[1] = windingCoordinates[1] + 6;
        }

        auto connectorJson = create_connector(connectorBottomCoordinates, connectorTopCoordinates, "Connector between winding " + std::to_string(windingIndex) + " and winding " + std::to_string(windingIndex + 1));
        device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);

        device["SubcircuitDefinition"]["Devices"].push_back(windingJson);
        device["SubcircuitDefinition"]["Devices"].push_back(topPinJson);
        device["SubcircuitDefinition"]["Devices"].push_back(bottomPinJson);
        device["SubcircuitDefinition"]["Devices"].push_back(acResistorJson);
        for (auto ladderJson : ladderJsons) {
            device["SubcircuitDefinition"]["Devices"].push_back(ladderJson);
        }
        for (auto ladderConnectorsJson : ladderConnectorsJsons) {
            device["SubcircuitDefinition"]["Connectors"].push_back(ladderConnectorsJson);
        }

        windingCoordinates[1] += 6;
        if (winding.get_isolation_side() == IsolationSide::PRIMARY) {
            numberLeftWindings++;
        }
        else {
            numberRightWindings++;
        }
    }


    // Magnetizing current and core losses
    {
        auto& settings_sb = Settings::GetInstance();
        auto coreLossTopology_sb = static_cast<CoreLossTopology>(settings_sb.get_circuit_simulator_core_loss_topology());
        auto coreResistanceCoefficients = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature, coreLossTopology_sb);
        // For Rosano, element count is 3 branches (R, RL, RLC) = ~4 visual rows; for Ridley, pairs of RL
        int numberCoreLadderPairElements = (coreLossTopology_sb == CoreLossTopology::ROSANO)
            ? 4
            : static_cast<int>(coreResistanceCoefficients.size() / 2);

        std::vector<int> coordinates = windingCoordinates;
        std::vector<ordered_json> ladderJsons;
        std::vector<ordered_json> ladderConnectorsJsons;
        coordinates[1] += numberCoreLadderPairElements * 6 - 5;

        {
            std::vector<int> connectorTopCoordinates = {windingCoordinates[0] + 2, windingCoordinates[1]};
            std::vector<int> connectorBottomCoordinates = {coordinates[0] + 2, coordinates[1] + 5};

            auto connectorJson = create_connector(connectorTopCoordinates, connectorBottomCoordinates, "Central column connector to core losses");
            device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
        }

        if (!coreResistanceCoefficients.empty()) {
            coordinates[0] -= 6;
            std::pair<std::vector<ordered_json>, std::vector<ordered_json>> aux;
            if (coreLossTopology_sb == CoreLossTopology::ROSANO) {
                aux = create_rosano_ladder(coreResistanceCoefficients, coordinates, "Core losses");
            } else {
                aux = create_ladder(coreResistanceCoefficients, coordinates, "Core losses");
            }
            if (coreResistanceCoefficients.size() > 0) {
                coordinates[0] -= 6;
            }
            ladderJsons = aux.first;
            ladderConnectorsJsons = aux.second;

            for (auto ladderJson : ladderJsons) {
                device["SubcircuitDefinition"]["Devices"].push_back(ladderJson);
            }
            for (auto ladderConnectorsJson : ladderConnectorsJsons) {
                device["SubcircuitDefinition"]["Connectors"].push_back(ladderConnectorsJson);
            }
        }
    }

    for (auto deviceJson : device["SubcircuitDefinition"]["Devices"]) {
        maximumLadderCoordinate = std::min(maximumLadderCoordinate, deviceJson["Top"].get<int>());
        maximumLeftCoordinate = std::max(maximumLeftCoordinate, deviceJson["Left"].get<int>());
        maximumRightCoordinate = std::min(maximumRightCoordinate, deviceJson["Left"].get<int>());
    }

    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        columnTopCoordinates[columnIndex][1] = maximumLadderCoordinate - 5;
    }
    if (columns.size() > 1) {
        columnTopCoordinates[1][0] = maximumLeftCoordinate + 6;
        columnBottomCoordinates[1][0] = maximumLeftCoordinate + 6;
    }
    if (columns.size() > 2) {
        columnTopCoordinates[2][0] = maximumRightCoordinate - 2;
        columnBottomCoordinates[2][0] = maximumRightCoordinate - 2;
    }


    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        auto column = columns[columnIndex];
        auto gapsInThisColumn = core.find_gaps_by_column(column);
        ordered_json coreChunkJson;
        std::vector<int> coordinates = columnBottomCoordinates[columnIndex];
        coordinates[1] -= (gapsInThisColumn.size() + 1) * 6 - 2;
        coordinates[0] -= 3;
        if (columnIndex == 0) {
            coreChunkJson = create_core(core.get_initial_permeability(), coordinates, coreEffectiveArea, coreEffectiveLengthMinusColumns, 90, "Core winding column and plates " + std::to_string(columnIndex));
        }
        else {
            coreChunkJson = create_core(core.get_initial_permeability(), coordinates, coreEffectiveArea, column.get_height(), 90, "Core lateral column " + std::to_string(columnIndex));

            std::vector<int> connectorTopCoordinates = {columnTopCoordinates[0][0], columnTopCoordinates[columnIndex][1]};
            std::vector<int> connectorBottomCoordinates = {columnTopCoordinates[columnIndex][0], coordinates[1] - 2};
            auto connectorJson = create_connector(connectorTopCoordinates, connectorBottomCoordinates, "Connector between column " + std::to_string(columnIndex) + " and top");
            device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
        }
        device["SubcircuitDefinition"]["Devices"].push_back(coreChunkJson);

        for (size_t gapIndex = 0; gapIndex < gapsInThisColumn.size(); ++gapIndex) {
            coordinates[1] += 6;
            auto gap = gapsInThisColumn[gapIndex];
            if (!gap.get_coordinates()) {
                throw GapException("Gap is not processed");
            }
            std::vector<int> gapCoordinates = {coordinates[0], coordinates[1]};

            if (gap.get_length() > 0) {

                auto gapJson = create_air_gap(gapCoordinates, gap.get_area().value(), gap.get_length(), 90, "Column " + std::to_string(columnIndex) + " gap " + std::to_string(gapIndex));
                device["SubcircuitDefinition"]["Devices"].push_back(gapJson);
            }
            else {
                std::vector<int> zeroGapConnectorTopCoordinates = {gapCoordinates[0] + 3, gapCoordinates[1]- 2};
                std::vector<int> zeroGapConnectorBottomCoordinates = {gapCoordinates[0] + 3, gapCoordinates[1] + 4};
                auto connectorJson = create_connector(zeroGapConnectorTopCoordinates, zeroGapConnectorBottomCoordinates, "Connector replacing gap of 0 length");
                device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
            }

        }
    }

    {
        std::vector<int> finalConnectorTopCoordinates = {columnTopCoordinates[0][0], columnTopCoordinates[0][1]};
        std::vector<int> finalConnectorBottomCoordinates = {columnTopCoordinates[0][0], maximumWindingCoordinate + 1};
        auto connectorJson = create_connector(finalConnectorTopCoordinates, finalConnectorBottomCoordinates, "Connector between winding 0 and top");
        device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
    }

    // Magnetic ground
    {
        std::vector<int> columnBottomCoordinatesAux = {0, columnTopCoordinates[0][1]};
        columnBottomCoordinatesAux[0] += 2;
        columnBottomCoordinatesAux[1] -= 2;
        auto groundJson = create_magnetic_ground(columnBottomCoordinatesAux, 180, "Magnetic ground");
        device["SubcircuitDefinition"]["Devices"].push_back(groundJson);
    }

    // Horizontal bottom connectors
    if (columns.size() == 1) {
        auto columnBottomCoordinatesAux = columnBottomCoordinates[0];
        columnBottomCoordinatesAux[1] = 0;
        columnBottomCoordinatesAux[0] += _modelSize / 2;
        auto bottomConnectorJson = create_connector(columnBottomCoordinates[0], columnBottomCoordinatesAux, "Bottom Connector between column " + std::to_string(0) + " and middle");
        device["SubcircuitDefinition"]["Connectors"].push_back(bottomConnectorJson);
        auto topConnectorJson = create_connector(columnTopCoordinates[0], columnBottomCoordinatesAux, "Top Connector between column " + std::to_string(0) + " and middle");
        device["SubcircuitDefinition"]["Connectors"].push_back(topConnectorJson);
        device["SubcircuitDefinition"]["Connectors"] = merge_connectors(device["SubcircuitDefinition"]["Connectors"]);
    }
    else {
        for (size_t columnIndex = 1; columnIndex < columns.size(); ++columnIndex) {
            auto connectorJson = create_connector(columnBottomCoordinates[0], columnBottomCoordinates[columnIndex], "Bottom Connector between column " + std::to_string(0) + " and columm " + std::to_string(columnIndex));
            device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
        }
    }
    device["SubcircuitDefinition"]["Connectors"] = merge_connectors(device["SubcircuitDefinition"]["Connectors"]);

    library["Devices"] = ordered_json::array();
    library["Devices"].push_back(device);
    simulation["Libraries"].push_back(library);
    return simulation.dump(2);
}

} // namespace OpenMagnetics
