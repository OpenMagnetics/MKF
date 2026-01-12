#include "support/Settings.h"
#include "TestingUtils.h"
#include "physical_models/MagnetizingInductance.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

bool verboseTests = false;


namespace OpenMagneticsTesting {


OpenMagnetics::Coil get_quick_coil(std::vector<int64_t> numberTurns,
                                                std::vector<int64_t> numberParallels,
                                                std::string shapeName,
                                                uint8_t interleavingLevel,
                                                WindingOrientation windingOrientation,
                                                WindingOrientation layersOrientation,
                                                CoilAlignment turnsAlignment,
                                                CoilAlignment sectionsAlignment,
                                                std::vector<OpenMagnetics::Wire> wires,
                                                bool useBobbin,
                                                int numberStacks){
    json coilJson;
    coilJson["functionalDescription"] = json::array();

    auto core = get_quick_core(shapeName, json::parse("[]"), numberStacks, "Dummy");
    bool auxUseBobbin = useBobbin;
    if (core.get_shape_family() == CoreShapeFamily::T) {
        auxUseBobbin = false;
    }
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, !auxUseBobbin);
    json bobbinJson;
    to_json(bobbinJson, bobbin);

    coilJson["bobbin"] = bobbinJson;

    for (size_t i = 0; i < numberTurns.size(); ++i){
        json isolationSideJson;
        to_json(isolationSideJson, get_isolation_side_from_index(i));
        json individualcoilJson;
        individualcoilJson["name"] = "winding " + std::to_string(i);
        individualcoilJson["numberTurns"] = numberTurns[i];
        individualcoilJson["numberParallels"] = numberParallels[i];
        individualcoilJson["isolationSide"] = isolationSideJson;
        if (i < wires.size()) {
            individualcoilJson["wire"] = wires[i];
        }
        else {
            individualcoilJson["wire"] = "Round 0.475 - Grade 1";
        }
        coilJson["functionalDescription"].push_back(individualcoilJson);
    }

    OpenMagnetics::Coil coil(coilJson, interleavingLevel, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    return coil;
}

OpenMagnetics::Coil get_quick_coil(std::vector<int64_t> numberTurns,
                                                std::vector<int64_t> numberParallels,
                                                double bobbinHeight,
                                                double bobbinWidth,
                                                std::vector<double> bobbinCenterCoodinates,
                                                uint8_t interleavingLevel,
                                                WindingOrientation windingOrientation,
                                                WindingOrientation layersOrientation,
                                                CoilAlignment turnsAlignment,
                                                CoilAlignment sectionsAlignment,
                                                std::vector<OpenMagnetics::Wire> wires){
    json coilJson;
    coilJson["functionalDescription"] = json::array();

    coilJson["bobbin"] = json();
    coilJson["bobbin"]["processedDescription"] = json();
    coilJson["bobbin"]["processedDescription"]["wallThickness"] = 0.001;
    coilJson["bobbin"]["processedDescription"]["columnThickness"] = 0.001;
    coilJson["bobbin"]["processedDescription"]["columnShape"] = ColumnShape::ROUND;
    coilJson["bobbin"]["processedDescription"]["columnDepth"] = bobbinCenterCoodinates[0] - bobbinWidth / 2;
    coilJson["bobbin"]["processedDescription"]["windingWindows"] = json::array();

    json windingWindow;
    windingWindow["height"] = bobbinHeight;
    windingWindow["width"] = bobbinWidth;
    windingWindow["coordinates"] = json::array();
    for (auto& coord : bobbinCenterCoodinates) {
        windingWindow["coordinates"].push_back(coord);
    }
    coilJson["bobbin"]["processedDescription"]["windingWindows"].push_back(windingWindow);

    for (size_t i = 0; i < numberTurns.size(); ++i){
        json isolationSideJson;
        to_json(isolationSideJson, get_isolation_side_from_index(i));
        json individualcoilJson;
        individualcoilJson["name"] = "winding " + std::to_string(i);
        individualcoilJson["numberTurns"] = numberTurns[i];
        individualcoilJson["numberParallels"] = numberParallels[i];
        individualcoilJson["isolationSide"] = isolationSideJson;
        if (i < wires.size()) {
            individualcoilJson["wire"] = wires[i];
        }
        else {
            individualcoilJson["wire"] = "Round 0.475 - Grade 1";
        }
        coilJson["functionalDescription"].push_back(individualcoilJson);
    }

    OpenMagnetics::Coil coil(coilJson, interleavingLevel, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    return coil;
}

OpenMagnetics::Coil get_quick_coil_no_compact(std::vector<int64_t> numberTurns,
                                                           std::vector<int64_t> numberParallels,
                                                           double bobbinHeight,
                                                           double bobbinWidth,
                                                           std::vector<double> bobbinCenterCoodinates,
                                                           uint8_t interleavingLevel,
                                                           WindingOrientation windingOrientation,
                                                           WindingOrientation layersOrientation,
                                                           CoilAlignment turnsAlignment,
                                                           CoilAlignment sectionsAlignment,
                                                           std::vector<OpenMagnetics::Wire> wires){
    json coilJson;
    coilJson["functionalDescription"] = json::array();

    coilJson["bobbin"] = json();
    coilJson["bobbin"]["processedDescription"] = json();
    coilJson["bobbin"]["processedDescription"]["wallThickness"] = 0.001;
    coilJson["bobbin"]["processedDescription"]["columnThickness"] = 0.001;
    coilJson["bobbin"]["processedDescription"]["columnShape"] = ColumnShape::ROUND;
    coilJson["bobbin"]["processedDescription"]["columnDepth"] = bobbinCenterCoodinates[0] - bobbinWidth / 2;
    coilJson["bobbin"]["processedDescription"]["windingWindows"] = json::array();

    json windingWindow;
    windingWindow["height"] = bobbinHeight;
    windingWindow["width"] = bobbinWidth;
    windingWindow["shape"] = WindingWindowShape::RECTANGULAR;
    windingWindow["coordinates"] = json::array();
    for (auto& coord : bobbinCenterCoodinates) {
        windingWindow["coordinates"].push_back(coord);
    }
    coilJson["bobbin"]["processedDescription"]["windingWindows"].push_back(windingWindow);

    for (size_t i = 0; i < numberTurns.size(); ++i){
        json isolationSideJson;
        to_json(isolationSideJson, get_isolation_side_from_index(i));
        json individualcoilJson;
        individualcoilJson["name"] = "winding " + std::to_string(i);
        individualcoilJson["numberTurns"] = numberTurns[i];
        individualcoilJson["numberParallels"] = numberParallels[i];
        individualcoilJson["isolationSide"] = "primary";
        if (i < wires.size()) {
            individualcoilJson["wire"] = wires[i];
        }
        else {
            individualcoilJson["wire"] = "Round 0.475 - Grade 1";
        }
        coilJson["functionalDescription"].push_back(individualcoilJson);
    }

    settings.set_coil_delimit_and_compact(false);

    OpenMagnetics::Coil coil(coilJson, interleavingLevel, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    return coil;
}

OpenMagnetics::Coil get_quick_toroidal_coil_no_compact(std::vector<int64_t> numberTurns,
                                                              std::vector<int64_t> numberParallels,
                                                              double bobbinRadialHeight,
                                                              double bobbinAngle,
                                                              double columnDepth,
                                                              uint8_t interleavingLevel,
                                                              WindingOrientation windingOrientation,
                                                              WindingOrientation layersOrientation,
                                                              CoilAlignment turnsAlignment,
                                                              CoilAlignment sectionsAlignment,
                                                              std::vector<OpenMagnetics::Wire> wires){
    json coilJson;
    coilJson["functionalDescription"] = json::array();

    coilJson["bobbin"] = json();
    coilJson["bobbin"]["processedDescription"] = json();
    coilJson["bobbin"]["processedDescription"]["wallThickness"] = 0.0;
    coilJson["bobbin"]["processedDescription"]["columnThickness"] = 0.0;
    coilJson["bobbin"]["processedDescription"]["columnShape"] = ColumnShape::ROUND;
    coilJson["bobbin"]["processedDescription"]["columnDepth"] = columnDepth;
    coilJson["bobbin"]["processedDescription"]["windingWindows"] = json::array();

    json windingWindow;
    windingWindow["radialHeight"] = bobbinRadialHeight;
    windingWindow["angle"] = bobbinAngle;
    windingWindow["shape"] = WindingWindowShape::ROUND;
    windingWindow["coordinates"] = json::array({0, 0, 0});
    coilJson["bobbin"]["processedDescription"]["windingWindows"].push_back(windingWindow);

    for (size_t i = 0; i < numberTurns.size(); ++i){
        json isolationSideJson;
        to_json(isolationSideJson, get_isolation_side_from_index(i));
        json individualcoilJson;
        individualcoilJson["name"] = "winding " + std::to_string(i);
        individualcoilJson["numberTurns"] = numberTurns[i];
        individualcoilJson["numberParallels"] = numberParallels[i];
        individualcoilJson["isolationSide"] = isolationSideJson;
        if (i < wires.size()) {
            individualcoilJson["wire"] = wires[i];
        }
        else {
            individualcoilJson["wire"] = "Round 0.475 - Grade 1";
        }
        coilJson["functionalDescription"].push_back(individualcoilJson);
    }

    settings.set_coil_delimit_and_compact(false);

    OpenMagnetics::Coil coil(coilJson, interleavingLevel, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    return coil;
}

Core get_quick_core(std::string shapeName,
                                          json basicGapping,
                                          int numberStacks,
                                          std::string materialName) {
    auto coreJson = json();

    std::string coreType;
    if (shapeName[0] == 'T' || (shapeName[0] == 'R' && shapeName[1] == ' ')) {
        coreType = "toroidal";
    }
    else {
        coreType = "two-piece set";
    }

    coreJson["functionalDescription"] = json();
    coreJson["functionalDescription"]["name"] = "GapReluctanceTest";
    coreJson["functionalDescription"]["type"] = coreType;
    coreJson["functionalDescription"]["material"] = materialName;
    coreJson["functionalDescription"]["shape"] = shapeName;
    coreJson["functionalDescription"]["gapping"] = basicGapping;
    coreJson["functionalDescription"]["numberStacks"] = numberStacks;

    Core core(coreJson);

    return core;
}

OpenMagnetics::Magnetic get_quick_magnetic(std::string shapeName,
                                                  json basicGapping,
                                                  std::vector<int64_t> numberTurns,
                                                  int numberStacks,
                                                  std::string materialName) {
    auto coreJson = json();

    std::string coreType;
    if (shapeName[0] == 'T') {
        coreType = "toroidal";
    }
    else {
        coreType = "two-piece set";
    }

    coreJson["functionalDescription"] = json();
    coreJson["functionalDescription"]["name"] = "GapReluctanceTest";
    coreJson["functionalDescription"]["type"] = coreType;
    coreJson["functionalDescription"]["material"] = materialName;
    coreJson["functionalDescription"]["shape"] = shapeName;
    coreJson["functionalDescription"]["gapping"] = basicGapping;
    coreJson["functionalDescription"]["numberStacks"] = numberStacks;

    Core core(coreJson);
    OpenMagnetics::Coil coil = get_quick_coil(numberTurns,
                                                     std::vector<int64_t>(numberTurns.size(), 1),
                                                     shapeName,
                                                     1);

    coil.set_sections_description(std::nullopt);
    coil.set_layers_description(std::nullopt);
    coil.set_turns_description(std::nullopt);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    return magnetic;
}

OpenMagnetics::Inputs get_quick_insulation_inputs(DimensionWithTolerance altitude,
                                                         Cti cti,
                                                         InsulationType insulationType,
                                                         DimensionWithTolerance mainSupplyVoltage,
                                                         OvervoltageCategory overvoltageCategory,
                                                         PollutionDegree pollutionDegree,
                                                         std::vector<InsulationStandards> standards,
                                                         double maximumVoltageRms,
                                                         double maximumVoltagePeak,
                                                         double frequency,
                                                         WiringTechnology wiringTechnology){
    OpenMagnetics::Inputs inputs;
    DesignRequirements designRequirements;
    InsulationRequirements insulationRequirements;
    OperatingPoint operatingPoint;
    OperatingPointExcitation excitation;
    SignalDescriptor voltage;
    Processed processedVoltage;

    processedVoltage.set_rms(maximumVoltageRms);
    processedVoltage.set_peak(maximumVoltagePeak);
    voltage.set_processed(processedVoltage);
    excitation.set_frequency(frequency);
    excitation.set_voltage(voltage);
    operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    inputs.get_mutable_operating_points().push_back(operatingPoint);
    insulationRequirements.set_altitude(altitude);
    insulationRequirements.set_cti(cti);
    insulationRequirements.set_insulation_type(insulationType);
    insulationRequirements.set_main_supply_voltage(mainSupplyVoltage);
    insulationRequirements.set_overvoltage_category(overvoltageCategory);
    insulationRequirements.set_pollution_degree(pollutionDegree);
    insulationRequirements.set_standards(standards);
    designRequirements.set_insulation(insulationRequirements);
    designRequirements.set_wiring_technology(wiringTechnology);
    inputs.set_design_requirements(designRequirements);
    return inputs;
}


InsulationRequirements get_quick_insulation_requirements(DimensionWithTolerance altitude,
                                                                        Cti cti,
                                                                        InsulationType insulationType,
                                                                        DimensionWithTolerance mainSupplyVoltage,
                                                                        OvervoltageCategory overvoltageCategory,
                                                                        PollutionDegree pollutionDegree,
                                                                        std::vector<InsulationStandards> standards){
    InsulationRequirements insulationRequirements;

    insulationRequirements.set_altitude(altitude);
    insulationRequirements.set_cti(cti);
    insulationRequirements.set_insulation_type(insulationType);
    insulationRequirements.set_main_supply_voltage(mainSupplyVoltage);
    insulationRequirements.set_overvoltage_category(overvoltageCategory);
    insulationRequirements.set_pollution_degree(pollutionDegree);
    insulationRequirements.set_standards(standards);
    return insulationRequirements;
}

json get_ground_gap(double gapLength) {
    auto constants = Constants();
    auto basicGapping = json::array();
    auto basicCentralGap = json();
    basicCentralGap["type"] = "subtractive";
    basicCentralGap["length"] = gapLength;
    auto basicLateralGap = json();
    basicLateralGap["type"] = "residual";
    basicLateralGap["length"] = constants.residualGap;
    basicGapping.push_back(basicCentralGap);
    basicGapping.push_back(basicLateralGap);
    basicGapping.push_back(basicLateralGap);

    return basicGapping;
}

json get_distributed_gap(double gapLength, int numberGaps) {
    auto constants = Constants();
    auto basicGapping = json::array();
    auto basicCentralGap = json();
    basicCentralGap["type"] = "subtractive";
    basicCentralGap["length"] = gapLength;
    auto basicLateralGap = json();
    basicLateralGap["type"] = "residual";
    basicLateralGap["length"] = constants.residualGap;
    for (int i = 0; i < numberGaps; ++i) {
        basicGapping.push_back(basicCentralGap);
    }
    basicGapping.push_back(basicLateralGap);
    basicGapping.push_back(basicLateralGap);

    return basicGapping;
}

json get_spacer_gap(double gapLength) {
    auto basicGapping = json::array();
    auto basicSpacerGap = json();
    basicSpacerGap["type"] = "additive";
    basicSpacerGap["length"] = gapLength;
    basicGapping.push_back(basicSpacerGap);
    basicGapping.push_back(basicSpacerGap);
    basicGapping.push_back(basicSpacerGap);

    return basicGapping;
}

json get_residual_gap() {
    auto constants = Constants();
    auto basicGapping = json::array();
    auto basicCentralGap = json();
    basicCentralGap["type"] = "residual";
    basicCentralGap["length"] = constants.residualGap;
    auto basicLateralGap = json();
    basicLateralGap["type"] = "residual";
    basicLateralGap["length"] = constants.residualGap;
    basicGapping.push_back(basicCentralGap);
    basicGapping.push_back(basicLateralGap);
    basicGapping.push_back(basicLateralGap);

    return basicGapping;
}

void print(std::vector<double> data) {
    for (auto i : data) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
}
void print(std::vector<std::vector<double>> data) {
    for (auto i : data) {
        print(i);
    }
}

void print(std::vector<int64_t> data) {
    for (auto i : data) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
}

void print(std::vector<uint64_t> data) {
    for (auto i : data) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
}

void print(std::vector<std::string> data) {
    for (auto i : data) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
}

void print(double data) {
    std::cout << data << std::endl;
}
void print(std::string data) {
    std::cout << data << std::endl;
}
void print_json(json data) {
    std::cout << data << std::endl;
}

void check_sections_description(OpenMagnetics::Coil coil,
                                std::vector<int64_t> numberTurns,
                                std::vector<int64_t> numberParallels,
                                uint8_t interleavingLevel,
                                WindingOrientation windingOrientation) {

    auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
    auto windingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
    auto bobbinArea = windingWindow.get_width().value() * windingWindow.get_height().value();
    auto sectionsDescription = coil.get_sections_description().value();
    std::vector<double> numberAssignedParallels(numberTurns.size(), 0);
    std::vector<int64_t> numberAssignedPhysicalTurns(numberTurns.size(), 0);
    std::map<std::string, std::vector<double>> dimensionsByName;
    std::map<std::string, std::vector<double>> coordinatesByName;
    double sectionsArea = 0;
    size_t numberInsulationSections = 0;

    for (auto& section : sectionsDescription){
        if(section.get_type() == ElectricalType::INSULATION) {
            numberInsulationSections++;
            sectionsArea += section.get_dimensions()[0] * section.get_dimensions()[1];
        }
        else {
            for (auto& partialWinding : section.get_partial_windings()) {
                auto currentIndividualWindingIndex = coil.get_winding_index_by_name(partialWinding.get_winding());
                auto currentIndividualWinding = coil.get_winding_by_name(partialWinding.get_winding());
                sectionsArea += section.get_dimensions()[0] * section.get_dimensions()[1];
                dimensionsByName[section.get_name()] = section.get_dimensions();
                coordinatesByName[section.get_name()] = section.get_coordinates();
                REQUIRE(roundFloat(section.get_coordinates()[0] - section.get_dimensions()[0] / 2, 6) >= roundFloat(windingWindow.get_coordinates().value()[0] - windingWindow.get_width().value() / 2, 6));
                REQUIRE(roundFloat(section.get_coordinates()[0] + section.get_dimensions()[0] / 2, 6) <= roundFloat(windingWindow.get_coordinates().value()[0] + windingWindow.get_width().value() / 2, 6));
                REQUIRE(roundFloat(section.get_coordinates()[1] - section.get_dimensions()[1] / 2, 6) >= roundFloat(windingWindow.get_coordinates().value()[1] - windingWindow.get_height().value() / 2, 6));
                REQUIRE(roundFloat(section.get_coordinates()[1] + section.get_dimensions()[1] / 2, 6) <= roundFloat(windingWindow.get_coordinates().value()[1] + windingWindow.get_height().value() / 2, 6));

                for (auto& parallelProportion : partialWinding.get_parallels_proportion()) {
                    numberAssignedParallels[currentIndividualWindingIndex] += parallelProportion;
                    numberAssignedPhysicalTurns[currentIndividualWindingIndex] += round(parallelProportion * currentIndividualWinding.get_number_turns());
                }
            }
            REQUIRE(section.get_filling_factor().value() > 0);
        }
    }
    for (size_t i = 0; i < sectionsDescription.size() - 1; ++i){
        if(sectionsDescription[i].get_type() == ElectricalType::INSULATION) {
        }
        else {
            if (windingOrientation == WindingOrientation::OVERLAPPING) {
                REQUIRE(sectionsDescription[i].get_coordinates()[0] < sectionsDescription[i + 1].get_coordinates()[0]);
                REQUIRE(sectionsDescription[i].get_coordinates()[1] == sectionsDescription[i + 1].get_coordinates()[1]);
            } 
            else if (windingOrientation == WindingOrientation::CONTIGUOUS) {
                REQUIRE(sectionsDescription[i].get_coordinates()[1] > sectionsDescription[i + 1].get_coordinates()[1]);
                REQUIRE(sectionsDescription[i].get_coordinates()[0] == sectionsDescription[i + 1].get_coordinates()[0]);
            }
        }
    }

    REQUIRE(roundFloat(bobbinArea, 6) == roundFloat(sectionsArea, 6));
    for (size_t i = 0; i < numberAssignedParallels.size() - 1; ++i){
        REQUIRE(round(numberAssignedParallels[i]) == round(numberParallels[i]));
        REQUIRE(round(numberAssignedPhysicalTurns[i]) == round(numberTurns[i] * numberParallels[i]));
    }
    REQUIRE(sectionsDescription.size() - numberInsulationSections == (interleavingLevel * numberTurns.size()));
    REQUIRE(!check_collisions(dimensionsByName, coordinatesByName));
}

void check_layers_description(OpenMagnetics::Coil coil,
                              WindingOrientation layersOrientation) {
    if (!coil.get_layers_description()) {
        return;
    }
    auto sections = coil.get_sections_description().value();
    std::map<std::string, std::vector<double>> dimensionsByName;
    std::map<std::string, std::vector<double>> coordinatesByName;

    for (auto& section : sections){
        auto layers = coil.get_layers_by_section(section.get_name());
        if(section.get_type() == ElectricalType::INSULATION) {
        }
        else {
            auto sectionParallelsProportionExpected = section.get_partial_windings()[0].get_parallels_proportion();
            std::vector<double> sectionParallelsProportion(sectionParallelsProportionExpected.size(), 0);
            for (auto& layer : layers){
                for (size_t i = 0; i < sectionParallelsProportion.size(); ++i){
                    sectionParallelsProportion[i] += layer.get_partial_windings()[0].get_parallels_proportion()[i];
                }
                REQUIRE(layer.get_filling_factor().value() > 0);

                dimensionsByName[layer.get_name()] = layer.get_dimensions();
                coordinatesByName[layer.get_name()] = layer.get_coordinates();
            }
            for (size_t i = 0; i < sectionParallelsProportion.size(); ++i){
                REQUIRE(roundFloat(sectionParallelsProportion[i], 9) == roundFloat(sectionParallelsProportionExpected[i], 9));
            }
            for (size_t i = 0; i < layers.size() - 1; ++i){
                if (layersOrientation == WindingOrientation::OVERLAPPING) {
                    REQUIRE(layers[i].get_coordinates()[0] < layers[i + 1].get_coordinates()[0]);
                    REQUIRE(layers[i].get_coordinates()[1] == layers[i + 1].get_coordinates()[1]);
                    REQUIRE(layers[i].get_coordinates()[2] == layers[i + 1].get_coordinates()[2]);
                } 
                else if (layersOrientation == WindingOrientation::CONTIGUOUS) {
                    REQUIRE(layers[i].get_coordinates()[1] > layers[i + 1].get_coordinates()[1]);
                    REQUIRE(layers[i].get_coordinates()[0] == layers[i + 1].get_coordinates()[0]);
                    REQUIRE(layers[i].get_coordinates()[2] == layers[i + 1].get_coordinates()[2]);
                }
            }
        }

    }

    REQUIRE(!check_collisions(dimensionsByName, coordinatesByName));
}


bool check_turns_description(OpenMagnetics::Coil coil) {
    if (!coil.get_turns_description()) {
        return true;
    }

    std::vector<std::vector<double>> parallelProportion;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        parallelProportion.push_back(std::vector<double>(coil.get_number_parallels(windingIndex), 0));
    }

    auto wires = coil.get_wires();

    auto turns = coil.get_turns_description().value();
    std::map<std::string, std::vector<double>> dimensionsByName;
    std::map<std::string, std::vector<double>> coordinatesByName;
    std::map<std::string, std::vector<double>> additionalCoordinatesByName;

    auto bobbin = coil.resolve_bobbin();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (bobbinWindingWindowShape == WindingWindowShape::ROUND) {
        coil.convert_turns_to_cartesian_coordinates();
    }


    int turnsIn0 = 0;
    for (auto& turn : turns){
        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        if (windingIndex == 0 && turn.get_parallel() == 0) {
            turnsIn0++;
        }
        parallelProportion[windingIndex][turn.get_parallel()] += 1.0 / coil.get_number_turns(windingIndex);

        if (bobbinWindingWindowShape != WindingWindowShape::ROUND || wires[windingIndex].get_type() != WireType::RECTANGULAR) {
            dimensionsByName[turn.get_name()] = turn.get_dimensions().value();
        }
        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            coordinatesByName[turn.get_name()] = turn.get_coordinates();
        }
        else {
            double xCoordinate = turn.get_coordinates()[0];
            double yCoordinate = turn.get_coordinates()[1];
            if (wires[windingIndex].get_type() != WireType::RECTANGULAR) {
                coordinatesByName[turn.get_name()] = {xCoordinate, yCoordinate};
            }
            if (turn.get_additional_coordinates()) {
                auto additionalCoordinates = turn.get_additional_coordinates().value();

                for (auto additionalCoordinate : additionalCoordinates){
                    double xAdditionalCoordinate = additionalCoordinate[0];
                    double yAdditionalCoordinate = additionalCoordinate[1];
                    additionalCoordinatesByName[turn.get_name()] = {xAdditionalCoordinate, yAdditionalCoordinate};
                }
            }
        }
    }

    bool equalToOne = true;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        for (size_t parallelIndex = 0; parallelIndex < coil.get_number_parallels(windingIndex); ++parallelIndex) {
            equalToOne &= roundFloat(parallelProportion[windingIndex][parallelIndex], 9) == 1;
        }
    }

    REQUIRE(equalToOne);
    bool collides = check_collisions(dimensionsByName, coordinatesByName, bobbinWindingWindowShape == WindingWindowShape::ROUND);
    REQUIRE(!collides);
    if (additionalCoordinatesByName.size() > 0) {
        collides |= check_collisions(dimensionsByName, additionalCoordinatesByName, bobbinWindingWindowShape == WindingWindowShape::ROUND);
    }
    REQUIRE(!collides);
    return !collides && equalToOne;
}

bool check_wire_standards(OpenMagnetics::Coil coil) {
    auto wires = coil.get_wires();
    std::optional<WireStandard> firstWireStandard = std::nullopt;
    for (auto wire : wires) {
        if (wire.get_standard()) {
            if (!firstWireStandard) {
                firstWireStandard = wire.get_standard().value();
            }
            else {
                REQUIRE(firstWireStandard.value() == wire.get_standard().value());
                if (firstWireStandard.value() != wire.get_standard().value()) {
                    return false;
                }
            }
        }
    }
    return true;
}

void check_winding_losses(OpenMagnetics::Mas mas) {
    for (size_t operatingPointIndex = 0; operatingPointIndex < mas.get_outputs().size(); ++operatingPointIndex) {
        auto output = mas.get_outputs()[operatingPointIndex];
        auto windingLosses = output.get_winding_losses().value();
        double totalWindingLosses = windingLosses.get_winding_losses();
        double totalWindingLossesByTurn = 0;
        double totalWindingLossesByLayer = 0;
        double totalWindingLossesBySection = 0;
        double totalWindingLossesByWinding = 0;

        auto windingLossesPerTurn = windingLosses.get_winding_losses_per_turn().value();
        auto windingLossesPerLayer = windingLosses.get_winding_losses_per_layer().value();
        auto windingLossesPerSection = windingLosses.get_winding_losses_per_section().value();
        auto windingLossesPerWinding = windingLosses.get_winding_losses_per_winding().value();

        {
            double ohmicLossesPerTurn = 0;
            double ohmicLossesPerLayer = 0;
            double ohmicLossesPerSection = 0;
            double ohmicLossesPerWinding = 0;
            for (auto lossesElement : windingLossesPerTurn) {
                ohmicLossesPerTurn += lossesElement.get_ohmic_losses()->get_losses();
            }
            for (auto lossesElement : windingLossesPerLayer) {
                ohmicLossesPerLayer += lossesElement.get_ohmic_losses()->get_losses();
            }
            for (auto lossesElement : windingLossesPerSection) {
                ohmicLossesPerSection += lossesElement.get_ohmic_losses()->get_losses();
            }
            for (auto lossesElement : windingLossesPerWinding) {
                ohmicLossesPerWinding += lossesElement.get_ohmic_losses()->get_losses();
            }
            REQUIRE_THAT(ohmicLossesPerTurn, Catch::Matchers::WithinAbs(ohmicLossesPerLayer, ohmicLossesPerTurn * 0.001));
            REQUIRE_THAT(ohmicLossesPerTurn, Catch::Matchers::WithinAbs(ohmicLossesPerSection, ohmicLossesPerTurn * 0.001));
            REQUIRE_THAT(ohmicLossesPerTurn, Catch::Matchers::WithinAbs(ohmicLossesPerWinding, ohmicLossesPerTurn * 0.001));
            totalWindingLossesByTurn += ohmicLossesPerTurn;
            totalWindingLossesByLayer += ohmicLossesPerLayer;
            totalWindingLossesBySection += ohmicLossesPerSection;
            totalWindingLossesByWinding += ohmicLossesPerWinding;
        }


        {
            double skinLossesPerTurn = 0;
            double skinLossesPerLayer = 0;
            double skinLossesPerSection = 0;
            double skinLossesPerWinding = 0;
            for (auto lossesElement : windingLossesPerTurn) {
                auto harmonicLosses = lossesElement.get_skin_effect_losses()->get_losses_per_harmonic();
                for(auto harmonicLoss : harmonicLosses) {
                    skinLossesPerTurn += harmonicLoss;
                }
            }
            for (auto lossesElement : windingLossesPerLayer) {
                auto harmonicLosses = lossesElement.get_skin_effect_losses()->get_losses_per_harmonic();
                for(auto harmonicLoss : harmonicLosses) {
                    skinLossesPerLayer += harmonicLoss;
                }
            }
            for (auto lossesElement : windingLossesPerSection) {
                auto harmonicLosses = lossesElement.get_skin_effect_losses()->get_losses_per_harmonic();
                for(auto harmonicLoss : harmonicLosses) {
                    skinLossesPerSection += harmonicLoss;
                }
            }
            for (auto lossesElement : windingLossesPerWinding) {
                auto harmonicLosses = lossesElement.get_skin_effect_losses()->get_losses_per_harmonic();
                for(auto harmonicLoss : harmonicLosses) {
                    skinLossesPerWinding += harmonicLoss;
                }
            }
            REQUIRE_THAT(skinLossesPerTurn, Catch::Matchers::WithinAbs(skinLossesPerLayer, skinLossesPerTurn * 0.001));
            REQUIRE_THAT(skinLossesPerTurn, Catch::Matchers::WithinAbs(skinLossesPerSection, skinLossesPerTurn * 0.001));
            REQUIRE_THAT(skinLossesPerTurn, Catch::Matchers::WithinAbs(skinLossesPerWinding, skinLossesPerTurn * 0.001));
            totalWindingLossesByTurn += skinLossesPerTurn;
            totalWindingLossesByLayer += skinLossesPerLayer;
            totalWindingLossesBySection += skinLossesPerSection;
            totalWindingLossesByWinding += skinLossesPerWinding;
        }


        {
            double proximityLossesPerTurn = 0;
            double proximityLossesPerLayer = 0;
            double proximityLossesPerSection = 0;
            double proximityLossesPerWinding = 0;
            for (auto lossesElement : windingLossesPerTurn) {
                auto harmonicLosses = lossesElement.get_proximity_effect_losses()->get_losses_per_harmonic();
                for(auto harmonicLoss : harmonicLosses) {
                    proximityLossesPerTurn += harmonicLoss;
                }
            }
            for (auto lossesElement : windingLossesPerLayer) {
                auto harmonicLosses = lossesElement.get_proximity_effect_losses()->get_losses_per_harmonic();
                for(auto harmonicLoss : harmonicLosses) {
                    proximityLossesPerLayer += harmonicLoss;
                }
            }
            for (auto lossesElement : windingLossesPerSection) {
                auto harmonicLosses = lossesElement.get_proximity_effect_losses()->get_losses_per_harmonic();
                for(auto harmonicLoss : harmonicLosses) {
                    proximityLossesPerSection += harmonicLoss;
                }
            }
            for (auto lossesElement : windingLossesPerWinding) {
                auto harmonicLosses = lossesElement.get_proximity_effect_losses()->get_losses_per_harmonic();
                for(auto harmonicLoss : harmonicLosses) {
                    proximityLossesPerWinding += harmonicLoss;
                }
            }
            REQUIRE_THAT(proximityLossesPerTurn, Catch::Matchers::WithinAbs(proximityLossesPerLayer, proximityLossesPerTurn * 0.001));
            REQUIRE_THAT(proximityLossesPerTurn, Catch::Matchers::WithinAbs(proximityLossesPerSection, proximityLossesPerTurn * 0.001));
            REQUIRE_THAT(proximityLossesPerTurn, Catch::Matchers::WithinAbs(proximityLossesPerWinding, proximityLossesPerTurn * 0.001));
            totalWindingLossesByTurn += proximityLossesPerTurn;
            totalWindingLossesByLayer += proximityLossesPerLayer;
            totalWindingLossesBySection += proximityLossesPerSection;
            totalWindingLossesByWinding += proximityLossesPerWinding;
        }

        REQUIRE_THAT(totalWindingLosses, Catch::Matchers::WithinAbs(totalWindingLossesByTurn, totalWindingLosses * 0.001));
        REQUIRE_THAT(totalWindingLosses, Catch::Matchers::WithinAbs(totalWindingLossesByLayer, totalWindingLosses * 0.001));
        REQUIRE_THAT(totalWindingLosses, Catch::Matchers::WithinAbs(totalWindingLossesBySection, totalWindingLosses * 0.001));
        REQUIRE_THAT(totalWindingLosses, Catch::Matchers::WithinAbs(totalWindingLossesByWinding, totalWindingLosses * 0.001));
    }
}

OpenMagnetics::Mas mas_loader(std::string path) {
    std::ifstream f(path);
    std::string data((std::istreambuf_iterator<char>(f)),
                     std::istreambuf_iterator<char>());
    auto masJson = json::parse(data);
    auto inputsJson = masJson["inputs"];
    auto magneticJson = masJson["magnetic"];
    auto outputsJson = masJson["outputs"];

    auto magnetic = OpenMagnetics::Magnetic(magneticJson);
    if (magneticJson["coil"]["bobbin"] == "Basic") {
        auto bobbinData = OpenMagnetics::Bobbin::create_quick_bobbin(magnetic.get_mutable_core(), false);
        to_json(magneticJson["coil"]["bobbin"], bobbinData);
    }
    auto coil = OpenMagnetics::Coil(magneticJson["coil"]);
    std::vector<OpenMagnetics::Outputs> outputs;
    if (outputsJson != nullptr) {
        outputs = std::vector<OpenMagnetics::Outputs>(outputsJson);
    }
    OpenMagnetics::Inputs inputs;
    std::vector<double> magnetizingInductancePerPoint;
    for (auto output : outputs) {
        if (output.get_magnetizing_inductance()) {
            auto magnetizingInductance = resolve_dimensional_values(output.get_magnetizing_inductance()->get_magnetizing_inductance());
            magnetizingInductancePerPoint.push_back(magnetizingInductance);
        }
    }

    if (magnetizingInductancePerPoint.size() > 0) {
        inputs = OpenMagnetics::Inputs(inputsJson, true, magnetizingInductancePerPoint);

    }
    else {
        try {
            MagnetizingInductance magnetizingInductanceModel;
            double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance().get_nominal().value();
            inputs = OpenMagnetics::Inputs(inputsJson, true, magnetizingInductance);
        }
        catch (const std::exception &e)
        {
            inputs = OpenMagnetics::Inputs(inputsJson, true);
        }
    }

    OpenMagnetics::Mas mas;
    mas.set_inputs(inputs);
    mas.set_magnetic(magnetic);
    mas.set_outputs(outputs);
    return mas;
}



} // namespace OpenMagneticsTesting