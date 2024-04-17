#include "Settings.h"
#include "Utils.h"
#include "TestingUtils.h"
#include <UnitTest++.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

bool verboseTests = true;


namespace OpenMagneticsTesting {


OpenMagnetics::CoilWrapper get_quick_coil(std::vector<int64_t> numberTurns,
                                                std::vector<int64_t> numberParallels,
                                                std::string shapeName,
                                                uint8_t interleavingLevel,
                                                OpenMagnetics::WindingOrientation windingOrientation,
                                                OpenMagnetics::WindingOrientation layersOrientation,
                                                OpenMagnetics::CoilAlignment turnsAlignment,
                                                OpenMagnetics::CoilAlignment sectionsAlignment,
                                                std::vector<OpenMagnetics::WireWrapper> wires,
                                                bool useBobbin){
    json coilJson;
    coilJson["functionalDescription"] = json::array();

    auto core = get_quick_core(shapeName, json::parse("[]"), 1, "Dummy");
    auto bobbin = OpenMagnetics::BobbinWrapper::create_quick_bobbin(core, !useBobbin);
    json bobbinJson;
    OpenMagnetics::to_json(bobbinJson, bobbin);

    coilJson["bobbin"] = bobbinJson;

    for (size_t i = 0; i < numberTurns.size(); ++i){
        json individualcoilJson;
        individualcoilJson["name"] = "winding " + std::to_string(i);
        individualcoilJson["numberTurns"] = numberTurns[i];
        individualcoilJson["numberParallels"] = numberParallels[i];
        individualcoilJson["isolationSide"] = "primary";
        if (i < wires.size()) {
            individualcoilJson["wire"] = wires[i];
        }
        else {
            individualcoilJson["wire"] = "0.475 - Grade 1";
        }
        coilJson["functionalDescription"].push_back(individualcoilJson);
    }

    OpenMagnetics::CoilWrapper coil(coilJson, interleavingLevel, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    return coil;
}

OpenMagnetics::CoilWrapper get_quick_coil(std::vector<int64_t> numberTurns,
                                                std::vector<int64_t> numberParallels,
                                                double bobbinHeight,
                                                double bobbinWidth,
                                                std::vector<double> bobbinCenterCoodinates,
                                                uint8_t interleavingLevel,
                                                OpenMagnetics::WindingOrientation windingOrientation,
                                                OpenMagnetics::WindingOrientation layersOrientation,
                                                OpenMagnetics::CoilAlignment turnsAlignment,
                                                OpenMagnetics::CoilAlignment sectionsAlignment,
                                                std::vector<OpenMagnetics::WireWrapper> wires){
    json coilJson;
    coilJson["functionalDescription"] = json::array();

    coilJson["bobbin"] = json();
    coilJson["bobbin"]["processedDescription"] = json();
    coilJson["bobbin"]["processedDescription"]["wallThickness"] = 0.001;
    coilJson["bobbin"]["processedDescription"]["columnThickness"] = 0.001;
    coilJson["bobbin"]["processedDescription"]["columnShape"] = OpenMagnetics::ColumnShape::ROUND;
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
        json individualcoilJson;
        individualcoilJson["name"] = "winding " + std::to_string(i);
        individualcoilJson["numberTurns"] = numberTurns[i];
        individualcoilJson["numberParallels"] = numberParallels[i];
        individualcoilJson["isolationSide"] = "primary";
        if (i < wires.size()) {
            individualcoilJson["wire"] = wires[i];
        }
        else {
            individualcoilJson["wire"] = "0.475 - Grade 1";
        }
        coilJson["functionalDescription"].push_back(individualcoilJson);
    }

    OpenMagnetics::CoilWrapper coil(coilJson, interleavingLevel, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    return coil;
}

OpenMagnetics::CoilWrapper get_quick_coil_no_compact(std::vector<int64_t> numberTurns,
                                                           std::vector<int64_t> numberParallels,
                                                           double bobbinHeight,
                                                           double bobbinWidth,
                                                           std::vector<double> bobbinCenterCoodinates,
                                                           uint8_t interleavingLevel,
                                                           OpenMagnetics::WindingOrientation windingOrientation,
                                                           OpenMagnetics::WindingOrientation layersOrientation,
                                                           OpenMagnetics::CoilAlignment turnsAlignment,
                                                           OpenMagnetics::CoilAlignment sectionsAlignment,
                                                           std::vector<OpenMagnetics::WireWrapper> wires){
    json coilJson;
    coilJson["functionalDescription"] = json::array();

    coilJson["bobbin"] = json();
    coilJson["bobbin"]["processedDescription"] = json();
    coilJson["bobbin"]["processedDescription"]["wallThickness"] = 0.001;
    coilJson["bobbin"]["processedDescription"]["columnThickness"] = 0.001;
    coilJson["bobbin"]["processedDescription"]["columnShape"] = OpenMagnetics::ColumnShape::ROUND;
    coilJson["bobbin"]["processedDescription"]["columnDepth"] = bobbinCenterCoodinates[0] - bobbinWidth / 2;
    coilJson["bobbin"]["processedDescription"]["windingWindows"] = json::array();

    json windingWindow;
    windingWindow["height"] = bobbinHeight;
    windingWindow["width"] = bobbinWidth;
    windingWindow["shape"] = OpenMagnetics::WindingWindowShape::RECTANGULAR;
    windingWindow["coordinates"] = json::array();
    for (auto& coord : bobbinCenterCoodinates) {
        windingWindow["coordinates"].push_back(coord);
    }
    coilJson["bobbin"]["processedDescription"]["windingWindows"].push_back(windingWindow);

    for (size_t i = 0; i < numberTurns.size(); ++i){
        json individualcoilJson;
        individualcoilJson["name"] = "winding " + std::to_string(i);
        individualcoilJson["numberTurns"] = numberTurns[i];
        individualcoilJson["numberParallels"] = numberParallels[i];
        individualcoilJson["isolationSide"] = "primary";
        if (i < wires.size()) {
            individualcoilJson["wire"] = wires[i];
        }
        else {
            individualcoilJson["wire"] = "0.475 - Grade 1";
        }
        coilJson["functionalDescription"].push_back(individualcoilJson);
    }

    auto settings = OpenMagnetics::Settings::GetInstance();
    settings->set_coil_delimit_and_compact(false);

    OpenMagnetics::CoilWrapper coil(coilJson, interleavingLevel, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    return coil;
}

OpenMagnetics::CoilWrapper get_quick_toroidal_coil_no_compact(std::vector<int64_t> numberTurns,
                                                              std::vector<int64_t> numberParallels,
                                                              double bobbinRadialHeight,
                                                              double bobbinAngle,
                                                              double columnDepth,
                                                              uint8_t interleavingLevel,
                                                              OpenMagnetics::WindingOrientation windingOrientation,
                                                              OpenMagnetics::WindingOrientation layersOrientation,
                                                              OpenMagnetics::CoilAlignment turnsAlignment,
                                                              OpenMagnetics::CoilAlignment sectionsAlignment,
                                                              std::vector<OpenMagnetics::WireWrapper> wires){
    json coilJson;
    coilJson["functionalDescription"] = json::array();

    coilJson["bobbin"] = json();
    coilJson["bobbin"]["processedDescription"] = json();
    coilJson["bobbin"]["processedDescription"]["wallThickness"] = 0.0;
    coilJson["bobbin"]["processedDescription"]["columnThickness"] = 0.0;
    coilJson["bobbin"]["processedDescription"]["columnShape"] = OpenMagnetics::ColumnShape::ROUND;
    coilJson["bobbin"]["processedDescription"]["columnDepth"] = columnDepth;
    coilJson["bobbin"]["processedDescription"]["windingWindows"] = json::array();

    json windingWindow;
    windingWindow["radialHeight"] = bobbinRadialHeight;
    windingWindow["angle"] = bobbinAngle;
    windingWindow["shape"] = OpenMagnetics::WindingWindowShape::ROUND;
    windingWindow["coordinates"] = json::array({0, 0, 0});
    coilJson["bobbin"]["processedDescription"]["windingWindows"].push_back(windingWindow);

    for (size_t i = 0; i < numberTurns.size(); ++i){
        json individualcoilJson;
        individualcoilJson["name"] = "winding " + std::to_string(i);
        individualcoilJson["numberTurns"] = numberTurns[i];
        individualcoilJson["numberParallels"] = numberParallels[i];
        individualcoilJson["isolationSide"] = "primary";
        if (i < wires.size()) {
            individualcoilJson["wire"] = wires[i];
        }
        else {
            individualcoilJson["wire"] = "0.475 - Grade 1";
        }
        coilJson["functionalDescription"].push_back(individualcoilJson);
    }

    auto settings = OpenMagnetics::Settings::GetInstance();
    settings->set_coil_delimit_and_compact(false);

    OpenMagnetics::CoilWrapper coil(coilJson, interleavingLevel, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
    return coil;
}

OpenMagnetics::CoreWrapper get_quick_core(std::string shapeName,
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

    OpenMagnetics::CoreWrapper core(coreJson);

    return core;
}

OpenMagnetics::MagneticWrapper get_quick_magnetic(std::string shapeName,
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

    OpenMagnetics::CoreWrapper core(coreJson);
    OpenMagnetics::CoilWrapper coil = get_quick_coil(numberTurns,
                                                     std::vector<int64_t>(numberTurns.size(), 1),
                                                     shapeName,
                                                     1);

    coil.set_sections_description(std::nullopt);
    coil.set_layers_description(std::nullopt);
    coil.set_turns_description(std::nullopt);
    OpenMagnetics::MagneticWrapper magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    return magnetic;
}

OpenMagnetics::InputsWrapper get_quick_insulation_inputs(OpenMagnetics::DimensionWithTolerance altitude,
                                                         OpenMagnetics::Cti cti,
                                                         OpenMagnetics::InsulationType insulationType,
                                                         OpenMagnetics::DimensionWithTolerance mainSupplyVoltage,
                                                         OpenMagnetics::OvervoltageCategory overvoltageCategory,
                                                         OpenMagnetics::PollutionDegree pollutionDegree,
                                                         std::vector<OpenMagnetics::InsulationStandards> standards,
                                                         double maximumVoltageRms,
                                                         double maximumVoltagePeak,
                                                         double frequency,
                                                         OpenMagnetics::WiringTechnology wiringTechnology){
    OpenMagnetics::InputsWrapper inputs;
    OpenMagnetics::DesignRequirements designRequirements;
    OpenMagnetics::InsulationRequirements insulationRequirements;
    OpenMagnetics::OperatingPoint operatingPoint;
    OpenMagnetics::OperatingPointExcitation excitation;
    OpenMagnetics::SignalDescriptor voltage;
    OpenMagnetics::Processed processedVoltage;

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


OpenMagnetics::InsulationRequirements get_quick_insulation_requirements(OpenMagnetics::DimensionWithTolerance altitude,
                                                                        OpenMagnetics::Cti cti,
                                                                        OpenMagnetics::InsulationType insulationType,
                                                                        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage,
                                                                        OpenMagnetics::OvervoltageCategory overvoltageCategory,
                                                                        OpenMagnetics::PollutionDegree pollutionDegree,
                                                                        std::vector<OpenMagnetics::InsulationStandards> standards){
    OpenMagnetics::InsulationRequirements insulationRequirements;

    insulationRequirements.set_altitude(altitude);
    insulationRequirements.set_cti(cti);
    insulationRequirements.set_insulation_type(insulationType);
    insulationRequirements.set_main_supply_voltage(mainSupplyVoltage);
    insulationRequirements.set_overvoltage_category(overvoltageCategory);
    insulationRequirements.set_pollution_degree(pollutionDegree);
    insulationRequirements.set_standards(standards);
    return insulationRequirements;
}

json get_grinded_gap(double gapLength) {
    auto constants = OpenMagnetics::Constants();
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
    auto constants = OpenMagnetics::Constants();
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
    auto constants = OpenMagnetics::Constants();
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

void check_sections_description(OpenMagnetics::CoilWrapper coil,
                                std::vector<int64_t> numberTurns,
                                std::vector<int64_t> numberParallels,
                                uint8_t interleavingLevel,
                                OpenMagnetics::WindingOrientation windingOrientation) {

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
        if(section.get_type() == OpenMagnetics::ElectricalType::INSULATION) {
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
                CHECK(OpenMagnetics::roundFloat(section.get_coordinates()[0] - section.get_dimensions()[0] / 2, 6) >= OpenMagnetics::roundFloat(windingWindow.get_coordinates().value()[0] - windingWindow.get_width().value() / 2, 6));
                CHECK(OpenMagnetics::roundFloat(section.get_coordinates()[0] + section.get_dimensions()[0] / 2, 6) <= OpenMagnetics::roundFloat(windingWindow.get_coordinates().value()[0] + windingWindow.get_width().value() / 2, 6));
                CHECK(OpenMagnetics::roundFloat(section.get_coordinates()[1] - section.get_dimensions()[1] / 2, 6) >= OpenMagnetics::roundFloat(windingWindow.get_coordinates().value()[1] - windingWindow.get_height().value() / 2, 6));
                CHECK(OpenMagnetics::roundFloat(section.get_coordinates()[1] + section.get_dimensions()[1] / 2, 6) <= OpenMagnetics::roundFloat(windingWindow.get_coordinates().value()[1] + windingWindow.get_height().value() / 2, 6));

                for (auto& parallelProportion : partialWinding.get_parallels_proportion()) {
                    numberAssignedParallels[currentIndividualWindingIndex] += parallelProportion;
                    numberAssignedPhysicalTurns[currentIndividualWindingIndex] += round(parallelProportion * currentIndividualWinding.get_number_turns());
                }
            }
            CHECK(section.get_filling_factor().value() > 0);
        }
    }
    for (size_t i = 0; i < sectionsDescription.size() - 1; ++i){
        if(sectionsDescription[i].get_type() == OpenMagnetics::ElectricalType::INSULATION) {
        }
        else {
            if (windingOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
                CHECK(sectionsDescription[i].get_coordinates()[0] < sectionsDescription[i + 1].get_coordinates()[0]);
                CHECK(sectionsDescription[i].get_coordinates()[1] == sectionsDescription[i + 1].get_coordinates()[1]);
            } 
            else if (windingOrientation == OpenMagnetics::WindingOrientation::CONTIGUOUS) {
                CHECK(sectionsDescription[i].get_coordinates()[1] > sectionsDescription[i + 1].get_coordinates()[1]);
                CHECK(sectionsDescription[i].get_coordinates()[0] == sectionsDescription[i + 1].get_coordinates()[0]);
            }
        }
    }

    CHECK(OpenMagnetics::roundFloat(bobbinArea, 6) == OpenMagnetics::roundFloat(sectionsArea, 6));
    for (size_t i = 0; i < numberAssignedParallels.size() - 1; ++i){
        CHECK(round(numberAssignedParallels[i]) == round(numberParallels[i]));
        CHECK(round(numberAssignedPhysicalTurns[i]) == round(numberTurns[i] * numberParallels[i]));
    }
    CHECK(sectionsDescription.size() - numberInsulationSections == (interleavingLevel * numberTurns.size()));
    CHECK(!OpenMagnetics::check_collisions(dimensionsByName, coordinatesByName));
}

void check_layers_description(OpenMagnetics::CoilWrapper coil,
                              OpenMagnetics::WindingOrientation layersOrientation) {
    if (!coil.get_layers_description()) {
        return;
    }
    auto sections = coil.get_sections_description().value();
    std::map<std::string, std::vector<double>> dimensionsByName;
    std::map<std::string, std::vector<double>> coordinatesByName;

    for (auto& section : sections){
        auto layers = coil.get_layers_by_section(section.get_name());
        if(section.get_type() == OpenMagnetics::ElectricalType::INSULATION) {
        }
        else {
            auto sectionParallelsProportionExpected = section.get_partial_windings()[0].get_parallels_proportion();
            std::vector<double> sectionParallelsProportion(sectionParallelsProportionExpected.size(), 0);
            for (auto& layer : layers){
                for (size_t i = 0; i < sectionParallelsProportion.size(); ++i){
                    sectionParallelsProportion[i] += layer.get_partial_windings()[0].get_parallels_proportion()[i];
                }
                CHECK(layer.get_filling_factor().value() > 0);

                dimensionsByName[layer.get_name()] = layer.get_dimensions();
                coordinatesByName[layer.get_name()] = layer.get_coordinates();
            }
            for (size_t i = 0; i < sectionParallelsProportion.size(); ++i){
                CHECK(OpenMagnetics::roundFloat(sectionParallelsProportion[i], 9) == OpenMagnetics::roundFloat(sectionParallelsProportionExpected[i], 9));
            }
            for (size_t i = 0; i < layers.size() - 1; ++i){
                if (layersOrientation == OpenMagnetics::WindingOrientation::OVERLAPPING) {
                    CHECK(layers[i].get_coordinates()[0] < layers[i + 1].get_coordinates()[0]);
                    CHECK(layers[i].get_coordinates()[1] == layers[i + 1].get_coordinates()[1]);
                    CHECK(layers[i].get_coordinates()[2] == layers[i + 1].get_coordinates()[2]);
                } 
                else if (layersOrientation == OpenMagnetics::WindingOrientation::CONTIGUOUS) {
                    CHECK(layers[i].get_coordinates()[1] > layers[i + 1].get_coordinates()[1]);
                    CHECK(layers[i].get_coordinates()[0] == layers[i + 1].get_coordinates()[0]);
                    CHECK(layers[i].get_coordinates()[2] == layers[i + 1].get_coordinates()[2]);
                }
            }
        }

    }

    CHECK(!OpenMagnetics::check_collisions(dimensionsByName, coordinatesByName));
}


bool check_turns_description(OpenMagnetics::CoilWrapper coil) {
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
    if (bobbinWindingWindowShape == OpenMagnetics::WindingWindowShape::ROUND) {
        coil.convert_turns_to_cartesian_coordinates();
    }


    int turnsIn0 = 0;
    for (auto& turn : turns){
        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        if (windingIndex == 0 && turn.get_parallel() == 0) {
            turnsIn0++;
        }
        parallelProportion[windingIndex][turn.get_parallel()] += 1.0 / coil.get_number_turns(windingIndex);

        if (bobbinWindingWindowShape != OpenMagnetics::WindingWindowShape::ROUND || wires[windingIndex].get_type() != OpenMagnetics::WireType::RECTANGULAR) {
            dimensionsByName[turn.get_name()] = turn.get_dimensions().value();
        }
        if (bobbinWindingWindowShape == OpenMagnetics::WindingWindowShape::RECTANGULAR) {
            coordinatesByName[turn.get_name()] = turn.get_coordinates();
        }
        else {
            double xCoordinate = turn.get_coordinates()[0];
            double yCoordinate = turn.get_coordinates()[1];
            if (wires[windingIndex].get_type() != OpenMagnetics::WireType::RECTANGULAR) {
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
            equalToOne &= OpenMagnetics::roundFloat(parallelProportion[windingIndex][parallelIndex], 9) == 1;
        }
    }

    CHECK(equalToOne);
    bool collides = OpenMagnetics::check_collisions(dimensionsByName, coordinatesByName, bobbinWindingWindowShape == OpenMagnetics::WindingWindowShape::ROUND);
    CHECK(!collides);
    if (additionalCoordinatesByName.size() > 0) {
        collides |= OpenMagnetics::check_collisions(dimensionsByName, additionalCoordinatesByName, bobbinWindingWindowShape == OpenMagnetics::WindingWindowShape::ROUND);
    }
    CHECK(!collides);
    return !collides && equalToOne;
}

} // namespace OpenMagneticsTesting