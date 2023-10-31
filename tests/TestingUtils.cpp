#include "TestingUtils.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace OpenMagneticsTesting {


OpenMagnetics::CoilWrapper get_quick_coil(std::vector<uint64_t> numberTurns,
                                                std::vector<uint64_t> numberParallels,
                                                std::string shapeName,
                                                uint64_t interleavingLevel,
                                                OpenMagnetics::WindingOrientation windingOrientation,
                                                OpenMagnetics::WindingOrientation layersOrientation,
                                                OpenMagnetics::CoilAlignment turnsAlignment,
                                                OpenMagnetics::CoilAlignment sectionsAlignment,
                                                std::vector<OpenMagnetics::WireWrapper> wires,
                                                bool useBobbin){
    json coilJson;
    coilJson["functionalDescription"] = json::array();

    auto core = get_core(shapeName, json::parse("[]"), 1, "Dummy");
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

OpenMagnetics::CoilWrapper get_quick_coil(std::vector<uint64_t> numberTurns,
                                                std::vector<uint64_t> numberParallels,
                                                double bobbinHeight,
                                                double bobbinWidth,
                                                std::vector<double> bobbinCenterCoodinates,
                                                uint64_t interleavingLevel,
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

OpenMagnetics::CoilWrapper get_quick_coil_no_compact(std::vector<uint64_t> numberTurns,
                                                           std::vector<uint64_t> numberParallels,
                                                           double bobbinHeight,
                                                           double bobbinWidth,
                                                           std::vector<double> bobbinCenterCoodinates,
                                                           uint64_t interleavingLevel,
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

    OpenMagnetics::CoilWrapper coil(coilJson, interleavingLevel, windingOrientation, layersOrientation, turnsAlignment, sectionsAlignment, false);
    return coil;
}

OpenMagnetics::CoreWrapper get_core(std::string shapeName,
                                    json basicGapping,
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

    return core;
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

void print(std::vector<uint64_t> data) {
    for (auto i : data) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
}

void print(std::vector<int64_t> data) {
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

void print(std::vector<SVG::Point> data) {
    for (auto i : data) {
        std::cout << '(' << i.first << ", " << i.second << "), ";
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

} // namespace OpenMagneticsTesting