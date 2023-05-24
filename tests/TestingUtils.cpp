#include "TestingUtils.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace OpenMagneticsTesting {


OpenMagnetics::WindingWrapper get_quick_winding(std::vector<uint64_t> numberTurns,
                                                std::vector<uint64_t> numberParallels,
                                                double bobbinHeight,
                                                double bobbinWidth,
                                                std::vector<double> bobbinCenterCoodinates,
                                                uint64_t interleavingLevel,
                                                OpenMagnetics::OrientationEnum windingOrientation,
                                                OpenMagnetics::OrientationEnum layersOrientation){
    json windingJson;
    windingJson["functionalDescription"] = json::array();

    windingJson["bobbin"] = json();
    windingJson["bobbin"]["processedDescription"] = json();
    windingJson["bobbin"]["processedDescription"]["wallThickness"] = 0.001;
    windingJson["bobbin"]["processedDescription"]["columnThickness"] = 0.001;
    windingJson["bobbin"]["processedDescription"]["columnShape"] = OpenMagnetics::ColumnShape::ROUND;
    windingJson["bobbin"]["processedDescription"]["columnDepth"] = bobbinCenterCoodinates[0] - bobbinWidth / 2;
    windingJson["bobbin"]["processedDescription"]["windingWindows"] = json::array();

    json windingWindow;
    windingWindow["height"] = bobbinHeight;
    windingWindow["width"] = bobbinWidth;
    windingWindow["coordinates"] = json::array();
    for (auto& coord : bobbinCenterCoodinates) {
        windingWindow["coordinates"].push_back(coord);
    }
    windingJson["bobbin"]["processedDescription"]["windingWindows"].push_back(windingWindow);

    for (size_t i = 0; i < numberTurns.size(); ++i){
        json individualWindingJson;
        individualWindingJson["name"] = "winding " + std::to_string(i);
        individualWindingJson["numberTurns"] = numberTurns[i];
        individualWindingJson["numberParallels"] = numberParallels[i];
        individualWindingJson["isolationSide"] = "primary";
        individualWindingJson["wire"] = "0.475 - Grade 1";
        windingJson["functionalDescription"].push_back(individualWindingJson);
    }

    OpenMagnetics::WindingWrapper winding(windingJson, interleavingLevel, windingOrientation, layersOrientation);
    return winding;
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