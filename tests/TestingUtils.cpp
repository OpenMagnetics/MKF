#include "TestingUtils.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace OpenMagneticsTesting {

OpenMagnetics::CoreWrapper get_core(std::string shapeName,
                                    json basicGapping,
                                    int numberStacks,
                                    std::string materialName) {
    auto coreJson = json();

    coreJson["functionalDescription"] = json();
    coreJson["functionalDescription"]["name"] = "GapReluctanceTest";
    coreJson["functionalDescription"]["type"] = "two-piece set";
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

} // namespace OpenMagneticsTesting