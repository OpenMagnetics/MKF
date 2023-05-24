#pragma once

#include "Constants.h"
#include "CoreWrapper.h"
#include "WindingWrapper.h"
#include "svg.hpp"

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
                                                OpenMagnetics::OrientationEnum windingOrientation = OpenMagnetics::OrientationEnum::HORIZONTAL,
                                                OpenMagnetics::OrientationEnum layersOrientation = OpenMagnetics::OrientationEnum::VERTICAL);

OpenMagnetics::CoreWrapper get_core(std::string shapeName,
                                    json basicGapping,
                                    int numberStacks = 1,
                                    std::string materialName = "N87");
json get_grinded_gap(double gapLength);
json get_distributed_gap(double gapLength, int numberGaps);
json get_spacer_gap(double gapLength);
json get_residual_gap();

void print(std::vector<double> data);
void print(std::vector<uint64_t> data);
void print(std::vector<std::string> data);
void print(std::vector<SVG::Point> data);
void print(double data);
void print(std::string data);
void print_json(json data);

} // namespace OpenMagneticsTesting