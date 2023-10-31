#pragma once

#include "Constants.h"
#include "CoreWrapper.h"
#include "CoilWrapper.h"
#include "svg.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace OpenMagneticsTesting {
OpenMagnetics::CoilWrapper get_quick_coil(std::vector<uint64_t> numberTurns,
                                          std::vector<uint64_t> numberParallels,
                                          std::string shapeName,
                                          uint64_t interleavingLevel,
                                          OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL,
                                          OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL,
                                          OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                          OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                          std::vector<OpenMagnetics::WireWrapper> wires = std::vector<OpenMagnetics::WireWrapper>({}),
                                          bool useBobbin = true);

OpenMagnetics::CoilWrapper get_quick_coil(std::vector<uint64_t> numberTurns,
                                          std::vector<uint64_t> numberParallels,
                                          double bobbinHeight,
                                          double bobbinWidth,
                                          std::vector<double> bobbinCenterCoodinates,
                                          uint64_t interleavingLevel,
                                          OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL,
                                          OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL,
                                          OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                          OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                          std::vector<OpenMagnetics::WireWrapper> wires = std::vector<OpenMagnetics::WireWrapper>({}));

OpenMagnetics::CoilWrapper get_quick_coil_no_compact(std::vector<uint64_t> numberTurns,
                                                     std::vector<uint64_t> numberParallels,
                                                     double bobbinHeight,
                                                     double bobbinWidth,
                                                     std::vector<double> bobbinCenterCoodinates,
                                                     uint64_t interleavingLevel,
                                                     OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL,
                                                     OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL,
                                                     OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                                     OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                                     std::vector<OpenMagnetics::WireWrapper> wires = std::vector<OpenMagnetics::WireWrapper>({}));


OpenMagnetics::InputsWrapper get_quick_insulation_inputs(OpenMagnetics::DimensionWithTolerance altitude,
                                                         OpenMagnetics::Cti cti,
                                                         OpenMagnetics::InsulationType insulation_type,
                                                         OpenMagnetics::DimensionWithTolerance main_supply_voltage,
                                                         OpenMagnetics::OvervoltageCategory overvoltage_category,
                                                         OpenMagnetics::PollutionDegree pollution_degree,
                                                         std::vector<OpenMagnetics::InsulationStandards> standards,
                                                         double maximumVoltageRms,
                                                         double maximumVoltagePeak,
                                                         double frequency,
                                                         OpenMagnetics::WiringTechnology wiringTechnology = OpenMagnetics::WiringTechnology::WOUND);

OpenMagnetics::CoreWrapper get_core(std::string shapeName,
                                    json basicGapping,
                                    int numberStacks = 1,
                                    std::string materialName = "N87");
json get_grinded_gap(double gapLength);
json get_distributed_gap(double gapLength, int numberGaps);
json get_spacer_gap(double gapLength);
json get_residual_gap();

void print(std::vector<double> data);
void print(std::vector<int64_t> data);
void print(std::vector<uint64_t> data);
void print(std::vector<std::string> data);
void print(std::vector<SVG::Point> data);
void print(double data);
void print(std::string data);
void print_json(json data);

} // namespace OpenMagneticsTesting