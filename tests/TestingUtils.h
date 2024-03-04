#pragma once

#include "Constants.h"
#include "CoreWrapper.h"
#include "CoilWrapper.h"
#include "MagneticWrapper.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

extern bool verboseTests;

namespace OpenMagneticsTesting {
OpenMagnetics::CoilWrapper get_quick_coil(std::vector<int64_t> numberTurns,
                                          std::vector<int64_t> numberParallels,
                                          std::string shapeName,
                                          uint8_t interleavingLevel,
                                          OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING,
                                          OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING,
                                          OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                          OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                          std::vector<OpenMagnetics::WireWrapper> wires = std::vector<OpenMagnetics::WireWrapper>({}),
                                          bool useBobbin = true);

OpenMagnetics::CoilWrapper get_quick_coil(std::vector<int64_t> numberTurns,
                                          std::vector<int64_t> numberParallels,
                                          double bobbinHeight,
                                          double bobbinWidth,
                                          std::vector<double> bobbinCenterCoodinates,
                                          uint8_t interleavingLevel,
                                          OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING,
                                          OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING,
                                          OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                          OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                          std::vector<OpenMagnetics::WireWrapper> wires = std::vector<OpenMagnetics::WireWrapper>({}));

OpenMagnetics::CoilWrapper get_quick_coil_no_compact(std::vector<int64_t> numberTurns,
                                                     std::vector<int64_t> numberParallels,
                                                     double bobbinHeight,
                                                     double bobbinWidth,
                                                     std::vector<double> bobbinCenterCoodinates,
                                                     uint8_t interleavingLevel,
                                                     OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING,
                                                     OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING,
                                                     OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                                     OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED,
                                                     std::vector<OpenMagnetics::WireWrapper> wires = std::vector<OpenMagnetics::WireWrapper>({}));

OpenMagnetics::CoilWrapper get_quick_toroidal_coil_no_compact(std::vector<int64_t> numberTurns,
                                                              std::vector<int64_t> numberParallels,
                                                              double bobbinRadialHeight,
                                                              double bobbinAngle,
                                                              double columnDepth,
                                                              uint8_t interleavingLevel,
                                                              OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING,
                                                              OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING,
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

OpenMagnetics::InsulationRequirements get_quick_insulation_requirements(OpenMagnetics::DimensionWithTolerance altitude,
                                                                        OpenMagnetics::Cti cti,
                                                                        OpenMagnetics::InsulationType insulation_type,
                                                                        OpenMagnetics::DimensionWithTolerance main_supply_voltage,
                                                                        OpenMagnetics::OvervoltageCategory overvoltage_category,
                                                                        OpenMagnetics::PollutionDegree pollution_degree,
                                                                        std::vector<OpenMagnetics::InsulationStandards> standards);

OpenMagnetics::CoreWrapper get_quick_core(std::string shapeName,
                                          json basicGapping,
                                          int numberStacks = 1,
                                          std::string materialName = "N87");

OpenMagnetics::MagneticWrapper get_quick_magnetic(std::string shapeName,
                                                  json basicGapping,
                                                  std::vector<int64_t> numberTurns,
                                                  int numberStacks = 1,
                                                  std::string materialName = "N87");
json get_grinded_gap(double gapLength);
json get_distributed_gap(double gapLength, int numberGaps);
json get_spacer_gap(double gapLength);
json get_residual_gap();

void print(std::vector<double> data);
void print(std::vector<std::vector<double>> data);
void print(std::vector<int64_t> data);
void print(std::vector<uint64_t> data);
void print(std::vector<std::string> data);
void print(double data);
void print(std::string data);
void print_json(json data);


void check_sections_description(OpenMagnetics::CoilWrapper coil,
                                std::vector<int64_t> numberTurns,
                                std::vector<int64_t> numberParallels,
                                uint8_t interleavingLevel,
                                OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING);

void check_layers_description(OpenMagnetics::CoilWrapper coil,
                                      OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING);


bool check_turns_description(OpenMagnetics::CoilWrapper coil);

} // namespace OpenMagneticsTesting