#pragma once

#include "Constants.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "constructive_models/Wire.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

extern bool verboseTests;

using namespace MAS;
using namespace OpenMagnetics;

namespace OpenMagneticsTesting {
OpenMagnetics::Coil get_quick_coil(std::vector<int64_t> numberTurns,
                                          std::vector<int64_t> numberParallels,
                                          std::string shapeName,
                                          uint8_t interleavingLevel = 1,
                                          MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING,
                                          MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING,
                                          MAS::CoilAlignment turnsAlignment = MAS::CoilAlignment::CENTERED,
                                          MAS::CoilAlignment sectionsAlignment = MAS::CoilAlignment::CENTERED,
                                          std::vector<OpenMagnetics::Wire> wires = std::vector<OpenMagnetics::Wire>({}),
                                          bool useBobbin = true,
                                          int numberStacks = 1);

OpenMagnetics::Coil get_quick_coil(std::vector<int64_t> numberTurns,
                                          std::vector<int64_t> numberParallels,
                                          double bobbinHeight,
                                          double bobbinWidth,
                                          std::vector<double> bobbinCenterCoodinates,
                                          uint8_t interleavingLevel = 1,
                                          MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING,
                                          MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING,
                                          MAS::CoilAlignment turnsAlignment = MAS::CoilAlignment::CENTERED,
                                          MAS::CoilAlignment sectionsAlignment = MAS::CoilAlignment::CENTERED,
                                          std::vector<OpenMagnetics::Wire> wires = std::vector<OpenMagnetics::Wire>({}));

OpenMagnetics::Coil get_quick_coil_no_compact(std::vector<int64_t> numberTurns,
                                                     std::vector<int64_t> numberParallels,
                                                     double bobbinHeight,
                                                     double bobbinWidth,
                                                     std::vector<double> bobbinCenterCoodinates,
                                                     uint8_t interleavingLevel = 1,
                                                     MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING,
                                                     MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING,
                                                     MAS::CoilAlignment turnsAlignment = MAS::CoilAlignment::CENTERED,
                                                     MAS::CoilAlignment sectionsAlignment = MAS::CoilAlignment::CENTERED,
                                                     std::vector<OpenMagnetics::Wire> wires = std::vector<OpenMagnetics::Wire>({}));

OpenMagnetics::Coil get_quick_toroidal_coil_no_compact(std::vector<int64_t> numberTurns,
                                                              std::vector<int64_t> numberParallels,
                                                              double bobbinRadialHeight,
                                                              double bobbinAngle,
                                                              double columnDepth,
                                                              uint8_t interleavingLevel = 1,
                                                              MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING,
                                                              MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING,
                                                              MAS::CoilAlignment turnsAlignment = MAS::CoilAlignment::CENTERED,
                                                              MAS::CoilAlignment sectionsAlignment = MAS::CoilAlignment::CENTERED,
                                                              std::vector<OpenMagnetics::Wire> wires = std::vector<OpenMagnetics::Wire>({}));


OpenMagnetics::Inputs get_quick_insulation_inputs(MAS::DimensionWithTolerance altitude,
                                                         MAS::Cti cti,
                                                         MAS::InsulationType insulation_type,
                                                         MAS::DimensionWithTolerance main_supply_voltage,
                                                         MAS::OvervoltageCategory overvoltage_category,
                                                         MAS::PollutionDegree pollution_degree,
                                                         std::vector<MAS::InsulationStandards> standards,
                                                         double maximumVoltageRms,
                                                         double maximumVoltagePeak,
                                                         double frequency,
                                                         MAS::WiringTechnology wiringTechnology = MAS::WiringTechnology::WOUND);

MAS::InsulationRequirements get_quick_insulation_requirements(MAS::DimensionWithTolerance altitude,
                                                                        MAS::Cti cti,
                                                                        MAS::InsulationType insulation_type,
                                                                        MAS::DimensionWithTolerance main_supply_voltage,
                                                                        MAS::OvervoltageCategory overvoltage_category,
                                                                        MAS::PollutionDegree pollution_degree,
                                                                        std::vector<MAS::InsulationStandards> standards);

Core get_quick_core(std::string shapeName,
                                          json basicGapping,
                                          int numberStacks = 1,
                                          std::string materialName = "N87");

OpenMagnetics::Magnetic get_quick_magnetic(std::string shapeName,
                                                  json basicGapping,
                                                  std::vector<int64_t> numberTurns,
                                                  int numberStacks = 1,
                                                  std::string materialName = "N87");
json get_ground_gap(double gapLength);
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


void check_sections_description(OpenMagnetics::Coil coil,
                                std::vector<int64_t> numberTurns,
                                std::vector<int64_t> numberParallels,
                                uint8_t interleavingLevel = 1,
                                MAS::WindingOrientation windingOrientation = MAS::WindingOrientation::OVERLAPPING);

void check_layers_description(OpenMagnetics::Coil coil,
                                      MAS::WindingOrientation layersOrientation = MAS::WindingOrientation::OVERLAPPING);


bool check_turns_description(OpenMagnetics::Coil coil);
bool check_wire_standards(OpenMagnetics::Coil coil);
void check_winding_losses(OpenMagnetics::Mas mas);

OpenMagnetics::Mas mas_loader(std::string path);

} // namespace OpenMagneticsTesting