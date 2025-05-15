#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "support/Utils.h"

#include "constructive_models/CoreWrapper.h"
#include "processors/InputsWrapper.h"
#include "constructive_models/MagneticWrapper.h"
#include <MAS.hpp>
#include "constructive_models/CoilWrapper.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

class MagnetizingInductance {
  private:
    std::map<std::string, std::string> _models;

  protected:
  public:
    MagnetizingInductance() {
        _models["gapReluctance"] = magic_enum::enum_name(Defaults().reluctanceModelDefault);
    }

    MagnetizingInductance(ReluctanceModels model) {
        _models["gapReluctance"] = magic_enum::enum_name(model);
    }

    MagnetizingInductance(std::string model) {
        _models["gapReluctance"] = model;
    }
    MagnetizingInductanceOutput calculate_inductance_from_number_turns_and_gapping(CoreWrapper core,
                                                                                   CoilWrapper coil,
                                                                                   OperatingPoint* operatingPoint = nullptr);

    MagnetizingInductanceOutput calculate_inductance_from_number_turns_and_gapping(MagneticWrapper magnetic,
                                                                                   OperatingPoint* operatingPoint = nullptr);

    std::vector<CoreGap> calculate_gapping_from_number_turns_and_inductance(CoreWrapper core,
                                                                      CoilWrapper coil,
                                                                      InputsWrapper* inputs,
                                                                      GappingType gappingType,
                                                                      size_t decimals = 4);

    int calculate_number_turns_from_gapping_and_inductance(CoreWrapper core, InputsWrapper* inputs, DimensionalValues preferredValue = DimensionalValues::NOMINAL);

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> calculate_inductance_and_magnetic_flux_density(
        CoreWrapper core,
        CoilWrapper coil,
        OperatingPoint* operatingPoint = nullptr);

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> calculate_inductance_and_magnetic_flux_density(
        MagneticWrapper magnetic,
        OperatingPoint* operatingPoint = nullptr);

    double calculate_inductance_air_solenoid(
        MagneticWrapper magnetic);

    double calculate_inductance_air_solenoid(
        CoreWrapper core,
        CoilWrapper coil);
};

} // namespace OpenMagnetics