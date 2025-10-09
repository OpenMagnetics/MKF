#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "support/Utils.h"

#include "constructive_models/Core.h"
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include <MAS.hpp>
#include "constructive_models/Coil.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

class MagnetizingInductance {
  private:
    std::map<std::string, std::string> _models;

  protected:
  public:
    MagnetizingInductance() {
        _models["gapReluctance"] = to_string(Defaults().reluctanceModelDefault);
    }

    MagnetizingInductance(ReluctanceModels model) {
        _models["gapReluctance"] = to_string(model);
    }

    MagnetizingInductance(std::string model) {
        _models["gapReluctance"] = model;
    }
    MagnetizingInductanceOutput calculate_inductance_from_number_turns_and_gapping(Core core,
                                                                                   Coil coil,
                                                                                   OperatingPoint* operatingPoint = nullptr);

    MagnetizingInductanceOutput calculate_inductance_from_number_turns_and_gapping(Magnetic magnetic,
                                                                                   OperatingPoint* operatingPoint = nullptr);

    std::vector<CoreGap> calculate_gapping_from_number_turns_and_inductance(Core core,
                                                                      Coil coil,
                                                                      Inputs* inputs,
                                                                      GappingType gappingType,
                                                                      size_t decimals = 4);

    int calculate_number_turns_from_gapping_and_inductance(Core core, Inputs* inputs, DimensionalValues preferredValue = DimensionalValues::NOMINAL);

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> calculate_inductance_and_magnetic_flux_density(
        Core core,
        Coil coil,
        OperatingPoint* operatingPoint = nullptr);

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> calculate_inductance_and_magnetic_flux_density(
        Magnetic magnetic,
        OperatingPoint* operatingPoint = nullptr);

    double calculate_inductance_air_solenoid(
        Magnetic magnetic);

    double calculate_inductance_air_solenoid(
        Core core,
        Coil coil);
};

} // namespace OpenMagnetics