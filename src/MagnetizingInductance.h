#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "Utils.h"

#include <CoreWrapper.h>
#include <InputsWrapper.h>
#include <MAS.hpp>
#include <CoilWrapper.h>
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
                                                                                   CoilWrapper winding,
                                                                                   OperatingPoint* operatingPoint);

    std::vector<CoreGap> calculate_gapping_from_number_turns_and_inductance(CoreWrapper core,
                                                                      CoilWrapper winding,
                                                                      InputsWrapper* inputs,
                                                                      GappingType gappingType,
                                                                      size_t decimals = 4);

    int calculate_number_turns_from_gapping_and_inductance(CoreWrapper core, InputsWrapper* inputs, DimensionalValues preferredValue = DimensionalValues::NOMINAL);

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> calculate_inductance_and_magnetic_flux_density(
        CoreWrapper core,
        CoilWrapper winding,
        OperatingPoint* operatingPoint);
};

} // namespace OpenMagnetics