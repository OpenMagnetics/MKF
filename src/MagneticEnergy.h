#pragma once
#include "Constants.h"
#include "Defaults.h"

#include <Reluctance.h>
#include <CoreWrapper.h>
#include <InputsWrapper.h>
#include <MAS.hpp>
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

class MagneticEnergy {
  private:
    std::map<std::string, std::string> _models;
  protected:
  public:
    MagneticEnergy(std::map<std::string, std::string> models) {
        auto defaults = OpenMagnetics::Defaults();
        _models = models;
        if (models.find("gapReluctance") == models.end()) {
            _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
        }
    }
    double get_ungapped_core_maximum_magnetic_energy(CoreWrapper core, OperatingPoint* operatingPoint = nullptr);
    double get_gap_maximum_magnetic_energy(CoreGap gapInfo, double magneticFluxDensitySaturation, double* frequency = nullptr);
    double get_core_maximum_magnetic_energy(CoreWrapper core, OperatingPoint* operatingPoint = nullptr);
    DimensionWithTolerance required_magnetic_energy(InputsWrapper inputs);
};

} // namespace OpenMagnetics