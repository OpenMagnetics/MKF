#pragma once
#include "Constants.h"
#include "Defaults.h"

#include "physical_models/Reluctance.h"
#include "constructive_models/CoreWrapper.h"
#include "processors/InputsWrapper.h"
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
    MagneticEnergy() {
        auto defaults = OpenMagnetics::Defaults();
        _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
    }
    static double get_ungapped_core_maximum_magnetic_energy(CoreWrapper core, std::optional<OperatingPoint> operatingPoint = std::nullopt);
    static double get_ungapped_core_maximum_magnetic_energy(CoreWrapper core, double temperature, std::optional<double> frequency = std::nullopt);
    double get_gap_maximum_magnetic_energy(CoreGap gapInfo, double magneticFluxDensitySaturation, std::optional<double> fringing_factor = std::nullopt);
    double calculate_core_maximum_magnetic_energy(CoreWrapper core, std::optional<OperatingPoint> operatingPoint = std::nullopt);
    double calculate_core_maximum_magnetic_energy(CoreWrapper core, double temperature, std::optional<double> frequency = std::nullopt);
    DimensionWithTolerance calculate_required_magnetic_energy(InputsWrapper inputs);
};

} // namespace OpenMagnetics