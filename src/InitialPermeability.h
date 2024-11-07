#pragma once
#include "Constants.h"

#include <CoreWrapper.h>
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

class InitialPermeability {
    private:
    protected:
    public:
        static double get_initial_permeability(std::string coreMaterialName,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt,
                                        std::optional<double> magneticFluxDensity = std::nullopt);

        static double get_initial_permeability(CoreMaterial coreMaterial,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt,
                                        std::optional<double> magneticFluxDensity = std::nullopt);
};

} // namespace OpenMagnetics