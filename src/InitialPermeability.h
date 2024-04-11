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
        double get_initial_permeability(std::string coreMaterialName,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt);

        double get_initial_permeability(CoreMaterial coreMaterial,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt);
};

} // namespace OpenMagnetics