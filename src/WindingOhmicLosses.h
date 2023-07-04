#pragma once
#include "MAS.hpp"
#include "Resistivity.h"
#include "WindingWrapper.h"
#include "WireWrapper.h"
#include "Utils.h"

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

class WindingOhmicLosses {
    private:
        std::shared_ptr<ResistivityModel> _resistivityModel;
    protected:
    public:
        WindingOhmicLosses() {
            _resistivityModel = OpenMagnetics::ResistivityModel::factory(OpenMagnetics::ResistivityModels::WIRE_MATERIAL);
        }
        double get_dc_resistance(Turn turn, WireWrapper wire, double temperature);
        double get_ohmic_losses(WindingWrapper winding, OperationPoint operationPoint, double temperature);
};

} // namespace OpenMagnetics