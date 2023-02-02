#pragma once
#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <filesystem>
#include <streambuf>
#include "Constants.h"
#include <magic_enum.hpp>

#include <CoreWrapper.h>
#include <WindingWrapper.h>
#include <InputsWrapper.h>
#include <MAS.hpp>

namespace OpenMagnetics {

    class MagnetizingInductance {
        private:
            std::map<std::string, std::string> _models;
        protected:
        public:
            MagnetizingInductance(std::map<std::string, std::string> models) {
                _models = models;
            }
            double get_inductance_from_number_turns_and_gapping(CoreWrapper core,
                                                                WindingWrapper winding,
                                                                OperationPoint operationPoint);

            std::vector<CoreGap> get_gapping_from_number_turns_and_inductance(CoreWrapper core,
                                                                              WindingWrapper winding,
                                                                              InputsWrapper inputs,
                                                                              GappingType gappingType, 
                                                                              size_t decimals = 4);

            double get_number_turns_from_gapping_and_inductance(CoreWrapper core,
                                                                InputsWrapper inputs);
    };

}