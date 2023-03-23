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
#include "Defaults.h"
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
                auto defaults = OpenMagnetics::Defaults();
                _models = models;
                if (models.find("gapReluctance") == models.end()) {
                    _models["gapReluctance"] = magic_enum::enum_name(defaults.reluctanceModelDefault);
                }
            }
            double get_inductance_from_number_turns_and_gapping(CoreWrapper core,
                                                                WindingWrapper winding,
                                                                OperationPoint* operationPoint);

            std::vector<CoreGap> get_gapping_from_number_turns_and_inductance(CoreWrapper core,
                                                                              WindingWrapper winding,
                                                                              InputsWrapper* inputs,
                                                                              GappingType gappingType, 
                                                                              size_t decimals = 4);

            int get_number_turns_from_gapping_and_inductance(CoreWrapper core,
                                                             InputsWrapper* inputs);

            std::pair<double, ElectromagneticParameter> get_inductance_and_magnetic_flux_density(CoreWrapper core,
                                                                                                 WindingWrapper winding,
                                                                                                 OperationPoint* operationPoint);
    };

}