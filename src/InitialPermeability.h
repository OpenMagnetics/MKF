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

namespace OpenMagnetics {

    class InitialPermeability {
        private:
        protected:
        public:
            double get_initial_permeability(CoreMaterialDataOrNameUnion material,
                                            double* temperature=nullptr,
                                            double* magneticFieldDcBias=nullptr,
                                            double* frequency=nullptr);
    };

}