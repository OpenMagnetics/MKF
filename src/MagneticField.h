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
#include <MAS.hpp>
#include <magic_enum.hpp>
#include <InputsWrapper.h>


namespace OpenMagnetics {

    class MagneticField {
        private:
        protected:
        public:
            static ElectromagneticParameter get_magnetic_flux(ElectromagneticParameter magnetizingCurrent, double reluctance, double numberTurns, double frequency);
            static ElectromagneticParameter get_magnetic_flux_density(ElectromagneticParameter magneticFlux, double area, double frequency);
            static ElectromagneticParameter get_magnetic_field_strength(ElectromagneticParameter magneticFluxDensity, double initialPermeability, double frequency);
    };

}