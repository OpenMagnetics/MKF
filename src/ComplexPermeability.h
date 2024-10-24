#pragma once
#include "Constants.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class ComplexPermeability {
    private:
    protected:
    public:
        std::pair<double, double> get_complex_permeability(std::string coreMaterialName, double frequency);

        std::pair<double, double> get_complex_permeability(CoreMaterial coreMaterial, double frequency);
};

} // namespace OpenMagnetics