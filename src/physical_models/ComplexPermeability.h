#pragma once
#include <MAS.hpp>
#include "spline.h"

using namespace MAS;

namespace OpenMagnetics {
inline std::map<std::string, tk::spline> complexPermeabilityRealInterps;
inline std::map<std::string, tk::spline> complexPermeabilityImaginaryInterps;

class ComplexPermeability {
    private:
    protected:
    public:
        std::pair<double, double> get_complex_permeability(std::string coreMaterialName, double frequency);
        std::pair<double, double> get_complex_permeability(CoreMaterial coreMaterial, double frequency);
        ComplexPermeabilityData calculate_complex_permeability_from_frequency_dependent_initial_permeability(CoreMaterial coreMaterial);
        ComplexPermeabilityData calculate_complex_permeability_from_frequency_dependent_initial_permeability(std::string coreMaterialName);
};

} // namespace OpenMagnetics