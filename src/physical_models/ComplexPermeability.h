#pragma once
#include <MAS.hpp>
#include "spline.h"

using namespace MAS;

namespace OpenMagnetics {
inline std::map<std::string, tk::spline> complexPermeabilityRealInterps;
inline std::map<std::string, tk::spline> complexPermeabilityImaginaryInterps;

class ComplexPermeability {
    private:
        // Mueller et al., "Novel Complex Permeability Model of Powder Magnetic
        // Materials" (PCIM 2024), Fig 1 — maps ΔF_{L,95-90} (eq 13) to F_µ
        // (the base-material / pressed-core permeability ratio). Inverse
        // lookup, monotonically increasing, clamped at the edges (F_µ → 1
        // for ΔF < 0.31, F_µ → 1000 for ΔF > 0.755). Linear in (ΔF, log F_µ).
        static double infer_F_mu_from_delta_FL(double deltaFL_95_90);
    protected:
    public:
        std::pair<double, double> get_complex_permeability(std::string coreMaterialName, double frequency);
        std::pair<double, double> get_complex_permeability(CoreMaterial coreMaterial, double frequency);
        ComplexPermeabilityData calculate_complex_permeability_from_frequency_dependent_initial_permeability(CoreMaterial coreMaterial);
        ComplexPermeabilityData calculate_complex_permeability_from_frequency_dependent_initial_permeability(std::string coreMaterialName);
};

} // namespace OpenMagnetics