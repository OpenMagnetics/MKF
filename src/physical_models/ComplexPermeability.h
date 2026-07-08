#pragma once
#include <MAS.hpp>
#include "spline.h"

using namespace MAS;

namespace OpenMagnetics {
// ABT #113: per-thread memo caches (lazily built from the frozen material
// catalog; per-thread copies are lock-free and semantically transparent).
inline thread_local std::map<std::string, tk::spline> complexPermeabilityRealInterps;
inline thread_local std::map<std::string, tk::spline> complexPermeabilityImaginaryInterps;
// ABT #167: measured frequency span of each interpolator's data. The splines are
// interpolators, not extrapolators — evaluated past the last measured point they
// diverge polynomially (e.g. µ'' → −2.5e6 at 1 GHz for materials whose data ends
// at 1.3 MHz, an active element that made impedance sweeps rise instead of roll
// off). Queries outside the span are clamped to the nearest measured endpoint.
inline thread_local std::map<std::string, std::pair<double, double>> complexPermeabilityRealFrequencySpans;
inline thread_local std::map<std::string, std::pair<double, double>> complexPermeabilityImaginaryFrequencySpans;

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