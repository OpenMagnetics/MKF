#pragma once
#include <MAS.hpp>
#include "spline.h"
#include "physical_models/ComplexPermeabilityEstimator.h"
#include <memory>

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
// off). Queries outside the span throw ComplexPermeabilityFrequencyOutOfRangeException;
// holding the endpoint instead hid stale data (K081 tabulated 1 Hz..1 MHz read its
// 1 MHz point at 350 MHz).
inline thread_local std::map<std::string, std::pair<double, double>> complexPermeabilityRealFrequencySpans;
inline thread_local std::map<std::string, std::pair<double, double>> complexPermeabilityImaginaryFrequencySpans;
// ABT #1861: the estimate of each material's untabulated range (ComplexPermeabilityEstimator), built
// the first time a frequency outside a part's table is asked for; and the reason, for the materials
// the estimator refuses, so the refusal is not recomputed per call.
inline thread_local std::map<std::string, std::shared_ptr<const ComplexPermeabilityEstimate>> complexPermeabilityEstimates;
inline thread_local std::map<std::string, std::string> complexPermeabilityEstimateRefusals;
// The catalogue generation (Utils.h) the maps above were filled under: the estimates read donor
// materials from the catalogue, so a cleared or reloaded catalogue empties them (ABT #710).
inline thread_local uint64_t complexPermeabilityCachesCatalogueGeneration = 0;

// A complex permeability and, per part, where it came from (TABULATED inside the part's table).
// realMethod / imaginaryMethod describe an estimated part (empty for a tabulated one).
struct ComplexPermeabilityValue {
    double real;
    double imaginary;
    ComplexPermeabilitySource realSource;
    ComplexPermeabilitySource imaginarySource;
    std::string realMethod;
    std::string imaginaryMethod;
    bool is_estimated() const {
        return realSource != ComplexPermeabilitySource::TABULATED || imaginarySource != ComplexPermeabilitySource::TABULATED;
    }
};

class ComplexPermeability {
    private:
        // Mueller et al., "Novel Complex Permeability Model of Powder Magnetic
        // Materials" (PCIM 2024), Fig 1 — maps ΔF_{L,95-90} (eq 13) to F_µ
        // (the base-material / pressed-core permeability ratio). Inverse
        // lookup, monotonically increasing, clamped at the edges (F_µ → 1
        // for ΔF < 0.31, F_µ → 1000 for ΔF > 0.755). Linear in (ΔF, log F_µ).
        static double infer_F_mu_from_delta_FL(double deltaFL_95_90);
        // Builds (once per material content) the real and imaginary interpolators and their spans,
        // and returns the key they are cached under (cache_key).
        std::string ensure_interpolators(const CoreMaterial& coreMaterial);
        // ABT #710: the caches are keyed by the material's name AND its permeability data. Keyed by
        // the name alone, an inline material that shared a name with one evaluated before (a database
        // record, or an edited copy) silently got the earlier material's curves.
        static std::string cache_key(const CoreMaterial& coreMaterial);
        // The tabulated range of the interpolators cached under cacheKey (both parts).
        static std::pair<double, double> frequency_range_for_key(const std::string& cacheKey, const std::string& materialName);
        // The cached estimate of the material's untabulated range; throws
        // ComplexPermeabilityEstimateUnavailableException (naming frequency) when there is none.
        const ComplexPermeabilityEstimate& estimate_for(const CoreMaterial& coreMaterial, const std::string& cacheKey, double frequency);
        // The cached estimate, built on first use; nullptr when the estimator refuses the material
        // (the reason is kept in complexPermeabilityEstimateRefusals).
        std::shared_ptr<const ComplexPermeabilityEstimate> find_or_build_estimate(const CoreMaterial& coreMaterial, const std::string& cacheKey);
    protected:
    public:
        // mu' and mu'' at a frequency, ALWAYS answered for a material the estimator covers (ABT #1861,
        // Alf 2026-10-14): inside each part's table, the table; outside it, the approved estimate
        // (ComplexPermeabilityEstimator: class neighbours, their tail, mu' held flat below its data),
        // regardless of Settings::allowMaterialDataExtrapolation, which governs the other material
        // data only. Throws ComplexPermeabilityEstimateUnavailableException outside the tabulated
        // range for the materials the estimator refuses (powders, MgZn, unclassified, derived tables).
        // For graphs and analysis of a chosen design: get_complex_permeability_with_source says
        // which parts are estimated. Advisers and filters use get_tabulated_complex_permeability.
        std::pair<double, double> get_complex_permeability(std::string coreMaterialName, double frequency);
        std::pair<double, double> get_complex_permeability(CoreMaterial coreMaterial, double frequency);
        ComplexPermeabilityValue get_complex_permeability_with_source(const CoreMaterial& coreMaterial, double frequency);
        // Inside the tabulated range only (get_frequency_range): throws
        // ComplexPermeabilityFrequencyOutOfRangeException outside it. The advisers and filters rank and
        // select parts on data, never on estimates (ABT #1861, proposal D5).
        std::pair<double, double> get_tabulated_complex_permeability(const CoreMaterial& coreMaterial, double frequency);
        // The stretches of frequency the material's complex permeability is estimated in (empty when
        // the estimator refuses the material: then nothing outside the table is answered at all).
        std::vector<ComplexPermeabilityEstimatedRange> get_estimated_ranges(const CoreMaterial& coreMaterial);
        // [minimum, maximum] frequency, in Hz, over which both mu' and mu'' are TABULATED. Outside it
        // get_complex_permeability estimates and get_tabulated_complex_permeability throws; a caller
        // that wants to stay on data reads this range and limits itself explicitly.
        std::pair<double, double> get_frequency_range(CoreMaterial coreMaterial);
        std::pair<double, double> get_frequency_range(std::string coreMaterialName);
        ComplexPermeabilityData calculate_complex_permeability_from_frequency_dependent_initial_permeability(CoreMaterial coreMaterial);
        ComplexPermeabilityData calculate_complex_permeability_from_frequency_dependent_initial_permeability(std::string coreMaterialName);
};

} // namespace OpenMagnetics
