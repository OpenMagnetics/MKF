#pragma once
#include <memory>
#include "Constants.h"

#include "constructive_models/Core.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <functional>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>
#include "spline.h"

using namespace MAS;

namespace OpenMagnetics {
// ABT #113: per-thread memo caches (lazily built from the frozen material
// catalog; per-thread copies are lock-free and semantically transparent).
inline thread_local std::map<std::string, std::variant<double, tk::spline>> initialPermeabilityMagneticFieldDcBiasInterps;
inline thread_local std::map<std::string, std::variant<double, tk::spline>> initialPermeabilityFrequencyInterps;
inline thread_local std::map<std::string, std::function<double(double)>> initialPermeabilityTemperatureInterps;

class InitialPermeability {

    private:
    protected:
    public:
        static double get_initial_permeability(std::string coreMaterialName,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt,
                                        std::optional<double> magneticFluxDensity = std::nullopt,
                                        std::optional<CoreShapeFamily> shapeFamily = std::nullopt);

        static double get_initial_permeability(const CoreMaterial& coreMaterial,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt,
                                        std::optional<double> magneticFluxDensity = std::nullopt,
                                        std::optional<CoreShapeFamily> shapeFamily = std::nullopt);
        // The material model's value as it stands, which past the saturation knee may fall BELOW 1: the
        // amplitude-permeability (Roshen loop) slope and the vendors' 1/(a + b H^c) fits tend to 0 there,
        // while a real material's reversible permeability never drops under the vacuum's (B = mu0 (H + M)).
        // get_initial_permeability gives 1 for it under a bias (the saturated material, ABT #1855) and
        // refuses it without one (the model outside its frequency or temperature range, ABT #1740); only
        // callers that model the vacuum asymptote themselves read it here -- the magnetisation-curve
        // builder, which carries B on the vacuum slope past that point, and the DC-bias fixed point.
        static double get_initial_permeability_model_value(const CoreMaterial& coreMaterial,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt,
                                        std::optional<double> magneticFluxDensity = std::nullopt,
                                        std::optional<CoreShapeFamily> shapeFamily = std::nullopt);
        static double get_initial_permeability(const CoreMaterial& coreMaterial, OperatingPoint operatingPoint,
                                        std::optional<CoreShapeFamily> shapeFamily = std::nullopt);
        static double get_initial_permeability(std::string coreMaterialName, OperatingPoint operatingPoint,
                                        std::optional<CoreShapeFamily> shapeFamily = std::nullopt);

        // ABT #358: vendors publish per-shape-family DC-bias fits under modifier keys like
        // "E", "EQ", "PQ", "E/ER/U", "EQ/LP" (a slash lists the families sharing one fit).
        // Those entries are PARTIAL — typically only the DC-bias factor — so the resolved
        // modifier is "default" with every factor the family entry provides overlaid on top.
        // Matching is exact per slash-separated token, never substring (a substring test would
        // make "E" hit "EQ/LP", the ABT #359 class of bug).
        static InitialPermeabilitModifier resolve_modifier(const PermeabilityPoint& permeabilityPoint,
                                        std::optional<CoreShapeFamily> shapeFamily = std::nullopt);

        static double has_temperature_dependency(const CoreMaterial& coreMaterial);
        static double has_frequency_dependency(const CoreMaterial& coreMaterial);
        static double has_magnetic_field_dc_bias_dependency(const CoreMaterial& coreMaterial);
        static std::optional<double> get_paramagnetic_temperature(const CoreMaterial& coreMaterial);
        // The highest temperature at which the material is ferromagnetic AND its mu_i(T) model has data:
        // the last point of the measured mu_i(T) table when the table ends below the Curie temperature
        // (between the two, mu_i(T) refuses, ABT #1739), otherwise the Curie temperature. Without a
        // Curie temperature: the table's last point; with neither: none. The loss-temperature loop
        // keeps the core below it (ABT #838).
        static std::optional<double> get_maximum_ferromagnetic_data_temperature(const CoreMaterial& coreMaterial);
        // The last temperature of the material's measured mu_i(T) table; none without such a table.
        static std::optional<double> get_temperature_table_end(const CoreMaterial& coreMaterial);
        static double get_initial_permeability_temperature_dependent(const CoreMaterial& coreMaterial, double temperature);
        static double get_initial_permeability_frequency_dependent(const CoreMaterial& coreMaterial, double frequency);
        static double get_initial_permeability_magnetic_field_dc_bias_dependent(const CoreMaterial& coreMaterial, double magneticFieldDcBias);
        // ABT #1093: the field strength H_dc that carries a given DC flux density in the material.
        // The datasheet bias curve is the REVERSIBLE (incremental) permeability, dB/dH along the
        // magnetisation curve, so B(H) = µ0·∫0^H µ_rev(h)·dh must be integrated and inverted; the
        // secant B = µ0·µ_rev(H)·H is wrong past the knee, and the fixed point µ ← µ_rev(B/(µ0·µ))
        // diverges there. With refuseAboveSaturation (gap sizing) it throws when biasFluxDensity
        // exceeds the material saturation at that temperature: no gap holds that flux. Without it
        // (analysing a given design) the march continues past the knee on the vacuum slope.
        static double get_magnetic_field_dc_bias_for_flux_density(const CoreMaterial& coreMaterial,
                                        double biasFluxDensity,
                                        double temperature,
                                        std::optional<double> frequency = std::nullopt,
                                        bool refuseAboveSaturation = true);
        // The inverse of get_magnetic_field_dc_bias_for_flux_density (ABT #1855): the DC flux density B(H) the
        // material carries at the field strength H, read on its magnetisation curve (on the vacuum slope past its
        // end). Throws for a negative H.
        static double get_flux_density_for_magnetic_field_dc_bias(const CoreMaterial& coreMaterial,
                                        double magneticFieldDcBias,
                                        double temperature,
                                        std::optional<double> frequency = std::nullopt);
        // The material's magnetisation curve B(H) = mu0 * integral of mu_rev, tabulated and cached per material,
        // temperature and frequency (ABT #1223). The cache keeps the most recently used curves only (ABT #838: the
        // loss-temperature loop asks for a new temperature on every pass, and an unbounded cache grew without end);
        // the caller shares ownership of the curve it gets, so an eviction never pulls it from under the caller.
        struct MagnetisationCurve {
            std::vector<double> fieldStrength;
            std::vector<double> fluxDensity;
        };
        static constexpr size_t magnetisationCurveCacheCapacity = 256;
        static std::shared_ptr<const MagnetisationCurve> get_magnetisation_curve(const CoreMaterial& coreMaterial, double temperature, std::optional<double> frequency);
        // Number of curves the cache holds now (bounded by magnetisationCurveCacheCapacity).
        static size_t get_magnetisation_curve_cache_size();
        // The DC operating point {H, B} of the material when the MMF N*I is imposed (current-driven, or a DC
        // set by the load currents): solves N*I = H*l_e + B(H)*A_e*R_gap on the magnetisation curve.
        static std::pair<double, double> get_dc_operating_point_for_magnetomotive_force(const CoreMaterial& coreMaterial,
                                        double magnetomotiveForce,
                                        double effectiveLength,
                                        double effectiveArea,
                                        double gapReluctance,
                                        double temperature,
                                        std::optional<double> frequency = std::nullopt);
        // The same crossing on a curve the caller already holds: get_magnetisation_curve builds its cache key
        // from the material's data on every call, which a caller solving many MMFs on one material (one per
        // switching cycle of a line-cycle record, ABT #1764) pays once by fetching the curve itself.
        static std::pair<double, double> get_dc_operating_point_for_magnetomotive_force(const MagnetisationCurve& curve,
                                        double magnetomotiveForce,
                                        double effectiveLength,
                                        double effectiveArea,
                                        double gapReluctance);
        static std::vector<PermeabilityPoint> sample_initial_permeability_by_frequency_modifier(PermeabilityPoint permeabilityPoint);
        static double calculate_frequency_for_initial_permeability_drop(const CoreMaterial& coreMaterial, double percentageDrop, double maximumError = 0.01);
        static std::vector<size_t> get_only_temperature_dependent_indexes(const CoreMaterial& coreMaterial);
        static std::vector<size_t> get_only_temperature_dependent_indexes(const std::vector<PermeabilityPoint>& permeabilityPoints);
        static std::vector<PermeabilityPoint> get_only_temperature_dependent_points(const CoreMaterial& coreMaterial);
        static std::vector<size_t> get_only_frequency_dependent_indexes(const CoreMaterial& coreMaterial);
        static std::vector<size_t> get_only_frequency_dependent_indexes(const std::vector<PermeabilityPoint>& permeabilityPoints);
        static std::vector<PermeabilityPoint> get_only_frequency_dependent_points(const CoreMaterial& coreMaterial);
        static std::vector<size_t> get_only_magnetic_field_dc_bias_dependent_indexes(const CoreMaterial& coreMaterial);
        static std::vector<size_t> get_only_magnetic_field_dc_bias_dependent_indexes(const std::vector<PermeabilityPoint>& permeabilityPoints);
        static std::vector<PermeabilityPoint> get_only_magnetic_field_dc_bias_dependent_points(const CoreMaterial& coreMaterial);
        static std::map<std::string, std::string> get_initial_permeability_equations(PermeabilityPoint permeabilityPoint);
        static double get_initial_permeability_formula(const CoreMaterial& coreMaterial,
                                                       std::optional<double> temperature,
                                                       std::optional<double> magneticFieldDcBias,
                                                       std::optional<double> frequency,
                                                       std::optional<double> magneticFluxDensity,
                                        std::optional<CoreShapeFamily> shapeFamily);
};

} // namespace OpenMagnetics
