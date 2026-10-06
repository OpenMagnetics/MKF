#pragma once
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

        static double get_initial_permeability(CoreMaterial coreMaterial,
                                        std::optional<double> temperature = std::nullopt,
                                        std::optional<double> magneticFieldDcBias = std::nullopt,
                                        std::optional<double> frequency = std::nullopt,
                                        std::optional<double> magneticFluxDensity = std::nullopt,
                                        std::optional<CoreShapeFamily> shapeFamily = std::nullopt);
        static double get_initial_permeability(CoreMaterial coreMaterial, OperatingPoint operatingPoint,
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

        static double has_temperature_dependency(CoreMaterial coreMaterial);
        static double has_frequency_dependency(CoreMaterial coreMaterial);
        static double has_magnetic_field_dc_bias_dependency(CoreMaterial coreMaterial);
        static double get_initial_permeability_temperature_dependent(CoreMaterial coreMaterial, double temperature);
        static double get_initial_permeability_frequency_dependent(CoreMaterial coreMaterial, double frequency);
        static double get_initial_permeability_magnetic_field_dc_bias_dependent(CoreMaterial coreMaterial, double magneticFieldDcBias);
        // ABT #1093: the field strength H_dc that carries a given DC flux density in the material.
        // The datasheet bias curve is the REVERSIBLE (incremental) permeability, dB/dH along the
        // magnetisation curve, so B(H) = µ0·∫0^H µ_rev(h)·dh must be integrated and inverted; the
        // secant B = µ0·µ_rev(H)·H is wrong past the knee, and the fixed point µ ← µ_rev(B/(µ0·µ))
        // diverges there. With refuseAboveSaturation (gap sizing) it throws when biasFluxDensity
        // exceeds the material saturation at that temperature: no gap holds that flux. Without it
        // (analysing a given design) the march continues past the knee on the vacuum slope.
        static double get_magnetic_field_dc_bias_for_flux_density(CoreMaterial coreMaterial,
                                        double biasFluxDensity,
                                        double temperature,
                                        std::optional<double> frequency = std::nullopt,
                                        bool refuseAboveSaturation = true);
        // The material's magnetisation curve B(H) = mu0 * integral of mu_rev, tabulated and cached per material,
        // temperature and frequency (ABT #1223).
        struct MagnetisationCurve {
            std::vector<double> fieldStrength;
            std::vector<double> fluxDensity;
        };
        static const MagnetisationCurve& get_magnetisation_curve(const CoreMaterial& coreMaterial, double temperature, std::optional<double> frequency);
        // The DC operating point {H, B} of the material when the MMF N*I is imposed (current-driven, or a DC
        // set by the load currents): solves N*I = H*l_e + B(H)*A_e*R_gap on the magnetisation curve.
        static std::pair<double, double> get_dc_operating_point_for_magnetomotive_force(CoreMaterial coreMaterial,
                                        double magnetomotiveForce,
                                        double effectiveLength,
                                        double effectiveArea,
                                        double gapReluctance,
                                        double temperature,
                                        std::optional<double> frequency = std::nullopt);
        static std::vector<PermeabilityPoint> sample_initial_permeability_by_frequency_modifier(PermeabilityPoint permeabilityPoint);
        static double calculate_frequency_for_initial_permeability_drop(CoreMaterial coreMaterial, double percentageDrop, double maximumError = 0.01);
        static std::vector<size_t> get_only_temperature_dependent_indexes(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_temperature_dependent_indexes(std::vector<PermeabilityPoint> permeabilityPoints);
        static std::vector<PermeabilityPoint> get_only_temperature_dependent_points(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_frequency_dependent_indexes(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_frequency_dependent_indexes(std::vector<PermeabilityPoint> permeabilityPoints);
        static std::vector<PermeabilityPoint> get_only_frequency_dependent_points(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_magnetic_field_dc_bias_dependent_indexes(CoreMaterial coreMaterial);
        static std::vector<size_t> get_only_magnetic_field_dc_bias_dependent_indexes(std::vector<PermeabilityPoint> permeabilityPoints);
        static std::vector<PermeabilityPoint> get_only_magnetic_field_dc_bias_dependent_points(CoreMaterial coreMaterial);
        static std::map<std::string, std::string> get_initial_permeability_equations(PermeabilityPoint permeabilityPoint);
        static double get_initial_permeability_formula(CoreMaterial coreMaterial,
                                                       std::optional<double> temperature,
                                                       std::optional<double> magneticFieldDcBias,
                                                       std::optional<double> frequency,
                                                       std::optional<double> magneticFluxDensity,
                                        std::optional<CoreShapeFamily> shapeFamily);
};

} // namespace OpenMagnetics
