#include "support/Utils.h"
#include "physical_models/ThermalResistance.h"
#include "support/Exceptions.h"
#include "physical_models/InitialPermeability.h"
#include <magic_enum.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace OpenMagnetics {
using namespace ThermalDefaults; // IMP-2


std::shared_ptr<CoreThermalResistanceModel> CoreThermalResistanceModel::factory(CoreThermalResistanceModels modelName) {
    if (modelName == CoreThermalResistanceModels::MANIKTALA) {
        return std::make_shared<CoreThermalResistanceManiktalaModel>();
    }
    else
        throw ModelNotAvailableException("Unknown Thermal Resistance mode, available options are: {MANIKTALA}");
}

std::shared_ptr<CoreThermalResistanceModel> CoreThermalResistanceModel::factory() {
    return factory(defaults.coreThermalResistanceModelDefault);
}

/**
 * @brief Calculates core thermal resistance using Maniktala's empirical formula.
 */
double CoreThermalResistanceManiktalaModel::get_core_thermal_resistance_reluctance(Core core) {
    double effectiveVolume = core.get_effective_volume();
    // Maniktala's correlation expects cm^3 (see the identical formula in
    // CoreTemperature.cpp); without the conversion Rth was ~1740x too high
    return 53 * pow(effectiveVolume * 1000000, -0.54);
}

//=============================================================================
// Static Thermal Resistance Calculation Methods
//=============================================================================

double ThermalResistance::calculateNaturalConvectionCoefficient(
    double surfaceTemperature,
    double ambientTemperature,
    double characteristicLength,
    SurfaceOrientation orientation) {
    
    // Get air properties at film temperature
    double filmTemp = (surfaceTemperature + ambientTemperature) / 2.0;
    FluidProperties air = FluidProperties::getAirProperties(filmTemp);
    
    // Temperature difference
    double deltaT = std::abs(surfaceTemperature - ambientTemperature);
    if (deltaT < 0.1) {
        deltaT = 0.1;  // Avoid division by zero, assume small gradient
    }
    
    // Grashof number: Gr = g * β * ΔT * L³ / ν²
    Constants constants;
    double Gr = constants.gravityAcceleration * air.thermalExpansionCoefficient * deltaT * 
                std::pow(characteristicLength, 3) / std::pow(air.kinematicViscosity, 2);
    
    // Rayleigh number: Ra = Gr * Pr
    double Ra = Gr * air.prandtlNumber;
    
    double Nu = 0.0;  // Nusselt number
    
    switch (orientation) {
        case SurfaceOrientation::VERTICAL: {
            // Churchill-Chu correlation for vertical surfaces (all Ra)
            double term = std::pow(1.0 + std::pow(0.492 / air.prandtlNumber, 9.0/16.0), 8.0/27.0);
            Nu = std::pow(0.825 + 0.387 * std::pow(Ra, 1.0/6.0) / term, 2);
            break;
        }
        case SurfaceOrientation::HORIZONTAL_TOP: {
            // Hot surface facing up (or cold facing down)
            if (Ra < 1e7) {
                Nu = 0.54 * std::pow(Ra, 0.25);
            } else {
                Nu = 0.15 * std::pow(Ra, 1.0/3.0);
            }
            break;
        }
        case SurfaceOrientation::HORIZONTAL_BOTTOM: {
            // Hot surface facing down (or cold facing up)
            Nu = 0.27 * std::pow(Ra, 0.25);
            break;
        }
    }
    
    // Conduction limit: as Ra -> 0 convection dies but heat still leaves the surface by pure
    // conduction into the surrounding air, whose plate limit is Nu -> O(1) (Nu = 1 recovers
    // h = k_air/L, the stagnant-air conduction film). This is a physical floor, not a tuned one.
    Nu = std::max(Nu, 1.0);

    // Heat transfer coefficient: h = Nu * k / L
    // ABT #527: the historical practical floors (h >= 5, later h >= 2 W/m2K) are removed — they
    // were engineering guards with no derivation, and with the component-scale characteristic
    // length (ABT #461) the correlation plus the conduction limit above needs no patching.
    return Nu * air.thermalConductivity / characteristicLength;
}

double ThermalResistance::calculateForcedConvectionCoefficient(
    double airVelocity,
    double characteristicLength,
    double temperature) {
    
    if (airVelocity <= 0) {
        return calculateNaturalConvectionCoefficient(temperature + 20, temperature, 
                                                     characteristicLength, SurfaceOrientation::VERTICAL);
    }
    
    FluidProperties air = FluidProperties::getAirProperties(temperature);
    
    // Reynolds number: Re = V * L / ν
    double Re = airVelocity * characteristicLength / air.kinematicViscosity;
    
    double Nu = 0.0;
    
    if (Re < 5e5) {
        // Laminar flow: Nu = 0.664 * Re^0.5 * Pr^(1/3)
        Nu = 0.664 * std::pow(Re, 0.5) * std::pow(air.prandtlNumber, 1.0/3.0);
    } else {
        // Turbulent flow: Nu = 0.037 * Re^0.8 * Pr^(1/3)
        Nu = 0.037 * std::pow(Re, 0.8) * std::pow(air.prandtlNumber, 1.0/3.0);
    }
    
    // Heat transfer coefficient
    double h = Nu * air.thermalConductivity / characteristicLength;
    
    // Typical range for forced convection: 25-250 W/(m²·K)
    // ABT #527: no practical floor — the correlation governs. At zero/low velocity the caller
    // uses the natural-convection path instead (see calculateForcedConvectionCoefficient's
    // velocity guard above), so this result is never asked to stand in for natural convection.
    return h;
}

double ThermalResistance::calculateRadiationCoefficient(
    double surfaceTemperature,
    double ambientTemperature,
    double emissivity) {
    
    // Convert to Kelvin
    Constants constants;
    double Ts = surfaceTemperature + constants.kelvinOffset;
    double Ta = ambientTemperature + constants.kelvinOffset;
    
    // Linearized radiation coefficient: h_rad = ε * σ * (Ts² + Ta²) * (Ts + Ta)
    double h_rad = emissivity * constants.stefanBoltzmannConstant * (Ts*Ts + Ta*Ta) * (Ts + Ta);
    
    return h_rad;
}

double ThermalResistance::calculateCylinderWallToOpeningsViewFactor(double radius, double height) {
    if (!(radius > 0) || !(height > 0)) {
        throw std::invalid_argument("ThermalResistance::calculateCylinderWallToOpeningsViewFactor: the cylinder needs "
                                    "a positive radius and height (radius " + std::to_string(radius) + " m, height " +
                                    std::to_string(height) + " m).");
    }
    const double heightOverRadius = height / radius;
    const double X = 2.0 + heightOverRadius * heightOverRadius;
    const double diskToDisk = 0.5 * (X - std::sqrt(X * X - 4.0));
    return (radius / height) * (1.0 - diskToDisk);
}

// IMP-NEW-08: View-factor-aware radiation coefficient
double ThermalResistance::calculateRadiationCoefficientWithViewFactor(
    double surfaceTemperature, double ambientTemperature,
    double emissivity, double viewFactor) {
    double h_rad = calculateRadiationCoefficient(surfaceTemperature, ambientTemperature, emissivity);
    return h_rad * std::clamp(viewFactor, 0.0, 1.0);
}

double ThermalResistance::getMaterialThermalConductivity(const std::string& materialName) {
    // Convert to lowercase for database lookup
    std::string lowerName = materialName;
    std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(), ::tolower);
    
    // Try wire material database (copper, aluminium have thermalConductivity as array)
    try {
        auto wireMaterial = find_wire_material_by_name(lowerName);
        // Use the interpolation function at room temperature (25°C)
        return getWireMaterialThermalConductivity(wireMaterial, 25.0);
    } catch (const std::exception& /* e */) { // IMP-8
        // Not a wire material, continue
    }
    
    // Try insulation material database. Its keys are the MAS names as written ("Kapton HN"), so
    // the name as given is tried before the lower-cased one: searching only the lower-cased name
    // never matched a mixed-case material, which fell through to the generic list and threw
    // MATERIAL_DATA_MISSING for a material MAS does carry a conductivity for.
    for (const auto& name : {materialName, lowerName}) {
        try {
            auto insulationMaterial = find_insulation_material_by_name(name);
            auto thermalCond = insulationMaterial.get_thermal_conductivity();
            if (thermalCond) {
                return *thermalCond;
            }
        } catch (const std::exception& /* e */) { // IMP-8
            // Not an insulation material under this spelling, continue
        }
    }
    
    // Try core material database. A core material resolves through the one core-material
    // resolver (MAS value, or a flagged class/family/all-sourced estimate) — not through get_nominal(),
    // which missed every min/max band and sent such materials to the generic list below.
    std::optional<CoreMaterial> coreMaterial;
    try {
        coreMaterial = find_core_material_by_name(materialName);  // Core materials may be case-sensitive
    } catch (const CoreMaterialNotFoundException& /* e */) {
        // Not a core material, continue
    }
    if (coreMaterial) {
        return resolve_core_material_heat_conductivity(coreMaterial.value()).value;
    }
    
    // Generic NON-core materials in W/(m·K) (core materials resolve above, through the resolver)
    static const std::map<std::string, double> conductivities = {
        {"copper", 385.0},
        {"aluminium", 237.0},
        {"aluminum", 237.0},
        {"air", 0.026},
        {"epoxy", 0.25},
        {"polyamide", 0.25},
        {"nylon", 0.25},
        {"pet", 0.15},
        {"pbt", 0.2},
        {"lcp", 0.3},
        {"mica", 0.5},
        {"kapton", 0.12},
        {"mylar", 0.15},
        {"teflon", 0.25},
        {"silicone", 0.2},
        {"thermal_compound", 1.0},
        {"solder", 50.0}
    };
    
    std::string name = materialName;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    
    auto it = conductivities.find(name);
    if (it != conductivities.end()) {
        return it->second;
    }

    // ABT #837: this used to `return 1.0` for ANY unrecognised name — a silent fallback of
    // exactly the kind the no-fallbacks rule forbids, and a badly chosen one: 1.0 W/(m*K) is
    // ~5x too conductive for the plastics this path actually sees (PET 0.15, PBT 0.2,
    // polyamide 0.25), so an unknown insulation material made the part look BETTER cooled
    // than it is. Nothing downstream could tell the guess from a real datasheet value.
    // Name the material and refuse instead; the map above is the list of what is known.
    throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
        "No thermal conductivity for material '" + materialName + "': it is not in the MAS "
        "material data and not one of the known generic materials. Add its thermal "
        "conductivity to MAS, or use a material that carries one — a thermal result derived "
        "from a guessed conductivity is not worth reporting.");
}

double ThermalResistance::getWireMaterialThermalConductivity(const WireMaterial& wireMaterial, double temperature) {
    // Get thermal conductivity from MAS database with temperature interpolation
    auto thermalConductivityOpt = wireMaterial.get_thermal_conductivity();
    
    if (!thermalConductivityOpt) {
        // Fallback to material name lookup  
        std::string materialName = wireMaterial.get_name();
        return getMaterialThermalConductivity(materialName);
    }
    
    auto& thermalCond = *thermalConductivityOpt;
    
    // Check if it's a constant value (just one entry)
    if (thermalCond.size() == 1) {
        return thermalCond[0].get_value();
    }
    
    // Interpolate based on temperature
    // thermalCond is array of {temperature, value} pairs
    if (thermalCond.size() == 0) {
        // An empty (but present) table is malformed data: silently treating it as
        // copper would be ~40% wrong for an aluminium wire
        throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING, "Wire material " + wireMaterial.get_name() + " has an empty thermal conductivity table");
    }
    
    // Find bracketing temperatures
    for (size_t i = 0; i < thermalCond.size() - 1; ++i) {
        double T1 = thermalCond[i].get_temperature();
        double T2 = thermalCond[i+1].get_temperature();
        double k1 = thermalCond[i].get_value();
        double k2 = thermalCond[i+1].get_value();
        
        if (temperature >= T1 && temperature <= T2) {
            // Linear interpolation
            double alpha = (temperature - T1) / (T2 - T1);
            return k1 + alpha * (k2 - k1);
        }
    }
    
    // If temperature is out of range, use nearest endpoint
    if (temperature < thermalCond[0].get_temperature()) {
        return thermalCond[0].get_value();
    } else {
        return thermalCond.back().get_value();
    }
}

namespace {

// THE APPROVED FALLBACK for a core material's HEAT CAPACITY (Alf, 2026-10-11): a material
// that carries no value of its own takes the MEDIAN of that value over every MAS core material of
// the SAME materialComposition that does carry one. Each contributing value is the manufacturer's
// own, with its provenance in MAS, so the fallback is sourced data rather than a constant, and it
// follows the database as data is added. The median, not the mean, so one outlying vendor table
// cannot drag a whole chemistry. (Heat CONDUCTIVITY does not use this: it goes through the tiered
// resolve_core_material_heat_conductivity below, ABT #1753, Alf's decision 2026-10-12.)
//
// It throws (MATERIAL_DATA_MISSING) when the material states no composition, when the composition
// is 'proprietary' (no chemistry to pool by), or when no material of that composition carries the
// value -- never a guessed number.
template <typename Getter>
double median_over_composition(const CoreMaterial& coreMaterial, Getter getValue, const std::string& quantity) {
    const auto composition = coreMaterial.get_material_composition();
    if (!composition || *composition == MaterialComposition::PROPRIETARY) {
        throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
            "Core material '" + coreMaterial.get_name() + "' has no " + quantity +
            " and no materialComposition to take the same-composition median from");
    }
    if (coreMaterialDatabase.empty()) {
        load_core_materials();
    }
    std::vector<double> values;
    for (const auto& [name, material] : coreMaterialDatabase) {
        if (material.get_material_composition() != composition) {
            continue;
        }
        if (const auto value = getValue(material)) {
            values.push_back(resolve_dimensional_values(*value));
        }
    }
    if (values.empty()) {
        throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
            "Core material '" + coreMaterial.get_name() + "' has no " + quantity + ", and no " +
            std::string(magic_enum::enum_name(*composition)) + " material in the database carries one");
    }
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    return values.size() % 2 == 1 ? values[middle] : 0.5 * (values[middle - 1] + values[middle]);
}

}  // namespace

double ThermalResistance::getCoreMaterialHeatCapacity(const CoreMaterial& coreMaterial) {
    if (const auto heatCapacity = coreMaterial.get_heat_capacity()) {
        return resolve_dimensional_values(*heatCapacity);
    }
    return median_over_composition(coreMaterial, [](const CoreMaterial& m) { return m.get_heat_capacity(); },
                                   "heatCapacity");
}

double ThermalResistance::getCoreMaterialThermalConductivity(const CoreMaterial& coreMaterial) {
    return resolve_core_material_heat_conductivity(coreMaterial).value;
}

// ============================================================================
// Core material heat conductivity: MAS value, or a flagged estimate from the sourced members of
// the material's class, family, or all sourced materials (ABT #1753). Alf, 2026-10-11: "create an interpolator for the ones that
// don't have data, from the ones that have data, as an approved fallback".
// ============================================================================
namespace {

// A heatConductivity is SOURCED when a provenance entry that is not retracted names the field.
// A bare value without one (e.g. copied between vendors) is used for its own material, but it is
// not evidence about any other material.
std::vector<std::string> heat_conductivity_sources(const CoreMaterial& material) {
    std::vector<std::string> sources;
    if (!material.get_heat_conductivity() || !material.get_provenance()) {
        return sources;
    }
    const auto provenance = material.get_provenance().value();
    for (const auto& entry : provenance) {
        if (entry.get_retracted() && entry.get_retracted().value()) {
            continue;
        }
        if (!entry.get_fields()) {
            continue;
        }
        const auto fields = entry.get_fields().value();
        if (std::find(fields.begin(), fields.end(), "heatConductivity") == fields.end()) {
            continue;
        }
        if (entry.get_source_url()) {
            sources.push_back(entry.get_source_url().value());
        }
        else if (entry.get_source_name()) {
            sources.push_back(entry.get_source_name().value());
        }
        // An entry naming the field but no document is not a citation.
    }
    return sources;
}

// The first document citing the value (the MATERIAL tier reports one source).
std::optional<std::string> heat_conductivity_source(const CoreMaterial& material) {
    const auto sources = heat_conductivity_sources(material);
    if (sources.empty()) {
        return std::nullopt;
    }
    return sources.front();
}

// Two-sided 95 % critical value of Student's t for `dof` degrees of freedom (standard table;
// statistics, not physics). Beyond 30 dof the normal value is within 4 %.
double student_t_975(size_t dof) {
    static const double table[] = {12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
                                   2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
                                   2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042};
    if (dof == 0) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Student t with zero degrees of freedom");
    }
    if (dof <= 30) {
        return table[dof - 1];
    }
    return 1.960;
}

struct LinearTrend {
    double intercept;
    double slope;
    double rSquared;
    double minimum;
    double maximum;
};

// Least squares k = a + b x; engaged only when the slope is significant at 95 % (|b| > t·SE(b)).
std::optional<LinearTrend> significant_linear_trend(const std::vector<std::pair<double, double>>& points) {
    const size_t n = points.size();
    if (n < 3) {
        return std::nullopt;
    }
    double meanX = 0, meanY = 0;
    for (const auto& [x, y] : points) {
        meanX += x;
        meanY += y;
    }
    meanX /= static_cast<double>(n);
    meanY /= static_cast<double>(n);
    double sxx = 0, sxy = 0, syy = 0;
    for (const auto& [x, y] : points) {
        sxx += (x - meanX) * (x - meanX);
        sxy += (x - meanX) * (y - meanY);
        syy += (y - meanY) * (y - meanY);
    }
    if (sxx <= 0 || syy <= 0) {
        return std::nullopt;
    }
    const double slope = sxy / sxx;
    const double intercept = meanY - slope * meanX;
    double residual = 0;
    for (const auto& [x, y] : points) {
        const double e = y - (intercept + slope * x);
        residual += e * e;
    }
    const double standardError = std::sqrt(residual / static_cast<double>(n - 2) / sxx);
    if (std::abs(slope) <= student_t_975(n - 2) * standardError) {
        return std::nullopt;
    }
    double minimum = std::numeric_limits<double>::max();
    double maximum = std::numeric_limits<double>::lowest();
    for (const auto& point : points) {
        minimum = std::min(minimum, point.first);
        maximum = std::max(maximum, point.first);
    }
    return LinearTrend{intercept, slope, sxy * sxy / (sxx * syy), minimum, maximum};
}

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    return n % 2 == 1 ? values[n / 2] : (values[n / 2 - 1] + values[n / 2]) / 2.0;
}

// The predictors a trend may follow, in the order they are tried: density first (it tracks the
// porosity, i.e. how much binder/air sits between conducting particles), then initial permeability
// (which tracks the same packing for powders).
std::optional<double> predictor_value(const std::string& predictor, const CoreMaterial& material) {
    if (predictor == "density") {
        return material.get_density();
    }
    return InitialPermeability::get_initial_permeability(material, Defaults().ambientTemperature);
}

// A heatConductivity that is present must be a finite positive number: anything else is malformed
// data, which throws (absence is not malformed; it walks the estimate tiers).
double checked_heat_conductivity(const CoreMaterial& material) {
    const double value = resolve_dimensional_values(material.get_heat_conductivity().value());
    if (!std::isfinite(value) || value <= 0) {
        std::ostringstream message;
        message << "Core material '" << material.get_name() << "' has a malformed heatConductivity in MAS (" << value
                << " W/(m K)): it must be a finite positive number. Fix the record in MAS.";
        throw MaterialException(ErrorCode::MATERIAL_INVALID_PROPERTY, message.str());
    }
    return value;
}

bool is_ferrite_family(MaterialType materialType) {
    return materialType == MaterialType::FERRITE;
}

std::string class_name(MaterialType materialType, MaterialComposition composition) {
    return std::string(magic_enum::enum_name(materialType)) + " " + std::string(magic_enum::enum_name(composition));
}

struct SourcedMember {
    const CoreMaterial* material;
    std::vector<std::string> sources;  // every document citing the value
    double value;
};

// Every core material whose heatConductivity a provenance entry cites. Unsourced values are not
// evidence about any other material and are left out.
std::vector<SourcedMember> sourced_core_materials() {
    if (coreMaterialDatabase.empty()) {
        load_core_materials();
    }
    std::vector<SourcedMember> members;
    for (const auto& [materialName, material] : coreMaterialDatabase) {
        auto sources = heat_conductivity_sources(material);
        if (sources.empty()) {
            continue;
        }
        members.push_back({&material, std::move(sources), checked_heat_conductivity(material)});
    }
    return members;
}

// Pools the members, collapsing family values: a document that gives one value to all its grades
// is ONE data point (else a vendor with many grades would outvote the rest); a document that gives
// its grades different values is per-grade data, one point per grade, and only those can show a
// trend. A value several documents cite is counted once, under the citing document that covers
// the most members of the pool (ties: the smallest identifier), so the grouping does not depend on
// the order of the provenance entries, nor split one table cited under two URLs.
// Fills the estimate's materials/documents/points and returns the data points.
std::vector<double> pool_data_points(const std::vector<SourcedMember>& members, CoreHeatConductivityEstimate& estimate,
                                     std::vector<const SourcedMember*>* perGradeMembers) {
    std::map<std::string, size_t> citationsInPool;
    for (const auto& member : members) {
        for (const auto& source : std::set<std::string>(member.sources.begin(), member.sources.end())) {
            ++citationsInPool[source];
        }
    }
    std::map<std::string, std::vector<const SourcedMember*>> membersBySource;
    for (const auto& member : members) {
        const std::string* widest = nullptr;
        for (const auto& source : member.sources) {
            if (!widest || citationsInPool.at(source) > citationsInPool.at(*widest) ||
                (citationsInPool.at(source) == citationsInPool.at(*widest) && source < *widest)) {
                widest = &source;
            }
        }
        membersBySource[*widest].push_back(&member);
    }
    std::vector<double> points;
    estimate.sourceMaterials.clear();
    for (const auto& [source, sourceMembers] : membersBySource) {
        std::set<double> distinct;
        for (const auto* member : sourceMembers) {
            distinct.insert(member->value);
            estimate.sourceMaterials.push_back(member->material->get_name());
        }
        if (distinct.size() == 1) {
            points.push_back(sourceMembers.front()->value);
            continue;
        }
        for (const auto* member : sourceMembers) {
            points.push_back(member->value);
            if (perGradeMembers) {
                perGradeMembers->push_back(member);
            }
        }
    }
    std::sort(estimate.sourceMaterials.begin(), estimate.sourceMaterials.end());
    estimate.numberDocuments = membersBySource.size();
    estimate.numberDataPoints = points.size();
    return points;
}

// Tier CLASS: in-span trend, else the class median. Disengaged when the class has no sourced member.
std::optional<CoreHeatConductivity> class_estimate(const CoreMaterial& coreMaterial, MaterialComposition composition,
                                                   const std::vector<SourcedMember>& sourced) {
    const auto materialType = coreMaterial.get_material();
    std::vector<SourcedMember> members;
    for (const auto& member : sourced) {
        if (member.material->get_material() == materialType && member.material->get_material_composition() &&
            member.material->get_material_composition().value() == composition) {
            members.push_back(member);
        }
    }
    if (members.empty()) {
        return std::nullopt;
    }
    CoreHeatConductivityEstimate estimate;
    estimate.tier = CoreHeatConductivityTier::CLASS;
    estimate.pool = "class " + class_name(materialType, composition);
    std::vector<const SourcedMember*> trendMembers;
    const auto points = pool_data_points(members, estimate, &trendMembers);

    double value = median(points);
    std::string method = "the median of its " + std::to_string(points.size()) + " sourced data points";

    std::optional<LinearTrend> bestTrend;
    std::string bestPredictor;
    double bestPredictorValue = 0;
    for (const std::string predictor : {"density", "initialPermeability"}) {
        if (trendMembers.size() < 3) {
            break;
        }
        const auto target = predictor_value(predictor, coreMaterial);
        if (!target) {
            continue;
        }
        std::vector<std::pair<double, double>> trendPoints;
        for (const auto* member : trendMembers) {
            const auto x = predictor_value(predictor, *member->material);
            if (x) {
                trendPoints.push_back({x.value(), member->value});
            }
        }
        const auto trend = significant_linear_trend(trendPoints);
        if (trend && (!bestTrend || trend->rSquared > bestTrend->rSquared)) {
            bestTrend = trend;
            bestPredictor = predictor;
            bestPredictorValue = target.value();
        }
    }
    if (bestTrend) {
        std::ostringstream span;
        span << bestPredictor << " = " << bestPredictorValue << " against the sourced span " << bestTrend->minimum
             << " to " << bestTrend->maximum;
        // Alf, 2026-10-11: follow the trend only inside the sourced span; outside it, the class median.
        if (bestPredictorValue >= bestTrend->minimum && bestPredictorValue <= bestTrend->maximum) {
            value = bestTrend->intercept + bestTrend->slope * bestPredictorValue;
            estimate.predictor = bestPredictor;
            estimate.predictorValue = bestPredictorValue;
            estimate.predictorMinimum = bestTrend->minimum;
            estimate.predictorMaximum = bestTrend->maximum;
            std::ostringstream trendText;
            trendText << "the least-squares trend on " << bestPredictor << " (r^2 = " << bestTrend->rSquared << ") at "
                      << span.str();
            method = trendText.str();
        }
        else {
            estimate.trendNotFollowedReason = "outside the trend's span: " + span.str();
            method += " (the class's " + bestPredictor + " trend is not followed: " + span.str() + ")";
        }
    }
    estimate.method = method;
    return CoreHeatConductivity{value, CoreHeatConductivityTier::CLASS, std::nullopt, estimate};
}

// Tiers FAMILY and ALL_SOURCED: the median of a broader pool. Disengaged when the pool is empty.
std::optional<CoreHeatConductivity> pooled_median_estimate(CoreHeatConductivityTier tier, const std::string& pool,
                                                           const std::vector<SourcedMember>& members) {
    if (members.empty()) {
        return std::nullopt;
    }
    CoreHeatConductivityEstimate estimate;
    estimate.tier = tier;
    estimate.pool = pool;
    const auto points = pool_data_points(members, estimate, nullptr);
    estimate.method = "the median of its " + std::to_string(points.size()) + " sourced data points";
    return CoreHeatConductivity{median(points), tier, std::nullopt, estimate};
}

CoreHeatConductivity estimate_core_material_heat_conductivity(const CoreMaterial& coreMaterial) {
    const auto materialType = coreMaterial.get_material();
    const auto sourced = sourced_core_materials();
    std::optional<CoreHeatConductivity> resolved;
    std::string skipped;

    // Tier CLASS: a physical composition only (a proprietary one is not a class).
    const auto composition = coreMaterial.get_material_composition();
    const bool physicalComposition = composition && composition.value() != MaterialComposition::PROPRIETARY;
    if (physicalComposition) {
        resolved = class_estimate(coreMaterial, composition.value(), sourced);
        if (!resolved) {
            skipped = "class " + class_name(materialType, composition.value()) + " has no sourced member";
        }
    }
    else {
        skipped = composition ? "its composition is PROPRIETARY" : "it has no materialComposition";
    }

    // Tier FAMILY: ferrites, or metal powders and alloys; a proprietary or missing composition only
    // says the material type, so its family is that type.
    if (!resolved) {
        std::vector<SourcedMember> members;
        std::string pool;
        if (physicalComposition) {
            const bool ferrite = is_ferrite_family(materialType);
            pool = ferrite ? "ferrite family (material type FERRITE)"
                           : "metal powder and alloy family (material types other than FERRITE)";
            for (const auto& member : sourced) {
                if (is_ferrite_family(member.material->get_material()) == ferrite) {
                    members.push_back(member);
                }
            }
        }
        else {
            pool = "material type " + std::string(magic_enum::enum_name(materialType));
            for (const auto& member : sourced) {
                if (member.material->get_material() == materialType) {
                    members.push_back(member);
                }
            }
        }
        resolved = pooled_median_estimate(CoreHeatConductivityTier::FAMILY, pool, members);
        if (!resolved) {
            skipped += "; " + pool + " has no sourced member";
        }
    }

    // Tier ALL_SOURCED: last resort.
    if (!resolved) {
        resolved = pooled_median_estimate(CoreHeatConductivityTier::ALL_SOURCED, "all sourced core materials", sourced);
    }
    if (!resolved) {
        throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
            "Core material '" + coreMaterial.get_name() + "' has no heatConductivity in MAS, and no core material in "
            "MAS carries a sourced heatConductivity to estimate it from: the MAS core material data is not loaded "
            "or has lost its provenance.");
    }
    if (!skipped.empty()) {
        resolved->estimate->method += " (" + skipped + ")";
    }
    return resolved.value();
}

// Where a material's conductivity came from is a note about the data, not about the user's design:
// it is logged at INFO, once per (kind of note, material) per session. Keying on the kind as well
// keeps a different note about the same material from being swallowed. The value, the result's own
// flags and estimate, and every exception are untouched.
void note_provenance_once(const std::string& kind, const std::string& material, const std::string& message) {
    static std::mutex notedMutex;
    static std::set<std::pair<std::string, std::string>> noted;
    std::lock_guard<std::mutex> lock(notedMutex);
    if (noted.insert({kind, material}).second) {
        OM_INFO_M("ThermalResistance", message);
    }
}

void note_once(const CoreMaterial& coreMaterial, const std::string& kind, const std::string& message) {
    note_provenance_once("core " + kind, coreMaterial.get_manufacturer_info().get_name() + "/" + coreMaterial.get_name(),
                         message);
}

} // namespace

CoreHeatConductivity resolve_core_material_heat_conductivity(const CoreMaterial& coreMaterial) {
    if (coreMaterial.get_heat_conductivity()) {
        CoreHeatConductivity resolved{checked_heat_conductivity(coreMaterial), CoreHeatConductivityTier::MATERIAL,
                                      heat_conductivity_source(coreMaterial), std::nullopt};
        if (resolved.is_unsourced()) {
            std::ostringstream message;
            message << "Core material '" << coreMaterial.get_name() << "': its MAS heatConductivity " << resolved.value
                    << " W/(m K) is UNSOURCED (no provenance entry cites it); used for this material, never pooled "
                       "for others.";
            note_once(coreMaterial, "unsourced", message.str());
        }
        return resolved;
    }
    auto resolved = estimate_core_material_heat_conductivity(coreMaterial);
    const auto& estimate = resolved.estimate.value();
    std::ostringstream message;
    message << "Core material '" << coreMaterial.get_name() << "' has no heatConductivity in MAS: ESTIMATED "
            << resolved.value << " W/(m K), tier " << magic_enum::enum_name(estimate.tier) << ", from the "
            << estimate.pool << ": " << estimate.sourceMaterials.size() << " sourced materials, "
            << estimate.numberDocuments << " documents, " << estimate.numberDataPoints
            << " independent data points, using " << estimate.method << ".";
    note_once(coreMaterial, "estimated", message.str());
    return resolved;
}


// ============================================================================
// Bobbin material thermal conductivity: the MAS value, or a flagged estimate from the sourced
// data points of the material's polymer class (ABT #1643, Alf 2026-10-18: "use a class-based
// estimator like the approved #1752/#1753 core one: k from the polymer class, flagged as an
// estimate in the results. MAS keeps only sourced values").
// ============================================================================
namespace {

// The documents a non-retracted provenance entry naming `field` cites (as the core-material
// heat_conductivity_sources does): an entry naming the field but no document is not a citation.
std::vector<std::string> provenance_sources(const std::optional<std::vector<MAS::Provenance>>& provenance,
                                            const std::string& field) {
    std::vector<std::string> sources;
    if (!provenance) {
        return sources;
    }
    for (const auto& entry : provenance.value()) {
        if (entry.get_retracted() && entry.get_retracted().value()) {
            continue;
        }
        if (!entry.get_fields()) {
            continue;
        }
        const auto fields = entry.get_fields().value();
        if (std::find(fields.begin(), fields.end(), field) == fields.end()) {
            continue;
        }
        if (entry.get_source_url()) {
            sources.push_back(entry.get_source_url().value());
        }
        else if (entry.get_source_name()) {
            sources.push_back(entry.get_source_name().value());
        }
    }
    return sources;
}


// A polymer class: the base polymer plus whether it is glass-fibre reinforced. Fibre content moves
// a moulding compound's conductivity well beyond the spread between grades of one polymer (DuPont
// Rynite PET: 0.29 at GF30, 0.33 at GF55), so unfilled and filled grades are separate classes.
struct BobbinPolymerClass {
    std::string polymer;  // "PA66", "PET/PBT", "LCP", "PPS", "PHENOLIC"
    bool glassFibre;
    bool operator==(const BobbinPolymerClass&) const = default;
    std::string name() const { return polymer + (glassFibre ? ", glass-fibre reinforced" : ", unreinforced"); }
};

// A thermal-conductivity data point MKF cites for a class MAS has no sourced member of. Each is
// read from a maker's document of the same polymer class; the document is named so it can be
// checked, and its value is never written into MAS (MAS keeps only the values of its own records).
struct BobbinDatasheetPoint {
    BobbinPolymerClass polymerClass;
    std::string grade;
    double value;  // W/(m K)
    std::string document;
};

// THE CITED DATASHEET POINTS (retrieved 2026-10-18; copies kept in the ABT #1643 run directory).
//
// PET/PBT, glass-fibre reinforced (for SKYT.5220FR = Skytra 5220FR, PET-GF30 FR; besides the
// sourced MAS polyesters, Pocan B4225 0.23 and Valox 420SE0 0.25): DuPont "Crastin PBT and Rynite PET thermoplastic polyester resins, Design
// Information, Module IV" (L-10978), section 2 property tables. One document giving each grade
// its own value is one point per grade. Rynite rows by ASTM C-177, Crastin rows by DIN 51046; the
// unreinforced, impact-modified and PBT-ASA blend columns of the same tables are left out (not the
// class). The FR530L row is the closest grade to Skytra 5220FR (PET-GF30 FR); the class median,
// 0.29, is what is used.
//
// LCP, glass-fibre reinforced (for Zen.6130L = Zenite 6130L, LCP-GF30; MAS has no sourced LCP):
// Celanese VECTRA E130i (LCP-GF30) data sheet, "Thermal conductivity, flow 0.348 W/(m K), ISO
// 22007-2", marked [OT] one-time tested. Caveats: it is the IN-FLOW value, and a moulded LCP wall
// conducts less across its thickness than along the flow, so this over-states the through-wall
// conductance of an LCP bobbin; and the copy read is a distributor-hosted one of the Celanese
// sheet. It is the only maker value of the class found.
//
// PA66-GF and PPS-GF take their points from MAS (BASF Ultramid A3X2G5 0.33 and A3X2G10 0.35, both
// DIN 52612-1, and the GFN (Natural) representative grade Ultramid A3EG7 0.35; Toray Torelina
// A504X90 0.3 through-thickness), so no datasheet point is cited here. A PA66-GF material without
// its own value (PA66 = Durethan AKV30H2.0, PA66-GF30) is therefore estimated at the median of
// those three BASF points, 0.35 (MAS c1b29c0). Phenolic (PF) has no sourced point in MAS or here,
// so a PF bobbin without its own value throws.
const std::vector<BobbinDatasheetPoint>& bobbin_datasheet_points() {
    static const std::string crastinRynite =
        "DuPont, Crastin PBT and Rynite PET Design Information, Module IV (L-10978), section 2: "
        "https://www.distrupol.com/Crastin_PBT_and_Rynite_PET_Design_Info_Module_IV.pdf";
    static const std::string vectra =
        "Celanese VECTRA E130i data sheet (LCP-GF30), thermal conductivity flow, ISO 22007-2: "
        "http://hongrunplastics.com/public/uploads/images/20250611/Celanese%20LCP%20VECTRA%20E130i.pdf";
    static const BobbinPolymerClass polyesterGf{"PET/PBT", true};
    static const BobbinPolymerClass lcpGf{"LCP", true};
    static const std::vector<BobbinDatasheetPoint> points = {
        {polyesterGf, "Rynite 530 (PET-GF30)", 0.29, crastinRynite},
        {polyesterGf, "Rynite 545 (PET-GF45)", 0.32, crastinRynite},
        {polyesterGf, "Rynite 555 (PET-GF55)", 0.33, crastinRynite},
        {polyesterGf, "Rynite FR530L (PET-GF30 FR)", 0.25, crastinRynite},
        {polyesterGf, "Crastin SK605 (PBT-GF30)", 0.28, crastinRynite},
        {polyesterGf, "Crastin SK645FR (PBT-GF30 FR)", 0.29, crastinRynite},
        {polyesterGf, "Crastin T845FR (co-PBT GF30 FR)", 0.29, crastinRynite},
        {lcpGf, "VECTRA E130i (LCP-GF30)", 0.348, vectra},
    };
    return points;
}

std::string lower_case(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) { return std::tolower(c); });
    return text;
}

// The class of a bobbin material, read from its MAS composition ("Polyamide 6.6, 30% glass
// fibre", "Polyethylene terephthalate, 30% glass fibre, flame retardant", "Liquid crystal polymer,
// 30% glass fibre, lubricated", ...). Disengaged when the composition does not say it.
std::optional<BobbinPolymerClass> bobbin_polymer_class(const InsulationMaterial& material, std::string& why) {
    if (!material.get_composition()) {
        why = "its MAS record states no composition";
        return std::nullopt;
    }
    const auto composition = lower_case(material.get_composition().value());
    std::optional<std::string> polymer;
    auto names = [&](std::initializer_list<const char*> keys) {
        for (auto key : keys) {
            if (composition.find(key) != std::string::npos) {
                return true;
            }
        }
        return false;
    };
    if (names({"polyamide 6.6", "polyamide 66", "pa66", "pa 66"})) {
        polymer = "PA66";
    }
    else if (names({"polyethylene terephthalate", "polybutylene terephthalate"})) {
        polymer = "PET/PBT";
    }
    else if (names({"liquid crystal polymer"})) {
        polymer = "LCP";
    }
    else if (names({"polyphenylene sulfide", "polyphenylene sulphide"})) {
        polymer = "PPS";
    }
    else if (names({"phenol"})) {
        polymer = "PHENOLIC";
    }
    if (!polymer) {
        why = "its composition '" + material.get_composition().value() + "' names no polymer class MKF estimates "
              "(PA66, PET/PBT, LCP, PPS, phenolic)";
        return std::nullopt;
    }
    const bool glass = composition.find("glass") != std::string::npos;
    const bool glassFibre = names({"glass fibre", "glass fiber"});
    if (glass && !glassFibre) {
        why = "its composition '" + material.get_composition().value() + "' states glass in a form other than "
              "fibre, which no class covers";
        return std::nullopt;
    }
    return BobbinPolymerClass{polymer.value(), glassFibre};
}

// A thermalConductivity that is present must be a finite positive number.
double checked_insulation_thermal_conductivity(const InsulationMaterial& material) {
    const double value = material.get_thermal_conductivity().value();
    if (!std::isfinite(value) || value <= 0) {
        std::ostringstream message;
        message << "Insulation material '" << material.get_name() << "' has a malformed thermalConductivity in MAS ("
                << value << " W/(m K)): it must be a finite positive number. Fix the record in MAS.";
        throw MaterialException(ErrorCode::MATERIAL_INVALID_PROPERTY, message.str());
    }
    return value;
}

void note_once_bobbin(const InsulationMaterial& material, const std::string& kind, const std::string& message) {
    note_provenance_once("bobbin " + kind, material.get_name(), message);
}

}  // namespace

BobbinThermalConductivity resolve_bobbin_material_thermal_conductivity(const InsulationMaterial& material) {
    if (material.get_thermal_conductivity()) {
        const auto sources = provenance_sources(material.get_provenance(), "thermalConductivity");
        BobbinThermalConductivity resolved{checked_insulation_thermal_conductivity(material),
                                           BobbinThermalConductivityTier::MATERIAL,
                                           sources.empty() ? std::nullopt : std::optional<std::string>(sources.front()),
                                           std::nullopt};
        if (resolved.is_unsourced()) {
            std::ostringstream message;
            message << "Bobbin material '" << material.get_name() << "': its MAS thermalConductivity " << resolved.value
                    << " W/(m K) is UNSOURCED (no provenance entry cites it); used for this material, never pooled "
                       "for others.";
            note_once_bobbin(material, "unsourced", message.str());
        }
        return resolved;
    }

    std::string why;
    const auto polymerClass = bobbin_polymer_class(material, why);
    if (!polymerClass) {
        throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
            "Bobbin material '" + material.get_name() + "' has no thermalConductivity in MAS, and its polymer class "
            "cannot be determined to estimate one: " + why + ".");
    }

    // The class's data points: one per (document, value) -- a document giving one value to several
    // grades or records is one point, as in the core-material estimator.
    BobbinThermalConductivityEstimate estimate;
    estimate.polymerClass = polymerClass->name();
    std::set<std::pair<std::string, double>> points;
    std::set<std::string> documents;
    if (insulationMaterialDatabase.empty()) {
        load_insulation_materials();
    }
    for (const auto& [name, member] : insulationMaterialDatabase) {
        if (!member.get_thermal_conductivity()) {
            continue;
        }
        const auto sources = provenance_sources(member.get_provenance(), "thermalConductivity");
        if (sources.empty()) {
            continue;  // an unsourced value is not evidence about another material
        }
        std::string ignored;
        const auto memberClass = bobbin_polymer_class(member, ignored);
        if (!memberClass || !(memberClass.value() == polymerClass.value())) {
            continue;
        }
        const double value = checked_insulation_thermal_conductivity(member);
        estimate.sourceMaterials.push_back(member.get_name());
        for (const auto& source : sources) {
            documents.insert(source);
        }
        points.insert({sources.front(), value});
    }
    for (const auto& point : bobbin_datasheet_points()) {
        if (!(point.polymerClass == polymerClass.value())) {
            continue;
        }
        estimate.datasheetGrades.push_back(point.grade);
        documents.insert(point.document);
        points.insert({point.document + " / " + point.grade, point.value});
    }
    if (points.empty()) {
        throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
            "Bobbin material '" + material.get_name() + "' has no thermalConductivity in MAS, and its class (" +
            polymerClass->name() + ") has no sourced data point to estimate it from: no MAS insulation material of "
            "the class carries a sourced value and MKF cites no datasheet of the class.");
    }
    std::vector<double> values;
    for (const auto& point : points) {
        values.push_back(point.second);
    }
    estimate.numberDocuments = documents.size();
    estimate.numberDataPoints = values.size();
    estimate.method = "the median of its " + std::to_string(values.size()) + " sourced data points";
    BobbinThermalConductivity resolved{median(values), BobbinThermalConductivityTier::CLASS, std::nullopt, estimate};

    std::ostringstream message;
    message << "Bobbin material '" << material.get_name() << "' has no thermalConductivity in MAS: ESTIMATED "
            << resolved.value << " W/(m K) from its polymer class (" << estimate.polymerClass << "): "
            << estimate.sourceMaterials.size() << " sourced MAS materials, " << estimate.datasheetGrades.size()
            << " cited datasheet grades, " << estimate.numberDocuments << " documents, using " << estimate.method << ".";
    note_once_bobbin(material, "estimated", message.str());
    return resolved;
}

double ThermalResistance::getInsulationMaterialThermalConductivity(const InsulationMaterial& material) {
    auto thermalCond = material.get_thermal_conductivity();
    if (thermalCond) {
        return *thermalCond;
    }
    throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
        "Insulation material '" + material.get_name() + "' has no thermal_conductivity defined");
}

} // namespace OpenMagnetics

