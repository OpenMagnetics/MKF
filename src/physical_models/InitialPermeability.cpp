#include "physical_models/InitialPermeability.h"
#include "physical_models/AmplitudePermeability.h"

#include "support/Utils.h"
#include "support/NumberFormat.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "spline.h"
#include <numbers>
#include <mutex>
#include <list>
#include <unordered_map>
#include <memory>
#include <map>
#include <sstream>
#include <streambuf>
#include <algorithm>
#include <vector>
#include <cfloat>
#include <magic_enum.hpp>
#include "support/Exceptions.h"
#include "support/Logger.h"



namespace OpenMagnetics {

double InitialPermeability::get_initial_permeability(std::string coreMaterialName,
                                                     std::optional<double> temperature,
                                                     std::optional<double> magneticFieldDcBias,
                                                     std::optional<double> frequency,
                                                     std::optional<double> magneticFluxDensity,
                                                     std::optional<CoreShapeFamily> shapeFamily) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_initial_permeability(coreMaterial, temperature, magneticFieldDcBias, frequency, magneticFluxDensity, shapeFamily);
}


// ABT #358: resolve the modifier to use for a given shape family. See the header for the
// key grammar; family entries are partial, so this returns "default" with the family's
// factors overlaid.
InitialPermeabilitModifier InitialPermeability::resolve_modifier(const PermeabilityPoint& permeabilityPoint,
                                                                std::optional<CoreShapeFamily> shapeFamily) {
    if (!permeabilityPoint.get_modifiers()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Permeability point has no modifiers to resolve");
    }
    auto modifiers = permeabilityPoint.get_modifiers().value();
    if (!modifiers.contains("default")) {
        // Every catalogued material carries "default"; without it there is no base fit to
        // overlay onto, and silently default-constructing one would ship mu = 1 physics.
        throw InvalidInputException(ErrorCode::MISSING_DATA,
                                   "Permeability modifiers carry no \"default\" entry for material fit");
    }
    InitialPermeabilitModifier resolved = modifiers["default"];
    if (!shapeFamily) {
        return resolved;
    }

    std::string familyName = std::string(magic_enum::enum_name(shapeFamily.value()));
    std::transform(familyName.begin(), familyName.end(), familyName.begin(), ::toupper);

    for (auto const& [key, familyModifier] : modifiers) {
        if (key == "default") {
            continue;
        }
        bool matches = false;
        std::stringstream keyStream(key);
        std::string token;
        while (std::getline(keyStream, token, '/')) {
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            std::transform(token.begin(), token.end(), token.begin(), ::toupper);
            if (token == familyName) {
                matches = true;
                break;
            }
        }
        if (!matches) {
            continue;
        }
        // Overlay only what the family entry actually provides.
        if (familyModifier.get_method()) {
            resolved.set_method(familyModifier.get_method());
        }
        if (familyModifier.get_temperature_factor()) {
            resolved.set_temperature_factor(familyModifier.get_temperature_factor());
        }
        if (familyModifier.get_frequency_factor()) {
            resolved.set_frequency_factor(familyModifier.get_frequency_factor());
        }
        if (familyModifier.get_magnetic_field_dc_bias_factor()) {
            resolved.set_magnetic_field_dc_bias_factor(familyModifier.get_magnetic_field_dc_bias_factor());
        }
        if (familyModifier.get_magnetic_flux_density_factor()) {
            resolved.set_magnetic_flux_density_factor(familyModifier.get_magnetic_flux_density_factor());
        }
        break;
    }
    return resolved;
}


double InitialPermeability::get_initial_permeability(const CoreMaterial& coreMaterial, OperatingPoint operatingPoint,
                                                     std::optional<CoreShapeFamily> shapeFamily) {
    if (operatingPoint.get_excitations_per_winding().size() == 0) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Operating point is missing excitations");
    }
    double temperature = operatingPoint.get_conditions().get_ambient_temperature();
    auto frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();

    std::optional<double> magneticFieldDcBias = std::nullopt;
    if (operatingPoint.get_excitations_per_winding()[0].get_magnetic_field_strength()) {
        if (operatingPoint.get_excitations_per_winding()[0].get_magnetic_field_strength()->get_processed()) {
            magneticFieldDcBias = operatingPoint.get_excitations_per_winding()[0].get_magnetic_field_strength()->get_processed()->get_offset();
        }
    }
    std::optional<double> magneticFluxDensity = std::nullopt;
    if (operatingPoint.get_excitations_per_winding()[0].get_magnetic_flux_density()) {
        if (operatingPoint.get_excitations_per_winding()[0].get_magnetic_flux_density()->get_processed()) {
            magneticFluxDensity = operatingPoint.get_excitations_per_winding()[0].get_magnetic_flux_density()->get_processed()->get_peak();
        }
    }
    return get_initial_permeability(coreMaterial, temperature, magneticFieldDcBias, frequency, magneticFluxDensity, shapeFamily);
}

double InitialPermeability::get_initial_permeability(std::string coreMaterialName, OperatingPoint operatingPoint,
                                                     std::optional<CoreShapeFamily> shapeFamily) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_initial_permeability(coreMaterial, operatingPoint, shapeFamily);
}

double InitialPermeability::has_frequency_dependency(const CoreMaterial& coreMaterial) {
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        const auto& permeabilityPoint = std::get<PermeabilityPoint>(initialPermeabilityData);
        if (permeabilityPoint.get_modifiers()) {
            InitialPermeabilitModifier modifiers = (*permeabilityPoint.get_modifiers())["default"];
            if (modifiers.get_frequency_factor()) {
                return true;
            }
            else {
                return false;
            }
        }
        return false;
    }
    else {
        const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
        for (const auto& point : permeabilityPoints) {
            if (point.get_frequency()) {
                if (point.get_frequency().value() > constants.quasiStaticFrequencyLimit) {
                    return true;
                }
            }
        }
        return false;
    }
}

double InitialPermeability::has_temperature_dependency(const CoreMaterial& coreMaterial) {
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    double firstTemperatureValue = DBL_MAX;
    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        return false;
    }
    else {
        const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
        for (const auto& point : permeabilityPoints) {
            if (point.get_temperature()) {
                if (firstTemperatureValue == DBL_MAX) {
                    firstTemperatureValue = point.get_temperature().value();
                }
                else if (fabs(firstTemperatureValue - point.get_temperature().value()) > 1e-9) {
                    return true;
                }
            }
        }
        return false;
    }
}

// The material's magnetisation curve B(H) = mu0 * integral_0^H mu_rev(h) dh, tabulated once per material,
// temperature and frequency. mu_rev is floored at 1: the vacuum term B = mu0 (H + M) never vanishes, and the
// datasheet fits (1/(a + b H^c) and kin) fall below it when extrapolated. The grid is geometric from 1e-6 of the
// saturation field strength to 1000 times it (0.8 % steps, trapezoid rule), so it resolves the low-field slope and
// the knee alike; a linear march on Hsat / 400 steps took up to 400 000 permeability evaluations per call, which
// made every biased inductance in an adviser cost seconds (ABT #1223).
namespace {
// Least-recently-used store of magnetisation curves. A curve is rebuilt from the same material data whenever it is
// asked for again after an eviction, so the bound changes memory, never a result.
struct MagnetisationCurveCache {
    std::mutex mutex;
    std::list<std::string> recency;  // front = most recently used
    std::unordered_map<std::string, std::pair<std::shared_ptr<const InitialPermeability::MagnetisationCurve>, std::list<std::string>::iterator>> curves;
};
MagnetisationCurveCache& magnetisation_curve_cache() {
    static MagnetisationCurveCache cache;
    return cache;
}
}

size_t InitialPermeability::get_magnetisation_curve_cache_size() {
    auto& cache = magnetisation_curve_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    return cache.curves.size();
}

std::shared_ptr<const InitialPermeability::MagnetisationCurve> InitialPermeability::get_magnetisation_curve(const CoreMaterial& coreMaterial, double temperature, std::optional<double> frequency) {
    // Keyed on the data the curve is built from, not the name alone (a custom material may reuse a catalogue
    // name): the initial permeability with its modifiers, and the saturation points that set Hsat.
    json initialPermeabilityJson = std::visit([](const auto& initial) { return json(initial); }, coreMaterial.get_permeability().get_initial());
    json saturationJson = json(coreMaterial.get_saturation());
    std::ostringstream keyStream;
    keyStream.precision(17);
    keyStream << coreMaterial.get_name() << '|' << initialPermeabilityJson.dump() << '|' << saturationJson.dump() << '|' << temperature << '|' << (frequency ? frequency.value() : -1.0);
    std::string key = keyStream.str();

    auto& cache = magnetisation_curve_cache();
    {
        std::lock_guard<std::mutex> lock(cache.mutex);
        auto found = cache.curves.find(key);
        if (found != cache.curves.end()) {
            cache.recency.splice(cache.recency.begin(), cache.recency, found->second.second);
            return found->second.first;
        }
    }

    InitialPermeability initialPermeability;
    double saturationFieldStrength = Core::get_magnetic_field_strength_saturation(coreMaterial, temperature);
    if (!(saturationFieldStrength > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Material " + coreMaterial.get_name() + " has no positive saturation field strength at " + std::to_string(temperature) + " C to build its magnetisation curve");
    }
    double vacuumPermeability = constants.vacuumPermeability;
    auto reversiblePermeabilityAt = [&](double fieldStrength) {
        return std::max(1.0, initialPermeability.get_initial_permeability_model_value(coreMaterial, temperature, fieldStrength, frequency));
    };
    auto curve = std::make_shared<MagnetisationCurve>();
    curve->materialName = coreMaterial.get_name();
    curve->saturationFluxDensity = Core::get_magnetic_flux_density_saturation(coreMaterial, temperature, false);
    curve->paramagnetic = is_paramagnetic(coreMaterial, temperature);
    const double firstFieldStrength = 1e-6 * saturationFieldStrength;
    const double lastFieldStrength = 1000 * saturationFieldStrength;
    const double ratio = 1.008;
    // From 0 to the first grid point mu_rev is its zero-field value.
    double previousPermeability = reversiblePermeabilityAt(0);
    curve->fieldStrength.push_back(0);
    curve->fluxDensity.push_back(0);
    double fieldStrength = firstFieldStrength;
    double permeability = reversiblePermeabilityAt(fieldStrength);
    curve->fieldStrength.push_back(fieldStrength);
    curve->fluxDensity.push_back(vacuumPermeability * 0.5 * (previousPermeability + permeability) * fieldStrength);
    previousPermeability = permeability;
    while (fieldStrength < lastFieldStrength) {
        double nextFieldStrength = fieldStrength * ratio;
        double nextPermeability = reversiblePermeabilityAt(nextFieldStrength);
        curve->fluxDensity.push_back(curve->fluxDensity.back() + vacuumPermeability * 0.5 * (previousPermeability + nextPermeability) * (nextFieldStrength - fieldStrength));
        curve->fieldStrength.push_back(nextFieldStrength);
        fieldStrength = nextFieldStrength;
        previousPermeability = nextPermeability;
    }

    std::lock_guard<std::mutex> lock(cache.mutex);
    auto found = cache.curves.find(key);
    if (found != cache.curves.end()) {
        // Another thread built the same curve meanwhile: both are the same numbers, keep the stored one.
        cache.recency.splice(cache.recency.begin(), cache.recency, found->second.second);
        return found->second.first;
    }
    cache.recency.push_front(key);
    std::shared_ptr<const MagnetisationCurve> stored = std::move(curve);
    cache.curves.emplace(std::move(key), std::make_pair(stored, cache.recency.begin()));
    while (cache.curves.size() > magnetisationCurveCacheCapacity) {
        cache.curves.erase(cache.recency.back());
        cache.recency.pop_back();
    }
    return stored;
}

double InitialPermeability::get_magnetic_field_dc_bias_for_flux_density(const CoreMaterial& coreMaterial, double biasFluxDensity, double temperature, std::optional<double> frequency, bool refuseAboveSaturation) {
    if (biasFluxDensity < 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "A DC bias flux density cannot be negative: " + std::to_string(biasFluxDensity) + " T");
    }
    if (biasFluxDensity == 0) {
        return 0;
    }
    double saturationFluxDensity = Core::get_magnetic_flux_density_saturation(coreMaterial, temperature, false);
    // Sizing a gap (refuseAboveSaturation) asks for a bias no gap can hold. Analysing a given design
    // does not: its operating point may sit past the knee, and the curve carries on there on the
    // vacuum slope, which is the material's B(H) beyond saturation.
    // A saturation verdict on the design, not malformed input: callers that rank candidates reject this
    // one and go on, while a genuine INVALID_INPUT still stops them (ABT #1534).
    if (refuseAboveSaturation && biasFluxDensity > saturationFluxDensity) {
        throw SaturationException(
            "The DC bias puts " + std::to_string(biasFluxDensity) + " T in material " + coreMaterial.get_name() +
            ", above its saturation of " + std::to_string(saturationFluxDensity) + " T at " + std::to_string(temperature) +
            " °C: no gap holds the target inductance at that bias, fewer turns or a larger core are needed");
    }

    auto curve = get_magnetisation_curve(coreMaterial, temperature, frequency);
    const auto& fluxDensities = curve->fluxDensity;
    const auto& fieldStrengths = curve->fieldStrength;
    if (biasFluxDensity <= fluxDensities.back()) {
        // B(H) is strictly increasing (mu_rev >= 1), so the crossing is found by bisection and interpolated linearly.
        size_t upper = std::lower_bound(fluxDensities.begin(), fluxDensities.end(), biasFluxDensity) - fluxDensities.begin();
        size_t lower = upper - 1;
        return fieldStrengths[lower] + (biasFluxDensity - fluxDensities[lower]) / (fluxDensities[upper] - fluxDensities[lower]) * (fieldStrengths[upper] - fieldStrengths[lower]);
    }
    if (!refuseAboveSaturation && biasFluxDensity > saturationFluxDensity) {
        // Analysing a design driven past saturation: a thousand times the saturation field strength
        // is far beyond the knee, the magnetisation no longer grows, and B rises on the vacuum slope.
        return fieldStrengths.back() + (biasFluxDensity - fluxDensities.back()) / constants.vacuumPermeability;
    }
    throw std::runtime_error("The DC-bias permeability curve of material " + coreMaterial.get_name() + " integrates to only " +
                             std::to_string(fluxDensities.back()) + " T at " + std::to_string(fieldStrengths.back()) + " A/m, below the " +
                             std::to_string(biasFluxDensity) + " T asked, although that is under its saturation of " +
                             std::to_string(saturationFluxDensity) + " T: the material's DC-bias factor is inconsistent with its saturation data");
}

double InitialPermeability::get_flux_density_for_magnetic_field_dc_bias(const CoreMaterial& coreMaterial, double magneticFieldDcBias, double temperature, std::optional<double> frequency) {
    if (!(magneticFieldDcBias >= 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "A DC bias field strength must be a non-negative number: " + std::to_string(magneticFieldDcBias) + " A/m");
    }
    // Ampere with no gap on a 1 m path: the MMF is the field strength itself, and the crossing is the curve's B(H).
    return get_dc_operating_point_for_magnetomotive_force(coreMaterial, magneticFieldDcBias, 1.0, 1.0, 0.0, temperature, frequency).second;
}

std::pair<double, double> InitialPermeability::get_dc_operating_point_for_magnetomotive_force(const CoreMaterial& coreMaterial, double magnetomotiveForce, double effectiveLength, double effectiveArea, double gapReluctance, double temperature, std::optional<double> frequency) {
    if (magnetomotiveForce < 0 || effectiveLength <= 0 || effectiveArea <= 0 || gapReluctance < 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "DC-bias load line needs a non-negative MMF and gap reluctance and positive core dimensions: MMF " + std::to_string(magnetomotiveForce) +
            " A, l_e " + std::to_string(effectiveLength) + " m, A_e " + std::to_string(effectiveArea) + " m2, R_gap " + std::to_string(gapReluctance) + " 1/H");
    }
    if (magnetomotiveForce == 0) {
        return {0, 0};
    }
    // Ampere around the magnetic circuit with the MMF imposed: N*I = H*l_e + B(H)*A_e*R_gap on the material's
    // magnetisation curve. The right side grows strictly with H, so the crossing on the tabulated curve is found by
    // bisection; past the table's end (1000 Hsat) B rises on the vacuum slope and the crossing is linear. Without a
    // gap this is H = N*I/l_e, the vendors' own rule for powder cores.
    // The caller's share of the cached curve keeps it alive for the call (the cache is bounded and may evict it).
    auto curve = get_magnetisation_curve(coreMaterial, temperature, frequency);
    return get_dc_operating_point_for_magnetomotive_force(*curve, magnetomotiveForce, effectiveLength, effectiveArea, gapReluctance);
}

std::pair<double, double> InitialPermeability::get_dc_operating_point_for_magnetomotive_force(const MagnetisationCurve& curve, double magnetomotiveForce, double effectiveLength, double effectiveArea, double gapReluctance) {
    if (magnetomotiveForce < 0 || effectiveLength <= 0 || effectiveArea <= 0 || gapReluctance < 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "DC-bias load line needs a non-negative MMF and gap reluctance and positive core dimensions: MMF " + std::to_string(magnetomotiveForce) +
            " A, l_e " + std::to_string(effectiveLength) + " m, A_e " + std::to_string(effectiveArea) + " m2, R_gap " + std::to_string(gapReluctance) + " 1/H");
    }
    if (magnetomotiveForce == 0) {
        return {0, 0};
    }
    const auto& fieldStrengths = curve.fieldStrength;
    const auto& fluxDensities = curve.fluxDensity;
    auto loadLine = [&](size_t index) {
        return fieldStrengths[index] * effectiveLength + fluxDensities[index] * effectiveArea * gapReluctance;
    };
    size_t last = fieldStrengths.size() - 1;
    if (magnetomotiveForce > loadLine(last)) {
        double vacuumPermeability = constants.vacuumPermeability;
        double slope = effectiveLength + vacuumPermeability * effectiveArea * gapReluctance;
        double extraFieldStrength = (magnetomotiveForce - loadLine(last)) / slope;
        return {fieldStrengths[last] + extraFieldStrength, fluxDensities[last] + vacuumPermeability * extraFieldStrength};
    }
    size_t lower = 0;
    size_t upper = last;
    while (upper - lower > 1) {
        size_t middle = (lower + upper) / 2;
        if (loadLine(middle) < magnetomotiveForce) {
            lower = middle;
        }
        else {
            upper = middle;
        }
    }
    double fraction = (magnetomotiveForce - loadLine(lower)) / (loadLine(upper) - loadLine(lower));
    return {fieldStrengths[lower] + fraction * (fieldStrengths[upper] - fieldStrengths[lower]),
            fluxDensities[lower] + fraction * (fluxDensities[upper] - fluxDensities[lower])};
}

std::pair<double, double> InitialPermeability::get_magnetomotive_force_for_flux_density(const MagnetisationCurve& curve, double fluxDensity, double effectiveLength, double effectiveArea, double gapReluctance) {
    if (!(fluxDensity >= 0) || effectiveLength <= 0 || effectiveArea <= 0 || gapReluctance < 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Inverse DC-bias load line needs a non-negative flux density and gap reluctance and positive core dimensions: B " + std::to_string(fluxDensity) +
            " T, l_e " + std::to_string(effectiveLength) + " m, A_e " + std::to_string(effectiveArea) + " m2, R_gap " + std::to_string(gapReluctance) + " 1/H");
    }
    if (fluxDensity == 0) {
        return {0, 0};
    }
    if (curve.paramagnetic) {
        // mu_r = 1 at every H: the inverse is mu0's at any flux density, the table's end is only where it was
        // sampled (1000 Hsat(T), a few A/m near the Curie point) (ABT #1942, Alf 2026-10-19).
        double fieldStrength = fluxDensity / constants.vacuumPermeability;
        return {fieldStrength, fieldStrength * effectiveLength + fluxDensity * effectiveArea * gapReluctance};
    }
    const auto& fieldStrengths = curve.fieldStrength;
    const auto& fluxDensities = curve.fluxDensity;
    if (fluxDensity > fluxDensities.back()) {
        // No extrapolation past the table (1000 Hsat, Alf 2026-10-18): the curve says nothing there.
        // The numbers as ResultsAnalyzer's findings print them (ABT #2001): these messages are the CORE_SATURATES
        // finding's reason for a point past the curve at ambient (ABT #2052).
        auto number = [](double value) {
            return format_significant(value);
        };
        if (fluxDensity > curve.saturationFluxDensity) {
            throw SaturationException(
                "The flux density " + number(fluxDensity) + " T asked of material " + curve.materialName +
                " lies past the end of its magnetisation curve, " + number(fluxDensities.back()) + " T at " +
                number(fieldStrengths.back()) + " A/m (1000 times its saturation field strength), and above its B_sat of " +
                number(curve.saturationFluxDensity) + " T: the core is driven far past saturation, and its curve has no field strength for that flux");
        }
        throw MagnetisationCurveExhaustedException("The magnetisation curve of material " + curve.materialName + " integrates to only " +
                                 number(fluxDensities.back()) + " T at " + number(fieldStrengths.back()) + " A/m, below the " +
                                 number(fluxDensity) + " T asked, although that is under its saturation of " +
                                 number(curve.saturationFluxDensity) + " T: the material's DC-bias factor is inconsistent with its saturation data");
    }
    // B(H) is strictly increasing (mu_rev >= 1); the segment and its linear interpolation are the forward
    // crossing's, so the two are exact inverses on the table.
    size_t upper = std::lower_bound(fluxDensities.begin(), fluxDensities.end(), fluxDensity) - fluxDensities.begin();
    size_t lower = upper - 1;
    double fieldStrength = fieldStrengths[lower] + (fluxDensity - fluxDensities[lower]) / (fluxDensities[upper] - fluxDensities[lower]) * (fieldStrengths[upper] - fieldStrengths[lower]);
    return {fieldStrength, fieldStrength * effectiveLength + fluxDensity * effectiveArea * gapReluctance};
}

double InitialPermeability::has_magnetic_field_dc_bias_dependency(const CoreMaterial& coreMaterial) {
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        return false;
    }
    else {
        const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
        for (const auto& point : permeabilityPoints) {
            if (point.get_magnetic_field_dc_bias()) {
                if (point.get_magnetic_field_dc_bias().value() > 0) {
                    return true;
                }
            }
        }
        return false;
    }
}

std::map<std::string, std::string> InitialPermeability::get_initial_permeability_equations(PermeabilityPoint permeabilityPoint) {
    std::map<std::string, std::string> equations;

    if (permeabilityPoint.get_modifiers()) {
        InitialPermeabilitModifier modifiers = (*permeabilityPoint.get_modifiers())["default"];
        if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MAGNETICS) {
            auto temperatureFactor = modifiers.get_temperature_factor();
            if (temperatureFactor) {
                equations["temperatureFactor"] = "1 + (a + b*T + c*T^2 + d*T^3 + e*T^4)";
            }

            auto frequencyFactor = modifiers.get_frequency_factor();
            if (frequencyFactor) {
                equations["frequencyFactor"] = "1 + (a + b*f + c*f^2 + d*f^3 + e*f^4)";
            }

            auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
            if (magneticFieldDcBiasFactor) {
                equations["magneticFieldDcBiasFactor"] = "0.01 / (a + b*(H^c))";
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MICROMETALS) {

            auto frequencyFactor = modifiers.get_frequency_factor();
            if (frequencyFactor) {
                equations["frequencyFactor"] = "(1.0 / (a + b * f^c) + d) / mu_ini";
            }

            auto temperatureFactor = modifiers.get_temperature_factor();
            if (temperatureFactor) {
                if (temperatureFactor->get_b()) {
                    equations["temperatureFactor"] = "1 + (a + c * T + e * T^2) / (1 + b * T + d * T^2) * 0.01";
                }
                else {
                    equations["temperatureFactor"] = "1 + (a * (T - 20) * 0.0001) * 0.01";
                }
            }

            auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
            if (magneticFieldDcBiasFactor) {
                equations["magneticFieldDcBiasFactor"] = "(1.0 / (a + b * H^c) + d) * 0.01";
            }

            auto magneticFluxDensityFactor = modifiers.get_magnetic_flux_density_factor();
            if (magneticFluxDensityFactor) {
                equations["magneticFluxDensityFactor"] = "(1.0 / (1.0 / (a + b * B^c) + 1.0 / (d * B^e) + 1.0 / f)) * 0.01";
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::POCO) {
            auto frequencyFactor = modifiers.get_frequency_factor();
            if (frequencyFactor) {
                equations["frequencyFactor"] = "(a / (1 + (f / b)^c) + d) * 0.01";
            }

            auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
            if (magneticFieldDcBiasFactor) {
                equations["magneticFieldDcBiasFactor"] = "a / (1 + (H / 79.57747 / b)^c) + d";
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::TDG) {
            auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
            if (magneticFieldDcBiasFactor) {
                equations["magneticFieldDcBiasFactor"] = "1 / (a + b * H^c)";
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::FAIR_RITE) {
            auto temperatureFactor = modifiers.get_temperature_factor();
            if (temperatureFactor) {
                equations["temperatureFactor"] = "(1 + a * T * 0.01)";
            }

        }
    }
    return equations;

}

double InitialPermeability::get_initial_permeability_formula(const CoreMaterial& coreMaterial,
                                                             std::optional<double> temperature,
                                                             std::optional<double> magneticFieldDcBias,
                                                             std::optional<double> frequency,
                                                             std::optional<double> magneticFluxDensity,
                                                             std::optional<CoreShapeFamily> shapeFamily) {
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    const auto& permeabilityPoint = std::get<PermeabilityPoint>(initialPermeabilityData);
    double initialPermeabilityValue = permeabilityPoint.get_value();

    if (permeabilityPoint.get_modifiers()) {
        InitialPermeabilitModifier modifiers = resolve_modifier(permeabilityPoint, shapeFamily);
        if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MAGNETICS) {
            auto temperatureFactor = modifiers.get_temperature_factor();
            if (temperature && temperatureFactor) {
                double permeabilityVariationDueToTemperature =
                    temperatureFactor->get_a() + temperatureFactor->get_b().value() * (temperature.value()) +
                    temperatureFactor->get_c().value() * pow(temperature.value(), 2) +
                    temperatureFactor->get_d().value() * pow(temperature.value(), 3) +
                    temperatureFactor->get_e().value() * pow(temperature.value(), 4);
                initialPermeabilityValue *= (1 + permeabilityVariationDueToTemperature);
            }

            auto frequencyFactor = modifiers.get_frequency_factor();
            if (frequency && frequencyFactor) {
                double permeabilityVariationDueToFrequency =
                    frequencyFactor->get_a() + frequencyFactor->get_b() * (frequency.value()) +
                    frequencyFactor->get_c() * pow(frequency.value(), 2) + frequencyFactor->get_d() * pow(frequency.value(), 3) +
                    frequencyFactor->get_e().value() * pow(frequency.value(), 4);

                initialPermeabilityValue *= (1 + permeabilityVariationDueToFrequency);
            }

            if (magneticFieldDcBias && modifiers.get_magnetic_field_dc_bias_factor()) {  // factor optional checked like the temperature/frequency siblings (deref of a disengaged optional is UB)
                auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
                double a = magneticFieldDcBiasFactor->get_a();
                double b = magneticFieldDcBiasFactor->get_b();
                double c = magneticFieldDcBiasFactor->get_c();
                double H = magneticFieldDcBias.value();
                double permeabilityVariationDueToMagneticFieldDcBias = 0.01 / (a + b * pow(roundFloat(fabs(H), 3), c));

                initialPermeabilityValue *= permeabilityVariationDueToMagneticFieldDcBias;
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MICROMETALS) {

            auto frequencyFactor = modifiers.get_frequency_factor();
            if (frequency && frequencyFactor) {
                double a = frequencyFactor->get_a();
                double b = frequencyFactor->get_b();
                double c = frequencyFactor->get_c();
                double d = frequencyFactor->get_d();
                double f = frequency.value();
                initialPermeabilityValue = 1.0 / (a + b * pow(f, c)) + d;
            }

            auto temperatureFactor = modifiers.get_temperature_factor();
            if (temperature && temperatureFactor) {
                double permeabilityVariationDueToTemperature;
                if (temperatureFactor->get_b()) {
                    double a = temperatureFactor->get_a();
                    double b = temperatureFactor->get_b().value();
                    double c = temperatureFactor->get_c().value();
                    double d = temperatureFactor->get_d().value();
                    double e = temperatureFactor->get_e().value();
                    double T = temperature.value();
                    permeabilityVariationDueToTemperature = (a + c * T + e * pow(T, 2)) / (1 + b * T + d * pow(T, 2));
                }
                else {
                    double a = temperatureFactor->get_a();
                    double T = temperature.value();
                    permeabilityVariationDueToTemperature = a * (T - 20) * 0.0001;
                }

                initialPermeabilityValue *= (1 + permeabilityVariationDueToTemperature * 0.01);
            }

            if (magneticFieldDcBias && modifiers.get_magnetic_field_dc_bias_factor()) {  // factor optional checked like the temperature/frequency siblings (deref of a disengaged optional is UB)
                auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
                double a = magneticFieldDcBiasFactor->get_a();
                double b = magneticFieldDcBiasFactor->get_b();
                double c = magneticFieldDcBiasFactor->get_c();
                double d = magneticFieldDcBiasFactor->get_d().value();
                double H = magneticFieldDcBias.value();
                double permeabilityVariationDueToMagneticFieldDcBias = 1.0 / (a + b * pow(roundFloat(fabs(H), 3), c)) + d;

                initialPermeabilityValue *= permeabilityVariationDueToMagneticFieldDcBias * 0.01;
            }

            if (magneticFluxDensity && modifiers.get_magnetic_flux_density_factor()) {  // factor optional checked like the siblings
                auto magneticFluxDensityFactor = modifiers.get_magnetic_flux_density_factor();
                double a = magneticFluxDensityFactor->get_a();
                double b = magneticFluxDensityFactor->get_b();
                double c = magneticFluxDensityFactor->get_c();
                double d = magneticFluxDensityFactor->get_d();
                double e = magneticFluxDensityFactor->get_e();
                double f = magneticFluxDensityFactor->get_f();
                double B = magneticFluxDensity.value();
                double permeabilityVariationDueToMagneticFluxDensity = 1.0 / (1.0 / (a + b * pow(B, c)) + 1.0 / (d * pow(B, e)) + 1.0 / f);

                initialPermeabilityValue *= permeabilityVariationDueToMagneticFluxDensity * 0.01;
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::POCO) {
            // ABT #169: percent-of-initial rolloff fitted to the POCO catalog
            // "Permeability vs Frequency" curves. Same logistic shape as the
            // POCO DC-bias factor: a is the rolling-off share, d the
            // high-frequency asymptote (a + d = 100 at DC), b the corner
            // frequency in Hz and c the steepness.
            auto frequencyFactor = modifiers.get_frequency_factor();
            if (frequency && frequencyFactor) {
                double a = frequencyFactor->get_a();
                double b = frequencyFactor->get_b();
                double c = frequencyFactor->get_c();
                double d = frequencyFactor->get_d();
                double f = frequency.value();
                double permeabilityVariationDueToFrequency = a / (1 + pow(f / b, c)) + d;

                initialPermeabilityValue *= permeabilityVariationDueToFrequency * 0.01;
            }

            if (magneticFieldDcBias && modifiers.get_magnetic_field_dc_bias_factor()) {  // factor optional checked like the temperature/frequency siblings (deref of a disengaged optional is UB)
                auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
                double a = magneticFieldDcBiasFactor->get_a();
                double b = magneticFieldDcBiasFactor->get_b();
                double c = magneticFieldDcBiasFactor->get_c();
                double d = magneticFieldDcBiasFactor->get_d().value();
                double H = magneticFieldDcBias.value();
                double permeabilityVariationDueToMagneticFieldDcBias = a / (1 + pow(fabs(H) / 79.57747 / b, c)) + d;


                initialPermeabilityValue *= permeabilityVariationDueToMagneticFieldDcBias * 0.01;
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::TDG) {
            if (magneticFieldDcBias && modifiers.get_magnetic_field_dc_bias_factor()) {  // factor optional checked like the temperature/frequency siblings (deref of a disengaged optional is UB)
                auto magneticFieldDcBiasFactor = modifiers.get_magnetic_field_dc_bias_factor();
                double a = magneticFieldDcBiasFactor->get_a();
                double b = magneticFieldDcBiasFactor->get_b();
                double c = magneticFieldDcBiasFactor->get_c();
                double H = magneticFieldDcBias.value();
                double permeabilityVariationDueToMagneticFieldDcBias = 1 / (a + b * pow(fabs(H) / 79.57747, c));


                initialPermeabilityValue *= permeabilityVariationDueToMagneticFieldDcBias * 0.01;
            }
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::FAIR_RITE) {
            auto temperatureFactor = modifiers.get_temperature_factor();
            if (temperature && temperatureFactor) {
                double permeabilityVariationDueToTemperature;
                double a = temperatureFactor->get_a();
                double T = temperature.value();
                permeabilityVariationDueToTemperature = a * T;

                initialPermeabilityValue *= (1 + permeabilityVariationDueToTemperature * 0.01);
            }

        }
    }
    return initialPermeabilityValue;
}

double get_closes_temperature_to_default_in_permeability_points(const std::vector<PermeabilityPoint>& points) {
    double closestDistance = DBL_MAX;
    double closestTemperature = DBL_MAX;
    for (const auto& point : points) {
        if (point.get_temperature()) {
            double distance = fabs(point.get_temperature().value() - defaults.ambientTemperature);
            if (distance < closestDistance) {
                closestDistance = distance;
                closestTemperature = point.get_temperature().value();
            }
        }
    }

    return closestTemperature;
}

double get_minimum_frequency_in_permeability_points(const std::vector<PermeabilityPoint>& points) {
    auto order_frequency = []( const auto &e1, const auto &e2 ) {
        if (!e1.get_frequency()) {
            return false;
        }
        if (!e2.get_frequency()) {
            return true;
        }
        return e1.get_frequency().value() < e2.get_frequency().value();
    };

    auto minimumFrequencyPoint = std::min_element( std::begin( points ), std::end( points ), order_frequency );
    if (minimumFrequencyPoint->get_frequency()) {
        return minimumFrequencyPoint->get_frequency().value();
    }
    else {
        return constants.quasiStaticFrequencyLimit;
    }
}

double get_maximum_frequency_in_permeability_points(const std::vector<PermeabilityPoint>& points) {
    auto order_frequency = []( const auto &e1, const auto &e2 ) {
        if (!e1.get_frequency()) {
            return false;
        }
        if (!e2.get_frequency()) {
            return true;
        }
        return e1.get_frequency().value() < e2.get_frequency().value();
    };

    auto maximumFrequencyPoint = std::max_element( std::begin( points ), std::end( points ), order_frequency );
    if (maximumFrequencyPoint->get_frequency()) {
        return maximumFrequencyPoint->get_frequency().value();
    }
    else {
        return DBL_MAX;
    }
}

double get_minimum_magnetic_field_dc_bias_in_permeability_points(const std::vector<PermeabilityPoint>& points) {
    auto minimum_magnetic_field_dc_bias = []( const auto &e1, const auto &e2 ) {
        if (!e1.get_magnetic_field_dc_bias()) {
            return false;
        }
        if (!e2.get_magnetic_field_dc_bias()) {
            return true;
        }
        return e1.get_magnetic_field_dc_bias().value() < e2.get_magnetic_field_dc_bias().value();
    };

    auto minimumMagneticFieldDcBiasPoint = std::min_element( std::begin( points ), std::end( points ), minimum_magnetic_field_dc_bias );
    if (minimumMagneticFieldDcBiasPoint->get_magnetic_field_dc_bias()) {
        return minimumMagneticFieldDcBiasPoint->get_magnetic_field_dc_bias().value();
    }
    else {
        return 0;
    }
}

std::vector<size_t> InitialPermeability::get_only_temperature_dependent_indexes(const std::vector<PermeabilityPoint>& permeabilityPoints) {
    std::vector<size_t> temperatureIndexes;
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
    for (size_t index = 0; index < permeabilityPoints.size(); ++index) {
        const auto& point = permeabilityPoints[index];
        if (point.get_frequency()) {
            if (point.get_frequency().value() > minimumFrequency) {
                continue;
            }
        }
        if (point.get_magnetic_field_dc_bias()) {
            if (point.get_magnetic_field_dc_bias().value() > minimumMagneticFieldDcBias) {
                continue;
            }
        }
        temperatureIndexes.push_back(index);
    }
    return temperatureIndexes;
}

std::vector<size_t> InitialPermeability::get_only_temperature_dependent_indexes(const CoreMaterial& coreMaterial) {
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    return get_only_temperature_dependent_indexes(permeabilityPoints);
}

std::vector<PermeabilityPoint> InitialPermeability::get_only_temperature_dependent_points(const CoreMaterial& coreMaterial) {
    std::vector<PermeabilityPoint> temperaturePoints;
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
    for (const auto& point : permeabilityPoints) {
        if (point.get_frequency()) {
            if (point.get_frequency().value() > minimumFrequency) {
                continue;
            }
        }
        if (point.get_magnetic_field_dc_bias()) {
            if (point.get_magnetic_field_dc_bias().value() > minimumMagneticFieldDcBias) {
                continue;
            }
        }
        temperaturePoints.push_back(point);
    }
    return temperaturePoints;
}

std::vector<size_t> InitialPermeability::get_only_frequency_dependent_indexes(const std::vector<PermeabilityPoint>& permeabilityPoints) {
    std::vector<size_t> frequencyIndexes;
    double defaultTemperature = get_closes_temperature_to_default_in_permeability_points(permeabilityPoints);
    double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
    for (size_t index = 0; index < permeabilityPoints.size(); ++index) {
        const auto& point = permeabilityPoints[index];
        if (point.get_temperature()) {
            if (point.get_temperature().value() != defaultTemperature) {
                continue;
            }
        }
        if (point.get_magnetic_field_dc_bias()) {
            if (point.get_magnetic_field_dc_bias().value() > minimumMagneticFieldDcBias) {
                continue;
            }
        }
        frequencyIndexes.push_back(index);
    }
    return frequencyIndexes;
}

std::vector<size_t> InitialPermeability::get_only_frequency_dependent_indexes(const CoreMaterial& coreMaterial) {
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    return get_only_frequency_dependent_indexes(permeabilityPoints);
}


std::vector<PermeabilityPoint> InitialPermeability::sample_initial_permeability_by_frequency_modifier(PermeabilityPoint permeabilityPoint) {
    std::vector<PermeabilityPoint> frequencyPoints;
    InitialPermeabilitModifier modifiers = (*permeabilityPoint.get_modifiers())["default"];
    auto frequencies = logarithmic_spaced_array(defaults.measurementFrequency, defaults.maximumFrequency, 100);
    for (auto frequency : frequencies) {
        double initialPermeabilityValue = 1;
        if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MAGNETICS) {
            double permeabilityVariationDueToFrequency =
                modifiers.get_frequency_factor()->get_a() + modifiers.get_frequency_factor()->get_b() * (frequency) +
                modifiers.get_frequency_factor()->get_c() * pow(frequency, 2) + modifiers.get_frequency_factor()->get_d() * pow(frequency, 3) +
                modifiers.get_frequency_factor()->get_e().value() * pow(frequency, 4);
            initialPermeabilityValue = permeabilityPoint.get_value() * (1 + permeabilityVariationDueToFrequency);
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::MICROMETALS) {
            double a = modifiers.get_frequency_factor()->get_a();
            double b = modifiers.get_frequency_factor()->get_b();
            double c = modifiers.get_frequency_factor()->get_c();
            double d = modifiers.get_frequency_factor()->get_d();
            double f = frequency;
            initialPermeabilityValue = 1.0 / (a + b * pow(f, c)) + d;
        }
        else if ((*modifiers.get_method()) == InitialPermeabilitModifierMethod::POCO) {
            // ABT #169: same percent-rolloff form as get_initial_permeability_formula.
            double a = modifiers.get_frequency_factor()->get_a();
            double b = modifiers.get_frequency_factor()->get_b();
            double c = modifiers.get_frequency_factor()->get_c();
            double d = modifiers.get_frequency_factor()->get_d();
            double f = frequency;
            initialPermeabilityValue = permeabilityPoint.get_value() * (a / (1 + pow(f / b, c)) + d) * 0.01;
        }

        PermeabilityPoint point;
        point.set_frequency(frequency);
        point.set_value(initialPermeabilityValue);
        frequencyPoints.push_back(point);
    }

    return frequencyPoints;
}


std::vector<PermeabilityPoint> InitialPermeability::get_only_frequency_dependent_points(const CoreMaterial& coreMaterial) {
    std::vector<PermeabilityPoint> frequencyPoints;
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();

    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        const auto& permeabilityPoint = std::get<PermeabilityPoint>(initialPermeabilityData);
        if (permeabilityPoint.get_modifiers()) {
            InitialPermeabilitModifier modifiers = (*permeabilityPoint.get_modifiers())["default"];
            if (modifiers.get_frequency_factor()) {
                frequencyPoints = sample_initial_permeability_by_frequency_modifier(permeabilityPoint);
            }
        }
    }
    else {
        const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
        double defaultTemperature = get_closes_temperature_to_default_in_permeability_points(permeabilityPoints);
        double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
        for (const auto& point : permeabilityPoints) {
            if (point.get_temperature()) {
                if (point.get_temperature().value() != defaultTemperature) {
                    continue;
                }
            }
            if (point.get_magnetic_field_dc_bias()) {
                if (point.get_magnetic_field_dc_bias().value() > minimumMagneticFieldDcBias) {
                    continue;
                }
            }
            frequencyPoints.push_back(point);
        }
    }

    return frequencyPoints;
}

std::vector<size_t> InitialPermeability::get_only_magnetic_field_dc_bias_dependent_indexes(const std::vector<PermeabilityPoint>& permeabilityPoints) {
    std::vector<size_t>  magneticFieldDcBiasIndexes;
    double defaultTemperature = get_closes_temperature_to_default_in_permeability_points(permeabilityPoints);
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    for (size_t index = 0; index < permeabilityPoints.size(); ++index) {
        const auto& point = permeabilityPoints[index];
        if (point.get_temperature()) {
            if (point.get_temperature().value() != defaultTemperature) {
                continue;
            }
        }
        if (point.get_frequency()) {
            if (point.get_frequency().value() > minimumFrequency) {
                continue;
            }
        }
         magneticFieldDcBiasIndexes.push_back(index);
    }
    return  magneticFieldDcBiasIndexes;
}


std::vector<size_t> InitialPermeability::get_only_magnetic_field_dc_bias_dependent_indexes(const CoreMaterial& coreMaterial) {
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    return get_only_magnetic_field_dc_bias_dependent_indexes(permeabilityPoints);
}

std::vector<PermeabilityPoint> InitialPermeability::get_only_magnetic_field_dc_bias_dependent_points(const CoreMaterial& coreMaterial) {
    std::vector<PermeabilityPoint>  magneticFieldDcBiasPoints;
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    double defaultTemperature = get_closes_temperature_to_default_in_permeability_points(permeabilityPoints);
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    for (const auto& point : permeabilityPoints) {
        if (point.get_temperature()) {
            if (point.get_temperature().value() != defaultTemperature) {
                continue;
            }
        }
        if (point.get_frequency()) {
            if (point.get_frequency().value() > minimumFrequency) {
                continue;
            }
        }
         magneticFieldDcBiasPoints.push_back(point);
    }
    return  magneticFieldDcBiasPoints;
}

std::optional<double> InitialPermeability::get_paramagnetic_temperature(const CoreMaterial& coreMaterial) {
    // ABT #1739: the temperature above which the material is paramagnetic (mu_i = 1): the stored Curie
    // temperature, or the last point of the measured mu_i(T) table if that runs above it (a Curie
    // temperature stored as the data sheet's guaranteed minimum). No Curie temperature: none.
    auto curieTemperature = coreMaterial.get_curie_temperature();
    if (!curieTemperature) {
        return std::nullopt;
    }
    double paramagneticTemperature = curieTemperature.value();
    auto tableEnd = get_temperature_table_end(coreMaterial);
    if (tableEnd) {
        paramagneticTemperature = std::max(paramagneticTemperature, tableEnd.value());
    }
    return paramagneticTemperature;
}

bool InitialPermeability::is_paramagnetic(const CoreMaterial& coreMaterial, double temperature) {
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        const auto& point = std::get<PermeabilityPoint>(initialPermeabilityData);
        if (point.get_value() == 1 && !point.get_modifiers()) {
            return true;
        }
    }
    // As get_initial_permeability_model_value: mu_i = 1 above the paramagnetic temperature.
    auto curieTemperature = coreMaterial.get_curie_temperature();
    if (!curieTemperature || temperature <= curieTemperature.value()) {
        return false;
    }
    auto paramagneticTemperature = get_paramagnetic_temperature(coreMaterial);
    return paramagneticTemperature && temperature > paramagneticTemperature.value();
}

std::optional<double> InitialPermeability::get_temperature_table_end(const CoreMaterial& coreMaterial) {
    if (!has_temperature_dependency(coreMaterial)) {
        return std::nullopt;
    }
    std::optional<double> tableEnd;
    // By index into the material's points: this runs on every permeability evaluation (each of the ~2600 grid
    // points of a magnetisation curve), where copying the points was most of the cost.
    const auto& allPoints = std::get<std::vector<PermeabilityPoint>>(coreMaterial.get_permeability().get_initial());
    for (size_t index : get_only_temperature_dependent_indexes(allPoints)) {
        const auto& point = allPoints[index];
        if (point.get_temperature() && (!tableEnd || point.get_temperature().value() > tableEnd.value())) {
            tableEnd = point.get_temperature().value();
        }
    }
    return tableEnd;
}

std::optional<double> InitialPermeability::get_maximum_ferromagnetic_data_temperature(const CoreMaterial& coreMaterial) {
    // ABT #838: get_initial_permeability_temperature_dependent answers up to the table's last point and
    // above get_paramagnetic_temperature (mu_i = 1); between the two it refuses. The ferromagnetic side
    // ends at the Curie temperature, so the range where the material is ferromagnetic and mu_i(T) answers
    // ends at the lower of the two.
    auto curieTemperature = coreMaterial.get_curie_temperature();
    auto tableEnd = get_temperature_table_end(coreMaterial);
    if (curieTemperature && tableEnd) {
        return std::min(curieTemperature.value(), tableEnd.value());
    }
    if (curieTemperature) {
        return curieTemperature;
    }
    return tableEnd;
}

namespace {
void combine_permeability_hash(std::size_t& seed, std::size_t value) {
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

void combine_permeability_optional_double(std::size_t& seed, const std::optional<double>& value) {
    combine_permeability_hash(seed, value ? std::hash<double>{}(value.value()) : 0x51ed270b);
}

void combine_permeability_point(std::size_t& seed, const PermeabilityPoint& point) {
    combine_permeability_hash(seed, std::hash<double>{}(point.get_value()));
    combine_permeability_optional_double(seed, point.get_frequency());
    combine_permeability_optional_double(seed, point.get_temperature());
    combine_permeability_optional_double(seed, point.get_magnetic_field_dc_bias());
    combine_permeability_optional_double(seed, point.get_magnetic_field_peak());
    combine_permeability_optional_double(seed, point.get_magnetic_flux_density_peak());
    combine_permeability_optional_double(seed, point.get_tolerance());
    if (point.get_modifiers()) {
        // the frequency/temperature/DC-bias formulas of a formula material: small, and rare
        combine_permeability_hash(seed, std::hash<std::string>{}(json(point.get_modifiers().value()).dump()));
    }
}
} // namespace

std::size_t InitialPermeability::hash_permeability(const Permeability& permeability) {
    std::size_t seed = 0;
    if (std::holds_alternative<PermeabilityPoint>(permeability)) {
        combine_permeability_point(seed, std::get<PermeabilityPoint>(permeability));
        return seed;
    }
    const auto& points = std::get<std::vector<PermeabilityPoint>>(permeability);
    combine_permeability_hash(seed, points.size());
    for (const auto& point : points) {
        combine_permeability_point(seed, point);
    }
    return seed;
}

std::string InitialPermeability::interpolator_cache_key(const CoreMaterial& coreMaterial) {
    return coreMaterial.get_name() + "#" + std::to_string(hash_permeability(coreMaterial.get_permeability().get_initial()));
}

double InitialPermeability::get_initial_permeability_temperature_dependent(const CoreMaterial& coreMaterial, double temperature) {
    const std::string cacheKey = interpolator_cache_key(coreMaterial);
    double initialPermeabilityValue = 1;
    // The points are read by index; they are copied only to build the interpolator, once per material.
    const auto& allPoints = std::get<std::vector<PermeabilityPoint>>(coreMaterial.get_permeability().get_initial());
    const auto temperatureIndexes = get_only_temperature_dependent_indexes(allPoints);
    if (temperatureIndexes.size() == 0) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "No temperature dependent points for material: " + coreMaterial.get_name());
    }

    if (!initialPermeabilityTemperatureInterps.contains(cacheKey)) {
        auto permeabilityPoints = get_only_temperature_dependent_points(coreMaterial);
        // Sort by temperature so spline x values are strictly increasing
        std::sort(permeabilityPoints.begin(), permeabilityPoints.end(),
            [](const auto& a, const auto& b) {
                return *a.get_temperature() < *b.get_temperature();
            });

        int n = permeabilityPoints.size();
        std::vector<double> x, y;

        for (int i = 0; i < n; i++) {
            if (x.size() == 0 || (*permeabilityPoints[i].get_temperature()) > x.back()) {
                x.push_back(*permeabilityPoints[i].get_temperature());
                y.push_back(permeabilityPoints[i].get_value());
            }
        }

        if (x.size() >= 3) {
            tk::spline interp(x, y, tk::spline::cspline_hermite);
            initialPermeabilityTemperatureInterps[cacheKey] =
                [interp](double t) { return interp(t); };
        }
        else if (x.size() == 2) {
            // tk::spline requires >= 3 points. Do manual linear interpolation.
            double x0 = x[0], x1 = x[1], y0 = y[0], y1 = y[1];
            initialPermeabilityTemperatureInterps[cacheKey] =
                [x0, x1, y0, y1](double t) {
                    return y0 + (y1 - y0) * (t - x0) / (x1 - x0);
                };
        }
        else {
            double c = y[0];
            initialPermeabilityTemperatureInterps[cacheKey] =
                [c](double) { return c; };
        }
    }

    auto byTemperature = [&](size_t a, size_t b) { return *allPoints[a].get_temperature() < *allPoints[b].get_temperature(); };
    double tabulatedMinimum = *allPoints[*std::min_element(temperatureIndexes.begin(), temperatureIndexes.end(), byTemperature)].get_temperature();
    double tabulatedMaximum = *allPoints[*std::max_element(temperatureIndexes.begin(), temperatureIndexes.end(), byTemperature)].get_temperature();

    // Above its Curie point a ferrite is paramagnetic: mu_i = 1 is the physics, not an extrapolation.
    // The stored Curie temperature can be a guaranteed minimum (Ferroxcube prints "TC >= 220 C") with
    // the measured mu_i(T) running above it (3C90 to 240 C), so mu_i = 1 starts only above both.
    auto curieTemperature = coreMaterial.get_curie_temperature();
    // Never below the Curie temperature: its table scan is needed only above the Curie point.
    if (curieTemperature && temperature > curieTemperature.value()) {
        auto paramagneticTemperature = get_paramagnetic_temperature(coreMaterial);
        if (paramagneticTemperature && temperature > paramagneticTemperature.value()) {
            return 1;
        }
    }

    // ABT #1739: below the Curie point, outside the tabulated mu_i(T) range, the material has no data.
    // The spline through the table is not a model of the material there, so it is not evaluated
    // unless Settings::allowMaterialDataExtrapolation asks for it (then with a WARNING).
    if (temperature < tabulatedMinimum || temperature > tabulatedMaximum) {
        MaterialTemperatureOutOfRangeException outOfRange(coreMaterial.get_name(), temperature, tabulatedMinimum,
                                                          tabulatedMaximum, curieTemperature);
        if (!material_data_extrapolation_allowed(outOfRange.what())) {
            throw outOfRange;
        }
    }

    initialPermeabilityValue = initialPermeabilityTemperatureInterps[cacheKey](temperature);

    // A ferromagnetic material below its Curie point has mu_i > 1. A value below 1 means the table (or its
    // extrapolation) says something unphysical there; it is reported, not floored to 1.
    if (initialPermeabilityValue < 1) {
        throw MaterialException(ErrorCode::MATERIAL_INVALID_PROPERTY,
                                "Material " + coreMaterial.get_name() + ": initial permeability at " +
                                std::to_string(temperature) + " C evaluates to " + std::to_string(initialPermeabilityValue) +
                                ", below 1, although the temperature is below the Curie point" +
                                (curieTemperature ? " (" + std::to_string(curieTemperature.value()) + " C)" : std::string(" (none stored)")) +
                                "; its mu_i(T) table runs from " + std::to_string(tabulatedMinimum) + " C to " +
                                std::to_string(tabulatedMaximum) + " C",
                                coreMaterial.get_name());
    }

    if (std::isnan(initialPermeabilityValue)) {
        throw NaNResultException("Initial Permeability must be a number, not NaN");
    }

    return initialPermeabilityValue;
}

double InitialPermeability::get_initial_permeability_frequency_dependent(const CoreMaterial& coreMaterial, double frequency) {
    const std::string cacheKey = interpolator_cache_key(coreMaterial);
    double initialPermeabilityValue = 1;
    // The points are copied only to build the interpolator, once per material: an interpolator is stored only
    // for a material that has them.
    if (!initialPermeabilityFrequencyInterps.contains(cacheKey)) {
        auto permeabilityPoints = get_only_frequency_dependent_points(coreMaterial);
        if (permeabilityPoints.size() == 0) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "No frequency dependent points for material: " + coreMaterial.get_name());
        }
        // Sort by frequency so spline x values are strictly increasing (the temperature
        // sibling sorts; this one relied on DB order — tk::spline asserts/UB on unsorted input).
        std::sort(permeabilityPoints.begin(), permeabilityPoints.end(),
            [](const auto& a, const auto& b) {
                return *a.get_frequency() < *b.get_frequency();
            });

        int n = permeabilityPoints.size();
        std::vector<double> x, y;

        for (int i = 0; i < n; i++) {
            if (x.size() == 0 || (*permeabilityPoints[i].get_frequency()) != x.back()) {
                x.push_back(*permeabilityPoints[i].get_frequency());
                y.push_back(permeabilityPoints[i].get_value());
            }
        }

        if (x.size() > 1) {
            tk::spline interp(x, y, tk::spline::cspline_hermite);
            initialPermeabilityFrequencyInterps[cacheKey] = interp;
        }
        else {
            initialPermeabilityFrequencyInterps[cacheKey] = permeabilityPoints[0].get_value();
        }
    }

    if (std::holds_alternative<double>(initialPermeabilityFrequencyInterps[cacheKey])) {
        initialPermeabilityValue = std::get<double>(initialPermeabilityFrequencyInterps[cacheKey]);
    }
    else {
        auto value = std::get<tk::spline>(initialPermeabilityFrequencyInterps[cacheKey])(frequency);
        initialPermeabilityValue = std::max(1., value);
    }


    if (std::isnan(initialPermeabilityValue)) {
        throw NaNResultException("Initial Permeability must be a number, not NaN");
    }

    return initialPermeabilityValue;
}

double InitialPermeability::get_initial_permeability_magnetic_field_dc_bias_dependent(const CoreMaterial& coreMaterial, double magneticFieldDcBias) {
    const std::string cacheKey = interpolator_cache_key(coreMaterial);
    double initialPermeabilityValue = 1;
    // The points are copied only to build the interpolator, once per material: an interpolator is stored only
    // for a material that has them.
    if (!initialPermeabilityMagneticFieldDcBiasInterps.contains(cacheKey)) {
        auto permeabilityPoints = get_only_magnetic_field_dc_bias_dependent_points(coreMaterial);
        if (permeabilityPoints.size() == 0) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "No magnetic field dc bias dependent points for material: " + coreMaterial.get_name());
        }
        // Sort by H bias so spline x values are strictly increasing (see the temperature
        // sibling; unsorted DB order would be UB in tk::spline).
        std::sort(permeabilityPoints.begin(), permeabilityPoints.end(),
            [](const auto& a, const auto& b) {
                return *a.get_magnetic_field_dc_bias() < *b.get_magnetic_field_dc_bias();
            });

        int n = permeabilityPoints.size();
        std::vector<double> x, y;

        for (int i = 0; i < n; i++) {
            if (x.size() == 0 || (*permeabilityPoints[i].get_magnetic_field_dc_bias()) != x.back()) {
                x.push_back(*permeabilityPoints[i].get_magnetic_field_dc_bias());
                y.push_back(permeabilityPoints[i].get_value());
            }
        }

        if (x.size() > 1) {
            tk::spline interp(x, y, tk::spline::cspline_hermite);
            initialPermeabilityMagneticFieldDcBiasInterps[cacheKey] = interp;
        }
        else {
            initialPermeabilityMagneticFieldDcBiasInterps[cacheKey] = permeabilityPoints[0].get_value();
        }
    }

    if (std::holds_alternative<double>(initialPermeabilityMagneticFieldDcBiasInterps[cacheKey])) {
        initialPermeabilityValue = std::get<double>(initialPermeabilityMagneticFieldDcBiasInterps[cacheKey]);
    }
    else {
        auto value = std::get<tk::spline>(initialPermeabilityMagneticFieldDcBiasInterps[cacheKey])(magneticFieldDcBias);
        initialPermeabilityValue = std::max(1., value);
    }

    if (std::isnan(initialPermeabilityValue)) {
        throw NaNResultException("Initial Permeability must be a number, not NaN");
    }

    return initialPermeabilityValue;
}

double InitialPermeability::calculate_frequency_for_initial_permeability_drop(const CoreMaterial& coreMaterial, double percentageDrop, double maximumError) {
    bool hasFrequencyDependency = has_frequency_dependency(coreMaterial);
    if (!hasFrequencyDependency) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::vector<PermeabilityPoint> permeabilityPoints;

    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();

    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        const auto& permeabilityPoint = std::get<PermeabilityPoint>(initialPermeabilityData);
        permeabilityPoints = sample_initial_permeability_by_frequency_modifier(permeabilityPoint);
        
    }
    else {
        permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
    }
    double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
    double maximumFrequency = get_maximum_frequency_in_permeability_points(permeabilityPoints);

    double initialPermeabilityValueReference = get_initial_permeability_frequency_dependent(coreMaterial, minimumFrequency);
    double initialPermeabilityAbsolute = initialPermeabilityValueReference * (1 - percentageDrop);
    double currentFrequency = (maximumFrequency + minimumFrequency) / 2;

    size_t timeout = 100;
    while (timeout > 0) {
        auto currentInitialPermeability = get_initial_permeability_frequency_dependent(coreMaterial, currentFrequency);
        if (fabs(currentInitialPermeability - initialPermeabilityAbsolute) / initialPermeabilityAbsolute < maximumError) {
            break;
        }
        else {
            if (currentInitialPermeability > initialPermeabilityAbsolute) {
                minimumFrequency = currentFrequency;
                currentFrequency = (maximumFrequency + currentFrequency) / 2;
            }
            else {
                maximumFrequency = currentFrequency;
                currentFrequency = (currentFrequency + minimumFrequency) / 2;
            }
        }

        timeout--;
    }

    return currentFrequency;
}


double InitialPermeability::get_initial_permeability(const CoreMaterial& coreMaterial,
                                                     std::optional<double> temperature,
                                                     std::optional<double> magneticFieldDcBias,
                                                     std::optional<double> frequency,
                                                     std::optional<double> magneticFluxDensity,
                                                     std::optional<CoreShapeFamily> shapeFamily) {
    double permeability = get_initial_permeability_model_value(coreMaterial, temperature, magneticFieldDcBias, frequency, magneticFluxDensity, shapeFamily);
    // A relative permeability below the vacuum's is not a material property: it is the model extrapolated
    // past the knee, where its B-H slope tends to 0 instead of to mu0 (ABT #1740: 3C96 read 0.499 at
    // 6.2 kA/m and 8e-13 at 5e9 A/m). The material there is saturated, with the vacuum's permeability.
    if (permeability < 1) {
        // Tell the two causes apart, as they need different answers: the same material and conditions
        // without the bias. If that is physical, the bias drove the model past the knee (the material is
        // saturated there); if it is not, the frequency or temperature dependence itself is outside its
        // valid range (e.g. the Magnetics frequencyFactor quartic, which goes negative), and calling that
        // "saturation" would send the reader after the wrong thing.
        bool biased = (magneticFieldDcBias && magneticFieldDcBias.value() != 0) || (magneticFluxDensity && magneticFluxDensity.value() != 0);
        std::optional<double> unbiasedPermeability;
        if (biased) {
            unbiasedPermeability = get_initial_permeability_model_value(coreMaterial, temperature, std::nullopt, frequency, std::nullopt, shapeFamily);
        }
        // The bias case is the material saturated at that bias: past the knee its reversible permeability is
        // the vacuum's, 1, the slope its own magnetisation curve takes there (get_magnetisation_curve) and the
        // saturated state of the DC-bias sweep (ABT #1790). Until ABT #1855 this was refused
        // (SaturationException); a saturated core now gets its numbers and ResultsAnalyzer::check_saturation
        // reports the saturation. A model outside its frequency or temperature range is not saturation: it
        // has no value to give, and must stop the run.
        if (biased && unbiasedPermeability.value() >= 1) {
            return 1;
        }
        std::ostringstream message;
        message << "Material " << coreMaterial.get_name() << ": its permeability model gives mu_r = " << permeability
                << " (below the vacuum's 1)";
        if (magneticFieldDcBias) {
            message << ", DC bias " << magneticFieldDcBias.value() << " A/m";
        }
        if (magneticFluxDensity) {
            message << ", flux density " << magneticFluxDensity.value() << " T";
        }
        message << ", frequency " << (frequency ? std::to_string(frequency.value()) + " Hz" : std::string("not given"));
        message << ", temperature " << (temperature ? std::to_string(temperature.value()) + " C" : std::string("not given"));
        if (biased) {
            message << ". Without the bias it is still " << unbiasedPermeability.value()
                    << ": the material's frequency or temperature dependence is outside its valid range here";
        }
        else {
            message << ": the material's frequency or temperature dependence is outside its valid range here";
        }
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, message.str());
    }
    return permeability;
}

double InitialPermeability::get_initial_permeability_model_value(const CoreMaterial& coreMaterial,
                                                     std::optional<double> temperature,
                                                     std::optional<double> magneticFieldDcBias,
                                                     std::optional<double> frequency,
                                                     std::optional<double> magneticFluxDensity,
                                                     std::optional<CoreShapeFamily> shapeFamily) {
    const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();

    double initialPermeabilityValue = 1;

    if (std::holds_alternative<PermeabilityPoint>(initialPermeabilityData)) {
        // Only the formula path carries vendor modifiers; the measured-points path below has
        // no per-family variants in MAS data, so shapeFamily is not consulted there.
        initialPermeabilityValue = get_initial_permeability_formula(coreMaterial, temperature, magneticFieldDcBias, frequency, magneticFluxDensity, shapeFamily);
    }
    else {
        double initialPermeabilityValueReference = 1;
        bool hasTemperatureRequirement = temperature.has_value();
        bool hasFrequencyRequirement = frequency.has_value();
        bool hasMagneticFieldDcBiasRequirement = magneticFieldDcBias.has_value();
        bool hasMagneticFluxDensityRequirement = magneticFluxDensity.has_value();
        bool hasTemperatureDependency = has_temperature_dependency(coreMaterial);
        bool hasFrequencyDependency = has_frequency_dependency(coreMaterial);
        bool hasMagneticFieldDcBiasDependency = has_magnetic_field_dc_bias_dependency(coreMaterial);
        double temperatureFactor = 1;
        double frequencyFactor = 1;
        double magneticFieldDcBiasFactor = 1;
        double saturationFactor = 1;

        // ABT #1739: the reference (mu_i at ambient for a temperature-dependent material) only
        // normalises the factors. A request at a temperature that engages no other factor is the
        // tabulated mu_i(T) itself and must not evaluate the ambient point: a material whose table
        // starts above 25 C (DN100H from 26.46 C) would otherwise throw at every temperature.
        bool temperatureOnly = hasTemperatureDependency && hasTemperatureRequirement &&
                               !(hasFrequencyDependency && hasFrequencyRequirement) &&
                               !(hasMagneticFieldDcBiasDependency && hasMagneticFieldDcBiasRequirement) &&
                               !(!hasMagneticFieldDcBiasDependency && (hasMagneticFieldDcBiasRequirement || hasMagneticFluxDensityRequirement));
        if (temperatureOnly) {
            initialPermeabilityValue = get_initial_permeability_temperature_dependent(coreMaterial, temperature.value());
        }
        else if (hasTemperatureDependency) {
            initialPermeabilityValueReference = get_initial_permeability_temperature_dependent(coreMaterial, defaults.ambientTemperature);
        }
        else if (hasFrequencyDependency) {
            const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
            const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
            double minimumFrequency = get_minimum_frequency_in_permeability_points(permeabilityPoints);
            initialPermeabilityValueReference = get_initial_permeability_frequency_dependent(coreMaterial, minimumFrequency);
        }
        else if (hasMagneticFieldDcBiasDependency) {
            const auto& initialPermeabilityData = coreMaterial.get_permeability().get_initial();
            const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
            double minimumMagneticFieldDcBias = get_minimum_magnetic_field_dc_bias_in_permeability_points(permeabilityPoints);
            initialPermeabilityValueReference = get_initial_permeability_magnetic_field_dc_bias_dependent(coreMaterial, minimumMagneticFieldDcBias);
        }
        else {
            // ABT #339: no dependency axis engaged — a single point, or points
            // all at one temperature with no frequency / DC-bias sweep. The
            // reference must come from the measured points (closest to ambient
            // temperature; points without a temperature count as ambient) —
            // leaving the seed value of 1 shipped mu_i = 1 as silent physics.
            const auto& permeabilityPoints = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
            if (permeabilityPoints.empty()) {
                if (coreMaterial.get_name() == DUMMY_SENTINEL_NAME) {
                    // The adviser's synthetic placeholder (Core::resolve_material)
                    // carries no permeability BY CONTRACT — material-aware filters
                    // skip the sentinel, and the exploration pipeline has always
                    // seen mu = 1 for it. Keep that documented behavior; every
                    // real material with an empty list still throws below.
                    return 1;
                }
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Initial permeability list is empty for material: " + coreMaterial.get_name());
            }
            double closestDistance = DBL_MAX;
            for (const auto& point : permeabilityPoints) {
                double pointTemperature = point.get_temperature() ? point.get_temperature().value() : defaults.ambientTemperature;
                double distance = fabs(pointTemperature - defaults.ambientTemperature);
                if (distance < closestDistance) {
                    closestDistance = distance;
                    initialPermeabilityValueReference = point.get_value();
                }
            }
        }

        if (temperatureOnly) {
            // initialPermeabilityValue is the tabulated mu_i(T); no factor applies.
        }
        else {
            if (hasTemperatureDependency && hasTemperatureRequirement) {
                double initialPermeabilityValueTemperatureDependent = get_initial_permeability_temperature_dependent(coreMaterial, temperature.value());
                temperatureFactor = initialPermeabilityValueTemperatureDependent / initialPermeabilityValueReference;
            }

            if (hasFrequencyDependency && hasFrequencyRequirement) {
                double initialPermeabilityValueFrequencyDependent = get_initial_permeability_frequency_dependent(coreMaterial, frequency.value());
                frequencyFactor = initialPermeabilityValueFrequencyDependent / initialPermeabilityValueReference;
            }

            if (hasMagneticFieldDcBiasDependency && hasMagneticFieldDcBiasRequirement) {
                double initialPermeabilityValueMagneticFieldDcBiasDependent = get_initial_permeability_magnetic_field_dc_bias_dependent(coreMaterial, magneticFieldDcBias.value());
                magneticFieldDcBiasFactor = initialPermeabilityValueMagneticFieldDcBiasDependent / initialPermeabilityValueReference;
            }

            // Add variation due to saturation
            if (!hasMagneticFieldDcBiasDependency && (hasMagneticFieldDcBiasRequirement || hasMagneticFluxDensityRequirement)) {
                double amplitudePermeability = initialPermeabilityValueReference;
                double auxTemperatude = defaults.ambientTemperature;
                if (temperature) {
                    auxTemperatude = temperature.value();
                }
                if(hasMagneticFieldDcBiasRequirement) {
                    // Use absolute value of DC bias - both positive and negative bias cause saturation
                    double absDcBias = std::abs(magneticFieldDcBias.value());
                    if (absDcBias <= 0) {
                        amplitudePermeability = initialPermeabilityValueReference;
                    }
                    else {
                        auto amplitudePermeabilityResult = AmplitudePermeability::get_amplitude_permeability(coreMaterial, std::nullopt, absDcBias, auxTemperatude);        
                        if (amplitudePermeabilityResult) {
                            amplitudePermeability = amplitudePermeabilityResult.value();
                        }
                    }
                }
                else {
                    if (magneticFluxDensity.value() <= 0) {
                        amplitudePermeability = initialPermeabilityValueReference;
                    }
                    else {
                        auto amplitudePermeabilityResult = AmplitudePermeability::get_amplitude_permeability(coreMaterial, magneticFluxDensity.value(), std::nullopt, auxTemperatude);        
                        if (amplitudePermeabilityResult) {
                            amplitudePermeability = amplitudePermeabilityResult.value();
                        }
                    }
                }
                if (amplitudePermeability < initialPermeabilityValueReference) {
                    saturationFactor = amplitudePermeability / initialPermeabilityValueReference;
                }
            }
            initialPermeabilityValue = initialPermeabilityValueReference * temperatureFactor * frequencyFactor * magneticFieldDcBiasFactor * saturationFactor;
        }
    }

    // The paramagnetic temperature is never below the Curie temperature, so its table scan is needed only above
    // the Curie point (this runs on every grid point of a magnetisation curve).
    auto curieTemperature = coreMaterial.get_curie_temperature();
    if (temperature && curieTemperature && temperature.value() > curieTemperature.value()) {
        auto paramagneticTemperature = get_paramagnetic_temperature(coreMaterial);
        if (paramagneticTemperature && temperature.value() > paramagneticTemperature.value()) {
            initialPermeabilityValue = 1;
        }
    }

    if (std::isnan(initialPermeabilityValue)) {
        throw NaNResultException("get_initial_permeability, Initial Permeability must be a number, not NaN");
    }

    return initialPermeabilityValue;
}

} // namespace OpenMagnetics
