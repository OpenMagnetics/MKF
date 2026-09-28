// SPDX-License-Identifier: MIT
// Magnetic Blade Runner: physics validation of MAS core-material records. See MaterialValidator.h.
#include "support/MaterialValidator.h"

#include "constructive_models/Core.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/InitialPermeability.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <type_traits>

#include <cmrc/cmrc.hpp>
#include <magic_enum.hpp>
CMRC_DECLARE(data);

namespace OpenMagnetics {

namespace {

// ---------------------------------------------------------------- rule codes
constexpr const char* MAT_PARSE = "MAT_PARSE";
constexpr const char* MAT_PROVENANCE = "MAT_PROVENANCE";
constexpr const char* MAT_LOSS_EVAL = "MAT_LOSS_EVAL";
constexpr const char* MAT_LOSS_OUT_OF_RANGE = "MAT_LOSS_OUT_OF_RANGE";
constexpr const char* MAT_LOSS_MONOTONIC = "MAT_LOSS_MONOTONIC";
constexpr const char* MAT_LOSS_ENVELOPE = "MAT_LOSS_ENVELOPE";
constexpr const char* MAT_LOSS_TEMPERATURE = "MAT_LOSS_TEMPERATURE";
constexpr const char* MAT_LOSS_HYSTERESIS_BOUND = "MAT_LOSS_HYSTERESIS_BOUND";
constexpr const char* MAT_LOSS_POINTS_INCONSISTENT = "MAT_LOSS_POINTS_INCONSISTENT";
constexpr const char* MAT_LOSS_FACTOR = "MAT_LOSS_FACTOR";
constexpr const char* MAT_BSAT_CEILING = "MAT_BSAT_CEILING";
constexpr const char* MAT_BSAT_CLASS = "MAT_BSAT_CLASS";
constexpr const char* MAT_CURIE = "MAT_CURIE";
constexpr const char* MAT_CURIE_SATURATION = "MAT_CURIE_SATURATION";
constexpr const char* MAT_PERM = "MAT_PERM";

// ---------------------------------------------------------------- physical constants and rule thresholds
// Saturation polarisation of pure iron, the highest of any material MAS can describe (its
// composition enum has no Fe-Co alloy): Ms = 1707 emu/cm^3 at 20 C and 1746 emu/cm^3 at 0 K,
// i.e. Js = 2.145 T and 2.194 T (Cullity & Graham, Introduction to Magnetic Materials, 2nd ed.,
// Table 4.2). B measured at a field H is Js + mu0*H at most.
constexpr double kIronPolarisationRoomTemperature = 2.16;  // T, used for T >= 0 C (2.145 T rounded up)
constexpr double kIronPolarisationZeroKelvin = 2.20;       // T, used below 0 C (2.194 T rounded up)
constexpr double kVacuumPermeability = 1.25663706212e-6;

// Envelope severities (Alf, 2026-09-28): below best/2 or above worst*2 SUSPICIOUS, below best/10
// or above worst*10 IMPOSSIBLE; WEAK classes stop at SUSPICIOUS.
constexpr double kEnvelopeSuspiciousFactor = 2.0;
constexpr double kEnvelopeImpossibleFactor = 10.0;
// Static class bounds: datasheet mu_i tolerance is +-25..30 %, so outside [min/2, max*2] is
// implausible; Bs and Tc outside the published extreme by more than 10 %.
constexpr double kPermeabilityMargin = 2.0;
constexpr double kSaturationMargin = 1.10;
constexpr double kCurieMargin = 0.10;
// Monotonicity: a loss that falls by more than 1 % when f or B rises is not a fit wobble.
constexpr double kMonotonicTolerance = 0.01;
// Hysteresis bound: at f <= 10 kHz eddy losses are negligible, so the energy per cycle Pv/f cannot
// exceed the rectangle 2Hc x 2Bpk that bounds any loop inside the major loop.
constexpr double kHysteresisMaximumFrequency = 10e3;
constexpr double kHysteresisSuspiciousFactor = 2.0;
constexpr double kHysteresisImpossibleFactor = 5.0;
// Measured points: same frequency (1 %), temperature within 2.5 C.
constexpr double kPointsSuspiciousFactor = 2.0;
constexpr double kPointsImpossibleFactor = 10.0;
constexpr double kPointsTemperatureWindow = 2.5;
// Temperature dependence of a ferrite: a jump of more than 2x per 5 C is not physical.
constexpr double kTemperatureStepJumpFactor = 2.0;

// ---------------------------------------------------------------- provenance detection
// MAS RFC 0011 (Draft) adds `provenance` to records. The generated CoreMaterial gets
// get_provenance() once the schema carries it; this switches the rule to per-record WARNING.
template <typename T, typename = void>
struct HasProvenance : std::false_type {};
template <typename T>
struct HasProvenance<T, std::void_t<decltype(std::declval<T&>().get_provenance())>> : std::true_type {};
constexpr bool kSchemaHasProvenance = HasProvenance<CoreMaterial>::value;
const char* kProvenanceSchemaGap =
    "MAT_PROVENANCE: the MAS core-material schema has no provenance field (RFC 0011 Draft)";

// ---------------------------------------------------------------- context
struct Ctx {
    std::string reference;
    std::vector<MaterialFinding>* out;
};

void emit(const Ctx& ctx, const std::string& code, MaterialFindingSeverity severity, double value, double threshold,
          const std::string& message) {
    // Revert harness (as Blade Runner's TAS_VALIDATOR_SUPPRESS): the listed codes behave as if the
    // check had never been written, so a test's must-fire case can be shown to go red per rule.
    if (const char* suppress = std::getenv("MKF_MATERIAL_VALIDATOR_SUPPRESS"); suppress != nullptr && *suppress != '\0') {
        const std::string list = std::string(",") + suppress + ",";
        if (list.find("," + code + ",") != std::string::npos) return;
    }
    MaterialFinding f;
    f.code = code;
    f.severity = severity;
    f.reference = ctx.reference;
    f.message = message;
    f.value = value;
    f.threshold = threshold;
    ctx.out->push_back(std::move(f));
}

std::string fmt(double v, int precision = 4) {
    std::ostringstream s;
    s.precision(precision);
    s << v;
    return s.str();
}

std::string point_label(double f, double b, double t) {
    return fmt(f / 1e3) + " kHz / " + fmt(b * 1e3) + " mT / " + fmt(t) + " C";
}

void erase_name_keyed_caches(const std::string& name) {
    initialPermeabilityTemperatureInterps.erase(name);
    initialPermeabilityFrequencyInterps.erase(name);
    initialPermeabilityMagneticFieldDcBiasInterps.erase(name);
    lossFactorInterps.erase(name);
}

struct CacheGuard {
    std::string name;
    explicit CacheGuard(std::string n) : name(std::move(n)) { erase_name_keyed_caches(name); }
    ~CacheGuard() { erase_name_keyed_caches(name); }
};

bool is_steinmetz_family(CoreLossesModels model) {
    switch (model) {
        case CoreLossesModels::STEINMETZ:
        case CoreLossesModels::IGSE:
        case CoreLossesModels::CIGSE:
        case CoreLossesModels::BARG:
        case CoreLossesModels::ALBACH:
        case CoreLossesModels::NSE:
        case CoreLossesModels::MSE:
            return true;
        default:
            return false;
    }
}

// ---------------------------------------------------------------- loss evaluator
class LossEvaluator {
  public:
    LossEvaluator(const CoreMaterial& material) : _material(material) {}

    // Picks the model MKF itself would use: the first of the settings' model order that the
    // material's data supports (CoreLosses::get_core_losses_model, without the by-name lookup).
    std::string select() {
        auto available = CoreLossesModel::get_methods(_material);
        for (auto model : settings.get_core_losses_model_names()) {
            if (std::find(available.begin(), available.end(), model) != available.end()) {
                _model = model;
                _modelInstance = CoreLossesModel::factory(model);
                _modelName = std::string(magic_enum::enum_name(model));
                return "";
            }
        }
        if (available.empty()) {
            return "the record has no loss model MKF can evaluate (no steinmetz/roshen/proprietary/lossFactor/magnetec data)";
        }
        return "none of the material's loss models is in the settings' model order";
    }

    bool selected() const { return _modelInstance != nullptr; }
    CoreLossesModels model() const { return _model; }
    const std::string& model_name() const { return _modelName; }

    // Frequency spans of the fitted Steinmetz ranges, when the model is a Steinmetz-family one.
    bool has_span() const { return is_steinmetz_family(_model); }
    std::vector<std::pair<double, double>> spans() const {
        std::vector<std::pair<double, double>> out;
        auto data = CoreLossesModel::get_method_data(_material, "steinmetz");
        if (!data.get_ranges()) {
            throw std::runtime_error("steinmetz method without ranges");
        }
        const auto ranges = data.get_ranges().value();
        for (const auto& range : ranges) {
            if (!range.get_minimum_frequency() || !range.get_maximum_frequency()) {
                throw std::runtime_error("a steinmetz range has no minimumFrequency/maximumFrequency");
            }
            out.emplace_back(range.get_minimum_frequency().value(), range.get_maximum_frequency().value());
        }
        return out;
    }

    double evaluate(double frequency, double peak, double temperature) {
        json excitationJson;
        excitationJson["frequency"] = frequency;
        excitationJson["magneticFluxDensity"]["processed"]["dutyCycle"] = 0.5;
        excitationJson["magneticFluxDensity"]["processed"]["label"] = "sinusoidal";
        excitationJson["magneticFluxDensity"]["processed"]["offset"] = 0;
        excitationJson["magneticFluxDensity"]["processed"]["peak"] = peak;
        excitationJson["magneticFluxDensity"]["processed"]["peakToPeak"] = 2 * peak;
        excitationJson["magneticFieldStrength"]["processed"]["offset"] = 0;
        excitationJson["magneticFieldStrength"]["processed"]["label"] = "sinusoidal";
        excitationJson["magneticFieldStrength"]["processed"]["peakToPeak"] = 0;
        OperatingPointExcitation excitation(excitationJson);
        if (_model == CoreLossesModels::PROPRIETARY && !has_volumetric_proprietary()) {
            // Magnetec publishes mass losses: W/kg x density.
            double density = Core::get_density(_material);
            if (!std::isfinite(density) || density <= 0) {
                throw std::runtime_error("mass-loss model needs the material density, which the record does not give");
            }
            return _modelInstance->get_core_mass_losses(_material, excitation, temperature) * density;
        }
        return _modelInstance->get_core_volumetric_losses(_material, excitation, temperature);
    }

  private:
    bool has_volumetric_proprietary() const {
        for (const auto& [key, methods] : _material.get_volumetric_losses()) {
            for (const auto& m : methods) {
                if (std::holds_alternative<CoreLossesMethodData>(m)) {
                    auto method = std::get<CoreLossesMethodData>(m).get_method();
                    if (method == VolumetricCoreLossesMethodType::MAGNETICS || method == VolumetricCoreLossesMethodType::MICROMETALS ||
                        method == VolumetricCoreLossesMethodType::POCO || method == VolumetricCoreLossesMethodType::TDG) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    const CoreMaterial& _material;
    CoreLossesModels _model = CoreLossesModels::STEINMETZ;
    std::shared_ptr<CoreLossesModel> _modelInstance;
    std::string _modelName;
};

bool in_spans(const std::vector<std::pair<double, double>>& spans, double f) {
    for (const auto& [lo, hi] : spans) {
        if (f >= lo && f <= hi) return true;
    }
    return false;
}

std::string spans_label(const std::vector<std::pair<double, double>>& spans) {
    std::string s;
    for (const auto& [lo, hi] : spans) {
        if (!s.empty()) s += ", ";
        s += fmt(lo / 1e3) + "-" + fmt(hi / 1e3) + " kHz";
    }
    return s;
}

// The in-span frequency nearest (in log) to `target`.
double nearest_in_span(const std::vector<std::pair<double, double>>& spans, double target) {
    double best = std::numeric_limits<double>::quiet_NaN();
    double bestDistance = std::numeric_limits<double>::max();
    for (const auto& [lo, hi] : spans) {
        double candidate = std::clamp(target, lo, hi);
        double distance = std::fabs(std::log(candidate / target));
        if (distance < bestDistance) {
            bestDistance = distance;
            best = candidate;
        }
    }
    return best;
}

std::vector<double> log_grid(double lo, double hi, int n) {
    std::vector<double> out;
    if (!(hi > lo)) {
        out.push_back(lo);
        return out;
    }
    for (int i = 0; i < n; ++i) {
        out.push_back(lo * std::pow(hi / lo, static_cast<double>(i) / (n - 1)));
    }
    return out;
}

bool is_ferrite_or_powder(const json& record) {
    std::string type = record.contains("material") && record["material"].is_string() ? record["material"].get<std::string>() : "";
    return type == "ferrite" || type == "powder";
}

void add_points(const json& pointsArray, std::vector<MaterialLossPoint>& out) {
    for (const auto& p : pointsArray) {
        if (!p.is_object() || !p.contains("value") || !p.contains("temperature") || !p.contains("magneticFluxDensity")) continue;
        const auto& excitation = p["magneticFluxDensity"];
        if (!excitation.contains("frequency") || !excitation.contains("magneticFluxDensity")) continue;
        const auto& signal = excitation["magneticFluxDensity"];
        if (!signal.contains("processed")) continue;
        const auto& processed = signal["processed"];
        if (!processed.contains("label") || processed["label"] != "sinusoidal") continue;
        if (!processed.contains("peak")) continue;
        if (processed.contains("offset") && std::fabs(processed["offset"].get<double>()) > 1e-12) continue;
        MaterialLossPoint point{excitation["frequency"].get<double>(), processed["peak"].get<double>(), p["temperature"].get<double>(),
                                p["value"].get<double>(), p.contains("origin") && p["origin"].is_string() ? p["origin"].get<std::string>() : ""};
        out.push_back(point);
    }
}

// ---------------------------------------------------------------- static-property checks
void check_saturation(const json& record, const CoreMaterial& material, const MaterialClassBounds* cls, const Ctx& ctx,
                      std::vector<std::string>& skipped) {
    if (!record.contains("saturation") || !record["saturation"].is_array() || record["saturation"].empty()) {
        skipped.push_back(std::string(MAT_BSAT_CEILING) + ": the record has no saturation data");
        skipped.push_back(std::string(MAT_BSAT_CLASS) + ": the record has no saturation data");
        return;
    }
    for (const auto& datum : material.get_saturation()) {
        double b = datum.get_magnetic_flux_density();
        double h = datum.get_magnetic_field();
        double t = datum.get_temperature();
        double ceiling = (t >= 0 ? kIronPolarisationRoomTemperature : kIronPolarisationZeroKelvin) + kVacuumPermeability * std::fabs(h);
        if (!std::isfinite(b) || b <= 0) {
            emit(ctx, MAT_BSAT_CEILING, MaterialFindingSeverity::IMPOSSIBLE, b, 0,
                 "saturation flux density " + fmt(b) + " T at " + fmt(t) + " C is not a positive number");
        }
        else if (b > ceiling) {
            emit(ctx, MAT_BSAT_CEILING, MaterialFindingSeverity::IMPOSSIBLE, b, ceiling,
                 "saturation flux density " + fmt(b) + " T at " + fmt(t) + " C, " + fmt(h) + " A/m exceeds pure iron's Js + mu0*H = " +
                     fmt(ceiling) + " T, the highest of any material MAS can describe");
        }
    }
    if (cls == nullptr) {
        skipped.push_back(std::string(MAT_BSAT_CLASS) + ": unclassified");
        return;
    }
    if (!cls->maximumSaturationFluxDensity) {
        skipped.push_back(std::string(MAT_BSAT_CLASS) + ": " + cls->name + ": no published class maximum");
        return;
    }
    double bs25 = Core::get_magnetic_flux_density_saturation(material, 25.0, false);
    double limit = *cls->maximumSaturationFluxDensity * kSaturationMargin;
    if (bs25 > limit) {
        emit(ctx, MAT_BSAT_CLASS, MaterialFindingSeverity::SUSPICIOUS, bs25, limit,
             "saturation flux density at 25 C " + fmt(bs25) + " T is above the " + cls->name + " maximum " +
                 fmt(*cls->maximumSaturationFluxDensity) + " T (+10 %); " + cls->saturationSource);
    }
}

void check_curie(const json& record, const CoreMaterial& material, const MaterialClassBounds* cls, const Ctx& ctx,
                 std::vector<std::string>& skipped) {
    double tc = Core::get_curie_temperature(material);
    if (!std::isfinite(tc)) {
        skipped.push_back(std::string(MAT_CURIE) + ": the record has no Curie temperature");
        skipped.push_back(std::string(MAT_CURIE_SATURATION) + ": the record has no Curie temperature");
        return;
    }
    // A ferromagnet has no spontaneous magnetisation above Tc: a saturation point well above it
    // contradicts the record's own Curie temperature.
    if (record.contains("saturation") && record["saturation"].is_array()) {
        for (const auto& datum : material.get_saturation()) {
            if (datum.get_temperature() > tc + 10 && datum.get_magnetic_flux_density() > 0.05) {
                emit(ctx, MAT_CURIE_SATURATION, MaterialFindingSeverity::IMPOSSIBLE, datum.get_temperature(), tc,
                     "saturation " + fmt(datum.get_magnetic_flux_density()) + " T reported at " + fmt(datum.get_temperature()) +
                         " C, above the record's own Curie temperature " + fmt(tc) + " C");
            }
        }
    }
    if (cls == nullptr) {
        skipped.push_back(std::string(MAT_CURIE) + ": unclassified");
        return;
    }
    if (!cls->minimumCurieTemperature && !cls->maximumCurieTemperature) {
        skipped.push_back(std::string(MAT_CURIE) + ": " + cls->name + ": no published class range");
        return;
    }
    if (cls->minimumCurieTemperature) {
        double limit = *cls->minimumCurieTemperature * (1 - kCurieMargin);
        if (tc < limit) {
            emit(ctx, MAT_CURIE, MaterialFindingSeverity::SUSPICIOUS, tc, limit,
                 "Curie temperature " + fmt(tc) + " C is below the " + cls->name + " minimum " + fmt(*cls->minimumCurieTemperature) +
                     " C (-10 %); " + cls->curieSource);
        }
    }
    if (cls->maximumCurieTemperature) {
        double limit = *cls->maximumCurieTemperature * (1 + kCurieMargin);
        if (tc > limit) {
            emit(ctx, MAT_CURIE, MaterialFindingSeverity::SUSPICIOUS, tc, limit,
                 "Curie temperature " + fmt(tc) + " C is above the " + cls->name + " maximum " + fmt(*cls->maximumCurieTemperature) +
                     " C (+10 %); " + cls->curieSource);
        }
    }
}

void check_permeability(const CoreMaterial& material, const MaterialClassBounds* cls, const Ctx& ctx, std::vector<std::string>& skipped) {
    double mu;
    try {
        mu = InitialPermeability::get_initial_permeability(material, 25.0);
    }
    catch (const std::exception& e) {
        emit(ctx, MAT_PERM, MaterialFindingSeverity::SUSPICIOUS, 0, 0,
             std::string("MKF cannot evaluate the initial permeability at 25 C: ") + e.what());
        return;
    }
    if (!std::isfinite(mu) || mu < 1) {
        emit(ctx, MAT_PERM, MaterialFindingSeverity::IMPOSSIBLE, mu, 1,
             "initial permeability at 25 C is " + fmt(mu) + ", below 1: no magnetic core material is diamagnetic");
        return;
    }
    if (cls == nullptr) {
        skipped.push_back(std::string(MAT_PERM) + ": unclassified");
        return;
    }
    if (!cls->minimumInitialPermeability && !cls->maximumInitialPermeability) {
        skipped.push_back(std::string(MAT_PERM) + ": " + cls->name + ": no published class range");
        return;
    }
    if (cls->minimumInitialPermeability && mu < *cls->minimumInitialPermeability / kPermeabilityMargin) {
        emit(ctx, MAT_PERM, MaterialFindingSeverity::SUSPICIOUS, mu, *cls->minimumInitialPermeability / kPermeabilityMargin,
             "initial permeability at 25 C " + fmt(mu) + " is below half the " + cls->name + " minimum " +
                 fmt(*cls->minimumInitialPermeability) + "; " + cls->initialPermeabilitySource);
    }
    if (cls->maximumInitialPermeability && mu > *cls->maximumInitialPermeability * kPermeabilityMargin) {
        emit(ctx, MAT_PERM, MaterialFindingSeverity::SUSPICIOUS, mu, *cls->maximumInitialPermeability * kPermeabilityMargin,
             "initial permeability at 25 C " + fmt(mu) + " is above twice the " + cls->name + " maximum " +
                 fmt(*cls->maximumInitialPermeability) + "; " + cls->initialPermeabilitySource);
    }
}

void check_loss_factor(const json& record, const Ctx& ctx) {
    if (!record.contains("volumetricLosses") || !record["volumetricLosses"].is_object()) return;
    for (const auto& [key, methods] : record["volumetricLosses"].items()) {
        if (!methods.is_array()) continue;
        for (const auto& m : methods) {
            if (!m.is_object() || !m.contains("method") || m["method"] != "lossFactor" || !m.contains("factors")) continue;
            for (const auto& factor : m["factors"]) {
                double value = factor.value("value", std::numeric_limits<double>::quiet_NaN());
                double frequency = factor.value("frequency", std::numeric_limits<double>::quiet_NaN());
                if (!std::isfinite(value) || value <= 0) {
                    emit(ctx, MAT_LOSS_FACTOR, MaterialFindingSeverity::IMPOSSIBLE, value, 0,
                         "loss factor tan(delta)/mu " + fmt(value) + " at " + fmt(frequency / 1e3) + " kHz is not positive: a passive core cannot generate energy");
                }
            }
        }
    }
}

// ---------------------------------------------------------------- measured points
void check_hysteresis_bound(const json& record, const CoreMaterial& material, const std::vector<MaterialLossPoint>& points, const Ctx& ctx,
                            std::vector<std::string>& skipped) {
    if (!is_ferrite_or_powder(record)) {
        skipped.push_back(std::string(MAT_LOSS_HYSTERESIS_BOUND) + ": applies to ferrites and powders only");
        return;
    }
    bool anyLowFrequency = false;
    for (const auto& p : points) {
        if (p.frequency > kHysteresisMaximumFrequency) continue;
        anyLowFrequency = true;
        double hc = Core::get_coercive_force(material, p.temperature, false);
        if (!std::isfinite(hc) || hc <= 0) {
            skipped.push_back(std::string(MAT_LOSS_HYSTERESIS_BOUND) + ": the record has no coercive force");
            return;
        }
        double energy = p.volumetricLosses / p.frequency;
        double bound = 4 * hc * p.magneticFluxDensityPeak;
        double ratio = energy / bound;
        if (ratio > kHysteresisSuspiciousFactor) {
            emit(ctx, MAT_LOSS_HYSTERESIS_BOUND,
                 ratio > kHysteresisImpossibleFactor ? MaterialFindingSeverity::IMPOSSIBLE : MaterialFindingSeverity::SUSPICIOUS, energy, bound,
                 "measured point " + point_label(p.frequency, p.magneticFluxDensityPeak, p.temperature) + " (" + p.origin + "): " +
                     fmt(p.volumetricLosses / 1e3) + " kW/m3 is " + fmt(energy) + " J/m3 per cycle, " + fmt(ratio, 3) +
                     "x the hysteresis bound 4*Hc*Bpk = " + fmt(bound) + " J/m3 (Hc " + fmt(hc) + " A/m)");
        }
    }
    if (!anyLowFrequency) {
        skipped.push_back(std::string(MAT_LOSS_HYSTERESIS_BOUND) + ": no measured point at or below 10 kHz");
    }
}

void check_points_consistency(const std::vector<MaterialLossPoint>& points, const Ctx& ctx, std::vector<std::string>& skipped) {
    if (points.size() < 3) {
        skipped.push_back(std::string(MAT_LOSS_POINTS_INCONSISTENT) + ": fewer than 3 measured points");
        return;
    }
    bool anyCompared = false;
    std::set<std::pair<size_t, size_t>> reportedPairs;
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        if (p.volumetricLosses <= 0 || p.magneticFluxDensityPeak <= 0) continue;
        // Neighbours at the same frequency and temperature, one below and one above p's B.
        std::optional<size_t> below, above;
        for (size_t j = 0; j < points.size(); ++j) {
            if (j == i) continue;
            const auto& q = points[j];
            if (q.volumetricLosses <= 0 || q.magneticFluxDensityPeak <= 0) continue;
            if (std::fabs(q.frequency - p.frequency) > 0.01 * p.frequency) continue;
            if (std::fabs(q.temperature - p.temperature) > kPointsTemperatureWindow) continue;
            if (q.magneticFluxDensityPeak < p.magneticFluxDensityPeak * 0.99) {
                if (!below || q.magneticFluxDensityPeak > points[*below].magneticFluxDensityPeak) below = j;
            }
            else if (q.magneticFluxDensityPeak > p.magneticFluxDensityPeak * 1.01) {
                if (!above || q.magneticFluxDensityPeak < points[*above].magneticFluxDensityPeak) above = j;
            }
        }
        if (!below || !above) continue;
        anyCompared = true;
        const auto& a = points[*below];
        const auto& b = points[*above];
        double x = std::log(p.magneticFluxDensityPeak / a.magneticFluxDensityPeak) / std::log(b.magneticFluxDensityPeak / a.magneticFluxDensityPeak);
        double interpolated = std::exp(std::log(a.volumetricLosses) + x * (std::log(b.volumetricLosses) - std::log(a.volumetricLosses)));
        double ratio = std::max(p.volumetricLosses / interpolated, interpolated / p.volumetricLosses);
        if (ratio >= kPointsSuspiciousFactor) {
            emit(ctx, MAT_LOSS_POINTS_INCONSISTENT,
                 ratio >= kPointsImpossibleFactor ? MaterialFindingSeverity::IMPOSSIBLE : MaterialFindingSeverity::SUSPICIOUS, p.volumetricLosses,
                 interpolated,
                 "measured point " + point_label(p.frequency, p.magneticFluxDensityPeak, p.temperature) + " (" + p.origin + ") is " +
                     fmt(p.volumetricLosses / 1e3) + " kW/m3, but the same-frequency, same-temperature points at " +
                     fmt(a.magneticFluxDensityPeak * 1e3) + " and " + fmt(b.magneticFluxDensityPeak * 1e3) + " mT interpolate to " +
                     fmt(interpolated / 1e3) + " kW/m3 (" + fmt(ratio, 3) + "x apart)");
        }
    }
    if (!anyCompared) {
        skipped.push_back(std::string(MAT_LOSS_POINTS_INCONSISTENT) + ": no measured point is bracketed in B by two others at the same f and T");
    }
}

// ---------------------------------------------------------------- model checks
struct EvalResult {
    bool ok;
    double value;
    std::string error;
};

EvalResult safe_evaluate(LossEvaluator& evaluator, double f, double b, double t) {
    try {
        double v = evaluator.evaluate(f, b, t);
        if (!std::isfinite(v) || v <= 0) {
            return {false, v, "the model returns " + fmt(v)};
        }
        return {true, v, ""};
    }
    catch (const std::exception& e) {
        return {false, std::numeric_limits<double>::quiet_NaN(), std::string("the model throws: ") + e.what()};
    }
}

void check_loss_model(const json& record, const CoreMaterial& material, const MaterialClassBounds* cls, const Ctx& ctx,
                      std::vector<std::string>& skipped, std::string& modelName) {
    const std::vector<std::string> modelRules = {MAT_LOSS_EVAL, MAT_LOSS_OUT_OF_RANGE, MAT_LOSS_MONOTONIC, MAT_LOSS_ENVELOPE, MAT_LOSS_TEMPERATURE};
    auto skipAll = [&](const std::string& reason) {
        for (const auto& rule : modelRules) skipped.push_back(rule + ": " + reason);
    };

    LossEvaluator evaluator(material);
    std::string why = evaluator.select();
    if (!evaluator.selected()) {
        skipAll(why);
        return;
    }
    modelName = evaluator.model_name();
    if (evaluator.model() == CoreLossesModels::LOSS_FACTOR) {
        skipAll("the loss-factor model needs a magnetizing inductance (a core), so it has no per-volume evaluation");
        return;
    }

    std::vector<std::pair<double, double>> spans;
    if (evaluator.has_span()) {
        try {
            spans = evaluator.spans();
        }
        catch (const std::exception& e) {
            emit(ctx, MAT_LOSS_EVAL, MaterialFindingSeverity::IMPOSSIBLE, 0, 0,
                 std::string("the Steinmetz fit's frequency span cannot be read: ") + e.what());
            return;
        }
        if (spans.empty()) {
            emit(ctx, MAT_LOSS_EVAL, MaterialFindingSeverity::IMPOSSIBLE, 0, 0, "the Steinmetz method has no ranges");
            return;
        }
    }
    auto inSpan = [&](double f) { return !evaluator.has_span() || in_spans(spans, f); };

    // Saturation bounds the flux densities the grids may use.
    double bsat = std::numeric_limits<double>::quiet_NaN();
    try {
        bsat = Core::get_magnetic_flux_density_saturation(material, 100.0, false);
    }
    catch (const std::exception&) {
    }

    // ---- reference points of the class: OUT_OF_RANGE, EVAL, ENVELOPE
    const MaterialLossReferencePoint* primary = nullptr;
    if (cls == nullptr) {
        skipped.push_back(std::string(MAT_LOSS_ENVELOPE) + ": unclassified");
        skipped.push_back(std::string(MAT_LOSS_OUT_OF_RANGE) + ": unclassified (no class reference point)");
    }
    else if (cls->lossReferencePoints.empty()) {
        skipped.push_back(std::string(MAT_LOSS_ENVELOPE) + ": " + cls->name + ": " + cls->lossAnchorGap);
        skipped.push_back(std::string(MAT_LOSS_OUT_OF_RANGE) + ": " + cls->name + ": no class reference point");
    }
    else {
        for (const auto& ref : cls->lossReferencePoints) {
            const std::string label = point_label(ref.frequency, ref.magneticFluxDensityPeak, ref.temperature);
            if (ref.primary) primary = &ref;
            if (!ref.primary && !evaluator.has_span()) {
                // A secondary point (e.g. the MHz points of MnZn power ferrite) belongs to the grades
                // fitted there; a model with no declared span says nothing about where it is valid.
                skipped.push_back(std::string(MAT_LOSS_ENVELOPE) + ": " + label + ": secondary point, and the " + modelName +
                                  " model declares no fitted span");
                continue;
            }
            if (!inSpan(ref.frequency)) {
                if (ref.primary) {
                    emit(ctx, MAT_LOSS_OUT_OF_RANGE, MaterialFindingSeverity::SUSPICIOUS, ref.frequency, 0,
                         "the " + cls->name + " reference point " + label + " is outside the fitted Steinmetz span (" + spans_label(spans) +
                             "): a consumer using this grade at its class's standard point would be extrapolating");
                }
                else {
                    skipped.push_back(std::string(MAT_LOSS_ENVELOPE) + ": " + label + " outside the fitted span");
                }
                continue;
            }
            auto r = safe_evaluate(evaluator, ref.frequency, ref.magneticFluxDensityPeak, ref.temperature);
            if (!r.ok) {
                emit(ctx, MAT_LOSS_EVAL, MaterialFindingSeverity::IMPOSSIBLE, r.value, 0,
                     "the " + modelName + " loss at " + label + " (inside the fitted span) is not a positive number: " + r.error);
                continue;
            }
            double lowSus = ref.best.volumetricLosses / kEnvelopeSuspiciousFactor;
            double lowImp = ref.best.volumetricLosses / kEnvelopeImpossibleFactor;
            double highSus = ref.worst.volumetricLosses * kEnvelopeSuspiciousFactor;
            double highImp = ref.worst.volumetricLosses * kEnvelopeImpossibleFactor;
            auto severity = [&](bool impossible) {
                return (impossible && !cls->weak) ? MaterialFindingSeverity::IMPOSSIBLE : MaterialFindingSeverity::SUSPICIOUS;
            };
            const std::string quality = cls->weak ? " [WEAK class: " + cls->anchorQuality + "]" : "";
            if (r.value < lowSus) {
                bool impossible = r.value < lowImp;
                emit(ctx, MAT_LOSS_ENVELOPE, severity(impossible), r.value, impossible ? lowImp : lowSus,
                     modelName + " loss at " + label + " is " + fmt(r.value / 1e3) + " kW/m3, below " + (impossible ? "a tenth" : "half") +
                         " of the best published " + cls->name + " grade (" + ref.best.manufacturer + " " + ref.best.grade + " " +
                         fmt(ref.best.volumetricLosses / 1e3) + " kW/m3; " + ref.best.source + ")" + quality);
            }
            else if (ref.primary && r.value > highSus) {
                // High side only at the primary point: its worst anchor covers every grade of the
                // class (legacy and high-Bs included), while a secondary point's worst anchor only
                // covers the grades marketed there (a 100 kHz grade at 1 MHz is legitimately far
                // above the worst MHz grade, and TDK/ACME publish nothing at 100 kHz / 100 mT).
                bool impossible = r.value > highImp;
                emit(ctx, MAT_LOSS_ENVELOPE, severity(impossible), r.value, impossible ? highImp : highSus,
                     modelName + " loss at " + label + " is " + fmt(r.value / 1e3) + " kW/m3, above " + (impossible ? "ten times" : "twice") +
                         " the worst published " + cls->name + " grade (" + ref.worst.manufacturer + " " + ref.worst.grade + " " +
                         fmt(ref.worst.volumetricLosses / 1e3) + " kW/m3; " + ref.worst.source + ")" + quality);
            }
        }
    }

    // ---- probe point for the shape checks
    double probeFrequency = primary ? primary->frequency : 100e3;
    double probePeak = primary ? primary->magneticFluxDensityPeak : 0.05;
    double probeTemperature = primary ? primary->temperature : (record.value("material", "") == "ferrite" ? 100.0 : 25.0);
    if (evaluator.has_span() && !in_spans(spans, probeFrequency)) {
        probeFrequency = nearest_in_span(spans, probeFrequency);
    }
    if (std::isfinite(bsat) && probePeak > 0.5 * bsat) {
        probePeak = 0.5 * bsat;
    }

    // ---- MAT_LOSS_EVAL + MAT_LOSS_MONOTONIC: grids in f and B inside the span
    double fLo = 10e3, fHi = 500e3;
    if (evaluator.has_span()) {
        fLo = std::numeric_limits<double>::max();
        fHi = 0;
        for (const auto& [lo, hi] : spans) {
            fLo = std::min(fLo, lo);
            fHi = std::max(fHi, hi);
        }
        fLo = std::max(fLo, 1e3);
        fHi = std::min(fHi, 10e6);
    }
    // The fits carry no flux-density span. Across the whole frequency span only a low flux
    // density is inside the data every range was fitted to (MHz ranges are measured at <= 50 mT),
    // so the frequency sweep runs at 50 mT (or the class point, if lower); higher B at MHz is an
    // extrapolation in B, where independently-fitted ranges legitimately disagree at their seams.
    const double frequencyGridPeak = std::min(0.05, probePeak);
    auto checkGrid = [&](const std::vector<double>& axis, bool isFrequency) {
        std::optional<std::pair<double, double>> previous;  // (axis value, loss)
        for (double x : axis) {
            double f = isFrequency ? x : probeFrequency;
            double b = isFrequency ? frequencyGridPeak : x;
            if (!inSpan(f)) continue;
            auto r = safe_evaluate(evaluator, f, b, probeTemperature);
            if (!r.ok) {
                emit(ctx, MAT_LOSS_EVAL, MaterialFindingSeverity::IMPOSSIBLE, r.value, 0,
                     "the " + modelName + " loss at " + point_label(f, b, probeTemperature) + " (inside the fitted span) is not a positive number: " + r.error);
                return;
            }
            if (previous && r.value < previous->second * (1 - kMonotonicTolerance)) {
                emit(ctx, MAT_LOSS_MONOTONIC, MaterialFindingSeverity::IMPOSSIBLE, r.value, previous->second,
                     modelName + " loss falls from " + fmt(previous->second / 1e3) + " to " + fmt(r.value / 1e3) + " kW/m3 as " +
                         (isFrequency ? "frequency rises from " + fmt(previous->first / 1e3) + " to " + fmt(x / 1e3) + " kHz (at " + fmt(b * 1e3) + " mT, "
                                      : "flux density rises from " + fmt(previous->first * 1e3) + " to " + fmt(x * 1e3) + " mT (at " + fmt(f / 1e3) + " kHz, ") +
                         fmt(probeTemperature) + " C): losses must grow with f and B");
                return;
            }
            previous = std::make_pair(x, r.value);
        }
    };
    checkGrid(log_grid(fLo, fHi, 12), true);
    double bHi = std::isfinite(bsat) ? std::min(0.3, 0.8 * bsat) : 0.2;
    checkGrid(log_grid(0.005, std::max(bHi, 0.006), 8), false);

    // ---- MAT_LOSS_TEMPERATURE: ferrites depend strongly on T (a valley near 80-100 C)
    if (record.value("material", "") != "ferrite") {
        skipped.push_back(std::string(MAT_LOSS_TEMPERATURE) + ": applies to ferrites only");
        return;
    }
    std::vector<std::pair<double, double>> curve;
    for (double t = 25; t <= 120 + 1e-9; t += 5) {
        auto r = safe_evaluate(evaluator, probeFrequency, probePeak, t);
        if (!r.ok) {
            emit(ctx, MAT_LOSS_EVAL, MaterialFindingSeverity::IMPOSSIBLE, r.value, 0,
                 "the " + modelName + " loss at " + point_label(probeFrequency, probePeak, t) + " (inside the fitted span) is not a positive number: " + r.error);
            return;
        }
        curve.emplace_back(t, r.value);
    }
    double minimum = curve.front().second, maximum = curve.front().second;
    for (size_t i = 0; i < curve.size(); ++i) {
        minimum = std::min(minimum, curve[i].second);
        maximum = std::max(maximum, curve[i].second);
        if (i > 0) {
            double jump = std::max(curve[i].second / curve[i - 1].second, curve[i - 1].second / curve[i].second);
            if (jump > kTemperatureStepJumpFactor) {
                emit(ctx, MAT_LOSS_TEMPERATURE, MaterialFindingSeverity::SUSPICIOUS, jump, kTemperatureStepJumpFactor,
                     modelName + " loss at " + fmt(probeFrequency / 1e3) + " kHz / " + fmt(probePeak * 1e3) + " mT changes " + fmt(jump, 3) +
                         "x between " + fmt(curve[i - 1].first) + " and " + fmt(curve[i].first) + " C");
                return;
            }
        }
    }
    if (maximum <= minimum * (1 + 1e-9)) {
        emit(ctx, MAT_LOSS_TEMPERATURE, MaterialFindingSeverity::SUSPICIOUS, maximum, minimum,
             modelName + " loss at " + fmt(probeFrequency / 1e3) + " kHz / " + fmt(probePeak * 1e3) +
                 " mT is identical from 25 to 120 C: a ferrite's loss depends strongly on temperature, so the fit's temperature coefficients are missing or ignored");
    }
}

} // namespace

// ---------------------------------------------------------------- public API
const char* to_string(MaterialFindingSeverity severity) {
    switch (severity) {
        case MaterialFindingSeverity::WARNING: return "WARNING";
        case MaterialFindingSeverity::SUSPICIOUS: return "SUSPICIOUS";
        case MaterialFindingSeverity::IMPOSSIBLE: return "IMPOSSIBLE";
    }
    return "UNKNOWN";
}

std::vector<std::string> MaterialValidator::check_codes() {
    return {MAT_PARSE, MAT_PROVENANCE, MAT_LOSS_EVAL, MAT_LOSS_OUT_OF_RANGE, MAT_LOSS_MONOTONIC, MAT_LOSS_ENVELOPE, MAT_LOSS_TEMPERATURE,
            MAT_LOSS_HYSTERESIS_BOUND, MAT_LOSS_POINTS_INCONSISTENT, MAT_LOSS_FACTOR, MAT_BSAT_CEILING, MAT_BSAT_CLASS, MAT_CURIE,
            MAT_CURIE_SATURATION, MAT_PERM};
}

std::vector<MaterialLossPoint> MaterialValidator::read_loss_points(const json& record) {
    std::vector<MaterialLossPoint> out;
    if (!record.contains("volumetricLosses") || !record["volumetricLosses"].is_object()) return out;
    for (const auto& [key, methods] : record["volumetricLosses"].items()) {
        if (!methods.is_array()) continue;
        for (const auto& m : methods) {
            if (m.is_array()) add_points(m, out);
        }
    }
    return out;
}

std::map<std::string, std::vector<MaterialLossPoint>> MaterialValidator::read_advanced_loss_points(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open advanced core materials file " + path);
    }
    std::map<std::string, std::vector<MaterialLossPoint>> out;
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        if (line.empty()) continue;
        if (lineNumber == 1 && line.rfind("version https://git-lfs", 0) == 0) {
            throw std::runtime_error(path + " is a git-LFS pointer, not the data: fetch it with git lfs pull");
        }
        json record = json::parse(line);
        if (!record.contains("name") || !record.contains("volumetricLosses")) continue;
        auto& points = out[record["name"].get<std::string>()];
        auto recordPoints = read_loss_points(record);
        points.insert(points.end(), recordPoints.begin(), recordPoints.end());
    }
    return out;
}

namespace {
MaterialVerdict validate_record(const json& record, const std::vector<MaterialLossPoint>& extraLossPoints, bool provenanceExpected);
}

MaterialVerdict MaterialValidator::validate(const json& record, const std::vector<MaterialLossPoint>& extraLossPoints) const {
    return validate_record(record, extraLossPoints, kSchemaHasProvenance);
}

namespace {
MaterialVerdict validate_record(const json& record, const std::vector<MaterialLossPoint>& extraLossPoints, bool provenanceExpected) {
    MaterialVerdict verdict;
    verdict.reference = record.is_object() && record.contains("name") && record["name"].is_string() ? record["name"].get<std::string>() : "<unnamed>";
    Ctx ctx{verdict.reference, &verdict.findings};

    // Provenance: per-record once the schema (generated CoreMaterial) or the record carries it.
    if (record.is_object() && record.contains("provenance")) {
        // present
    }
    else if (provenanceExpected) {
        emit(ctx, MAT_PROVENANCE, MaterialFindingSeverity::WARNING, 0, 0, "the record has no provenance (MAS RFC 0011)");
    }
    else {
        verdict.skipped.push_back(kProvenanceSchemaGap);
    }

    CoreMaterial material;
    try {
        material = CoreMaterial(record);
    }
    catch (const std::exception& e) {
        emit(ctx, MAT_PARSE, MaterialFindingSeverity::IMPOSSIBLE, 0, 0, std::string("MKF cannot read the record: ") + e.what());
        verdict.valid = false;
        return verdict;
    }

    CacheGuard guard(material.get_name());
    const MaterialClassBounds* cls = classify_material(record);
    if (cls) verdict.materialClass = cls->name;

    auto guarded = [&](const char* rule, auto&& check) {
        try {
            check();
        }
        catch (const std::exception& e) {
            verdict.skipped.push_back(std::string(rule) + ": MKF could not evaluate it: " + e.what());
        }
    };

    guarded(MAT_BSAT_CEILING, [&] { check_saturation(record, material, cls, ctx, verdict.skipped); });
    guarded(MAT_CURIE, [&] { check_curie(record, material, cls, ctx, verdict.skipped); });
    guarded(MAT_PERM, [&] { check_permeability(material, cls, ctx, verdict.skipped); });
    check_loss_factor(record, ctx);

    std::vector<MaterialLossPoint> points = MaterialValidator::read_loss_points(record);
    points.insert(points.end(), extraLossPoints.begin(), extraLossPoints.end());
    if (points.empty()) {
        verdict.skipped.push_back(std::string(MAT_LOSS_HYSTERESIS_BOUND) + ": no measured loss points");
        verdict.skipped.push_back(std::string(MAT_LOSS_POINTS_INCONSISTENT) + ": no measured loss points");
    }
    else {
        guarded(MAT_LOSS_HYSTERESIS_BOUND, [&] { check_hysteresis_bound(record, material, points, ctx, verdict.skipped); });
        check_points_consistency(points, ctx, verdict.skipped);
    }

    std::string modelName;
    guarded(MAT_LOSS_EVAL, [&] { check_loss_model(record, material, cls, ctx, verdict.skipped, modelName); });

    for (const auto& f : verdict.findings) {
        if (f.severity == MaterialFindingSeverity::IMPOSSIBLE) verdict.valid = false;
    }
    return verdict;
}
} // namespace

std::vector<MaterialVerdict> MaterialValidator::validate_catalogue_text(
    const std::string& ndjson, const std::map<std::string, std::vector<MaterialLossPoint>>& extraLossPoints) const {
    std::vector<MaterialVerdict> verdicts;
    // Once any record carries provenance, the catalogue has adopted it: a record without it is
    // then a per-record WARNING, even before the schema itself declares the field.
    const bool provenanceExpected = kSchemaHasProvenance || ndjson.find("\"provenance\"") != std::string::npos;
    std::istringstream in(ndjson);
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        if (lineNumber == 1 && line.rfind("version https://git-lfs", 0) == 0) {
            throw std::runtime_error("the core materials file is a git-LFS pointer, not the data");
        }
        json record;
        try {
            record = json::parse(line);
        }
        catch (const std::exception& e) {
            MaterialVerdict v;
            v.reference = "<line " + std::to_string(lineNumber) + ">";
            v.valid = false;
            v.findings.push_back(MaterialFinding{MAT_PARSE, MaterialFindingSeverity::IMPOSSIBLE, "coreMaterial", v.reference,
                                                 std::string("the line is not JSON: ") + e.what(), 0, 0});
            verdicts.push_back(v);
            continue;
        }
        std::vector<MaterialLossPoint> extra;
        if (record.contains("name") && record["name"].is_string()) {
            auto it = extraLossPoints.find(record["name"].get<std::string>());
            if (it != extraLossPoints.end()) extra = it->second;
        }
        verdicts.push_back(validate_record(record, extra, provenanceExpected));
    }
    return verdicts;
}

std::vector<MaterialVerdict> MaterialValidator::validate_catalogue(const std::string& coreMaterialsPath,
                                                                   const std::optional<std::string>& advancedCoreMaterialsPath) const {
    std::ifstream in(coreMaterialsPath);
    if (!in) {
        throw std::runtime_error("cannot open core materials file " + coreMaterialsPath);
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::map<std::string, std::vector<MaterialLossPoint>> extra;
    if (advancedCoreMaterialsPath) {
        extra = read_advanced_loss_points(*advancedCoreMaterialsPath);
    }
    return validate_catalogue_text(buffer.str(), extra);
}

std::vector<MaterialVerdict> MaterialValidator::validate_embedded_catalogue() const {
    auto fs = cmrc::data::get_filesystem();
    auto data = fs.open("MAS/data/core_materials.ndjson");
    return validate_catalogue_text(std::string(data.begin(), data.end()));
}

void to_json(json& j, const MaterialFinding& f) {
    j = json{{"code", f.code}, {"severity", to_string(f.severity)}, {"component", f.component}, {"reference", f.reference},
             {"message", f.message}, {"value", f.value}, {"threshold", f.threshold}};
}

void to_json(json& j, const MaterialVerdict& v) {
    j = json{{"reference", v.reference}, {"valid", v.valid}, {"materialClass", v.materialClass}, {"findings", v.findings}, {"skipped", v.skipped}};
}

json summarize_material_verdicts(const std::vector<MaterialVerdict>& verdicts) {
    json summary;
    size_t invalid = 0, suspicious = 0, withoutProvenanceField = 0, unclassified = 0;
    json byCode = json::object();
    json skippedByCode = json::object();
    json impossible = json::array();
    json suspiciousList = json::array();
    for (const auto& v : verdicts) {
        if (!v.valid) ++invalid;
        if (v.materialClass.empty()) ++unclassified;
        bool anySuspicious = false;
        for (const auto& f : v.findings) {
            std::string sev = to_string(f.severity);
            auto& counts = byCode[f.code];
            counts[sev] = (counts.contains(sev) ? counts[sev].get<size_t>() : 0) + 1;
            if (f.severity == MaterialFindingSeverity::IMPOSSIBLE) impossible.push_back(f);
            if (f.severity == MaterialFindingSeverity::SUSPICIOUS) {
                suspiciousList.push_back(f);
                anySuspicious = true;
            }
        }
        if (anySuspicious) ++suspicious;
        for (const auto& s : v.skipped) {
            if (s == kProvenanceSchemaGap) ++withoutProvenanceField;
            std::string code = s.substr(0, s.find(':'));
            skippedByCode[code] = skippedByCode.value(code, 0) + 1;
        }
    }
    summary["records"] = verdicts.size();
    summary["invalid"] = invalid;
    summary["recordsWithSuspicious"] = suspicious;
    summary["unclassified"] = unclassified;
    summary["findingsByCode"] = byCode;
    summary["skippedByCode"] = skippedByCode;
    if (withoutProvenanceField > 0) {
        summary["provenance"] = std::to_string(withoutProvenanceField) + "/" + std::to_string(verdicts.size()) +
                                " materials without provenance; schema has no field, RFC 0011";
    }
    summary["impossible"] = impossible;
    summary["suspicious"] = suspiciousList;
    return summary;
}

} // namespace OpenMagnetics
