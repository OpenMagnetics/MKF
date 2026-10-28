#include "physical_models/CoreMaterialCheck.h"

#include "constructive_models/Core.h"
#include "physical_models/AmplitudePermeability.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/ThermalResistance.h"
#include "processors/Inputs.h"
#include "support/NumberFormat.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <set>

namespace OpenMagnetics {

namespace {

// The refusal a design-path check raises: its code and its whole message.
struct Refusal {
    ErrorCode code;
    std::string what;
};

// Runs one design-path check. An MKF exception is the check's refusal and is returned; anything
// else (a bug, not a verdict on the data) propagates.
std::optional<Refusal> refusal_of(const std::function<void()>& check) {
    try {
        check();
        return std::nullopt;
    }
    catch (const OpenMagneticsException& exception) {
        return Refusal{exception.code(), exception.what()};
    }
}

struct RefusedInterval {
    double from;
    double to;
    Refusal refusal;
};

// The temperature where a check starts (or stops) refusing, between a temperature it accepts and one it
// refuses with `code`: bisection on the check itself, to a hundredth of the gap.
double refusal_boundary(const std::function<void(double)>& check, double accepted, double refused, ErrorCode code) {
    const double tolerance = std::abs(refused - accepted) * 1e-6;
    while (std::abs(refused - accepted) > tolerance) {
        double middle = (accepted + refused) / 2;
        auto refusal = refusal_of([&] { check(middle); });
        if (refusal && refusal->code == code) {
            refused = middle;
        }
        else {
            accepted = middle;
        }
    }
    return refused;
}

// Where inside [samples.front(), samples.back()] the check refuses, as intervals of one refusal code.
// The check runs at every sample and at the middle of every pair of neighbours; the edge of each
// refused interval next to an accepted sample is found by bisection.
std::vector<RefusedInterval> refused_intervals(std::vector<double> samples, const std::function<void(double)>& check) {
    std::sort(samples.begin(), samples.end());
    samples.erase(std::unique(samples.begin(), samples.end()), samples.end());
    std::vector<double> temperatures;
    for (size_t index = 0; index < samples.size(); ++index) {
        if (index > 0) {
            temperatures.push_back((samples[index - 1] + samples[index]) / 2);
        }
        temperatures.push_back(samples[index]);
    }
    std::vector<std::optional<Refusal>> outcomes;
    for (double temperature : temperatures) {
        outcomes.push_back(refusal_of([&] { check(temperature); }));
    }

    std::vector<RefusedInterval> intervals;
    size_t index = 0;
    while (index < temperatures.size()) {
        if (!outcomes[index]) {
            ++index;
            continue;
        }
        const ErrorCode code = outcomes[index]->code;
        size_t last = index;
        while (last + 1 < temperatures.size() && outcomes[last + 1] && outcomes[last + 1]->code == code) {
            ++last;
        }
        double from = temperatures[index];
        double to = temperatures[last];
        if (index > 0 && !outcomes[index - 1]) {
            from = refusal_boundary(check, temperatures[index - 1], temperatures[index], code);
        }
        if (last + 1 < temperatures.size() && !outcomes[last + 1]) {
            to = refusal_boundary(check, temperatures[last + 1], temperatures[last], code);
        }
        // The message of the refusal where the interval starts.
        auto refusal = refusal_of([&] { check(from); });
        intervals.push_back({from, to, refusal && refusal->code == code ? refusal.value() : outcomes[index].value()});
        index = last + 1;
    }
    return intervals;
}

std::vector<double> temperatures_of(const std::vector<BhCycleElement>& points) {
    std::vector<double> temperatures;
    for (const auto& point : points) {
        temperatures.push_back(point.get_temperature());
    }
    return temperatures;
}

std::vector<double> initial_permeability_temperatures(const CoreMaterial& material) {
    std::vector<double> temperatures;
    if (!std::holds_alternative<std::vector<PermeabilityPoint>>(material.get_permeability().get_initial())) {
        return temperatures;
    }
    for (const auto& point : InitialPermeability::get_only_temperature_dependent_points(material)) {
        if (point.get_temperature()) {
            temperatures.push_back(point.get_temperature().value());
        }
    }
    return temperatures;
}

std::optional<CoreMaterialSpan> span_of(const std::vector<double>& values) {
    if (values.empty()) {
        return std::nullopt;
    }
    auto [minimum, maximum] = std::minmax_element(values.begin(), values.end());
    return CoreMaterialSpan{*minimum, *maximum};
}

// The samples of `temperatures` inside [span.minimum, span.maximum], with the span's ends.
std::vector<double> samples_within(const std::vector<double>& temperatures, const CoreMaterialSpan& span) {
    std::vector<double> samples{span.minimum, span.maximum};
    for (double temperature : temperatures) {
        if (temperature > span.minimum && temperature < span.maximum) {
            samples.push_back(temperature);
        }
    }
    return samples;
}

// A span cut to below the material's Curie temperature, where the design prices a ferromagnetic material
// (above it every loss and permeability model is refused or extrapolated: CoreLosses ABT #1485).
std::optional<CoreMaterialSpan> below_curie(CoreMaterialSpan span, const CoreMaterial& material) {
    if (material.get_curie_temperature()) {
        const double curie = material.get_curie_temperature().value();
        if (span.minimum >= curie) {
            return std::nullopt;
        }
        // The design refuses at the Curie temperature itself; the span stops just short of it.
        span.maximum = std::min(span.maximum, std::nextafter(curie, span.minimum));
    }
    return span;
}

struct VolumetricLossesPointReading {
    double frequency;
    double magneticFluxDensityPeak;
    double temperature;
};

std::vector<VolumetricLossesPointReading> read_volumetric_losses_points(const CoreMaterial& material) {
    std::vector<VolumetricLossesPointReading> readings;
    for (const auto& point : CoreLossesModel::get_volumetric_losses_data(material)) {
        const auto& excitation = point.get_magnetic_flux_density();
        if (!excitation.get_magnetic_flux_density() || !excitation.get_magnetic_flux_density()->get_processed() ||
            !excitation.get_magnetic_flux_density()->get_processed()->get_peak()) {
            throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
                "Material " + material.get_name() + " has a volumetric-loss point at " + format_quantity(excitation.get_frequency(), "Hz") +
                " and " + format_significant(point.get_temperature()) + " C that states no peak flux density");
        }
        readings.push_back({excitation.get_frequency(), excitation.get_magnetic_flux_density()->get_processed()->get_peak().value(),
                            point.get_temperature()});
    }
    return readings;
}

std::string temperature_range_text(const RefusedInterval& interval) {
    if (interval.from == interval.to) {
        return "at " + format_significant(interval.from) + " C";
    }
    return "from " + format_significant(interval.from) + " C to " + format_significant(interval.to) + " C";
}

std::string span_text(const CoreMaterialSpan& span) {
    if (span.minimum == span.maximum) {
        return format_significant(span.minimum) + " C";
    }
    return format_significant(span.minimum) + " C to " + format_significant(span.maximum) + " C";
}

const std::string replaceAdvice = "Replace this copy of the material with the library's, or correct its data.";

void add_refusals(CoreMaterialCheck& check, const std::string& property, const std::string& subject, const CoreMaterialSpan& judged,
                  const std::vector<RefusedInterval>& intervals, std::optional<CoreMaterialSpan> frequency) {
    for (const auto& interval : intervals) {
        CoreMaterialProblem problem;
        problem.code = interval.refusal.code;
        problem.property = property;
        problem.engineMessage = interval.refusal.what;
        problem.temperature = CoreMaterialSpan{interval.from, interval.to};
        problem.frequency = frequency;
        const bool whole = interval.from == judged.minimum && interval.to == judged.maximum;
        problem.message = "The design calculation refuses " + subject + " of material " + check.material + " " +
                          temperature_range_text(interval) +
                          (whole ? ", all of the " + span_text(judged) + " it is judged over" : " (judged over " + span_text(judged) + ")") +
                          ". " + interval.refusal.what + ". " + replaceAdvice;
        check.problems.push_back(problem);
    }
}

void add_refusal(CoreMaterialCheck& check, const std::string& property, const std::string& subject, const Refusal& refusal) {
    CoreMaterialProblem problem;
    problem.code = refusal.code;
    problem.property = property;
    problem.engineMessage = refusal.what;
    problem.message = "The design calculation refuses " + subject + " of material " + check.material + ". " + refusal.what + ". " + replaceAdvice;
    check.problems.push_back(problem);
}

std::string frequency_span_text(const SteinmetzCoreLossesMethodRangeDatum& range) {
    return format_quantity(range.get_minimum_frequency().value(), "Hz") + " to " + format_quantity(range.get_maximum_frequency().value(), "Hz");
}

void check_steinmetz_ranges(CoreMaterialCheck& check, const CoreMaterial& material, const std::optional<CoreMaterialSpan>& lossSpan,
                            const std::vector<double>& tabulatedTemperatures) {
    const auto methods = CoreLossesModel::get_methods(material);
    if (std::find(methods.begin(), methods.end(), CoreLossesModels::STEINMETZ) == methods.end()) {
        check.notChecked.push_back("volumetricLosses: the material carries no Steinmetz ranges; its other loss data is judged only "
                                   "when a design evaluates it at an excitation");
        return;
    }
    // The span check every Steinmetz evaluation runs first (missing range frequencies are refused there).
    if (auto refusal = refusal_of([&] { CoreLossesModel::get_steinmetz_fitted_span(material); })) {
        add_refusal(check, "volumetricLosses", "the Steinmetz core-loss data", refusal.value());
        return;
    }
    if (!lossSpan) {
        check.notChecked.push_back("volumetricLosses: all the temperatures the material's data covers are at or above its Curie "
                                   "temperature, where no design prices its losses");
        return;
    }
    const auto ranges = CoreLossesModel::get_method_data(material, "Steinmetz").get_ranges().value();
    for (const auto& range : ranges) {
        auto samples = samples_within(tabulatedTemperatures, lossSpan.value());
        for (double root : CoreLossesModel::get_temperature_factor_roots(range)) {
            if (root > lossSpan->minimum && root < lossSpan->maximum) {
                samples.push_back(root);
            }
        }
        auto intervals = refused_intervals(samples, [&](double temperature) { CoreLossesModel::get_temperature_factor(range, temperature); });
        add_refusals(check, "volumetricLosses", "the core-loss data (Steinmetz range " + frequency_span_text(range) + ")", lossSpan.value(),
                     intervals, CoreMaterialSpan{range.get_minimum_frequency().value(), range.get_maximum_frequency().value()});
    }
}

} // namespace

CoreMaterialCheck check_core_material(const CoreMaterial& material) {
    CoreMaterialCheck check;
    check.material = material.get_name();

    // What the data covers.
    const auto permeabilityTemperatures = initial_permeability_temperatures(material);
    const auto saturationTemperatures = temperatures_of(material.get_saturation());
    const auto remanenceTemperatures = material.get_remanence() ? temperatures_of(material.get_remanence().value()) : std::vector<double>{};
    const auto coerciveForceTemperatures = material.get_coercive_force() ? temperatures_of(material.get_coercive_force().value()) : std::vector<double>{};
    std::vector<VolumetricLossesPointReading> lossPoints;
    if (auto refusal = refusal_of([&] { lossPoints = read_volumetric_losses_points(material); })) {
        add_refusal(check, "volumetricLosses", "the volumetric-loss points", refusal.value());
    }
    std::vector<double> lossPointTemperatures, lossPointFrequencies, lossPointFluxDensities;
    for (const auto& point : lossPoints) {
        lossPointTemperatures.push_back(point.temperature);
        lossPointFrequencies.push_back(point.frequency);
        lossPointFluxDensities.push_back(point.magneticFluxDensityPeak);
    }

    std::vector<double> tabulatedTemperatures;
    for (const std::vector<double>* temperatures : std::initializer_list<const std::vector<double>*>{&permeabilityTemperatures, &saturationTemperatures, &remanenceTemperatures, &coerciveForceTemperatures, &lossPointTemperatures}) {
        tabulatedTemperatures.insert(tabulatedTemperatures.end(), temperatures->begin(), temperatures->end());
    }
    auto temperatureSpan = span_of(tabulatedTemperatures);
    if (!temperatureSpan) {
        throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
            "Material " + material.get_name() + " tabulates no property at any temperature (initial permeability, saturation, remanence, "
            "coercive force or volumetric-loss points), so there is no temperature span to judge it over. Add its datasheet data");
    }
    check.coverage.temperature = temperatureSpan.value();

    std::vector<double> lossFrequencies = lossPointFrequencies;
    const auto methods = CoreLossesModel::get_methods(material);
    if (std::find(methods.begin(), methods.end(), CoreLossesModels::STEINMETZ) != methods.end()) {
        // A Steinmetz range without its frequencies has no span to cover: check_steinmetz_ranges reports
        // that refusal as a problem, so it is not recorded twice here.
        auto spanRefusal = refusal_of([&] {
            auto [minimum, maximum] = CoreLossesModel::get_steinmetz_fitted_span(material);
            lossFrequencies.push_back(minimum);
            lossFrequencies.push_back(maximum);
        });
        static_cast<void>(spanRefusal);
    }
    check.coverage.lossFrequency = span_of(lossFrequencies);
    check.coverage.lossMagneticFluxDensity = span_of(lossPointFluxDensities);

    // Core losses: a fitted model states no temperature span of its own, so it is judged over every
    // temperature the material's data covers, below its Curie point.
    check_steinmetz_ranges(check, material, below_curie(temperatureSpan.value(), material), tabulatedTemperatures);

    // Initial permeability, over its own table (outside it the design extrapolates or refuses by setting).
    if (auto span = span_of(permeabilityTemperatures)) {
        if (auto judged = below_curie(span.value(), material)) {
            auto intervals = refused_intervals(samples_within(permeabilityTemperatures, judged.value()),
                                               [&](double temperature) { InitialPermeability::get_initial_permeability(material, temperature); });
            add_refusals(check, "permeability", "the initial permeability", judged.value(), intervals, std::nullopt);
        }
    }
    else {
        check.notChecked.push_back("permeability: the material tabulates no temperature-dependent initial permeability");
    }

    // Saturation flux density, over its own table.
    if (auto span = span_of(saturationTemperatures)) {
        if (auto judged = below_curie(span.value(), material)) {
            auto intervals = refused_intervals(samples_within(saturationTemperatures, judged.value()),
                                               [&](double temperature) { Core::get_magnetic_flux_density_saturation(material, temperature, false); });
            add_refusals(check, "saturation", "the saturation flux density", judged.value(), intervals, std::nullopt);
        }
    }
    else {
        check.notChecked.push_back("saturation: the material carries no saturation data");
    }

    // Remanence: the hysteresis-loop resolver judges it against the saturation flux density, so over the
    // temperatures both tables cover.
    if (remanenceTemperatures.empty()) {
        check.notChecked.push_back("remanence: the material carries none; a design estimates it from its class");
    }
    else if (auto remanenceSpan = span_of(remanenceTemperatures), saturationSpan = span_of(saturationTemperatures); saturationSpan) {
        CoreMaterialSpan both{std::max(remanenceSpan->minimum, saturationSpan->minimum), std::min(remanenceSpan->maximum, saturationSpan->maximum)};
        std::optional<CoreMaterialSpan> judged = both.minimum <= both.maximum ? below_curie(both, material) : std::nullopt;
        if (judged) {
            std::vector<double> temperatures = remanenceTemperatures;
            temperatures.insert(temperatures.end(), saturationTemperatures.begin(), saturationTemperatures.end());
            auto intervals = refused_intervals(samples_within(temperatures, judged.value()),
                                               [&](double temperature) { resolve_core_material_remanence(material, temperature); });
            add_refusals(check, "remanence", "the remanence", judged.value(), intervals, std::nullopt);
        }
        else {
            check.notChecked.push_back("remanence: its table and the saturation table share no temperature below the Curie point");
        }
    }
    else {
        check.notChecked.push_back("remanence: the material carries no saturation data to judge it against");
    }

    // Coercive force, over its own table.
    if (auto span = span_of(coerciveForceTemperatures)) {
        if (auto judged = below_curie(span.value(), material)) {
            auto intervals = refused_intervals(samples_within(coerciveForceTemperatures, judged.value()),
                                               [&](double temperature) { resolve_core_material_coercive_force(material, temperature); });
            add_refusals(check, "coerciveForce", "the coercive force", judged.value(), intervals, std::nullopt);
        }
    }
    else {
        check.notChecked.push_back("coerciveForce: the material carries none; a design estimates it from its class");
    }

    // Heat conductivity: a value the material carries must be usable (absence walks the estimate tiers).
    if (material.get_heat_conductivity()) {
        if (auto refusal = refusal_of([&] { resolve_core_material_heat_conductivity(material); })) {
            add_refusal(check, "heatConductivity", "the heat conductivity", refusal.value());
        }
    }
    else {
        check.notChecked.push_back("heatConductivity: the material carries none; a design estimates it from its class");
    }

    check.valid = check.problems.empty();
    return check;
}

namespace {

CoreMaterialLossValue price_volumetric_losses(CoreLosses& coreLosses, const CoreMaterial& material, double frequency, double magneticFluxDensityPeak,
                                              double temperature) {
    auto waveform = Inputs::create_waveform(WaveformLabel::SINUSOIDAL, 2 * magneticFluxDensityPeak, frequency);
    SignalDescriptor magneticFluxDensity;
    magneticFluxDensity.set_waveform(waveform);
    magneticFluxDensity.set_processed(Inputs::calculate_processed_data(waveform, frequency));
    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    excitation.set_magnetic_flux_density(magneticFluxDensity);

    CoreMaterialLossValue value;
    if (auto refusal = refusal_of([&] { value.volumetricLosses = coreLosses.get_core_volumetric_losses(material, excitation, temperature); })) {
        value.volumetricLosses = std::nullopt;
        value.code = refusal->code;
        value.message = refusal->what;
    }
    return value;
}

} // namespace

CoreMaterialLossComparison compare_core_material_volumetric_losses(const CoreMaterial& reference, const CoreMaterial& other) {
    CoreMaterialLossComparison comparison;
    comparison.reference = reference.get_name();
    comparison.other = other.get_name();

    std::vector<VolumetricLossesPointReading> points = read_volumetric_losses_points(reference);
    if (!points.empty()) {
        comparison.pointsSource = "the " + std::to_string(points.size()) + " volumetric-loss points of " + reference.get_name();
    }
    else {
        const auto methods = CoreLossesModel::get_methods(reference);
        if (std::find(methods.begin(), methods.end(), CoreLossesModels::STEINMETZ) == methods.end()) {
            throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
                "Material " + reference.get_name() + " carries neither volumetric-loss points nor Steinmetz ranges, so its data states no "
                "frequency to compare its losses at");
        }
        if (reference.get_saturation().empty()) {
            throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
                "Material " + reference.get_name() + " carries no saturation data, so its data states no flux density to compare its losses at");
        }
        // Validates every range's frequencies (missing ones throw).
        CoreLossesModel::get_steinmetz_fitted_span(reference);
        const auto ranges = CoreLossesModel::get_method_data(reference, "Steinmetz").get_ranges().value();
        std::set<double> temperatures;
        for (const auto& saturation : reference.get_saturation()) {
            if (!reference.get_curie_temperature() || saturation.get_temperature() < reference.get_curie_temperature().value()) {
                temperatures.insert(saturation.get_temperature());
            }
        }
        if (temperatures.empty()) {
            throw MaterialException(ErrorCode::MATERIAL_DATA_MISSING,
                "Material " + reference.get_name() + " tabulates its saturation flux density only at or above its Curie temperature, so its "
                "data states no temperature to compare its losses at");
        }
        for (double temperature : temperatures) {
            const double saturation = Core::get_magnetic_flux_density_saturation(reference, temperature, false);
            for (const auto& range : ranges) {
                const double frequency = std::sqrt(range.get_minimum_frequency().value() * range.get_maximum_frequency().value());
                for (double fraction : {0.5, 0.25, 0.125}) {
                    points.push_back({frequency, saturation * fraction, temperature});
                }
            }
        }
        // Steinmetz coefficients state a frequency span but no flux-density span, so the flux densities are a
        // chosen spread (1/2, 1/4, 1/8 of the saturation flux density), not points of the data; pointsSource says so.
        comparison.pointsSource = "the temperatures " + reference.get_name() + " tabulates its saturation flux density at (below its Curie "
                                  "point), at the geometric centre of each of its " + std::to_string(ranges.size()) +
                                  " Steinmetz ranges, and at 1/2, 1/4 and 1/8 of that saturation flux density: a chosen spread, because "
                                  "Steinmetz coefficients state no flux-density span";
    }

    CoreLosses coreLosses;
    for (const auto& point : points) {
        CoreMaterialLossComparisonPoint comparisonPoint;
        comparisonPoint.frequency = point.frequency;
        comparisonPoint.magneticFluxDensityPeak = point.magneticFluxDensityPeak;
        comparisonPoint.temperature = point.temperature;
        comparisonPoint.reference = price_volumetric_losses(coreLosses, reference, point.frequency, point.magneticFluxDensityPeak, point.temperature);
        comparisonPoint.other = price_volumetric_losses(coreLosses, other, point.frequency, point.magneticFluxDensityPeak, point.temperature);
        if (!comparisonPoint.reference.volumetricLosses) {
            ++comparison.pointsRefusedByReference;
        }
        if (!comparisonPoint.other.volumetricLosses) {
            ++comparison.pointsRefusedByOther;
        }
        if (comparisonPoint.reference.volumetricLosses && comparisonPoint.other.volumetricLosses && comparisonPoint.reference.volumetricLosses.value() > 0) {
            const double difference = (comparisonPoint.other.volumetricLosses.value() - comparisonPoint.reference.volumetricLosses.value()) /
                                      comparisonPoint.reference.volumetricLosses.value();
            comparisonPoint.relativeDifference = difference;
            if (!comparison.maximumAbsoluteRelativeDifference || std::abs(difference) > comparison.maximumAbsoluteRelativeDifference.value()) {
                comparison.maximumAbsoluteRelativeDifference = std::abs(difference);
            }
        }
        comparison.points.push_back(comparisonPoint);
    }
    return comparison;
}

void to_json(json& j, const CoreMaterialSpan& span) {
    j = json{{"minimum", span.minimum}, {"maximum", span.maximum}};
}

void to_json(json& j, const CoreMaterialProblem& problem) {
    j = json{{"code", to_string(problem.code)}, {"property", problem.property}, {"message", problem.message}, {"engineMessage", problem.engineMessage}};
    j["temperature"] = problem.temperature ? json(problem.temperature.value()) : json(nullptr);
    j["frequency"] = problem.frequency ? json(problem.frequency.value()) : json(nullptr);
}

void to_json(json& j, const CoreMaterialCoverage& coverage) {
    j = json{{"temperature", coverage.temperature}};
    j["lossFrequency"] = coverage.lossFrequency ? json(coverage.lossFrequency.value()) : json(nullptr);
    j["lossMagneticFluxDensity"] = coverage.lossMagneticFluxDensity ? json(coverage.lossMagneticFluxDensity.value()) : json(nullptr);
}

void to_json(json& j, const CoreMaterialCheck& check) {
    j = json{{"material", check.material}, {"valid", check.valid}, {"problems", check.problems}, {"coverage", check.coverage},
             {"notChecked", check.notChecked}};
}

void to_json(json& j, const CoreMaterialLossValue& value) {
    j = json::object();
    j["volumetricLosses"] = value.volumetricLosses ? json(value.volumetricLosses.value()) : json(nullptr);
    j["code"] = value.code ? json(to_string(value.code.value())) : json(nullptr);
    j["message"] = value.message ? json(value.message.value()) : json(nullptr);
}

void to_json(json& j, const CoreMaterialLossComparisonPoint& point) {
    j = json{{"frequency", point.frequency}, {"magneticFluxDensityPeak", point.magneticFluxDensityPeak}, {"temperature", point.temperature},
             {"reference", point.reference}, {"other", point.other}};
    j["relativeDifference"] = point.relativeDifference ? json(point.relativeDifference.value()) : json(nullptr);
}

void to_json(json& j, const CoreMaterialLossComparison& comparison) {
    j = json{{"reference", comparison.reference}, {"other", comparison.other}, {"pointsSource", comparison.pointsSource},
             {"points", comparison.points}, {"pointsRefusedByReference", comparison.pointsRefusedByReference},
             {"pointsRefusedByOther", comparison.pointsRefusedByOther}};
    j["maximumAbsoluteRelativeDifference"] =
        comparison.maximumAbsoluteRelativeDifference ? json(comparison.maximumAbsoluteRelativeDifference.value()) : json(nullptr);
}

} // namespace OpenMagnetics
