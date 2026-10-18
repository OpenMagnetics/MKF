#include "physical_models/ComplexPermeability.h"
#include "physical_models/InitialPermeability.h"

#include "support/Utils.h"
#include <algorithm>
#include <limits>
#include <math.h>
#include "support/Exceptions.h"
#include "support/Logger.h"
#include <tuple>


namespace OpenMagnetics {

std::pair<double, double> ComplexPermeability::get_complex_permeability(std::string coreMaterialName, double frequency) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_complex_permeability(coreMaterial, frequency);
}

ComplexPermeabilityData ComplexPermeability::calculate_complex_permeability_from_frequency_dependent_initial_permeability(std::string coreMaterialName) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return calculate_complex_permeability_from_frequency_dependent_initial_permeability(coreMaterial);
}
double ComplexPermeability::infer_F_mu_from_delta_FL(double deltaFL_95_90) {
    // Visual readings from Fig 1 (Mueller PCIM 2024). The three plotted
    // curves (4µ / 14µ / 205µ initial) overlap closely; this is the
    // representative middle curve. Values in (ΔF, F_µ).
    static const std::vector<std::pair<double, double>> table = {
        {0.31,    1.0},
        {0.32,    2.0},
        {0.36,    5.0},
        {0.43,   10.0},
        {0.55,   25.0},
        {0.65,   50.0},
        {0.69,  100.0},
        {0.73,  300.0},
        {0.755, 1000.0},
    };
    if (!std::isfinite(deltaFL_95_90)) return 1.0;
    if (deltaFL_95_90 <= table.front().first) return table.front().second;
    if (deltaFL_95_90 >= table.back().first)  return table.back().second;
    for (size_t i = 1; i < table.size(); ++i) {
        if (deltaFL_95_90 <= table[i].first) {
            double x0 = table[i - 1].first;
            double x1 = table[i].first;
            double y0 = std::log10(table[i - 1].second);
            double y1 = std::log10(table[i].second);
            double t = (deltaFL_95_90 - x0) / (x1 - x0);
            return std::pow(10.0, y0 + t * (y1 - y0));
        }
    }
    return 1.0;
}

ComplexPermeabilityData ComplexPermeability::calculate_complex_permeability_from_frequency_dependent_initial_permeability(CoreMaterial coreMaterial) {
    // Closed-form sheet expression (paper eq. 9 with hysteresis phase θh = 0)
    // expanded into real and imaginary parts, anchored at the frequency
    // where the bare-material µ' has dropped to ~67.78 % of its DC value
    // (the value of the closed-form at normalized frequency = 1).
    // calculate_frequency_for_initial_permeability_drop takes the drop FRACTION
    // (targets reference * (1 - drop)), so anchoring where u' has dropped TO
    // 67.78% means passing 1 - 0.6778 = 0.3222. Passing 0.6778 anchored at
    // u' = 0.32*u_DC, far too high in frequency.
    double frequencyFor67Point78Drop = InitialPermeability::calculate_frequency_for_initial_permeability_drop(coreMaterial, 0.3222);
    double initialPermeability = InitialPermeability::get_initial_permeability(coreMaterial);

    // Distributed-air-gap correction (paper eqs 10-12). Powder cores have
    // the magnetic particles suspended in a non-magnetic binder; the bare
    // particle permeability µ_material is F_µ × the pressed-core µ_core.
    // Eq 10 wraps the sheet expression in a gap divider so the predicted
    // µ_core(f) falloff is much more gradual than the homogeneous case.
    // For ferrite (F_µ → 1, F_g → 1) the divider is identity and the
    // result matches the previous homogeneous behavior exactly.
    double F_mu = 1.0;
    double F_g = 1.0;
    double f_L90 = InitialPermeability::calculate_frequency_for_initial_permeability_drop(coreMaterial, 0.10);
    double f_L95 = InitialPermeability::calculate_frequency_for_initial_permeability_drop(coreMaterial, 0.05);
    if (std::isfinite(f_L90) && std::isfinite(f_L95) && f_L90 > 0 && f_L90 > f_L95) {
        double deltaFL = (f_L90 - f_L95) / f_L90;                     // eq 13
        F_mu = infer_F_mu_from_delta_FL(deltaFL);                     // Fig 1
        if (F_mu > 1.0 && initialPermeability > 1.0) {
            // Eq 12 simplified: F_g = F_µ·(µ′_core − 1) / (F_µ·µ′_core − 1)
            F_g = F_mu * (initialPermeability - 1.0)
                / (F_mu * initialPermeability - 1.0);
        }
    }
    double oneMinusFg = 1.0 - F_g;
    double materialPermeabilityScale = F_mu * initialPermeability;

    // Normalized span, and why it reaches so far past the anchor: get_complex_permeability
    // answers only inside this table's frequency range (it throws outside it), and it used to
    // hold the last tabulated point instead — a constant permeability, not a falling one.
    // At 0.01..100x the anchor, a material anchored at 18.8 kHz
    // (Nanoperm 80000) tabulated only 188 Hz..1.88 MHz, and every sweep above 1.88 MHz
    // read back the same mu' and mu''. That reported 23,565 at 25.6 MHz where the
    // material's own initial-permeability table says 156 — 151x high — and pushed the
    // modelled self-resonance of a common-mode choke from ~50 MHz down to 4.33 MHz
    // (ABT #843). The closed form is defined for any normalized frequency; only the
    // tabulation was short. Keep the same ~10 points per decade.
    auto normalizedFrequencies = logarithmic_spaced_array(0.01, 1e5, 70);
    // The table must also reach down to where the material's own initial-permeability data
    // starts. get_complex_permeability answers only inside this table and throws outside it, and
    // 0.01x the anchor lies above the data for materials anchored at MHz (Edge 26: 340 kHz), so
    // a 100 kHz request that the data covers would have no table. Extend the grid downwards at
    // the same logarithmic step, ending exactly at the data's lowest frequency; the points above
    // 0.01x the anchor are unchanged.
    {
        double lowestDataFrequency = std::numeric_limits<double>::max();
        for (const auto& point : InitialPermeability::get_only_frequency_dependent_points(coreMaterial)) {
            if (point.get_frequency()) {
                lowestDataFrequency = std::min(lowestDataFrequency, point.get_frequency().value());
            }
        }
        // A material whose mu(f) is a FORMULA (Poco/Magnetics/Micrometals frequency factors: no
        // tabulated points) is defined down to DC, so its table must reach low frequencies too.
        // Stopping at 0.01x the anchor refused GPC 26 (anchor ~1 MHz) at 1 kHz with "complex
        // permeability data only from 10 kHz" — data it never had; the limit was this grid's
        // (user report 2026-10-09). Down there the closed form gives mu' -> mu_i and mu'' -> 0.
        // A formula material stores ONE initial-permeability point (plus its modifiers);
        // get_only_frequency_dependent_points samples that formula from the measurement frequency
        // (10 kHz) up, so those samples are not data and must not set the table's lower end.
        constexpr double formulaLowestFrequency = 1;  // Hz
        bool isFormulaMaterial = std::holds_alternative<PermeabilityPoint>(coreMaterial.get_permeability().get_initial());
        bool hasTabulatedPoints = !isFormulaMaterial && lowestDataFrequency < std::numeric_limits<double>::max();
        double tableLowestFrequency = hasTabulatedPoints ? lowestDataFrequency : formulaLowestFrequency;
        if (frequencyFor67Point78Drop > 0) {
            double lowestNormalizedFrequency = tableLowestFrequency / frequencyFor67Point78Drop;
            double logarithmicStep = normalizedFrequencies[1] / normalizedFrequencies[0];
            std::vector<double> lowerNormalizedFrequencies;
            double normalizedFrequency = normalizedFrequencies.front() / logarithmicStep;
            while (normalizedFrequency > lowestNormalizedFrequency) {
                lowerNormalizedFrequencies.push_back(normalizedFrequency);
                normalizedFrequency /= logarithmicStep;
            }
            if (lowestNormalizedFrequency < normalizedFrequencies.front()) {
                lowerNormalizedFrequencies.push_back(lowestNormalizedFrequency);
            }
            std::reverse(lowerNormalizedFrequencies.begin(), lowerNormalizedFrequencies.end());
            normalizedFrequencies.insert(normalizedFrequencies.begin(), lowerNormalizedFrequencies.begin(), lowerNormalizedFrequencies.end());
        }
    }
    // Where the material's own initial-permeability curve actually ends. Past that point
    // the interpolator turns back UP (Nanoperm 80000: 156 at its last point, 25.6 MHz, then
    // 451 at 50 MHz and 1,057 at 100 MHz), and permeability does not recover after roll-off.
    // Knowing the real end lets us follow the data up to it and extrapolate beyond it,
    // instead of guessing from the shape where the data stopped (ABT #843).
    double lastTabulatedFrequency = 0;
    double lastTabulatedPermeability = 0;
    double beyondTheDataLogSlope = 0;
    {
        auto initialPermeabilityData = coreMaterial.get_permeability().get_initial();
        if (std::holds_alternative<std::vector<PermeabilityPoint>>(initialPermeabilityData)) {
            auto points = std::get<std::vector<PermeabilityPoint>>(initialPermeabilityData);
            std::sort(points.begin(), points.end(), [](const PermeabilityPoint& a, const PermeabilityPoint& b) {
                return a.get_frequency().value_or(0) < b.get_frequency().value_or(0);
            });
            // Take the slope from the material's own LAST TWO points, not from the sampled
            // grid: near the end the interpolator has already begun to turn, so a slope read
            // off the samples can come out positive and send the extrapolation back up.
            if (points.size() >= 2) {
                const auto& last = points[points.size() - 1];
                const auto& secondToLast = points[points.size() - 2];
                if (last.get_frequency() && secondToLast.get_frequency() && last.get_value() > 0 && secondToLast.get_value() > 0 && *last.get_frequency() > *secondToLast.get_frequency()) {
                    lastTabulatedFrequency = *last.get_frequency();
                    lastTabulatedPermeability = last.get_value();
                    beyondTheDataLogSlope = log(last.get_value() / secondToLast.get_value()) / log(*last.get_frequency() / *secondToLast.get_frequency());
                    // Permeability past roll-off falls; never let the tail climb.
                    beyondTheDataLogSlope = std::min(beyondTheDataLogSlope, 0.0);
                }
            }
        }
    }
    std::vector<PermeabilityPoint> real;
    std::vector<PermeabilityPoint> imaginary;

    for (auto normalizedFrequency : normalizedFrequencies) {
        double sqrtFn = sqrt(normalizedFrequency);
        double s2 = sin(2 * sqrtFn);
        double sh2 = sinh(2 * sqrtFn);
        double c2 = cos(2 * sqrtFn);
        double ch2 = cosh(2 * sqrtFn);
        double denom = 2 * sqrtFn * (c2 + ch2);
        double muMatRealNormalized = (s2 + sh2) / denom;
        double muMatImagNormalized = -(s2 - sh2) / denom;

        // Bare magnetic material (paper eq. 9, real and imaginary parts).
        double muMatReal = materialPermeabilityScale * muMatRealNormalized;
        double muMatImag = materialPermeabilityScale * muMatImagNormalized;

        // Apply the gap divider (paper eq. 10) as a complex division:
        //   µ_core = µ_material / (F_g + (1 − F_g) · µ_material).
        // For F_g = 1 (homogeneous case) the divider collapses to 1 and
        // µ_core = µ_material directly.
        double denomReal = F_g + oneMinusFg * muMatReal;
        double denomImag = oneMinusFg * muMatImag;
        double denomMagnitudeSquared = denomReal * denomReal + denomImag * denomImag;
        double muCoreReal = (muMatReal * denomReal + muMatImag * denomImag) / denomMagnitudeSquared;
        double muCoreImag = (muMatImag * denomReal - muMatReal * denomImag) / denomMagnitudeSquared;

        double permeabilityPointFrequency = normalizedFrequency * frequencyFor67Point78Drop;

        // Put mu' back on the material's OWN measured curve. The closed form above supplies
        // the SHAPE of the loss (the mu''/mu' angle) and the distributed-gap correction, but
        // its roll-off is the sheet solution's ~1/sqrt(f), while a real nanocrystalline
        // material falls closer to ~1/f. Left alone it reported 7,375 at 25.6 MHz where
        // Nanoperm 80000's own initial-permeability table says 156 — still 47x high after
        // the tabulation range was widened. Rescale both parts by the ratio of the measured
        // mu'(f) to the modelled one, which pins mu' to the data and carries mu'' along at
        // the modelled loss angle (ABT #843). Materials without a frequency-dependent table
        // never reach this function, so there is always a curve to read.
        // The table ENDS where the material's own frequency model stops describing a material.
        // A formula mu(f) is a fit with a valid range: XFlux 60's Magnetics frequencyFactor quartic
        // reads 59.9 at 100 kHz, 58.5 at 1 MHz, 30.7 at 20 MHz and 0.83 at 26.6 MHz, and goes
        // negative after that. A relative permeability below the vacuum's is not a material value
        // (get_initial_permeability refuses it, ABT #1740), so there is no table point there: the
        // points before it stand, and get_complex_permeability answers inside them and throws
        // above them, like any other end of the tabulated range. Asking the checked getter here
        // instead made that one out-of-range point throw away the whole table, and the material
        // had no complex permeability at 100 kHz either (ABT #1790). The model value is read as
        // it stands only to find that end; it is never used below 1, clamped or floored.
        // A point also needs the model 5 % above it (the local slope of the Kramers-Kronig loss
        // angle below), so the table ends at the first frequency where either read is below 1.
        // Past the last TABULATED point the model is not read at all (see below).
        constexpr double localSlopeRelativeStep = 1.05;
        bool beyondTheTabulatedData = lastTabulatedFrequency > 0 && permeabilityPointFrequency > lastTabulatedFrequency;
        double measuredPermeability = 0;
        double permeabilityAbove = 0;
        if (!beyondTheTabulatedData) {
            measuredPermeability = InitialPermeability::get_initial_permeability_model_value(coreMaterial, std::nullopt, std::nullopt, permeabilityPointFrequency);
            permeabilityAbove = InitialPermeability::get_initial_permeability_model_value(coreMaterial, std::nullopt, std::nullopt, permeabilityPointFrequency * localSlopeRelativeStep);
            if (measuredPermeability < 1 || permeabilityAbove < 1) {
                break;
            }
        }

        // Past its last tabulated point the initial-permeability interpolator turns back UP:
        // Nanoperm 80000 reads 156 at 25.6 MHz (its last point), then 451 at 50 MHz and
        // 1,057 at 100 MHz. Permeability does not recover after roll-off, so a rise means we
        // have run off the end of the data, not that the material got more permeable. Carry
        // on with the last real log-log slope instead of following the rise (ABT #843).
        if (beyondTheTabulatedData) {
            measuredPermeability = lastTabulatedPermeability * pow(permeabilityPointFrequency / lastTabulatedFrequency, beyondTheDataLogSlope);
        }

        if (measuredPermeability > 0 && muCoreReal > 0) {
            // ABT #847: the sheet solution's loss angle is anchored by resistivity and
            // dimensions, which for a high-permeability nanocrystalline material lands the
            // eddy knee orders of magnitude above the real relaxation knee — so mu'' came
            // out ~38x too small right where the material is relaxing (Nanoperm 8000 at
            // 100 kHz: model 105, Kramers-Kronig from its own table ~3,700). For a local
            // power law |mu| ~ f^-n, Kramers-Kronig fixes the loss angle at delta = n*pi/2:
            // read n from the measured curve itself and take the LARGER of the two angles —
            // relaxation loss where the table rolls off, sheet eddy loss where it is flat
            // (a flat table says nothing about eddy loss, which the closed form does).
            double sheetAngle = std::atan2(std::max(muCoreImag, 0.0), muCoreReal);
            double localSlope = 0;
            if (beyondTheTabulatedData) {
                localSlope = beyondTheDataLogSlope;
            }
            else {
                double permeabilityBelow = InitialPermeability::get_initial_permeability(coreMaterial, std::nullopt, std::nullopt, permeabilityPointFrequency / localSlopeRelativeStep);
                if (permeabilityAbove > 0 && permeabilityBelow > 0) {
                    localSlope = log(permeabilityAbove / permeabilityBelow) / (2 * log(localSlopeRelativeStep));
                }
            }
            // n in [0, 0.95]: negative slopes only (a rising read is interpolator artefact),
            // capped below 1 so the angle stays short of pi/2 and mu' stays positive.
            double kramersKronigAngle = std::clamp(-localSlope, 0.0, 0.95) * std::numbers::pi / 2.0;
            double lossAngle = std::max(sheetAngle, kramersKronigAngle);
            // The table is the measured |mu|; split it by the loss angle so the impedance
            // magnitude stays pinned to the data.
            muCoreReal = measuredPermeability * std::cos(lossAngle);
            muCoreImag = measuredPermeability * std::sin(lossAngle);
        }

        PermeabilityPoint realPermeabilityPoint;
        realPermeabilityPoint.set_frequency(permeabilityPointFrequency);
        realPermeabilityPoint.set_value(muCoreReal);
        PermeabilityPoint imaginaryPermeabilityPoint;
        imaginaryPermeabilityPoint.set_frequency(permeabilityPointFrequency);
        imaginaryPermeabilityPoint.set_value(muCoreImag);
        real.push_back(realPermeabilityPoint);
        imaginary.push_back(imaginaryPermeabilityPoint);
    }

    ComplexPermeabilityData complexPermeabilityData;
    complexPermeabilityData.set_real(real);
    complexPermeabilityData.set_imaginary(imaginary);
    return complexPermeabilityData;
}


namespace {
void combine_hash(std::size_t& seed, std::size_t value) {
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

// The caches are per thread; a catalogue clear or reload on any thread empties them on this one.
void drop_complex_permeability_caches_if_the_catalogue_changed() {
    const uint64_t currentCatalogueGeneration = catalogueGeneration.load(std::memory_order_relaxed);
    if (complexPermeabilityCachesCatalogueGeneration == currentCatalogueGeneration) {
        return;
    }
    complexPermeabilityRealInterps.clear();
    complexPermeabilityImaginaryInterps.clear();
    complexPermeabilityRealFrequencySpans.clear();
    complexPermeabilityImaginaryFrequencySpans.clear();
    complexPermeabilityEstimates.clear();
    complexPermeabilityEstimateRefusals.clear();
    complexPermeabilityCachesCatalogueGeneration = currentCatalogueGeneration;
}
} // namespace

std::string ComplexPermeability::cache_key(const CoreMaterial& coreMaterial) {
    // Everything the curves are built from: the tabulated complex permeability when there is one,
    // otherwise the initial permeability they are derived from (the initial permeability is part
    // of the key either way: the estimates of the untabulated range read it too).
    std::size_t seed = 0;
    const auto& permeability = coreMaterial.get_permeability();
    combine_hash(seed, InitialPermeability::hash_permeability(permeability.get_initial()));
    const auto complex = permeability.get_complex();
    if (complex) {
        combine_hash(seed, 1);
        combine_hash(seed, InitialPermeability::hash_permeability(complex->get_real()));
        combine_hash(seed, InitialPermeability::hash_permeability(complex->get_imaginary()));
    }
    return coreMaterial.get_name() + "#" + std::to_string(seed);
}

std::string ComplexPermeability::ensure_interpolators(const CoreMaterial& coreMaterial) {
    // Cached per material: advisers that scan thousands of cores share only a few materials,
    // and rebuilding the derived curves (calculate_frequency_for_initial_permeability_drop runs
    // O(40) pow() per call, three times per material) dominated DMC/CMC core selection.
    drop_complex_permeability_caches_if_the_catalogue_changed();
    const std::string& materialName = coreMaterial.get_name();
    std::string cacheKey = cache_key(coreMaterial);
    if (complexPermeabilityRealInterps.contains(cacheKey) &&
        complexPermeabilityImaginaryInterps.contains(cacheKey)) {
        return cacheKey;
    }

    ComplexPermeabilityData complexPermeabilityData;
    if (!coreMaterial.get_permeability().get_complex()) {
        if (InitialPermeability::has_frequency_dependency(coreMaterial)) {
            complexPermeabilityData = calculate_complex_permeability_from_frequency_dependent_initial_permeability(coreMaterial);
        }
        else {
            throw MaterialDataMissingException(materialName, "Complex permeability");
        }
    }
    else {
        complexPermeabilityData = coreMaterial.get_permeability().get_complex().value();
    }

    auto realPart = complexPermeabilityData.get_real();
    auto imaginaryPart = complexPermeabilityData.get_imaginary();

    if (!std::holds_alternative<std::vector<PermeabilityPoint>>(realPart) ||
        !std::holds_alternative<std::vector<PermeabilityPoint>>(imaginaryPart)) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Complex permeability data is not in expected format for " + materialName);
    }

    auto buildInterpolator = [&](std::vector<PermeabilityPoint> points, const std::string& part,
                                 std::map<std::string, tk::spline>& interps,
                                 std::map<std::string, std::pair<double, double>>& spans) {
        for (const auto& point : points) {
            if (!point.get_frequency()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Complex permeability " + part + " point without frequency in " + materialName);
            }
        }
        std::sort(points.begin(), points.end(), [](const PermeabilityPoint& b1, const PermeabilityPoint& b2) {
            return b1.get_frequency().value() < b2.get_frequency().value();
        });
        std::vector<double> x, y;
        for (const auto& point : points) {
            if (x.empty() || fabs(*point.get_frequency() - x.back()) > 1e-9) {
                x.push_back(*point.get_frequency());
                y.push_back(point.get_value());
            }
        }
        if (x.size() < 2) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Not enough complex permeability " + part + " data for " + materialName);
        }
        interps[cacheKey] = tk::spline(x, y, tk::spline::cspline_hermite);
        spans[cacheKey] = {x.front(), x.back()};
    };

    buildInterpolator(std::get<std::vector<PermeabilityPoint>>(realPart), "real", complexPermeabilityRealInterps, complexPermeabilityRealFrequencySpans);
    buildInterpolator(std::get<std::vector<PermeabilityPoint>>(imaginaryPart), "imaginary", complexPermeabilityImaginaryInterps, complexPermeabilityImaginaryFrequencySpans);
    return cacheKey;
}

std::pair<double, double> ComplexPermeability::get_frequency_range(std::string coreMaterialName) {
    return get_frequency_range(Core::resolve_material(coreMaterialName));
}

std::pair<double, double> ComplexPermeability::get_frequency_range(CoreMaterial coreMaterial) {
    return frequency_range_for_key(ensure_interpolators(coreMaterial), coreMaterial.get_name());
}

std::pair<double, double> ComplexPermeability::frequency_range_for_key(const std::string& cacheKey, const std::string& materialName) {
    const auto& realSpan = complexPermeabilityRealFrequencySpans.at(cacheKey);
    const auto& imaginarySpan = complexPermeabilityImaginaryFrequencySpans.at(cacheKey);
    // Both parts are needed for every value returned, so the usable range is where both exist.
    double minimum = std::max(realSpan.first, imaginarySpan.first);
    double maximum = std::min(realSpan.second, imaginarySpan.second);
    if (minimum > maximum) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Complex permeability of " + materialName +
                                    ": the real part is tabulated over " + std::to_string(realSpan.first) + "-" + std::to_string(realSpan.second) +
                                    " Hz and the imaginary part over " + std::to_string(imaginarySpan.first) + "-" + std::to_string(imaginarySpan.second) +
                                    " Hz; the two do not overlap");
    }
    return {minimum, maximum};
}

std::pair<double, double> ComplexPermeability::get_tabulated_complex_permeability(const CoreMaterial& coreMaterial, double frequency) {
    // Inside the tabulated range only. Outside it there is no data: holding the last point
    // (as this used to) read K081's stale 1 MHz point at 350 MHz, and a mu' floor of 1 hid a
    // tabulated value below vacuum instead of reporting it. Both are gone; a mu' below 1 in
    // the data is returned as it is.
    const std::string cacheKey = ensure_interpolators(coreMaterial);
    auto [minimumFrequency, maximumFrequency] = frequency_range_for_key(cacheKey, coreMaterial.get_name());
    if (!(frequency >= minimumFrequency && frequency <= maximumFrequency)) {
        throw ComplexPermeabilityFrequencyOutOfRangeException(coreMaterial.get_name(), frequency, minimumFrequency, maximumFrequency);
    }
    double complexPermeabilityRealValue = complexPermeabilityRealInterps.at(cacheKey)(frequency);
    if (std::isnan(complexPermeabilityRealValue)) {
        throw NaNResultException("complex Permeability real part must be a number, not NaN");
    }
    double complexPermeabilityImaginaryValue = complexPermeabilityImaginaryInterps.at(cacheKey)(frequency);
    if (std::isnan(complexPermeabilityImaginaryValue)) {
        throw NaNResultException("complex Permeability imaginary part must be a number, not NaN");
    }
    return {complexPermeabilityRealValue, complexPermeabilityImaginaryValue};
}

std::shared_ptr<const ComplexPermeabilityEstimate> ComplexPermeability::find_or_build_estimate(const CoreMaterial& coreMaterial, const std::string& cacheKey) {
    if (complexPermeabilityEstimateRefusals.contains(cacheKey)) {
        return nullptr;
    }
    auto cached = complexPermeabilityEstimates.find(cacheKey);
    if (cached != complexPermeabilityEstimates.end()) {
        return cached->second;
    }
    std::shared_ptr<const ComplexPermeabilityEstimate> estimate;
    try {
        estimate = std::make_shared<const ComplexPermeabilityEstimate>(estimate_complex_permeability(coreMaterial));
    }
    catch (const ComplexPermeabilityEstimatorRefusal& exception) {
        // Only the estimator's own refusal (unsupported class, too few donors, ...) is recorded; the
        // caller turns it into ComplexPermeabilityEstimateUnavailableException. Anything else (e.g.
        // malformed MAS data) propagates.
        complexPermeabilityEstimateRefusals.emplace(cacheKey, exception.reason());
        return nullptr;
    }
    // One WARNING per material, as the heat-conductivity estimate (ABT #1753).
    OM_WARNING(estimate->summary());
    complexPermeabilityEstimates[cacheKey] = estimate;
    return estimate;
}

const ComplexPermeabilityEstimate& ComplexPermeability::estimate_for(const CoreMaterial& coreMaterial, const std::string& cacheKey, double frequency) {
    auto estimate = find_or_build_estimate(coreMaterial, cacheKey);
    if (!estimate) {
        auto [minimumFrequency, maximumFrequency] = frequency_range_for_key(cacheKey, coreMaterial.get_name());
        throw ComplexPermeabilityEstimateUnavailableException(coreMaterial.get_name(), frequency, minimumFrequency, maximumFrequency,
                                                              complexPermeabilityEstimateRefusals.at(cacheKey));
    }
    return *estimate;
}

ComplexPermeabilityValue ComplexPermeability::get_complex_permeability_with_source(const CoreMaterial& coreMaterial, double frequency) {
    // The range first: building the tables throws for parts that do not overlap.
    const std::string cacheKey = ensure_interpolators(coreMaterial);
    auto [minimumFrequency, maximumFrequency] = frequency_range_for_key(cacheKey, coreMaterial.get_name());
    if (frequency >= minimumFrequency && frequency <= maximumFrequency) {
        auto [real, imaginary] = get_tabulated_complex_permeability(coreMaterial, frequency);
        return {real, imaginary, ComplexPermeabilitySource::TABULATED, ComplexPermeabilitySource::TABULATED, "", ""};
    }
    // A frequency that is not a finite positive number has no complex permeability, tabulated or not.
    if (!std::isfinite(frequency) || frequency <= 0) {
        throw ComplexPermeabilityFrequencyOutOfRangeException(coreMaterial.get_name(), frequency, minimumFrequency, maximumFrequency);
    }
    // Outside the range where both parts are tabulated: each part from its own table where it has
    // one, the approved estimate where it has not (ABT #1861).
    const std::string& materialName = coreMaterial.get_name();
    const auto& estimate = estimate_for(coreMaterial, cacheKey, frequency);
    ComplexPermeabilityValue value;
    auto evaluatePart = [&](ComplexPermeabilityPart part, std::map<std::string, tk::spline>& interps,
                            std::map<std::string, std::pair<double, double>>& spans,
                            double& result, ComplexPermeabilitySource& source, std::string& method) {
        const auto& span = spans.at(cacheKey);
        if (frequency >= span.first && frequency <= span.second) {
            result = interps.at(cacheKey)(frequency);
            if (std::isnan(result)) {
                throw NaNResultException("complex permeability " + std::string(part == ComplexPermeabilityPart::REAL ? "real" : "imaginary") +
                                         " part must be a number, not NaN");
            }
            source = ComplexPermeabilitySource::TABULATED;
            return;
        }
        auto side = frequency < span.first ? ComplexPermeabilityGapSide::BELOW : ComplexPermeabilityGapSide::ABOVE;
        const auto& gap = estimate.gap(part, side);
        try {
            std::tie(result, source) = gap.evaluate(frequency);
        }
        catch (const ComplexPermeabilityEstimatorRefusal& exception) {
            throw ComplexPermeabilityEstimateUnavailableException(materialName, frequency, minimumFrequency, maximumFrequency, exception.reason());
        }
        method = gap.method(source);
    };
    evaluatePart(ComplexPermeabilityPart::REAL, complexPermeabilityRealInterps, complexPermeabilityRealFrequencySpans,
                 value.real, value.realSource, value.realMethod);
    evaluatePart(ComplexPermeabilityPart::IMAGINARY, complexPermeabilityImaginaryInterps, complexPermeabilityImaginaryFrequencySpans,
                 value.imaginary, value.imaginarySource, value.imaginaryMethod);
    return value;
}

std::pair<double, double> ComplexPermeability::get_complex_permeability(CoreMaterial coreMaterial, double frequency) {
    auto value = get_complex_permeability_with_source(coreMaterial, frequency);
    return {value.real, value.imaginary};
}

std::vector<ComplexPermeabilityEstimatedRange> ComplexPermeability::get_estimated_ranges(const CoreMaterial& coreMaterial) {
    auto estimate = find_or_build_estimate(coreMaterial, ensure_interpolators(coreMaterial));
    if (!estimate) {
        // Refused: nothing outside the table is answered, so nothing is estimated.
        return {};
    }
    return OpenMagnetics::get_estimated_ranges(*estimate);
}

} // namespace OpenMagnetics
