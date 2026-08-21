#include "physical_models/ComplexPermeability.h"
#include "physical_models/InitialPermeability.h"

#include "support/Utils.h"
#include <algorithm>
#include <math.h>
#include "support/Exceptions.h"


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
    // CLAMPS the interpolation to this table's frequency range, so whatever the last
    // tabulated point holds is what every higher frequency gets — a constant permeability,
    // not a falling one. At 0.01..100x the anchor, a material anchored at 18.8 kHz
    // (Nanoperm 80000) tabulated only 188 Hz..1.88 MHz, and every sweep above 1.88 MHz
    // read back the same mu' and mu''. That reported 23,565 at 25.6 MHz where the
    // material's own initial-permeability table says 156 — 151x high — and pushed the
    // modelled self-resonance of a common-mode choke from ~50 MHz down to 4.33 MHz
    // (ABT #843). The closed form is defined for any normalized frequency; only the
    // tabulation was short. Keep the same ~10 points per decade.
    auto normalizedFrequencies = logarithmic_spaced_array(0.01, 1e5, 70);
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
        double measuredPermeability = InitialPermeability::get_initial_permeability(coreMaterial, std::nullopt, std::nullopt, permeabilityPointFrequency);

        // Past its last tabulated point the initial-permeability interpolator turns back UP:
        // Nanoperm 80000 reads 156 at 25.6 MHz (its last point), then 451 at 50 MHz and
        // 1,057 at 100 MHz. Permeability does not recover after roll-off, so a rise means we
        // have run off the end of the data, not that the material got more permeable. Carry
        // on with the last real log-log slope instead of following the rise (ABT #843).
        if (lastTabulatedFrequency > 0 && permeabilityPointFrequency > lastTabulatedFrequency) {
            measuredPermeability = lastTabulatedPermeability * pow(permeabilityPointFrequency / lastTabulatedFrequency, beyondTheDataLogSlope);
        }

        if (measuredPermeability > 0 && muCoreReal > 0) {
            double scaleToMeasured = measuredPermeability / muCoreReal;
            muCoreReal *= scaleToMeasured;
            muCoreImag *= scaleToMeasured;
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


std::pair<double, double> ComplexPermeability::get_complex_permeability(CoreMaterial coreMaterial, double frequency) {
    // Fast path: if both interpolators are already cached for this material,
    // skip the expensive recomputation of the complex permeability data
    // (calculate_frequency_for_initial_permeability_drop runs O(40) pow() per
    // call and is invoked 3 times per material). Without this guard, advisers
    // that scan thousands of cores share only a few materials and end up
    // rebuilding the same per-material curves thousands of times, which
    // dominates DMC/CMC core selection wall time.
    const std::string& materialNameForCache = coreMaterial.get_name();
    if (complexPermeabilityRealInterps.contains(materialNameForCache) &&
        complexPermeabilityImaginaryInterps.contains(materialNameForCache)) {
        auto realSpan = complexPermeabilityRealFrequencySpans[materialNameForCache];
        double cachedReal = std::max(1., complexPermeabilityRealInterps[materialNameForCache](std::clamp(frequency, realSpan.first, realSpan.second)));
        if (std::isnan(cachedReal)) {
            throw NaNResultException("complex Permeability real part must be a number, not NaN");
        }
        auto imaginarySpan = complexPermeabilityImaginaryFrequencySpans[materialNameForCache];
        double cachedImag = complexPermeabilityImaginaryInterps[materialNameForCache](std::clamp(frequency, imaginarySpan.first, imaginarySpan.second));
        if (std::isnan(cachedImag)) {
            throw NaNResultException("complex Permeability imaginary part must be a number, not NaN");
        }
        return {cachedReal, cachedImag};
    }

    ComplexPermeabilityData complexPermeabilityData;
    if (!coreMaterial.get_permeability().get_complex()) {
        if (InitialPermeability::has_frequency_dependency(coreMaterial)) {
            complexPermeabilityData = calculate_complex_permeability_from_frequency_dependent_initial_permeability(coreMaterial);
        }
        else {
            throw MaterialDataMissingException(coreMaterial.get_name(), "Complex permeability");
        }
    }
    else {
        complexPermeabilityData = coreMaterial.get_permeability().get_complex().value();
    }

    auto realPart = complexPermeabilityData.get_real();
    auto imaginaryPart = complexPermeabilityData.get_imaginary();

    if (!std::holds_alternative<std::vector<PermeabilityPoint>>(realPart) ||
        !std::holds_alternative<std::vector<PermeabilityPoint>>(imaginaryPart)) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Complex permeability data is not in expected format for " + coreMaterial.get_name());
    }
    auto realPermeabilityPoints = std::get<std::vector<PermeabilityPoint>>(realPart);
    auto imaginaryPermeabilityPoints = std::get<std::vector<PermeabilityPoint>>(imaginaryPart);

    if (realPermeabilityPoints.size() < 2) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Not enough complex permeability data for  " + coreMaterial.get_name());
    }

    if (!complexPermeabilityRealInterps.contains(coreMaterial.get_name()))
    {
        int n = realPermeabilityPoints.size();
        std::vector<double> x, y;


        std::sort(realPermeabilityPoints.begin(), realPermeabilityPoints.end(), [](const PermeabilityPoint& b1, const PermeabilityPoint& b2) {
            return b1.get_frequency().value() < b2.get_frequency().value();
        });


        for (int i = 0; i < n; i++) {
            if (x.empty() || fabs(*realPermeabilityPoints[i].get_frequency() - x.back()) > 1e-9) {
                x.push_back(*realPermeabilityPoints[i].get_frequency());
                y.push_back(realPermeabilityPoints[i].get_value());
            }
        }

        tk::spline interp(x, y, tk::spline::cspline_hermite);
        complexPermeabilityRealInterps[coreMaterial.get_name()] = interp;
        complexPermeabilityRealFrequencySpans[coreMaterial.get_name()] = {x.front(), x.back()};
    }
    auto realSpan = complexPermeabilityRealFrequencySpans[coreMaterial.get_name()];
    double complexPermeabilityRealValue = std::max(1., complexPermeabilityRealInterps[coreMaterial.get_name()](std::clamp(frequency, realSpan.first, realSpan.second)));

    if (std::isnan(complexPermeabilityRealValue)) {
        throw NaNResultException("complex Permeability real part must be a number, not NaN");
    }

    if (!complexPermeabilityImaginaryInterps.contains(coreMaterial.get_name()))
    {
        int n = imaginaryPermeabilityPoints.size();
        std::vector<double> x, y;

        std::sort(imaginaryPermeabilityPoints.begin(), imaginaryPermeabilityPoints.end(), [](const PermeabilityPoint& b1, const PermeabilityPoint& b2) {
            return b1.get_frequency().value() < b2.get_frequency().value();
        });


        for (int i = 0; i < n; i++) {
            if (x.empty() || fabs(*imaginaryPermeabilityPoints[i].get_frequency() - x.back()) > 1e-9) {
                x.push_back(*imaginaryPermeabilityPoints[i].get_frequency());
                y.push_back(imaginaryPermeabilityPoints[i].get_value());
            }
        }

        tk::spline interp(x, y, tk::spline::cspline_hermite);
        complexPermeabilityImaginaryInterps[coreMaterial.get_name()] = interp;
        complexPermeabilityImaginaryFrequencySpans[coreMaterial.get_name()] = {x.front(), x.back()};
    }
    auto imaginarySpan = complexPermeabilityImaginaryFrequencySpans[coreMaterial.get_name()];
    double complexPermeabilityImaginaryValue = complexPermeabilityImaginaryInterps[coreMaterial.get_name()](std::clamp(frequency, imaginarySpan.first, imaginarySpan.second));

    if (std::isnan(complexPermeabilityImaginaryValue)) {
        throw NaNResultException("complex Permeability imaginary part must be a number, not NaN");
    }

    return {complexPermeabilityRealValue, complexPermeabilityImaginaryValue};
}

} // namespace OpenMagnetics
