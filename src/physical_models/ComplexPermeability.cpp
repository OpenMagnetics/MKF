#include "physical_models/ComplexPermeability.h"
#include "physical_models/InitialPermeability.h"

#include "support/Utils.h"
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
    double frequencyFor67Point78Drop = InitialPermeability::calculate_frequency_for_initial_permeability_drop(coreMaterial, 0.6778);
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

    auto normalizedFrequencies = logarithmic_spaced_array(0.01, 100, 40);
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
    }
    double complexPermeabilityRealValue = std::max(1., complexPermeabilityRealInterps[coreMaterial.get_name()](frequency));

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
    }
    double complexPermeabilityImaginaryValue = complexPermeabilityImaginaryInterps[coreMaterial.get_name()](frequency);

    if (std::isnan(complexPermeabilityImaginaryValue)) {
        throw NaNResultException("complex Permeability imaginary part must be a number, not NaN");
    }

    return {complexPermeabilityRealValue, complexPermeabilityImaginaryValue};
}

} // namespace OpenMagnetics
