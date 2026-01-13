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
ComplexPermeabilityData ComplexPermeability::calculate_complex_permeability_from_frequency_dependent_initial_permeability(CoreMaterial coreMaterial) {
    double frequencyFor67Point78Drop = InitialPermeability::calculate_frequency_for_initial_permeability_drop(coreMaterial, 0.6778);
    double initialPermeability = InitialPermeability::get_initial_permeability(coreMaterial);
    auto normalizedFrequencies = logarithmic_spaced_array(0.01, 100, 40);
    std::vector<PermeabilityPoint> real;
    std::vector<PermeabilityPoint> imaginary;

    for (auto normalizedFrequency : normalizedFrequencies) {
        double complexPermeabilityRealNormalizedValue = (sin(2 * sqrt(normalizedFrequency)) + sinh(2 * sqrt(normalizedFrequency))) / (2 * sqrt(normalizedFrequency) * (cos(2 * sqrt(normalizedFrequency)) + cosh(2 * sqrt(normalizedFrequency))));
        double complexPermeabilityImaginaryNormalizedValue =-(sin(2 * sqrt(normalizedFrequency)) - sinh(2 * sqrt(normalizedFrequency))) / (2 * sqrt(normalizedFrequency) * (cos(2 * sqrt(normalizedFrequency)) + cosh(2 * sqrt(normalizedFrequency))));
        double permeabilityPointFrequency = normalizedFrequency * frequencyFor67Point78Drop;
        double complexPermeabilityRealValue = initialPermeability * complexPermeabilityRealNormalizedValue;
        double complexPermeabilityImaginaryValue = initialPermeability * complexPermeabilityImaginaryNormalizedValue;
        PermeabilityPoint realPermeabilityPoint;
        realPermeabilityPoint.set_frequency(permeabilityPointFrequency);
        realPermeabilityPoint.set_value(complexPermeabilityRealValue);
        PermeabilityPoint imaginaryPermeabilityPoint;
        imaginaryPermeabilityPoint.set_frequency(permeabilityPointFrequency);
        imaginaryPermeabilityPoint.set_value(complexPermeabilityImaginaryValue);
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
            if (x.size() == 0 || (*realPermeabilityPoints[i].get_frequency()) != x.back()) {
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
            if (x.size() == 0 || (*imaginaryPermeabilityPoints[i].get_frequency()) != x.back()) {
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
