#include "ComplexPermeability.h"

#include "Utils.h"

std::map<std::string, tk::spline> complexPermeabilityRealInterps;
std::map<std::string, tk::spline> complexPermeabilityImaginaryInterps;

namespace OpenMagnetics {

std::pair<double, double> ComplexPermeability::get_complex_permeability(std::string coreMaterialName, double frequency) {
    CoreMaterial coreMaterial = CoreWrapper::resolve_material(coreMaterialName);
    return get_complex_permeability(coreMaterial, frequency);
}

std::pair<double, double> ComplexPermeability::get_complex_permeability(CoreMaterial coreMaterial, double frequency) {
    if (!coreMaterial.get_permeability().get_complex()) {
        throw missing_material_data_exception("Missing complex data in material " + coreMaterial.get_name());
    }
    auto complexPermeabilityData = coreMaterial.get_permeability().get_complex().value();

    double complexPermeabilityValue = 1;

    auto realPart = complexPermeabilityData.get_real();
    auto imaginaryPart = complexPermeabilityData.get_imaginary();

    auto realPermeabilityPoints = std::get<std::vector<PermeabilityPoint>>(realPart);
    auto imaginaryPermeabilityPoints = std::get<std::vector<PermeabilityPoint>>(imaginaryPart);

    if (realPermeabilityPoints.size() < 2) {
        throw std::runtime_error("Not enough complex permeability data for  " + coreMaterial.get_name());
    }

    if (!complexPermeabilityRealInterps.contains(coreMaterial.get_name()))
    {
        int n = realPermeabilityPoints.size();
        std::vector<double> x, y;


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
        throw std::runtime_error("complex Permeability real part must be a number, not NaN");
    }

    if (!complexPermeabilityImaginaryInterps.contains(coreMaterial.get_name()))
    {
        int n = imaginaryPermeabilityPoints.size();
        std::vector<double> x, y;


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
        throw std::runtime_error("complex Permeability imaginary part must be a number, not NaN");
    }

    return {complexPermeabilityRealValue, complexPermeabilityImaginaryValue};
}

} // namespace OpenMagnetics
