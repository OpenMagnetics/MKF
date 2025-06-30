#include "physical_models/AmplitudePermeability.h"


namespace OpenMagnetics {

std::optional<double> AmplitudePermeability::get_amplitude_permeability(std::string coreMaterialName, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak, double temperature) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_amplitude_permeability(coreMaterial, magneticFluxDensityPeak, magneticFieldStrengthPeak, temperature);
}

std::optional<double> AmplitudePermeability::get_amplitude_permeability(CoreMaterial coreMaterial, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak, double temperature) {
    // TODO: Add more ways of getting this, like data from the manufacturer
    double amplitudePermeability = 1;
    if (!(magneticFieldStrengthPeak && !magneticFluxDensityPeak)) {
        throw std::runtime_error("Either H or B must be specified");
    }
    if (magneticFieldStrengthPeak) {
        auto [upperPath, lowerPath] = BHLoopRoshenModel().get_hysteresis_loop(coreMaterial, temperature, std::nullopt, magneticFieldStrengthPeak.value());
        if (upperPath.get_x_points().size() > 1){
            amplitudePermeability = fabs(upperPath.get_y_points()[1] - upperPath.get_y_points()[0]) / fabs(upperPath.get_x_points()[1] - upperPath.get_x_points()[0]) / constants.vacuumPermeability;
        }
        else {
            return std::nullopt;
        }
    }
    else {
        auto [upperPath, lowerPath] = BHLoopRoshenModel().get_hysteresis_loop(coreMaterial, temperature, magneticFluxDensityPeak.value(), std::nullopt);
        if (upperPath.get_x_points().size() > 1){
            amplitudePermeability = fabs(upperPath.get_y_points()[1] - upperPath.get_y_points()[0]) / fabs(upperPath.get_x_points()[1] - upperPath.get_x_points()[0]) / constants.vacuumPermeability;
        }
        else {
            return std::nullopt;
        }
    }
    return amplitudePermeability;
}

std::map<std::string, double> BHLoopRoshenModel::get_major_loop_parameters(double saturationMagneticFieldStrength,
                                                                           double saturationMagneticFluxDensity,
                                                                           double coerciveForce,
                                                                           double remanence) {
    std::map<std::string, double> majorLoopParameters;
    double a1 = 0;
    double b1 = 0;
    double b2 = 0;
    double Hc = coerciveForce;
    double H0 = saturationMagneticFieldStrength;
    double B0 = saturationMagneticFluxDensity;
    double H1 = 0;
    double B1 = remanence;
    double H2 = -saturationMagneticFieldStrength;
    double B2 = -saturationMagneticFluxDensity;
    b1 = (H0 / B0 + Hc / B0 - H1 / B1 - Hc / B1) / (H0 - H1);
    a1 = (Hc - B1 * b1 * Hc) / B1;
    b2 = (H2 + Hc - B2 * a1) / (B2 * fabs(H2 + Hc));
    majorLoopParameters["a1"] = a1;
    majorLoopParameters["b1"] = b1;
    majorLoopParameters["b2"] = b2;
    return majorLoopParameters;
}

std::pair<Curve2D, Curve2D> BHLoopRoshenModel::get_hysteresis_loop(std::string coreMaterialName, double temperature, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_hysteresis_loop(coreMaterial, temperature, magneticFluxDensityPeak, magneticFieldStrengthPeak);
}
std::pair<Curve2D, Curve2D> BHLoopRoshenModel::get_hysteresis_loop(CoreMaterial coreMaterial, double temperature, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak) {

    double saturationMagneticFieldStrength = Core::get_magnetic_field_strength_saturation(coreMaterial, temperature);
    double saturationMagneticFluxDensity = Core::get_magnetic_flux_density_saturation(coreMaterial, temperature, false);
    double coerciveForce = Core::get_coercive_force(coreMaterial, temperature);
    double remanence = Core::get_remanence(coreMaterial, temperature);

    auto majorLoopParameters = get_major_loop_parameters(saturationMagneticFieldStrength, saturationMagneticFluxDensity,
                                                         coerciveForce, remanence);

    double a1 = majorLoopParameters["a1"];
    double b1 = majorLoopParameters["b1"];
    double b2 = majorLoopParameters["b2"];

    std::vector<double> magneticFieldStrengthPoints;
    for (double i = -saturationMagneticFieldStrength; i <= saturationMagneticFieldStrength; i += constants.roshenMagneticFieldStrengthStep) {
        magneticFieldStrengthPoints.push_back(i);
    }

    auto bh_curve_half_loop = [&](double magneticFieldStrength, double a, double b) {
        return (magneticFieldStrength + coerciveForce) / (a + b * fabs(magneticFieldStrength + coerciveForce));
    };

    auto calculate_magnetic_flux_density = [&](double magneticFieldStrength, bool loop_is_upper = true) {
        double magneticFluxDensity;
        if (loop_is_upper) {
            if (-saturationMagneticFieldStrength <= magneticFieldStrength && magneticFieldStrength < -coerciveForce) {
                magneticFluxDensity = bh_curve_half_loop(magneticFieldStrength, a1, b2);
            }
            else {
                magneticFluxDensity = bh_curve_half_loop(magneticFieldStrength, a1, b1);
            }
        }
        else {
            if (-saturationMagneticFieldStrength <= magneticFieldStrength && magneticFieldStrength < coerciveForce) {
                magneticFluxDensity = -bh_curve_half_loop(-magneticFieldStrength, a1, b1);
            }
            else {
                magneticFluxDensity = -bh_curve_half_loop(-magneticFieldStrength, a1, b2);
            }
        }

        return magneticFluxDensity;
    };

    auto calculate_magnetic_flux_density_waveform = [&](std::vector<double> magneticFieldStrengthWaveform,
                                                  bool loop_is_upper = true) {
        std::vector<double> magneticFluxDensityWaveform;
        for (auto& magneticFieldStrength : magneticFieldStrengthWaveform) {
            double magneticFluxDensity = calculate_magnetic_flux_density(magneticFieldStrength, loop_is_upper);
            magneticFluxDensityWaveform.push_back(magneticFluxDensity);
        }

        return magneticFluxDensityWaveform;
    };

    std::vector<double> upperMagneticFluxDensityWaveform =
        calculate_magnetic_flux_density_waveform(magneticFieldStrengthPoints, true);
    std::vector<double> lowerMagneticFluxDensityWaveform =
        calculate_magnetic_flux_density_waveform(magneticFieldStrengthPoints, false);
    std::vector<double> difference;

    for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
        difference.push_back(fabs(upperMagneticFluxDensityWaveform[i] - lowerMagneticFluxDensityWaveform[i]));
    }

    Curve2D upperPath;
    Curve2D lowerPath;

    if (magneticFieldStrengthPeak) {
        auto it = std::min_element(magneticFieldStrengthPoints.begin(), magneticFieldStrengthPoints.end(), [magneticFieldStrengthPeak] (double a, double b) {
            return std::abs(magneticFieldStrengthPeak.value() - a) < std::abs(magneticFieldStrengthPeak.value() - b);
        });
        double indexOfDesiredHPeak = std::distance(std::begin(magneticFieldStrengthPoints), it);
        double differenceAtDesiredHPeak = fabs(upperMagneticFluxDensityWaveform[indexOfDesiredHPeak] - lowerMagneticFluxDensityWaveform[indexOfDesiredHPeak]);

        for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
            upperMagneticFluxDensityWaveform[i] -= differenceAtDesiredHPeak / 2;
        }
        for (size_t i = 0; i < lowerMagneticFluxDensityWaveform.size(); i++) {
            lowerMagneticFluxDensityWaveform[i] += differenceAtDesiredHPeak / 2;
        }

        for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
            auto elemB = upperMagneticFluxDensityWaveform[i];
            auto elemH = magneticFieldStrengthPoints[i];
            if (elemH <= magneticFieldStrengthPeak.value() && elemH >= -magneticFieldStrengthPeak.value()) {
                upperPath.get_mutable_x_points().push_back(elemH);
                upperPath.get_mutable_y_points().push_back(elemB);
            }
        }

        for (size_t i = 0; i < lowerMagneticFluxDensityWaveform.size(); i++) {
            auto elemB = lowerMagneticFluxDensityWaveform[i];
            auto elemH = magneticFieldStrengthPoints[i];
            if (elemH <= magneticFieldStrengthPeak.value() && elemH >= -magneticFieldStrengthPeak.value()) {
                lowerPath.get_mutable_x_points().push_back(elemH);
                lowerPath.get_mutable_y_points().push_back(elemB);
            }
        }
    }
    else if (magneticFluxDensityPeak) {
        double magneticFluxDensityDifference = magneticFluxDensityPeak.value();
        size_t timeout = 0;
        double abs_tol = 0.001;
        while (fabs(magneticFluxDensityDifference) > abs_tol && timeout < 10) {
            size_t mininumPosition =
                std::distance(std::begin(difference), std::min_element(std::begin(difference), std::end(difference)));
            magneticFluxDensityDifference =
                fabs(upperMagneticFluxDensityWaveform[mininumPosition]) - magneticFluxDensityPeak.value();

            for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
                upperMagneticFluxDensityWaveform[i] -= magneticFluxDensityDifference / 16;
            }
            for (size_t i = 0; i < lowerMagneticFluxDensityWaveform.size(); i++) {
                lowerMagneticFluxDensityWaveform[i] += magneticFluxDensityDifference / 16;
            }

            difference.clear();
            for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
                difference.push_back(fabs(upperMagneticFluxDensityWaveform[i] - lowerMagneticFluxDensityWaveform[i]));
            }
            timeout++;
            abs_tol += timeout * 0.0001;
        }

        for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
            auto elemB = upperMagneticFluxDensityWaveform[i];
            auto elemH = magneticFieldStrengthPoints[i];
            if (elemB <= magneticFluxDensityPeak.value() && elemB >= -magneticFluxDensityPeak.value()) {
                upperPath.get_mutable_x_points().push_back(elemH);
                upperPath.get_mutable_y_points().push_back(elemB);
            }
        }

        for (size_t i = 0; i < lowerMagneticFluxDensityWaveform.size(); i++) {
            auto elemB = lowerMagneticFluxDensityWaveform[i];
            auto elemH = magneticFieldStrengthPoints[i];
            if (elemB <= magneticFluxDensityPeak.value() && elemB >= -magneticFluxDensityPeak.value()) {
                lowerPath.get_mutable_x_points().push_back(elemH);
                lowerPath.get_mutable_y_points().push_back(elemB);
            }
        }
    }

    return {upperPath, lowerPath};
}
} // namespace OpenMagnetics
