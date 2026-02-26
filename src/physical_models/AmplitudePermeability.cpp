#include "physical_models/AmplitudePermeability.h"
#include "support/Exceptions.h"


namespace OpenMagnetics {

std::optional<double> AmplitudePermeability::get_amplitude_permeability(std::string coreMaterialName, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak, double temperature) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_amplitude_permeability(coreMaterial, magneticFluxDensityPeak, magneticFieldStrengthPeak, temperature);
}

std::optional<double> AmplitudePermeability::get_amplitude_permeability(CoreMaterial coreMaterial, std::optional<double> magneticFluxDensityPeak, std::optional<double> magneticFieldStrengthPeak, double temperature) {
    // Calculate amplitude permeability using B-H loop slope
    // This uses the average slope of the loop, not the instantaneous derivative
    
    double amplitudePermeability = 1;
    if (!(magneticFieldStrengthPeak && !magneticFluxDensityPeak)) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Either H or B must be specified");
    }
    if (magneticFieldStrengthPeak) {
        auto [upperPath, lowerPath] = BHLoopRoshenModel().get_hysteresis_loop(coreMaterial, temperature, std::nullopt, magneticFieldStrengthPeak.value());
        if (upperPath.get_x_points().size() > 1){
            // Calculate slope between first two points on upper branch
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

std::pair<Curve2D, Curve2D> BHLoopRoshenModel::get_hysteresis_loop(std::string coreMaterialName, 
                                                                    double temperature, 
                                                                    std::optional<double> magneticFluxDensityPeak, 
                                                                    std::optional<double> magneticFieldStrengthPeak) {
    CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);
    return get_hysteresis_loop(coreMaterial, temperature, magneticFluxDensityPeak, magneticFieldStrengthPeak);
}

std::pair<Curve2D, Curve2D> BHLoopRoshenModel::get_hysteresis_loop(CoreMaterial coreMaterial, 
                                                                    double temperature, 
                                                                    std::optional<double> magneticFluxDensityPeak, 
                                                                    std::optional<double> magneticFieldStrengthPeak) {

    double saturationMagneticFieldStrength = Core::get_magnetic_field_strength_saturation(coreMaterial, temperature);
    double saturationMagneticFluxDensity = Core::get_magnetic_flux_density_saturation(coreMaterial, temperature, false);
    double coerciveForce = Core::get_coercive_force(coreMaterial, temperature);
    double remanence = Core::get_remanence(coreMaterial, temperature);

    auto majorLoopParameters = get_major_loop_parameters(saturationMagneticFieldStrength, saturationMagneticFluxDensity,
                                                         coerciveForce, remanence);

    double a1 = majorLoopParameters["a1"];
    double b1 = majorLoopParameters["b1"];
    double b2 = majorLoopParameters["b2"];

    // OPTIMIZED: Use fewer points with adaptive spacing instead of fixed step
    std::vector<double> magneticFieldStrengthPoints;
    
    // Critical points that must be included
    magneticFieldStrengthPoints.push_back(-saturationMagneticFieldStrength);
    magneticFieldStrengthPoints.push_back(-coerciveForce);
    magneticFieldStrengthPoints.push_back(0);
    magneticFieldStrengthPoints.push_back(coerciveForce);
    magneticFieldStrengthPoints.push_back(saturationMagneticFieldStrength);
    
    // Add points near origin for accurate permeability calculation
    // These are critical for slope calculation at low H
    double H_near_zero = saturationMagneticFieldStrength * 0.001;  // Very close to origin
    magneticFieldStrengthPoints.push_back(H_near_zero);
    magneticFieldStrengthPoints.push_back(-H_near_zero);
    magneticFieldStrengthPoints.push_back(coerciveForce * 0.1);
    magneticFieldStrengthPoints.push_back(-coerciveForce * 0.1);
    
    // Add points around the operating region if specified
    if (magneticFieldStrengthPeak) {
        double H_peak = magneticFieldStrengthPeak.value();
        // Add points near the peak for better resolution
        if (H_peak > 0) {
            magneticFieldStrengthPoints.push_back(H_peak * 0.001);  // Very close to 0
            magneticFieldStrengthPoints.push_back(H_peak * 0.01);   // Close to 0
            magneticFieldStrengthPoints.push_back(H_peak * 0.1);
            magneticFieldStrengthPoints.push_back(H_peak * 0.5);
            magneticFieldStrengthPoints.push_back(H_peak * 0.9);
            magneticFieldStrengthPoints.push_back(H_peak * 0.95);
            magneticFieldStrengthPoints.push_back(H_peak);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.001);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.01);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.1);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.5);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.9);
            magneticFieldStrengthPoints.push_back(-H_peak * 0.95);
            magneticFieldStrengthPoints.push_back(-H_peak);
        }
    }
    
    // Sort and remove duplicates
    std::sort(magneticFieldStrengthPoints.begin(), magneticFieldStrengthPoints.end());
    auto last = std::unique(magneticFieldStrengthPoints.begin(), magneticFieldStrengthPoints.end(),
        [](double a, double b) { return fabs(a - b) < 1e-6; });
    magneticFieldStrengthPoints.erase(last, magneticFieldStrengthPoints.end());
    
    // Ensure bounds
    if (magneticFieldStrengthPoints.front() > -saturationMagneticFieldStrength) {
        magneticFieldStrengthPoints.insert(magneticFieldStrengthPoints.begin(), -saturationMagneticFieldStrength);
    }
    if (magneticFieldStrengthPoints.back() < saturationMagneticFieldStrength) {
        magneticFieldStrengthPoints.push_back(saturationMagneticFieldStrength);
    }

    auto bh_curve_half_loop = [&](double magneticFieldStrength, double a, double b) {
        return (magneticFieldStrength + coerciveForce) / (a + b * fabs(magneticFieldStrength + coerciveForce));
    };

    auto calculate_magnetic_flux_density = [&](double magneticFieldStrength, bool loop_is_upper = true) {
        double magneticFluxDensity;
        if (loop_is_upper) {
            if (magneticFieldStrength < -coerciveForce) {
                magneticFluxDensity = bh_curve_half_loop(magneticFieldStrength, a1, b2);
            }
            else {
                magneticFluxDensity = bh_curve_half_loop(magneticFieldStrength, a1, b1);
            }
        }
        else {
            if (magneticFieldStrength < coerciveForce) {
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
        magneticFluxDensityWaveform.reserve(magneticFieldStrengthWaveform.size());
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

    Curve2D upperPath;
    Curve2D lowerPath;

    if (magneticFieldStrengthPeak) {
        // OPTIMIZED: Find index with binary search instead of min_element
        double H_peak = magneticFieldStrengthPeak.value();
        size_t indexOfDesiredHPeak = 0;
        double minDiff = fabs(magneticFieldStrengthPoints[0] - H_peak);
        for (size_t i = 1; i < magneticFieldStrengthPoints.size(); ++i) {
            double diff = fabs(magneticFieldStrengthPoints[i] - H_peak);
            if (diff < minDiff) {
                minDiff = diff;
                indexOfDesiredHPeak = i;
            }
        }
        
        double differenceAtDesiredHPeak = fabs(upperMagneticFluxDensityWaveform[indexOfDesiredHPeak] - 
                                                lowerMagneticFluxDensityWaveform[indexOfDesiredHPeak]);

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
        
        // OPTIMIZED: Pre-allocate difference vector
        std::vector<double> difference;
        difference.reserve(upperMagneticFluxDensityWaveform.size());
        
        while (fabs(magneticFluxDensityDifference) > abs_tol && timeout < 10) {
            difference.clear();
            for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
                difference.push_back(fabs(upperMagneticFluxDensityWaveform[i] - lowerMagneticFluxDensityWaveform[i]));
            }
            
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
