#include "CoreLosses.h"
#include "Resistivity.h"

#include "Constants.h"
#include "InputsWrapper.h"

#include <cmath>
#include <cfloat>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "spline.h"
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

std::shared_ptr<CoreLossesModel> CoreLossesModel::factory(std::map<std::string, std::string> models) {
    return factory(magic_enum::enum_cast<OpenMagnetics::CoreLossesModels>(models["coreLosses"]).value());
}

std::shared_ptr<CoreLossesModel> CoreLossesModel::factory(json models) {
    std::string model = models["coreLosses"];
    return factory(magic_enum::enum_cast<OpenMagnetics::CoreLossesModels>(model).value());
}

std::shared_ptr<CoreLossesModel> CoreLossesModel::factory(CoreLossesModels modelName) {
    if (modelName == CoreLossesModels::STEINMETZ) {
        std::shared_ptr<CoreLossesModel> coreLossesModel(new CoreLossesSteinmetzModel);
        return coreLossesModel;
    }
    else if (modelName == CoreLossesModels::IGSE) {
        std::shared_ptr<CoreLossesModel> coreLossesModel(new CoreLossesIGSEModel);
        return coreLossesModel;
    }
    else if (modelName == CoreLossesModels::MSE) {
        std::shared_ptr<CoreLossesModel> coreLossesModel(new CoreLossesMSEModel);
        return coreLossesModel;
    }
    else if (modelName == CoreLossesModels::NSE) {
        std::shared_ptr<CoreLossesModel> coreLossesModel(new CoreLossesNSEModel);
        return coreLossesModel;
    }
    else if (modelName == CoreLossesModels::ALBACH) {
        std::shared_ptr<CoreLossesModel> coreLossesModel(new CoreLossesAlbachModel);
        return coreLossesModel;
    }
    else if (modelName == CoreLossesModels::BARG) {
        std::shared_ptr<CoreLossesModel> coreLossesModel(new CoreLossesBargModel);
        return coreLossesModel;
    }
    else if (modelName == CoreLossesModels::ROSHEN) {
        std::shared_ptr<CoreLossesModel> coreLossesModel(new CoreLossesRoshenModel);
        return coreLossesModel;
    }
    else if (modelName == CoreLossesModels::PROPRIETARY) {
        std::shared_ptr<CoreLossesModel> coreLossesModel(new CoreLossesProprietaryModel);
        return coreLossesModel;
    }

    else
        throw std::runtime_error("Unknown Core losses mode, available options are: {STEINMETZ, IGSE, BARG, ALBACH, "
                                 "ROSHEN, OUYANG, NSE, MSE, PROPRIETARY}");
}

CoreLossesMethodData get_method_data(CoreMaterial materialData, std::string method) {
    auto volumetricLossesMethodsVariants = materialData.get_volumetric_losses();
    for (auto& volumetricLossesMethodVariant : volumetricLossesMethodsVariants) {
        auto volumetricLossesMethods = volumetricLossesMethodVariant.second;
        for (auto& volumetricLossesMethod : volumetricLossesMethods) {
            if (std::holds_alternative<OpenMagnetics::CoreLossesMethodData>(volumetricLossesMethod)) {
                auto methodData = std::get<OpenMagnetics::CoreLossesMethodData>(volumetricLossesMethod);
                if (methodData.get_method() == method) {
                    return methodData;
                }
            }
        }
    }
    throw std::runtime_error("Material " + materialData.get_name() + " does not have method:" + method);
}

OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum CoreLossesModel::get_steinmetz_coefficients(
    CoreMaterialDataOrNameUnion material,
    double frequency) {
    OpenMagnetics::CoreMaterial materialData;
    // If the material is a string, we have to load its data from the database, unless it is dummy (in order to avoid
    // long loading operatings)
    if (std::holds_alternative<std::string>(material) && std::get<std::string>(material) != "dummy") {
        materialData = OpenMagnetics::find_core_material_by_name(std::get<std::string>(material));
    }
    else {
        materialData = std::get<OpenMagnetics::CoreMaterial>(material);
    }

    auto volumetricLossesMethodsVariants = materialData.get_volumetric_losses();

    auto steinmetzData = get_method_data(materialData, "steinmetz");
    auto ranges = steinmetzData.get_ranges().value();
    double minimumMaterialFrequency = 100000000;
    double minimumMaterialFrequencyIndex = -1;
    double maximumMaterialFrequency = 0;
    double maximumMaterialFrequencyIndex = -1;
    for (size_t i = 0; i < ranges.size(); ++i) {
        if (!ranges[i].get_minimum_frequency()) {
            throw std::runtime_error("Missing minimum frequency in material");
        }
        if (!ranges[i].get_maximum_frequency()) {
            throw std::runtime_error("Missing maximum frequency in material");
        }

        if (frequency >= ranges[i].get_minimum_frequency().value() &&
            frequency <= ranges[i].get_maximum_frequency().value()) {
            return ranges[i];
        }

        if (minimumMaterialFrequency > ranges[i].get_minimum_frequency().value()) {
            minimumMaterialFrequency = ranges[i].get_minimum_frequency().value();
            minimumMaterialFrequencyIndex = i;
        }
        if (maximumMaterialFrequency < ranges[i].get_maximum_frequency().value()) {
            maximumMaterialFrequency = ranges[i].get_minimum_frequency().value();
            maximumMaterialFrequencyIndex = i;
        }
    }

    if (frequency < minimumMaterialFrequency) {
        return ranges[minimumMaterialFrequencyIndex];
    }
    if (frequency > maximumMaterialFrequency) {
        return ranges[maximumMaterialFrequencyIndex];
    }

    throw std::runtime_error("Error getting Steinmetz coefficients");
}

std::map<std::string, double> CoreLossesSteinmetzModel::get_core_losses(CoreWrapper core,
                                                                        OperatingPointExcitation excitation,
                                                                        double temperature) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();

    magneticFluxDensity = InputsWrapper::standarize_waveform(magneticFluxDensity, frequency);

    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() -
                                       magneticFluxDensity.get_processed().value().get_offset();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();

    OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(core.get_functional_description().get_material(), frequency);
    }
    double k = steinmetzDatum.get_k();
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double volumetricLosses = k * pow(frequency, alpha) * pow(magneticFluxDensityAcPeak, beta);

    volumetricLosses = CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);

    std::map<std::string, double> result;
    result["totalLosses"] = volumetricLosses * effectiveVolume;
    result["totalVolumetricLosses"] = volumetricLosses;

    return result;
};

double CoreLossesIGSEModel::get_ki(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum) {
    auto constants = Constants();
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double k = steinmetzDatum.get_k();
    std::vector<double> theta_vect(constants.number_points_samples_waveforms);
    std::generate(theta_vect.begin(), theta_vect.end(), [n = 0, &constants]() mutable {
        return n++ * 2 * std::numbers::pi / constants.number_points_samples_waveforms;
    });
    double theta_int = 0;

    for (auto& theta : theta_vect) {
        theta_int += pow(fabs(cos(theta)), alpha) * pow(2, beta - alpha) * 2 * std::numbers::pi /
                     constants.number_points_samples_waveforms;
    }

    double ki = k / (pow(2 * std::numbers::pi, alpha - 1) * theta_int);
    return ki;
}

std::map<std::string, double> CoreLossesIGSEModel::get_core_losses(CoreWrapper core,
                                                                   OperatingPointExcitation excitation,
                                                                   double temperature) {
    auto constants = Constants();
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    magneticFluxDensity = InputsWrapper::standarize_waveform(magneticFluxDensity, frequency);
    double magneticFluxDensityPeakToPeak = magneticFluxDensity.get_processed().value().get_peak_to_peak().value();
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(core.get_functional_description().get_material(), frequency);
    }

    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double ki = get_ki(steinmetzDatum);

    double volumetricLossesSum = 0;
    double timeDifference;

    for (size_t i = 0; i < magneticFluxDensityWaveform.size() - 1; ++i) {
        if (magneticFluxDensity.get_waveform().value().get_time()) {
            timeDifference = magneticFluxDensityTime[i + 1] - magneticFluxDensityTime[i];
        }
        else {
            timeDifference = 1 / frequency / constants.number_points_samples_waveforms;
        }
        volumetricLossesSum +=
            pow(fabs((magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / timeDifference), alpha) *
            timeDifference;
    }

    double volumetricLosses = ki * pow(magneticFluxDensityPeakToPeak, beta - alpha) * frequency * volumetricLossesSum;
    volumetricLosses = CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);

    std::map<std::string, double> result;
    result["totalLosses"] = volumetricLosses * effectiveVolume;
    result["totalVolumetricLosses"] = volumetricLosses;

    return result;
}

std::map<std::string, double> CoreLossesAlbachModel::get_core_losses(CoreWrapper core,
                                                                     OperatingPointExcitation excitation,
                                                                     double temperature) {
    auto constants = Constants();
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    magneticFluxDensity = InputsWrapper::standarize_waveform(magneticFluxDensity, frequency);
    double magneticFluxDensityPeakToPeak = magneticFluxDensity.get_processed().value().get_peak_to_peak().value();
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() -
                                       magneticFluxDensity.get_processed().value().get_offset();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(core.get_functional_description().get_material(), frequency);
    }

    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    double k = steinmetzDatum.get_k();
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();

    double equivalentSinusoidalFrequency = 0;
    double timeDifference;

    for (size_t i = 0; i < magneticFluxDensityWaveform.size() - 1; ++i) {
        if (magneticFluxDensity.get_waveform().value().get_time()) {
            timeDifference = magneticFluxDensityTime[i + 1] - magneticFluxDensityTime[i];
        }
        else {
            timeDifference = 1 / frequency / constants.number_points_samples_waveforms;
        }
        equivalentSinusoidalFrequency +=
            pow((magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / magneticFluxDensityPeakToPeak,
                2) /
            timeDifference;
    }

    equivalentSinusoidalFrequency = 2 / pow(std::numbers::pi, 2) * equivalentSinusoidalFrequency;

    double volumetricLosses =
        k * frequency * pow(equivalentSinusoidalFrequency, alpha - 1) * pow(magneticFluxDensityAcPeak, beta);
    volumetricLosses = CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);

    std::map<std::string, double> result;
    result["totalLosses"] = volumetricLosses * effectiveVolume;
    result["totalVolumetricLosses"] = volumetricLosses;

    return result;
}

std::map<std::string, double> CoreLossesMSEModel::get_core_losses(CoreWrapper core,
                                                                  OperatingPointExcitation excitation,
                                                                  double temperature) {
    auto constants = Constants();
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    magneticFluxDensity = InputsWrapper::standarize_waveform(magneticFluxDensity, frequency);
    double magneticFluxDensityPeakToPeak = magneticFluxDensity.get_processed().value().get_peak_to_peak().value();
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() -
                                       magneticFluxDensity.get_processed().value().get_offset();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(core.get_functional_description().get_material(), frequency);
    }

    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    double k = steinmetzDatum.get_k();
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();

    double equivalentSinusoidalFrequency = 0;
    double timeDifference;

    for (size_t i = 0; i < magneticFluxDensityWaveform.size() - 1; ++i) {
        if (magneticFluxDensity.get_waveform().value().get_time()) {
            timeDifference = magneticFluxDensityTime[i + 1] - magneticFluxDensityTime[i];
        }
        else {
            timeDifference = 1 / frequency / constants.number_points_samples_waveforms;
        }
        equivalentSinusoidalFrequency +=
            pow((magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / timeDifference, 2) *
            timeDifference;
    }

    equivalentSinusoidalFrequency =
        2 / pow(std::numbers::pi, 2) / pow(magneticFluxDensityPeakToPeak, 2) * equivalentSinusoidalFrequency;

    double volumetricLosses =
        k * frequency * pow(equivalentSinusoidalFrequency, alpha - 1) * pow(magneticFluxDensityAcPeak, beta);
    volumetricLosses = CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);

    std::map<std::string, double> result;
    result["totalLosses"] = volumetricLosses * effectiveVolume;
    result["totalVolumetricLosses"] = volumetricLosses;

    return result;
}

double CoreLossesNSEModel::get_kn(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum) {
    auto constants = Constants();
    double alpha = steinmetzDatum.get_alpha();
    double k = steinmetzDatum.get_k();
    std::vector<double> theta_vect(constants.number_points_samples_waveforms);
    std::generate(theta_vect.begin(), theta_vect.end(), [n = 0, &constants]() mutable {
        return n++ * 2 * std::numbers::pi / constants.number_points_samples_waveforms;
    });
    double theta_int = 0;

    for (auto& theta : theta_vect) {
        theta_int += pow(fabs(cos(theta)), alpha) * 2 * std::numbers::pi / constants.number_points_samples_waveforms;
    }

    double ki = k / (pow(2 * std::numbers::pi, alpha - 1) * theta_int);
    return ki;
}

std::map<std::string, double> CoreLossesNSEModel::get_core_losses(CoreWrapper core,
                                                                  OperatingPointExcitation excitation,
                                                                  double temperature) {
    auto constants = Constants();
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    magneticFluxDensity = InputsWrapper::standarize_waveform(magneticFluxDensity, frequency);
    double magneticFluxDensityPeakToPeak = magneticFluxDensity.get_processed().value().get_peak_to_peak().value();
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(core.get_functional_description().get_material(), frequency);
    }

    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double kn = get_kn(steinmetzDatum);

    double volumetricLossesSum = 0;
    double timeDifference;

    for (size_t i = 0; i < magneticFluxDensityWaveform.size() - 1; ++i) {
        if (magneticFluxDensity.get_waveform().value().get_time()) {
            timeDifference = magneticFluxDensityTime[i + 1] - magneticFluxDensityTime[i];
        }
        else {
            timeDifference = 1 / frequency / constants.number_points_samples_waveforms;
        }
        volumetricLossesSum +=
            pow(abs((magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / timeDifference), alpha) *
            timeDifference;
    }

    double volumetricLosses =
        kn * pow(magneticFluxDensityPeakToPeak / 2, beta - alpha) * frequency * volumetricLossesSum;
    volumetricLosses = CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);

    std::map<std::string, double> result;
    result["totalLosses"] = volumetricLosses * effectiveVolume;
    result["totalVolumetricLosses"] = volumetricLosses;

    return result;
}

std::map<std::string, double> CoreLossesBargModel::get_core_losses(CoreWrapper core,
                                                                   OperatingPointExcitation excitation,
                                                                   double temperature) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    magneticFluxDensity = InputsWrapper::standarize_waveform(magneticFluxDensity, frequency);
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(core.get_functional_description().get_material(), frequency);
    }

    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double k = steinmetzDatum.get_k();
    double dutyCycle = 0.5;
    // double dutyCycle = magneticFluxDensity.get_processed().value().get_duty_cycle().value();
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() -
                                       magneticFluxDensity.get_processed().value().get_offset();

    double lossesFrameT1 = std::numbers::pi / 4 * k * pow(frequency, alpha) * pow(magneticFluxDensityAcPeak, beta);
    lossesFrameT1 = CoreLossesModel::apply_temperature_coefficients(lossesFrameT1, steinmetzDatum, temperature);
    std::vector<double> dutyCycleValues = {0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,
                                           0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9};
    std::vector<double> factorValues = {1.45,  1.4,  1.35, 1.275, 1.25,  1.2,  1.15, 1.075, 1,
                                        1.075, 1.15, 1.2,  1.25,  1.275, 1.35, 1.4,  1.45};

    size_t n = dutyCycleValues.size();

    tk::spline interp(dutyCycleValues, factorValues, tk::spline::cspline, true);
    double dutyCycleFactor = std::max(1., interp(dutyCycle));

    double volumetricLosses = dutyCycleFactor * lossesFrameT1;

    std::map<std::string, double> result;
    result["totalLosses"] = volumetricLosses * effectiveVolume;
    result["totalVolumetricLosses"] = volumetricLosses;

    return result;
}

std::map<std::string, double> CoreLossesRoshenModel::get_core_losses(CoreWrapper core,
                                                                     OperatingPointExcitation excitation,
                                                                     double temperature) {
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    auto parameters = get_roshen_parameters(core, excitation, temperature);
    double hysteresisVolumetricLosses = get_hysteresis_losses_density(parameters, excitation);
    double eddyCurrentsVolumetricLosses = get_eddy_current_losses_density(core, excitation, parameters["resistivity"]);
    double excessEddyCurrentsVolumetricLosses = 0;
    if (parameters.count("excessLossesCoefficient")) {
        excessEddyCurrentsVolumetricLosses = get_excess_eddy_current_losses_density(
            excitation, parameters["resistivity"], parameters["excessLossesCoefficient"]);
    }
    double volumetricLosses =
        hysteresisVolumetricLosses + eddyCurrentsVolumetricLosses + excessEddyCurrentsVolumetricLosses;

    std::map<std::string, double> result;
    result["totalLosses"] = volumetricLosses * effectiveVolume;
    result["hysteresisLosses"] = hysteresisVolumetricLosses * effectiveVolume;
    result["eddyCurrentLosses"] = (eddyCurrentsVolumetricLosses + excessEddyCurrentsVolumetricLosses) * effectiveVolume;
    result["totalVolumetricLosses"] = volumetricLosses;

    return result;
}

std::map<std::string, double> CoreLossesRoshenModel::get_roshen_parameters(CoreWrapper core,
                                                                           OperatingPointExcitation excitation,
                                                                           double temperature) {
    std::map<std::string, double> roshenParameters;
    auto materialData =  core.get_material();

    auto roshenData = get_method_data(materialData, "roshen");

    roshenParameters["coerciveForce"] = core.get_coercive_force(temperature);
    roshenParameters["remanence"] = core.get_remanence(temperature);
    roshenParameters["saturationMagneticFluxDensity"] = core.get_magnetic_flux_density_saturation(temperature, false);
    roshenParameters["saturationMagneticFieldStrength"] = core.get_magnetic_field_strength_saturation(temperature);

    if (roshenData.get_coefficients()) {
        auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
        double frequency = excitation.get_frequency();
        magneticFluxDensity = InputsWrapper::standarize_waveform(magneticFluxDensity, frequency);
        double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() -
                                           magneticFluxDensity.get_processed().value().get_offset();

        auto roshenCoefficients = roshenData.get_coefficients().value();
        roshenParameters["excessLossesCoefficient"] = roshenCoefficients.get_excess_losses_coefficient();
        roshenParameters["resistivityFrequencyCoefficient"] =
            roshenCoefficients.get_resistivity_frequency_coefficient();
        roshenParameters["resistivityMagneticFluxDensityCoefficient"] =
            roshenCoefficients.get_resistivity_magnetic_flux_density_coefficient();
        roshenParameters["resistivityOffset"] = roshenCoefficients.get_resistivity_offset();
        roshenParameters["resistivityTemperatureCoefficient"] =
            roshenCoefficients.get_resistivity_temperature_coefficient();
        double resistivity = roshenParameters["resistivityOffset"] +
                             roshenParameters["resistivityTemperatureCoefficient"] * (temperature - 25) +
                             roshenParameters["resistivityMagneticFluxDensityCoefficient"] * magneticFluxDensityAcPeak +
                             roshenParameters["resistivityFrequencyCoefficient"] * frequency;
        roshenParameters["resistivity"] = resistivity;
    }
    else {
        auto resistivityModel = OpenMagnetics::ResistivityModel::factory(OpenMagnetics::ResistivityModels::CORE_MATERIAL);
        auto resistivity = (*resistivityModel).get_resistivity(materialData, temperature);
        roshenParameters["resistivity"] = resistivity;
    }

    return roshenParameters;
}

std::map<std::string, double> get_major_loop_parameters(double saturationMagneticFieldStrength,
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

double CoreLossesRoshenModel::get_hysteresis_losses_density(std::map<std::string, double> parameters,
                                                            OperatingPointExcitation excitation) {
    auto constants = Constants();
    double saturationMagneticFieldStrength = parameters["saturationMagneticFieldStrength"];
    double saturationMagneticFluxDensity = parameters["saturationMagneticFluxDensity"];
    double coerciveForce = parameters["coerciveForce"];
    double remanence = parameters["remanence"];
    double frequency = excitation.get_frequency();
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() -
                                       magneticFluxDensity.get_processed().value().get_offset();
    auto majorLoopParameters = get_major_loop_parameters(saturationMagneticFieldStrength, saturationMagneticFluxDensity,
                                                         coerciveForce, remanence);

    double a1 = majorLoopParameters["a1"];
    double b1 = majorLoopParameters["b1"];
    double b2 = majorLoopParameters["b2"];

    std::vector<double> magneticFieldStrengthPoints;
    for (double i = -saturationMagneticFieldStrength; i <= saturationMagneticFieldStrength;
         i += constants.roshen_magnetic_field_strength_step) {
        magneticFieldStrengthPoints.push_back(i);
    }

    auto bh_curve_half_loop = [&](double magneticFieldStrength, double a, double b) {
        return (magneticFieldStrength + coerciveForce) / (a + b * fabs(magneticFieldStrength + coerciveForce));
    };

    auto get_magnetic_flux_density = [&](double magneticFieldStrength, bool loop_is_upper = true) {
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

    auto get_magnetic_flux_density_waveform = [&](std::vector<double> magneticFieldStrength_waveform,
                                                  bool loop_is_upper = true) {
        std::vector<double> magneticFluxDensityWaveform;
        for (auto& magneticFieldStrength : magneticFieldStrength_waveform) {
            double magneticFluxDensity = get_magnetic_flux_density(magneticFieldStrength, loop_is_upper);
            magneticFluxDensityWaveform.push_back(magneticFluxDensity);
        }

        return magneticFluxDensityWaveform;
    };

    std::vector<double> upperMagneticFluxDensityWaveform =
        get_magnetic_flux_density_waveform(magneticFieldStrengthPoints, true);
    std::vector<double> lowerMagneticFluxDensityWaveform =
        get_magnetic_flux_density_waveform(magneticFieldStrengthPoints, false);
    std::vector<double> difference;

    _hysteresisMajorH = magneticFieldStrengthPoints;
    _hysteresisMajorLoopTop = upperMagneticFluxDensityWaveform;
    _hysteresisMajorLoopBottom = lowerMagneticFluxDensityWaveform;
    for (size_t i = 0; i < upperMagneticFluxDensityWaveform.size(); i++) {
        difference.push_back(fabs(upperMagneticFluxDensityWaveform[i] - lowerMagneticFluxDensityWaveform[i]));
    }

    double magneticFluxDensityDifference = magneticFluxDensityAcPeak;
    size_t timeout = 0;
    double abs_tol = 0.001;
    while (fabs(magneticFluxDensityDifference) > abs_tol && timeout < 10) {
        size_t mininumPosition =
            std::distance(std::begin(difference), std::min_element(std::begin(difference), std::end(difference)));
        magneticFluxDensityDifference =
            fabs(upperMagneticFluxDensityWaveform[mininumPosition]) - magneticFluxDensityAcPeak;

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

    std::vector<double> cutUpperMagneticFluxDensityWaveform;
    std::vector<double> cutLowerMagneticFluxDensityWaveform;
    for (auto& elem : upperMagneticFluxDensityWaveform) {
        if (elem <= magneticFluxDensityAcPeak && elem >= -magneticFluxDensityAcPeak) {
            cutUpperMagneticFluxDensityWaveform.push_back(elem);
        }
    }

    for (auto& elem : lowerMagneticFluxDensityWaveform) {
        if (elem <= magneticFluxDensityAcPeak && elem >= -magneticFluxDensityAcPeak) {
            cutLowerMagneticFluxDensityWaveform.push_back(elem);
        }
    }

    _hysteresisMinorLoopTop = cutUpperMagneticFluxDensityWaveform;
    _hysteresisMinorLoopBottom = cutLowerMagneticFluxDensityWaveform;

    size_t minimum_length =
        std::min(cutUpperMagneticFluxDensityWaveform.size(), cutLowerMagneticFluxDensityWaveform.size());
    double bhArea = 0;
    for (size_t i = 0; i < minimum_length; i++) {
        bhArea += fabs(cutUpperMagneticFluxDensityWaveform[i] - cutLowerMagneticFluxDensityWaveform[i]) *
                  constants.roshen_magnetic_field_strength_step;
    }

    double hysteresisLossesDensity = bhArea * frequency;

    if (bhArea < 0) {
        throw std::runtime_error("Negative hysteresis losses");
    }

    return hysteresisLossesDensity;
}

double CoreLossesRoshenModel::get_eddy_current_losses_density(CoreWrapper core,
                                                              OperatingPointExcitation excitation,
                                                              double resistivity) {
    auto constants = Constants();
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    magneticFluxDensity = InputsWrapper::standarize_waveform(magneticFluxDensity, frequency);
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }
    double centralColumnArea = core.get_processed_description().value().get_columns()[0].get_area();

    double volumetricLossesIntegration = 0;
    double timeDifference;

    for (size_t i = 0; i < magneticFluxDensityWaveform.size() - 1; ++i) {
        if (magneticFluxDensity.get_waveform().value().get_time()) {
            timeDifference = magneticFluxDensityTime[i + 1] - magneticFluxDensityTime[i];
        }
        else {
            timeDifference = 1 / frequency / constants.number_points_samples_waveforms;
        }
        volumetricLossesIntegration +=
            pow((magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / timeDifference, 2) *
            timeDifference;
    }

    double eddyCurrentLossesDensity =
        centralColumnArea / 8 / std::numbers::pi / resistivity * frequency * volumetricLossesIntegration;

    return eddyCurrentLossesDensity;
}

double CoreLossesRoshenModel::get_excess_eddy_current_losses_density(OperatingPointExcitation excitation,
                                                                     double resistivity,
                                                                     double alphaTimesN0) {
    auto constants = Constants();
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    magneticFluxDensity = InputsWrapper::standarize_waveform(magneticFluxDensity, frequency);
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    double volumetricLossesIntegration = 0;
    double timeDifference;

    for (size_t i = 0; i < magneticFluxDensityWaveform.size() - 1; ++i) {
        if (magneticFluxDensity.get_waveform().value().get_time()) {
            timeDifference = magneticFluxDensityTime[i + 1] - magneticFluxDensityTime[i];
        }
        else {
            timeDifference = 1 / frequency / constants.number_points_samples_waveforms;
        }
        volumetricLossesIntegration +=
            pow(fabs(magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / timeDifference, 3 / 2) *
            timeDifference;
    }

    double excessEddyCurrentLossesDensity = sqrt(alphaTimesN0 / resistivity) * frequency * volumetricLossesIntegration;

    return excessEddyCurrentLossesDensity;
}

std::map<std::string, double> CoreLossesProprietaryModel::get_core_losses(CoreWrapper core,
                                                                   OperatingPointExcitation excitation,
                                                                   double temperature) {

    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    std::map<std::string, double> result;
    auto materialData = core.get_material();
    double volumetricLosses = -1;

    if (materialData.get_manufacturer_info().get_name() == "Micrometals") {
        auto micrometalsData = get_method_data(materialData, "micrometals");
        double a = micrometalsData.get_a().value();
        double b = micrometalsData.get_b().value();
        double c = micrometalsData.get_c().value();
        double d = micrometalsData.get_d().value();
        volumetricLosses = frequency / (a / pow(magneticFluxDensityAcPeak, 3) + b / pow(magneticFluxDensityAcPeak, 2.3) + c / pow(magneticFluxDensityAcPeak, 1.65)) + d * pow(magneticFluxDensityAcPeak, 2) * pow(frequency, 2);
    }

    if (materialData.get_manufacturer_info().get_name() == "Magnetics") {
        auto micrometalsData = get_method_data(materialData, "magnetics");
        double a = micrometalsData.get_a().value();
        double b = micrometalsData.get_b().value();
        double c = micrometalsData.get_c().value();
        volumetricLosses = a * pow(magneticFluxDensityAcPeak, b) * pow(frequency, c);
    }
    result["totalLosses"] = volumetricLosses * effectiveVolume;
    result["totalVolumetricLosses"] = volumetricLosses;

    return result;
}

double CoreLossesSteinmetzModel::get_frequency_from_core_losses(CoreWrapper core,
                                                                SignalDescriptor magneticFluxDensity,
                                                                double temperature,
                                                                double coreLosses) {
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() -
                                       magneticFluxDensity.get_processed().value().get_offset();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();

    OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    double frequency = 100000;

    steinmetzDatum = get_steinmetz_coefficients(core.get_functional_description().get_material(), frequency);
    double alpha = steinmetzDatum.get_alpha();
    double convergeAlpha = 0;

    do {
        double k = steinmetzDatum.get_k();
        alpha = steinmetzDatum.get_alpha();
        double beta = steinmetzDatum.get_beta();
        double volumetricLosses = coreLosses / effectiveVolume / CoreLossesModel::apply_temperature_coefficients(1, steinmetzDatum, temperature);

        frequency = pow(volumetricLosses / (k * pow(magneticFluxDensityAcPeak, beta)), 1 / alpha);
        steinmetzDatum = get_steinmetz_coefficients(core.get_functional_description().get_material(), frequency);
        convergeAlpha = steinmetzDatum.get_alpha();
    }
    while (convergeAlpha != alpha);

    return frequency;
}

SignalDescriptor CoreLossesSteinmetzModel::get_magnetic_flux_density_from_core_losses(CoreWrapper core,
                                                                                              double frequency,
                                                                                              double temperature,
                                                                                              double coreLosses) {
    SignalDescriptor magneticFluxDensity;
    Processed processed;
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(0);
    magneticFluxDensity.set_processed(processed);

    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();

    OpenMagnetics::SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;

    steinmetzDatum = get_steinmetz_coefficients(core.get_functional_description().get_material(), frequency);

    double k = steinmetzDatum.get_k();
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double volumetricLosses = coreLosses / effectiveVolume / CoreLossesModel::apply_temperature_coefficients(1, steinmetzDatum, temperature);

    double magneticFluxDensityAcPeak = pow(volumetricLosses / (k * pow(frequency, alpha)), 1 / beta);

    processed.set_peak(magneticFluxDensityAcPeak);
    processed.set_peak_to_peak(magneticFluxDensityAcPeak * 2);
    magneticFluxDensity.set_processed(processed);
    return magneticFluxDensity;
}

double CoreLossesProprietaryModel::get_frequency_from_core_losses(CoreWrapper core,
                                                                    SignalDescriptor magneticFluxDensity,
                                                                    double temperature,
                                                                    double coreLosses) {

    double frequency = -1;
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    auto materialData = core.get_material();
    double volumetricLosses = coreLosses / effectiveVolume;

    if (materialData.get_manufacturer_info().get_name() == "Micrometals") {
        auto micrometalsData = get_method_data(materialData, "micrometals");
        double a = micrometalsData.get_a().value();
        double b = micrometalsData.get_b().value();
        double c = micrometalsData.get_c().value();
        double d = micrometalsData.get_d().value();
        double equation_a = d * pow(magneticFluxDensityAcPeak, 2);
        double equation_b = 1 / (a / pow(magneticFluxDensityAcPeak, 3) + b / pow(magneticFluxDensityAcPeak, 2.3) + c / pow(magneticFluxDensityAcPeak, 1.65));
        double equation_c = -volumetricLosses;
        frequency = (-equation_b + sqrt(pow(equation_b, 2) - 4 * equation_a * equation_c)) / (2 * equation_a);
    }

    if (materialData.get_manufacturer_info().get_name() == "Magnetics") {
        auto micrometalsData = get_method_data(materialData, "magnetics");
        double a = micrometalsData.get_a().value();
        double b = micrometalsData.get_b().value();
        double c = micrometalsData.get_c().value();
        frequency = pow(volumetricLosses / (a * pow(magneticFluxDensityAcPeak, b)), 1 / c);
    }

    return frequency;
}

double CoreLossesModel::_get_frequency_from_core_losses(CoreWrapper core,
                                                       SignalDescriptor magneticFluxDensity,
                                                       double temperature,
                                                       double coreLosses) {
    double minimumError = DBL_MAX;
    double frequencyMinimumError = -1;
    OperatingPointExcitation operatingPointExcitation;
    operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);


    for (int frequency = 10000; frequency < 2000000; frequency+=5000)
    {
        operatingPointExcitation.set_frequency(frequency);

        auto coreLossesCalculated = get_core_losses(core, operatingPointExcitation, temperature);
        double error = fabs(coreLossesCalculated["totalLosses"] - coreLosses) / coreLosses;
        if (error < minimumError) {
            minimumError = error;
            frequencyMinimumError = frequency;
        }
    }
    return frequencyMinimumError;
}

SignalDescriptor CoreLossesModel::_get_magnetic_flux_density_from_core_losses(CoreWrapper core,
                                                                                     double frequency,
                                                                                     double temperature,
                                                                                     double coreLosses) {

    OperatingPointExcitation operatingPointExcitation;
    SignalDescriptor magneticFluxDensity;
    Processed processed;
    operatingPointExcitation.set_frequency(frequency);
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(0);

    double previousMinimumError = DBL_MAX;
    double minimumError = DBL_MAX;
    SignalDescriptor magneticFluxDensityMinimumError;


    for (int i = 5; i < 1000; i+=5)
    {
        processed.set_peak(double(i) / 1000);
        processed.set_peak_to_peak(2 * double(i) / 1000);
        magneticFluxDensity.set_processed(processed);
        operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);

        auto coreLossesCalculated = get_core_losses(core, operatingPointExcitation, temperature);
        double error = fabs(coreLossesCalculated["totalLosses"] - coreLosses) / coreLosses;
        if (error < minimumError) {
            minimumError = error;
            magneticFluxDensityMinimumError = magneticFluxDensity;
        }
        if (previousMinimumError < error) {
            break;
        }
        else {
            previousMinimumError = error;
        }
    }
    return magneticFluxDensityMinimumError;
}

} // namespace OpenMagnetics
