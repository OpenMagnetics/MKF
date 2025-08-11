#include "physical_models/CoreLosses.h"
#include "physical_models/Resistivity.h"
#include "physical_models/InitialPermeability.h"

#include "processors/Inputs.h"
#include "physical_models/Reluctance.h"
#include "physical_models/MagneticField.h"

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
#include "levmar.h"

std::map<std::string, tk::spline> lossFactorInterps;

namespace OpenMagnetics {

std::vector<CoreLossesModels> CoreLossesModel::get_methods(CoreMaterialDataOrNameUnion material) {
    CoreMaterial materialData;
    // If the material is a string, we have to load its data from the database, unless it is dummy (in order to
    // avoid long loading operatings)
    if (std::holds_alternative<std::string>(material) && std::get<std::string>(material) != "dummy") {
        materialData =
            find_core_material_by_name(std::get<std::string>(material));
    }
    else {
        materialData = std::get<CoreMaterial>(material);
    }

    std::vector<CoreLossesModels> models;
    {
        std::vector<VolumetricCoreLossesMethodType> methods;
        auto volumetricLossesMethodsVariants = materialData.get_volumetric_losses();
        for (auto& volumetricLossesMethodVariant : volumetricLossesMethodsVariants) {
            auto volumetricLossesMethods = volumetricLossesMethodVariant.second;
            for (auto& volumetricLossesMethod : volumetricLossesMethods) {
                if (std::holds_alternative<CoreLossesMethodData>(volumetricLossesMethod)) {
                    auto methodData = std::get<CoreLossesMethodData>(volumetricLossesMethod);
                    methods.push_back(methodData.get_method());
                }
            }
        }

        if (std::count(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::STEINMETZ)) {
            models.push_back(CoreLossesModels::STEINMETZ);
            models.push_back(CoreLossesModels::IGSE);
            models.push_back(CoreLossesModels::BARG);
            models.push_back(CoreLossesModels::ALBACH);
            models.push_back(CoreLossesModels::MSE);
        }
        if (std::count(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::ROSHEN)) {
            models.push_back(CoreLossesModels::ROSHEN);
        }
        if (std::count(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::MAGNETICS) || std::count(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::MICROMETALS) || std::count(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::POCO) || std::count(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::TDG)) {
            models.push_back(CoreLossesModels::PROPRIETARY);
        }
        if (std::count(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::LOSS_FACTOR)) {
            models.push_back(CoreLossesModels::LOSS_FACTOR);
        }
    }

    if (materialData.get_mass_losses()) {
        std::vector<MassCoreLossesMethodType> methods;
        auto massLossesMethodsVariants = materialData.get_mass_losses().value();
        for (auto& massLossesMethodVariant : massLossesMethodsVariants) {
            auto massLossesMethods = massLossesMethodVariant.second;
            for (auto& massLossesMethod : massLossesMethods) {
                if (std::holds_alternative<MagneticsCoreLossesMethodData>(massLossesMethod)) {
                    auto methodData = std::get<MagneticsCoreLossesMethodData>(massLossesMethod);
                    methods.push_back(methodData.get_method());
                }
            }
        }

        if (std::count(methods.begin(), methods.end(), MassCoreLossesMethodType::MAGNETEC)) {
            models.push_back(CoreLossesModels::PROPRIETARY);
        }
    }

    return models;
}

std::shared_ptr<CoreLossesModel> CoreLosses::get_core_losses_model(std::string materialName) {
    std::shared_ptr<CoreLossesModel> coreLossesModelForMaterial = nullptr;

    auto availableMethodsForMaterial = CoreLossesModel::get_methods(materialName);
    for (auto& [modelName, coreLossesModel] : _coreLossesModels) {
        if (std::find(availableMethodsForMaterial.begin(), availableMethodsForMaterial.end(), modelName) != availableMethodsForMaterial.end()) {
            coreLossesModelForMaterial = coreLossesModel;
            break;
        }
    }
    if (coreLossesModelForMaterial == nullptr) {
        throw std::runtime_error("No model found for material: " + materialName);
    }

    return coreLossesModelForMaterial;
}


CoreLossesOutput CoreLosses::calculate_core_losses(Core core, OperatingPointExcitation excitation, double temperature) {
    auto coreLossesModelForMaterial = get_core_losses_model(core.get_material_name());

    CoreLossesOutput coreLossesOutput = coreLossesModelForMaterial->get_core_losses(core, excitation, temperature);
    return coreLossesOutput;
}
double CoreLosses::get_core_volumetric_losses(CoreMaterial coreMaterial, OperatingPointExcitation excitation, double temperature){
    auto coreLossesModelForMaterial = get_core_losses_model(coreMaterial.get_name());

    double coreVolumetricLosses = coreLossesModelForMaterial->get_core_volumetric_losses(coreMaterial, excitation, temperature);
    return coreVolumetricLosses;
}

double CoreLosses::get_core_losses_series_resistance(Core core, double frequency, double temperature, double magnetizingInductance){
    auto coreLossesModelForMaterial = get_core_losses_model(core.get_material_name());

    double coreLossesSeriesResistance = coreLossesModelForMaterial->get_core_losses_series_resistance(core, frequency, temperature, magnetizingInductance);
    return coreLossesSeriesResistance;
}


std::shared_ptr<CoreLossesModel> CoreLossesModel::factory(std::map<std::string, std::string> models) {
    return factory(magic_enum::enum_cast<CoreLossesModels>(models["coreLosses"]).value());
}

std::shared_ptr<CoreLossesModel> CoreLossesModel::factory(json models) {
    std::string model = models["coreLosses"];
    return factory(magic_enum::enum_cast<CoreLossesModels>(model).value());
}

std::shared_ptr<CoreLossesModel> CoreLossesModel::factory(CoreLossesModels modelName) {
    if (modelName == CoreLossesModels::STEINMETZ) {
        return std::make_shared<CoreLossesSteinmetzModel>();
    }
    else if (modelName == CoreLossesModels::IGSE) {
        return std::make_shared<CoreLossesIGSEModel>();
    }
    else if (modelName == CoreLossesModels::MSE) {
        return std::make_shared<CoreLossesMSEModel>();
    }
    else if (modelName == CoreLossesModels::NSE) {
        return std::make_shared<CoreLossesNSEModel>();
    }
    else if (modelName == CoreLossesModels::ALBACH) {
        return std::make_shared<CoreLossesAlbachModel>();
    }
    else if (modelName == CoreLossesModels::BARG) {
        return std::make_shared<CoreLossesBargModel>();
    }
    else if (modelName == CoreLossesModels::ROSHEN) {
        return std::make_shared<CoreLossesRoshenModel>();
    }
    else if (modelName == CoreLossesModels::PROPRIETARY) {
        return std::make_shared<CoreLossesProprietaryModel>();
    }
    else if (modelName == CoreLossesModels::LOSS_FACTOR) {
        return std::make_shared<CoreLossesLossFactorModel>();
    }

    else
        throw std::runtime_error("Unknown Core losses mode, available options are: {STEINMETZ, IGSE, BARG, ALBACH, "
                                 "ROSHEN, OUYANG, NSE, MSE, PROPRIETARY, LOSS_FACTOR}");
}

std::vector<VolumetricLossesPoint> CoreLossesModel::get_volumetric_losses_data(CoreMaterial materialData) {
    auto volumetricLossesMethodsVariants = materialData.get_volumetric_losses();

    for (auto& volumetricLossesMethodVariant : volumetricLossesMethodsVariants) {
        if (volumetricLossesMethodVariant.first != "default") {
            continue;
        }
        auto volumetricLossesMethods = volumetricLossesMethodVariant.second;
        for (auto& volumetricLossesMethod : volumetricLossesMethods) {
            if (std::holds_alternative<std::vector<VolumetricLossesPoint>>(volumetricLossesMethod)) {
                return std::get<std::vector<VolumetricLossesPoint>>(volumetricLossesMethod);
            }
        }
    }
    return std::vector<VolumetricLossesPoint>();
}


CoreLossesMethodData CoreLossesModel::get_method_data(CoreMaterial materialData, std::string method) {
    auto volumetricLossesMethodsVariants = materialData.get_volumetric_losses();
    std::string methodUpper = method;
    std::transform(methodUpper.begin(), methodUpper.end(), methodUpper.begin(), ::toupper);
    for (auto& volumetricLossesMethodVariant : volumetricLossesMethodsVariants) {
        if (volumetricLossesMethodVariant.first != "default") {
            continue;
        }
        auto volumetricLossesMethods = volumetricLossesMethodVariant.second;
        for (auto& volumetricLossesMethod : volumetricLossesMethods) {
            if (std::holds_alternative<CoreLossesMethodData>(volumetricLossesMethod)) {
                auto methodData = std::get<CoreLossesMethodData>(volumetricLossesMethod);
                std::string methodDataNameString = std::string{magic_enum::enum_name(methodData.get_method())};

                if (methodDataNameString == methodUpper) {
                    return methodData;
                }
            }
        }
    }
    throw std::runtime_error("Material " + materialData.get_name() + " does not have method:" + method);
}

SteinmetzCoreLossesMethodRangeDatum CoreLossesModel::get_steinmetz_coefficients(CoreMaterialDataOrNameUnion material, double frequency) {
    CoreMaterial materialData;
    // If the material is a string, we have to load its data from the database, unless it is dummy (in order to avoid
    // long loading operatings)
    if (std::holds_alternative<std::string>(material) && std::get<std::string>(material) != "dummy") {
        materialData = find_core_material_by_name(std::get<std::string>(material));
    }
    else {
        materialData = std::get<CoreMaterial>(material);
    }

    auto volumetricLossesMethodsVariants = materialData.get_volumetric_losses();

    auto steinmetzData = CoreLossesModel::get_method_data(materialData, "steinmetz");
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


double steinmetz_equation_with_temperature_and_log(double x[], double frequency, double magneticFluxDensityAcPeak, double temperature) {
    double temperatureCoefficient = x[3] - x[4] * temperature + x[5] * pow(temperature, 2);
    if (temperatureCoefficient < 0) {
        return x[0] + frequency * x[1] + magneticFluxDensityAcPeak * x[2];
    }
    else {
        return x[0] + frequency * x[1] + magneticFluxDensityAcPeak * x[2] + log10(temperatureCoefficient);
    }
}
double steinmetz_equation_with_temperature_and_log(double x[], double logK, double alpha, double beta, double frequency, double magneticFluxDensityAcPeak, double temperature) {
    double temperatureCoefficient = x[0] - x[1] * temperature + x[2] * pow(temperature, 2);
    if (temperatureCoefficient < 0) {
        return logK + frequency * alpha + magneticFluxDensityAcPeak * beta;
    }
    else {
        return logK + frequency * alpha + magneticFluxDensityAcPeak * beta + log10(temperatureCoefficient);
    }
}


double steinmetz_equation_log(double x[], double frequency, double magneticFluxDensityAcPeak) {
    if (x[1] < 0 || x[2] < 0) {
        return 0;
    }
    return x[0] + frequency * x[1] + magneticFluxDensityAcPeak * x[2];
}

double steinmetz_equation_with_temperature(double x[], double frequency, double magneticFluxDensityAcPeak, double temperature) {
    if (pow(10, x[0]) < 0) {
            return DBL_MAX;
    }
    if (pow(10, x[1]) < 1) {
            return DBL_MAX;
    }
    if (pow(10, x[2]) < 2) {
            return DBL_MAX;
    }

    return x[0] * pow(frequency, x[1]) * pow(magneticFluxDensityAcPeak, x[2]) * (x[3] - x[4] * temperature + x[5] * pow(temperature, 2));
}

double steinmetz_equation(double x[], double frequency, double magneticFluxDensityAcPeak) {
    // for(int i=0; i<10; ++i) {
    //     if (x[i] < 0) {
    //         return 0;
    //     }
    // }

    return x[0] * pow(frequency, x[1]) * pow(magneticFluxDensityAcPeak, x[2]);
}

void steinmetz_equation_func(double *p, double *x, int m, int n, void *data) {
    double* aux = static_cast <double*> (data);

    for(int i=0; i<n; ++i) {
        auto frequency = aux[3 + 2 * i];
        auto magneticFluxDensityAcPeak = aux[3 + 2 * i + 1];
        x[i]=steinmetz_equation_log(p, frequency, magneticFluxDensityAcPeak);
    }

    // x[0] = pow(10, x[0]);
}

void steinmetz_equation_with_temperature_func(double *p, double *x, int m, int n, void *data) {
    double* aux = static_cast <double*> (data);

    for(int i=0; i<n; ++i) {
        auto frequency = aux[3 + 3 * i];
        auto magneticFluxDensityAcPeak = aux[3 + 3 * i + 1];
        auto temperature = aux[3 + 3 * i + 2];
        x[i]=steinmetz_equation_with_temperature_and_log(p, frequency, magneticFluxDensityAcPeak, temperature);
    }
}


void steinmetz_equation_first_no_temperature_func(double *p, double *x, int m, int n, void *data) {
    double* aux = static_cast <double*> (data);

    for(int i=0; i<n; ++i) {
        auto frequency = aux[2 * i];
        auto magneticFluxDensityAcPeak = aux[2 * i + 1];
        x[i]=steinmetz_equation_log(p, frequency, magneticFluxDensityAcPeak);
    }
}

void steinmetz_equation_first_only_temperature_func(double *p, double *x, int m, int n, void *data) {
    double* aux = static_cast <double*> (data);
    double logK = aux[0];
    double alpha = aux[1];
    double beta = aux[2];

    for(int i=0; i<n; ++i) {
        auto frequency = aux[3 + 3 * i];
        auto magneticFluxDensityAcPeak = aux[3 + 3 * i + 1];
        auto temperature = aux[3 + 3 * i + 2];
        x[i]=steinmetz_equation_with_temperature_and_log(p, logK, alpha, beta, frequency, magneticFluxDensityAcPeak, temperature);
    }
}

std::pair<std::vector<SteinmetzCoreLossesMethodRangeDatum>, std::vector<double>> CoreLossesSteinmetzModel::calculate_steinmetz_coefficients(std::vector<VolumetricLossesPoint> volumetricLosses, std::vector<std::pair<double, double>> ranges) {
    std::vector<double> bestErrorPerRange;
    std::vector<SteinmetzCoreLossesMethodRangeDatum> steinmetzCoefficientsPerRange;

    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

    double lmInitMu = 1e-03;
    double lmStopThresh = 1e-25;
    double lmDiffDelta = 1e-19;

    opts[0]= lmInitMu;
    opts[1]= lmStopThresh;
    opts[2]= lmStopThresh;
    opts[3]= lmStopThresh;
    opts[4]= lmDiffDelta;

    size_t loopIterations = 6;

    std::vector<double> distinctTemperatures;
    for (size_t index = 0; index < volumetricLosses.size(); ++index) {
        auto temperature = volumetricLosses[index].get_temperature();
        if (std::find(distinctTemperatures.begin(), distinctTemperatures.end(), temperature) == distinctTemperatures.end()) {
            distinctTemperatures.push_back(temperature);
        }
    }
    size_t numberInputs = distinctTemperatures.size() > 1? 3 : 2;

    size_t numberUnknowns = numberInputs == 3? 6 : 3;

    std::vector<std::vector<VolumetricLossesPoint>> volumetricLossesChunks;

    if (ranges.size() > 1) {
        for (size_t index = 0; index < ranges.size(); ++index) {
            volumetricLossesChunks.push_back(std::vector<VolumetricLossesPoint>());
        }
        for (auto volumetricLoss : volumetricLosses) {
            auto frequency = volumetricLoss.get_magnetic_flux_density().get_frequency();
            for (size_t chunkIndex = 0; chunkIndex < ranges.size(); ++chunkIndex) {
                if (ranges[chunkIndex].first * 0.8 <= frequency && frequency <= ranges[chunkIndex].second * 1.2) {
                    volumetricLossesChunks[chunkIndex].push_back(volumetricLoss);
                }
            }
        }
    }
    else {
        volumetricLossesChunks = {volumetricLosses};
    }

    // We remove empty or too small chunks
    bool continueCleaning = true;
    while (continueCleaning) {
        continueCleaning = false;        
        for (size_t chunkIndex = 0; chunkIndex < volumetricLossesChunks.size(); ++chunkIndex) {
            auto volumetricLossesChunk = volumetricLossesChunks[chunkIndex];

            if (volumetricLossesChunk.size() <= numberUnknowns) {
                if (chunkIndex == 0 && volumetricLossesChunks.size() == 1) {
                    if (volumetricLossesChunks[0].size() > 3) {
                        numberInputs = 2;
                        numberUnknowns = 3;
                        break;
                    }
                    else {
                        throw std::runtime_error("Too few points");
                    }
                }
                else if (chunkIndex == 0) {
                    for (auto volumetricLosses : volumetricLossesChunk) {
                        volumetricLossesChunks[1].push_back(volumetricLosses);
                        ranges[1] = {ranges[0].first, ranges[1].second};
                    }
                    ranges.erase(ranges.begin() + chunkIndex);
                    volumetricLossesChunks.erase(volumetricLossesChunks.begin() + chunkIndex);
                    continueCleaning = true;        
                    break;
                }
                else {
                    for (auto volumetricLosses : volumetricLossesChunk) {
                        volumetricLossesChunks[chunkIndex - 1].push_back(volumetricLosses);
                        ranges[chunkIndex - 1] = {ranges[chunkIndex - 1].first, ranges[chunkIndex].second};
                    }
                    ranges.erase(ranges.begin() + chunkIndex);
                    volumetricLossesChunks.erase(volumetricLossesChunks.begin() + chunkIndex);
                    continueCleaning = true;        
                    break;
                }
            }
        }
    }

    for (size_t chunkIndex = 0; chunkIndex < volumetricLossesChunks.size(); ++chunkIndex) {
        double bestError = DBL_MAX;
        double initialState = 10;
        std::vector<double> bestCoefficients;
        auto volumetricLossesChunk = volumetricLossesChunks[chunkIndex];
        for (size_t loopIndex = 0; loopIndex < loopIterations; ++loopIndex) {
            size_t numberElements = volumetricLossesChunk.size();
            size_t numberElements100C = 0;

            double* volumetricLossesArray = new double[numberElements];
            // double volumetricLossesArray[numberElements];

            for (size_t index = 0; index < numberElements; ++index) {
                volumetricLossesArray[index] = log10(volumetricLossesChunk[index].get_value());
                auto temperature = volumetricLossesChunk[index].get_temperature();
                if (temperature >= 90 && temperature <= 110) {
                    numberElements100C++;
                }
                if (std::find(distinctTemperatures.begin(), distinctTemperatures.end(), temperature) == distinctTemperatures.end()) {
                    distinctTemperatures.push_back(temperature);
                }
            }

            double* coefficients = new double[numberUnknowns];
            // double coefficients[numberUnknowns];
            for (size_t index = 0; index < numberUnknowns; ++index) {
                coefficients[index] = initialState;
            }

            double* volumetricLossesInputs = new double[3 + numberElements * numberInputs];
            // double volumetricLossesInputs[numberElements * numberInputs];
            for (size_t index = 0; index < numberElements; ++index) {
                volumetricLossesInputs[3 + numberInputs * index] = log10(volumetricLossesChunk[index].get_magnetic_flux_density().get_frequency());
                volumetricLossesInputs[3 + numberInputs * index + 1] = log10(volumetricLossesChunk[index].get_magnetic_flux_density().get_magnetic_flux_density()->get_processed()->get_peak().value());
                if (numberInputs == 3) {
                    volumetricLossesInputs[3 + numberInputs * index + 2] = volumetricLossesChunk[index].get_temperature();
                }
            }

            if (numberInputs == 3 && numberElements100C >= 3) {
                double* tempCoefficients = new double[3];
                for (size_t index = 0; index < 3; ++index) {
                    tempCoefficients[index] = initialState;
                }

                double* tempVolumetricLossesInputs = new double[numberElements100C * 2];
                double* tempVolumetricLossesArray = new double[numberElements100C];
                // double tempVolumetricLossesInputs[numberElements * numberInputs];
                size_t index100C = 0;
                for (size_t index = 0; index < numberElements; ++index) {
                    if (volumetricLossesChunk[index].get_temperature() >= 90 && volumetricLossesChunk[index].get_temperature() <= 110) {
                        tempVolumetricLossesArray[index100C] = log10(volumetricLossesChunk[index].get_value());
                        tempVolumetricLossesInputs[2 * index100C] = log10(volumetricLossesChunk[index].get_magnetic_flux_density().get_frequency());
                        tempVolumetricLossesInputs[2 * index100C + 1] = log10(volumetricLossesChunk[index].get_magnetic_flux_density().get_magnetic_flux_density()->get_processed()->get_peak().value());
                        index100C++;
                    }
                }

                dlevmar_dif(steinmetz_equation_first_no_temperature_func, tempCoefficients, tempVolumetricLossesArray, 3, numberElements100C, 10000, opts, info, NULL, NULL, static_cast<void*>(tempVolumetricLossesInputs));
                coefficients[0] = tempCoefficients[0];
                coefficients[1] = tempCoefficients[1];
                coefficients[2] = tempCoefficients[2];
                volumetricLossesInputs[0] = tempCoefficients[0];
                volumetricLossesInputs[1] = tempCoefficients[1];
                volumetricLossesInputs[2] = tempCoefficients[2];

                for (size_t index = 0; index < 3; ++index) {
                    tempCoefficients[index] = initialState;
                }
                dlevmar_dif(steinmetz_equation_first_only_temperature_func, tempCoefficients, volumetricLossesArray, 3, numberElements, 10000, opts, info, NULL, NULL, static_cast<void*>(volumetricLossesInputs));
                coefficients[3] = tempCoefficients[0];
                coefficients[4] = tempCoefficients[1];
                coefficients[5] = tempCoefficients[2];
            }
            else {
                dlevmar_dif(steinmetz_equation_with_temperature_func, coefficients, volumetricLossesArray, numberUnknowns, numberElements, 10000, opts, info, NULL, NULL, static_cast<void*>(volumetricLossesInputs));
            }

            double errorAverage = 0;
            for (size_t index = 0; index < numberElements; ++index) {
                auto frequency = volumetricLossesInputs[3 + numberInputs * index];
                auto magneticFluxDensityAcPeak = volumetricLossesInputs[3 + numberInputs * index + 1];
                double modeledVolumetricLosses = 0;
                if (numberInputs == 3) {
                    auto temperature = volumetricLossesInputs[3 + numberInputs * index + 2];
                    modeledVolumetricLosses = steinmetz_equation_with_temperature_and_log(coefficients, frequency, magneticFluxDensityAcPeak, temperature);
                }
                else {
                    modeledVolumetricLosses = steinmetz_equation_log(coefficients, frequency, magneticFluxDensityAcPeak);
                }
                double error = fabs(pow(10, volumetricLossesArray[index]) - pow(10, modeledVolumetricLosses)) / pow(10, volumetricLossesArray[index]);
                errorAverage += error;
            }


            errorAverage /= numberElements;
            initialState /= 10;

            if (errorAverage < bestError) {
                bestError = errorAverage;
                bestCoefficients.clear();
                for (size_t index = 0; index < numberUnknowns; ++index) {
                    bestCoefficients.push_back(coefficients[index]);
                }
                bestCoefficients[0] = pow(10, bestCoefficients[0]); 
            }

            delete[] volumetricLossesArray;
            delete[] coefficients;
            delete[] volumetricLossesInputs;
        }
        SteinmetzCoreLossesMethodRangeDatum steinmetzCoreLossesMethodRangeDatum;
        steinmetzCoreLossesMethodRangeDatum.set_k(bestCoefficients[0]);
        steinmetzCoreLossesMethodRangeDatum.set_alpha(bestCoefficients[1]);
        steinmetzCoreLossesMethodRangeDatum.set_beta(bestCoefficients[2]);
        if (numberInputs == 3) {
            steinmetzCoreLossesMethodRangeDatum.set_ct0(bestCoefficients[3]);
            steinmetzCoreLossesMethodRangeDatum.set_ct1(bestCoefficients[4]);
            steinmetzCoreLossesMethodRangeDatum.set_ct2(bestCoefficients[5]);
        }
        steinmetzCoreLossesMethodRangeDatum.set_minimum_frequency(ranges[chunkIndex].first);
        steinmetzCoreLossesMethodRangeDatum.set_maximum_frequency(ranges[chunkIndex].second);
        steinmetzCoefficientsPerRange.push_back(steinmetzCoreLossesMethodRangeDatum);
        bestErrorPerRange.push_back(bestError);
    }

    return {steinmetzCoefficientsPerRange, bestErrorPerRange};
};

CoreLossesOutput CoreLossesSteinmetzModel::get_core_losses(Core core,
                                                  OperatingPointExcitation excitation,
                                                  double temperature) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    auto material = core.resolve_material();

    CoreLossesOutput result;
    result.set_magnetic_flux_density(magneticFluxDensity);
    result.set_method_used(_modelName);
    result.set_origin(ResultOrigin::SIMULATION);
    result.set_temperature(temperature);

    if (usesVolumetricLosses(material)) {
        auto volumetricLosses = get_core_volumetric_losses(material, excitation, temperature);
        result.set_core_losses(volumetricLosses * effectiveVolume);
        result.set_volumetric_losses(volumetricLosses);
    }
    else {
        auto massLosses = get_core_mass_losses(material, excitation, temperature);
        result.set_core_losses(massLosses * effectiveVolume);
        result.set_mass_losses(massLosses);
    }


    return result;
};

double CoreLossesSteinmetzModel::get_core_volumetric_losses(CoreMaterial coreMaterial,
                                                            OperatingPointExcitation excitation,
                                                            double temperature) {
    if (!excitation.get_magnetic_flux_density()) {
        throw std::runtime_error("Missing magnetizing flux density in excitation");
    }
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = Inputs::get_switching_frequency(excitation);
    double mainHarmonicMagneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() - magneticFluxDensity.get_processed().value().get_offset();
    double magneticFluxDensityAcPeak = Inputs::get_magnetic_flux_density_peak(excitation, frequency) - magneticFluxDensity.get_processed().value().get_offset();

    magneticFluxDensity = Inputs::standardize_waveform(magneticFluxDensity, frequency);


    SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(coreMaterial, frequency);
    }
    double k = steinmetzDatum.get_k();
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double volumetricLosses;


    if (beta > 2) {
        volumetricLosses = k * pow(frequency, alpha) * pow(mainHarmonicMagneticFluxDensityAcPeak, beta - 2) * pow(magneticFluxDensityAcPeak, 2);
    }
    else {
        volumetricLosses = k * pow(frequency, alpha) * pow(magneticFluxDensityAcPeak, beta);
    }

    return CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);
};

double CoreLossesModel::get_core_losses_series_resistance(Core core,
                                                          double frequency,
                                                          double temperature,
                                                          double magnetizingInductance) {

    double virtualCurrentRms = 1;
    auto coreMaterial = core.resolve_material();
    double effectiveArea = core.get_processed_description().value().get_effective_parameters().get_effective_area();

    double initialPermeability = InitialPermeability::get_initial_permeability(coreMaterial, temperature);
    auto reluctanceModel = ReluctanceModel::factory();
    auto reluctance = reluctanceModel->get_core_reluctance(core, initialPermeability).get_core_reluctance();

    int64_t numberTurnsPrimary = sqrt(magnetizingInductance * reluctance);
    auto operatingPoint = Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, magnetizingInductance, temperature, {}, {virtualCurrentRms * sqrt(2)});
    operatingPoint = Inputs::process_operating_point(operatingPoint, magnetizingInductance);
    auto excitation = operatingPoint.get_excitations_per_winding()[0];
    auto magneticFlux = MagneticField::calculate_magnetic_flux(excitation.get_magnetizing_current().value(), reluctance, numberTurnsPrimary);
    auto magneticFluxDensity = MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
    excitation.set_magnetic_flux_density(magneticFluxDensity);

    double coreLosses = get_core_losses(core, excitation, temperature).get_core_losses();

    auto seriesResistance = coreLosses / pow(virtualCurrentRms, 2);

    return seriesResistance;
}

double CoreLossesIGSEModel::get_ki(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum) {
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double k = steinmetzDatum.get_k();
    std::vector<double> theta_vect(settings->get_inputs_number_points_sampled_waveforms());
    std::generate(theta_vect.begin(), theta_vect.end(), [n = 0]() mutable {
        return n++ * 2 * std::numbers::pi / settings->get_inputs_number_points_sampled_waveforms();
    });
    double theta_int = 0;

    for (auto& theta : theta_vect) {
        theta_int += pow(fabs(cos(theta)), alpha) * pow(2, beta - alpha) * 2 * std::numbers::pi /
                     settings->get_inputs_number_points_sampled_waveforms();
    }

    double ki = k / (pow(2 * std::numbers::pi, alpha - 1) * theta_int);
    return ki;
}

double CoreLossesIGSEModel::get_core_volumetric_losses(CoreMaterial coreMaterial,
                                                       OperatingPointExcitation excitation,
                                                       double temperature) {

    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = Inputs::get_switching_frequency(excitation);
    double magneticFluxDensityAcPeakToPeak = Inputs::get_magnetic_flux_density_peak_to_peak(excitation, frequency);

    magneticFluxDensity = Inputs::standardize_waveform(magneticFluxDensity, frequency);
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(coreMaterial, frequency);
    }

    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double ki = get_ki(steinmetzDatum);

    double volumetricLossesSum = 0;
    double timeDifference;
    size_t numberPoints = magneticFluxDensityWaveform.size();

    if (frequency / excitation.get_frequency() > 1) {
        numberPoints = round(numberPoints / (frequency / excitation.get_frequency()));
    }

    for (size_t i = 0; i < numberPoints - 1; ++i) {
        if (magneticFluxDensity.get_waveform().value().get_time()) {
            timeDifference = magneticFluxDensityTime[i + 1] - magneticFluxDensityTime[i];
        }
        else {
            timeDifference = 1 / frequency / settings->get_inputs_number_points_sampled_waveforms();
        }
        volumetricLossesSum +=
            pow(fabs((magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / timeDifference), alpha) *
            timeDifference;
    }
 
    // double volumetricLosses = ki * pow(mainHarmonicMagneticFluxDensityPeakToPeak, beta - alpha) * frequency * volumetricLossesSum;
    double volumetricLosses = ki * pow(magneticFluxDensityAcPeakToPeak, beta - alpha) * frequency * volumetricLossesSum;
    return CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);
}

double CoreLossesAlbachModel::get_core_volumetric_losses(CoreMaterial coreMaterial,
                                                         OperatingPointExcitation excitation,
                                                         double temperature) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = Inputs::get_switching_frequency(excitation);
    double mainHarmonicMagneticFluxDensityPeak = magneticFluxDensity.get_processed().value().get_peak().value();
    double magneticFluxDensityPeakToPeak = Inputs::get_magnetic_flux_density_peak_to_peak(excitation, frequency);
    double magneticFluxDensityAcPeak = Inputs::get_magnetic_flux_density_peak(excitation, frequency) - magneticFluxDensity.get_processed().value().get_offset();

    magneticFluxDensity = Inputs::standardize_waveform(magneticFluxDensity, frequency);
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(coreMaterial, frequency);
    }

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
            timeDifference = 1 / frequency / settings->get_inputs_number_points_sampled_waveforms();
        }
        equivalentSinusoidalFrequency +=
            pow((magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / magneticFluxDensityPeakToPeak,
                2) /
            timeDifference;
    }

    equivalentSinusoidalFrequency = 2 / pow(std::numbers::pi, 2) * equivalentSinusoidalFrequency;


    double volumetricLosses;
    if (beta > 2) {
        volumetricLosses = k * frequency * pow(equivalentSinusoidalFrequency, alpha - 1) * pow(mainHarmonicMagneticFluxDensityPeak, beta - 2) * pow(magneticFluxDensityAcPeak, 2);
    }
    else {
        volumetricLosses = k * frequency * pow(equivalentSinusoidalFrequency, alpha - 1) * pow(magneticFluxDensityAcPeak, beta);
    }

    return CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);
}

double CoreLossesMSEModel::get_core_volumetric_losses(CoreMaterial coreMaterial,
                                                      OperatingPointExcitation excitation,
                                                      double temperature) {
    double frequency = Inputs::get_switching_frequency(excitation);
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    magneticFluxDensity = Inputs::standardize_waveform(magneticFluxDensity, frequency);
    double mainHarmonicMagneticFluxDensityPeak = magneticFluxDensity.get_processed().value().get_peak().value();
    double magneticFluxDensityPeakToPeak = Inputs::get_magnetic_flux_density_peak_to_peak(excitation, frequency);
    double magneticFluxDensityAcPeak = Inputs::get_magnetic_flux_density_peak(excitation, frequency) -
                                       magneticFluxDensity.get_processed().value().get_offset();


    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(coreMaterial, frequency);
    }

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
            timeDifference = 1 / frequency / settings->get_inputs_number_points_sampled_waveforms();
        }
        equivalentSinusoidalFrequency +=
            pow((magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / timeDifference, 2) *
            timeDifference;
    }

    equivalentSinusoidalFrequency =
        2 / pow(std::numbers::pi, 2) / pow(magneticFluxDensityPeakToPeak, 2) * equivalentSinusoidalFrequency;

    double volumetricLosses;
    if (beta > 2) {
        volumetricLosses = k * frequency * pow(equivalentSinusoidalFrequency, alpha - 1) * pow(mainHarmonicMagneticFluxDensityPeak, beta - 2) * pow(magneticFluxDensityAcPeak, 2);
    }
    else {
        volumetricLosses = k * frequency * pow(equivalentSinusoidalFrequency, alpha - 1) * pow(magneticFluxDensityAcPeak, beta);
    }

    return CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);
}

double CoreLossesNSEModel::get_kn(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum) {
    double alpha = steinmetzDatum.get_alpha();
    double k = steinmetzDatum.get_k();
    std::vector<double> theta_vect(settings->get_inputs_number_points_sampled_waveforms());
    std::generate(theta_vect.begin(), theta_vect.end(), [n = 0]() mutable {
        return n++ * 2 * std::numbers::pi / settings->get_inputs_number_points_sampled_waveforms();
    });
    double theta_int = 0;

    for (auto& theta : theta_vect) {
        theta_int += pow(fabs(cos(theta)), alpha) * 2 * std::numbers::pi / settings->get_inputs_number_points_sampled_waveforms();
    }

    double ki = k / (pow(2 * std::numbers::pi, alpha - 1) * theta_int);
    return ki;
}

double CoreLossesNSEModel::get_core_volumetric_losses(CoreMaterial coreMaterial,
                                                      OperatingPointExcitation excitation,
                                                      double temperature) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = Inputs::get_switching_frequency(excitation);
    double mainHarmonicMagneticFluxDensityPeak = magneticFluxDensity.get_processed().value().get_peak().value();

    magneticFluxDensity = Inputs::standardize_waveform(magneticFluxDensity, frequency);
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(coreMaterial, frequency);
    }

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
            timeDifference = 1 / frequency / settings->get_inputs_number_points_sampled_waveforms();
        }
        volumetricLossesSum +=
            pow(abs((magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / timeDifference), alpha) *
            timeDifference;
    }

    double volumetricLosses = kn * pow(mainHarmonicMagneticFluxDensityPeak, beta - alpha) * frequency * volumetricLossesSum;
    return CoreLossesModel::apply_temperature_coefficients(volumetricLosses, steinmetzDatum, temperature);
}

double get_plateau_duty_cycle(std::vector<double> data) {
    double maxValue = *std::max_element(data.begin(), data.end());
    int numberPlateauPoints = 0;
    for (auto datum : data) {
        if (fabs(maxValue - datum) / maxValue < 0.02) {
            numberPlateauPoints++;
        }
    }
    auto onPoints = (data.size() / 2 - numberPlateauPoints) / data.size();
    return onPoints;
}

double CoreLossesBargModel::get_core_volumetric_losses(CoreMaterial coreMaterial,
                                                      OperatingPointExcitation excitation,
                                                      double temperature) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = Inputs::get_switching_frequency(excitation);
    magneticFluxDensity = Inputs::standardize_waveform(magneticFluxDensity, frequency);
    double mainHarmonicMagneticFluxDensityPeak = magneticFluxDensity.get_processed().value().get_peak().value();
    double magneticFluxDensityAcPeak = Inputs::get_magnetic_flux_density_peak(excitation, frequency) -
                                       magneticFluxDensity.get_processed().value().get_offset();


    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }

    SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    if (is_steinmetz_datum_loaded()) {
        steinmetzDatum = get_steinmetz_datum();
    }
    else {
        steinmetzDatum = get_steinmetz_coefficients(coreMaterial, frequency);
    }

    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double k = steinmetzDatum.get_k();
    double dutyCycle = get_plateau_duty_cycle(magneticFluxDensity.get_waveform().value().get_data());
    // double dutyCycle = magneticFluxDensity.get_processed().value().get_duty_cycle().value();

    double lossesFrameT1;
    if (beta > 2) {
        lossesFrameT1 = std::numbers::pi / 4 * k * pow(frequency, alpha) * pow(mainHarmonicMagneticFluxDensityPeak, beta - 2) * pow(magneticFluxDensityAcPeak, 2);
    }
    else {
        lossesFrameT1 = std::numbers::pi / 4 * k * pow(frequency, alpha) * pow(magneticFluxDensityAcPeak, beta);
    }

    lossesFrameT1 = CoreLossesModel::apply_temperature_coefficients(lossesFrameT1, steinmetzDatum, temperature);
    std::vector<double> plateauDutyCycleValues = {0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5};
    std::vector<double> factorValues = {1.45,  1.4,  1.35, 1.275, 1.25,  1.2,  1.15, 1.075, 1};

    tk::spline interp(plateauDutyCycleValues, factorValues, tk::spline::cspline_hermite, true);
    double dutyCycleFactor = std::max(1., interp(dutyCycle));

    return dutyCycleFactor * lossesFrameT1;
} 
 
CoreLossesOutput CoreLossesRoshenModel::get_core_losses(Core core,
                                                        OperatingPointExcitation excitation,
                                                        double temperature) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    auto parameters = get_roshen_parameters(core, excitation, temperature);
    double hysteresisVolumetricLosses = get_hysteresis_losses_density(parameters, excitation);
    double eddyCurrentsVolumetricLosses = get_eddy_current_losses_density(core, excitation, parameters["resistivity"]);
    double excessEddyCurrentsVolumetricLosses = 0;
    if (parameters.count("excessLossesCoefficient")) {
        excessEddyCurrentsVolumetricLosses = get_excess_eddy_current_losses_density(
            excitation, parameters["resistivity"], parameters["excessLossesCoefficient"]);
    }
    double volumetricLosses = hysteresisVolumetricLosses + eddyCurrentsVolumetricLosses + excessEddyCurrentsVolumetricLosses;

    CoreLossesOutput result;
    result.set_core_losses(volumetricLosses * effectiveVolume);
    result.set_eddy_current_core_losses((eddyCurrentsVolumetricLosses + excessEddyCurrentsVolumetricLosses) * effectiveVolume);
    result.set_hysteresis_core_losses(hysteresisVolumetricLosses * effectiveVolume);
    result.set_magnetic_flux_density(magneticFluxDensity);
    result.set_method_used("Roshen");
    result.set_origin(ResultOrigin::SIMULATION);
    result.set_temperature(temperature);
    result.set_volumetric_losses(volumetricLosses);

    return result;
} 
 
double CoreLossesRoshenModel::get_core_volumetric_losses(CoreMaterial coreMaterial,
                                                         OperatingPointExcitation excitation,
                                                         double temperature) {
    Core ringCore;
    ringCore.set_name("Dummy Ring Core");
    ringCore.get_mutable_functional_description().set_material(coreMaterial);
    ringCore.get_mutable_functional_description().set_shape("T 10/6/4");
    ringCore.get_mutable_functional_description().set_number_stacks(1);
    ringCore.get_mutable_functional_description().set_type(CoreType::TOROIDAL);
    ringCore.get_mutable_functional_description().set_gapping({});
    ringCore.process_data();

    double volumetricLosses = 0;
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    auto parameters = get_roshen_parameters(ringCore, excitation, temperature);
    double hysteresisVolumetricLosses = get_hysteresis_losses_density(parameters, excitation);
    double eddyCurrentsVolumetricLosses = get_eddy_current_losses_density(ringCore, excitation, parameters["resistivity"]);
    double excessEddyCurrentsVolumetricLosses = 0;
    if (parameters.count("excessLossesCoefficient")) {
        excessEddyCurrentsVolumetricLosses = get_excess_eddy_current_losses_density(
            excitation, parameters["resistivity"], parameters["excessLossesCoefficient"]);
    }
    volumetricLosses = hysteresisVolumetricLosses + eddyCurrentsVolumetricLosses + excessEddyCurrentsVolumetricLosses;
    return volumetricLosses;
} 
 
std::map<std::string, double> CoreLossesRoshenModel::get_roshen_parameters(Core core,
                                                                           OperatingPointExcitation excitation,
                                                                           double temperature) {
    std::map<std::string, double> roshenParameters;
    auto materialData =  core.resolve_material();

    auto roshenData = CoreLossesModel::get_method_data(materialData, "roshen");

    roshenParameters["coerciveForce"] = core.get_coercive_force(temperature);
    roshenParameters["remanence"] = core.get_remanence(temperature);
    roshenParameters["saturationMagneticFluxDensity"] = core.get_magnetic_flux_density_saturation(temperature, false);
    roshenParameters["saturationMagneticFieldStrength"] = core.get_magnetic_field_strength_saturation(temperature);

    if (roshenData.get_coefficients()) {
        auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
        double frequency = excitation.get_frequency();
        magneticFluxDensity = Inputs::standardize_waveform(magneticFluxDensity, frequency);
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
        auto resistivityModel = ResistivityModel::factory(ResistivityModels::CORE_MATERIAL);
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

    auto calculate_magnetic_flux_density_waveform = [&](std::vector<double> magneticFieldStrength_waveform,
                                                  bool loop_is_upper = true) {
        std::vector<double> magneticFluxDensityWaveform;
        for (auto& magneticFieldStrength : magneticFieldStrength_waveform) {
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

    // _hysteresisMajorH = magneticFieldStrengthPoints;
    // _hysteresisMajorLoopTop = upperMagneticFluxDensityWaveform;
    // _hysteresisMajorLoopBottom = lowerMagneticFluxDensityWaveform;
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

    // _hysteresisMinorLoopTop = cutUpperMagneticFluxDensityWaveform;
    // _hysteresisMinorLoopBottom = cutLowerMagneticFluxDensityWaveform;

    size_t minimum_length =
        std::min(cutUpperMagneticFluxDensityWaveform.size(), cutLowerMagneticFluxDensityWaveform.size());
    double bhArea = 0;
    for (size_t i = 0; i < minimum_length; i++) {
        bhArea += fabs(cutUpperMagneticFluxDensityWaveform[i] - cutLowerMagneticFluxDensityWaveform[i]) *
                  constants.roshenMagneticFieldStrengthStep;
    }

    double hysteresisLossesDensity = bhArea * frequency;

    if (bhArea < 0) {
        throw std::runtime_error("Negative hysteresis losses");
    }

    return hysteresisLossesDensity;
}

double CoreLossesRoshenModel::get_eddy_current_losses_density(Core core,
                                                              OperatingPointExcitation excitation,
                                                              double resistivity) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    magneticFluxDensity = Inputs::standardize_waveform(magneticFluxDensity, frequency);
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value().get_data();
    std::vector<double> magneticFluxDensityTime;
    if (magneticFluxDensity.get_waveform().value().get_time()) {
        magneticFluxDensityTime = magneticFluxDensity.get_waveform().value().get_time().value();
    }
    if (!core.get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }

    double centralColumnArea = core.get_processed_description().value().get_columns()[0].get_area();

    double volumetricLossesIntegration = 0;
    double timeDifference;

    for (size_t i = 0; i < magneticFluxDensityWaveform.size() - 1; ++i) {
        if (magneticFluxDensity.get_waveform().value().get_time()) {
            timeDifference = magneticFluxDensityTime[i + 1] - magneticFluxDensityTime[i];
        }
        else {
            timeDifference = 1 / frequency / settings->get_inputs_number_points_sampled_waveforms();
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
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = excitation.get_frequency();
    magneticFluxDensity = Inputs::standardize_waveform(magneticFluxDensity, frequency);
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
            timeDifference = 1 / frequency / settings->get_inputs_number_points_sampled_waveforms();
        }
        volumetricLossesIntegration +=
            pow(fabs(magneticFluxDensityWaveform[i + 1] - magneticFluxDensityWaveform[i]) / timeDifference, 3 / 2) *
            timeDifference;
    }

    double excessEddyCurrentLossesDensity = sqrt(alphaTimesN0 / resistivity) * frequency * volumetricLossesIntegration;

    return excessEddyCurrentLossesDensity;
}

double CoreLossesProprietaryModel::get_core_volumetric_losses(CoreMaterial coreMaterial,
                                                             OperatingPointExcitation excitation,
                                                             double temperature) {

    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = Inputs::get_switching_frequency(excitation);
    double mainHarmonicMagneticFluxDensityPeak = magneticFluxDensity.get_processed().value().get_peak().value();

    double magneticFluxDensityAcPeak = Inputs::get_magnetic_flux_density_peak(excitation, frequency) - magneticFluxDensity.get_processed().value().get_offset();

    CoreLossesOutput result;
    double volumetricLosses = -1;

    if (coreMaterial.get_manufacturer_info().get_name() == "Micrometals") {
        auto micrometalsData = CoreLossesModel::get_method_data(coreMaterial, "micrometals");
        double a = micrometalsData.get_a().value();
        double b = micrometalsData.get_b().value();
        double c = micrometalsData.get_c().value();
        double d = micrometalsData.get_d().value();
        volumetricLosses = frequency / (a / pow(magneticFluxDensityAcPeak, 3) + b / pow(magneticFluxDensityAcPeak, 2.3) + c / pow(magneticFluxDensityAcPeak, 1.65)) + d * pow(magneticFluxDensityAcPeak, 2) * pow(frequency, 2);
    }

    else if (coreMaterial.get_manufacturer_info().get_name() == "Magnetics") {
        auto magneticsData = CoreLossesModel::get_method_data(coreMaterial, "magnetics");
        double a = magneticsData.get_a().value();
        double b = magneticsData.get_b().value();
        double c = magneticsData.get_c().value();
        if (b > 2) {
            volumetricLosses = a * pow(mainHarmonicMagneticFluxDensityPeak, b - 2) * pow(frequency, c) * pow(magneticFluxDensityAcPeak, 2);
        }
        else {
            volumetricLosses = a * pow(magneticFluxDensityAcPeak, b) * pow(frequency, c);
        }
    }
    else if (coreMaterial.get_manufacturer_info().get_name() == "Poco") {
        auto magneticsData = CoreLossesModel::get_method_data(coreMaterial, "poco");
        double a = magneticsData.get_a().value();
        double b = magneticsData.get_b().value();
        double c = magneticsData.get_c().value();
        volumetricLosses = 1000 * (a * pow(magneticFluxDensityAcPeak * 10, b) * frequency / 1000 + c * pow(magneticFluxDensityAcPeak * 10 * frequency / 1000, 2));
    }

    else if (coreMaterial.get_manufacturer_info().get_name() == "TDG") {
        auto magneticsData = CoreLossesModel::get_method_data(coreMaterial, "tdg");
        double a = magneticsData.get_a().value();
        double b = magneticsData.get_b().value();
        double c = magneticsData.get_c().value();
        double d = magneticsData.get_d().value();
        volumetricLosses = 1000 * pow(magneticFluxDensityAcPeak * 10, a) * (b * frequency / 1000 + c * pow(frequency / 1000, d));
    }
    else {
        // throw std::invalid_argument("No volumetric losses method for manufacturer: " + coreMaterial.get_manufacturer_info().get_name());
    }

    return volumetricLosses;
}

std::map<std::string, std::string> CoreLossesProprietaryModel::get_core_volumetric_losses_equations(CoreMaterial coreMaterial) {
    std::map<std::string, std::string> equations;
    if (coreMaterial.get_manufacturer_info().get_name() == "Micrometals") {
        equations["volumetricCoreLosses"] = "f / (a / B^3 + b / B^2.3 + c / B^1.65) + d * B^2 * f^2";
    }
    else if (coreMaterial.get_manufacturer_info().get_name() == "Magnetics") {
        equations["volumetricCoreLosses"] = "a * B^b * f^c";
    }
    else if (coreMaterial.get_manufacturer_info().get_name() == "Poco") {
        equations["volumetricCoreLosses"] = "1000 * (a * f / 1000 * (B * 10)^b + c * (B * 10 * f / 1000)^2)";
    }
    else if (coreMaterial.get_manufacturer_info().get_name() == "TDG") {
        equations["volumetricCoreLosses"] = "1000 * ((B * 10)^a) * (b * f / 1000 + c * (f / 1000)^d)";
    }
    else {
        throw std::invalid_argument("No volumetric losses method for manufacturer: " + coreMaterial.get_manufacturer_info().get_name());
    }

    return equations;
}

std::map<std::string, std::string> CoreLossesProprietaryModel::get_core_volumetric_losses_equations(CoreLossesMethodData coreLossesMethodData) {
    std::map<std::string, std::string> equations;
    if (coreLossesMethodData.get_method() == VolumetricCoreLossesMethodType::MICROMETALS) {
        equations["volumetricCoreLosses"] = "f / (a / B^3 + b / B^2.3 + c / B^1.65) + d * B^2 * f^2";
    }
    else if (coreLossesMethodData.get_method() == VolumetricCoreLossesMethodType::MAGNETICS) {
        equations["volumetricCoreLosses"] = "a * B^b * f^c";
    }
    else if (coreLossesMethodData.get_method() == VolumetricCoreLossesMethodType::POCO) {
        equations["volumetricCoreLosses"] = "1000 * (a * f / 1000 * (B * 10)^b + c * (B * 10 * f / 1000)^2)";
    }
    else if (coreLossesMethodData.get_method() == VolumetricCoreLossesMethodType::TDG) {
        equations["volumetricCoreLosses"] = "1000 * ((B * 10)^a) * (b * f / 1000 + c * (f / 1000)^d)";
    }
    else {
        throw std::invalid_argument("No volumetric losses method for method: " + std::string{magic_enum::enum_name(coreLossesMethodData.get_method())});
    }

    return equations;
}

double CoreLossesProprietaryModel::get_core_mass_losses(CoreMaterial coreMaterial,
                                                             OperatingPointExcitation excitation,
                                                             double temperature) {

    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double frequency = Inputs::get_switching_frequency(excitation);
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() - magneticFluxDensity.get_processed().value().get_offset();

    CoreLossesOutput result;
    double massLosses = -1;

    if (coreMaterial.get_manufacturer_info().get_name() == "Magnetec") {
        massLosses = 80 * pow(frequency / 100000, 1.8) * pow(magneticFluxDensityAcPeak * 2 / 0.3, 2);
    }
    return massLosses;
}
 
double CoreLossesSteinmetzModel::get_frequency_from_core_losses(Core core,
                                                                SignalDescriptor magneticFluxDensity,
                                                                double temperature,
                                                                double coreLosses) {
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value() - magneticFluxDensity.get_processed().value().get_offset();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();

    SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    double frequency = 100000;

    steinmetzDatum = get_steinmetz_coefficients(core.resolve_material(), frequency);
    double alpha = steinmetzDatum.get_alpha();
    double convergeAlpha = 0;

    do {
        double k = steinmetzDatum.get_k();
        alpha = steinmetzDatum.get_alpha();
        double beta = steinmetzDatum.get_beta();
        double volumetricLosses = coreLosses / effectiveVolume / CoreLossesModel::apply_temperature_coefficients(1, steinmetzDatum, temperature);

        frequency = pow(volumetricLosses / (k * pow(magneticFluxDensityAcPeak, beta)), 1 / alpha);
        steinmetzDatum = get_steinmetz_coefficients(core.resolve_material(), frequency);
        convergeAlpha = steinmetzDatum.get_alpha();
    }
    while (convergeAlpha != alpha);

    return frequency;
}

SignalDescriptor CoreLossesSteinmetzModel::get_magnetic_flux_density_from_core_losses(Core core,
                                                                                              double frequency,
                                                                                              double temperature,
                                                                                              double coreLosses) {
    SignalDescriptor magneticFluxDensity;
    Processed processed;
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(0);
    magneticFluxDensity.set_processed(processed);

    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();

    SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;

    steinmetzDatum = get_steinmetz_coefficients(core.resolve_material(), frequency);

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

double CoreLossesLossFactorModel::calculate_magnetizing_inductance_from_excitation(Core core, OperatingPointExcitation excitation, double temperature) {
    auto currentPeak = excitation.get_magnetizing_current()->get_processed()->get_peak().value();
    auto magneticFluxDensityPeak = excitation.get_magnetic_flux_density()->get_processed()->get_peak().value();
    auto coreMaterial = core.resolve_material();
    double effectiveArea = core.get_processed_description().value().get_effective_parameters().get_effective_area();

    double initialPermeability = InitialPermeability::get_initial_permeability(coreMaterial, temperature);
    auto reluctanceModel = ReluctanceModel::factory();
    auto reluctance = reluctanceModel->get_core_reluctance(core, initialPermeability).get_core_reluctance();
    int64_t numberTurns = round(ceilFloat(magneticFluxDensityPeak / currentPeak * reluctance * effectiveArea, 0));

    double magnetizingInductance = numberTurns * numberTurns / reluctance;
    return magnetizingInductance;
}

CoreLossesOutput CoreLossesLossFactorModel::get_core_losses(Core core,
                                                        OperatingPointExcitation excitation,
                                                        double temperature) {
    if (!excitation.get_magnetizing_current()) {
        throw std::runtime_error("Missing magnetizing current in excitation");
    }
    if (!excitation.get_magnetizing_current()->get_processed()) {
        throw std::runtime_error("Magnetizing current not processed");
    }
    if (!excitation.get_magnetizing_current()->get_processed()->get_rms()) {
        auto magnetizingCurrent = excitation.get_magnetizing_current().value();
        magnetizingCurrent.set_processed(Inputs::calculate_processed_data(excitation.get_magnetizing_current()->get_waveform().value(), excitation.get_frequency()));
        excitation.set_magnetizing_current(magnetizingCurrent);

        // throw std::runtime_error("Missing RMS value in magnetizing current");
    }
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    auto coreMaterial = core.resolve_material();
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
    double magnetizingInductance = calculate_magnetizing_inductance_from_excitation(core, excitation, temperature);

    auto volumetricLosses = get_core_volumetric_losses(coreMaterial, excitation, temperature, magnetizingInductance);

    CoreLossesOutput result;
    result.set_core_losses(volumetricLosses * effectiveVolume);
    result.set_magnetic_flux_density(magneticFluxDensity);
    result.set_method_used(_modelName);
    result.set_origin(ResultOrigin::SIMULATION);
    result.set_temperature(temperature);
    result.set_volumetric_losses(volumetricLosses);

    return result;
} 

double CoreLossesLossFactorModel::get_core_volumetric_losses(CoreMaterial coreMaterial,
                                                             OperatingPointExcitation excitation,
                                                             double temperature,
                                                             double magnetizingInductance) {
    
    if (!excitation.get_magnetizing_current()) {
        throw std::runtime_error("Missing magnetizing current in excitation");
    }
    if (!excitation.get_magnetizing_current()->get_processed()) {
        throw std::runtime_error("Magnetizing current not processed");
    }
    if (!excitation.get_magnetizing_current()->get_processed()->get_rms()) {
        auto magnetizingCurrent = excitation.get_magnetizing_current().value();
        magnetizingCurrent.set_processed(Inputs::calculate_processed_data(excitation.get_magnetizing_current()->get_waveform().value(), excitation.get_frequency()));
        excitation.set_magnetizing_current(magnetizingCurrent);
        // throw std::runtime_error("Missing RMS value in magnetizing current");
    }
    auto currentRms = excitation.get_magnetizing_current()->get_processed()->get_rms().value();
    double frequency = Inputs::get_switching_frequency(excitation);

    auto seriesResistance = get_core_losses_series_resistance(coreMaterial, frequency, temperature, magnetizingInductance);
    double volumetricLosses = seriesResistance * pow(currentRms, 2);

    return volumetricLosses;
}

double CoreLossesLossFactorModel::get_core_losses_series_resistance(CoreMaterial coreMaterial,
                                                                    double frequency,
                                                                    double temperature,
                                                                    double magnetizingInductance) {

    CoreLossesOutput result;
    double initialPermeability = InitialPermeability::get_initial_permeability(coreMaterial, temperature);

    if (!lossFactorInterps.contains(coreMaterial.get_name())) {

        auto lossFactorData = CoreLossesModel::get_method_data(coreMaterial, "loss_factor");
        auto lossFactorPoints = lossFactorData.get_factors().value();

        int n = lossFactorPoints.size();
        std::vector<double> x, y;


        for (int i = 0; i < n; i++) {
            if (x.size() == 0 || (*lossFactorPoints[i].get_frequency()) != x.back()) {
                x.push_back(*lossFactorPoints[i].get_frequency());
                y.push_back(lossFactorPoints[i].get_value());
            }
        }

        tk::spline interp(x, y, tk::spline::cspline_hermite);
        lossFactorInterps[coreMaterial.get_name()] = interp;
    }
    double lossFactorValue = lossFactorInterps[coreMaterial.get_name()](frequency);

    auto lossTangent = lossFactorValue * initialPermeability;
    auto seriesResistance = lossTangent * 2 * std::numbers::pi * frequency * magnetizingInductance;

    return seriesResistance;
}

double CoreLossesProprietaryModel::get_frequency_from_core_losses(Core core,
                                                                    SignalDescriptor magneticFluxDensity,
                                                                    [[maybe_unused]]double temperature,
                                                                    double coreLosses) {

    double frequency = -1;
    double magneticFluxDensityAcPeak = magneticFluxDensity.get_processed().value().get_peak().value();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    auto materialData = core.resolve_material();
    double volumetricLosses = coreLosses / effectiveVolume;

    if (materialData.get_manufacturer_info().get_name() == "Micrometals") {
        auto micrometalsData = CoreLossesModel::get_method_data(materialData, "micrometals");
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
        auto magneticsData = CoreLossesModel::get_method_data(materialData, "magnetics");
        double a = magneticsData.get_a().value();
        double b = magneticsData.get_b().value();
        double c = magneticsData.get_c().value();
        frequency = pow(volumetricLosses / (a * pow(magneticFluxDensityAcPeak, b)), 1. / c);
    }

    if (materialData.get_manufacturer_info().get_name() == "Poco") {
        auto magneticsData = CoreLossesModel::get_method_data(materialData, "poco");
        double a = magneticsData.get_a().value();
        double b = magneticsData.get_b().value();
        double c = magneticsData.get_c().value();
        double auxA = c * pow(magneticFluxDensityAcPeak * 10, 2);
        double auxB = a * pow(magneticFluxDensityAcPeak * 10, b);
        double auxC = -volumetricLosses / 1000;
        frequency = 1000 * (-auxB + sqrt(pow(auxB, 2) - 4 * auxA * auxC)) / (2 * auxA);
    }

    if (materialData.get_manufacturer_info().get_name() == "TDG") {
        throw std::runtime_error("Not implemented fot TDG");
    }

    if (materialData.get_manufacturer_info().get_name() == "Magnetec") {
        double mass = core.get_mass();
        double massLosses = coreLosses / mass;
        frequency = pow(massLosses / 80 / pow(magneticFluxDensityAcPeak * 2 / 0.3, 2), 1.0 / 1.8) * 100000;
    }

    return frequency;
}

double CoreLossesModel::_get_frequency_from_core_losses(Core core,
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
        double error = fabs(coreLossesCalculated.get_core_losses() - coreLosses) / coreLosses;
        if (error < minimumError) {
            minimumError = error;
            frequencyMinimumError = frequency;
        }
    }
    return frequencyMinimumError;
}

SignalDescriptor CoreLossesModel::_get_magnetic_flux_density_from_core_losses(Core core,
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
        double error = fabs(coreLossesCalculated.get_core_losses() - coreLosses) / coreLosses;
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
