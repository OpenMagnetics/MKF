#include "CoreTemperature.h"

#include "../tests/TestingUtils.h"
#include "Constants.h"
#include "InputsWrapper.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {
std::shared_ptr<CoreTemperatureModel> CoreTemperatureModel::factory(CoreTemperatureModels modelName) {
    if (modelName == CoreTemperatureModels::KAZIMIERCZUK) {
        std::shared_ptr<CoreTemperatureModel> CoreTemperatureModel(new CoreTemperatureKazimierczukModel);
        return CoreTemperatureModel;
    }
    else if (modelName == CoreTemperatureModels::MANIKTALA) {
        std::shared_ptr<CoreTemperatureModel> CoreTemperatureModel(new CoreTemperatureManiktalaModel);
        return CoreTemperatureModel;
    }
    else if (modelName == CoreTemperatureModels::TDK) {
        std::shared_ptr<CoreTemperatureModel> CoreTemperatureModel(new CoreTemperatureTdkModel);
        return CoreTemperatureModel;
    }
    else if (modelName == CoreTemperatureModels::DIXON) {
        std::shared_ptr<CoreTemperatureModel> CoreTemperatureModel(new CoreTemperatureDixonModel);
        return CoreTemperatureModel;
    }
    else if (modelName == CoreTemperatureModels::AMIDON) {
        std::shared_ptr<CoreTemperatureModel> CoreTemperatureModel(new CoreTemperatureAmidonModel);
        return CoreTemperatureModel;
    }

    else
        throw std::runtime_error(
            "Unknown Core losses mode, available options are: {KAZIMIERCZUK, MANIKTALA, TDK, DIXON, AMIDON}");
}

std::map<std::string, double> CoreTemperatureManiktalaModel::get_core_temperature(CoreWrapper core,
                                                                                  double coreLosses,
                                                                                  double ambientTemperature) {
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    double thermalResistance = 53 * pow(effectiveVolume * 1000000, -0.54);
    double temperatureRise = coreLosses * thermalResistance;
    double maximumTemperature = ambientTemperature + temperatureRise;

    std::map<std::string, double> result;
    result["maximumTemperature"] = maximumTemperature;

    return result;
}

std::map<std::string, double> CoreTemperatureKazimierczukModel::get_core_temperature(CoreWrapper core,
                                                                                     double coreLosses,
                                                                                     double ambientTemperature) {
    double width = core.get_processed_description().value().get_width();
    double height = core.get_processed_description().value().get_height();
    double depth = core.get_processed_description().value().get_depth();
    double coreSurface = 2 * depth * height + 2 * height * width + 2 * width * depth;
    double temperatureRise = pow(0.1 * coreLosses / coreSurface, 0.826);
    double maximumTemperature = ambientTemperature + temperatureRise;

    std::map<std::string, double> result;
    result["maximumTemperature"] = maximumTemperature;

    return result;
}

std::map<std::string, double> CoreTemperatureTdkModel::get_core_temperature(CoreWrapper core,
                                                                            double coreLosses,
                                                                            double ambientTemperature) {
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    double thermalResistance = 1 / sqrt(effectiveVolume * 1000000);
    double temperatureRise = coreLosses * thermalResistance;
    double maximumTemperature = ambientTemperature + temperatureRise;

    std::map<std::string, double> result;
    result["maximumTemperature"] = maximumTemperature;

    return result;
}

std::map<std::string, double> CoreTemperatureDixonModel::get_core_temperature(CoreWrapper core,
                                                                              double coreLosses,
                                                                              double ambientTemperature) {
    double centralColumnArea = core.get_processed_description().value().get_columns()[0].get_area();
    double windingWindowArea = core.get_processed_description().value().get_winding_windows()[0].get_area().value();
    double areaProduct = centralColumnArea * windingWindowArea * 100000000;
    double thermalResistance = 23 * pow(areaProduct, -0.37);
    double temperatureRise = coreLosses * thermalResistance;
    double maximumTemperature = ambientTemperature + temperatureRise;

    std::map<std::string, double> result;
    result["maximumTemperature"] = maximumTemperature;

    return result;
}

std::map<std::string, double> CoreTemperatureAmidonModel::get_core_temperature(CoreWrapper core,
                                                                               double coreLosses,
                                                                               double ambientTemperature) {
    double width = core.get_processed_description().value().get_width();
    double height = core.get_processed_description().value().get_height();
    double depth = core.get_processed_description().value().get_depth();
    double coreSurface = 2 * depth * height + 2 * height * width + 2 * width * depth;
    double temperatureRise = pow(coreLosses * 1000 / (coreSurface * 10000), 0.833);
    double maximumTemperature = ambientTemperature + temperatureRise;

    std::map<std::string, double> result;
    result["maximumTemperature"] = maximumTemperature;

    return result;
}

} // namespace OpenMagnetics
