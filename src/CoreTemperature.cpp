#include "CoreTemperature.h"

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
        return std::make_shared<CoreTemperatureKazimierczukModel>();
    }
    else if (modelName == CoreTemperatureModels::MANIKTALA) {
        return std::make_shared<CoreTemperatureManiktalaModel>();
    }
    else if (modelName == CoreTemperatureModels::TDK) {
        return std::make_shared<CoreTemperatureTdkModel>();
    }
    else if (modelName == CoreTemperatureModels::DIXON) {
        return std::make_shared<CoreTemperatureDixonModel>();
    }
    else if (modelName == CoreTemperatureModels::AMIDON) {
        return std::make_shared<CoreTemperatureAmidonModel>();
    }

    else
        throw std::runtime_error(
            "Unknown Core losses mode, available options are: {KAZIMIERCZUK, MANIKTALA, TDK, DIXON, AMIDON}");
}

TemperatureOutput CoreTemperatureManiktalaModel::get_core_temperature(CoreWrapper core,
                                                                      double coreLosses,
                                                                      double ambientTemperature) {
    double cubeVolume = core.get_processed_description().value().get_depth() * core.get_processed_description().value().get_width() * core.get_processed_description().value().get_height();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    double thermalResistance = 53 * pow(effectiveVolume * 1000000, -0.54);
    double temperatureRise = coreLosses * thermalResistance;
    double maximumTemperature = ambientTemperature + temperatureRise;

    TemperatureOutput result;
    result.set_bulk_thermal_resistance((maximumTemperature - ambientTemperature) / cubeVolume);
    result.set_initial_temperature(ambientTemperature);
    result.set_maximum_temperature(maximumTemperature);
    result.set_method_used("Maniktala");
    result.set_origin(ResultOrigin::SIMULATION);
    return result;
}

TemperatureOutput CoreTemperatureKazimierczukModel::get_core_temperature(CoreWrapper core,
                                                                         double coreLosses,
                                                                         double ambientTemperature) {
    double width = core.get_processed_description().value().get_width();
    double height = core.get_processed_description().value().get_height();
    double depth = core.get_processed_description().value().get_depth();
    double cubeVolume = depth * width * height;
    double coreSurface = 2 * depth * height + 2 * height * width + 2 * width * depth;
    double temperatureRise = pow(0.1 * coreLosses / coreSurface, 0.826);
    double maximumTemperature = ambientTemperature + temperatureRise;

    TemperatureOutput result;
    result.set_bulk_thermal_resistance((maximumTemperature - ambientTemperature) / cubeVolume);
    result.set_initial_temperature(ambientTemperature);
    result.set_maximum_temperature(maximumTemperature);
    result.set_method_used("Kazimierczuk");
    result.set_origin(ResultOrigin::SIMULATION);

    return result;
}

TemperatureOutput CoreTemperatureTdkModel::get_core_temperature(CoreWrapper core,
                                                                double coreLosses,
                                                                double ambientTemperature) {
    double cubeVolume = core.get_processed_description().value().get_depth() * core.get_processed_description().value().get_width() * core.get_processed_description().value().get_height();
    double effectiveVolume = core.get_processed_description().value().get_effective_parameters().get_effective_volume();
    double thermalResistance = 1 / sqrt(effectiveVolume * 1000000);
    double temperatureRise = coreLosses * thermalResistance;
    double maximumTemperature = ambientTemperature + temperatureRise;

    TemperatureOutput result;
    result.set_bulk_thermal_resistance((maximumTemperature - ambientTemperature) / cubeVolume);
    result.set_initial_temperature(ambientTemperature);
    result.set_maximum_temperature(maximumTemperature);
    result.set_method_used("Tdk");
    result.set_origin(ResultOrigin::SIMULATION);

    return result;
}

TemperatureOutput CoreTemperatureDixonModel::get_core_temperature(CoreWrapper core,
                                                                  double coreLosses,
                                                                  double ambientTemperature) {
    double cubeVolume = core.get_processed_description().value().get_depth() * core.get_processed_description().value().get_width() * core.get_processed_description().value().get_height();
    double centralColumnArea = core.get_processed_description().value().get_columns()[0].get_area();
    double windingWindowArea = core.get_processed_description().value().get_winding_windows()[0].get_area().value();
    double areaProduct = centralColumnArea * windingWindowArea * 100000000;
    double thermalResistance = 23 * pow(areaProduct, -0.37);
    double temperatureRise = coreLosses * thermalResistance;
    double maximumTemperature = ambientTemperature + temperatureRise;

    TemperatureOutput result;
    result.set_bulk_thermal_resistance((maximumTemperature - ambientTemperature) / cubeVolume);
    result.set_initial_temperature(ambientTemperature);
    result.set_maximum_temperature(maximumTemperature);
    result.set_method_used("Dixon");
    result.set_origin(ResultOrigin::SIMULATION);

    return result;
}

TemperatureOutput CoreTemperatureAmidonModel::get_core_temperature(CoreWrapper core,
                                                                   double coreLosses,
                                                                   double ambientTemperature) {
    double width = core.get_processed_description().value().get_width();
    double height = core.get_processed_description().value().get_height();
    double depth = core.get_processed_description().value().get_depth();
    double cubeVolume = depth * width * height;
    double coreSurface = 2 * depth * height + 2 * height * width + 2 * width * depth;
    double temperatureRise = pow(coreLosses * 1000 / (coreSurface * 10000), 0.833);
    double maximumTemperature = ambientTemperature + temperatureRise;

    TemperatureOutput result;
    result.set_bulk_thermal_resistance((maximumTemperature - ambientTemperature) / cubeVolume);
    result.set_initial_temperature(ambientTemperature);
    result.set_maximum_temperature(maximumTemperature);
    result.set_method_used("Amidon");
    result.set_origin(ResultOrigin::SIMULATION);

    return result;
}

} // namespace OpenMagnetics
