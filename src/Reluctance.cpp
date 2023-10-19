#include "Reluctance.h"
#include "Defaults.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {


double ReluctanceModel::get_ungapped_core_reluctance(CoreWrapper core, double initialPermeability) {
    auto constants = Constants();
    double absolutePermeability = constants.vacuumPermeability * initialPermeability;
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();

    double reluctanceCore = effectiveLength / (absolutePermeability * effectiveArea);
    return reluctanceCore;
}

double ReluctanceModel::get_ungapped_core_reluctance(CoreWrapper core, OperatingPoint* operatingPoint) {
    OpenMagnetics::InitialPermeability initialPermeability;

    auto coreMaterial = core.get_functional_description().get_material();

    double initialPermeabilityValue;
    if (operatingPoint != nullptr) {
        double temperature =
            operatingPoint->get_conditions().get_ambient_temperature(); // TODO: Use a future calculated temperature
        _magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(temperature, true);
        auto frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();
        initialPermeabilityValue = initialPermeability.get_initial_permeability(
            coreMaterial, &temperature, nullptr, &frequency);
    }
    else {
        initialPermeabilityValue =
            initialPermeability.get_initial_permeability(coreMaterial);
        _magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(true);
    }
    return get_ungapped_core_reluctance(core, initialPermeabilityValue);
}

std::map<std::string, double> ReluctanceZhangModel::get_gap_reluctance(CoreGap gapInfo) {
    double perimeter = 0;
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw std::runtime_error("Gap Distance Closest Normal Surface is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gap_shape = *(gapInfo.get_shape());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    auto reluctance_internal = gapLength / (constants.vacuumPermeability * gapArea);
    double reluctance_fringing = 0;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];
    auto gapSectionDepth = gapSectionDimensions[1];

    if (gap_shape == ColumnShape::ROUND) {
        perimeter = std::numbers::pi * gapSectionWidth;
    }
    else { // TODO: Properly calculate perimeter for all shapes
        perimeter = gapSectionWidth * 2 + gapSectionDepth * 2;
    }

    if (gapLength > 0) {
        reluctance_fringing = std::numbers::pi / (constants.vacuumPermeability * perimeter *
                                                  log((2 * distanceClosestNormalSurface + gapLength) / gapLength));
    }

    if (std::isnan(reluctance_internal) || reluctance_internal == 0) {
        throw std::runtime_error("reluctance_internal cannot be 0 or NaN");
    }

    if (std::isnan(reluctance_fringing) || reluctance_fringing == 0) {
        throw std::runtime_error("reluctance_fringing cannot be 0 or NaN");
    }

    double reluctance = 1. / (1. / reluctance_internal + 1. / reluctance_fringing);

    if (gapLength > 0) {
        fringingFactor = gapLength / (constants.vacuumPermeability * gapArea * reluctance);
    }
    std::map<std::string, double> result;
    result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringingFactor);
    result["reluctance"] = reluctance;
    if (reluctance > 0) {
        result["permeance"] = 1 / reluctance;
    }
    else {
        result["permeance"] = std::numeric_limits<double>::infinity();
    }
    result["fringing_factor"] = fringingFactor;

    return result;
};

std::map<std::string, double> ReluctanceMuehlethalerModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw std::runtime_error("Gap Distance Closest Normal Surface is not set");
    }
    auto gap_shape = *(gapInfo.get_shape());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    double reluctance;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];
    auto gapSectionDepth = gapSectionDimensions[1];

    if (gap_shape == ColumnShape::ROUND) {
        double gamma_r = get_reluctance_type_1(gapLength / 2, gapSectionWidth / 2, distanceClosestNormalSurface) /
                         (gapLength / constants.vacuumPermeability / (gapSectionWidth / 2));
        reluctance = pow(gamma_r, 2) * gapLength /
                     (constants.vacuumPermeability * std::numbers::pi * pow(gapSectionWidth / 2, 2));
        fringingFactor = 1 / gamma_r;
    }
    else {
        double gamma_x = get_reluctance_type_1(gapLength / 2, gapSectionWidth, distanceClosestNormalSurface) /
                         (gapLength / constants.vacuumPermeability / gapSectionWidth);
        double gamma_y = get_reluctance_type_1(gapLength / 2, gapSectionDepth, distanceClosestNormalSurface) /
                         (gapLength / constants.vacuumPermeability / gapSectionDepth);
        double gamma = gamma_x * gamma_y;
        reluctance = gamma * gapLength / (constants.vacuumPermeability * gapSectionDepth * gapSectionWidth);
        fringingFactor = 1 / gamma;
    }

    std::map<std::string, double> result;
    result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringingFactor);
    result["reluctance"] = reluctance;
    if (reluctance > 0) {
        result["permeance"] = 1 / reluctance;
    }
    else {
        result["permeance"] = std::numeric_limits<double>::infinity();
    }
    result["fringing_factor"] = fringingFactor;

    return result;
};

std::map<std::string, double> ReluctanceEffectiveAreaModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gap_shape = *(gapInfo.get_shape());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    double reluctance;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];
    auto gapSectionDepth = gapSectionDimensions[1];

    if (gapLength > 0) {
        if (gap_shape == ColumnShape::ROUND) {
            fringingFactor = pow(1 + gapLength / gapSectionWidth, 2);
        }
        else {
            fringingFactor = (gapSectionDepth + gapLength) * (gapSectionWidth + gapLength) /
                              (gapSectionDepth * gapSectionWidth);
        }
    }

    reluctance = gapLength / (constants.vacuumPermeability * gapArea * fringingFactor);

    std::map<std::string, double> result;
    result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringingFactor);
    result["reluctance"] = reluctance;
    if (reluctance > 0) {
        result["permeance"] = 1 / reluctance;
    }
    else {
        result["permeance"] = std::numeric_limits<double>::infinity();
    }
    result["fringing_factor"] = fringingFactor;

    return result;
};

std::map<std::string, double> ReluctanceEffectiveLengthModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gap_shape = *(gapInfo.get_shape());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    double reluctance;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];
    auto gapSectionDepth = gapSectionDimensions[1];

    if (gapLength > 0) {
        if (gap_shape == ColumnShape::ROUND) {
            fringingFactor = pow(1 + gapLength / gapSectionWidth, 2);
        }
        else {
            fringingFactor = (1 + gapLength / gapSectionDepth) * (1 + gapLength / gapSectionWidth);
        }
    }

    reluctance = gapLength / (constants.vacuumPermeability * gapArea * fringingFactor);

    std::map<std::string, double> result;
    result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringingFactor);
    result["reluctance"] = reluctance;
    if (reluctance > 0) {
        result["permeance"] = 1 / reluctance;
    }
    else {
        result["permeance"] = std::numeric_limits<double>::infinity();
    }
    result["fringing_factor"] = fringingFactor;

    return result;
};

std::map<std::string, double> ReluctancePartridgeModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw std::runtime_error("Gap Distance Closest Normal Surface is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    double reluctance;
    double fringingFactor = 1;

    if (gapLength > 0) {
        fringingFactor = 1 + gapLength / sqrt(gapArea) * log(2 * 2 * distanceClosestNormalSurface / gapLength);
    }

    reluctance = gapLength / (constants.vacuumPermeability * gapArea * fringingFactor);

    std::map<std::string, double> result;
    result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringingFactor);
    result["reluctance"] = reluctance;
    if (reluctance > 0) {
        result["permeance"] = 1 / reluctance;
    }
    else {
        result["permeance"] = std::numeric_limits<double>::infinity();
    }
    result["fringing_factor"] = fringingFactor;

    return result;
};

std::map<std::string, double> ReluctanceStengleinModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw std::runtime_error("Gap Distance Closest Normal Surface is not set");
    }
    if (!gapInfo.get_coordinates()) {
        throw std::runtime_error("Gap Corrdinates are not set");
    }
    if (!gapInfo.get_distance_closest_parallel_surface()) {
        throw std::runtime_error("Gap Distance Closest Parallel Surface is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gapCoordinates = *(gapInfo.get_coordinates());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    auto distanceClosestParallelSurface = *(gapInfo.get_distance_closest_parallel_surface());
    double reluctance;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];

    if (gapLength > 0) {
        double c = gapSectionWidth / 2 + distanceClosestParallelSurface;
        double b = gapSectionWidth / 2 + 0.001;
        double l1 = distanceClosestNormalSurface * 2;
        double lg = gapLength;
        double rc = gapSectionWidth / 2;
        double rx = gapSectionWidth / 2;
        double aux1 = 1 + 2. / sqrt(std::numbers::pi) * lg / (2 * rc) * log(2.1 * rx / lg);
        double aux2 = 1. / 6. * (pow(c, 2) + 2 * c * b + pow(b, 2)) / pow(b, 2);

        double gamma = aux1 + (aux2 - aux1) * pow(lg / l1, 2 * std::numbers::pi);

        fringingFactor = alpha(rx, l1, lg) * pow(gapCoordinates[1] / l1, 2) + gamma;
    }

    reluctance = gapLength / (constants.vacuumPermeability * gapArea * fringingFactor);

    std::map<std::string, double> result;
    result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringingFactor);
    result["reluctance"] = reluctance;
    if (reluctance > 0) {
        result["permeance"] = 1 / reluctance;
    }
    else {
        result["permeance"] = std::numeric_limits<double>::infinity();
    }
    result["fringing_factor"] = fringingFactor;

    return result;
};

std::map<std::string, double> ReluctanceClassicModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    double reluctance;
    double fringingFactor = 1;

    reluctance = gapLength / (constants.vacuumPermeability * gapArea);

    std::map<std::string, double> result;
    result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringingFactor);
    result["reluctance"] = reluctance;
    if (reluctance > 0) {
        result["permeance"] = 1 / reluctance;
    }
    else {
        result["permeance"] = std::numeric_limits<double>::infinity();
    }
    result["fringing_factor"] = fringingFactor;

    return result;
};

std::map<std::string, double> ReluctanceBalakrishnanModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    double reluctance;
    double fringingFactor = 1;
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto gapSectionDepth = gapSectionDimensions[1];

    reluctance = 1. / (constants.vacuumPermeability *
                       (gapArea / gapLength +
                        2. * gapSectionDepth / std::numbers::pi *
                            (1 + log(std::numbers::pi * distanceClosestNormalSurface / (2 * gapLength)))));

    if (gapLength > 0) {
        fringingFactor = gapLength / (constants.vacuumPermeability * gapArea * reluctance);
    }

    std::map<std::string, double> result;
    result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringingFactor);
    result["reluctance"] = reluctance;
    if (reluctance > 0) {
        result["permeance"] = 1 / reluctance;
    }
    else {
        result["permeance"] = std::numeric_limits<double>::infinity();
    }
    result["fringing_factor"] = fringingFactor;

    return result;
};

std::shared_ptr<ReluctanceModel> ReluctanceModel::factory(ReluctanceModels modelName) {
    if (modelName == ReluctanceModels::ZHANG) {
        return std::make_shared<ReluctanceZhangModel>();
    }
    else if (modelName == ReluctanceModels::PARTRIDGE) {
        return std::make_shared<ReluctancePartridgeModel>();
    }
    else if (modelName == ReluctanceModels::EFFECTIVE_AREA) {
        return std::make_shared<ReluctanceEffectiveAreaModel>();
    }
    else if (modelName == ReluctanceModels::EFFECTIVE_LENGTH) {
        return std::make_shared<ReluctanceEffectiveLengthModel>();
    }
    else if (modelName == ReluctanceModels::MUEHLETHALER) {
        return std::make_shared<ReluctanceMuehlethalerModel>();
    }
    else if (modelName == ReluctanceModels::STENGLEIN) {
        return std::make_shared<ReluctanceStengleinModel>();
    }
    else if (modelName == ReluctanceModels::BALAKRISHNAN) {
        return std::make_shared<ReluctanceBalakrishnanModel>();
    }
    else if (modelName == ReluctanceModels::CLASSIC) {
        return std::make_shared<ReluctanceClassicModel>();
    }

    else
        throw std::runtime_error("Unknown Reluctance mode, available options are: {ZHANG, PARTRIDGE, EFFECTIVE_AREA, "
                                 "EFFECTIVE_LENGTH, MUEHLETHALER, STENGLEIN, BALAKRISHNAN, CLASSIC}");
}
std::shared_ptr<ReluctanceModel> ReluctanceModel::factory() {
    auto defaults = Defaults();
    return factory(defaults.reluctanceModelDefault);
}
} // namespace OpenMagnetics
