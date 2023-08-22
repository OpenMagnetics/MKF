#include "Reluctance.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

std::map<std::string, double> ReluctanceZhangModel::get_gap_reluctance(CoreGap gapInfo) {
    double perimeter = 0;
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
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
    else { // TODO: Properly calcualte perimeter for all shapes
        perimeter = gapSectionWidth * 2 + gapSectionDepth * 2;
    }

    if (gapLength > 0) {
        reluctance_fringing = std::numbers::pi / (constants.vacuumPermeability * perimeter *
                                                  log((2 * distanceClosestNormalSurface + gapLength) / gapLength));
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
        std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctanceZhangModel);
        return reluctanceModel;
    }
    else if (modelName == ReluctanceModels::PARTRIDGE) {
        std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctancePartridgeModel);
        return reluctanceModel;
    }
    else if (modelName == ReluctanceModels::EFFECTIVE_AREA) {
        std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctanceEffectiveAreaModel);
        return reluctanceModel;
    }
    else if (modelName == ReluctanceModels::EFFECTIVE_LENGTH) {
        std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctanceEffectiveLengthModel);
        return reluctanceModel;
    }
    else if (modelName == ReluctanceModels::MUEHLETHALER) {
        std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctanceMuehlethalerModel);
        return reluctanceModel;
    }
    else if (modelName == ReluctanceModels::STENGLEIN) {
        std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctanceStengleinModel);
        return reluctanceModel;
    }
    else if (modelName == ReluctanceModels::BALAKRISHNAN) {
        std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctanceBalakrishnanModel);
        return reluctanceModel;
    }
    else if (modelName == ReluctanceModels::CLASSIC) {
        std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctanceClassicModel);
        return reluctanceModel;
    }

    else
        throw std::runtime_error("Unknown Reluctance mode, available options are: {ZHANG, PARTRIDGE, EFFECTIVE_AREA, "
                                 "EFFECTIVE_LENGTH, MUEHLETHALER, STENGLEIN, BALAKRISHNAN, CLASSIC}");
}
} // namespace OpenMagnetics
