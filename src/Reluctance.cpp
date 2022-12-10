#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <streambuf>
#include <magic_enum.hpp>
#include "Reluctance.h"


namespace OpenMagnetics {

    std::map<std::string, double> ReluctanceZhangModel::get_gap_reluctance(CoreGap gapInfo)
    {
        double perimeter = 0;
        auto constants = Constants();
        auto gap_length = gapInfo.get_length();
        auto gap_area = *(gapInfo.get_area());
        auto gap_shape = *(gapInfo.get_shape());
        auto gap_section_dimensions = *(gapInfo.get_section_dimensions());
        auto distance_closest_normal_surface = *(gapInfo.get_distance_closest_normal_surface());
        auto reluctance_internal = gap_length / (constants.vacuum_permeability * gap_area);
        double reluctance_fringing = 0;
        double fringing_factor = 1;
        auto gap_section_width = gap_section_dimensions[0];
        auto gap_section_depth = gap_section_dimensions[1];

        if (gap_shape == ShapeEnum::ROUND) {
            perimeter = std::numbers::pi * gap_section_width;
        }
        else { // TODO: Properly calcualte perimeter for all shapes
            perimeter = gap_section_width * 2 + gap_section_depth * 2;
        }

        if (gap_length > 0) {
            reluctance_fringing = std::numbers::pi / (constants.vacuum_permeability * perimeter * log((2* distance_closest_normal_surface + gap_length) / gap_length));
        }

        double reluctance = 1. / (1. / reluctance_internal + 1. / reluctance_fringing);

        if (gap_length > 0) {
            fringing_factor = gap_length / (constants.vacuum_permeability * gap_area * reluctance);
        }
        std::map<std::string, double> result;
        result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringing_factor);
        result["reluctance"] = reluctance;
        result["fringing_factor"] = fringing_factor;

        return result;
    };

    std::map<std::string, double> ReluctanceMuehlethalerModel::get_gap_reluctance(CoreGap gapInfo)
    {
        auto constants = Constants();
        auto gap_length = gapInfo.get_length();
        auto gap_shape = *(gapInfo.get_shape());
        auto gap_section_dimensions = *(gapInfo.get_section_dimensions());
        auto distance_closest_normal_surface = *(gapInfo.get_distance_closest_normal_surface());
        double reluctance;
        double fringing_factor = 1;
        auto gap_section_width = gap_section_dimensions[0];
        auto gap_section_depth = gap_section_dimensions[1];

        if (gap_shape == ShapeEnum::ROUND) {
            double gamma_r = get_reluctance_type_1(gap_length / 2, gap_section_width / 2, distance_closest_normal_surface) / (gap_length / constants.vacuum_permeability / (gap_section_width / 2));
            reluctance = pow(gamma_r, 2) * gap_length / (constants.vacuum_permeability * std::numbers::pi * pow(gap_section_width / 2, 2));
            fringing_factor = 1 / gamma_r;
        }
        else {
            double gamma_x = get_reluctance_type_1(gap_length / 2, gap_section_width, distance_closest_normal_surface) / (gap_length / constants.vacuum_permeability / gap_section_width);
            double gamma_y = get_reluctance_type_1(gap_length / 2, gap_section_depth, distance_closest_normal_surface) / (gap_length / constants.vacuum_permeability / gap_section_depth);
            double gamma = gamma_x * gamma_y;
            reluctance = gamma * gap_length / (constants.vacuum_permeability * gap_section_depth * gap_section_width);
            fringing_factor = 1 / gamma;
        }

        std::map<std::string, double> result;
        result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringing_factor);
        result["reluctance"] = reluctance;
        result["fringing_factor"] = fringing_factor;

        return result;
    };

    std::map<std::string, double> ReluctanceEffectiveAreaModel::get_gap_reluctance(CoreGap gapInfo)
    {
        auto constants = Constants();
        auto gap_length = gapInfo.get_length();
        auto gap_area = *(gapInfo.get_area());
        auto gap_shape = *(gapInfo.get_shape());
        auto gap_section_dimensions = *(gapInfo.get_section_dimensions());
        double reluctance;
        double fringing_factor = 1;
        auto gap_section_width = gap_section_dimensions[0];
        auto gap_section_depth = gap_section_dimensions[1];

        if (gap_length > 0) {
            if (gap_shape == ShapeEnum::ROUND) {
                fringing_factor = pow(1 + gap_length / gap_section_width, 2);
            }
            else{
                fringing_factor = (gap_section_depth + gap_length) * (gap_section_width + gap_length) / (gap_section_depth * gap_section_width);
            }
        }

        reluctance = gap_length / (constants.vacuum_permeability * gap_area * fringing_factor);

        std::map<std::string, double> result;
        result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringing_factor);
        result["reluctance"] = reluctance;
        result["fringing_factor"] = fringing_factor;

        return result;
    };

    std::map<std::string, double> ReluctanceEffectiveLengthModel::get_gap_reluctance(CoreGap gapInfo)
    {
        auto constants = Constants();
        auto gap_length = gapInfo.get_length();
        auto gap_area = *(gapInfo.get_area());
        auto gap_shape = *(gapInfo.get_shape());
        auto gap_section_dimensions = *(gapInfo.get_section_dimensions());
        double reluctance;
        double fringing_factor = 1;
        auto gap_section_width = gap_section_dimensions[0];
        auto gap_section_depth = gap_section_dimensions[1];

        if (gap_length > 0) {
            if (gap_shape == ShapeEnum::ROUND) {
                fringing_factor = pow(1 + gap_length / gap_section_width, 2);
            }
            else{
                fringing_factor = (1 + gap_length / gap_section_depth) * (1 + gap_length / gap_section_width);
            }
        }

        reluctance = gap_length / (constants.vacuum_permeability * gap_area * fringing_factor);

        std::map<std::string, double> result;
        result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringing_factor);
        result["reluctance"] = reluctance;
        result["fringing_factor"] = fringing_factor;

        return result;
    };


    std::map<std::string, double> ReluctanceMcLymanModel::get_gap_reluctance(CoreGap gapInfo)
    {
        auto constants = Constants();
        auto gap_length = gapInfo.get_length();
        auto gap_area = *(gapInfo.get_area());
        auto gap_section_dimensions = *(gapInfo.get_section_dimensions());
        auto distance_closest_normal_surface = *(gapInfo.get_distance_closest_normal_surface());
        double reluctance;
        double fringing_factor = 1;

        if (gap_length > 0) {
            fringing_factor = 1 + gap_length / sqrt(gap_area) * log(2 * 2 * distance_closest_normal_surface / gap_length);
        }

        reluctance = gap_length / (constants.vacuum_permeability * gap_area * fringing_factor);

        std::map<std::string, double> result;
        result["maximum_storable_energy"] = get_gap_maximum_storable_energy(gapInfo, fringing_factor);
        result["reluctance"] = reluctance;
        result["fringing_factor"] = fringing_factor;

        return result;
    };


    std::shared_ptr<ReluctanceModel> ReluctanceModel::factory(ReluctanceModels modelName)
    {
        if (modelName == ReluctanceModels::ZHANG) {
            std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctanceZhangModel);
            return reluctanceModel;
        }
        else if (modelName == ReluctanceModels::MCLYMAN) {
            std::shared_ptr<ReluctanceModel> reluctanceModel(new ReluctanceMcLymanModel);
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

        else throw std::runtime_error("Unknown Reluctance mode, available options are: {ZHANG, MCLYMAN, EFFECTIVE_AREA, EFFECTIVE_LENGTH, MUEHLETHALER}");
    }
}
