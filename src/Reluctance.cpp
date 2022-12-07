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

    class ZhangModel: public Reluctance
    {
        // Based on Improved Calculation Method for Inductance Value of the Air-Gap Inductor by Xinsheng Zhang
        // https://sci-hub.wf/https://ieeexplore.ieee.org/document/9332553
        public:
            std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo)
            {
                double perimeter = 0;
                double gap_length = gapInfo.get_length();
                double gap_area = *(gapInfo.get_area());
                double gap_shape = *(gapInfo.get_shape());
                double gap_section_dimensions = *(gapInfo.get_section_dimensions());
                double distance_closest_normal_surface = *(gapInfo.get_distance_closest_normal_surface());
                double reluctance_internal = gap_length / (constants.vacuum_permeability * gap_area);
                double reluctance_fringing = 0;
                double fringing_factor = 1;
                double gap_section_width = gap_section_dimensions[0];
                double gap_section_depth = gap_section_dimensions[1];

                if (gap_shape == ShapeEnum::ROUND) {
                    perimeter = std::numbers::pi * gap_section_width;
                }
                else { // TODO: Properly calcualte perimeter for all shapes
                    perimeter = column_width * 2 + gap_section_depth * 2;
                }

                if (gap_length > 0) {
                    reluctance_fringing = std::numbers::pi / (constants.vacuum_permeability * perimeter * log((distance_closest_normal_surface + gap_length) / gap_length));
                }

                double reluctance = 1 / (1 / reluctance_internal + 1 / reluctance_fringing);

                if (gap_length > 0) {
                    fringing_factor = gap_length / (constants.vacuum_permeability * gap_area * reluctance)
                }
                std::map<std::string, double> result;
                result["reluctance"] = reluctance
                result["fringing_factor"] = fringing_factor

                return result;
            }
    };
