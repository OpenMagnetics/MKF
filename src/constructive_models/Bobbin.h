#pragma once

#include "json.hpp"

#include "constructive_models/Core.h"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>

using json = nlohmann::json;

namespace OpenMagnetics {
class Bobbin : public MAS::Bobbin {
  private:
  public:
    Bobbin(const json& j, bool includeProcessedDescription = true) {
        from_json(j, *this);
        if (get_functional_description()) {
            if (includeProcessedDescription) {
                process_data();
            }
        }
        // else {
        //     if (!get_processed_description()) {
        //         throw std::runtime_error("Missing data in bobbin");
        //     }
        // }
    }

    Bobbin(const MAS::Bobbin bobbin) {
        set_functional_description(bobbin.get_functional_description());

        if (bobbin.get_processed_description()) {
            set_processed_description(bobbin.get_processed_description());
        }
        if (bobbin.get_distributors_info()) {
            set_distributors_info(bobbin.get_distributors_info());
        }
        if (bobbin.get_manufacturer_info()) {
            set_manufacturer_info(bobbin.get_manufacturer_info());
        }
        if (bobbin.get_name()) {
            set_name(bobbin.get_name());
        }
    }


    Bobbin() = default;
    virtual ~Bobbin() = default;
    void process_data();
    
    static double get_filling_factor(double windingWindowWidth, double windingWindowHeight);
    static std::vector<double> get_winding_window_dimensions(double coreWindingWindowWidth, double coreWindingWindowHeight);
    std::vector<double> get_winding_window_dimensions(size_t windingWindowIndex = 0);
    double get_winding_window_area(size_t windingWindowIndex = 0);
    WindingWindowShape get_winding_window_shape(size_t windingWindowIndex = 0);
    std::vector<double> get_winding_window_coordinates(size_t windingWindowIndex = 0);
    WindingOrientation get_winding_window_sections_orientation(size_t windingWindowIndex = 0);
    CoilAlignment get_winding_window_sections_alignment(size_t windingWindowIndex = 0);
    static Bobbin create_quick_bobbin(double windingWindowHeight, double windingWindowWidth);
    static Bobbin create_quick_bobbin(Core core, bool nullDimensions = false);
    static Bobbin create_quick_bobbin(Core core, double thickness);
    static Bobbin create_quick_bobbin(Core core, double wallThickness, double columnThickness);
    bool check_if_fits(double dimension, bool isHorizontalOrRadial = true, size_t windingWindowIndex = 0);
    void set_winding_orientation(WindingOrientation windingOrientation, size_t windingWindowIndex = 0);
    std::optional<WindingOrientation> get_winding_orientation(size_t windingWindowIndex = 0);
    std::vector<double> get_maximum_dimensions();
};



class BobbinDataProcessor{
    public:
        BobbinDataProcessor() = default;
        virtual ~BobbinDataProcessor() = default;
        virtual CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) = 0;

        static std::shared_ptr<BobbinDataProcessor> factory(OpenMagnetics::Bobbin bobbin);
};

} // namespace OpenMagnetics