#pragma once

#include "json.hpp"

#include "CoreWrapper.h"
#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
using json = nlohmann::json;

namespace OpenMagnetics {

class BobbinDataProcessor{
    public:
        BobbinDataProcessor() = default;
        virtual ~BobbinDataProcessor() = default;
        virtual CoreBobbinProcessedDescription process_data(Bobbin bobbin) = 0;

        static std::shared_ptr<BobbinDataProcessor> factory(Bobbin bobbin);
};

class BobbinWrapper : public Bobbin {
  private:
  public:
    BobbinWrapper(const json& j) {
        from_json(j, *this);
        if (get_functional_description()) {
            auto processor = BobbinDataProcessor::factory(*this);
            set_processed_description((*processor).process_data(*this));
        }
        else {
            if (!get_processed_description()) {
                throw std::runtime_error("Missing data in bobbin");
            }
        }
    }
    BobbinWrapper() = default;
    virtual ~BobbinWrapper() = default;
    
    static double get_filling_factor(double windingWindowWidth, double windingWindowHeight);
    static std::vector<double> get_winding_window_dimensions(double coreWindingWindowWidth, double coreWindingWindowHeight);
    std::vector<double> get_winding_window_dimensions(size_t windingWindowIndex = 0);
    std::vector<double> get_winding_window_coordinates(size_t windingWindowIndex = 0);
    WindingOrientation get_winding_window_sections_orientation(size_t windingWindowIndex = 0);
    static BobbinWrapper create_quick_bobbin(double windingWindowHeight, double windingWindowWidth);
    static BobbinWrapper create_quick_bobbin(CoreWrapper core, bool nullDimensions=false);

};
} // namespace OpenMagnetics