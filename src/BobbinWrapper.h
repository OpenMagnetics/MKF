#pragma once

#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

namespace OpenMagnetics {
using nlohmann::json;


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
    static std::pair<double, double> get_winding_window_dimensions(double coreWindingWindowWidth, double coreWindingWindowHeight);
    static BobbinWrapper create_quick_bobbin(double windingWindowHeight, double windingWindowWidth);
    static BobbinWrapper create_quick_bobbin(MagneticCore core);

};
} // namespace OpenMagnetics