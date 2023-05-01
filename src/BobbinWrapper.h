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
        auto processor = BobbinDataProcessor::factory(*this);
        
        set_processed_description((*processor).process_data(*this));
    }
    BobbinWrapper() = default;
    virtual ~BobbinWrapper() = default;
    static double get_filling_factor(double windingWindowWidth, double windingWindowHeight);

};
} // namespace OpenMagnetics