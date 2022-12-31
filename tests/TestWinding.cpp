#include <UnitTest++.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>
#include <nlohmann/json-schema.hpp>
#include "json.hpp"
#include "WindingWrapper.h"
#include <magic_enum.hpp>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;
#include <typeinfo>


SUITE(WindingProcessedDescription)
{
    std::string filePath = __FILE__;
    auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");

}

SUITE(WindingGeometricalDescription)
{
    std::string filePath = __FILE__;
    auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");
}

SUITE(WindingFunctionalDescription)
{
    std::string filePath = __FILE__;
    auto masPath = filePath.substr(0, filePath.rfind("/")).append("/../../MAS/");

    TEST(inductor_42_turns)
    {
        auto windingFilePath = masPath + "samples/magnetic/winding/inductor_42_turns.json";
        std::ifstream json_file(windingFilePath);

        auto windingJson = json::parse(json_file);

        OpenMagnetics::WindingWrapper winding(windingJson);

        auto function_description = winding.get_functional_description()[0];

        json windingWrapperJson;
        to_json(windingWrapperJson, function_description);

        CHECK(windingWrapperJson == windingJson["functionalDescription"][0]);
    }

}