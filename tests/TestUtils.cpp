#include "InitialPermeability.h"
#include "Utils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>

SUITE(Utils) {
    TEST(LoadDatabaseJson) {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/")).append("/masData.json");
        std::cout << masPath << std::endl;

        std::ifstream ifs(masPath);
        json masData = json::parse(ifs);

        // for (auto& element : masData["coreShapes"].items()) {
        //     json jf = element.value();
        // std::cout << "shapes: " << element.key() << std::endl;
        // std::cout << jf << std::endl;
        //     int familySubtype = 1;
        //     if (jf["familySubtype"] != nullptr) {
        //         familySubtype = jf["familySubtype"];
        //     }
        //     jf["familySubtype"] = std::to_string(familySubtype);
        //     std::string family = jf["family"];
        //     std::transform(family.begin(), family.end(), family.begin(), ::toupper);
        //     jf["family"] = family;
            
        //     OpenMagnetics::CoreShape coreShape(jf);

        // }

        // for (auto& element : masData["wires"].items()) {
        //     json jf = element.value();
        // std::cout << "wires: " << element.key() << std::endl;
        // std::cout << jf << std::endl;
        //     std::string standardName;
        //     if (jf["standardName"].is_string()){
        //         standardName = jf["standardName"];
        //     }
        //     else {
        //         standardName = std::to_string(jf["standardName"].get<int>());
        //     }

        //     jf["standardName"] = standardName;

        //     OpenMagnetics::WireS wire(jf);

        //     wireDatabase[jf["name"]] = wire;
        // }
        // for (auto& element : masData["bobbins"].items()) {
        //     json jf = element.value();
        // std::cout << "bobbins: " << element.key() << std::endl;
        // std::cout << jf << std::endl;

        //     OpenMagnetics::BobbinWrapper bobbin(jf);
        //     bobbinDatabase[jf["name"]] = bobbin;
        // }

        // for (auto& element : masData["insulationMaterials"].items()) {
        //     json jf = element.value();
        //     std::cout << "insulationMaterials: " << element.key() << std::endl;
        //     std::cout << jf << std::endl;
        //     OpenMagnetics::InsulationMaterialWrapper insulationMaterial(jf);
        //     insulationMaterialDatabase[jf["name"]] = insulationMaterial;
        // }

        // for (auto& element : masData["wireMaterials"].items()) {
        //     json jf = element.value();
        //     std::cout << "wireMaterials: " << element.key() << std::endl;
        //     std::cout << jf << std::endl;
        //     OpenMagnetics::WireMaterial wireMaterial(jf);
        //     wireMaterialDatabase[jf["name"]] = wireMaterial;
        // }


        OpenMagnetics::load_databases(masData, true);
    }

}