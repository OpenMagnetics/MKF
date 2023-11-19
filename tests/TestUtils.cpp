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

    // TEST(Modified_Bessel) {
    //     double calculatedValue = OpenMagnetics::modified_bessel_first_kind(0.0, std::complex<double>{1.0, 0.0}).real();
    //     double expectedValue = 1.2660658777520084;
    //     std::cout << "calculatedValue: " << calculatedValue << std::endl;
    //     std::cout << "expectedValue: " << expectedValue << std::endl;
    //     CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    // }


    TEST(Bessel) {
        double calculatedValue = OpenMagnetics::bessel_first_kind(0.0, std::complex<double>{1.0, 0.0}).real();
        double expectedValue = 0.7651976865579666;
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);

        double calculatedBerValue = OpenMagnetics::kelvin_function_real(0.0, 1.0);
        double expectedBerValue = 0.98438178;
        CHECK_CLOSE(expectedBerValue, calculatedBerValue, expectedBerValue * 0.001);

        double calculatedBeiValue = OpenMagnetics::kelvin_function_imaginary(0.0, 1.0);
        double expectedBeiValue = 0.24956604;
        CHECK_CLOSE(expectedBeiValue, calculatedBeiValue, expectedBeiValue * 0.001);

        double calculatedBerpValue = OpenMagnetics::derivative_kelvin_function_real(0.0, 1.0);
        double expectedBerpValue = -0.06244575217903096;
        CHECK_CLOSE(expectedBerpValue, calculatedBerpValue, fabs(expectedBerpValue) * 0.001);

        double calculatedBeipValue = OpenMagnetics::derivative_kelvin_function_imaginary(0.0, 1.0);
        double expectedBeipValue = 0.49739651146809727;
        CHECK_CLOSE(expectedBeipValue, calculatedBeipValue, expectedBeipValue * 0.001);

    }
}