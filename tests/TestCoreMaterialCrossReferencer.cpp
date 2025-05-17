#include "support/Settings.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "support/Utils.h"
#include "processors/Inputs.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;


SUITE(CoreMaterialCrossReferencer) {
    auto settings = Settings::GetInstance();

    TEST(Test_All_Core_Materials) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "PC95");


        auto scorings = coreMaterialCrossReferencer.get_scorings();
        auto scoredValues = coreMaterialCrossReferencer.get_scored_values();
        json results;
        results["cores"] = json::array();
        results["scorings"] = json::array();
        for (auto& [coreMaterial, scoring] : crossReferencedCoreMaterials) {
            std::string name = coreMaterial.get_name();

            json coreMaterialJson;
            to_json(coreMaterialJson, coreMaterial);
            results["cores"].push_back(coreMaterialJson);
            results["scorings"].push_back(scoring);

            json result;
            result["scoringPerFilter"] = json();
            result["scoredValuePerFilter"] = json();
            for (auto& filter : magic_enum::enum_names<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>()) {
                std::string filterString(filter);

                result["scoringPerFilter"][filterString] = scorings[name][magic_enum::enum_cast<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()];
                result["scoredValuePerFilter"][filterString] = scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()];
                CHECK(!std::isnan(scorings[name][magic_enum::enum_cast<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()]));
                CHECK(!std::isnan(scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()]));
            };
            results["data"].push_back(result);
        }
    }

    TEST(Test_All_Core_Materials_Only_TDK) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("TDK");

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "PC95");
    }

    TEST(Test_All_Core_Materials_Powder) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::string coreMaterialName = "Kool Mµ MAX 26";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "Kool Mµ Hƒ 26");
    }

    TEST(Test_All_Core_Materials_Powder_Only_Micrometals) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Micrometals");

        std::string coreMaterialName = "Kool Mµ MAX 26";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "SM 40");
    }

    TEST(Test_All_Core_Materials_Powder_Only_Micrometals_Ferrite) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Micrometals");

        std::string coreMaterialName = "3C95";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "MS 160");
    }

    TEST(Test_All_Core_Materials_Only_Volumetric_Losses) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::map<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters, double> weights;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 1;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::REMANENCE] = 0;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::SATURATION] = 0;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 0.5;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0;

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "JNP96A");
    }

    TEST(Test_All_Core_Materials_Only_Volumetric_Losses_Powder) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::map<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters, double> weights;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 0;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::REMANENCE] = 0;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::SATURATION] = 0;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 1;
        weights[CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0;

        std::string coreMaterialName = "Kool Mµ MAX 26";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "Mix 40");
    }

    TEST(Test_All_Core_Materials_Only_Fair_Rite) {
        settings->reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Fair-Rite");
        settings->set_use_only_cores_in_stock(false);

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 50, 20);

        auto scorings = coreMaterialCrossReferencer.get_scorings();
        auto scoredValues = coreMaterialCrossReferencer.get_scored_values();
        json results;
        results["coreMaterials"] = json::array();
        results["scorings"] = json::array();
        results["data"] = json::array();

        for (auto& [coreMaterial, scoring] : crossReferencedCoreMaterials) {
            std::string name = coreMaterial.get_name();

            json coreMaterialJson;
            to_json(coreMaterialJson, coreMaterial);
            results["coreMaterials"].push_back(coreMaterialJson);
            results["scorings"].push_back(scoring);

            json result;
            result["scoringPerFilter"] = json();
            result["scoredValuePerFilter"] = json();
            for (auto& filter : magic_enum::enum_names<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>()) {
                std::string filterString(filter);
                result["scoringPerFilter"][filterString] = scorings[name][magic_enum::enum_cast<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()];
                result["scoredValuePerFilter"][filterString] = scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()];
                std::cout << "name: " << name << std::endl;
                std::cout << "scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()]: " << scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()] << std::endl;
            };
            results["data"].push_back(result);
        }

        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "95");
    }

}