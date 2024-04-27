#include "Settings.h"
#include "CoreMaterialCrossReferencer.h"
#include "Utils.h"
#include "InputsWrapper.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>


SUITE(CoreMaterialCrossReferencer) {
    auto settings = OpenMagnetics::Settings::GetInstance();

    TEST(Test_All_Core_Materials) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::string coreMaterialName = "3C97";
        OpenMagnetics::CoreMaterial coreMaterial = OpenMagnetics::CoreWrapper::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "JNP96A");


        auto scorings = coreMaterialCrossReferencer.get_scorings();
        auto scoredValues = coreMaterialCrossReferencer.get_scored_values();
        json results;
        results["cores"] = json::array();
        results["scorings"] = json::array();
        for (auto& [coreMaterial, scoring] : crossReferencedCoreMaterials) {
            std::string name = coreMaterial.get_name();

            json coreMaterialJson;
            OpenMagnetics::to_json(coreMaterialJson, coreMaterial);
            results["cores"].push_back(coreMaterialJson);
            results["scorings"].push_back(scoring);

            json result;
            result["scoringPerFilter"] = json();
            result["scoredValuePerFilter"] = json();
            for (auto& filter : magic_enum::enum_names<OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>()) {
                std::string filterString(filter);

                result["scoringPerFilter"][filterString] = scorings[name][magic_enum::enum_cast<OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()];
                result["scoredValuePerFilter"][filterString] = scoredValues[name][magic_enum::enum_cast<OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()];
                CHECK(!std::isnan(scorings[name][magic_enum::enum_cast<OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()]));
                CHECK(!std::isnan(scoredValues[name][magic_enum::enum_cast<OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters>(filterString).value()]));
            };
            results["data"].push_back(result);
        }
    }

    TEST(Test_All_Core_Materials_Only_TDK) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("TDK");

        std::string coreMaterialName = "3C97";
        OpenMagnetics::CoreMaterial coreMaterial = OpenMagnetics::CoreWrapper::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "N95");
    }

    TEST(Test_All_Core_Materials_Powder) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::string coreMaterialName = "Kool Mµ MAX 26";
        OpenMagnetics::CoreMaterial coreMaterial = OpenMagnetics::CoreWrapper::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "Kool Mµ Hƒ 26");
    }

    TEST(Test_All_Core_Materials_Powder_Only_Micrometals) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Micrometals");

        std::string coreMaterialName = "Kool Mµ MAX 26";
        OpenMagnetics::CoreMaterial coreMaterial = OpenMagnetics::CoreWrapper::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "SM 26");
    }

    TEST(Test_All_Core_Materials_Only_Volumetric_Losses) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::map<OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters, double> weights;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 1;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::REMANENCE] = 0;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::SATURATION] = 0;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 0.5;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0;

        std::string coreMaterialName = "3C97";
        OpenMagnetics::CoreMaterial coreMaterial = OpenMagnetics::CoreWrapper::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "JNP96A");
    }

    TEST(Test_All_Core_Materials_Only_Volumetric_Losses_Powder) {
        settings->reset();
        OpenMagnetics::clear_databases();
        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::map<OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters, double> weights;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 0;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::REMANENCE] = 0;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::SATURATION] = 0;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 1;
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0;

        std::string coreMaterialName = "Kool Mµ MAX 26";
        OpenMagnetics::CoreMaterial coreMaterial = OpenMagnetics::CoreWrapper::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        CHECK(crossReferencedCoreMaterials.size() > 0);

        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "Kool Mµ Hƒ 125");
    }
}