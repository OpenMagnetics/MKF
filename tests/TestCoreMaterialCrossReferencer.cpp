#include "support/Settings.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "support/Utils.h"
#include "processors/Inputs.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace { 
    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer(std::map<std::string, std::string>{{"coreLosses", "STEINMETZ"}});

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "P45");

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
            for (auto& filter : magic_enum::enum_names<CoreMaterialCrossReferencerFilters>()) {
                std::string filterString(filter);

                result["scoringPerFilter"][filterString] = scorings[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()];
                result["scoredValuePerFilter"][filterString] = scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()];
                REQUIRE(!std::isnan(scorings[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()]));
                REQUIRE(!std::isnan(scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()]));
            };
            results["data"].push_back(result);
        }
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Only_TDK", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("TDK");

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "N95");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Powder", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::string coreMaterialName = "Kool Mµ MAX 26";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "Kool Mµ Hƒ 26");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Powder_Only_Micrometals", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Micrometals");

        std::string coreMaterialName = "Kool Mµ MAX 26";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "SM 40");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Powder_Only_Micrometals_Ferrite", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        settings.set_core_cross_referencer_allow_different_core_material_type(true);
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Micrometals");

        std::string coreMaterialName = "3C95";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "MS 160");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Only_Volumetric_Losses", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::map<CoreMaterialCrossReferencerFilters, double> weights;
        weights[CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 1;
        weights[CoreMaterialCrossReferencerFilters::REMANENCE] = 0;
        weights[CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0;
        weights[CoreMaterialCrossReferencerFilters::SATURATION] = 0;
        weights[CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0;
        weights[CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 0.5;
        weights[CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0;

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "JNP96A");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Only_Volumetric_Losses_Powder", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::map<CoreMaterialCrossReferencerFilters, double> weights;
        weights[CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 0;
        weights[CoreMaterialCrossReferencerFilters::REMANENCE] = 0;
        weights[CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0;
        weights[CoreMaterialCrossReferencerFilters::SATURATION] = 0;
        weights[CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0;
        weights[CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 1;
        weights[CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0;

        std::string coreMaterialName = "Kool Mµ MAX 26";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "Kool Mµ MAX 40");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Only_Fair_Rite", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Fair-Rite");
        settings.set_use_only_cores_in_stock(false);

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
            for (auto& filter : magic_enum::enum_names<CoreMaterialCrossReferencerFilters>()) {
                std::string filterString(filter);
                result["scoringPerFilter"][filterString] = scorings[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()];
                result["scoredValuePerFilter"][filterString] = scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()];
            };
            results["data"].push_back(result);
        }

        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "95");
    }

}  // namespace
