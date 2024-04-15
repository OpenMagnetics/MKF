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
        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "3C92");
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
        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "PC47");
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
        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "Kool Mµ MAX 40");
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
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::DENSITY] = 0;

        std::string coreMaterialName = "3C97";
        OpenMagnetics::CoreMaterial coreMaterial = OpenMagnetics::CoreWrapper::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        CHECK(crossReferencedCoreMaterials.size() > 0);
        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "High DC Bias XFlux 40");
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
        weights[OpenMagnetics::CoreMaterialCrossReferencer::CoreMaterialCrossReferencerFilters::DENSITY] = 0;

        std::string coreMaterialName = "Kool Mµ MAX 26";
        OpenMagnetics::CoreMaterial coreMaterial = OpenMagnetics::CoreWrapper::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        CHECK(crossReferencedCoreMaterials.size() > 0);
        CHECK(crossReferencedCoreMaterials[0].first.get_name() == "Mix 35");
    }
}