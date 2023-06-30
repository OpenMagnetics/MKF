#include "WindingSkinEffectLosses.h"
#include "Utils.h"
#include "CoreWrapper.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

SUITE(WindingSkinEffectLosses) {
    double maximumError = 0.01;

    TEST(Test_Skin_Depth_Wire_Material_Data_20C) {
        auto materialData = OpenMagnetics::find_wire_material_by_name("copper");
        double frequency = 123000;
        double temperature = 20;

        auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLossesModel();  // TODO change to factory
        auto skinDepth = windingSkinEffectLossesModel.get_skin_depth(materialData, frequency, temperature);
        double expectedSkinDepth = 186.09e-6;
        CHECK_CLOSE(skinDepth, expectedSkinDepth, expectedSkinDepth * maximumError);
    }

    TEST(Test_Skin_Depth_Wire_Name_20C) {
        double frequency = 123000;
        double temperature = 20;

        auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLossesModel();  // TODO change to factory
        auto skinDepth = windingSkinEffectLossesModel.get_skin_depth("copper", frequency, temperature);
        double expectedSkinDepth = 186.09e-6;
        CHECK_CLOSE(skinDepth, expectedSkinDepth, expectedSkinDepth * maximumError);
    }

    TEST(Test_Skin_Depth_Wire_Material_Data_200C) {
        auto materialData = OpenMagnetics::find_wire_material_by_name("copper");
        double frequency = 123000;
        double temperature = 120;

        auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLossesModel();  // TODO change to factory
        auto skinDepth = windingSkinEffectLossesModel.get_skin_depth(materialData, frequency, temperature);
        double expectedSkinDepth = 220e-6;
        CHECK_CLOSE(skinDepth, expectedSkinDepth, expectedSkinDepth * maximumError);
    }
}