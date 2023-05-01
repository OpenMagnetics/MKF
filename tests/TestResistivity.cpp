#include "Resistivity.h"
#include "Utils.h"
#include "CoreWrapper.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

SUITE(Resistivity) {
    double maximumError = 0.05;

    TEST(Test_Core_Material_Resistivity_1_Point) {
        auto material_data = OpenMagnetics::find_core_material_by_name("3C97");
        double temperature = 42;

        auto resistivityModel = OpenMagnetics::ResistivityModel::factory(OpenMagnetics::ResistivityModels::CORE_MATERIAL);

        double expectedResistivity = 5;
        auto resistivity = (*resistivityModel).get_resistivity(material_data, temperature);
        CHECK_CLOSE(resistivity, expectedResistivity, expectedResistivity * maximumError);

    }
    TEST(Test_Core_Material_Resistivity_Several_Points) {
        auto material_data = OpenMagnetics::find_core_material_by_name("3C94");
        double temperature = 42;

        auto resistivityModel = OpenMagnetics::ResistivityModel::factory(OpenMagnetics::ResistivityModels::CORE_MATERIAL);

        double expectedResistivity = 2.3;
        auto resistivity = (*resistivityModel).get_resistivity(material_data, temperature);
        CHECK_CLOSE(resistivity, expectedResistivity, expectedResistivity * maximumError);

    }
}