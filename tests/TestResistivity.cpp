#include "physical_models/Resistivity.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

SUITE(Resistivity) {
    double maximumError = 0.05;

    TEST(Test_Core_Material_Resistivity_1_Point) {
        auto materialData = find_core_material_by_name("3C97");
        double temperature = 42;

        auto resistivityModel = ResistivityModel::factory(ResistivityModels::CORE_MATERIAL);

        double expectedResistivity = 5;
        auto resistivity = (*resistivityModel).get_resistivity(materialData, temperature);
        CHECK_CLOSE(resistivity, expectedResistivity, expectedResistivity * maximumError);

    }

    TEST(Test_Core_Material_Resistivity_Several_Points) {
        auto materialData = find_core_material_by_name("3C94");
        double temperature = 42;

        auto resistivityModel = ResistivityModel::factory(ResistivityModels::CORE_MATERIAL);

        double expectedResistivity = 2.3;
        auto resistivity = (*resistivityModel).get_resistivity(materialData, temperature);
        CHECK_CLOSE(resistivity, expectedResistivity, expectedResistivity * maximumError);
    }

    TEST(Test_Wire_Material_Resistivity_Copper_20) {
        auto materialData = find_wire_material_by_name("copper");
        double temperature = 20;

        auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);

        double expectedResistivity = 0.00000001678;
        auto resistivity = (*resistivityModel).get_resistivity(materialData, temperature);
        CHECK_CLOSE(resistivity, expectedResistivity, expectedResistivity * maximumError);
    }

    TEST(Test_Wire_Material_Resistivity_Copper_200) {
        auto materialData = find_wire_material_by_name("copper");
        double temperature = 200;

        auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);

        double expectedResistivity = 2.9e-08;
        auto resistivity = (*resistivityModel).get_resistivity(materialData, temperature);
        CHECK_CLOSE(resistivity, expectedResistivity, expectedResistivity * maximumError);
    }

    TEST(Test_Wire_Material_Resistivity_Aluminium_20) {
        auto materialData = find_wire_material_by_name("aluminium");
        double temperature = 20;

        auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);

        double expectedResistivity = 0.0000000265;
        auto resistivity = (*resistivityModel).get_resistivity(materialData, temperature);
        CHECK_CLOSE(resistivity, expectedResistivity, expectedResistivity * maximumError);
    }

    TEST(Test_Wire_Material_Resistivity_Aluminium_200) {
        auto materialData = find_wire_material_by_name("aluminium");
        double temperature = 200;

        auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);

        double expectedResistivity = 4.69e-08;
        auto resistivity = (*resistivityModel).get_resistivity(materialData, temperature);
        CHECK_CLOSE(resistivity, expectedResistivity, expectedResistivity * maximumError);
    }
}