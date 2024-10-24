#include "InitialPermeability.h"
#include "ComplexPermeability.h"
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

SUITE(InitialPermeability) {
    TEST(Test_Initial_Permeability_51) {
        OpenMagnetics::InitialPermeability initialPermeability;
        std::string materialName = "51";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        CHECK(initialPermeabilityValue == 350);
        {
            double temperature = 120;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 896;
            CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, 0.01 * expected);
        }
        {
            double temperature = -20;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 259;
            CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, 0.01 * expected);
        }
    }

    TEST(Test_Initial_Permeability_Mix_3) {
        OpenMagnetics::InitialPermeability initialPermeability;
        std::string materialName = "Mix 3";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        CHECK(initialPermeabilityValue == 35);
        {
            double temperature = 89;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 35.6158;
            CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, 0.01 * expected);
        }

        {
            double magneticFieldDcBias = 7957.7471546;
            double initialPermeabilityValueWithMagneticFieldDcBias = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 28.0527;
            CHECK_CLOSE(initialPermeabilityValueWithMagneticFieldDcBias, expected, 0.01 * expected);
        }

        {
            double magneticFluxDensity = 0.100;
            double initialPermeabilityValueWithMagneticFluxDensity = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, std::nullopt, magneticFluxDensity);
            double expected = 41.5;
            CHECK_CLOSE(initialPermeabilityValueWithMagneticFluxDensity, expected, 0.01 * expected);
        }

        {
            double frequency = 1250000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 34.8;
            CHECK_CLOSE(expected, initialPermeabilityValueWithFrequency, 0.01 * expected);
        }
        {
            double frequency = 125000000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 32.8;
            CHECK_CLOSE(expected, initialPermeabilityValueWithFrequency, 0.01 * expected);
        }
    }

    TEST(Test_Initial_Permeability_XFlux_60) {
        OpenMagnetics::InitialPermeability initialPermeability;
        std::string materialName = "XFlux 60";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        CHECK(initialPermeabilityValue == 60);
        {
            double temperature = 89;
            double initialPermeabilityValueWithTemperature =
                initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 60 * 1.0073;
            CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, 0.01 * expected);
        }

        {
            double magneticFieldDcBias = 3978.873577;
            double initialPermeabilityValueWithMagneticFieldDcBias =
                initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 60 * 0.9601;
            CHECK_CLOSE(initialPermeabilityValueWithMagneticFieldDcBias, expected, 0.01 * expected);
        }

        {
            double frequency = 1250000;
            double initialPermeabilityValueWithFrequency =
                initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 60 * 0.968;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, 0.01 * expected);
        }
    }

    TEST(Test_Initial_Permeability_N88) {
        OpenMagnetics::InitialPermeability initialPermeability;
        std::string materialName = "N88";
        auto materialData = OpenMagnetics::find_core_material_by_name(materialName);
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);

        double expected = 1900;
        double manufacturerTolerance = 0.25;
        CHECK_CLOSE(initialPermeabilityValue, expected, manufacturerTolerance * expected);

        double temperature = 80;
        double initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
        expected = 3200;
        manufacturerTolerance = 0.25;
        CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);

        temperature = 200;
        initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);

        expected = 4500;
        manufacturerTolerance = 0.25;
        CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);

        temperature = 300;
        initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);

        expected = 1;
        manufacturerTolerance = 0.25;
        CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);
    }

    TEST(Test_Initial_Permeability_N30) {
        OpenMagnetics::InitialPermeability initialPermeability;
        std::string materialName = "N30";
        auto materialData = OpenMagnetics::find_core_material_by_name(materialName);
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);

        double expected = 4300;
        double manufacturerTolerance = 0.25;
        CHECK_CLOSE(initialPermeabilityValue, expected, manufacturerTolerance * expected);

        double temperature = 80;
        double initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
        expected = 4300;
        manufacturerTolerance = 0.25;
        CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);

        temperature = 200;
        initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);

        expected = 1;
        manufacturerTolerance = 0.25;
        CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);

        temperature = 300;
        initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);

        expected = 1;
        manufacturerTolerance = 0.25;
        CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);
    }
}


SUITE(ComplexPermeability) {
    TEST(Test_Complex_Permeability_N22) {
        OpenMagnetics::ComplexPermeability complexPermeability;
        std::string materialName = "N22";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST(Test_Complex_Permeability_3C97) {
        OpenMagnetics::ComplexPermeability complexPermeability;
        std::string materialName = "3C97";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST(Test_Complex_Permeability_N49) {
        OpenMagnetics::ComplexPermeability complexPermeability;
        std::string materialName = "N49";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST(Test_Complex_Permeability_67) {
        OpenMagnetics::ComplexPermeability complexPermeability;
        std::string materialName = "67";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }
}