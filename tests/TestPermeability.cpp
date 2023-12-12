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

SUITE(InitialPermeability) {
    TEST(Test_XFlux_60) {
        OpenMagnetics::InitialPermeability initialPermeability;
        std::string materialName = "XFlux 60";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        CHECK(initialPermeabilityValue == 60);
        double temperature = 89;
        double initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
        double expected = 60 * 1.0073;
        CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, 0.01 * expected);

        double magneticFieldDcBias = 3978.873577;
        double initialPermeabilityValueWithMagneticFieldDcBias =
            initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
        expected = 60 * 0.9601;
        CHECK_CLOSE(initialPermeabilityValueWithMagneticFieldDcBias, expected, 0.01 * expected);

        double frequency = 1250000;
        double initialPermeabilityValueWithFrequency =
            initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
        expected = 60 * 0.968;
        CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, 0.01 * expected);
    }

    TEST(Test_N88) {
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

    TEST(Test_N30) {
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