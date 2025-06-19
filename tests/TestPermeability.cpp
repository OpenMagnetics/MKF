#include "physical_models/InitialPermeability.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/AmplitudePermeability.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <random>
using json = nlohmann::json;
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

SUITE(InitialPermeability) {
    TEST(Test_Initial_Permeability_3C97) {
        InitialPermeability initialPermeability;
        std::string materialName = "3C97";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        CHECK(initialPermeabilityValue == 3341.5);
        {
            double magneticFieldDcBias = 60;
            double temperature = 25;
            double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData, temperature, magneticFieldDcBias, std::nullopt);
            double expected = 2310;
            CHECK_CLOSE(initialPermeabilityValue, expected, 0.01 * expected);
        }
    }

    TEST(Test_Initial_Permeability_51) {
        InitialPermeability initialPermeability;
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
        InitialPermeability initialPermeability;
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
        InitialPermeability initialPermeability;
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
        InitialPermeability initialPermeability;
        std::string materialName = "N88";
        auto materialData = find_core_material_by_name(materialName);
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
        InitialPermeability initialPermeability;
        std::string materialName = "N30";
        auto materialData = find_core_material_by_name(materialName);
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

    TEST(Test_Initial_Permeability_X_Indmix_A) {
        InitialPermeability initialPermeability;
        std::string materialName = "X-Indmix A";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData, 85);
        double expected = 60;
        CHECK_CLOSE(expected, initialPermeabilityValue, 0.01 * expected);
    }

    TEST(Test_Initial_Permeability_Nanoperm_1000) {
        InitialPermeability initialPermeability;
        std::string materialName = "Nanoperm 1000";
        auto materialData = find_core_material_by_name(materialName);
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        double manufacturerTolerance = 0.05;

        {
            double expected = 1000;
            CHECK_CLOSE(initialPermeabilityValue, expected, manufacturerTolerance * expected);
        }

        {
            double temperature = 80;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 1000;
            CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);
        }

        {
            double temperature = 2000;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 1;
            CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);
        }

        {
            double frequency = 100000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 1000;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }

        {
            double frequency = 3000000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 570;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }

        {
            double magneticFieldDcBias = 100;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 950;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }

        {
            double magneticFieldDcBias = 1000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 580;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }

        {
            double frequency = 3000000;
            double magneticFieldDcBias = 1000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, frequency);
            double expected = 336;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }
    }

    TEST(Test_Initial_Permeability_Nanoperm_80000) {
        InitialPermeability initialPermeability;
        std::string materialName = "Nanoperm 80000";
        auto materialData = find_core_material_by_name(materialName);
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        double manufacturerTolerance = 0.05;

        {
            double expected = 80000;
            CHECK_CLOSE(initialPermeabilityValue, expected, manufacturerTolerance * expected);
        }

        {
            double temperature = 80;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 80000;
            CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);
        }

        {
            double temperature = 2000;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 1;
            CHECK_CLOSE(initialPermeabilityValueWithTemperature, expected, manufacturerTolerance * expected);
        }

        {
            double frequency = 1000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 80000;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }

        {
            double frequency = 30000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 40000;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }

        {
            double magneticFieldDcBias = 1;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 80000;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }

        {
            double magneticFieldDcBias = 10;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 43000;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }

        {
            double frequency = 30000;
            double magneticFieldDcBias = 10;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, frequency);
            double expected = 21000;
            CHECK_CLOSE(initialPermeabilityValueWithFrequency, expected, manufacturerTolerance * expected);
        }
    }

    TEST(Test_Frequency_For_Initial_Permeability_Drop_Nanoperm_80000) {
        srand (time(NULL));
        InitialPermeability initialPermeability;
        std::string materialName = "Nanoperm 80000";
        auto materialData = find_core_material_by_name(materialName);
        double manufacturerTolerance = 0.05;
        std::uniform_real_distribution<double> unif(0, 1);
        std::default_random_engine re;
        for (size_t i = 0; i < 1000; ++i)
        {
            double percentageDrop = unif(re);
            double frequencyForDrop = initialPermeability.calculate_frequency_for_initial_permeability_drop(materialData, percentageDrop);
            double expectedInitialPermeability = 80000 * (1 - percentageDrop);
            double initialPermeabilityValueWithFrequencyDrop = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequencyForDrop);
            CHECK_CLOSE(initialPermeabilityValueWithFrequencyDrop, expectedInitialPermeability, manufacturerTolerance * 80000);
        }
    }

    TEST(Test_Frequency_For_Initial_Permeability_Drop_XFlux_60) {
        srand (time(NULL));
        InitialPermeability initialPermeability;
        std::string materialName = "XFlux 60";
        auto materialData = find_core_material_by_name(materialName);
        double manufacturerTolerance = 0.05;
        std::uniform_real_distribution<double> unif(0, 1);
        std::default_random_engine re;
        for (size_t i = 0; i < 20; ++i)
        {
            double percentageDrop = unif(re);
            double frequencyForDrop = initialPermeability.calculate_frequency_for_initial_permeability_drop(materialData, percentageDrop);
            double expectedInitialPermeability = 60 * (1 - percentageDrop);
            double initialPermeabilityValueWithFrequencyDrop = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequencyForDrop);
            CHECK_CLOSE(initialPermeabilityValueWithFrequencyDrop, expectedInitialPermeability, manufacturerTolerance * expectedInitialPermeability);
        }
    }
}


SUITE(ComplexPermeability) {
    TEST(Test_Complex_Permeability_N22) {
        ComplexPermeability complexPermeability;
        std::string materialName = "N22";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST(Test_Complex_Permeability_3C97) {
        ComplexPermeability complexPermeability;
        std::string materialName = "3C97";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST(Test_Complex_Permeability_N49) {
        ComplexPermeability complexPermeability;
        std::string materialName = "N49";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST(Test_Complex_Permeability_67) {
        ComplexPermeability complexPermeability;
        std::string materialName = "67";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST(Test_Complex_Permeability_Nanoperm_8000) {
        ComplexPermeability complexPermeability;
        std::string materialName = "Nanoperm 8000";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 1000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST(Test_Complex_Permeability_XFlux_60) {
        ComplexPermeability complexPermeability;
        std::string materialName = "XFlux 60";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 1000000);
        CHECK(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        CHECK(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);

        auto complexPermeabilityValues = complexPermeability.calculate_complex_permeability_from_frequency_dependent_initial_permeability(materialData);
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        {
            OpenMagnetics::Curve2D curve;
            for (auto point : std::get<std::vector<PermeabilityPoint>>(complexPermeabilityValues.get_real())) {
                curve.get_mutable_x_points().push_back(point.get_frequency().value());
                curve.get_mutable_y_points().push_back(point.get_value());
            }

            auto outFile = outputFilePath;
            outFile.append("Test_Complex_Permeability_XFlux_60_Real.svg");

            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_curve(curve, true);
            painter.export_svg();
            CHECK(std::filesystem::exists(outFile));
        }
        {
            OpenMagnetics::Curve2D curve;
            for (auto point : std::get<std::vector<PermeabilityPoint>>(complexPermeabilityValues.get_imaginary())) {
                curve.get_mutable_x_points().push_back(point.get_frequency().value());
                curve.get_mutable_y_points().push_back(point.get_value());
            }

            auto outFile = outputFilePath;
            outFile.append("Test_Complex_Permeability_XFlux_60_Imaginary.svg");

            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_curve(curve, true);
            painter.export_svg();
            CHECK(std::filesystem::exists(outFile));
        }
    }
}

SUITE(AmplitudePermeability) {
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
    TEST(Test_BH_Loop_3C97) {
        std::string materialName = "3C97";
        auto curves = OpenMagnetics::BHLoopRoshenModel().get_hysteresis_loop(materialName, 25, 0.2, std::nullopt);

        auto outFile = outputFilePath;
        outFile.append("Test_BH_Loop_3C97_Upper.svg");

        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(curves.first);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
    }
}