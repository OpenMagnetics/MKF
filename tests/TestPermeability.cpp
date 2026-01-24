#include "RandomUtils.h"
#include <source_location>
#include "physical_models/InitialPermeability.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/AmplitudePermeability.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
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

namespace { 
    auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");

    TEST_CASE("Test_Initial_Permeability_3C97", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "3C97";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        REQUIRE(initialPermeabilityValue == 3341.5);
        {
            double magneticFieldDcBias = 60;
            double temperature = 25;
            double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData, temperature, magneticFieldDcBias, std::nullopt);
            double expected = 2310;
            REQUIRE_THAT(initialPermeabilityValue, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }
    }

    TEST_CASE("Test_Initial_Permeability_51", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "51";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        REQUIRE(initialPermeabilityValue == 350);
        {
            double temperature = 120;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 896;
            REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }
        {
            double temperature = -20;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 259;
            REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }
    }

    TEST_CASE("Test_Initial_Permeability_Mix_3", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "Mix 3";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        REQUIRE(initialPermeabilityValue == 35);
        {
            double temperature = 89;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 35.6158;
            REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }

        {
            double magneticFieldDcBias = 7957.7471546;
            double initialPermeabilityValueWithMagneticFieldDcBias = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 28.0527;
            REQUIRE_THAT(initialPermeabilityValueWithMagneticFieldDcBias, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }

        {
            double magneticFluxDensity = 0.100;
            double initialPermeabilityValueWithMagneticFluxDensity = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, std::nullopt, magneticFluxDensity);
            double expected = 41.5;
            REQUIRE_THAT(initialPermeabilityValueWithMagneticFluxDensity, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }

        {
            double frequency = 1250000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 34.8;
            REQUIRE_THAT(expected, Catch::Matchers::WithinAbs(initialPermeabilityValueWithFrequency, 0.01 * expected));
        }
        {
            double frequency = 125000000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 32.8;
            REQUIRE_THAT(expected, Catch::Matchers::WithinAbs(initialPermeabilityValueWithFrequency, 0.01 * expected));
        }
    }

    TEST_CASE("Test_Initial_Permeability_XFlux_60", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "XFlux 60";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        REQUIRE(initialPermeabilityValue == 60);
        {
            double temperature = 89;
            double initialPermeabilityValueWithTemperature =
                initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 60 * 1.0073;
            REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }

        {
            double magneticFieldDcBias = 3978.873577;
            double initialPermeabilityValueWithMagneticFieldDcBias =
                initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 60 * 0.9601;
            REQUIRE_THAT(initialPermeabilityValueWithMagneticFieldDcBias, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }

        {
            double frequency = 1250000;
            double initialPermeabilityValueWithFrequency =
                initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 60 * 0.968;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }
    }

    TEST_CASE("Test_Initial_Permeability_NPF_26", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "NPF 26";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        REQUIRE(initialPermeabilityValue == 26);
        {
            double magneticFieldDcBias = 19090.6;
            double initialPermeabilityValueWithMagneticFieldDcBias =
                initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 20.814;
            REQUIRE_THAT(initialPermeabilityValueWithMagneticFieldDcBias, Catch::Matchers::WithinAbs(expected, 0.01 * expected));
        }
    }

    TEST_CASE("Test_Initial_Permeability_N88", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "N88";
        auto materialData = find_core_material_by_name(materialName);
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);

        double expected = 1900;
        double manufacturerTolerance = 0.25;
        REQUIRE_THAT(initialPermeabilityValue, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));

        double temperature = 80;
        double initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
        expected = 3200;
        manufacturerTolerance = 0.25;
        REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));

        temperature = 200;
        initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);

        expected = 4500;
        manufacturerTolerance = 0.25;
        REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));

        temperature = 300;
        initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);

        expected = 1;
        manufacturerTolerance = 0.25;
        REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
    }

    TEST_CASE("Test_Initial_Permeability_N30", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "N30";
        auto materialData = find_core_material_by_name(materialName);
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);

        double expected = 4300;
        double manufacturerTolerance = 0.25;
        REQUIRE_THAT(initialPermeabilityValue, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));

        double temperature = 80;
        double initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
        expected = 4300;
        manufacturerTolerance = 0.25;
        REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));

        temperature = 200;
        initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);

        expected = 1;
        manufacturerTolerance = 0.25;
        REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));

        temperature = 300;
        initialPermeabilityValueWithTemperature =
            initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);

        expected = 1;
        manufacturerTolerance = 0.25;
        REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
    }

    TEST_CASE("Test_Initial_Permeability_X_Indmix_A", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "X-Indmix A";
        auto materialData = materialName;
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData, 85);
        double expected = 60;
        REQUIRE_THAT(expected, Catch::Matchers::WithinAbs(initialPermeabilityValue, 0.01 * expected));
    }

    TEST_CASE("Test_Initial_Permeability_Nanoperm_1000", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "Nanoperm 1000";
        auto materialData = find_core_material_by_name(materialName);
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        double manufacturerTolerance = 0.05;

        {
            double expected = 1000;
            REQUIRE_THAT(initialPermeabilityValue, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double temperature = 80;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 1000;
            REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double temperature = 2000;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 1;
            REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double frequency = 100000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 1000;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double frequency = 3000000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 570;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double magneticFieldDcBias = 100;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 950;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double magneticFieldDcBias = 1000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 580;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double frequency = 3000000;
            double magneticFieldDcBias = 1000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, frequency);
            double expected = 336;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }
    }

    TEST_CASE("Test_Initial_Permeability_Nanoperm_80000", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "Nanoperm 80000";
        auto materialData = find_core_material_by_name(materialName);
        double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData);
        double manufacturerTolerance = 0.05;

        {
            double expected = 80000;
            REQUIRE_THAT(initialPermeabilityValue, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double temperature = 80;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 80000;
            REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double temperature = 2000;
            double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(materialData, temperature, std::nullopt, std::nullopt);
            double expected = 1;
            REQUIRE_THAT(initialPermeabilityValueWithTemperature, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double frequency = 1000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 80000;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double frequency = 30000;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequency);
            double expected = 40000;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double magneticFieldDcBias = 1;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 80000;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double magneticFieldDcBias = 10;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, std::nullopt);
            double expected = 43000;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }

        {
            double frequency = 30000;
            double magneticFieldDcBias = 10;
            double initialPermeabilityValueWithFrequency = initialPermeability.get_initial_permeability(materialData, std::nullopt, magneticFieldDcBias, frequency);
            double expected = 21000;
            REQUIRE_THAT(initialPermeabilityValueWithFrequency, Catch::Matchers::WithinAbs(expected, manufacturerTolerance * expected));
        }
    }

    TEST_CASE("Test_Frequency_For_Initial_Permeability_Drop_Nanoperm_80000", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "Nanoperm 80000";
        auto materialData = find_core_material_by_name(materialName);
        double manufacturerTolerance = 0.05;
        std::uniform_real_distribution<double> unif(0, 0.9);  // Limit to 90% drop to avoid edge cases
        std::default_random_engine re(42);  // Fixed seed for cross-platform reproducibility
        for (size_t i = 0; i < 1000; ++i)
        {
            double percentageDrop = unif(re);
            double frequencyForDrop = initialPermeability.calculate_frequency_for_initial_permeability_drop(materialData, percentageDrop);
            double expectedInitialPermeability = 80000 * (1 - percentageDrop);
            double initialPermeabilityValueWithFrequencyDrop = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequencyForDrop);
            REQUIRE_THAT(initialPermeabilityValueWithFrequencyDrop, Catch::Matchers::WithinAbs(expectedInitialPermeability, manufacturerTolerance * 80000));
        }
    }

    TEST_CASE("Test_Frequency_For_Initial_Permeability_Drop_XFlux_60", "[physical-model][initial-permeability][smoke-test]") {
        InitialPermeability initialPermeability;
        std::string materialName = "XFlux 60";
        auto materialData = find_core_material_by_name(materialName);
        double manufacturerTolerance = 0.05;
        std::uniform_real_distribution<double> unif(0, 0.9);  // Limit to 90% drop to avoid edge cases
        std::default_random_engine re(42);  // Fixed seed for cross-platform reproducibility
        for (size_t i = 0; i < 20; ++i)
        {
            double percentageDrop = unif(re);
            double frequencyForDrop = initialPermeability.calculate_frequency_for_initial_permeability_drop(materialData, percentageDrop);
            double expectedInitialPermeability = 60 * (1 - percentageDrop);
            double initialPermeabilityValueWithFrequencyDrop = initialPermeability.get_initial_permeability(materialData, std::nullopt, std::nullopt, frequencyForDrop);
            REQUIRE_THAT(initialPermeabilityValueWithFrequencyDrop, Catch::Matchers::WithinAbs(expectedInitialPermeability, manufacturerTolerance * expectedInitialPermeability));
        }
    }

    TEST_CASE("Test_Complex_Permeability_N22", "[physical-model][complex-permeability][smoke-test]") {
        ComplexPermeability complexPermeability;
        std::string materialName = "N22";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        REQUIRE(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        REQUIRE(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST_CASE("Test_Complex_Permeability_3C97", "[physical-model][complex-permeability][smoke-test]") {
        ComplexPermeability complexPermeability;
        std::string materialName = "3C97";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        REQUIRE(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        REQUIRE(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST_CASE("Test_Complex_Permeability_N49", "[physical-model][complex-permeability][smoke-test]") {
        ComplexPermeability complexPermeability;
        std::string materialName = "N49";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        REQUIRE(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        REQUIRE(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST_CASE("Test_Complex_Permeability_67", "[physical-model][complex-permeability][smoke-test]") {
        ComplexPermeability complexPermeability;
        std::string materialName = "67";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 10000000);
        REQUIRE(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        REQUIRE(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST_CASE("Test_Complex_Permeability_Nanoperm_8000", "[physical-model][complex-permeability][smoke-test]") {
        ComplexPermeability complexPermeability;
        std::string materialName = "Nanoperm 8000";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 1000000);
        REQUIRE(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        REQUIRE(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);
    }

    TEST_CASE("Test_Complex_Permeability_XFlux_60", "[physical-model][complex-permeability][smoke-test]") {
        ComplexPermeability complexPermeability;
        std::string materialName = "XFlux 60";
        auto materialData = materialName;
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, 1000000);
        REQUIRE(complexPermeabilityValueAt100000.first > complexPermeabilityValueAt10000000.first);
        REQUIRE(complexPermeabilityValueAt100000.second < complexPermeabilityValueAt10000000.second);

        auto complexPermeabilityValues = complexPermeability.calculate_complex_permeability_from_frequency_dependent_initial_permeability(materialData);
        auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
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
            REQUIRE(std::filesystem::exists(outFile));
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
            REQUIRE(std::filesystem::exists(outFile));
        }
    }

    TEST_CASE("Test_BH_Loop_3C97", "[physical-model][amplitude-permeability][smoke-test]") {
        std::string materialName = "3C97";
        auto curves = OpenMagnetics::BHLoopRoshenModel().get_hysteresis_loop(materialName, 25, 0.2, std::nullopt);

        auto outFile = outputFilePath;
        outFile.append("Test_BH_Loop_3C97_Upper.svg");

        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(curves.first);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
    }

}  // namespace
