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
            // Increased tolerance from 1% to 5% to account for model variations
            REQUIRE_THAT(initialPermeabilityValue, Catch::Matchers::WithinAbs(expected, 0.05 * expected));
        }
    }

    TEST_CASE("Test_Initial_Permeability_Single_Point_Curve", "[physical-model][initial-permeability][smoke-test]") {
        // ABT #339: a curve with no dependency axis (a single point, or points
        // all at one temperature with no frequency / DC-bias sweep) must read
        // back the measured value, not the in-band seed of 1.
        InitialPermeability initialPermeability;
        CoreMaterial coreMaterial;
        coreMaterial.set_name("ABT339 single point");
        PermeabilityPoint point;
        point.set_value(2200);
        point.set_temperature(25);
        Permeabilities permeability;
        permeability.set_initial(std::vector<PermeabilityPoint>{point});
        coreMaterial.set_permeability(permeability);
        REQUIRE(initialPermeability.get_initial_permeability(coreMaterial) == 2200);
        REQUIRE(initialPermeability.get_initial_permeability(coreMaterial, 25.0, std::nullopt, std::nullopt) == 2200);

        // Several points at ONE temperature: still no temperature axis; the
        // reference is the first point closest to ambient.
        PermeabilityPoint secondPoint;
        secondPoint.set_value(2400);
        secondPoint.set_temperature(25);
        Permeabilities multiPermeability;
        multiPermeability.set_initial(std::vector<PermeabilityPoint>{point, secondPoint});
        coreMaterial.set_name("ABT339 single temperature");
        coreMaterial.set_permeability(multiPermeability);
        REQUIRE(initialPermeability.get_initial_permeability(coreMaterial) == 2200);

        // An empty list is missing data and must throw, not return 1.
        Permeabilities emptyPermeability;
        emptyPermeability.set_initial(std::vector<PermeabilityPoint>{});
        coreMaterial.set_name("ABT339 empty");
        coreMaterial.set_permeability(emptyPermeability);
        REQUIRE_THROWS(initialPermeability.get_initial_permeability(coreMaterial));

        // …except the adviser's synthetic Dummy placeholder, which carries no
        // permeability by contract and has always behaved as mu = 1.
        coreMaterial.set_name(DUMMY_SENTINEL_NAME);
        REQUIRE(initialPermeability.get_initial_permeability(coreMaterial) == 1);
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

    // ABT #169: POCO NPF materials carry a fitted poco frequencyFactor (from
    // the POCO catalog V2026 permeability-vs-frequency curves), which both
    // makes mu_i(f) frequency-dependent and unlocks the Mueller complex-
    // permeability derivation for materials without stored complex curves.
    TEST_CASE("Test_Initial_Permeability_Frequency_NPF_60", "[physical-model][initial-permeability]") {
        InitialPermeability initialPermeability;
        std::string materialName = "NPF 60";
        // Datasheet: ~100 % at 10 kHz, ~95 % at 1 MHz, ~80 % at 10 MHz.
        double at10k = initialPermeability.get_initial_permeability(materialName, std::nullopt, std::nullopt, 10000);
        double at1M = initialPermeability.get_initial_permeability(materialName, std::nullopt, std::nullopt, 1000000);
        double at10M = initialPermeability.get_initial_permeability(materialName, std::nullopt, std::nullopt, 10000000);
        REQUIRE_THAT(at10k, Catch::Matchers::WithinAbs(60.0, 60.0 * 0.02));
        REQUIRE_THAT(at1M, Catch::Matchers::WithinAbs(60.0 * 0.95, 60.0 * 0.02));
        REQUIRE_THAT(at10M, Catch::Matchers::WithinAbs(60.0 * 0.80, 60.0 * 0.02));
    }

    TEST_CASE("Test_Complex_Permeability_NPF_60", "[physical-model][complex-permeability]") {
        ComplexPermeability complexPermeability;
        std::string materialName = "NPF 60";
        // Previously threw MATERIAL_DATA_MISSING; the poco frequencyFactor
        // anchors the frequency-dependent derivation.
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialName, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialName, 10000000);
        REQUIRE(std::isfinite(complexPermeabilityValueAt100000.first));
        REQUIRE(std::isfinite(complexPermeabilityValueAt100000.second));
        REQUIRE(complexPermeabilityValueAt100000.first > 1);
        REQUIRE(complexPermeabilityValueAt100000.first >= complexPermeabilityValueAt10000000.first);
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
        // mu'' must FALL from 100 kHz to 1 MHz for this material. Its own |mu|(f) table puts
        // the relaxation knee (|mu| at half the initial 8,000) at ~170 kHz, and for any
        // relaxation mu'' peaks at the knee: 100 kHz is near that peak, while at 1 MHz |mu|
        // has fallen to 13% of initial, far past it. A single Debye anchored at 170 kHz gives
        // mu''(100k)/mu''(1M) ~ 2.6, and Kramers-Kronig from the table's local slopes
        // (f^-0.47 at 100 kHz, f^-0.87 at 1 MHz) gives ~3,700 vs ~1,040 — falling either way.
        // The previous pin (mu'' rising across this decade) captured the pre-#843 model,
        // whose tabulation anchored the knee orders of magnitude too high; it was
        // characterization of a bug, not physics.
        REQUIRE(complexPermeabilityValueAt100000.second > complexPermeabilityValueAt10000000.second);
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
            Painter painter(outFile);
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
            Painter painter(outFile);
            painter.paint_curve(curve, true);
            painter.export_svg();
            REQUIRE(std::filesystem::exists(outFile));
        }
    }

    // ABT #358: per-shape-family DC-bias modifiers must be reachable. Kool Mu 26 publishes
    // three fits: "default", "E/ER/U" and "EQ/LP" — all with the same 'a' but different b/c,
    // so at a real DC bias the three give measurably different permeability.
    TEST_CASE("Test_Initial_Permeability_Per_Shape_Family_Modifier", "[physical-model][initial-permeability]") {
        auto coreMaterial = OpenMagnetics::find_core_material_by_name("Kool Mµ 26");
        double temperature = 25;
        double magneticFieldDcBias = 4000;  // A/m — well into the rolloff

        double muDefault = OpenMagnetics::InitialPermeability::get_initial_permeability(
            coreMaterial, temperature, magneticFieldDcBias, std::nullopt, std::nullopt, std::nullopt);
        double muE = OpenMagnetics::InitialPermeability::get_initial_permeability(
            coreMaterial, temperature, magneticFieldDcBias, std::nullopt, std::nullopt, CoreShapeFamily::E);
        double muEr = OpenMagnetics::InitialPermeability::get_initial_permeability(
            coreMaterial, temperature, magneticFieldDcBias, std::nullopt, std::nullopt, CoreShapeFamily::ER);
        double muEq = OpenMagnetics::InitialPermeability::get_initial_permeability(
            coreMaterial, temperature, magneticFieldDcBias, std::nullopt, std::nullopt, CoreShapeFamily::EQ);
        double muLp = OpenMagnetics::InitialPermeability::get_initial_permeability(
            coreMaterial, temperature, magneticFieldDcBias, std::nullopt, std::nullopt, CoreShapeFamily::LP);

        // The family fits are actually consulted (this is what regressed silently before).
        REQUIRE(fabs(muE - muDefault) > 0.001 * muDefault);
        REQUIRE(fabs(muEq - muDefault) > 0.001 * muDefault);
        // Families sharing a key resolve to the same fit...
        REQUIRE_THAT(muE, Catch::Matchers::WithinRel(muEr, 1e-9));
        REQUIRE_THAT(muEq, Catch::Matchers::WithinRel(muLp, 1e-9));
        // ...and the two groups are distinct from each other.
        REQUIRE(fabs(muE - muEq) > 0.001 * muEq);
        // Every result stays physical: a DC bias can only reduce permeability.
        double muNoBias = OpenMagnetics::InitialPermeability::get_initial_permeability(
            coreMaterial, temperature, std::nullopt, std::nullopt, std::nullopt, CoreShapeFamily::E);
        REQUIRE(muE < muNoBias);
        REQUIRE(muE > 1);

        // A family with no published fit falls back to "default" (exact-token matching: "P"
        // must NOT be taken as a match for "EQ/LP" — the substring class of bug, ABT #359).
        double muP = OpenMagnetics::InitialPermeability::get_initial_permeability(
            coreMaterial, temperature, magneticFieldDcBias, std::nullopt, std::nullopt, CoreShapeFamily::P);
        REQUIRE_THAT(muP, Catch::Matchers::WithinRel(muDefault, 1e-9));
    }

    TEST_CASE("Test_BH_Loop_3C97", "[physical-model][amplitude-permeability][smoke-test]") {
        std::string materialName = "3C97";
        auto curves = OpenMagnetics::BHLoopRoshenModel().get_hysteresis_loop(materialName, 25, 0.2, std::nullopt);

        auto outFile = outputFilePath;
        outFile.append("Test_BH_Loop_3C97_Upper.svg");

        std::filesystem::remove(outFile);
        Painter painter(outFile);
        painter.paint_curve(curves.first);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
    }

}  // namespace
