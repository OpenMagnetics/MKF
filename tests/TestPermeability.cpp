#include "RandomUtils.h"
#include <source_location>
#include "physical_models/InitialPermeability.h"
#include "physical_models/ComplexPermeability.h"
#include "support/Exceptions.h"
#include "physical_models/AmplitudePermeability.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <algorithm>
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
        // ABT #1739: 3C97's mu_i(T) table is Ferroxcube's own measurement, unrounded ("FXC Materials
        // Power.xlsx", 10 kHz, 0.1 mT); at the default 25 C it gives 3341.24. The old pin 3341.5 was
        // the same measurement rounded in the table.
        REQUIRE_THAT(initialPermeabilityValue, Catch::Matchers::WithinAbs(3341.24, 0.005));
        {
            double magneticFieldDcBias = 60;
            double temperature = 25;
            double initialPermeabilityValue = initialPermeability.get_initial_permeability(materialData, temperature, magneticFieldDcBias, std::nullopt);
            double expected = 2310;
            // Increased tolerance from 1% to 5% to account for model variations
            REQUIRE_THAT(initialPermeabilityValue, Catch::Matchers::WithinAbs(expected, 0.05 * expected));
        }
    }

    // ABT #1739: 3C99's mu_i(T) table is Ferroxcube's measurement, -40..180 C; its Curie point is 300 C.
    // Between 180 C and 300 C (and below -40 C) there is no data: MKF throws instead of running the
    // spline there. Above the Curie point the material is paramagnetic and mu_i = 1.
    TEST_CASE("Test_Initial_Permeability_Outside_Temperature_Table_Throws", "[physical-model][initial-permeability][initial-permeability-temperature-range]") {
        settings.reset();
        auto material = find_core_material_by_name("3C99");
        REQUIRE(material.get_curie_temperature());
        REQUIRE(material.get_curie_temperature().value() == 300);

        CHECK(InitialPermeability::get_initial_permeability_temperature_dependent(material, 100) > 1);
        CHECK_THROWS_AS(InitialPermeability::get_initial_permeability_temperature_dependent(material, 200), MaterialTemperatureOutOfRangeException);
        CHECK_THROWS_AS(InitialPermeability::get_initial_permeability_temperature_dependent(material, -60), MaterialTemperatureOutOfRangeException);
        try {
            InitialPermeability::get_initial_permeability_temperature_dependent(material, 200);
            FAIL("no exception at 200 C");
        }
        catch (const MaterialTemperatureOutOfRangeException& exception) {
            CHECK(exception.code() == ErrorCode::MATERIAL_TEMPERATURE_OUT_OF_RANGE);
            CHECK(static_cast<int>(exception.code()) == 306);
            CHECK(exception.table_minimum() == -40);
            CHECK(exception.table_maximum() == 180);
        }
        // Above the Curie point: paramagnetic, exactly 1, no throw.
        CHECK(InitialPermeability::get_initial_permeability_temperature_dependent(material, 350) == 1);

        // The same query through the public entry point throws too (temperature factor at 200 C).
        InitialPermeability initialPermeability;
        CHECK_THROWS_AS(initialPermeability.get_initial_permeability(material, 200.0, std::nullopt, std::nullopt), MaterialTemperatureOutOfRangeException);
    }

    TEST_CASE("Test_Initial_Permeability_Table_Starting_Above_Ambient_Evaluates_Inside_It", "[physical-model][initial-permeability][initial-permeability-temperature-range]") {
        // ABT #1739: DN100H's mu_i(T) table runs 26.46..139.82 C (Curie 130 C). A request at a
        // temperature inside the table is the table's own value; the ambient (25 C) reference that
        // only normalises frequency / DC-bias factors must not be evaluated for it, or the material
        // throws at every temperature.
        settings.reset();
        auto material = find_core_material_by_name("DN100H");
        InitialPermeability initialPermeability;
        double tabulated = InitialPermeability::get_initial_permeability_temperature_dependent(material, 60);
        REQUIRE(tabulated > 1);
        CHECK(initialPermeability.get_initial_permeability(material, 60.0, std::nullopt, std::nullopt) == tabulated);
        // Above both its Curie temperature (130 C) and the last measured point (139.82 C): paramagnetic, 1.
        CHECK(initialPermeability.get_initial_permeability(material, 140.0, std::nullopt, std::nullopt) == 1);
        // At 25 C there is no data: it still throws.
        CHECK_THROWS_AS(initialPermeability.get_initial_permeability(material, 25.0, std::nullopt, std::nullopt), MaterialTemperatureOutOfRangeException);
        CHECK_THROWS_AS(initialPermeability.get_initial_permeability(material, std::nullopt, std::nullopt, std::nullopt), MaterialTemperatureOutOfRangeException);
    }

    TEST_CASE("Test_Initial_Permeability_Measured_Table_Above_Minimum_Curie_Is_Evaluated", "[physical-model][initial-permeability][initial-permeability-temperature-range]") {
        // ABT #1739: Ferroxcube prints 3C90's Curie temperature as "TC >= 220 C", a guaranteed minimum, and
        // MAS stores 220 C; the measured mu_i(T) runs on to 240 C (225..240 C: 636, 91, 46, 34). mu_i = 1
        // starts only above both the Curie temperature and the last measured point.
        settings.reset();
        auto material = find_core_material_by_name("3C90");
        REQUIRE(material.get_curie_temperature());
        REQUIRE(material.get_curie_temperature().value() == 220);
        REQUIRE(InitialPermeability::get_paramagnetic_temperature(material).value() == 240);
        InitialPermeability initialPermeability;
        CHECK_THAT(InitialPermeability::get_initial_permeability_temperature_dependent(material, 230), Catch::Matchers::WithinAbs(91, 1e-9));
        CHECK_THAT(initialPermeability.get_initial_permeability(material, 230.0, std::nullopt, std::nullopt), Catch::Matchers::WithinAbs(91, 1e-9));
        CHECK(InitialPermeability::get_initial_permeability_temperature_dependent(material, 245) == 1);
        CHECK(initialPermeability.get_initial_permeability(material, 245.0, std::nullopt, std::nullopt) == 1);

        // 3C94: the same 220 C minimum, but its measured table ends at 180 C: no data from 180 to 220 C.
        auto material3C94 = find_core_material_by_name("3C94");
        REQUIRE(material3C94.get_curie_temperature().value() == 220);
        CHECK_THROWS_AS(InitialPermeability::get_initial_permeability_temperature_dependent(material3C94, 200), MaterialTemperatureOutOfRangeException);
        CHECK(InitialPermeability::get_initial_permeability_temperature_dependent(material3C94, 230) == 1);
    }

    TEST_CASE("Test_Initial_Permeability_Outside_Temperature_Table_Extrapolation_Flag", "[physical-model][initial-permeability][initial-permeability-temperature-range]") {
        settings.reset();
        settings.set_allow_material_data_extrapolation(true);
        auto material = find_core_material_by_name("3C99");
        Logger::getInstance().disableCollector();
        Logger::getInstance().enableCollector(LogLevel::WARNING);
        double value = InitialPermeability::get_initial_permeability_temperature_dependent(material, 200);
        auto collected = Logger::getInstance().drainCollected();
        Logger::getInstance().disableCollector();
        settings.reset();
        CHECK(value > 1);
        size_t warnings = 0;
        for (const auto& record : collected) {
            if (record.moduleOfOrigin == kMaterialDataExtrapolationModule && record.level == LogLevel::WARNING &&
                record.message.find("3C99") != std::string::npos && record.message.find("200.000000 C") != std::string::npos &&
                record.message.find("180.000000 C") != std::string::npos) {
                ++warnings;
            }
        }
        CHECK(warnings == 1);
    }

    TEST_CASE("Test_Initial_Permeability_Below_One_Below_Curie_Throws", "[physical-model][initial-permeability][initial-permeability-temperature-range]") {
        // A table whose own point says mu_i < 1 below the Curie point is unphysical: it is reported, not
        // floored to 1. (3C99's 180 C point replaced by 0.5; the rest of the record is untouched.)
        settings.reset();
        json materialJson;
        to_json(materialJson, find_core_material_by_name("3C99"));
        materialJson["name"] = "ABT1739 3C99 with mu 0.5 at 180 C";
        bool replaced = false;
        for (auto& point : materialJson["permeability"]["initial"]) {
            if (point.contains("temperature") && point["temperature"].get<double>() == 180) {
                point["value"] = 0.5;
                replaced = true;
            }
        }
        REQUIRE(replaced);
        CoreMaterial material;
        from_json(materialJson, material);
        CHECK(InitialPermeability::get_initial_permeability_temperature_dependent(material, 25) > 1);
        CHECK_THROWS_AS(InitialPermeability::get_initial_permeability_temperature_dependent(material, 180), MaterialException);
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
        // The high point is 10 MHz or the top of N22's complex-permeability table, whichever is
        // lower: above the table there is no data and get_complex_permeability throws.
        double highFrequency = std::min(10000000.0, complexPermeability.get_frequency_range(materialData).second);
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, highFrequency);
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
        // The high point is 10 MHz or the top of N49's complex-permeability table, whichever is
        // lower: above the table there is no data and get_complex_permeability throws.
        double highFrequency = std::min(10000000.0, complexPermeability.get_frequency_range(materialData).second);
        auto complexPermeabilityValueAt100000 = complexPermeability.get_complex_permeability(materialData, 100000);
        auto complexPermeabilityValueAt10000000 = complexPermeability.get_complex_permeability(materialData, highFrequency);
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

    // Outside its tabulated range a material's complex permeability used to hold the last point:
    // stale inline K081 data tabulated 1 Hz..1 MHz returned its 1 MHz point at 350 MHz and put a
    // common-mode choke's resonance there. It then threw (86c514b3). ABT #1861 (Alf approved
    // 2026-10-14): it is now ESTIMATED from the material's class neighbours
    // (ComplexPermeabilityEstimator) and marked as estimated, part by part; the tabulated-only getter
    // keeps throwing, naming the material, the frequency and the range. K081's table (MAS) runs from
    // 113335.906 Hz (its mu'' below that is a winding-resistance tail, removed) to 1 GHz; its mu' is
    // tabulated further down, so at 10 kHz only mu'' is estimated.
    TEST_CASE("Test_Complex_Permeability_Outside_Range_Estimated", "[physical-model][complex-permeability][complex-permeability-range]") {
        ComplexPermeability complexPermeability;
        auto k081 = find_core_material_by_name("K081");
        auto [minimumFrequency, maximumFrequency] = complexPermeability.get_frequency_range("K081");
        REQUIRE_THAT(minimumFrequency, Catch::Matchers::WithinRel(113335.906, 1e-9));
        REQUIRE_THAT(maximumFrequency, Catch::Matchers::WithinRel(1e9, 1e-9));

        for (double frequency : {minimumFrequency, maximumFrequency}) {
            auto value = complexPermeability.get_complex_permeability_with_source(k081, frequency);
            CHECK_FALSE(value.is_estimated());
            CHECK_NOTHROW(complexPermeability.get_tabulated_complex_permeability(k081, frequency));
        }

        auto below = complexPermeability.get_complex_permeability_with_source(k081, 10000);
        CHECK(below.is_estimated());
        CHECK(below.realSource == ComplexPermeabilitySource::TABULATED);
        CHECK(below.imaginarySource == ComplexPermeabilitySource::ESTIMATED_NEIGHBOURS);
        CHECK(std::isfinite(below.real));
        CHECK(below.imaginary > 0);
        CHECK_THAT(below.imaginaryMethod, Catch::Matchers::ContainsSubstring("FERRITE NI_ZN"));

        auto above = complexPermeability.get_complex_permeability_with_source(k081, maximumFrequency * 2);
        CHECK(above.realSource != ComplexPermeabilitySource::TABULATED);
        CHECK(above.imaginarySource != ComplexPermeabilitySource::TABULATED);
        CHECK(std::isfinite(above.real));
        CHECK(above.imaginary > 0);

        // The plain getter answers the same estimate.
        auto [real, imaginary] = complexPermeability.get_complex_permeability("K081", 10000);
        CHECK(real == below.real);
        CHECK(imaginary == below.imaginary);

        // The estimated stretches are labelled: below the mu'' table, and above both tables.
        auto ranges = complexPermeability.get_estimated_ranges(k081);
        CHECK(std::any_of(ranges.begin(), ranges.end(), [&](const ComplexPermeabilityEstimatedRange& range) {
            return range.part == ComplexPermeabilityPart::IMAGINARY && range.fromFrequency <= 10000 && range.toFrequency >= 10000;
        }));
        CHECK(std::any_of(ranges.begin(), ranges.end(), [&](const ComplexPermeabilityEstimatedRange& range) {
            return range.part == ComplexPermeabilityPart::REAL && range.fromFrequency <= 2 * maximumFrequency && range.toFrequency >= 2 * maximumFrequency;
        }));

        // Selection stays on data: the tabulated-only getter throws outside the range.
        CHECK_THROWS_AS(complexPermeability.get_tabulated_complex_permeability(k081, minimumFrequency / 2), ComplexPermeabilityFrequencyOutOfRangeException);
        CHECK_THROWS_AS(complexPermeability.get_tabulated_complex_permeability(k081, maximumFrequency * 2), ComplexPermeabilityFrequencyOutOfRangeException);
        try {
            complexPermeability.get_tabulated_complex_permeability(k081, 10000);
            FAIL("10 kHz is below K081's tabulated range and the tabulated-only getter must throw");
        }
        catch (const ComplexPermeabilityFrequencyOutOfRangeException& e) {
            CHECK(e.material_name() == "K081");
            CHECK(e.frequency() == 10000);
            CHECK(e.range_minimum() == minimumFrequency);
            CHECK(e.range_maximum() == maximumFrequency);
        }

        // A material whose complex permeability is derived from its initial-permeability table has
        // no MAS table to estimate from: it still throws outside the derived table's range, as
        // ComplexPermeabilityEstimateUnavailableException (a ComplexPermeabilityFrequencyOutOfRangeException).
        auto [derivedMinimum, derivedMaximum] = complexPermeability.get_frequency_range("XFlux 60");
        CHECK_THROWS_AS(complexPermeability.get_complex_permeability("XFlux 60", derivedMinimum / 2), ComplexPermeabilityEstimateUnavailableException);
        CHECK_THROWS_AS(complexPermeability.get_complex_permeability("XFlux 60", derivedMaximum * 2), ComplexPermeabilityEstimateUnavailableException);
    }

    // ABT #1861 guard: the approved estimator reproduces hidden data. 3C97 (Ferroxcube MnZn, tabulated
    // 100 kHz..30 MHz) with everything above 1 MHz removed: the estimate of |mu| at its own tabulated
    // frequencies above 1 MHz must be within the hold-out's class median for that pattern (MnZn,
    // nothing above 1 MHz: 19 %; 3C97 itself: 9 %). The truncated copy keeps the name, so 3C97 is
    // its own donor neither way; it is estimated directly, not through ComplexPermeability, which
    // caches per material name.
    // It must also beat holding mu' and mu'' at their last tabulated values: an estimate that only meets
    // the edge (D7) and then stays flat passes the 19 % bound on its first points but not this.
    TEST_CASE("Test_Complex_Permeability_Estimator_Holdout_3C97_Above_1MHz", "[physical-model][complex-permeability][abt-1861]") {
        auto full = find_core_material_by_name("3C97");
        auto complexData = full.get_permeability().get_complex().value();
        auto realPoints = std::get<std::vector<PermeabilityPoint>>(complexData.get_real());
        auto imaginaryPoints = std::get<std::vector<PermeabilityPoint>>(complexData.get_imaginary());
        const double cut = 1e6;
        std::vector<PermeabilityPoint> keptReal, keptImaginary;
        std::vector<double> hiddenFrequencies;
        for (const auto& point : realPoints) {
            if (point.get_frequency().value() <= cut) keptReal.push_back(point);
            else hiddenFrequencies.push_back(point.get_frequency().value());
        }
        for (const auto& point : imaginaryPoints) {
            if (point.get_frequency().value() <= cut) keptImaginary.push_back(point);
            else hiddenFrequencies.push_back(point.get_frequency().value());
        }
        REQUIRE(keptReal.size() >= 3);
        REQUIRE(keptImaginary.size() >= 3);
        REQUIRE(hiddenFrequencies.size() >= 10);
        auto truncated = full;
        ComplexPermeabilityData truncatedData;
        truncatedData.set_real(keptReal);
        truncatedData.set_imaginary(keptImaginary);
        auto permeability = truncated.get_permeability();
        permeability.set_complex(truncatedData);
        truncated.set_permeability(permeability);

        auto estimate = estimate_complex_permeability(truncated);
        const auto& realAbove = estimate.gap(ComplexPermeabilityPart::REAL, ComplexPermeabilityGapSide::ABOVE);
        const auto& imaginaryAbove = estimate.gap(ComplexPermeabilityPart::IMAGINARY, ComplexPermeabilityGapSide::ABOVE);
        REQUIRE(realAbove.refusal.empty());
        REQUIRE(imaginaryAbove.refusal.empty());
        for (const auto* gap : {&realAbove, &imaginaryAbove}) {
            REQUIRE(gap->donors.size() >= 3);
            for (const auto& donor : gap->donors) {
                CHECK(donor.name != "3C97");
            }
        }

        // The baseline the estimate has to beat: mu' and mu'' held at the truncated table's last points.
        ComplexPermeability tabulated;
        // (Read from the kept points: ComplexPermeability caches its tables per material name.)
        auto lastValue = [](const std::vector<PermeabilityPoint>& points) {
            double lastFrequency = 0, sum = 0;
            size_t count = 0;
            for (const auto& point : points) {
                double frequency = point.get_frequency().value();
                if (frequency > lastFrequency) { lastFrequency = frequency; sum = 0; count = 0; }
                if (frequency == lastFrequency) { sum += point.get_value(); ++count; }
            }
            return sum / static_cast<double>(count);
        };
        double edgeReal = lastValue(keptReal);
        double edgeImaginary = lastValue(keptImaginary);
        std::vector<double> errors, heldEdgeErrors;
        for (double frequency : hiddenFrequencies) {
            auto [trueReal, trueImaginary] = tabulated.get_tabulated_complex_permeability(full, frequency);
            double predictedReal = realAbove.evaluate(frequency).first;
            double predictedImaginary = imaginaryAbove.evaluate(frequency).first;
            errors.push_back(std::fabs(std::log10(std::hypot(predictedReal, predictedImaginary) / std::hypot(trueReal, trueImaginary))));
            heldEdgeErrors.push_back(std::fabs(std::log10(std::hypot(edgeReal, edgeImaginary) / std::hypot(trueReal, trueImaginary))));
        }
        auto median = [](std::vector<double> values) {
            std::sort(values.begin(), values.end());
            return values.size() % 2 ? values[values.size() / 2] : (values[values.size() / 2 - 1] + values[values.size() / 2]) / 2;
        };
        double medianError = median(errors);
        double heldEdgeMedianError = median(heldEdgeErrors);
        INFO("median |Z| error " << (std::pow(10.0, medianError) - 1) * 100 << " % over " << errors.size() << " hidden points; holding the "
             << "edge values: " << (std::pow(10.0, heldEdgeMedianError) - 1) * 100 << " %");
        CHECK(medianError <= std::log10(1.19));
        CHECK(medianError < heldEdgeMedianError);
        // D7 (Alf, 2026-10-14): the estimate meets the truncated table's own last value at the edge.
        CHECK_THAT(realAbove.evaluate(realAbove.edgeFrequency).first, Catch::Matchers::WithinRel(edgeReal, 1e-9));
        CHECK_THAT(imaginaryAbove.evaluate(imaginaryAbove.edgeFrequency).first, Catch::Matchers::WithinRel(edgeImaginary, 1e-9));
    }

    // ABT #1861 guard: the estimator covers FERRITE MnZn and NiZn only (Alf, 2026-10-14: powders, MgZn
    // and unclassified keep throwing). Mix 26 (iron powder, tabulated 10 kHz..70 MHz) and 46 (MgZn)
    // outside their tables throw ComplexPermeabilityEstimateUnavailableException, saying the class is
    // not covered.
    TEST_CASE("Test_Complex_Permeability_Estimator_Refuses_Powder_And_MgZn", "[physical-model][complex-permeability][abt-1861]") {
        ComplexPermeability complexPermeability;
        for (std::string name : {"Mix 26", "46"}) {
            INFO(name);
            auto material = find_core_material_by_name(name);
            REQUIRE(material.get_permeability().get_complex());
            auto [minimumFrequency, maximumFrequency] = complexPermeability.get_frequency_range(material);
            double outside = minimumFrequency > 1000 ? 1000 : maximumFrequency * 2;
            CHECK_NOTHROW(complexPermeability.get_complex_permeability(material, minimumFrequency));
            try {
                complexPermeability.get_complex_permeability(material, outside);
                FAIL(name + " outside its tabulated range must throw: the estimator does not cover its class");
            }
            catch (const ComplexPermeabilityEstimateUnavailableException& e) {
                CHECK(e.material_name() == name);
                CHECK(e.frequency() == outside);
                CHECK_THAT(e.reason(), Catch::Matchers::ContainsSubstring("is not covered"));
            }
            CHECK(complexPermeability.get_estimated_ranges(material).empty());
        }
    }

    // mu' is returned as tabulated, not floored at 1: past its ferromagnetic resonance K081's
    // table goes negative (-3.4 at 337 MHz), and that is data to report, not to clamp.
    TEST_CASE("Test_Complex_Permeability_Real_Part_Not_Floored", "[physical-model][complex-permeability][complex-permeability-range]") {
        ComplexPermeability complexPermeability;
        auto [realPart, imaginaryPart] = complexPermeability.get_complex_permeability("K081", 3.37e8);
        CHECK(realPart < 0);
        CHECK(imaginaryPart > 0);
    }

    // A07's mu' is tabulated (MAS Kramers-Kronig extension) to the top of its mu'', and it is
    // negative at 10 MHz (about -132). It is returned signed: flooring it at 1 is what this guards.
    TEST_CASE("Test_Complex_Permeability_A07_Negative_Real_Part_At_10_MHz", "[physical-model][complex-permeability][complex-permeability-sign]") {
        ComplexPermeability complexPermeability;
        auto [minimumFrequency, maximumFrequency] = complexPermeability.get_frequency_range("A07");
        REQUIRE(maximumFrequency > 1e7);
        auto [realPart, imaginaryPart] = complexPermeability.get_complex_permeability("A07", 1e7);
        CHECK(realPart < 0);
        CHECK(imaginaryPart > 0);
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

    // ABT #1790: XFlux 60's mu(f) is the Magnetics frequencyFactor quartic, a fit that reads 59.9 at
    // 100 kHz, 58.5 at 1 MHz and falls below the vacuum's 1 near 26.7 MHz (then goes negative). The
    // derived complex-permeability table ends where that model stops giving a material value: the
    // material answers inside the table and throws above its end. Reading the refused point with the
    // checked getter threw the whole table away, and XFlux 60 had no complex permeability at 100 kHz.
    TEST_CASE("Test_Complex_Permeability_Derived_Table_Ends_Where_Frequency_Model_Drops_Below_Vacuum", "[physical-model][complex-permeability][complex-permeability-range][abt-1790]") {
        ComplexPermeability complexPermeability;
        auto coreMaterial = OpenMagnetics::find_core_material_by_name("XFlux 60");

        // The premise: the material's own model is below the vacuum's permeability at 30 MHz and is
        // refused there, while it is a material value at 20 MHz.
        CHECK(OpenMagnetics::InitialPermeability::get_initial_permeability_model_value(coreMaterial, std::nullopt, std::nullopt, 3e7) < 1);
        CHECK_THROWS_AS(OpenMagnetics::InitialPermeability::get_initial_permeability(coreMaterial, std::nullopt, std::nullopt, 3e7), CalculationException);
        CHECK(OpenMagnetics::InitialPermeability::get_initial_permeability(coreMaterial, std::nullopt, std::nullopt, 2e7) > 1);

        // Inside the table: answers.
        auto [realAt100kHz, imaginaryAt100kHz] = complexPermeability.get_complex_permeability(coreMaterial, 1e5);
        auto [realAt1MHz, imaginaryAt1MHz] = complexPermeability.get_complex_permeability(coreMaterial, 1e6);
        CHECK(realAt100kHz > 50);
        CHECK(realAt100kHz < 60);
        CHECK(realAt1MHz > 50);
        CHECK(realAt1MHz < realAt100kHz);
        CHECK(imaginaryAt100kHz > 0);
        CHECK(imaginaryAt1MHz > imaginaryAt100kHz);

        // The table ends past 1 MHz and before the model's crossing of 1; its last point is a
        // frequency the checked initial-permeability getter accepts, and every tabulated |mu| is a
        // material value (not clamped to 1, not below it).
        auto [minimumFrequency, maximumFrequency] = complexPermeability.get_frequency_range(coreMaterial);
        CHECK(maximumFrequency > 1e7);
        CHECK(maximumFrequency < 2.67e7);
        CHECK(OpenMagnetics::InitialPermeability::get_initial_permeability(coreMaterial, std::nullopt, std::nullopt, maximumFrequency) > 1);
        auto derived = complexPermeability.calculate_complex_permeability_from_frequency_dependent_initial_permeability(coreMaterial);
        auto realPoints = std::get<std::vector<PermeabilityPoint>>(derived.get_real());
        auto imaginaryPoints = std::get<std::vector<PermeabilityPoint>>(derived.get_imaginary());
        REQUIRE(realPoints.size() == imaginaryPoints.size());
        REQUIRE(realPoints.size() > 10);
        for (size_t index = 0; index < realPoints.size(); ++index) {
            INFO("f = " << realPoints[index].get_frequency().value() << " Hz");
            CHECK(std::hypot(realPoints[index].get_value(), imaginaryPoints[index].get_value()) > 1);
        }

        // Above the table: throws, at its edge and where the model is below 1.
        CHECK_NOTHROW(complexPermeability.get_complex_permeability(coreMaterial, maximumFrequency));
        CHECK_THROWS_AS(complexPermeability.get_complex_permeability(coreMaterial, maximumFrequency * 1.001), ComplexPermeabilityFrequencyOutOfRangeException);
        CHECK_THROWS_AS(complexPermeability.get_complex_permeability(coreMaterial, 3e7), ComplexPermeabilityFrequencyOutOfRangeException);
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

    // ABT #838: the loss-temperature loop asks for the magnetisation curve at a new temperature on every pass. The
    // cache was keyed on the exact temperature and never evicted, so it grew for the life of the process. It is
    // bounded now, and a curve rebuilt after its eviction is the same numbers, bit for bit.
    TEST_CASE("Magnetisation curve cache stays bounded and an evicted curve rebuilds identically", "[physical-model][initial-permeability]") {
        auto coreMaterial = OpenMagnetics::find_core_material_by_name("3C97");
        double frequency = 100000;
        auto first = OpenMagnetics::InitialPermeability::get_magnetisation_curve(coreMaterial, 25.125, frequency);
        std::vector<double> firstFieldStrength = first->fieldStrength;
        std::vector<double> firstFluxDensity = first->fluxDensity;
        REQUIRE(firstFieldStrength.size() > 1000);

        const size_t capacity = OpenMagnetics::InitialPermeability::magnetisationCurveCacheCapacity;
        for (size_t index = 0; index < capacity + 20; ++index) {
            OpenMagnetics::InitialPermeability::get_magnetisation_curve(coreMaterial, 30 + 0.25 * index, frequency);
        }
        CHECK(OpenMagnetics::InitialPermeability::get_magnetisation_curve_cache_size() <= capacity);

        // The first curve was evicted; the caller's copy is still alive, and the rebuild matches it exactly.
        REQUIRE(first->fieldStrength == firstFieldStrength);
        auto rebuilt = OpenMagnetics::InitialPermeability::get_magnetisation_curve(coreMaterial, 25.125, frequency);
        CHECK(rebuilt.get() != first.get());
        CHECK(rebuilt->fieldStrength == firstFieldStrength);
        CHECK(rebuilt->fluxDensity == firstFluxDensity);
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


    TEST_CASE("Test_Initial_Permeability_Below_Vacuum_Is_Refused", "[physical-model][initial-permeability][abt-1740]") {
        // Past the knee the amplitude permeability of a ferrite (the slope of its Roshen loop) tends to 0, not
        // to the vacuum's 1: 3C96 gave mu_r = 0.499 at 6.2 kA/m and 8e-13 at 5e9 A/m, which the DC-bias fixed
        // point then used as a material property. A relative permeability below 1 is the model extrapolated
        // past saturation, and is refused.
        auto material = OpenMagnetics::Core::resolve_material("3C96");
        double modelValue = OpenMagnetics::InitialPermeability::get_initial_permeability_model_value(material, 25.0, 6.2e3, 100e3);
        INFO("3C96 model permeability at 6.2 kA/m: " << modelValue);
        REQUIRE(modelValue < 1);
        CHECK_THROWS_WITH(OpenMagnetics::InitialPermeability::get_initial_permeability(material, 25.0, 6.2e3, 100e3),
                          Catch::Matchers::ContainsSubstring("is saturated") && Catch::Matchers::ContainsSubstring("3C96"));
        // It is a SaturationException, the type an adviser rejects a candidate on (ABT #1740, #1534).
        CHECK_THROWS_AS(OpenMagnetics::InitialPermeability::get_initial_permeability(material, 25.0, 6.2e3, 100e3),
                        OpenMagnetics::SaturationException);

        // Without a bias, a value below 1 is not a saturation and must not be reported as one: the Magnetics
        // frequencyFactor quartic of Kool Mu 60 goes negative near 30 MHz. In the CoreAdviser's suppression
        // run these refusals once all read "is saturated", sending the reader after the wrong cause.
        auto koolMu = OpenMagnetics::Core::resolve_material("Kool Mµ 60");
        double highFrequencyModelValue = OpenMagnetics::InitialPermeability::get_initial_permeability_model_value(koolMu, 25.0, std::nullopt, 30e6);
        INFO("Kool Mu 60 model permeability at 30 MHz: " << highFrequencyModelValue);
        REQUIRE(highFrequencyModelValue < 1);
        CHECK_THROWS_WITH(OpenMagnetics::InitialPermeability::get_initial_permeability(koolMu, 25.0, std::nullopt, 30e6),
                          Catch::Matchers::ContainsSubstring("frequency or temperature dependence is outside its valid range") &&
                          !Catch::Matchers::ContainsSubstring("saturated"));
        // ... and it is not a SaturationException: a model out of its range must stop an adviser, not be
        // taken as a verdict on the candidate.
        bool refusedAsSaturation = false;
        bool refused = false;
        try {
            OpenMagnetics::InitialPermeability::get_initial_permeability(koolMu, 25.0, std::nullopt, 30e6);
        }
        catch (const OpenMagnetics::SaturationException&) {
            refusedAsSaturation = true;
        }
        catch (const OpenMagnetics::CalculationException&) {
            refused = true;
        }
        CHECK_FALSE(refusedAsSaturation);
        CHECK(refused);

        // Below the knee the value is the model's, unchanged.
        double biased = OpenMagnetics::InitialPermeability::get_initial_permeability(material, 25.0, 100.0, 100e3);
        CHECK(biased == OpenMagnetics::InitialPermeability::get_initial_permeability_model_value(material, 25.0, 100.0, 100e3));
        CHECK(biased > 1);

        // The magnetisation curve, which carries B on the vacuum slope past that point, is still built and
        // still inverts a flux density beyond saturation.
        double saturationFluxDensity = OpenMagnetics::Core::get_magnetic_flux_density_saturation(material, 25.0, false);
        double fieldStrength = 0;
        CHECK_NOTHROW(fieldStrength = OpenMagnetics::InitialPermeability::get_magnetic_field_dc_bias_for_flux_density(material, 1.2 * saturationFluxDensity, 25.0, 100e3, false));
        CHECK(fieldStrength > 0);

        // Sizing a gap for a bias above saturation is refused as a saturation too (ABT #1534 asked for a type
        // a caller can catch exactly, instead of INVALID_INPUT shared with malformed input).
        CHECK_THROWS_AS(OpenMagnetics::InitialPermeability::get_magnetic_field_dc_bias_for_flux_density(material, 1.2 * saturationFluxDensity, 25.0, 100e3, true),
                        OpenMagnetics::SaturationException);
    }

    TEST_CASE("Test_DC_Operating_Point_On_A_Held_Magnetisation_Curve_Is_The_Material_One", "[physical-model][initial-permeability][abt-1764]") {
        // The per-cycle line bias of a line-cycle record (ABT #1764) solves one load line per switching cycle,
        // ~2000 per record. Looking the material's curve up for each of them (a cache key built from the
        // material's JSON on every call) cost 14 us a cycle, 1.0 s of an 8.5 s DMC CoreAdviser run; the record
        // now fetches the curve once and solves on it. Both entries must give the same operating point, to the
        // bit, and that point must sit on the load line N*I = H*l_e + B*A_e*R_gap.
        auto material = OpenMagnetics::Core::resolve_material("Kool Mµ 60");
        double temperature = 25.0;
        double frequency = 50.0;
        double effectiveLength = 0.1;
        double effectiveArea = 1e-4;
        // The caller's share of the cached curve (the cache is bounded and may evict it, ABT #838).
        const auto curveShare = OpenMagnetics::InitialPermeability::get_magnetisation_curve(material, temperature, frequency);
        const auto& curve = *curveShare;
        for (double gapReluctance : {0.0, 2e6}) {
            for (double magnetomotiveForce : {0.0, 1e-3, 0.5, 7.3, 144.0, 238.1, 1500.0, 4e4, 1e9}) {
                auto fromMaterial = OpenMagnetics::InitialPermeability::get_dc_operating_point_for_magnetomotive_force(
                    material, magnetomotiveForce, effectiveLength, effectiveArea, gapReluctance, temperature, frequency);
                auto fromCurve = OpenMagnetics::InitialPermeability::get_dc_operating_point_for_magnetomotive_force(
                    curve, magnetomotiveForce, effectiveLength, effectiveArea, gapReluctance);
                INFO("MMF " << magnetomotiveForce << " A, R_gap " << gapReluctance << " 1/H: H " << fromCurve.first << " A/m, B " << fromCurve.second << " T");
                CHECK(fromCurve.first == fromMaterial.first);
                CHECK(fromCurve.second == fromMaterial.second);
                double loadLine = fromCurve.first * effectiveLength + fromCurve.second * effectiveArea * gapReluctance;
                CHECK_THAT(loadLine, Catch::Matchers::WithinRel(magnetomotiveForce, 1e-9));
            }
        }
        CHECK_THROWS_AS(OpenMagnetics::InitialPermeability::get_dc_operating_point_for_magnetomotive_force(curve, -1.0, effectiveLength, effectiveArea, 0.0),
                        OpenMagnetics::InvalidInputException);
    }

    // ABT #1736: remanence and coercive force are tabulated at 25 and 100 C for most ferrites, and the
    // BH-loop permeability model is asked up to the Curie temperature. Extrapolating the table's slope
    // made them cross zero below Curie (N87: Br at 149.8 C, Hc at 207 C, Curie 210 C; Fair-Rite 97: Br
    // at 146 C, Hc at 175 C, Curie 220 C), where the loop is singular: N87 under ANY bias gave mu = 44
    // at 150 C (3896 unbiased), and 97 gave mu = 0 at exactly 175 C, which the magnetizing model then
    // reported as a saturated core at 0.108 T. They now go linearly from the last data point to zero AT
    // the Curie temperature. The loop's saturation point is not anchored: it stays the saturation
    // data the saturation gates read (Alf, 2026-10-12).
    TEST_CASE("BH-loop permeability has no singular temperature between the remanence data and Curie", "[physical-model][initial-permeability][amplitude-permeability][abt-1736]") {
        struct Case { std::string material; double temperature; };
        for (auto& [materialName, temperature] : std::vector<Case>{{"N87", 150.0}, {"N87", 149.8}, {"N87", 207.0}, {"97", 175.0}, {"97", 146.0}}) {
            INFO(materialName << " at " << temperature << " C");
            auto material = OpenMagnetics::Core::resolve_material(materialName);
            double curieTemperature = OpenMagnetics::Core::get_curie_temperature(material);
            REQUIRE(temperature < curieTemperature);

            // Below Curie the material is ferromagnetic: it has a remanence and a coercive force.
            CHECK(OpenMagnetics::Core::get_remanence(material, temperature) > 0);
            CHECK(OpenMagnetics::Core::get_coercive_force(material, temperature) > 0);

            // A bias far below the knee leaves the permeability where it is without bias, here and
            // a fraction of a degree to either side.
            for (double probedTemperature : {temperature - 0.2, temperature, temperature + 0.2}) {
                INFO("probed at " << probedTemperature << " C");
                double unbiased = OpenMagnetics::InitialPermeability::get_initial_permeability_model_value(material, probedTemperature, std::nullopt, 100e3);
                for (double magneticFieldDcBias : {0.078, 1.0}) {
                    INFO("H = " << magneticFieldDcBias << " A/m");
                    double biased = 0;
                    REQUIRE_NOTHROW(biased = OpenMagnetics::InitialPermeability::get_initial_permeability(material, probedTemperature, magneticFieldDcBias, 100e3));
                    CHECK_THAT(biased, Catch::Matchers::WithinRel(unbiased, 0.05));
                }
            }
        }

        // The anchor itself: zero at the Curie temperature and above, the data inside the data.
        auto n87 = OpenMagnetics::Core::resolve_material("N87");
        double n87Curie = OpenMagnetics::Core::get_curie_temperature(n87);
        CHECK(OpenMagnetics::Core::get_remanence(n87, n87Curie) == 0);
        CHECK(OpenMagnetics::Core::get_coercive_force(n87, n87Curie) == 0);
        CHECK(OpenMagnetics::Core::get_remanence(n87, n87Curie + 20) == 0);
        // No loop at or above Curie: refused, by name.
        CHECK_THROWS_WITH(OpenMagnetics::BHLoopRoshenModel().get_hysteresis_loop(n87, n87Curie, std::nullopt, 10.0),
                          Catch::Matchers::ContainsSubstring("not ferromagnetic at or above its Curie temperature"));

        // A material without a Curie temperature cannot be evaluated past its data: no guess.
        auto withoutCurie = n87;
        withoutCurie.set_curie_temperature(std::nullopt);
        CHECK_THROWS_WITH(OpenMagnetics::Core::get_remanence(withoutCurie, 150.0), Catch::Matchers::ContainsSubstring("has no Curie temperature"));
        CHECK_THROWS_WITH(OpenMagnetics::Core::get_coercive_force(withoutCurie, 150.0), Catch::Matchers::ContainsSubstring("has no Curie temperature"));
        CHECK_NOTHROW(OpenMagnetics::Core::get_remanence(withoutCurie, 60.0));
    }

}  // namespace
