#include "InitialPermeability.h"
#include "Utils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <vector>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;
#include <typeinfo>

SUITE(InitialPermeability) {
    TEST(Test_XFlux_60) {
        OpenMagnetics::InitialPermeability initial_permeability;
        std::string material_name = "XFlux 60";
        auto material_data = material_name;
        double initial_permeability_value = initial_permeability.get_initial_permeability(material_data);
        CHECK(initial_permeability_value == 60);
        double temperature = 89;
        double initial_permeability_value_with_temperature =
            initial_permeability.get_initial_permeability(material_data, &temperature, nullptr, nullptr);
        double expected = 60 * 1.0073;
        CHECK_CLOSE(initial_permeability_value_with_temperature, expected, 0.01 * expected);

        double magnetic_field_dc_bias = 3978.873577;
        double initial_permeability_value_with_magnetic_field_dc_bias =
            initial_permeability.get_initial_permeability(material_data, nullptr, &magnetic_field_dc_bias, nullptr);
        expected = 60 * 0.9601;
        CHECK_CLOSE(initial_permeability_value_with_magnetic_field_dc_bias, expected, 0.01 * expected);

        double frequency = 1250000;
        double initial_permeability_value_with_frequency =
            initial_permeability.get_initial_permeability(material_data, nullptr, nullptr, &frequency);
        expected = 60 * 0.968;
        CHECK_CLOSE(initial_permeability_value_with_frequency, expected, 0.01 * expected);
    }

    TEST(Test_N88) {
        OpenMagnetics::InitialPermeability initial_permeability;
        std::string material_name = "N88";
        auto material_data = OpenMagnetics::find_core_material_by_name(material_name);
        double initial_permeability_value = initial_permeability.get_initial_permeability(material_data);

        double expected = 1900;
        double manufacturer_tolerance = 0.25;
        CHECK_CLOSE(initial_permeability_value, expected, manufacturer_tolerance * expected);

        double temperature = 80;
        double initial_permeability_value_with_temperature =
            initial_permeability.get_initial_permeability(material_data, &temperature, nullptr, nullptr);
        expected = 3200;
        manufacturer_tolerance = 0.25;
        CHECK_CLOSE(initial_permeability_value_with_temperature, expected, manufacturer_tolerance * expected);

        temperature = 200;
        initial_permeability_value_with_temperature =
            initial_permeability.get_initial_permeability(material_data, &temperature, nullptr, nullptr);

        expected = 4500;
        manufacturer_tolerance = 0.25;
        CHECK_CLOSE(initial_permeability_value_with_temperature, expected, manufacturer_tolerance * expected);

        temperature = 300;
        initial_permeability_value_with_temperature =
            initial_permeability.get_initial_permeability(material_data, &temperature, nullptr, nullptr);

        expected = 1;
        manufacturer_tolerance = 0.25;
        CHECK_CLOSE(initial_permeability_value_with_temperature, expected, manufacturer_tolerance * expected);
    }

    TEST(Test_N30) {
        OpenMagnetics::InitialPermeability initial_permeability;
        std::string material_name = "N30";
        auto material_data = OpenMagnetics::find_core_material_by_name(material_name);
        double initial_permeability_value = initial_permeability.get_initial_permeability(material_data);

        double expected = 4300;
        double manufacturer_tolerance = 0.25;
        CHECK_CLOSE(initial_permeability_value, expected, manufacturer_tolerance * expected);

        double temperature = 80;
        double initial_permeability_value_with_temperature =
            initial_permeability.get_initial_permeability(material_data, &temperature, nullptr, nullptr);
        expected = 4300;
        manufacturer_tolerance = 0.25;
        CHECK_CLOSE(initial_permeability_value_with_temperature, expected, manufacturer_tolerance * expected);

        temperature = 200;
        initial_permeability_value_with_temperature =
            initial_permeability.get_initial_permeability(material_data, &temperature, nullptr, nullptr);

        expected = 1;
        manufacturer_tolerance = 0.25;
        CHECK_CLOSE(initial_permeability_value_with_temperature, expected, manufacturer_tolerance * expected);

        temperature = 300;
        initial_permeability_value_with_temperature =
            initial_permeability.get_initial_permeability(material_data, &temperature, nullptr, nullptr);

        expected = 1;
        manufacturer_tolerance = 0.25;
        CHECK_CLOSE(initial_permeability_value_with_temperature, expected, manufacturer_tolerance * expected);
    }
}