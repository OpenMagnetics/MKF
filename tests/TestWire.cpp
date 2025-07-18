#include "constructive_models/Wire.h"
#include "support/Utils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

SUITE(Wire) {
    auto masPath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("MAS/").string();
    double max_error = 0.05;

    TEST(Test_Sample_Wire) {
        auto wireFilePath = masPath + "samples/magnetic/wire/round/0.000140.json";
        std::ifstream json_file(wireFilePath);
        auto wireJson = json::parse(json_file);

        OpenMagnetics::Wire wire(wireJson);

        auto conductingDiameter = wire.get_conducting_diameter().value().get_nominal().value();

        CHECK(conductingDiameter == wireJson["conductingDiameter"]["nominal"]);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_1) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05);
        double expectedValue = 0.755;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_1) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.1e-05);
        double expectedValue = 0.64;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_1) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048);
        double expectedValue = 0.87;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_2) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 2);
        double expectedValue = 0.616;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_2) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.1e-05, 2);
        double expectedValue = 0.455;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_2) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 2);
        double expectedValue = 0.8;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_3) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 3);
        double expectedValue = 0.523;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_3) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.1e-05, 3);
        double expectedValue = 0.334;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_3) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 3);
        double expectedValue = 0.741;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_1_Nema) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.79;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_1_Nema) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.3e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.71;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_1_Nema) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 1, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.89;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_2_Nema) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 2, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.65;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_2_Nema) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.3e-05, 2, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.52;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_2_Nema) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 2, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.81;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_3_Nema) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 3, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.55;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_3_Nema) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(4e-05, 3, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.51;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_3_Nema) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 3, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.74;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Thin_Round_Insulated_Wire_1_Layer_Thin_Layer_Thickness) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.000101, 1, 3.81e-05, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.321961;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Filling_Factors_Thick_Round_Insulated_Wire_1_Layer_Thick_Layer_Thickness) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00129, 1, 7.62e-05, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.799184;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Outer_Diameter_Thick_Round_Insulated_Wire_1_Layer_Thick_Layer_Thickness) {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_round(0.00129, 1, 7.62e-05, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00144;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Filling_Factors_Thick_Round_Insulated_Wire_3_Layer_Thick_Layer_Thickness) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00129, 3, 7.62e-05, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.5446;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Outer_Diameter_Thick_Round_Insulated_Wire_3_Layers_Thick_Layer_Thickness) {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_round(0.00129, 3, 7.62e-05, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.001748;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Filling_Factor_Thick_Litz_Wire_Served_1_Layer_Few_Strands) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_served_litz(0.0001, 66, 1, 1, WireStandard::IEC_60317, false);
        double expectedValue = 0.458122;

        CHECK_CLOSE(expectedValue, fillingFactor, max_error * expectedValue);
    }

    TEST(Test_Outer_Diameter_Litz_Wire_Grade_1_Served_1_Layer_Few_Strands) {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_served_litz(2e-05, 10, 1, 1, WireStandard::IEC_60317);
        double expectedMinimumValue = 0.000112;
        double expectedMaximumValue = 0.000142;

        CHECK(outerDiameter > expectedMinimumValue);
        CHECK(outerDiameter < expectedMaximumValue);
    }

    TEST(Test_Outer_Diameter_Litz_Wire_Grade_2_Served_1_Layer_1_Few_Strands) {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_served_litz(2e-05, 10, 2, 1, WireStandard::IEC_60317);
        double expectedOuterDiameter = 0.000153;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Outer_Diameter_Thick_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Few_Strands) {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_insulated_litz(0.000102, 66, 3, 7.62e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00152908;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Outer_Diameter_Thin_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Many_Strands) {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_insulated_litz(5.1e-05, 825, 3, 7.62e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00253238;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Outer_Diameter_Thin_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Many_Strands_Diameter_Not_In_Db) {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_insulated_litz(5.42e-05, 825, 3, 7.62e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00253238;

        CHECK(expectedOuterDiameter < outerDiameter);
    }

    TEST(Test_Filling_Factor_Litz_Wire_Grade_2_Served_1_Layer_1_Few_Strands) {
        auto outerDiameter = OpenMagnetics::Wire::get_filling_factor_served_litz(2e-05, 10, 2, 1, WireStandard::IEC_60317);
        double expectedOuterDiameter = 0.235;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Filling_Factor_Thick_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Few_Strands) {
        auto outerDiameter = OpenMagnetics::Wire::get_filling_factor_insulated_litz(0.000102, 66, 3, 7.62e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.3449;

        CHECK_CLOSE(expectedOuterDiameter, outerDiameter, max_error * expectedOuterDiameter);
    }

    TEST(Test_Outer_Width_Small_Rectangular_Grade_1) {
        auto outerWidth = OpenMagnetics::Wire::get_outer_width_rectangular(0.002, 1, WireStandard::IEC_60317);
        double expectedOuterWidth = 0.00206;

        CHECK_CLOSE(expectedOuterWidth, outerWidth, max_error * expectedOuterWidth);
    }

    TEST(Test_Outer_Width_Small_Rectangular_Grade_2) {
        auto outerWidth = OpenMagnetics::Wire::get_outer_width_rectangular(0.002, 2, WireStandard::IEC_60317);
        double expectedOuterWidth = 0.00217;

        CHECK_CLOSE(expectedOuterWidth, outerWidth, max_error * expectedOuterWidth);
    }

    TEST(Test_Outer_Width_Large_Rectangular_Grade_1) {
        auto outerWidth = OpenMagnetics::Wire::get_outer_width_rectangular(0.016, 1, WireStandard::IEC_60317);
        double expectedOuterWidth = 0.01608;

        CHECK_CLOSE(expectedOuterWidth, outerWidth, max_error * expectedOuterWidth);
    }

    TEST(Test_Outer_Width_Large_Rectangular_Grade_2) {
        auto outerWidth = OpenMagnetics::Wire::get_outer_width_rectangular(0.016, 2, WireStandard::IEC_60317);
        double expectedOuterWidth = 0.01614;

        CHECK_CLOSE(expectedOuterWidth, outerWidth, max_error * expectedOuterWidth);
    }

    TEST(Test_Outer_Height_Small_Rectangular_Grade_1) {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(0.0008, 1, WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00088;

        CHECK_CLOSE(expectedOuterHeight, outerHeight, max_error * expectedOuterHeight);
    }

    TEST(Test_Outer_Height_Small_Rectangular_Grade_2) {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(0.0008, 2, WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00092;

        CHECK_CLOSE(expectedOuterHeight, outerHeight, max_error * expectedOuterHeight);
    }

    TEST(Test_Outer_Height_Large_Rectangular_Grade_1) {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(0.0045, 1, WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00456;

        CHECK_CLOSE(expectedOuterHeight, outerHeight, max_error * expectedOuterHeight);
    }

    TEST(Test_Outer_Height_Large_Rectangular_Grade_2) {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(0.0045, 2, WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00467;

        CHECK_CLOSE(expectedOuterHeight, outerHeight, max_error * expectedOuterHeight);
    }

    TEST(Test_Filling_Factor_Small_Rectangular_Grade_2) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_rectangular(0.002, 0.0008, 2, WireStandard::IEC_60317);
        double expectedFillingFactor = 0.720267;

        CHECK_CLOSE(expectedFillingFactor, fillingFactor, max_error * expectedFillingFactor);
    }

    TEST(Test_Filling_Factor_Large_Rectangular_Grade_2) {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_rectangular(0.016, 0.0045, 2, WireStandard::IEC_60317);
        double expectedFillingFactor = 0.948615;

        CHECK_CLOSE(expectedFillingFactor, fillingFactor, max_error * expectedFillingFactor);
    }

    TEST(Test_Conducting_Area_Small_Rectangular) {
        auto conductingArea = OpenMagnetics::Wire::get_conducting_area_rectangular(0.002, 0.0008, WireStandard::IEC_60317);
        double expectedConductingArea = 0.00000146;

        CHECK_CLOSE(expectedConductingArea, conductingArea, max_error * expectedConductingArea);
    }

    TEST(Test_Conducting_Area_Large_Rectangular) {
        auto conductingArea = OpenMagnetics::Wire::get_conducting_area_rectangular(0.016, 0.0045, WireStandard::IEC_60317);
        double expectedConductingArea = 0.00007114;

        CHECK_CLOSE(expectedConductingArea, conductingArea, max_error * expectedConductingArea);
    }

    TEST(Test_Outer_Height_Tiny_Rectangular_Grade_2) {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(1e-9, 2, WireStandard::IEC_60317);
        double expectedOuterHeight = 1.2e-9;

        CHECK_CLOSE(expectedOuterHeight, outerHeight, max_error * expectedOuterHeight);
    }

    TEST(Test_Outer_Diameter_Litz_Wire_Unserved_Medium_Strands) {
        auto diameter = OpenMagnetics::Wire::get_outer_diameter_served_litz(0.000071, 270, 1, 0, WireStandard::IEC_60317);
        double expectedMaximumValue = 0.001767;
        double expectedMinimumValue = 0.001641;

        CHECK(diameter > expectedMinimumValue);
        CHECK(diameter < expectedMaximumValue);
    }

    TEST(Test_Outer_Diameter_Litz_Wire_Served_1_Layer_Medium_Strands) {
        auto diameter = OpenMagnetics::Wire::get_outer_diameter_served_litz(0.000071, 270, 1, 1, WireStandard::IEC_60317);
        double expectedMaximumValue = 0.001807;
        double expectedMinimumValue = 0.001666;

        CHECK(diameter > expectedMinimumValue);
        CHECK(diameter < expectedMaximumValue);
    }
}

SUITE(Wire_Effective_Current_Density) {
    auto label = WaveformLabel::TRIANGULAR;
    double offset = 0;
    double peakToPeak = 2 * 1.73205;
    double dutyCycle = 0.5;
    double magnetizingInductance = 1e-3;
    double max_error = 0.01;
    double temperature = 22;

    TEST(Test_Effective_Current_Density_Medium_Frequency_Round_Operation_Point) {
        double frequency = 100000;
        Processed processed;
        processed.set_label(label);
        processed.set_offset(offset);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(inputs.get_primary_excitation(), temperature);
        double expectedEffectiveCurrentDensity = 5.33e6;

        CHECK_CLOSE(expectedEffectiveCurrentDensity, effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity);
    }

    TEST(Test_Effective_Current_Density_Low_Frequency_Round) {
        double frequency = 10;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 5.093e6;

        CHECK_CLOSE(expectedEffectiveCurrentDensity, effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity);
    }

    TEST(Test_Effective_Current_Density_Medium_Frequency_Round) {
        double frequency = 100000;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 5.283e6;

        CHECK_CLOSE(expectedEffectiveCurrentDensity, effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity);
    }

    TEST(Test_Effective_Current_Density_High_Frequency_Round) {
        double frequency = 1000000;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 11.19e6;

        CHECK_CLOSE(expectedEffectiveCurrentDensity, effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity);
    }

    TEST(Test_Effective_Current_Density_Low_Frequency_Litz) {
        double frequency = 10;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 5.093e5;

        CHECK_CLOSE(expectedEffectiveCurrentDensity, effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity);
    }

    TEST(Test_Effective_Current_Density_High_Frequency_Litz) {
        double frequency = 10000000;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 5.24e5;

        CHECK_CLOSE(expectedEffectiveCurrentDensity, effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity);
    }

    TEST(Test_Effective_Current_Density_Low_Frequency_Rectangular) {
        double frequency = 10;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 3.96e5;

        CHECK_CLOSE(expectedEffectiveCurrentDensity, effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity);
    }

    TEST(Test_Effective_Current_Density_High_Frequency_Rectangular) {
        double frequency = 1000000;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 2.09e6;

        CHECK_CLOSE(expectedEffectiveCurrentDensity, effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity);
    }

    TEST(Test_Conducting_Area_Large_Rectangular) {
        double frequency = 10;

        Processed processed;
        processed.set_label(label);
        processed.set_offset(offset);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              label,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              offset);
        
        auto wire = find_wire_by_name("Round 0.5 - Grade 1");
        // double expectedNumberParallels = 69;

        // CHECK_CLOSE(expectedConductingArea, conductingArea, max_error * expectedConductingArea);
    }
}

SUITE(Wire_Number_Parallels_Calculation) {
    auto label = WaveformLabel::TRIANGULAR;
    double offset = 0;
    double peakToPeak = 2 * 1.73205;
    double dutyCycle = 0.5;
    double magnetizingInductance = 1e-3;
    double max_error = 0.01;
    double temperature = 22;

    OpenMagnetics::Inputs setup_inputs(double frequency) {

        Processed processed;
        processed.set_label(label);
        processed.set_offset(offset);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        return OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                       magnetizingInductance,
                                                                                       temperature,
                                                                                       label,
                                                                                       peakToPeak,
                                                                                       dutyCycle,
                                                                                       offset);
    }

    TEST(Test_Number_Parallels_Low_Frequency_Round_1_Parallel) {
        double frequency = 10;

        auto inputs = setup_inputs(frequency);
        
        auto wire = find_wire_by_name("Round 0.5 - Grade 1");
        double maximumEffectiveCurrentDensity = 5.5e6;
        double expectedNumberParallels = 1;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        CHECK_EQUAL(expectedNumberParallels, numberParallels);
    }

    TEST(Test_Number_Parallels_Low_Frequency_Round_2_Parallels) {
        double frequency = 10;

        auto inputs = setup_inputs(frequency);
        
        auto wire = find_wire_by_name("Round 0.5 - Grade 1");
        double maximumEffectiveCurrentDensity = 5e6;
        double expectedNumberParallels = 2;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        CHECK_EQUAL(expectedNumberParallels, numberParallels);
    }

    TEST(Test_Number_Parallels_High_Frequency_Round_3_Parallels) {
        double frequency = 1000000;

        auto inputs = setup_inputs(frequency);
        
        auto wire = find_wire_by_name("Round 0.5 - Grade 1");
        double maximumEffectiveCurrentDensity = 5e6;
        double expectedNumberParallels = 3;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        CHECK_EQUAL(expectedNumberParallels, numberParallels);
    }

    TEST(Test_Number_Parallels_High_Frequency_Litz_2_Parallels) {
        double frequency = 1000000;

        auto inputs = setup_inputs(frequency);
        
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        double maximumEffectiveCurrentDensity = 5e5;
        double expectedNumberParallels = 2;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        CHECK_EQUAL(expectedNumberParallels, numberParallels);
    }

    TEST(Test_Number_Parallels_Low_Frequency_Rectangular_1_Parallels) {
        double frequency = 10;

        auto inputs = setup_inputs(frequency);
        
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        double maximumEffectiveCurrentDensity = 5e6;
        double expectedNumberParallels = 1;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        CHECK_EQUAL(expectedNumberParallels, numberParallels);
    }

    TEST(Test_Number_Parallels_High_Frequency_Rectangular_4_Parallels) {
        double frequency = 1000000;

        auto inputs = setup_inputs(frequency);
        
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        double maximumEffectiveCurrentDensity = 1e6;
        double expectedNumberParallels = 3;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        CHECK_EQUAL(expectedNumberParallels, numberParallels);
    }
}

SUITE(Wire_Coating) {
    double max_error = 0.05;
    TEST(Test_Coating_Label_Uniqueness) {
        auto wires = get_wires();
        std::vector<std::string> coatingLabels;
        for (auto wire : wires) {
            auto coatingLabel = wire.encode_coating_label();
            coatingLabels.push_back(coatingLabel);
        }

        CHECK(std::find(coatingLabels.begin(), coatingLabels.end(), "Bare") != coatingLabels.end());
    }

    TEST(Test_Coating_Decoding) {
        auto wires = get_wires();
        for (auto wire : wires) {
            auto coatingLabel = wire.encode_coating_label();
            auto coating = wire.resolve_coating();
            if (coating) {
                auto decodedCoating = OpenMagnetics::Wire::decode_coating_label(coatingLabel);
                CHECK(coating->get_type().value() == decodedCoating->get_type().value());
                if (coating->get_number_layers()) {
                    CHECK(coating->get_number_layers().value() == decodedCoating->get_number_layers().value());
                }
                if (coating->get_temperature_rating()) {
                    CHECK(coating->get_temperature_rating().value() == decodedCoating->get_temperature_rating().value());
                }
                if (coating->get_breakdown_voltage() && coating->get_type().value() == InsulationWireCoatingType::INSULATED) {
                    CHECK(coating->get_breakdown_voltage().value() == decodedCoating->get_breakdown_voltage().value());
                }
                if (coating->get_grade()) {
                    CHECK(coating->get_grade().value() == decodedCoating->get_grade().value());
                }
            }
        }
    }


    TEST(Test_Coating_Relative_Permittivity) {
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Round 0.80 - Grade 1"));
        auto relativePermittivity = wire.get_coating_relative_permittivity();
        CHECK_CLOSE(3.7, relativePermittivity, max_error * 3.7);
    }

    TEST(Test_Coating_Relative_Permittivity_Web0) {
        auto wire = OpenMagnetics::Wire(json::parse(R"({"type": "round", "conductingDiameter": {"nominal": 0.001}, "material": "copper", "outerDiameter": {"nominal": 0.001062}, "coating": {"breakdownVoltage": 2700.0, "grade": 1, "type": "enamelled"}, "name": "Round 1.00 - Grade 1", "numberConductors": 1, "standard": "IEC 60317", "standardName": "1.00 mm"})"));
        auto relativePermittivity = wire.get_coating_relative_permittivity();
        CHECK_CLOSE(3.7, relativePermittivity, max_error * 3.7);
    }

    TEST(Test_Coating_Material) {
        auto material = find_insulation_material_by_name("ETFE");
        json mierda;
        to_json(mierda, material);
        std::cout << mierda << std::endl;
    }

}

SUITE(Wire_Equivalents) {
    double max_error = 0.05;

    TEST(Test_Find_Round_By_Dimension_European) {
        auto wire = find_wire_by_dimension(0.00072, WireType::ROUND, WireStandard::IEC_60317, false);
        CHECK_EQUAL(wire.get_standard_name().value(), "0.71 mm");
    }

    TEST(Test_Find_Round_By_Dimension_American) {
        auto wire = find_wire_by_dimension(0.00072, WireType::ROUND, WireStandard::NEMA_MW_1000_C);
        CHECK_EQUAL(wire.get_standard_name().value(), "21 AWG");
    }

    TEST(Test_Find_Among_All_By_Dimension) {
        auto wire = find_wire_by_dimension(0.00072);
        CHECK_EQUAL(wire.get_standard_name().value(), "21 AWG");
    }

    TEST(Test_Find_Rectangular_By_Dimension) {
        auto wire = find_wire_by_dimension(0.00072, WireType::RECTANGULAR);
        auto conductingHeight = resolve_dimensional_values(wire.get_conducting_height().value());
        CHECK_EQUAL(conductingHeight, 0.0008);
    }

    TEST(Test_Find_Foil_By_Dimension) {
        auto wire = find_wire_by_dimension(0.00072, WireType::FOIL);
        auto conductingWidth = resolve_dimensional_values(wire.get_conducting_width().value());
        CHECK_EQUAL(conductingWidth, 0.0007);
    }

    TEST(Test_Litz_To_Litz_Equivalent) {
        double effectiveFrequency = 1234981;
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::LITZ, effectiveFrequency);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::LITZ);
        CHECK_EQUAL(numberConductors, 1000);
        CHECK_CLOSE(conductingDimension, 0.00005, max_error * 0.00005);
    }

    TEST(Test_Round_To_Litz_Equivalent) {
        double effectiveFrequency = 1234981;
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::LITZ, effectiveFrequency);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::LITZ);
        CHECK_EQUAL(numberConductors, 71);
        CHECK_CLOSE(conductingDimension, 0.00006, max_error * 0.00006);
    }

    TEST(Test_Rectangular_To_Litz_Equivalent) {
        double effectiveFrequency = 1234981;
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::LITZ, effectiveFrequency);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::LITZ);
        CHECK_EQUAL(numberConductors, 914);
        CHECK_CLOSE(conductingDimension, 0.00006, max_error * 0.00006);
    }

    TEST(Test_Foil_To_Litz_Equivalent) {
        double effectiveFrequency = 1234981;
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Foil 0.2"));
        oldWire.set_nominal_value_conducting_height(0.01);

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::LITZ, effectiveFrequency);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::LITZ);
        CHECK_EQUAL(numberConductors, 725);
        CHECK_CLOSE(conductingDimension, 0.00006, max_error * 0.00006);
    }

    TEST(Test_Litz_To_Round_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::ROUND);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::ROUND);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0016, max_error * 0.0016);
    }

    TEST(Test_Round_To_Round_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::ROUND);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::ROUND);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0005, max_error * 0.0005);
    }

    TEST(Test_Rectangular_To_Round_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::ROUND);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::ROUND);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0009, max_error * 0.0009);
    }

    TEST(Test_Foil_To_Round_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Foil 0.2"));
        oldWire.set_nominal_value_conducting_height(0.001);

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::ROUND);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::ROUND);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0002, max_error * 0.0002);
    }

    TEST(Test_Litz_To_Rectangular_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::RECTANGULAR);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::RECTANGULAR);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0016, max_error * 0.0016);
    }

    TEST(Test_Round_To_Rectangular_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Round 0.80 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::RECTANGULAR);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::RECTANGULAR);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0008, max_error * 0.0008);
    }

    TEST(Test_Rectangular_To_Rectangular_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::RECTANGULAR);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::RECTANGULAR);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.00085, max_error * 0.00085);
    }

    TEST(Test_Foil_To_Rectangular_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Foil 0.2"));
        oldWire.set_nominal_value_conducting_height(0.001);

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::RECTANGULAR);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::RECTANGULAR);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0008, max_error * 0.0008);
    }

    TEST(Test_Litz_To_Foil_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::FOIL);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::FOIL);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0016, max_error * 0.0016);
    }

    TEST(Test_Round_To_Foil_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Round 0.80 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::FOIL);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::FOIL);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0008, max_error * 0.0008);
    }

    TEST(Test_Rectangular_To_Foil_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::FOIL);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::FOIL);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0008, max_error * 0.0008);
    }

    TEST(Test_Foil_To_Foil_Equivalent) {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Foil 0.2"));
        oldWire.set_nominal_value_conducting_height(0.001);

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::FOIL);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        CHECK(newWire.get_type() == WireType::FOIL);
        CHECK_EQUAL(numberConductors, 1);
        CHECK_CLOSE(conductingDimension, 0.0002, max_error * 0.0002);
    }

    TEST(Test_Equivalent_Web_0) {
        auto oldWire = OpenMagnetics::Wire(json::parse(R"({"coating":{"breakdownVoltage":13500,"grade":null,"material":"FEP","numberLayers":3,"temperatureRating":155,"thickness":null,"thicknessLayers":0.0000762,"type":"insulated"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001024},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":null,"numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001095},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"18 AWG","strand":null,"type":"round"})"));
        WireType newWireType;
        from_json("litz", newWireType);
        double effectivefrequency = 110746;

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, newWireType, effectivefrequency);
        CHECK(newWire.get_type() == WireType::LITZ);
    }
}
