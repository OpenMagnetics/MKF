#include <source_location>
#include "constructive_models/Wire.h"
#include "support/Utils.h"
#include "json.hpp"

#include <numbers>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    auto masPath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("MAS/").string();
    double max_error = 0.07;
    auto label = WaveformLabel::TRIANGULAR;
    double offset = 0;
    double peakToPeak = 2 * 1.73205;
    double dutyCycle = 0.5;
    double magnetizingInductance = 1e-3;
    double temperature = 22;

    OpenMagnetics::Inputs setup_inputs(double frequency) {

        ProcessedWaveform processed;
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

    TEST_CASE("Test_Sample_Wire", "[constructive-model][wire][smoke-test]") {
        auto wireFilePath = masPath + "samples/magnetic/wire/round/0.000140.json";
        std::ifstream json_file(wireFilePath);
        auto wireJson = json::parse(json_file);

        OpenMagnetics::Wire wire(wireJson);

        auto conductingDiameter = wire.get_conducting_diameter().value().get_nominal().value();

        REQUIRE(conductingDiameter == wireJson["conductingDiameter"]["nominal"]);
    }

    // A litz bundle's own coating describes the SERVING; the insulation that separates one
    // conductor from another is the STRAND enamel. Reading the bundle used to report exactly
    // zero for every catalogue litz — no entry carries a serving thickness and none sets an
    // outerDiameter, so the derived (outer - conducting)/2 collapsed to (strand - strand)/2.
    // A user hit this as COIL_SHORTED_TURNS on a perfectly valid close-wound litz coil.
    TEST_CASE("Litz reports the strand enamel as its turn-to-turn insulation", "[constructive-model][wire][litz]") {
        for (const auto& name : {std::string("Litz 12x0.03 - Grade 1 - Double Served"),
                                 std::string("Litz 12x0.03 - Grade 1 - Unserved")}) {
            auto wire = OpenMagnetics::find_wire_by_name(name);
            auto strandConductingDiameter =
                resolve_dimensional_values(wire.resolve_strand().get_conducting_diameter());

            INFO("wire: " << name);
            // The whole point: not zero, and physically sane (a fraction of the strand).
            REQUIRE(wire.get_coating_thickness() > 0);
            REQUIRE(wire.get_coating_thickness() < strandConductingDiameter);
        }
    }

    // ABT #857. The thermal model asks a litz wire what its coating is made of, in order to
    // get the inter-strand matrix conductivity -- Temperature.cpp documents that value as
    // "the strand enamel/insulation that heat must cross between conductors". It used to ask
    // the wire's OWN coating, which for litz is the bundle SERVING: {type: served,
    // material: null} in every catalogue and synthesized litz. resolve_coating_insulation_material
    // then threw "Coating is missing material information" and the whole temperature field
    // died. Measured on the production bug queue: 11 of 29 recent designs failed
    // plot_temperature_field and every single one was litz, while every round, foil and
    // planar design succeeded.
    TEST_CASE("Litz resolves its coating material from the strand enamel rather than the bundle serving", "[constructive-model][wire][litz][coating]") {
        for (const auto& name : {std::string("Litz 12x0.03 - Grade 1 - Double Served"),
                                 std::string("Litz 12x0.03 - Grade 1 - Unserved")}) {
            auto wire = OpenMagnetics::find_wire_by_name(name);
            INFO("wire: " << name);

            // Before the fix the SERVED bundle threw here rather than returning a number.
            double thermalConductivity = wire.get_coating_thermal_conductivity();
            REQUIRE(thermalConductivity > 0);
            // Enamel, not copper: polyurethane/polyesterimide sits around 0.2-0.3 W/(m K).
            // A value anywhere near copper's 400 would mean we had resolved the conductor.
            REQUIRE(thermalConductivity < 1.0);

            // The dielectric properties travel through the same resolution point, so they
            // must answer too rather than throwing on the serving's absent material.
            REQUIRE(wire.get_coating_relative_permittivity() > 1.0);
        }
    }

    // The outer envelope must be the BUNDLE, not one strand. get_maximum_outer_width() used to
    // fall back to get_maximum_conducting_width(), which for litz returns the strand diameter —
    // so a 12-strand bundle claimed to be as wide as one of its own strands.
    TEST_CASE("Litz outer width is the bundle envelope rather than a single strand", "[constructive-model][wire][litz]") {
        auto wire = OpenMagnetics::find_wire_by_name("Litz 12x0.03 - Grade 1 - Double Served");
        auto strandConductingDiameter =
            resolve_dimensional_values(wire.resolve_strand().get_conducting_diameter());

        REQUIRE(wire.get_maximum_outer_width() > strandConductingDiameter);
        // 12 strands cannot fit inside less than ~sqrt(12) strand diameters.
        REQUIRE(wire.get_maximum_outer_width() > 3 * strandConductingDiameter);
        REQUIRE(wire.get_maximum_outer_width() == Catch::Approx(wire.get_maximum_outer_height()));
    }

    // The guard that catches a REAL fault must survive: strands with no insulation of their own
    // genuinely leave copper exposed, and close-wound turns of that are shorted.
    TEST_CASE("Litz with bare strands still reports no insulation", "[constructive-model][wire][litz]") {
        auto wire = OpenMagnetics::find_wire_by_name("Litz 12x0.03 - Grade 1 - Double Served");
        auto strand = wire.resolve_strand();
        InsulationWireCoating bareCoating;
        bareCoating.set_type(InsulationWireCoatingType::BARE);
        strand.set_coating(bareCoating);
        wire.set_strand(strand);

        REQUIRE(wire.get_coating_thickness() == 0);
    }

    // ABT #853: a wire that DECLARES a coating but gives neither its thickness nor an outer
    // dimension used to derive (outer - conducting) / 2 with the outer dimension defaulted to the
    // conducting one -- "absent" read as "equal to the conductor", i.e. as bare metal -- and that
    // zero then declared close-wound turns a short circuit. Enamelled wire with a grade resolves
    // from the standard's table; anything else unresolvable is an incomplete wire and throws.
    // A wire with NO coating at all stays bare (0): foil and planar conductors legitimately rely
    // on separate insulation layers, and that declaration is explicit, not missing.
    TEST_CASE("Declared coating without thickness or outer dimension resolves from the table or throws, never 0", "[constructive-model][wire][coating]") {
        auto rectangular = OpenMagnetics::Wire(json::parse(R"({"type": "rectangular", "material": "copper", "numberConductors": 1,
            "conductingWidth": {"nominal": 0.004}, "conductingHeight": {"nominal": 0.001},
            "coating": {"type": "enamelled", "grade": 1}})"));
        double rectangularCoating = rectangular.get_coating_thickness();
        INFO("rectangular grade-1 coating: " << rectangularCoating);
        REQUIRE(rectangularCoating > 0);
        REQUIRE(rectangularCoating < 0.001);  // a fraction of the conductor, not the conductor
        REQUIRE(rectangularCoating == Catch::Approx(std::min(
            (OpenMagnetics::Wire::get_outer_width_rectangular(0.004, 1, WireStandard::IEC_60317) - 0.004) / 2,
            (OpenMagnetics::Wire::get_outer_height_rectangular(0.001, 1, WireStandard::IEC_60317) - 0.001) / 2)));

        auto round = OpenMagnetics::Wire(json::parse(R"({"type": "round", "material": "copper", "numberConductors": 1,
            "conductingDiameter": {"nominal": 0.001}, "coating": {"type": "enamelled", "grade": 2}})"));
        double roundCoating = round.get_coating_thickness();
        INFO("round grade-2 coating: " << roundCoating);
        REQUIRE(roundCoating > 0);
        REQUIRE(roundCoating < 0.001);
        REQUIRE(roundCoating == Catch::Approx(
            (OpenMagnetics::Wire::get_outer_diameter_round(0.001, 2, WireStandard::IEC_60317) - 0.001) / 2));

        // Declares an insulation but nothing that sizes it: cannot be derived, must not be 0.
        auto foil = OpenMagnetics::Wire(json::parse(R"({"type": "foil", "material": "copper", "numberConductors": 1,
            "conductingWidth": {"nominal": 0.0002}, "conductingHeight": {"nominal": 0.01},
            "coating": {"type": "insulated"}})"));
        REQUIRE_THROWS_AS(foil.get_coating_thickness(), OpenMagnetics::InvalidInputException);

        // Explicitly bare: no coating key at all. Still 0 -- that is a declaration, not a gap.
        auto bareFoil = OpenMagnetics::Wire(json::parse(R"({"type": "foil", "material": "copper", "numberConductors": 1,
            "conductingWidth": {"nominal": 0.0002}, "conductingHeight": {"nominal": 0.01}})"));
        REQUIRE(bareFoil.get_coating_thickness() == 0);

        // Outer dimensions given: the derived route is still the answer.
        auto foilWithOuter = OpenMagnetics::Wire(json::parse(R"({"type": "foil", "material": "copper", "numberConductors": 1,
            "conductingWidth": {"nominal": 0.0002}, "conductingHeight": {"nominal": 0.01},
            "outerWidth": {"nominal": 0.00025}, "outerHeight": {"nominal": 0.01005},
            "coating": {"type": "insulated"}})"));
        REQUIRE(foilWithOuter.get_coating_thickness() == Catch::Approx(0.000025));
    }

    TEST_CASE("Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_1", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05);
        double expectedValue = 0.755;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_1", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.1e-05);
        double expectedValue = 0.64;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_1", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048);
        double expectedValue = 0.87;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 2);
        double expectedValue = 0.616;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.1e-05, 2);
        double expectedValue = 0.455;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 2);
        double expectedValue = 0.8;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_3", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 3);
        double expectedValue = 0.523;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_3", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.1e-05, 3);
        double expectedValue = 0.334;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_3", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 3);
        double expectedValue = 0.741;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_1_Nema", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.79;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_1_Nema", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.3e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.71;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_1_Nema", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 1, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.89;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_2_Nema", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 2, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.65;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_2_Nema", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(1.3e-05, 2, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.52;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_2_Nema", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 2, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.81;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Medium_Round_Enamelled_Wire_Grade_3_Nema", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(5.4e-05, 3, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.55;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Small_Round_Enamelled_Wire_Grade_3_Nema", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(4e-05, 3, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.51;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Large_Round_Enamelled_Wire_Grade_3_Nema", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00048, 3, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.74;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Thin_Round_Insulated_Wire_1_Layer_Thin_Layer_Thickness", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.000101, 1, 3.81e-05, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.321961;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Filling_Factors_Thick_Round_Insulated_Wire_1_Layer_Thick_Layer_Thickness", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00129, 1, 7.62e-05, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.799184;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Outer_Diameter_Thick_Round_Insulated_Wire_1_Layer_Thick_Layer_Thickness", "[constructive-model][wire][smoke-test]") {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_round(0.00129, 1, 7.62e-05, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00144;

        REQUIRE_THAT(expectedOuterDiameter, Catch::Matchers::WithinAbs(outerDiameter, max_error * expectedOuterDiameter));
    }

    TEST_CASE("Test_Filling_Factors_Thick_Round_Insulated_Wire_3_Layer_Thick_Layer_Thickness", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_round(0.00129, 3, 7.62e-05, WireStandard::NEMA_MW_1000_C);
        double expectedValue = 0.5446;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Outer_Diameter_Thick_Round_Insulated_Wire_3_Layers_Thick_Layer_Thickness", "[constructive-model][wire][smoke-test]") {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_round(0.00129, 3, 7.62e-05, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.001748;

        REQUIRE_THAT(expectedOuterDiameter, Catch::Matchers::WithinAbs(outerDiameter, max_error * expectedOuterDiameter));
    }

    TEST_CASE("Test_Filling_Factor_Thick_Litz_Wire_Served_1_Layer_Few_Strands", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_served_litz(0.0001, 66, 1, 1, WireStandard::IEC_60317, false);
        double expectedValue = 0.458122;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedValue));
    }

    TEST_CASE("Test_Outer_Diameter_Litz_Wire_Grade_1_Served_1_Layer_Few_Strands", "[constructive-model][wire][smoke-test]") {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_served_litz(2e-05, 10, 1, 1, WireStandard::IEC_60317);
        double expectedMinimumValue = 0.000112;
        double expectedMaximumValue = 0.000142;

        REQUIRE(outerDiameter > expectedMinimumValue);
        REQUIRE(outerDiameter < expectedMaximumValue);
    }

    TEST_CASE("Test_Outer_Diameter_Litz_Wire_Grade_2_Served_1_Layer_1_Few_Strands", "[constructive-model][wire][smoke-test]") {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_served_litz(2e-05, 10, 2, 1, WireStandard::IEC_60317);
        double expectedOuterDiameter = 0.000153;

        REQUIRE_THAT(expectedOuterDiameter, Catch::Matchers::WithinAbs(outerDiameter, max_error * expectedOuterDiameter));
    }

    TEST_CASE("Test_Outer_Diameter_Thick_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Few_Strands", "[constructive-model][wire][smoke-test]") {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_insulated_litz(0.000102, 66, 3, 7.62e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00152908;

        REQUIRE_THAT(expectedOuterDiameter, Catch::Matchers::WithinAbs(outerDiameter, max_error * expectedOuterDiameter));
    }

    TEST_CASE("Test_Outer_Diameter_Thin_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Many_Strands", "[constructive-model][wire][smoke-test]") {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_insulated_litz(5.1e-05, 825, 3, 7.62e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00253238;

        REQUIRE_THAT(expectedOuterDiameter, Catch::Matchers::WithinAbs(outerDiameter, max_error * expectedOuterDiameter));
    }

    TEST_CASE("Test_Outer_Diameter_Thin_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Many_Strands_Diameter_Not_In_Db", "[constructive-model][wire][smoke-test]") {
        auto outerDiameter = OpenMagnetics::Wire::get_outer_diameter_insulated_litz(5.42e-05, 825, 3, 7.62e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.00253238;

        REQUIRE(expectedOuterDiameter < outerDiameter);
    }

    TEST_CASE("Test_Filling_Factor_Litz_Wire_Grade_2_Served_1_Layer_1_Few_Strands", "[constructive-model][wire][smoke-test]") {
        load_wires();
        auto outerDiameter = OpenMagnetics::Wire::get_filling_factor_served_litz(2e-05, 10, 2, 1, WireStandard::IEC_60317);
        double expectedOuterDiameter = 0.33978;

        REQUIRE_THAT(expectedOuterDiameter, Catch::Matchers::WithinAbs(outerDiameter, max_error * expectedOuterDiameter));
    }

    TEST_CASE("Test_Filling_Factor_Thick_Litz_Insulated_Wire_3_Layers_Thick_Layer_Thickness_Few_Strands", "[constructive-model][wire][smoke-test]") {
        load_wires();
        auto outerDiameter = OpenMagnetics::Wire::get_filling_factor_insulated_litz(0.000102, 66, 3, 7.62e-05, 1, WireStandard::NEMA_MW_1000_C);
        double expectedOuterDiameter = 0.45201;

        REQUIRE_THAT(expectedOuterDiameter, Catch::Matchers::WithinAbs(outerDiameter, max_error * expectedOuterDiameter));
    }

    TEST_CASE("Test_Outer_Width_Small_Rectangular_Grade_1", "[constructive-model][wire][smoke-test]") {
        auto outerWidth = OpenMagnetics::Wire::get_outer_width_rectangular(0.002, 1, WireStandard::IEC_60317);
        double expectedOuterWidth = 0.00206;

        REQUIRE_THAT(expectedOuterWidth, Catch::Matchers::WithinAbs(outerWidth, max_error * expectedOuterWidth));
    }

    TEST_CASE("Test_Outer_Width_Small_Rectangular_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto outerWidth = OpenMagnetics::Wire::get_outer_width_rectangular(0.002, 2, WireStandard::IEC_60317);
        double expectedOuterWidth = 0.00217;

        REQUIRE_THAT(expectedOuterWidth, Catch::Matchers::WithinAbs(outerWidth, max_error * expectedOuterWidth));
    }

    TEST_CASE("Test_Outer_Width_Large_Rectangular_Grade_1", "[constructive-model][wire][smoke-test]") {
        auto outerWidth = OpenMagnetics::Wire::get_outer_width_rectangular(0.016, 1, WireStandard::IEC_60317);
        double expectedOuterWidth = 0.01608;

        REQUIRE_THAT(expectedOuterWidth, Catch::Matchers::WithinAbs(outerWidth, max_error * expectedOuterWidth));
    }

    TEST_CASE("Test_Outer_Width_Large_Rectangular_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto outerWidth = OpenMagnetics::Wire::get_outer_width_rectangular(0.016, 2, WireStandard::IEC_60317);
        double expectedOuterWidth = 0.01614;

        REQUIRE_THAT(expectedOuterWidth, Catch::Matchers::WithinAbs(outerWidth, max_error * expectedOuterWidth));
    }

    TEST_CASE("Test_Outer_Height_Small_Rectangular_Grade_1", "[constructive-model][wire][smoke-test]") {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(0.0008, 1, WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00088;

        REQUIRE_THAT(expectedOuterHeight, Catch::Matchers::WithinAbs(outerHeight, max_error * expectedOuterHeight));
    }

    TEST_CASE("Test_Outer_Height_Small_Rectangular_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(0.0008, 2, WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00092;

        REQUIRE_THAT(expectedOuterHeight, Catch::Matchers::WithinAbs(outerHeight, max_error * expectedOuterHeight));
    }

    TEST_CASE("Test_Outer_Height_Large_Rectangular_Grade_1", "[constructive-model][wire][smoke-test]") {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(0.0045, 1, WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00456;

        REQUIRE_THAT(expectedOuterHeight, Catch::Matchers::WithinAbs(outerHeight, max_error * expectedOuterHeight));
    }

    TEST_CASE("Test_Outer_Height_Large_Rectangular_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(0.0045, 2, WireStandard::IEC_60317);
        double expectedOuterHeight = 0.00467;

        REQUIRE_THAT(expectedOuterHeight, Catch::Matchers::WithinAbs(outerHeight, max_error * expectedOuterHeight));
    }

    TEST_CASE("Test_Filling_Factor_Small_Rectangular_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_rectangular(0.002, 0.0008, 2, WireStandard::IEC_60317);
        double expectedFillingFactor = 0.720267;

        REQUIRE_THAT(expectedFillingFactor, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedFillingFactor));
    }

    TEST_CASE("Test_Filling_Factor_Large_Rectangular_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto fillingFactor = OpenMagnetics::Wire::get_filling_factor_rectangular(0.016, 0.0045, 2, WireStandard::IEC_60317);
        double expectedFillingFactor = 0.948615;

        REQUIRE_THAT(expectedFillingFactor, Catch::Matchers::WithinAbs(fillingFactor, max_error * expectedFillingFactor));
    }

    TEST_CASE("Test_Conducting_Area_Small_Rectangular", "[constructive-model][wire][smoke-test]") {
        auto conductingArea = OpenMagnetics::Wire::get_conducting_area_rectangular(0.002, 0.0008, WireStandard::IEC_60317);
        double expectedConductingArea = 0.00000146;

        REQUIRE_THAT(expectedConductingArea, Catch::Matchers::WithinAbs(conductingArea, max_error * expectedConductingArea));
    }

    TEST_CASE("Test_Conducting_Area_Large_Rectangular", "[constructive-model][wire][smoke-test]") {
        auto conductingArea = OpenMagnetics::Wire::get_conducting_area_rectangular(0.016, 0.0045, WireStandard::IEC_60317);
        double expectedConductingArea = 0.00007114;

        REQUIRE_THAT(expectedConductingArea, Catch::Matchers::WithinAbs(conductingArea, max_error * expectedConductingArea));
    }

    TEST_CASE("Test_Conducting_Area_Very_Large_Rectangular_BugFix", "[constructive-model][wire][smoke-test]") {
        // Regression test for bug where spline interpolator returned crazy values
        // for large wires (e.g., 13.5mm x 2.2mm), causing conducting area to be
        // ~4000x too large and winding losses ~4000x too low.
        // The fix ensures we fall back to theoretical area when proportion is out of bounds.
        double conductingWidth = 0.0135;   // 13.5mm
        double conductingHeight = 0.0022;  // 2.2mm
        double theoreticalArea = conductingWidth * conductingHeight;  // 29.7 mm²
        
        auto conductingArea = OpenMagnetics::Wire::get_conducting_area_rectangular(conductingWidth, conductingHeight, WireStandard::IEC_60317);
        
        // Area should be close to theoretical (within ~5% for enamel coating)
        // Before fix: area was ~122,000 mm² (proportion ~4118)
        // After fix: area should be ~29.7 mm² (theoretical)
        REQUIRE(conductingArea > theoreticalArea * 0.8);  // At least 80% of theoretical
        REQUIRE(conductingArea < theoreticalArea * 1.1);  // At most 110% of theoretical
    }

    TEST_CASE("Test_Outer_Height_Tiny_Rectangular_Grade_2", "[constructive-model][wire][smoke-test]") {
        auto outerHeight = OpenMagnetics::Wire::get_outer_height_rectangular(1e-9, 2, WireStandard::IEC_60317);
        double expectedOuterHeight = 1.2e-9;

        REQUIRE_THAT(expectedOuterHeight, Catch::Matchers::WithinAbs(outerHeight, max_error * expectedOuterHeight));
    }

    TEST_CASE("Test_Outer_Diameter_Litz_Wire_Unserved_Medium_Strands", "[constructive-model][wire][smoke-test]") {
        auto diameter = OpenMagnetics::Wire::get_outer_diameter_served_litz(0.000071, 270, 1, 0, WireStandard::IEC_60317);
        double expectedMaximumValue = 0.001767;
        double expectedMinimumValue = 0.001641;

        REQUIRE(diameter > expectedMinimumValue);
        REQUIRE(diameter < expectedMaximumValue);
    }

    TEST_CASE("Test_Outer_Diameter_Litz_Wire_Served_1_Layer_Medium_Strands", "[constructive-model][wire][smoke-test]") {
        auto diameter = OpenMagnetics::Wire::get_outer_diameter_served_litz(0.000071, 270, 1, 1, WireStandard::IEC_60317);
        double expectedMaximumValue = 0.001807;
        double expectedMinimumValue = 0.001666;

        REQUIRE(diameter > expectedMinimumValue);
        REQUIRE(diameter < expectedMaximumValue);
    }

    TEST_CASE("Test_Effective_Current_Density_Medium_Frequency_Round_Operation_Point", "[constructive-model][wire][smoke-test]") {
        double frequency = 100000;
        ProcessedWaveform processed;
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

        REQUIRE_THAT(expectedEffectiveCurrentDensity, Catch::Matchers::WithinAbs(effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity));
    }

    TEST_CASE("Test_Effective_Current_Density_Low_Frequency_Round", "[constructive-model][wire][smoke-test]") {
        double frequency = 10;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 5.093e6;

        REQUIRE_THAT(expectedEffectiveCurrentDensity, Catch::Matchers::WithinAbs(effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity));
    }

    TEST_CASE("Test_Effective_Current_Density_Medium_Frequency_Round", "[constructive-model][wire][smoke-test]") {
        double frequency = 100000;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 5.283e6;

        REQUIRE_THAT(expectedEffectiveCurrentDensity, Catch::Matchers::WithinAbs(effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity));
    }

    TEST_CASE("Test_Effective_Current_Density_High_Frequency_Round", "[constructive-model][wire][smoke-test]") {
        double frequency = 1000000;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 11.19e6;

        REQUIRE_THAT(expectedEffectiveCurrentDensity, Catch::Matchers::WithinAbs(effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity));
    }

    TEST_CASE("Test_Effective_Current_Density_Low_Frequency_Litz", "[constructive-model][wire][smoke-test]") {
        double frequency = 10;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 5.093e5;

        REQUIRE_THAT(expectedEffectiveCurrentDensity, Catch::Matchers::WithinAbs(effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity));
    }

    TEST_CASE("Test_Effective_Current_Density_High_Frequency_Litz", "[constructive-model][wire][smoke-test]") {
        double frequency = 10000000;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 5.24e5;

        REQUIRE_THAT(expectedEffectiveCurrentDensity, Catch::Matchers::WithinAbs(effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity));
    }

    TEST_CASE("Test_Effective_Current_Density_Low_Frequency_Rectangular", "[constructive-model][wire][smoke-test]") {
        double frequency = 10;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 3.96e5;

        REQUIRE_THAT(expectedEffectiveCurrentDensity, Catch::Matchers::WithinAbs(effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity));
    }

    TEST_CASE("Test_Effective_Current_Density_High_Frequency_Rectangular", "[constructive-model][wire][smoke-test]") {
        double frequency = 1000000;
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(1, frequency, temperature);
        double expectedEffectiveCurrentDensity = 2.09e6;

        REQUIRE_THAT(expectedEffectiveCurrentDensity, Catch::Matchers::WithinAbs(effectiveCurrentDensity, max_error * expectedEffectiveCurrentDensity));
    }

    TEST_CASE("Test_Conducting_Area_Large_Rectangular_2", "[constructive-model][wire][smoke-test]") {
        double frequency = 10;

        ProcessedWaveform processed;
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
        // NOTE: the original conducting-area computation this test was written around was
        // deleted at some point, leaving only commented-out expectations. Pin what is left:
        // the wire lookup and its standard conducting area (pi * (0.25 mm)^2).
        REQUIRE(wire.get_conducting_diameter());
        double conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
        REQUIRE_THAT(conductingDiameter, Catch::Matchers::WithinAbs(0.0005, 0.0005 * 0.01));
        // Database wires carry the diameter, not a materialized conducting area; pin the
        // area implied by the standard diameter instead (pi * (d/2)^2).
        double conductingArea = std::numbers::pi * conductingDiameter * conductingDiameter / 4;
        REQUIRE_THAT(conductingArea, Catch::Matchers::WithinAbs(1.9635e-7, 1.9635e-7 * 0.02));
        auto& excitation = inputs.get_mutable_operating_points()[0].get_mutable_excitations_per_winding()[0];
        REQUIRE(excitation.get_current()->get_processed());
        CHECK(excitation.get_current()->get_processed()->get_rms().value() > 0);
    }

    TEST_CASE("Test_Number_Parallels_Low_Frequency_Round_1_Parallel", "[constructive-model][wire][smoke-test]") {
        double frequency = 10;

        auto inputs = setup_inputs(frequency);
        
        auto wire = find_wire_by_name("Round 0.5 - Grade 1");
        double maximumEffectiveCurrentDensity = 5.5e6;
        double expectedNumberParallels = 1;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        REQUIRE(expectedNumberParallels == numberParallels);
    }

    TEST_CASE("Test_Number_Parallels_Low_Frequency_Round_2_Parallels", "[constructive-model][wire][smoke-test]") {
        double frequency = 10;

        auto inputs = setup_inputs(frequency);
        
        auto wire = find_wire_by_name("Round 0.5 - Grade 1");
        double maximumEffectiveCurrentDensity = 5e6;
        double expectedNumberParallels = 2;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        REQUIRE(expectedNumberParallels == numberParallels);
    }

    TEST_CASE("Test_Number_Parallels_High_Frequency_Round_3_Parallels", "[constructive-model][wire][smoke-test]") {
        double frequency = 1000000;

        auto inputs = setup_inputs(frequency);
        
        auto wire = find_wire_by_name("Round 0.5 - Grade 1");
        double maximumEffectiveCurrentDensity = 5e6;
        double expectedNumberParallels = 3;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        REQUIRE(expectedNumberParallels == numberParallels);
    }

    TEST_CASE("Test_Number_Parallels_High_Frequency_Litz_2_Parallels", "[constructive-model][wire][smoke-test]") {
        double frequency = 1000000;

        auto inputs = setup_inputs(frequency);
        
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        double maximumEffectiveCurrentDensity = 5e5;
        double expectedNumberParallels = 2;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        REQUIRE(expectedNumberParallels == numberParallels);
    }

    TEST_CASE("Test_Number_Parallels_Low_Frequency_Rectangular_1_Parallels", "[constructive-model][wire][smoke-test]") {
        double frequency = 10;

        auto inputs = setup_inputs(frequency);
        
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        double maximumEffectiveCurrentDensity = 5e6;
        double expectedNumberParallels = 1;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        REQUIRE(expectedNumberParallels == numberParallels);
    }

    TEST_CASE("Test_Number_Parallels_High_Frequency_Rectangular_4_Parallels", "[constructive-model][wire][smoke-test]") {
        double frequency = 1000000;

        auto inputs = setup_inputs(frequency);
        
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        double maximumEffectiveCurrentDensity = 1e6;
        double expectedNumberParallels = 3;
        auto numberParallels = OpenMagnetics::Wire::calculate_number_parallels_needed(inputs, wire, maximumEffectiveCurrentDensity);

        REQUIRE(expectedNumberParallels == numberParallels);
    }

    TEST_CASE("Test_Coating_Label_Uniqueness", "[constructive-model][wire][smoke-test]") {
        auto wires = get_wires();
        std::vector<std::string> coatingLabels;
        for (auto wire : wires) {
            auto coatingLabel = wire.encode_coating_label();
            coatingLabels.push_back(coatingLabel);
        }

        REQUIRE(std::find(coatingLabels.begin(), coatingLabels.end(), "Bare") != coatingLabels.end());
    }

    TEST_CASE("Test_Coating_Decoding", "[constructive-model][wire][smoke-test]") {
        auto wires = get_wires();
        for (auto wire : wires) {
            auto coatingLabel = wire.encode_coating_label();
            auto coating = wire.resolve_coating();
            if (coating) {
                auto decodedCoating = OpenMagnetics::Wire::decode_coating_label(coatingLabel);
                REQUIRE(coating->get_type().value() == decodedCoating->get_type().value());
                if (coating->get_number_layers()) {
                    REQUIRE(coating->get_number_layers().value() == decodedCoating->get_number_layers().value());
                }
                if (coating->get_temperature_rating()) {
                    REQUIRE(coating->get_temperature_rating().value() == decodedCoating->get_temperature_rating().value());
                }
                if (coating->get_breakdown_voltage() && coating->get_type().value() == InsulationWireCoatingType::INSULATED) {
                    REQUIRE(coating->get_breakdown_voltage().value() == decodedCoating->get_breakdown_voltage().value());
                }
                if (coating->get_grade()) {
                    REQUIRE(coating->get_grade().value() == decodedCoating->get_grade().value());
                }
            }
        }
    }

    TEST_CASE("Test_Coating_Relative_Permittivity", "[constructive-model][wire][smoke-test]") {
        auto wire = OpenMagnetics::Wire(find_wire_by_name("Round 0.80 - Grade 1"));
        auto relativePermittivity = wire.get_coating_relative_permittivity();
        REQUIRE_THAT(3.7, Catch::Matchers::WithinAbs(relativePermittivity, max_error * 3.7));
    }

    TEST_CASE("Test_Coating_Relative_Permittivity_Web0", "[constructive-model][wire][smoke-test]") {
        auto wire = OpenMagnetics::Wire(json::parse(R"({"type": "round", "conductingDiameter": {"nominal": 0.001}, "material": "copper", "outerDiameter": {"nominal": 0.001062}, "coating": {"breakdownVoltage": 2700.0, "grade": 1, "type": "enamelled"}, "name": "Round 1.00 - Grade 1", "numberConductors": 1, "standard": "IEC 60317", "standardName": "1.00 mm"})"));
        auto relativePermittivity = wire.get_coating_relative_permittivity();
        REQUIRE_THAT(3.7, Catch::Matchers::WithinAbs(relativePermittivity, max_error * 3.7));
    }

    TEST_CASE("Test_Coating_Material", "[constructive-model][wire][smoke-test]") {
        auto material = find_insulation_material_by_name("ETFE");
        json jsonDta;
        to_json(jsonDta, material);
        // The ETFE insulation material must round-trip with its key dielectric data.
        REQUIRE(!jsonDta.empty());
        CHECK(jsonDta["name"] == "ETFE");
        REQUIRE(jsonDta.contains("dielectricStrength"));
        CHECK(jsonDta["dielectricStrength"].size() > 0);
    }

    TEST_CASE("Test_Find_Round_By_Dimension_European", "[constructive-model][wire][smoke-test]") {
        auto wire = find_wire_by_dimension(0.00072, WireType::ROUND, WireStandard::IEC_60317, false);
        REQUIRE(wire.get_standard_name().value() == "0.71 mm");
    }

    TEST_CASE("Test_Find_Round_By_Dimension_American", "[constructive-model][wire][smoke-test]") {
        auto wire = find_wire_by_dimension(0.00072, WireType::ROUND, WireStandard::NEMA_MW_1000_C);
        REQUIRE(wire.get_standard_name().value() == "21 AWG");
    }

    TEST_CASE("Test_Find_Among_All_By_Dimension", "[constructive-model][wire][smoke-test]") {
        auto wire = find_wire_by_dimension(0.00072);
        REQUIRE(wire.get_standard_name().value() == "21 AWG");
    }

    TEST_CASE("Test_Find_Rectangular_By_Dimension", "[constructive-model][wire][smoke-test]") {
        auto wire = find_wire_by_dimension(0.00072, WireType::RECTANGULAR);
        auto conductingHeight = resolve_dimensional_values(wire.get_conducting_height().value());
        REQUIRE(conductingHeight == 0.0008);
    }

    TEST_CASE("Test_Find_Foil_By_Dimension", "[constructive-model][wire][smoke-test]") {
        auto wire = find_wire_by_dimension(0.00072, WireType::FOIL);
        auto conductingWidth = resolve_dimensional_values(wire.get_conducting_width().value());
        REQUIRE(conductingWidth == 0.0007);
    }

    TEST_CASE("Test_Litz_To_Litz_Equivalent", "[constructive-model][wire][smoke-test]") {
        double effectiveFrequency = 1234981;
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::LITZ, effectiveFrequency);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::LITZ);
        REQUIRE(numberConductors == 1000);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.00005, max_error * 0.00005));
    }

    TEST_CASE("Test_Round_To_Litz_Equivalent", "[constructive-model][wire][smoke-test]") {
        double effectiveFrequency = 1234981;
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::LITZ, effectiveFrequency);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::LITZ);
        REQUIRE(numberConductors == 71);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.00006, max_error * 0.00006));
    }

    TEST_CASE("Test_Rectangular_To_Litz_Equivalent", "[constructive-model][wire][smoke-test]") {
        double effectiveFrequency = 1234981;
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::LITZ, effectiveFrequency);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::LITZ);
        REQUIRE(numberConductors == 914);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.00006, max_error * 0.00006));
    }

    TEST_CASE("Test_Foil_To_Litz_Equivalent", "[constructive-model][wire][smoke-test]") {
        double effectiveFrequency = 1234981;
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Foil 0.2"));
        oldWire.set_nominal_value_conducting_height(0.01);

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::LITZ, effectiveFrequency);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::LITZ);
        REQUIRE(numberConductors == 725);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.00006, max_error * 0.00006));
    }

    TEST_CASE("Test_Litz_To_Round_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::ROUND);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::ROUND);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0016, max_error * 0.0016));
    }

    TEST_CASE("Test_Round_To_Round_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Round 0.5 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::ROUND);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::ROUND);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0005, max_error * 0.0005));
    }

    TEST_CASE("Test_Rectangular_To_Round_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::ROUND);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::ROUND);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0009, max_error * 0.0009));
    }

    TEST_CASE("Test_Foil_To_Round_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Foil 0.2"));
        oldWire.set_nominal_value_conducting_height(0.001);

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::ROUND);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::ROUND);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0002, max_error * 0.0002));
    }

    TEST_CASE("Test_Litz_To_Rectangular_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::RECTANGULAR);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::RECTANGULAR);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0016, max_error * 0.0016));
    }

    TEST_CASE("Test_Round_To_Rectangular_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Round 0.80 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::RECTANGULAR);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::RECTANGULAR);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0008, max_error * 0.0008));
    }

    TEST_CASE("Test_Rectangular_To_Rectangular_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::RECTANGULAR);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::RECTANGULAR);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.00085, max_error * 0.00085));
    }

    TEST_CASE("Test_Foil_To_Rectangular_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Foil 0.2"));
        oldWire.set_nominal_value_conducting_height(0.001);

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::RECTANGULAR);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::RECTANGULAR);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0008, max_error * 0.0008));
    }

    TEST_CASE("Test_Litz_To_Foil_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Litz 1000x0.05 - Grade 1 - Single Served"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::FOIL);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::FOIL);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0016, max_error * 0.0016));
    }

    TEST_CASE("Test_Round_To_Foil_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Round 0.80 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::FOIL);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::FOIL);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0008, max_error * 0.0008));
    }

    TEST_CASE("Test_Rectangular_To_Foil_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Rectangular 3.15x0.85 - Grade 1"));
        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::FOIL);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::FOIL);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0008, max_error * 0.0008));
    }

    TEST_CASE("Test_Foil_To_Foil_Equivalent", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(find_wire_by_name("Foil 0.2"));
        oldWire.set_nominal_value_conducting_height(0.001);

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, WireType::FOIL);

        auto conductingDimension = resolve_dimensional_values(newWire.get_minimum_conducting_dimension());
        auto numberConductors = newWire.get_number_conductors().value();
        REQUIRE(newWire.get_type() == WireType::FOIL);
        REQUIRE(numberConductors == 1);
        REQUIRE_THAT(conductingDimension, Catch::Matchers::WithinAbs(0.0002, max_error * 0.0002));
    }

    TEST_CASE("Test_Equivalent_Web_0", "[constructive-model][wire][smoke-test]") {
        auto oldWire = OpenMagnetics::Wire(json::parse(R"({"coating":{"breakdownVoltage":13500,"grade":null,"material":"FEP","numberLayers":3,"temperatureRating":155,"thickness":null,"thicknessLayers":0.0000762,"type":"insulated"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001024},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":null,"numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001095},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"18 AWG","strand":null,"type":"round"})"));
        WireType newWireType;
        from_json("litz", newWireType);
        double effectivefrequency = 110746;

        auto newWire = OpenMagnetics::Wire::get_equivalent_wire(oldWire, newWireType, effectivefrequency);
        REQUIRE(newWire.get_type() == WireType::LITZ);
    }


    // ABT #898: the wire the ADVISERS synthesise must carry the coating whose table
    // produced its outer diameter. get_wire_for_frequency/get_wire_for_conducting_area in
    // `exact` mode size the outer diameter with get_outer_diameter_round(d, grade 1,
    // IEC 60317) -- an ENAMELLED grade-1 wire -- but used to stamp no coating at all, so
    // Wire::get_coating_thickness() read "no coating object" and answered 0. Every turn of
    // the fast adviser's dummy coil was then bare copper, and StrayCapacitance rightly
    // called close-wound bare turns a short circuit / an infinite capacitance.
    TEST_CASE("Test_Synthesised_Adviser_Wire_Is_Insulated", "[constructive-model][wire][abt898]") {
        for (auto wire : {OpenMagnetics::Wire::get_wire_for_frequency(100000, 25, true),
                          OpenMagnetics::Wire::get_wire_for_conducting_area(1e-7, 25, true)}) {
            auto conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
            auto outerDiameter = resolve_dimensional_values(wire.get_outer_diameter().value());
            REQUIRE(outerDiameter > conductingDiameter);
            // The coating must exist, must not claim to be bare, and must measure the very
            // gap the outer diameter already describes.
            auto coating = wire.resolve_coating();
            REQUIRE(coating);
            REQUIRE(coating->get_type().value() != InsulationWireCoatingType::BARE);
            REQUIRE(wire.get_coating_thickness() > 0);
            REQUIRE_THAT(wire.get_coating_thickness(),
                         Catch::Matchers::WithinRel((outerDiameter - conductingDiameter) / 2, 1e-9));
        }
    }

}  // namespace
