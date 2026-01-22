#include <source_location>
#include "support/Painter.h"
#include "physical_models/MagneticField.h"
#include "json.hpp"
#include "TestingUtils.h"
#include "Models.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fstream>
#include <string>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace MAS;
using namespace OpenMagnetics;

namespace { 
    double maximumError = 0.05;
    auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
    bool plot = false;

    std::vector<int64_t> numberTurns = {1};
    std::vector<int64_t> numberParallels = {1};
    std::vector<double> turnsRatios = {};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    double voltagePeakToPeak = 2000;
    double frequency = 125000;
    std::string coreShape = "PQ 26/25";
    std::string coreMaterial = "3C97";
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
    WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
    WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
    CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
    CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

    OpenMagnetics::Coil coil;
    Core core;
    OpenMagnetics::Inputs inputs;


    void setup() {
        coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        inputs = OpenMagnetics::Inputs::create_quick_operating_point(frequency, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();
    }


    TEST_CASE("Test_Magnetic_Field_Frequencies", "[physical-model][magnetic-field][smoke-test]") {
        numberTurns = {2};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 2;
        settings.set_harmonic_amplitude_threshold(0.01);

        sectionsAlignment = CoilAlignment::SPREAD;
        turnsAlignment = CoilAlignment::CENTERED;
        setup();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON);
        settings.set_magnetic_field_include_fringing(false);

        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic);
        auto field_0 = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
        auto field_1 = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[1];
        auto field_2 = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[2];

        REQUIRE_THAT(frequency, Catch::Matchers::WithinAbs(windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0].get_frequency(), maximumError));
        REQUIRE_THAT(3* frequency, Catch::Matchers::WithinAbs(windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[1].get_frequency(), maximumError));
        REQUIRE_THAT(5 * frequency, Catch::Matchers::WithinAbs(windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[2].get_frequency(), maximumError));

        REQUIRE(windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[1].get_data()[0].get_imaginary() <  windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0].get_data()[0].get_imaginary());
        REQUIRE(windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[2].get_data()[0].get_imaginary() <  windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0].get_data()[0].get_imaginary());
        settings.reset();
    }

    TEST_CASE("Test_Magnetic_Field_One_Turn_Round", "[physical-model][magnetic-field][smoke-test]") {

        numberTurns = {1};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 1;
        numberStacks = 1;
        voltagePeakToPeak = 2000;
        frequency = 125000;
        setup();
        auto turn = coil.get_turns_description().value()[0];

        std::vector<FieldPoint> points;
        FieldPoint fieldPoint;
        double maximumWidth = coil.resolve_wire(0).get_maximum_outer_width();
        double maximumHeight = coil.resolve_wire(0).get_maximum_outer_height();
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] - (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1]});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] + (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1]});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0], turn.get_coordinates()[1] - (maximumHeight / 2) * 1.0001});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0], turn.get_coordinates()[1] + (maximumHeight / 2) * 1.0001});
        points.push_back(fieldPoint);
        Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON);
        settings.set_magnetic_field_mirroring_dimension(0);
        settings.set_magnetic_field_include_fringing(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double expectedValue = harmonicAmplitude / (2 * M_PI * maximumWidth / 2);

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(field.get_data()[0].get_imaginary(), expectedValue * maximumError));
        REQUIRE_THAT(field.get_data()[0].get_real(), Catch::Matchers::WithinAbs(-field.get_data()[1].get_real(), expectedValue * maximumError));
        REQUIRE_THAT(field.get_data()[0].get_imaginary(), Catch::Matchers::WithinAbs(-field.get_data()[1].get_imaginary(), expectedValue * maximumError));
        REQUIRE_THAT(field.get_data()[0].get_real(), Catch::Matchers::WithinAbs(-field.get_data()[2].get_imaginary(), expectedValue * maximumError));
        REQUIRE_THAT(field.get_data()[0].get_imaginary(), Catch::Matchers::WithinAbs(-field.get_data()[2].get_real(), expectedValue * maximumError));
        REQUIRE_THAT(field.get_data()[0].get_real(), Catch::Matchers::WithinAbs(field.get_data()[3].get_imaginary(), expectedValue * maximumError));
        REQUIRE_THAT(field.get_data()[0].get_imaginary(), Catch::Matchers::WithinAbs(field.get_data()[3].get_real(), expectedValue * maximumError));
        if (plot) {

            settings.set_painter_mode(PainterModes::QUIVER);

            settings.set_magnetic_field_number_points_x(50);
            settings.set_magnetic_field_number_points_y(50);

            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Magnetic_Field_One_Turn_Round.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

 
        settings.reset();
    }

    TEST_CASE("Test_Magnetic_Field_Two_Turns_Round_Same_Current", "[physical-model][magnetic-field][smoke-test]") {
        numberTurns = {2};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 2;

        sectionsAlignment = CoilAlignment::SPREAD;
        turnsAlignment = CoilAlignment::CENTERED;
        setup();

        auto turn_0 = coil.get_turns_description().value()[0];
        auto turn_1 = coil.get_turns_description().value()[1];

        std::vector<FieldPoint> points;
        FieldPoint fieldPoint;
        fieldPoint.set_point(std::vector<double>{(turn_0.get_coordinates()[0] + turn_1.get_coordinates()[0]) / 2, turn_0.get_coordinates()[1]});
        points.push_back(fieldPoint);
        Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON);
        // MagneticField magneticField(MagneticFieldStrengthModels::LAMMERANER);
        settings.set_magnetic_field_include_fringing(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double expectedValue = 0;

        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(field.get_data()[0].get_real(), maximumError));

        settings.reset();
    }

    TEST_CASE("Test_Magnetic_Field_Two_Turns_Round_Opposite_Current", "[physical-model][magnetic-field][smoke-test]") {
        numberTurns = {1, 1};
        numberParallels = {1, 1};
        turnsRatios = {1};
        interleavingLevel = 1;

        sectionsAlignment = CoilAlignment::SPREAD;
        turnsAlignment = CoilAlignment::CENTERED;
        setup();

        auto turn_0 = coil.get_turns_description().value()[0];
        auto turn_1 = coil.get_turns_description().value()[1];

        std::vector<FieldPoint> points;
        FieldPoint fieldPoint;
        fieldPoint.set_point(std::vector<double>{(turn_0.get_coordinates()[0] + turn_1.get_coordinates()[0]) / 2, turn_0.get_coordinates()[1]});
        points.push_back(fieldPoint);
        Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON);
        settings.set_magnetic_field_mirroring_dimension(0);
        settings.set_magnetic_field_include_fringing(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double distanceCenterPoint = fieldPoint.get_point()[0] - turn_0.get_coordinates()[0];
        double expectedValue = -2 * harmonicAmplitude / (2 * M_PI * distanceCenterPoint);
        REQUIRE((expectedValue - field.get_data()[0].get_imaginary()) / expectedValue < maximumError);
        settings.reset();
    }

    TEST_CASE("Test_Magnetic_Field_Two_Turns_Round_Opposite_Current_Lammeraner", "[physical-model][magnetic-field][smoke-test]") {
        numberTurns = {1, 1};
        numberParallels = {1, 1};
        turnsRatios = {1};
        interleavingLevel = 1;

        sectionsAlignment = CoilAlignment::SPREAD;
        turnsAlignment = CoilAlignment::CENTERED;
        setup();

        auto turn_0 = coil.get_turns_description().value()[0];
        auto turn_1 = coil.get_turns_description().value()[1];

        std::vector<FieldPoint> points;
        FieldPoint fieldPoint;
        fieldPoint.set_point(std::vector<double>{(turn_0.get_coordinates()[0] + turn_1.get_coordinates()[0]) / 2, turn_0.get_coordinates()[1]});
        points.push_back(fieldPoint);
        Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);


        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        MagneticField magneticField(MagneticFieldStrengthModels::LAMMERANER);
        settings.set_magnetic_field_include_fringing(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double distanceCenterPoint = fieldPoint.get_point()[0] - turn_0.get_coordinates()[0];
        double expectedValue = -2 * harmonicAmplitude / (2 * M_PI * distanceCenterPoint);

        REQUIRE((expectedValue - field.get_data()[0].get_real()) / expectedValue < maximumError);
        settings.reset();
    }

    TEST_CASE("Test_Magnetic_Image_Method", "[physical-model][magnetic-field][smoke-test]") {
        gapping = OpenMagneticsTesting::get_residual_gap();
        coreShape = "P 9/5";
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        setup();

        WindingWindowElement windingWindow = core.get_processed_description().value().get_winding_windows()[0];

        double coreColumnHeight = core.get_columns()[0].get_height();

        auto turn_0 = coil.get_turns_description().value()[0];
        std::vector<FieldPoint> points;
        FieldPoint fieldPoint;
        fieldPoint.set_point(std::vector<double>{turn_0.get_coordinates()[0], coreColumnHeight / 2});
        points.push_back(fieldPoint);
        Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);


        OpenMagnetics::Wire wire = find_wire_by_name("Round 0.475 - Grade 1");
        auto wires = std::vector<OpenMagnetics::Wire>({wire});
        coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        coil.delimit_and_compact();
        auto turn = coil.get_turns_description().value()[0];

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON);
        settings.set_magnetic_field_mirroring_dimension(1);
        settings.set_magnetic_field_include_fringing(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
 
        auto outFile = outputFilePath;
        outFile.append("Test_Magnetic_Image_Method.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_mirroring_dimension(1);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(field.get_data().size() == 1);
        settings.reset();
    }


    TEST_CASE("Test_Magnetic_Field_One_Turn_Rectangular", "[physical-model][magnetic-field][smoke-test]") {
        numberTurns = {1};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 1;
        numberStacks = 1;
        voltagePeakToPeak = 2000;
        frequency = 125000;
        coreShape = "PQ 26/25";
        coreMaterial = "3C97";
        gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        sectionOrientation = WindingOrientation::OVERLAPPING;
        layersOrientation = WindingOrientation::OVERLAPPING;
        sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        turnsAlignment = CoilAlignment::CENTERED;
        setup();


        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0028);
        wire.set_nominal_value_conducting_height(0.00076);
        wire.set_nominal_value_outer_height(0.0007676);
        wire.set_nominal_value_outer_width(0.002838);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);
        coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);

        coil.delimit_and_compact();

        auto turn = coil.get_turns_description().value()[0];

        std::vector<FieldPoint> points;
        FieldPoint fieldPoint;
        double maximumWidth = coil.resolve_wire(0).get_maximum_outer_width();
        double maximumHeight = coil.resolve_wire(0).get_maximum_outer_height();
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] - (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1]});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] + (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1]});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0], turn.get_coordinates()[1] - (maximumHeight / 2) * 1.0001});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0], turn.get_coordinates()[1] + (maximumHeight / 2) * 1.0001});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] + (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1] + 0.00001});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] + (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1] - 0.00001});
        points.push_back(fieldPoint);
        Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON);
        settings.set_magnetic_field_include_fringing(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double expectedValue = harmonicAmplitude / (2 * M_PI * maximumWidth / 2);

        REQUIRE_THAT(field.get_data()[4].get_real(), Catch::Matchers::WithinAbs(-field.get_data()[5].get_real(), expectedValue * maximumError));
        REQUIRE_THAT(field.get_data()[4].get_imaginary(), Catch::Matchers::WithinAbs(field.get_data()[5].get_imaginary(), expectedValue * maximumError));
        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Magnetic_Field_One_Turn_Rectangular.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            settings.set_painter_mode(PainterModes::QUIVER);
            settings.set_painter_logarithmic_scale(false);
            settings.set_painter_mirroring_dimension(0);
            settings.set_painter_include_fringing(false);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            settings.reset();
        }
    }


    TEST_CASE("Test_Magnetic_Field_One_Turn_Foil", "[physical-model][magnetic-field][smoke-test]") {
        numberTurns = {1};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 1;
        numberStacks = 1;
        voltagePeakToPeak = 2000;
        frequency = 125000;
        coreShape = "PQ 26/25";
        coreMaterial = "3C97";
        gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        sectionOrientation = WindingOrientation::OVERLAPPING;
        layersOrientation = WindingOrientation::OVERLAPPING;
        sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        turnsAlignment = CoilAlignment::CENTERED;
        setup();


        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.00001);
        wire.set_nominal_value_conducting_height(0.0076);
        wire.set_nominal_value_outer_height(0.007676);
        wire.set_nominal_value_outer_width(0.000011);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::FOIL);
        wires.push_back(wire);
        coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);

        coil.delimit_and_compact();

        auto turn = coil.get_turns_description().value()[0];

        std::vector<FieldPoint> points;
        FieldPoint fieldPoint;
        double maximumWidth = coil.resolve_wire(0).get_maximum_outer_width();
        double maximumHeight = coil.resolve_wire(0).get_maximum_outer_height();
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] - (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1]});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] + (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1]});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0], turn.get_coordinates()[1] - (maximumHeight / 2) * 1.0001});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0], turn.get_coordinates()[1] + (maximumHeight / 2) * 1.0001});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] + (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1] + 0.00001});
        points.push_back(fieldPoint);
        fieldPoint.set_point(std::vector<double>{turn.get_coordinates()[0] + (maximumWidth / 2) * 1.0001, turn.get_coordinates()[1] - 0.00001});
        points.push_back(fieldPoint);
        Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON);
        settings.set_magnetic_field_include_fringing(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double expectedValue = harmonicAmplitude / (2 * M_PI * maximumWidth / 2);

        REQUIRE_THAT(field.get_data()[4].get_real(), Catch::Matchers::WithinAbs(-field.get_data()[5].get_real(), expectedValue * maximumError));
        REQUIRE_THAT(field.get_data()[4].get_imaginary(), Catch::Matchers::WithinAbs(field.get_data()[5].get_imaginary(), expectedValue * maximumError));
        settings.reset();
        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Magnetic_Field_One_Turn_Foil.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            settings.set_painter_mode(PainterModes::QUIVER);
            settings.set_painter_logarithmic_scale(false);
            settings.set_painter_mirroring_dimension(0);
            settings.set_painter_include_fringing(false);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            settings.reset();
        }
    }

}  // namespace
