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
            Painter painter(outFile);
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
        // Disable core-wall image currents so the test checks the free-space two-wire field against
        // the closed-form -2*I/(2*pi*d), exactly like the Binns-Lawrenson sibling above. Without this
        // the default mirroring adds image contributions (~1.1x) and the closed form no longer applies.
        settings.set_magnetic_field_mirroring_dimension(0);
        settings.set_magnetic_field_include_fringing(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double distanceCenterPoint = fieldPoint.get_point()[0] - turn_0.get_coordinates()[0];
        double expectedValue = -2 * harmonicAmplitude / (2 * M_PI * distanceCenterPoint);

        // July 2026: after the Lammeraner atan2(dy,dx) fix the field is correctly AZIMUTHAL — for two
        // horizontally-separated turns probed on their shared axis it lands in the imaginary (Y)
        // component, like the Binns-Lawrenson sibling above. The old swapped atan2 rotated it 90 deg
        // into the real (X) component, which this test used to check.
        REQUIRE((expectedValue - field.get_data()[0].get_imaginary()) / expectedValue < maximumError);
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
        Painter painter(outFile);
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
            Painter painter(outFile);
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
            Painter painter(outFile);
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


    // ABT #376: DOWELL's field model. It was declared in MagneticFieldStrengthModels but the
    // factory had no branch for it, so selecting it threw MODEL_NOT_AVAILABLE — and the
    // model-comparison sweeps that named it failed by the thousand.
    //
    // What defines Dowell is the MMF STAIRCASE: layers are assumed to span the winding breadth b,
    // so the field is parallel to them and set only by the ampere-turns enclosed — H = MMF / b,
    // with no distance dependence at all. These assertions pin exactly that, because it is what
    // separates his model from the 2D ones next to it.
    TEST_CASE("Test_Magnetic_Field_Dowell_Staircase", "[physical-model][magnetic-field][dowell]") {
        settings.reset();
        auto model = MagneticField::factory(MagneticFieldStrengthModels::DOWELL);
        REQUIRE(model != nullptr);

        double windingBreadth = 0.02;
        model->_windingWindowBreadth = windingBreadth;

        double current = 3.0;
        FieldPoint inducingFieldPoint;
        inducingFieldPoint.set_point({0.010, 0.0});
        inducingFieldPoint.set_value(current);

        auto fieldAt = [&](double radialPosition, double axialPosition) {
            FieldPoint inducedFieldPoint;
            inducedFieldPoint.set_point({radialPosition, axialPosition});
            return model->get_magnetic_field_strength_between_two_points(inducingFieldPoint, inducedFieldPoint);
        };

        // One conductor is one step of height I/b, felt only OUTSIDE it.
        double expectedStep = current / windingBreadth;
        CHECK_THAT(fieldAt(0.012, 0.0).get_imaginary(), Catch::Matchers::WithinRel(expectedStep, 1e-9));
        CHECK(fieldAt(0.008, 0.0).get_imaginary() == 0);
        // On the conductor itself, half the step — Dowell's treatment of the layer carrying it.
        CHECK_THAT(fieldAt(0.010, 0.0).get_imaginary(), Catch::Matchers::WithinRel(expectedStep / 2, 1e-9));

        // The field is PARALLEL to the layers: no radial component anywhere.
        CHECK(fieldAt(0.012, 0.0).get_real() == 0);
        CHECK(fieldAt(0.030, 0.005).get_real() == 0);

        // NO distance dependence — the property that makes this a 1D model. Doubling the
        // separation, or moving along the layer, leaves the step untouched (a 2D model would
        // fall off with distance; that difference is the point of offering both).
        CHECK_THAT(fieldAt(0.014, 0.0).get_imaginary(), Catch::Matchers::WithinRel(expectedStep, 1e-9));
        CHECK_THAT(fieldAt(0.012, 0.009).get_imaginary(), Catch::Matchers::WithinRel(expectedStep, 1e-9));

        // Staircase: N conductors enclosed give N steps, i.e. the enclosed ampere-turns over b.
        double accumulatedField = 0;
        size_t numberEnclosedTurns = 5;
        for (size_t turnIndex = 0; turnIndex < numberEnclosedTurns; ++turnIndex) {
            FieldPoint turnFieldPoint;
            turnFieldPoint.set_point({0.006 + 0.001 * double(turnIndex), 0.0});
            turnFieldPoint.set_value(current);
            FieldPoint inducedFieldPoint;
            inducedFieldPoint.set_point({0.015, 0.0});
            accumulatedField += model->get_magnetic_field_strength_between_two_points(turnFieldPoint, inducedFieldPoint).get_imaginary();
        }
        CHECK_THAT(accumulatedField,
                   Catch::Matchers::WithinRel(double(numberEnclosedTurns) * current / windingBreadth, 1e-9));

        // Without a breadth there is no Dowell field, and the model says so rather than guessing.
        auto modelWithoutBreadth = MagneticField::factory(MagneticFieldStrengthModels::DOWELL);
        FieldPoint anyInducedPoint;
        anyInducedPoint.set_point({0.012, 0.0});
        CHECK_THROWS(modelWithoutBreadth->get_magnetic_field_strength_between_two_points(inducingFieldPoint, anyInducedPoint));
        settings.reset();
    }

    // ABT #1304: the closed-form wall-sheet kernel against the thing it integrates -- 20000 Binns-Lawrenson filaments
    // along the same segment -- beside the sheet, beyond its ends, off its line, and 1 um from it.
    TEST_CASE("Test_Magnetic_Field_Mmf_Boundary_Sheet_Kernel_Matches_Filaments", "[physical-model][magnetic-field][mmf-boundary]") {
        const double length = 1e-3;
        WindowWallSheet sheet{{0.0, 0.0}, {0.0, length}, 1.0};
        const size_t numberFilaments = 20000;
        MagneticFieldStrengthBinnsLawrensonModel filamentModel;
        std::vector<std::vector<double>> probes = {{0.3e-3, 0.5e-3}, {-0.05e-3, 0.2e-3}, {0.2e-3, 1.5e-3},
                                                   {1e-6, 0.5e-3}, {-1e-6, 0.5e-3}, {0.0, 1.3e-3}, {2e-3, -1e-3}};
        for (auto& probe : probes) {
            double filamentX = 0;
            double filamentY = 0;
            FieldPoint induced;
            induced.set_point(probe);
            for (size_t index = 0; index < numberFilaments; ++index) {
                FieldPoint filament;
                filament.set_point({0.0, (index + 0.5) * length / numberFilaments});
                filament.set_value(1.0 / numberFilaments);
                auto [fieldX, fieldY] = filamentModel.get_magnetic_field_strength_components_between_two_points(filament, induced);
                filamentX += fieldX;
                filamentY += fieldY;
            }
            auto [sheetX, sheetY] = MagneticFieldStrengthMmfBoundaryModel::calculate_sheet_field(sheet, 1.0, probe[0], probe[1]);
            double scale = hypot(filamentX, filamentY);
            CHECK_THAT(sheetX, Catch::Matchers::WithinAbs(filamentX, 1e-3 * scale));
            CHECK_THAT(sheetY, Catch::Matchers::WithinAbs(filamentY, 1e-3 * scale));
        }
        // Across the sheet the tangential field jumps by the sheet density I / L, as a current sheet must.
        auto [insideX, insideY] = MagneticFieldStrengthMmfBoundaryModel::calculate_sheet_field(sheet, 1.0, 1e-9, 0.5e-3);
        auto [outsideX, outsideY] = MagneticFieldStrengthMmfBoundaryModel::calculate_sheet_field(sheet, 1.0, -1e-9, 0.5e-3);
        CHECK_THAT(std::abs(insideY - outsideY), Catch::Matchers::WithinRel(1.0 / length, 1e-5));
    }

    // ABT #1304: a mu -> inf wall with a ferrite MMF drop along it has tangential H equal to that drop's density. The
    // MMF_BOUNDARY puts exactly that on each wall; the historical lattice (net current, no MMF return) cannot, and
    // this test fails with the ROSHEN fringing model in its place. Checked on the column and lateral walls, away from the gap.
    TEST_CASE("Test_Magnetic_Field_Mmf_Boundary_Honours_The_Core_Walls", "[physical-model][magnetic-field][mmf-boundary]") {
        numberTurns = {8};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 1;
        numberStacks = 1;
        voltagePeakToPeak = 2000;
        frequency = 125000;
        setup();
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        auto operatingPoint = inputs.get_operating_point(0);

        settings.reset();
        settings.set_magnetic_field_mirroring_dimension(1);

        auto processed = core.get_processed_description().value();
        auto window = processed.get_winding_windows()[0];
        const double A = window.get_width().value();
        const double B = window.get_height().value();
        const double x0 = core.get_columns()[core.get_main_column_index()].get_width() / 2;
        const double x1 = x0 + A;

        // The ferrite sheets of the two vertical walls (the originals: full-height, largest image factor).
        auto sheets = CoilMesher::generate_window_wall_sheets(magnetic, operatingPoint);
        double columnSheet = 0;
        double lateralSheet = 0;
        for (auto& sheet : sheets) {
            bool fullHeight = std::abs(sheet.start[1] + B / 2) < 1e-12 && std::abs(sheet.end[1] - B / 2) < 1e-12;
            if (fullHeight && sheet.start[0] == x0 && sheet.end[0] == x0 && std::abs(sheet.current) > std::abs(columnSheet)) {
                columnSheet = sheet.current;
            }
            if (fullHeight && sheet.start[0] == x1 && sheet.end[0] == x1 && std::abs(sheet.current) > std::abs(lateralSheet)) {
                lateralSheet = sheet.current;
            }
        }
        REQUIRE(columnSheet != 0);
        REQUIRE(lateralSheet != 0);

        double harmonicAmplitude = operatingPoint.get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        const double ampereTurns = harmonicAmplitude * numberTurns[0];
        const double windowFieldScale = ampereTurns / B;

        std::vector<FieldPoint> points;
        for (double fraction : {-0.4, -0.3, 0.3, 0.4}) {   // the 1 mm ground gap sits at y = 0
            FieldPoint point;
            point.set_point({x0 + 1e-6, fraction * B});
            points.push_back(point);
            point.set_point({x1 - 1e-6, fraction * B});
            points.push_back(point);
        }
        Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);
        MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON, MagneticFieldStrengthFringingEffectModels::MMF_BOUNDARY);
        auto field = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField).get_field_per_frequency()[0];

        for (size_t index = 0; index < field.get_data().size(); ++index) {
            bool onColumn = index % 2 == 0;
            double expectedTangential = std::abs((onColumn ? columnSheet : lateralSheet) * ampereTurns / B);
            CHECK_THAT(std::abs(field.get_data()[index].get_imaginary()), Catch::Matchers::WithinAbs(expectedTangential, 0.03 * windowFieldScale));
        }
        settings.reset();
    }

    // ABT #1304: with MMF_BOUNDARY every lattice cell is neutral, so the image sum converges -- one and two mirror cells agree.
    TEST_CASE("Test_Magnetic_Field_Mmf_Boundary_Lattice_Converges", "[physical-model][magnetic-field][mmf-boundary]") {
        numberTurns = {8};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 1;
        numberStacks = 1;
        voltagePeakToPeak = 2000;
        frequency = 125000;
        setup();
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        auto operatingPoint = inputs.get_operating_point(0);

        auto fieldWithCells = [&](int cells) {
            settings.reset();
            settings.set_magnetic_field_mirroring_dimension(cells);
            MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON, MagneticFieldStrengthFringingEffectModels::MMF_BOUNDARY);
            return magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic).get_field_per_frequency()[0];
        };
        auto one = fieldWithCells(1);
        auto two = fieldWithCells(2);
        REQUIRE(one.get_data().size() == two.get_data().size());
        double scale = 0;
        for (auto& point : two.get_data()) {
            scale = std::max(scale, hypot(point.get_real(), point.get_imaginary()));
        }
        for (size_t index = 0; index < one.get_data().size(); ++index) {
            CHECK_THAT(one.get_data()[index].get_real(), Catch::Matchers::WithinAbs(two.get_data()[index].get_real(), 0.01 * scale));
            CHECK_THAT(one.get_data()[index].get_imaginary(), Catch::Matchers::WithinAbs(two.get_data()[index].get_imaginary(), 0.01 * scale));
        }
        settings.reset();
    }

    // ABT #1304: MMF_BOUNDARY is refused -- not quietly replaced by another model -- where its sources cannot be consistent.
    TEST_CASE("Test_Magnetic_Field_Mmf_Boundary_Refuses_What_It_Cannot_Model", "[physical-model][magnetic-field][mmf-boundary]") {
        numberTurns = {8};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 1;
        numberStacks = 1;
        voltagePeakToPeak = 2000;
        frequency = 125000;
        setup();
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        auto operatingPoint = inputs.get_operating_point(0);
        settings.reset();
        for (auto fieldModel : {MagneticFieldStrengthModels::LAMMERANER, MagneticFieldStrengthModels::DOWELL, MagneticFieldStrengthModels::WANG, MagneticFieldStrengthModels::ALBACH}) {
            MagneticField magneticField(fieldModel, MagneticFieldStrengthFringingEffectModels::MMF_BOUNDARY);
            CHECK_THROWS_AS(magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic), ModelNotAvailableException);
        }
        settings.set_magnetic_field_include_fringing(false);
        MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON, MagneticFieldStrengthFringingEffectModels::MMF_BOUNDARY);
        CHECK_THROWS_AS(magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic), InvalidInputException);
        settings.reset();
    }

}  // namespace
