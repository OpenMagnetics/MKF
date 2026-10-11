#include <source_location>
#include "support/Painter.h"
#include "physical_models/MagneticField.h"
#include "physical_models/Reluctance.h"
#include "json.hpp"
#include "TestingUtils.h"
#include "Models.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <fstream>
#include <map>
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

    // MMF_BOUNDARY: Ampere round the winding window. Just inside a mu -> inf wall the tangential H is the MMF the wall
    // carries, so the circulation along the column face over a gap's mouth is that gap's MMF and along the rest of the
    // walls the ferrite's. MMF_BOUNDARY counts every gap once: the mouth carries the gap's reluctance share R_gap / R_total
    // of NI (MKF's own reluctance model), the remaining walls the rest, the whole contour NI, and every wall exactly its
    // own sources. (ROSHEN puts Roshen's gap MMF 0.9 H_g g at the mouth ON TOP of a lattice whose walls already carry the
    // whole NI -- a double count. A count-once Roshen variant, mouth = 0.9 H_g g with H_g = NI / (R_total A_e mu0), was
    // tried and dropped: on real gaps 0.9 H_g g exceeds NI -- PQ 26/25 1 mm 1.054 NI, E42 1 mm 1.122 NI -- because
    // H_g = B / mu0 ignores the fringing that lowers the gap's reluctance.) Ground gaps 0.1 mm and 1 mm, and the
    // ungapped core (residual gaps), on PQ 26/25. Also: Roshen's Hx IS the mouth sheet's with its coincident image.
    TEST_CASE("Test_Magnetic_Field_Mmf_Boundary_Counts_Each_Gap_Once", "[physical-model][magnetic-field][mmf-boundary]") {
        for (double groundGap : {1e-4, 1e-3, 0.0}) {
            const bool gapped = groundGap > 0;
            numberTurns = {8};
            numberParallels = {1};
            turnsRatios = {};
            interleavingLevel = 1;
            numberStacks = 1;
            voltagePeakToPeak = 2000;
            frequency = 125000;
            gapping = gapped ? OpenMagneticsTesting::get_ground_gap(groundGap) : OpenMagneticsTesting::get_residual_gap();
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
            const double y0 = -B / 2;
            const double y1 = B / 2;

            // The main-column gap and its reluctance share of the MMF, from MKF's reluctance model.
            std::optional<CoreGap> mainGap;
            std::optional<size_t> mainGapIndex;
            auto gappingNow = core.get_functional_description().get_gapping();
            for (size_t gapIndex = 0; gapIndex < gappingNow.size(); ++gapIndex) {
                if (static_cast<size_t>(core.find_closest_column_index_by_coordinates(gappingNow[gapIndex].get_coordinates().value())) == core.get_main_column_index()) {
                    mainGap = gappingNow[gapIndex];
                    mainGapIndex = gapIndex;
                }
            }
            REQUIRE(mainGap);
            auto reluctanceModel = ReluctanceModel::factory();
            auto reluctanceOutput = reluctanceModel->get_core_reluctance(core, operatingPoint);
            const double totalReluctance = reluctanceOutput.get_core_reluctance();
            REQUIRE(reluctanceOutput.get_reluctance_per_gap());
            const double gapMmfPerAmpereTurn = reluctanceOutput.get_reluctance_per_gap().value()[mainGapIndex.value()].get_reluctance() / totalReluctance;
            const double gapLength = mainGap->get_length();
            const double gapY = mainGap->get_coordinates().value()[1];
            if (gapped) {
                CHECK(gapMmfPerAmpereTurn > 0.3);  // a gap whose MMF a double count would show
            }

            // Roshen's closed form = 0.9 x (mouth sheet of density H_g + its coincident image).
            for (std::vector<double> probe : {std::vector<double>{x0 + 0.3e-3, gapY + 0.2e-3}, {x0 + 2e-3, gapY - 1e-3}, {x0 + 0.05e-3, gapY}}) {
                FieldPoint point;
                point.set_point(probe);
                auto roshen = MagneticFieldStrengthRoshenModel().get_magnetic_field_strength_between_gap_and_point(mainGap.value(), 1.0, point);
                WindowWallSheet mouth{{x0, gapY - gapLength / 2}, {x0, gapY + gapLength / 2}, 0};
                auto [sheetX, sheetY] = MagneticFieldStrengthMmfBoundaryModel::calculate_sheet_field(mouth, -0.9 * gapLength, probe[0], probe[1]);
                double scale = hypot(roshen.get_real(), roshen.get_imaginary());
                CHECK_THAT(2 * sheetX, Catch::Matchers::WithinAbs(roshen.get_real(), 1e-9 * scale));
            }

            // Contour 5 um inside the walls, counter-clockwise, midpoint rule with steps of a quarter of that inset (the
            // sheets' ends make the field near the corners vary on the inset's scale). The mouth window reaches 0.5 mm past
            // the gap's ends: at 5 um from the face less than 1 % of the mouth's MMF lies beyond, and the column
            // ferrite's own density over the window is added to the expectation.
            const double inset = 5e-6;
            const double mouthHalfWindow = gapLength / 2 + 0.5e-3;
            std::vector<FieldPoint> points;
            std::vector<std::array<double, 2>> steps;
            std::vector<bool> onMouth;
            std::vector<size_t> wallOf;
            auto addWall = [&](double ax, double ay, double bx, double by, size_t wall) {
                const bool columnWall = wall == 3;
                size_t samples = static_cast<size_t>(std::ceil(hypot(bx - ax, by - ay) / (inset / 4)));
                for (size_t index = 0; index < samples; ++index) {
                    double t = (index + 0.5) / samples;
                    FieldPoint point;
                    point.set_point({ax + t * (bx - ax), ay + t * (by - ay)});
                    points.push_back(point);
                    steps.push_back({(bx - ax) / samples, (by - ay) / samples});
                    onMouth.push_back(columnWall && std::abs(ay + t * (by - ay) - gapY) < mouthHalfWindow);
                    wallOf.push_back(wall);
                }
            };
            addWall(x0 + inset, y0 + inset, x1 - inset, y0 + inset, 0);
            addWall(x1 - inset, y0 + inset, x1 - inset, y1 - inset, 1);
            addWall(x1 - inset, y1 - inset, x0 + inset, y1 - inset, 2);
            addWall(x0 + inset, y1 - inset, x0 + inset, y0 + inset, 3);

            // The window's own sources (each boundary segment once, at its largest image factor), per wall: bottom,
            // lateral, top, column. They must return exactly NI, and the column's mouth sheet must be the gap's reluctance share.
            std::map<std::array<long, 4>, double> originals;
            for (auto& sheet : CoilMesher::generate_window_wall_sheets(magnetic, operatingPoint)) {
                std::array<long, 4> key{std::lround(sheet.start[0] * 1e9), std::lround(sheet.start[1] * 1e9), std::lround(sheet.end[0] * 1e9), std::lround(sheet.end[1] * 1e9)};
                if (std::abs(sheet.current) > std::abs(originals[key])) {
                    originals[key] = sheet.current;
                }
            }
            std::array<double, 4> expectedPerWall{0, 0, 0, 0};
            double columnSheet = 0;
            double mouthSheet = 0;
            const long nx0 = std::lround(x0 * 1e9), nx1 = std::lround(x1 * 1e9), ny0 = std::lround(y0 * 1e9), ny1 = std::lround(y1 * 1e9);
            for (auto& [key, current] : originals) {
                const bool inWindowX = key[0] >= nx0 && key[2] <= nx1 && key[0] <= nx1 && key[2] >= nx0;
                const bool inWindowY = key[1] >= ny0 && key[3] <= ny1 && key[1] <= ny1 && key[3] >= ny0;
                if (key[0] == key[2] && inWindowY && (key[0] == nx0 || key[0] == nx1)) {
                    expectedPerWall[key[0] == nx0 ? 3 : 1] += current;
                    if (key[0] == nx0 && key[1] == ny0 && key[3] == ny1) {
                        columnSheet = current;
                    }
                    else if (key[0] == nx0) {
                        mouthSheet += current;
                    }
                }
                else if (key[1] == key[3] && inWindowX && (key[1] == ny0 || key[1] == ny1)) {
                    expectedPerWall[key[1] == ny0 ? 0 : 2] += current;
                }
            }
            CHECK_THAT(expectedPerWall[0] + expectedPerWall[1] + expectedPerWall[2] + expectedPerWall[3], Catch::Matchers::WithinAbs(-1.0, 1e-6));
            CHECK_THAT(-mouthSheet, Catch::Matchers::WithinRel(gapMmfPerAmpereTurn, 1e-9));
            REQUIRE(columnSheet != 0);
            const double expectedMouth = gapMmfPerAmpereTurn + std::abs(columnSheet) * 2 * mouthHalfWindow / B;

            Field inducedField;
            inducedField.set_data(points);
            inducedField.set_frequency(frequency);
            MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON, MagneticFieldStrengthFringingEffectModels::MMF_BOUNDARY);
            auto field = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField).get_field_per_frequency()[0];
            REQUIRE(field.get_data().size() == points.size());

            double harmonicAmplitude = operatingPoint.get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
            const double ampereTurns = harmonicAmplitude * numberTurns[0];
            double mouthCirculation = 0;
            double wallCirculation = 0;
            std::array<double, 4> circulationPerWall{0, 0, 0, 0};
            for (size_t index = 0; index < points.size(); ++index) {
                double contribution = field.get_data()[index].get_real() * steps[index][0] + field.get_data()[index].get_imaginary() * steps[index][1];
                (onMouth[index] ? mouthCirculation : wallCirculation) += contribution;
                circulationPerWall[wallOf[index]] += contribution;
            }
            // Total: NI. Mouth: the gap's reluctance share, once (plus the column ferrite over the window). Rest: what the gap leaves.
            // Every wall: exactly its own sources -- a double count cannot satisfy all four, since the total is NI.
            CHECK_THAT(std::abs(mouthCirculation + wallCirculation), Catch::Matchers::WithinRel(ampereTurns, 0.005));
            CHECK_THAT(std::abs(mouthCirculation), Catch::Matchers::WithinAbs(expectedMouth * ampereTurns, 0.01 * ampereTurns));
            CHECK_THAT(std::abs(wallCirculation), Catch::Matchers::WithinAbs((1 - expectedMouth) * ampereTurns, 0.01 * ampereTurns));
            for (size_t wall = 0; wall < 4; ++wall) {
                CHECK_THAT(std::abs(circulationPerWall[wall]), Catch::Matchers::WithinAbs(std::abs(expectedPerWall[wall]) * ampereTurns, 0.01 * ampereTurns));
            }
        }
        gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        settings.reset();
    }

    // ABT #1759: the EP family's back leg is displaced only in z (x = 0); it is the window's whole return limb, which the
    // 2-D frame represents by the outer wall x1. Its lateral (residual) gap at x = 0 used to make the wall sources throw
    // ("lateral gaps but none on the x > 0 side"), and with them the default parallel current split, the winding losses
    // and the Painter field plot of every EP with a paralleled winding (Test_Coil_Painter_Web_3). The back leg's gap now
    // opens on x1 with the lateral MMF share, the ferrite walls carry the rest, and the sources return exactly NI.
    TEST_CASE("Test_Magnetic_Field_Mmf_Boundary_Ep_Back_Leg_Gap_Opens_On_The_Outer_Wall", "[physical-model][magnetic-field][mmf-boundary]") {
        settings.reset();
        clear_databases();
        auto jsonPath = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_painter_web_3_1473.json");
        std::ifstream jsonFile(jsonPath);
        OpenMagnetics::Magnetic magnetic(json::parse(jsonFile));
        auto core = magnetic.get_core();
        if (!core.get_processed_description()) {
            core.process_data();
        }
        core.process_gap_or_throw();
        magnetic.set_core(core);
        REQUIRE(core.get_shape_family() == CoreShapeFamily::EP);
        auto columns = core.get_columns();
        const size_t mainColumnIndex = core.get_main_column_index();
        const double x0 = columns[mainColumnIndex].get_width() / 2;
        const double x1 = x0 + core.get_processed_description()->get_winding_windows()[0].get_width().value();
        // The fixture is the case: one lateral limb, at x = 0, with a residual gap.
        size_t lateralColumns = 0;
        for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
            if (columnIndex != mainColumnIndex) {
                ++lateralColumns;
                REQUIRE(columns[columnIndex].get_coordinates()[0] == 0);
            }
        }
        REQUIRE(lateralColumns == 1);

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 10, 0.5, 0, {32});
        auto operatingPoint = inputs.get_operating_point(0);
        std::vector<WindowWallSheet> sheets;
        REQUIRE_NOTHROW(sheets = CoilMesher::generate_window_wall_sheets(magnetic, operatingPoint));

        // The function itself throws unless the window's original sources return exactly NI; here: where they sit. Sheets
        // are filtered to the window's height (the lattice's y images fall outside it); along x the cell's own source
        // and its nearest x image coincide on each wall, with the same factor for both walls.
        const double halfHeight = core.get_processed_description()->get_winding_windows()[0].get_height().value() / 2;
        double backLegGapSheet = 0;
        double mainGapSheet = 0;
        double backLegGapOnColumnWall = 0;
        for (auto& sheet : sheets) {
            const double length = std::hypot(sheet.end[0] - sheet.start[0], sheet.end[1] - sheet.start[1]);
            if (std::abs(sheet.start[1]) > halfHeight + 1e-12 || std::abs(sheet.end[1]) > halfHeight + 1e-12) {
                continue;
            }
            const bool residualMouth = std::abs(length - 5e-6) < 1e-9;
            if (sheet.start[0] == x1 && sheet.end[0] == x1 && residualMouth) {
                backLegGapSheet += sheet.current;
            }
            if (sheet.start[0] == x0 && sheet.end[0] == x0 && residualMouth) {
                backLegGapOnColumnWall += sheet.current;
            }
            if (sheet.start[0] == x0 && sheet.end[0] == x0 && std::abs(length - 0.000414) < 1e-9) {
                mainGapSheet += sheet.current;
            }
        }
        CHECK(backLegGapSheet < 0);
        CHECK(backLegGapOnColumnWall == 0);
        CHECK(mainGapSheet < backLegGapSheet);   // the 0.414 mm ground gap drops far more MMF than the 5 um residual one

        // End to end: the default parallel current split (2-parallel secondary) builds on these sources.
        auto twoWindingInputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 10, 0.5, 0, {32});
        REQUIRE_NOTHROW(MagneticField::calculate_parallel_current_split(magnetic, twoWindingInputs.get_operating_point(0), 25));
        settings.reset();
    }

    // Roshen (IEEE Trans. Magn. 43(8), 2007, as reproduced in Electronics 10, 1683, eqs. 6-7), l_g the half gap, xi the
    // distance from the leg face, y along it: Hx = H_g / (2 pi) ln[(xi^2 + (y - l_g)^2) / (xi^2 + (y + l_g)^2)],
    // Hy = -H_g / pi [atan(2 xi l_g / (xi^2 + y^2 - l_g^2)) + m pi] -- together the field of the mouth as a uniform
    // current sheet of density H_g on the mu -> inf face plus its coincident image, curl-free in the window. MKF had
    // atan(xi l_g / ...) (the factor 2 lost): Hx was that sheet's, Hy was not any field's partner, and the window
    // circulation of the gap field alone was not zero. Checked: both components equal the sheet's (x 0.9); round a
    // closed contour 5 um inside the window walls the gap field circulates 0 in total and 0.9 H_g g along the face over
    // its mouth; the whole ROSHEN window field (lattice + gap) circulates exactly NI. PQ 26/25, 0.1 mm and 1 mm gaps.
    TEST_CASE("Test_Magnetic_Field_Roshen_Gap_Field_Is_Curl_Free", "[physical-model][magnetic-field]") {
        for (double groundGap : {1e-4, 1e-3}) {
            numberTurns = {8};
            numberParallels = {1};
            turnsRatios = {};
            interleavingLevel = 1;
            numberStacks = 1;
            voltagePeakToPeak = 2000;
            frequency = 125000;
            gapping = OpenMagneticsTesting::get_ground_gap(groundGap);
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
            const double y0 = -B / 2;
            const double y1 = B / 2;
            std::optional<CoreGap> mainGap;
            for (auto& gap : core.get_functional_description().get_gapping()) {
                if (static_cast<size_t>(core.find_closest_column_index_by_coordinates(gap.get_coordinates().value())) == core.get_main_column_index()) {
                    mainGap = gap;
                }
            }
            REQUIRE(mainGap);
            const double gapLength = mainGap->get_length();
            const double gapY = mainGap->get_coordinates().value()[1];
            MagneticFieldStrengthRoshenModel roshenModel;

            // Both components: 0.9 x (mouth sheet of density H_g = 1 + its coincident image).
            for (std::vector<double> probe : {std::vector<double>{x0 + 0.3e-3, gapY + 0.2e-3}, {x0 + 2e-3, gapY - 1e-3}, {x0 + 0.05e-3, gapY},
                                              {x0 + 1e-3, gapY + 0.4e-3}, {x0 + 0.2e-3, gapY - 0.6 * gapLength}}) {
                FieldPoint point;
                point.set_point(probe);
                auto roshen = roshenModel.get_magnetic_field_strength_between_gap_and_point(mainGap.value(), 1.0, point);
                WindowWallSheet mouth{{x0, gapY - gapLength / 2}, {x0, gapY + gapLength / 2}, 0};
                auto [sheetX, sheetY] = MagneticFieldStrengthMmfBoundaryModel::calculate_sheet_field(mouth, -0.9 * gapLength, probe[0], probe[1]);
                double scale = hypot(roshen.get_real(), roshen.get_imaginary());
                CHECK_THAT(2 * sheetX, Catch::Matchers::WithinAbs(roshen.get_real(), 1e-9 * scale));
                CHECK_THAT(2 * sheetY, Catch::Matchers::WithinAbs(roshen.get_imaginary(), 1e-9 * scale));
            }

            // Closed contour 5 um inside the walls, counter-clockwise, midpoint rule (steps of a quarter of the inset).
            const double inset = 5e-6;
            const double mouthHalfWindow = gapLength / 2 + 0.5e-3;
            std::vector<FieldPoint> points;
            std::vector<std::array<double, 2>> steps;
            std::vector<bool> onMouth;
            auto addWall = [&](double ax, double ay, double bx, double by, bool columnWall) {
                size_t samples = static_cast<size_t>(std::ceil(hypot(bx - ax, by - ay) / (inset / 4)));
                for (size_t index = 0; index < samples; ++index) {
                    double t = (index + 0.5) / samples;
                    FieldPoint point;
                    point.set_point({ax + t * (bx - ax), ay + t * (by - ay)});
                    points.push_back(point);
                    steps.push_back({(bx - ax) / samples, (by - ay) / samples});
                    onMouth.push_back(columnWall && std::abs(ay + t * (by - ay) - gapY) < mouthHalfWindow);
                }
            };
            addWall(x0 + inset, y0 + inset, x1 - inset, y0 + inset, false);
            addWall(x1 - inset, y0 + inset, x1 - inset, y1 - inset, false);
            addWall(x1 - inset, y1 - inset, x0 + inset, y1 - inset, false);
            addWall(x0 + inset, y1 - inset, x0 + inset, y0 + inset, true);

            // The gap field alone (H_g = 1): curl-free, so 0 round the contour, and 0.9 g along the face over its mouth.
            double gapTotal = 0;
            double gapMouth = 0;
            for (size_t index = 0; index < points.size(); ++index) {
                auto h = roshenModel.get_magnetic_field_strength_between_gap_and_point(mainGap.value(), 1.0, points[index]);
                double contribution = h.get_real() * steps[index][0] + h.get_imaginary() * steps[index][1];
                gapTotal += contribution;
                if (onMouth[index]) {
                    gapMouth += contribution;
                }
            }
            CHECK_THAT(gapTotal, Catch::Matchers::WithinAbs(0, 0.002 * 0.9 * gapLength));
            CHECK_THAT(std::abs(gapMouth), Catch::Matchers::WithinRel(0.9 * gapLength, 0.02));

            // The whole ROSHEN window field: the lattice carries NI and the curl-free gap field adds nothing round the contour.
            Field inducedField;
            inducedField.set_data(points);
            inducedField.set_frequency(frequency);
            MagneticField magneticField(MagneticFieldStrengthModels::BINNS_LAWRENSON, MagneticFieldStrengthFringingEffectModels::ROSHEN);
            auto field = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField).get_field_per_frequency()[0];
            REQUIRE(field.get_data().size() == points.size());
            double harmonicAmplitude = operatingPoint.get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
            const double ampereTurns = harmonicAmplitude * numberTurns[0];
            double windowCirculation = 0;
            for (size_t index = 0; index < points.size(); ++index) {
                windowCirculation += field.get_data()[index].get_real() * steps[index][0] + field.get_data()[index].get_imaginary() * steps[index][1];
            }
            CHECK_THAT(std::abs(windowCirculation), Catch::Matchers::WithinRel(ampereTurns, 0.005));
        }
        gapping = OpenMagneticsTesting::get_ground_gap(0.001);
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
