#include "Painter.h"
#include "MagneticField.h"
#include "json.hpp"
#include "TestingUtils.h"
#include "Models.h"

#include <UnitTest++.h>

#include <fstream>
#include <string>
#include <matplot/matplot.h>


SUITE(MagneticField) {
    std::string filePath = __FILE__;
    double maximumError = 0.05;
    auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");

    std::vector<int64_t> numberTurns = {1};
    std::vector<int64_t> numberParallels = {1};
    std::vector<double> turnsRatios = {};
    uint8_t interleavingLevel = 1;
    int64_t numberStacks = 1;
    double voltagePeakToPeak = 2000;
    double frequency = 125000;
    std::string coreShape = "PQ 26/25";
    std::string coreMaterial = "3C97";
    auto gapping = OpenMagneticsTesting::get_grinded_gap(0.001);
    OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
    OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
    OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
    OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

    OpenMagnetics::CoilWrapper coil;
    OpenMagnetics::CoreWrapper core;
    OpenMagnetics::InputsWrapper inputs;


    void setup() {
        coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(frequency, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();
    }


    TEST(Test_Magnetic_Field_Frequencies) {
        numberTurns = {2};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 2;

        sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        setup();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::MagneticField magneticField(OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON);
        magneticField.set_fringing_effect(false);
        magneticField.set_winding_losses_harmonic_amplitude_threshold(0.01);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic);
        auto field_0 = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
        auto field_1 = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[1];
        auto field_2 = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[2];

        CHECK_CLOSE(frequency, windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0].get_frequency(), maximumError);
        CHECK_CLOSE(3* frequency, windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[1].get_frequency(), maximumError);
        CHECK_CLOSE(5 * frequency, windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[2].get_frequency(), maximumError);

        CHECK(windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[1].get_data()[0].get_real() <  windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0].get_data()[0].get_real());
        CHECK(windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[2].get_data()[0].get_real() <  windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0].get_data()[0].get_real());

    }

    TEST(Test_Magnetic_Field_One_Turn_Round) {

        numberTurns = {1};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 1;
        numberStacks = 1;
        voltagePeakToPeak = 2000;
        frequency = 125000;
        setup();
        auto turn = coil.get_turns_description().value()[0];

        std::vector<OpenMagnetics::FieldPoint> points;
        OpenMagnetics::FieldPoint fieldPoint;
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
        OpenMagnetics::Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::MagneticField magneticField(OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON);
        magneticField.set_mirroring_dimension(0);
        magneticField.set_fringing_effect(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double expectedValue = harmonicAmplitude / (2 * std::numbers::pi * maximumWidth / 2);

        CHECK_CLOSE(expectedValue, field.get_data()[0].get_real(), expectedValue * maximumError);
        CHECK_CLOSE(field.get_data()[0].get_real(), -field.get_data()[1].get_real(), expectedValue * maximumError);
        CHECK_CLOSE(field.get_data()[0].get_imaginary(), field.get_data()[1].get_imaginary(), expectedValue * maximumError);
        CHECK_CLOSE(field.get_data()[0].get_real(), -field.get_data()[2].get_imaginary(), expectedValue * maximumError);
        CHECK_CLOSE(field.get_data()[0].get_imaginary(), field.get_data()[2].get_real(), expectedValue * maximumError);
        CHECK_CLOSE(field.get_data()[0].get_real(), field.get_data()[3].get_imaginary(), expectedValue * maximumError);
        CHECK_CLOSE(field.get_data()[0].get_imaginary(), field.get_data()[3].get_real(), expectedValue * maximumError);
    }

    TEST(Test_Magnetic_Field_Two_Turns_Round_Same_Current) {
        numberTurns = {2};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 2;

        sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        setup();

        auto turn_0 = coil.get_turns_description().value()[0];
        auto turn_1 = coil.get_turns_description().value()[1];

        std::vector<OpenMagnetics::FieldPoint> points;
        OpenMagnetics::FieldPoint fieldPoint;
        fieldPoint.set_point(std::vector<double>{(turn_0.get_coordinates()[0] + turn_1.get_coordinates()[0]) / 2, turn_0.get_coordinates()[1]});
        points.push_back(fieldPoint);
        OpenMagnetics::Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::MagneticField magneticField(OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON);
        // OpenMagnetics::MagneticField magneticField(OpenMagnetics::MagneticFieldStrengthModels::LAMMERANER);
        magneticField.set_fringing_effect(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double expectedValue = 0;

        CHECK_CLOSE(expectedValue, field.get_data()[0].get_real(), maximumError);

    }

    TEST(Test_Magnetic_Field_Two_Turns_Round_Opposite_Current) {
        numberTurns = {1, 1};
        numberParallels = {1, 1};
        turnsRatios = {1};
        interleavingLevel = 1;

        sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        setup();

        auto turn_0 = coil.get_turns_description().value()[0];
        auto turn_1 = coil.get_turns_description().value()[1];

        std::vector<OpenMagnetics::FieldPoint> points;
        OpenMagnetics::FieldPoint fieldPoint;
        fieldPoint.set_point(std::vector<double>{(turn_0.get_coordinates()[0] + turn_1.get_coordinates()[0]) / 2, turn_0.get_coordinates()[1]});
        points.push_back(fieldPoint);
        OpenMagnetics::Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::MagneticField magneticField(OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON);
        magneticField.set_fringing_effect(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double distanceCenterPoint = fieldPoint.get_point()[0] - turn_0.get_coordinates()[0];
        double expectedValue = -2 * harmonicAmplitude / (2 * std::numbers::pi * distanceCenterPoint);

        CHECK((expectedValue - field.get_data()[0].get_real()) / expectedValue < maximumError);
    }

    TEST(Test_Magnetic_Field_Two_Turns_Round_Opposite_Current_Lammeraner) {
        numberTurns = {1, 1};
        numberParallels = {1, 1};
        turnsRatios = {1};
        interleavingLevel = 1;

        sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        setup();

        auto turn_0 = coil.get_turns_description().value()[0];
        auto turn_1 = coil.get_turns_description().value()[1];

        std::vector<OpenMagnetics::FieldPoint> points;
        OpenMagnetics::FieldPoint fieldPoint;
        fieldPoint.set_point(std::vector<double>{(turn_0.get_coordinates()[0] + turn_1.get_coordinates()[0]) / 2, turn_0.get_coordinates()[1]});
        points.push_back(fieldPoint);
        OpenMagnetics::Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);


        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::MagneticField magneticField(OpenMagnetics::MagneticFieldStrengthModels::LAMMERANER);
        magneticField.set_fringing_effect(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double distanceCenterPoint = fieldPoint.get_point()[0] - turn_0.get_coordinates()[0];
        double expectedValue = -2 * harmonicAmplitude / (2 * std::numbers::pi * distanceCenterPoint);

        CHECK((expectedValue - field.get_data()[0].get_real()) / expectedValue < maximumError);
    }

    TEST(Test_Magnetic_Image_Method) {
        gapping = OpenMagneticsTesting::get_residual_gap();
        coreShape = "P 9/5";
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        setup();

        OpenMagnetics::WindingWindowElement windingWindow = core.get_processed_description().value().get_winding_windows()[0];

        double coreColumnHeight = core.get_columns()[0].get_height();

        auto turn_0 = coil.get_turns_description().value()[0];
        std::vector<OpenMagnetics::FieldPoint> points;
        OpenMagnetics::FieldPoint fieldPoint;
        fieldPoint.set_point(std::vector<double>{turn_0.get_coordinates()[0], coreColumnHeight / 2});
        points.push_back(fieldPoint);
        OpenMagnetics::Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);


        OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("0.475 - Grade 1");
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});
        coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        coil.delimit_and_compact();
        auto turn = coil.get_turns_description().value()[0];

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::MagneticField magneticField(OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON);
        magneticField.set_fringing_effect(false);
        magneticField.set_mirroring_dimension(1);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
 
        auto outFile = outputFilePath;
        outFile.append("Test_Magnetic_Image_Method.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::QUIVER);
        painter.set_logarithmic_scale(false);
        painter.set_mirroring_dimension(1);
        painter.set_fringing_effect(false);
        painter.set_maximum_scale_value(std::nullopt);
        painter.set_minimum_scale_value(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(field.get_data().size() == 1);
    }


    TEST(Test_Magnetic_Field_One_Turn_Rectangular) {
        numberTurns = {1};
        numberParallels = {1};
        turnsRatios = {};
        interleavingLevel = 1;
        numberStacks = 1;
        voltagePeakToPeak = 2000;
        frequency = 125000;
        coreShape = "PQ 26/25";
        coreMaterial = "3C97";
        gapping = OpenMagneticsTesting::get_grinded_gap(0.001);
        sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        setup();


        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0028);
        wire.set_nominal_value_conducting_height(0.00076);
        wire.set_nominal_value_outer_height(0.0007676);
        wire.set_nominal_value_outer_width(0.002838);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::RECTANGULAR);
        wires.push_back(wire);
        coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);

        coil.delimit_and_compact();

        auto turn = coil.get_turns_description().value()[0];

        std::vector<OpenMagnetics::FieldPoint> points;
        OpenMagnetics::FieldPoint fieldPoint;
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
        OpenMagnetics::Field inducedField;
        inducedField.set_data(points);
        inducedField.set_frequency(frequency);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::MagneticField magneticField(OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON);
        magneticField.set_fringing_effect(false);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(inputs.get_operating_point(0), magnetic, inducedField);
        auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

        double harmonicAmplitude = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_current().value().get_harmonics().value().get_amplitudes()[1];
        double expectedValue = harmonicAmplitude / (2 * std::numbers::pi * maximumWidth / 2);

        CHECK_CLOSE(field.get_data()[4].get_real(), -field.get_data()[5].get_real(), expectedValue * maximumError);
        CHECK_CLOSE(field.get_data()[4].get_imaginary(), field.get_data()[5].get_imaginary(), expectedValue * maximumError);
        // std::cout << "_______________________________________________________________" << std::endl;
        // std::cout << "field.get_data()[4].get_real(): " << field.get_data()[4].get_real() << std::endl;
        // std::cout << "field.get_data()[4].get_imaginary(): " << field.get_data()[4].get_imaginary() << std::endl;
        // std::cout << "field.get_data()[5].get_real(): " << field.get_data()[5].get_real() << std::endl;
        // std::cout << "field.get_data()[5].get_imaginary(): " << field.get_data()[5].get_imaginary() << std::endl;
        // auto outFile = outputFilePath;
        // outFile.append("Test_Magnetic_Field_One_Turn_Rectangular.svg");
        // std::filesystem::remove(outFile);
        // OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::QUIVER);
        // painter.set_number_points_x(3);
        // painter.set_number_points_y(3);
        // painter.set_logarithmic_scale(false);
        // painter.set_fringing_effect(false);
        // painter.set_mirroring_dimension(0);
        // painter.set_maximum_scale_value(std::nullopt);
        // painter.set_minimum_scale_value(std::nullopt);
        // painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        // painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        // painter.paint_coil_turns(magnetic);
        // painter.export_svg();
    }
}