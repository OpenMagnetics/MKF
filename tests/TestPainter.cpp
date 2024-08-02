#include "Painter.h"
#include "json.hpp"
#include "TestingUtils.h"
#include "Settings.h"

#include <UnitTest++.h>

#include <fstream>
#include <string>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <chrono>
#include <thread>

SUITE(FieldPainter) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

    TEST(Test_Painter_Contour_Many_Turns) {
        std::vector<int64_t> numberTurns = {23, 13};
        std::vector<int64_t> numberParallels = {2, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(true);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Contour_Many_Turns_Logarithmic_Scale) {
        std::vector<int64_t> numberTurns = {23, 13};
        std::vector<int64_t> numberParallels = {2, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_Logarithmic_Scale.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(true);
        settings->set_painter_include_fringing(true);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Contour_Many_Turns_No_Fringing) {
        std::vector<int64_t> numberTurns = {23, 13};
        std::vector<int64_t> numberParallels = {2, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
 
        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_No_Fringing.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Contour_Many_Turns_Limit_Scale) {
        std::vector<int64_t> numberTurns = {23, 13};
        std::vector<int64_t> numberParallels = {2, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_Limit_Scale.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(true);
        settings->set_painter_maximum_value_colorbar(25500);
        settings->set_painter_minimum_value_colorbar(0);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Contour_One_Turn) {
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_One_Turn.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(true);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Quiver_Many_Turns) {
        std::vector<int64_t> numberTurns = {23, 13};
        std::vector<int64_t> numberParallels = {2, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(true);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Quiver_One_Turn) {
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_One_Turn.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(true);
        settings->set_painter_include_fringing(true);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Quiver_Many_Turns_No_Fringing) {
        std::vector<int64_t> numberTurns = {23, 13};
        std::vector<int64_t> numberParallels = {2, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_No_Fringing.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_include_fringing(false);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Quiver_Many_Turns_Logarithmic_Scale) {
        std::vector<int64_t> numberTurns = {23, 13};
        std::vector<int64_t> numberParallels = {2, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_Logarithmic_Scale.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(true);
        settings->set_painter_include_fringing(true);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Quiver_Many_Turns_Limit_Scale) {
        std::vector<int64_t> numberTurns = {23, 13};
        std::vector<int64_t> numberParallels = {2, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_Limit_Scale.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(true);
        settings->set_painter_maximum_value_colorbar(2500);
        settings->set_painter_minimum_value_colorbar(0);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Quiver_One_Turn_Rectangular) {
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        // std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Rectangular 2.36x1.12 - Grade 1");
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_One_Turn_Rectangular.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Contour_One_Turn_Rectangular) {
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        // std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Rectangular 2.36x1.12 - Grade 1");
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_One_Turn_Rectangular.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Quiver_Many_Turns_Rectangular) {
        std::vector<int64_t> numberTurns = {10};
        std::vector<int64_t> numberParallels = {1};
        // std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Rectangular 4.50x0.90 - Grade 1");
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires, false);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_Rectangular.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Contour_Many_Turns_Rectangular) {
        std::vector<int64_t> numberTurns = {10};
        std::vector<int64_t> numberParallels = {1};
        // std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Rectangular 4.50x0.90 - Grade 1");
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires, false);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_Rectangular.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Quiver_One_Turn_Foil) {
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        // std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Foil 0.15");
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithTolerance);
        wire.set_outer_width(wire.get_conducting_width().value());
        wire.set_outer_height(dimensionWithTolerance);
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_One_Turn_Foil.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Contour_One_Turn_Foil) {
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        // std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Foil 0.15");
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithTolerance);
        wire.set_outer_width(wire.get_conducting_width().value());
        wire.set_outer_height(dimensionWithTolerance);
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_One_Turn_Foil.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Quiver_Many_Turns_Foil) {
        std::vector<int64_t> numberTurns = {10};
        std::vector<int64_t> numberParallels = {1};
        // std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Foil 0.15");
        OpenMagnetics::DimensionWithTolerance dimensionWithToleranceHeight;
        dimensionWithToleranceHeight.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithToleranceHeight);
        OpenMagnetics::DimensionWithTolerance dimensionWithToleranceWidth;
        dimensionWithToleranceWidth.set_nominal(0.2e-3);
        wire.set_outer_width(dimensionWithToleranceWidth);
        wire.set_outer_height(dimensionWithToleranceHeight);
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_Foil.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Contour_Many_Turns_Foil) {
        std::vector<int64_t> numberTurns = {10};
        std::vector<int64_t> numberParallels = {1};
        // std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        // auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Foil 0.15");
        OpenMagnetics::DimensionWithTolerance dimensionWithToleranceHeight;
        dimensionWithToleranceHeight.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithToleranceHeight);
        OpenMagnetics::DimensionWithTolerance dimensionWithToleranceWidth;
        dimensionWithToleranceWidth.set_nominal(0.2e-3);
        wire.set_outer_width(dimensionWithToleranceWidth);
        wire.set_outer_height(dimensionWithToleranceHeight);
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_Foil.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }
}


SUITE(ToridalFieldPainter) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

    TEST(Test_Painter_Toroid_Round_Wires) {
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({100, 5});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 1e-3;
        std::string shapeName = "T 20/10/7";

        OpenMagnetics::Processed processed;
        processed.set_label(label);
        processed.set_offset(offset);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset,
                                                                                         turnsRatios);

        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;


        settings->reset();

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_Toroid_Round_Wires.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile, true);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            // settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
            settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
            settings->set_painter_logarithmic_scale(false);
            settings->set_painter_include_fringing(true);
            settings->set_painter_number_points_x(50);
            settings->set_painter_number_points_y(50);
            settings->set_painter_maximum_value_colorbar(std::nullopt);
            settings->set_painter_minimum_value_colorbar(std::nullopt);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }

    TEST(Test_Painter_Toroid_Quiver_One_Turn_Rectangular) {
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        std::vector<OpenMagnetics::WireWrapper> wires;
        wires.push_back({OpenMagnetics::find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_One_Turn_Rectangular.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        settings->set_painter_number_points_x(100);
        settings->set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Toroid_Quiver_One_Turn_Rectangular_Inner) {
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        std::vector<OpenMagnetics::WireWrapper> wires;
        wires.push_back({OpenMagnetics::find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_One_Turn_Rectangular_Inner.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        settings->set_painter_number_points_x(100);
        settings->set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Toroid_Quiver_Four_Turns_Rectangular_Inner) {
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        std::vector<OpenMagnetics::WireWrapper> wires;
        wires.push_back({OpenMagnetics::find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_Four_Turns_Rectangular_Inner.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        settings->set_painter_number_points_x(100);
        settings->set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Toroid_Quiver_Four_Turns_Rectangular_Spread) {
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        std::vector<OpenMagnetics::WireWrapper> wires;
        wires.push_back({OpenMagnetics::find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_Four_Turns_Rectangular_Spread.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        settings->set_painter_number_points_x(100);
        settings->set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Toroid_Quiver_Two_Turn_Rectangular) {
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        std::vector<OpenMagnetics::WireWrapper> wires;
        wires.push_back({OpenMagnetics::find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_Two_Turn_Rectangular.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        settings->set_painter_number_points_x(100);
        settings->set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Toroid_Rectangular_Wires) {
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns = {11, 90};
        std::vector<int64_t> numberParallels = {1, 1};
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        uint8_t interleavingLevel = 1;
        std::string coreShape = "T 20/10/7";
        auto emptyGapping = json::array();
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::vector<OpenMagnetics::WireWrapper> wires;

        wires.push_back({OpenMagnetics::find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        wires.push_back({OpenMagnetics::find_wire_by_name("Round 0.335 - Grade 1")});


        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 1e-3;
        std::string shapeName = "T 20/10/7";

        OpenMagnetics::Processed processed;
        processed.set_label(label);
        processed.set_offset(offset);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset,
                                                                                         turnsRatios);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         sectionOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;


        settings->reset();

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_Toroid_Rectangular_Wires.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile, true);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            // settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
            settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
            settings->set_painter_logarithmic_scale(false);
            settings->set_painter_include_fringing(true);
            settings->set_painter_number_points_x(50);
            settings->set_painter_number_points_y(50);
            settings->set_painter_maximum_value_colorbar(std::nullopt);
            settings->set_painter_minimum_value_colorbar(std::nullopt);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }
}

SUITE(CoilPainterToroid) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");

    TEST(Test_Painter_T_Core) { 
        OpenMagnetics::clear_databases();
        settings->set_coil_try_rewind(false);
        std::vector<int64_t> numberTurns = {12, 12};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 80/20/50";
        std::string coreMaterial = "3C97"; 
        settings->set_coil_delimit_and_compact(false);
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel,
                                                        OpenMagnetics::WindingOrientation::OVERLAPPING,
                                                        OpenMagnetics::WindingOrientation::OVERLAPPING,
                                                        OpenMagnetics::CoilAlignment::SPREAD,
                                                        OpenMagnetics::CoilAlignment::SPREAD
                                                        );
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_T_Core.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_T_Core_Overlapping) {
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns = {2, 2};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "T 80/20/50";
        std::string coreMaterial = "3C97"; 
        settings->set_coil_delimit_and_compact(false);
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, OpenMagnetics::WindingOrientation::OVERLAPPING);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_T_Core_Overlapping.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }


    TEST(Test_Painter_T_Core_Contiguous) {
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns = {72, 72};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "T 80/20/50";
        std::string coreMaterial = "3C97";
        // settings->set_coil_delimit_and_compact(false);
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_T_Core_Contiguous_turns.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            // painter.paint_coil_layers(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_T_Core_Contiguous_layers.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_layers(magnetic);
            // painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
            settings->reset();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_T_Core_Contiguous_sections.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_sections(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
            settings->reset();
        }
    }

    TEST(Test_Painter_T_Core_Contiguous_Sections_With_Margin) {
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns = {2, 2};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "T 80/20/50";
        std::string coreMaterial = "3C97"; 
        settings->set_coil_delimit_and_compact(false);
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, OpenMagnetics::WindingOrientation::CONTIGUOUS);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        double margin = 0.005;
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 1.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        coil.add_margin_to_section_by_index(3, std::vector<double>{margin * 3.5, margin * 3.5});

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_T_Core_Contiguous_Sections_With_Margin.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_sections(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }
}

SUITE(CoilPainter) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");

    TEST(Test_Painter_Pq_Core_Distributed_Gap) {
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Distributed_Gap.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Distributed_Gap_Many) {
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.001, 9);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Distributed_Gap_Many.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap) {
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_U_Core_Distributed_Gap) {
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "U 10/8/3";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.001, 3);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_U_Core_Distributed_Gap.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_U_Core_Grinded_Gap) {
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "U 10/8/3";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_U_Core_Grinded_Gap.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Bobbin) {
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Section) {
        std::vector<int64_t> numberTurns = {42};
        std::vector<int64_t> numberParallels = {3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Sections.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Bobbin_And_Section) {
        std::vector<int64_t> numberTurns = {54};
        std::vector<int64_t> numberParallels = {2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);


        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_And_Section.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Bobbin_And_Two_Sections) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_And_Sections.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Epx_Core_Grinded_Gap) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "EPX 9/9";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Epx_Core_Grinded_Gap.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Epx_Core_Spacer_Gap) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "EPX 9/9";
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto basicSpacerGap = json();
        basicSpacerGap["type"] = "additive";
        basicSpacerGap["length"] = 0.0003;
        gapping.push_back(basicSpacerGap);
        gapping.push_back(basicSpacerGap);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Epx_Core_Spacer_Gap.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_P_Core_Grinded_Gap) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "P 3.3/2.6";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_P_Core_Grinded_Gap.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_U80_Core_Grinded_Gap) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "U 80/65/32";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_U80_Core_Grinded_Gap.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Ep_Core_Grinded_Gap) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "EP 10";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Ep_Core_Grinded_Gap.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_All_Cores) {

        for (std::string& shapeName : OpenMagnetics::get_shape_names()){
            if (shapeName.contains("PQI") || shapeName.contains("R ") || shapeName.contains("T ") || shapeName.contains("UI ")) {
                continue;
            }
            std::vector<int64_t> numberTurns = {42, 42};
            std::vector<int64_t> numberParallels = {3, 3};
            uint8_t interleavingLevel = 2;
            int64_t numberStacks = 1;
            std::string coreShape = shapeName;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
            auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

            std::replace( shapeName.begin(), shapeName.end(), '.', '_');
            std::replace( shapeName.begin(), shapeName.end(), '/', '_');
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_Core_" + shapeName + ".svg");
        std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_sections(magnetic);

            painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Layers_No_Interleaving) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Layers_No_Interleaving.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_No_Interleaving) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_No_Interleaving.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Top_Alignment) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Top_Alignment.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Bottom_Alignment) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Bottom_Alignment.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Spread_Alignment) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Spread_Alignment.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Vertical_Sections) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Vertical_Sections_Vectical_Layers) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {3, 3};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Vectical_Layers.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers_Spread_Turns) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Spread_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers_Inner_Turns) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Inner_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers_Outer_Turns) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Outer_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers_Centered_Turns) {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Centered_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Foil_Centered) {
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_outer_height(0.014);
        wire.set_nominal_value_outer_width(0.0002);
        wire.set_type(OpenMagnetics::WireType::FOIL);
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Foil_Centered.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Foil_Top) {
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_outer_height(0.014);
        wire.set_nominal_value_outer_width(0.0002);
        wire.set_type(OpenMagnetics::WireType::FOIL);
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Foil_Top.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Litz) {
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto wire = OpenMagnetics::find_wire_by_name("Litz TXXL07/28FXXX-2(MWXX)");
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Litz.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Foil_With_Insulation_Centered) {
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_outer_height(0.014);
        wire.set_nominal_value_outer_width(0.0002);
        wire.set_nominal_value_conducting_height(0.0139);
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_type(OpenMagnetics::WireType::FOIL);
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Foil_With_Insulation_Centered.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Centered) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Centered) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }


    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Horizontal_Centered) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Top) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Top.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Inner) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Inner.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Outer) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Outer.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Bottom) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Bottom.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Spread) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Spread_Two_Sections) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Spread_Two_Sections.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Spread_One_Section) {
        std::vector<int64_t> numberTurns = {23};
        std::vector<int64_t> numberParallels = {2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Spread_One_Section.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Spread) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_Two_Sections) {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_Two_Sections.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_One_Section) {
        std::vector<int64_t> numberTurns = {23};
        std::vector<int64_t> numberParallels = {2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_One_Section.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Bobbin_Vertical_Sections_And_Insulation) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Vertical_Sections_And_Insulation.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Bobbin_Horizontal_Sections_And_Insulation) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Horizontal_Sections_And_Insulation.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Bobbin_Layers_And_Insulation) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Layers_And_Insulation.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Painter_Pq_Core_Bobbin_Turns_And_Insulation) {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Turns_And_Insulation.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Test_Turns_Not_Fitting) {

        settings->set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {6, 6};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_outer_diameter(0.0015);
        wire.set_type(OpenMagnetics::WireType::ROUND);
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire, wire});

        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Turns_Not_Fitting.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        try
        {
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(false);
        }
        catch (const std::exception &e)
        {
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(true);
            std::cerr << e.what() << std::endl;
        }
        settings->reset();
    }

    TEST(Test_Painter_Planar) {
        std::vector<int64_t> numberTurns = {12, 6};
        std::vector<int64_t> numberParallels = {1, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "ER 18/3/10";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;


        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_height(0.000034);
        wire.set_nominal_value_conducting_width(0.001);
        wire.set_nominal_value_outer_height(0.000134);
        wire.set_nominal_value_outer_width(0.0011);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::RECTANGULAR);
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire, wire});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires, false);
        // auto bobbin = coil.resolve_bobbin();
        // auto processedDescription = bobbin.get_processed_description().value();
        // processedDescription.get_mutable_winding_windows()[0].set_height(processedDescription.get_mutable_winding_windows()[0].get_height().value() / 2);
        // bobbin.set_processed_description(processedDescription);
        // coil.set_bobbin(bobbin);


        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.wind();
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Planar.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, true);
        settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
        settings->set_painter_logarithmic_scale(false);
        settings->set_painter_include_fringing(false);
        settings->set_painter_maximum_value_colorbar(std::nullopt);
        settings->set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));
        settings->reset();
    }

    TEST(Course0) {
        for (int64_t index = 1; index < 20; index++) {

            std::vector<int64_t> numberTurns = {index, index * 10};
            std::vector<int64_t> numberParallels = {1, 1};
            uint8_t interleavingLevel = 1;
            int64_t numberStacks = 1;
            std::string coreShape = "PQ 26/25";
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
            auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

            auto outFile = outputFilePath;
            outFile.append("Course" + std::to_string(index) + ".svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);

            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
            settings->reset();
        }
    }
}

SUITE(WirePainter) {
    auto settings = OpenMagnetics::Settings::GetInstance();

    TEST(Test_Wire_Painter_Round_Enamelled_Grade_1) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Round 0.335 - Grade 1");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_Enamelled_Grade_1.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }
    TEST(Test_Wire_Painter_Round_Enamelled_Grade_2) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Round 0.335 - Grade 2");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_Enamelled_Grade_2.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Round_Enamelled_Grade_9) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Round 0.071 - FIW 9");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_Enamelled_Grade_9.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Round_TIW_PFA) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Round T17A01PXXX-1.5");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_TIW_PFA.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Round_DIW_ETFE) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Round D32A01TXX-3");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_DIW_ETFE.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Round_SIW_PEF) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Round S24A01FX-2");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_SIW_PEF.svg");
            std::filesystem::remove(outFile)
;            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Round_SIW_TCA) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Round TCA1 31 AWG");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_SIW_TCA.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Rectangular_Enamelled_Grade_1) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Rectangular 2x0.80 - Grade 1");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Rectangular_Enamelled_Grade_1.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Rectangular_Enamelled_Grade_2) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Rectangular 2x0.80 - Grade 2");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Rectangular_Enamelled_Grade_2.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Foil_Bare) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Foil 0.7");
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithTolerance);
        wire.set_outer_width(wire.get_conducting_width().value());
        wire.set_outer_height(dimensionWithTolerance);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Foil_Bare.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Litz_Simple_TIW) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Litz TXXL07/28FXXX-2(MWXX)");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_Insulation.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Litz_Simple_Bare) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Litz 75x0.3 - Grade 1 - Unserved");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_Bare.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Litz_Simple_Single_Served) {
        OpenMagnetics::clear_databases();

        auto wire = OpenMagnetics::find_wire_by_name("Litz 75x0.3 - Grade 1 - Single Served");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_Single_Served.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Litz_Simple_Double_Served) {
        OpenMagnetics::clear_databases(); 
        settings->set_painter_simple_litz(false);
        settings->set_painter_advanced_litz(true);

        auto wire = OpenMagnetics::find_wire_by_name("Litz 270x0.02 - Grade 1 - Single Served");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_Double_Served.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Litz_Simple_Insulated) {
        OpenMagnetics::clear_databases(); 
        settings->set_painter_simple_litz(false);
        settings->set_painter_advanced_litz(true);

        auto wire = OpenMagnetics::find_wire_by_name("Litz DXXL550/44TXX-3(MWXX)");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_sInsulated.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Litz_Few_Strands) {
        OpenMagnetics::clear_databases(); 
        settings->set_painter_simple_litz(false);
        settings->set_painter_advanced_litz(true);

        auto wire = OpenMagnetics::find_wire_by_name("Litz 4x0.071 - Grade 2 - Double Served");

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Few_Strands.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    // TEST(Test_Wire_Painter_All_Litz) {
    //     OpenMagnetics::clear_databases(); 
    //     settings->set_painter_simple_litz(false);
    //     settings->set_painter_advanced_litz(true);

    //     auto wires = OpenMagnetics::get_wires(OpenMagnetics::WireType::LITZ);

    //     for (auto wire : wires) {
    //         std::string wireName = wire.get_name().value();
    //         std::replace( wireName.begin(), wireName.end(), '/', '_');
    //         std::cout << wireName << std::endl;
    //         auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
    //         auto outFile = outputFilePath;
    //         outFile.append("Test_Wire_Painter_Litz_" + wireName + ".svg");
    //         std::filesystem::remove(outFile);
    //         OpenMagnetics::Painter painter(outFile);
    //         painter.paint_wire(wire);
    //         painter.export_svg();
    //     }
    //     settings->reset();
    // }
}

SUITE(WirePainterCurrentDensity) {
    auto settings = OpenMagnetics::Settings::GetInstance();

    TEST(Test_Wire_Painter_Current_Density_Round_Enamelled_Grade_1) {
        OpenMagnetics::clear_databases();

        double currentPeakToPeak = 1000;
        double currentRms = currentPeakToPeak / sqrt(3);
        double frequency = 1250000;
        auto wire = OpenMagnetics::find_wire_by_name("Round 3.55 - Grade 1");
        double conductingArea = wire.calculate_conducting_area();
        double currentDensity = currentRms / conductingArea;
        std::cout << "currentDensity: " << currentDensity << std::endl;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Round_Enamelled_Grade_1.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.paint_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings->reset();
    }
}