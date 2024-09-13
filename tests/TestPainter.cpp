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


SUITE(ToroidalFieldPainter) {
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


    TEST(Test_Coil_Painter_Web_0) {
        OpenMagnetics::clear_databases();

        OpenMagnetics::MagneticWrapper magnetic(json::parse(R"({"coil": {"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.005974999999999999, "columnShape": "rectangular", "columnThickness": 0.0013999999999999993, "columnWidth": 0.005999999999999999, "coordinates": [0, 0, 0], "pins": null, "wallThickness": 0.0011204000000000006, "windingWindows": [{"angle": null, "area": 0.00011625151999999999, "coordinates": [0.008799999999999999, 0, 0], "height": 0.0207592, "radialHeight": null, "sectionsOrientation": "overlapping", "shape": "rectangular", "width": 0.0056}]}}, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "primary", "numberParallels": 1, "numberTurns": 20, "wire": {"coating": {"breakdownVoltage": 18000, "grade": null, "material": "PFA", "numberLayers": 3, "temperatureRating": 180, "thickness": null, "thicknessLayers": 7.62e-05, "type": "insulated"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 6.532502100168472e-07}, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.000912}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Rubadue", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round T19A01PXXX-3", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.001369}, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": "19 AWG", "strand": null, "type": "round"}}], "layersDescription": [{"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 8.673617379884035e-19], "dimensions": [0.0045850000000000005, 0.020700000000000003], "fillingFactor": 0.9971482523411308, "insulationMaterial": null, "name": "primary section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [1], "winding": "primary"}], "section": "primary section 0", "turnsAlignment": "centered", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}], "sectionsDescription": [{"coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0], "dimensions": [0.0045850000000000005, 0.020700000000000003], "fillingFactor": 0.9624513055163947, "layersAlignment": null, "layersOrientation": "overlapping", "margin": [0, 0], "name": "primary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [1], "winding": "primary"}], "type": "conduction", "windingStyle": "windByConsecutiveTurns"}], "turnsDescription": [{"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.009832500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 0", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.008797500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 1", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.007762500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 2", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.006727500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 3", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.005692500000000002], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 4", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.004657500000000002], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 5", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.0036225000000000016], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 6", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.0025875000000000013], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 7", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.0015525000000000011], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 8", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.000517500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 9", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0005174999999999991], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 10", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0015524999999999992], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 11", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0025874999999999995], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 12", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0036225], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 13", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0046575], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 14", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0056925000000000005], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 15", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.006727500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 16", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.007762500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 17", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.008797500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 18", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.009832500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 19", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}]}, "core": {"distributorsInfo": [{"cost": 1.87, "country": "USA", "distributedArea": "International", "email": null, "link": "https://www.digikey.com/en/products/detail/epcos-tdk-electronics/B66229G0500X127/3914560", "name": "Digi-Key", "phone": null, "quantity": 660, "reference": "495-B66229G0500X127-ND", "updatedAt": "16/04/2024"}, {"cost": 1.86, "country": "USA", "distributedArea": "International", "email": null, "link": "https://www.mouser.es/ProductDetail/EPCOS-TDK/B66229G0500X127?qs=%2FsLciWRBLmC6djd1SV6djQ%3D%3D", "name": "Mouser", "phone": null, "quantity": 355, "reference": "871-B66229G0500X127", "updatedAt": "16/04/2024"}], "functionalDescription": {"coating": null, "gapping": [{"area": 8.5e-05, "coordinates": [0, 0.00025, 0], "distanceClosestNormalSurface": 0.011, "distanceClosestParallelSurface": 0.006999999999999999, "length": 0.0005, "sectionDimensions": [0.0092, 0.00915], "shape": "rectangular", "type": "subtractive"}, {"area": 4.1e-05, "coordinates": [0.013825, 0, 0], "distanceClosestNormalSurface": 0.011498, "distanceClosestParallelSurface": 0.006999999999999999, "length": 5e-06, "sectionDimensions": [0.004451, 0.00915], "shape": "rectangular", "type": "residual"}, {"area": 4.1e-05, "coordinates": [-0.013825, 0, 0], "distanceClosestNormalSurface": 0.011498, "distanceClosestParallelSurface": 0.006999999999999999, "length": 5e-06, "sectionDimensions": [0.004451, 0.00915], "shape": "rectangular", "type": "residual"}], "material": "N27", "numberStacks": 1, "shape": {"aliases": ["E 32/9", "EF 32"], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0329, "minimum": 0.0313, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0164, "minimum": 0.0158, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0095, "minimum": 0.0088, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0118, "minimum": 0.0112, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0237, "minimum": 0.0227, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0095, "minimum": 0.0089, "nominal": null}}, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 32/16/9", "type": "standard"}, "type": "two-piece set"}, "geometricalDescription": null, "manufacturerInfo": {"cost": null, "datasheetUrl": "https://product.tdk.com/system/files/dam/doc/product/ferrite/ferrite/ferrite-acc/data_sheet/80/db/fer/e_32_16_9.pdf", "family": null, "name": "TDK", "orderCode": null, "reference": "B66229G0500X127 (N27", "status": "production"}, "name": "E 32/16/9 - N27 - Gapped 0.500 mm", "processedDescription": {"columns": [{"area": 8.5e-05, "coordinates": [0, 0, 0], "depth": 0.00915, "height": 0.023, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "central", "width": 0.0092}, {"area": 4.1e-05, "coordinates": [0.013825, 0, 0], "depth": 0.00915, "height": 0.023, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.004451}, {"area": 4.1e-05, "coordinates": [-0.013825, 0, 0], "depth": 0.00915, "height": 0.023, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.004451}], "depth": 0.00915, "effectiveParameters": {"effectiveArea": 8.316165619067432e-05, "effectiveLength": 0.07431657444154591, "effectiveVolume": 6.180289412976495e-06, "minimumArea": 8.143500000000005e-05}, "height": 0.032200000000000006, "width": 0.032100000000000004, "windingWindows": [{"angle": null, "area": 0.00016099999999999998, "coordinates": [0.0046, 0], "height": 0.023, "radialHeight": null, "sectionsOrientation": null, "shape": null, "width": 0.006999999999999999}]}}, "distributorsInfo": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "", "recommendations": null, "reference": "E 32/16/9 - N27 - 0.5 mm, Turns: 20, Order: 0, Non-Interleaved, Margin Taped 00", "status": null}, "rotation": null})"));
        // settings->set_painter_simple_litz(false);
        // settings->set_painter_advanced_litz(false);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Painter_Web_0.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            CHECK(std::filesystem::exists(outFile));
        }
        settings->reset();
    }


    TEST(Test_Coil_Painter_Web_1) {
        OpenMagnetics::clear_databases();

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_1.svg");
        OpenMagnetics::MagneticWrapper magnetic(json::parse(R"({"coil": {"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.005974999999999999, "columnShape": "rectangular", "columnThickness": 0.0013999999999999993, "columnWidth": 0.005999999999999999, "coordinates": [0, 0, 0], "pins": null, "wallThickness": 0.0011204000000000006, "windingWindows": [{"angle": null, "area": 0.00011625151999999999, "coordinates": [0.008799999999999999, 0, 0], "height": 0.0207592, "radialHeight": null, "sectionsOrientation": "overlapping", "shape": "rectangular", "width": 0.0056}]}}, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "primary", "numberParallels": 1, "numberTurns": 18, "wire": {"coating": {"breakdownVoltage": 18000, "grade": null, "material": "PFA", "numberLayers": 3, "temperatureRating": 180, "thickness": null, "thicknessLayers": 7.62e-05, "type": "insulated"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 6.532502100168472e-07}, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.000912}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Rubadue", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round T19A01PXXX-3", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.001369}, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": "19 AWG", "strand": null, "type": "round"}}], "layersDescription": [{"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 8.673617379884035e-19], "dimensions": [0.0045850000000000005, 0.020700000000000003], "fillingFactor": 0.9971482523411308, "insulationMaterial": null, "name": "primary section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [1], "winding": "primary"}], "section": "primary section 0", "turnsAlignment": "centered", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}], "sectionsDescription": [{"coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0], "dimensions": [0.0045850000000000005, 0.020700000000000003], "fillingFactor": 0.9624513055163947, "layersAlignment": null, "layersOrientation": "overlapping", "margin": [0, 0], "name": "primary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [1], "winding": "primary"}], "type": "conduction", "windingStyle": "windByConsecutiveTurns"}], "turnsDescription": [{"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.009832500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 0", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.008797500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 1", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.007762500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 2", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.006727500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 3", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.005692500000000002], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 4", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.004657500000000002], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 5", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.0036225000000000016], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 6", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.0025875000000000013], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 7", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.0015525000000000011], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 8", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.000517500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 9", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0005174999999999991], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 10", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0015524999999999992], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 11", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0025874999999999995], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 12", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0036225], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 13", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0046575], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 14", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0056925000000000005], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 15", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.006727500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 16", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.007762500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 17", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.008797500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 18", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.009832500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 19", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}]}, "core": {"distributorsInfo": [], "functionalDescription": {"coating": null, "gapping": [{"length": 0.00019999999999999998, "type": "subtractive"}, {"length": 5e-06, "type": "residual"}, {"length": 5e-06, "type": "residual"}], "material": "3F36", "numberStacks": 1, "shape": {"aliases": ["E 34.6/9"], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0353, "minimum": 0.0339, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01464, "minimum": 0.0139, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00972, "minimum": 0.0089, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01005, "minimum": 0.00951, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0262, "minimum": 0.025, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0097, "minimum": 0.0091, "nominal": null}}, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 34/14/9", "type": "standard"}, "type": "two-piece set"}, "geometricalDescription": [{"coordinates": [0, 0, 0], "dimensions": null, "insulationMaterial": null, "machining": [{"coordinates": [0, 0.0001, 0], "length": 0.0002}], "material": "3F36", "rotation": [3.141592653589793, 3.141592653589793, 0], "shape": {"aliases": ["E 34.6/9"], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0353, "minimum": 0.0339, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01464, "minimum": 0.0139, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00972, "minimum": 0.0089, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01005, "minimum": 0.00951, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0262, "minimum": 0.025, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0097, "minimum": 0.0091, "nominal": null}}, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 34/14/9", "type": "standard"}, "type": "half set"}, {"coordinates": [0, 0, 0], "dimensions": null, "insulationMaterial": null, "machining": null, "material": "3F36", "rotation": [0, 0, 0], "shape": {"aliases": ["E 34.6/9"], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0353, "minimum": 0.0339, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01464, "minimum": 0.0139, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00972, "minimum": 0.0089, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01005, "minimum": 0.00951, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0262, "minimum": 0.025, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0097, "minimum": 0.0091, "nominal": null}}, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 34/14/9", "type": "standard"}, "type": "half set"}], "manufacturerInfo": {"cost": null, "datasheetUrl": "https://ferroxcube.com/upload/media/product/file/Pr_ds/E34_14_9.pdf", "family": null, "name": "Ferroxcube", "orderCode": null, "reference": "E34/14/9-3F36-G200", "status": "production"}, "name": "E 34/14/9 - 3F36 - Gapped 0.200 mm", "processedDescription": {"columns": [{"area": 8.8e-05, "coordinates": [0, 0, 0], "depth": 0.00931, "height": 0.01956, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "central", "width": 0.0094}, {"area": 4.2e-05, "coordinates": [0.01505, 0, 0], "depth": 0.00931, "height": 0.01956, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.0045}, {"area": 4.2e-05, "coordinates": [-0.01505, 0, 0], "depth": 0.00931, "height": 0.01956, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.0045}], "depth": 0.009309999999999999, "effectiveParameters": {"effectiveArea": 8.490169410503223e-05, "effectiveLength": 0.06957187369722902, "effectiveVolume": 5.9067699389560765e-06, "minimumArea": 8.360379999999998e-05}, "height": 0.02854, "width": 0.0346, "windingWindows": [{"angle": null, "area": 0.000158436, "coordinates": [0.0047, 0], "height": 0.01956, "radialHeight": null, "sectionsOrientation": null, "shape": null, "width": 0.0081}]}}, "distributorsInfo": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "", "recommendations": null, "reference": "E 32/16/9 - N27 - 0.5 mm, Turns: 20, Order: 0, Non-Interleaved, Margin Taped 00", "status": null}, "rotation": null})"));

        OpenMagnetics::Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings->reset();
    }

    TEST(Test_Coil_Painter_Web_2) {
        OpenMagnetics::clear_databases();
        // settings->set_coil_delimit_and_compact(false);

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_2.svg");
        OpenMagnetics::MagneticWrapper magnetic(json::parse(R"({"coil": {"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.003912540921738127, "columnShape": "round", "columnThickness": 0.0010750409217381266, "columnWidth": 0.003912540921738127, "coordinates": [0, 0, 0], "pins": null, "wallThickness": 0.0007636652757989004, "windingWindows": [{"angle": null, "area": 2.0400047558919626e-05, "coordinates": [0.004956270460869064, 0, 0], "height": 0.009772669448402199, "radialHeight": null, "sectionsAlignment": "spread", "sectionsOrientation": "contiguous", "shape": "rectangular", "width": 0.002087459078261874}]}}, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 32, "wire": {"coating": {"breakdownVoltage": 2300, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 1.5904312808798332e-07}, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000455, "minimum": 0.00044500000000000003, "nominal": 0.00045000000000000004}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.45 - Grade 1", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000491, "minimum": 0.000472, "nominal": null}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": "0.45 mm", "strand": null, "type": "round"}}, {"connections": null, "isolationSide": "secondary", "name": "Secondary", "numberParallels": 2, "numberTurns": 1, "wire": {"coating": {"breakdownVoltage": 1350, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 2.0508395382450515e-07}, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000513, "minimum": 0.000505, "nominal": 0.0005110000000000001}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 24.0 - Single Build", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0005510000000000001, "minimum": 0.000531, "nominal": 0.0005409999999990001}, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": "24 AWG", "strand": null, "type": "round"}}], "layersDescription": [{"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.0029603347242010986], "dimensions": [0.0004815, 0.003852000000000001], "fillingFactor": 0.7903427177982021, "insulationMaterial": null, "name": "Primary section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.25], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "inner or top", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.0029603347242010986], "dimensions": [0.0004815, 0.003852000000000001], "fillingFactor": 0.7903427177982021, "insulationMaterial": null, "name": "Primary section 0 layer 1", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.25], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "inner or top", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.0029603347242010986], "dimensions": [0.0004815, 0.003852000000000001], "fillingFactor": 0.7903427177982021, "insulationMaterial": null, "name": "Primary section 0 layer 2", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.25], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "inner or top", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.0029603347242010986], "dimensions": [0.0004815, 0.003852000000000001], "fillingFactor": 0.7903427177982021, "insulationMaterial": null, "name": "Primary section 0 layer 3", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.25], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "inner or top", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.00495627, -0.0009517214828007338], "dimensions": [0.002087459078261874, 2.5e-05], "fillingFactor": 1, "insulationMaterial": null, "name": "Insulation between Primary and Primary section 1 layer 0", "orientation": "contiguous", "partialWindings": [], "section": "Insulation between Primary and Primary section 1", "turnsAlignment": "spread", "type": "insulation", "windingStyle": null}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.0057295000000005016, -0.0034787782413993667], "dimensions": [0.0005409999999990001, 0.0010819999999979997], "fillingFactor": 0.22200177067914698, "insulationMaterial": null, "name": "Secondary section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [1, 1], "winding": "Secondary"}], "section": "Secondary section 0", "turnsAlignment": "outer or bottom", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.00495627, -0.006005834448400199], "dimensions": [0.002087459078261874, 2.5e-05], "fillingFactor": 1, "insulationMaterial": null, "name": "Insulation between Secondary and Secondary section 3 layer 0", "orientation": "contiguous", "partialWindings": [], "section": "Insulation between Secondary and Secondary section 3", "turnsAlignment": "spread", "type": "insulation", "windingStyle": null}], "sectionsDescription": [{"coordinateSystem": "cartesian", "coordinates": [0.004875540921738127, 0.0029603347242010986], "dimensions": [0.0019260000000000006, 0.0038520000000000013], "fillingFactor": 0.7292119353768599, "layersAlignment": null, "layersOrientation": "overlapping", "margin": [0, 0], "name": "Primary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [1], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"coordinateSystem": "cartesian", "coordinates": [0.004956270460869064, -0.0009517217585996346], "dimensions": [0.002087459078261874, 2.5e-05], "fillingFactor": 1, "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Primary section 1", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.0057295000000005016, -0.0034787782413993667], "dimensions": [0.0005409999999990002, 0.0010819999999979995], "fillingFactor": 0.05753547898874283, "layersAlignment": null, "layersOrientation": "overlapping", "margin": [0, 0], "name": "Secondary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [1, 1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}, {"coordinateSystem": "cartesian", "coordinates": [0.004956270460869064, -0.0060058347241991], "dimensions": [0.002087459078261874, 2.5e-05], "fillingFactor": 1, "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Secondary and Secondary section 3", "partialWindings": [], "type": "insulation", "windingStyle": null}], "turnsDescription": [{"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.004645584724201099], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 0", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.004164084724201099], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 1", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.0036825847242010984], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 2", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.0032010847242010983], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 3", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.002719584724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 4", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.002238084724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 5", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.001756584724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 6", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.001275084724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 7", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.004645584724201099], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 8", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.004164084724201099], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 9", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.0036825847242010984], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 10", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.0032010847242010983], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 11", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.002719584724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 12", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.002238084724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 13", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.001756584724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 14", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.001275084724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 15", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.004645584724201099], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 16", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.004164084724201099], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 17", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.0036825847242010984], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 18", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.0032010847242010983], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 19", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.002719584724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 20", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.002238084724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 21", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.001756584724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 22", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.001275084724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 23", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.004645584724201099], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 24", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.004164084724201099], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 25", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.0036825847242010984], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 26", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.0032010847242010983], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 27", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.002719584724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 28", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.002238084724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 29", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.001756584724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 30", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.001275084724201098], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 31", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0057295000000005016, -0.0032082782413998673], "dimensions": [0.0005409999999990001, 0.0005409999999990001], "layer": "Secondary section 0 layer 0", "length": 0.03599951021748859, "name": "Secondary parallel 0 turn 0", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Secondary section 0", "winding": "Secondary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0057295000000005016, -0.003749278241398867], "dimensions": [0.0005409999999990001, 0.0005409999999990001], "layer": "Secondary section 0 layer 0", "length": 0.03599951021748859, "name": "Secondary parallel 1 turn 0", "orientation": "clockwise", "parallel": 1, "rotation": 0, "section": "Secondary section 0", "winding": "Secondary"}]}, "core": {"distributorsInfo": [], "functionalDescription": {"coating": null, "gapping": [{"area": 2.6e-05, "coordinates": [0, 0, 0], "distanceClosestNormalSurface": 0.005443, "distanceClosestParallelSurface": 0.0031625000000000004, "length": 0.000414, "sectionDimensions": [0.005675, 0.005675], "shape": "round", "type": "subtractive"}, {"area": 0.000163, "coordinates": [0, 0, -0.007038], "distanceClosestNormalSurface": 0.005648, "distanceClosestParallelSurface": 0.0031625000000000004, "length": 5e-06, "sectionDimensions": [0.078555, 0.002075], "shape": "irregular", "type": "residual"}], "material": "3C91", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0184, "minimum": 0.0176, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0085, "minimum": 0.0083, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01125, "minimum": 0.01075, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0058, "minimum": 0.0055, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0124, "minimum": 0.0116, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00585, "minimum": 0.0055, "nominal": null}, "K": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.002925, "minimum": null, "nominal": null}}, "family": "ep", "familySubtype": null, "magneticCircuit": "open", "name": "EP 17", "type": "standard"}, "type": "two-piece set"}, "geometricalDescription": [{"coordinates": [0, 0, 0], "dimensions": null, "insulationMaterial": null, "machining": [{"coordinates": [0, 0.0001035, 0], "length": 0.000207}], "material": "3C91", "rotation": [3.141592653589793, 3.141592653589793, 0], "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0184, "minimum": 0.0176, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0085, "minimum": 0.0083, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01125, "minimum": 0.01075, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0058, "minimum": 0.0055, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0124, "minimum": 0.0116, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00585, "minimum": 0.0055, "nominal": null}, "K": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.002925, "minimum": null, "nominal": null}}, "family": "ep", "familySubtype": null, "magneticCircuit": "open", "name": "EP 17", "type": "standard"}, "type": "half set"}, {"coordinates": [0, 0, 0], "dimensions": null, "insulationMaterial": null, "machining": [{"coordinates": [0, -0.0001035, 0], "length": 0.000207}], "material": "3C91", "rotation": [0, 0, 0], "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0184, "minimum": 0.0176, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0085, "minimum": 0.0083, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01125, "minimum": 0.01075, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0058, "minimum": 0.0055, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0124, "minimum": 0.0116, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00585, "minimum": 0.0055, "nominal": null}, "K": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.002925, "minimum": null, "nominal": null}}, "family": "ep", "familySubtype": null, "magneticCircuit": "open", "name": "EP 17", "type": "standard"}, "type": "half set"}], "manufacturerInfo": {"cost": null, "datasheetUrl": "https://ferroxcube.com/upload/media/product/file/Pr_ds/EP17.pdf", "family": null, "name": "Ferroxcube", "orderCode": null, "reference": "EP17-3C91-A100", "status": "production"}, "name": "EP 17 - 3C91 - Gapped 0.414 mm", "processedDescription": {"columns": [{"area": 2.6e-05, "coordinates": [0, 0, 0], "depth": 0.005675, "height": 0.0113, "minimumDepth": null, "minimumWidth": null, "shape": "round", "type": "central", "width": 0.005675}, {"area": 0.000163, "coordinates": [0, 0, -0.007038], "depth": 0.002075, "height": 0.0113, "minimumDepth": null, "minimumWidth": 0.003001, "shape": "irregular", "type": "lateral", "width": 0.078555}], "depth": 0.011, "effectiveParameters": {"effectiveArea": 3.4461049794381896e-05, "effectiveLength": 0.02849715978271948, "effectiveVolume": 9.820420422707531e-07, "minimumArea": 2.5790801226066944e-05}, "height": 0.016800000000000002, "width": 0.018000000000000002, "windingWindows": [{"angle": null, "area": 3.573625e-05, "coordinates": [0.0028374999999999997, 0], "height": 0.0113, "radialHeight": null, "sectionsAlignment": null, "sectionsOrientation": null, "shape": null, "width": 0.0031625000000000004}]}}})"));

        magnetic.get_mutable_coil().set_layers_orientation(OpenMagnetics::WindingOrientation::OVERLAPPING, "Primary section 0");
        magnetic.get_mutable_coil().set_layers_orientation(OpenMagnetics::WindingOrientation::CONTIGUOUS, "Secondary section 0");
        magnetic.get_mutable_coil().set_turns_alignment(OpenMagnetics::CoilAlignment::INNER_OR_TOP, "Primary section 0");
        magnetic.get_mutable_coil().set_turns_alignment(OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM, "Secondary section 0");

        magnetic.get_mutable_coil().clear();
        magnetic.get_mutable_coil().wind();
        // magnetic.get_mutable_coil().delimit_and_compact();

        OpenMagnetics::Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings->reset();
    }

    TEST(Test_Coil_Painter_Web_3) {
        OpenMagnetics::clear_databases();
        // settings->set_coil_delimit_and_compact(false);

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_3.svg");
        OpenMagnetics::MagneticWrapper magnetic(json::parse(R"({"coil": {"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.003912540921738127, "columnShape": "round", "columnThickness": 0.0010750409217381266, "columnWidth": 0.003912540921738127, "coordinates": [0, 0, 0], "pins": null, "wallThickness": 0.0007636652757989004, "windingWindows": [{"angle": null, "area": 2.0400047558919626e-05, "coordinates": [0.004956270460869064, 0, 0], "height": 0.009772669448402199, "radialHeight": null, "sectionsAlignment": "outer or bottom", "sectionsOrientation": "contiguous", "shape": "rectangular", "width": 0.002087459078261874}]}}, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 32, "wire": {"coating": {"breakdownVoltage": 2300, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 1.5904312808798332e-07}, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000455, "minimum": 0.00044500000000000003, "nominal": 0.00045000000000000004}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.45 - Grade 1", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000491, "minimum": 0.000472, "nominal": null}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": "0.45 mm", "strand": null, "type": "round"}}, {"connections": null, "isolationSide": "secondary", "name": "Secondary", "numberParallels": 2, "numberTurns": 1, "wire": {"coating": {"breakdownVoltage": 1350, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 2.0508395382450515e-07}, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000513, "minimum": 0.000505, "nominal": 0.0005110000000000001}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 24.0 - Single Build", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0005510000000000001, "minimum": 0.000531, "nominal": 0.0005409999999990001}, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": "24 AWG", "strand": null, "type": "round"}}], "layersDescription": [{"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, -0.0014062832242030996], "dimensions": [0.0004815, 0.004746103000000001], "fillingFactor": 0.7903427177982021, "insulationMaterial": null, "name": "Primary section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.25], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, -0.0014062832242030996], "dimensions": [0.0004815, 0.004746103000000001], "fillingFactor": 0.7903427177982021, "insulationMaterial": null, "name": "Primary section 0 layer 1", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.25], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, -0.0014062832242030996], "dimensions": [0.0004815, 0.004746103000000001], "fillingFactor": 0.7903427177982021, "insulationMaterial": null, "name": "Primary section 0 layer 2", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.25], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, -0.0014062832242030996], "dimensions": [0.0004815, 0.004746103000000001], "fillingFactor": 0.7903427177982021, "insulationMaterial": null, "name": "Primary section 0 layer 3", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.25], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.00495627, -0.003791834448404199], "dimensions": [0.002087459078261874, 2.5e-05], "fillingFactor": 1, "insulationMaterial": null, "name": "Insulation between Primary and Primary section 1 layer 0", "orientation": "contiguous", "partialWindings": [], "section": "Insulation between Primary and Primary section 1", "turnsAlignment": "spread", "type": "insulation", "windingStyle": null}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.0057295000000005016, -0.0043453347242020995], "dimensions": [0.0005409999999990001, 0.0010819999999979997], "fillingFactor": 0.22200177067914698, "insulationMaterial": null, "name": "Secondary section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [1, 1], "winding": "Secondary"}], "section": "Secondary section 0", "turnsAlignment": "outer or bottom", "type": "conduction", "windingStyle": "windByConsecutiveParallels"}], "sectionsDescription": [{"coordinateSystem": "cartesian", "coordinates": [0.004875540921738127, -0.0014062832242030996], "dimensions": [0.0019260000000000006, 0.004746103000000001], "fillingFactor": 0.7292119353768599, "layersAlignment": null, "layersOrientation": "overlapping", "margin": [0, 0], "name": "Primary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [1], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"coordinateSystem": "cartesian", "coordinates": [0.004956270460869064, -0.0037918347242031], "dimensions": [0.002087459078261874, 2.5e-05], "fillingFactor": 1, "layersAlignment": null, "layersOrientation": "contiguous", "margin": null, "name": "Insulation between Primary and Primary section 1", "partialWindings": [], "type": "insulation", "windingStyle": null}, {"coordinateSystem": "cartesian", "coordinates": [0.0057295000000005016, -0.0043453347242020995], "dimensions": [0.0005409999999990002, 0.0010819999999979995], "fillingFactor": 0.05753547898874283, "layersAlignment": null, "layersOrientation": "overlapping", "margin": [0, 0], "name": "Secondary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [1, 1], "winding": "Secondary"}], "type": "conduction", "windingStyle": "windByConsecutiveParallels"}], "turnsDescription": [{"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.0007260182757969007], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 0", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, 0.00011678927579690046], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 1", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, -0.0004924397242030997], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 2", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, -0.0011016687242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 3", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, -0.0017108977242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 4", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, -0.0023201267242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 5", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, -0.0029293557242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 6", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004153290921738127, -0.0035385847242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 0", "length": 0.02609589649590736, "name": "Primary parallel 0 turn 7", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.0007260182757969007], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 8", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, 0.00011678927579690046], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 9", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, -0.0004924397242030997], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 10", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, -0.0011016687242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 11", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, -0.0017108977242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 12", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, -0.0023201267242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 13", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, -0.0029293557242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 14", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.004634790921738127, -0.0035385847242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 1", "length": 0.029121250221314333, "name": "Primary parallel 0 turn 15", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.0007260182757969007], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 16", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, 0.00011678927579690046], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 17", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, -0.0004924397242030997], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 18", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, -0.0011016687242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 19", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, -0.0017108977242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 20", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, -0.0023201267242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 21", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, -0.0029293557242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 22", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0051162909217381276, -0.0035385847242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 2", "length": 0.03214660394672131, "name": "Primary parallel 0 turn 23", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.0007260182757969007], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 24", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, 0.00011678927579690046], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 25", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, -0.0004924397242030997], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 26", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, -0.0011016687242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 27", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, -0.0017108977242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 28", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, -0.0023201267242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 29", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, -0.0029293557242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 30", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005597790921738128, -0.0035385847242031], "dimensions": [0.0004815, 0.0004815], "layer": "Primary section 0 layer 3", "length": 0.03517195767212828, "name": "Primary parallel 0 turn 31", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0057295000000005016, -0.0040748347242026], "dimensions": [0.0005409999999990001, 0.0005409999999990001], "layer": "Secondary section 0 layer 0", "length": 0.03599951021748859, "name": "Secondary parallel 0 turn 0", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Secondary section 0", "winding": "Secondary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0057295000000005016, -0.0046158347242016], "dimensions": [0.0005409999999990001, 0.0005409999999990001], "layer": "Secondary section 0 layer 0", "length": 0.03599951021748859, "name": "Secondary parallel 1 turn 0", "orientation": "clockwise", "parallel": 1, "rotation": 0, "section": "Secondary section 0", "winding": "Secondary"}]}, "core": {"distributorsInfo": [], "functionalDescription": {"coating": null, "gapping": [{"area": 2.6e-05, "coordinates": [0, 0, 0], "distanceClosestNormalSurface": 0.005443, "distanceClosestParallelSurface": 0.0031625000000000004, "length": 0.000414, "sectionDimensions": [0.005675, 0.005675], "shape": "round", "type": "subtractive"}, {"area": 0.000163, "coordinates": [0, 0, -0.007038], "distanceClosestNormalSurface": 0.005648, "distanceClosestParallelSurface": 0.0031625000000000004, "length": 5e-06, "sectionDimensions": [0.078555, 0.002075], "shape": "irregular", "type": "residual"}], "material": "3C91", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0184, "minimum": 0.0176, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0085, "minimum": 0.0083, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01125, "minimum": 0.01075, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0058, "minimum": 0.0055, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0124, "minimum": 0.0116, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00585, "minimum": 0.0055, "nominal": null}, "K": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.002925, "minimum": null, "nominal": null}}, "family": "ep", "familySubtype": null, "magneticCircuit": "open", "name": "EP 17", "type": "standard"}, "type": "two-piece set"}, "geometricalDescription": [{"coordinates": [0, 0, 0], "dimensions": null, "insulationMaterial": null, "machining": [{"coordinates": [0, 0.0001035, 0], "length": 0.000207}], "material": "3C91", "rotation": [3.141592653589793, 3.141592653589793, 0], "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0184, "minimum": 0.0176, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0085, "minimum": 0.0083, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01125, "minimum": 0.01075, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0058, "minimum": 0.0055, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0124, "minimum": 0.0116, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00585, "minimum": 0.0055, "nominal": null}, "K": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.002925, "minimum": null, "nominal": null}}, "family": "ep", "familySubtype": null, "magneticCircuit": "open", "name": "EP 17", "type": "standard"}, "type": "half set"}, {"coordinates": [0, 0, 0], "dimensions": null, "insulationMaterial": null, "machining": [{"coordinates": [0, -0.0001035, 0], "length": 0.000207}], "material": "3C91", "rotation": [0, 0, 0], "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0184, "minimum": 0.0176, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0085, "minimum": 0.0083, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.01125, "minimum": 0.01075, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0058, "minimum": 0.0055, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0124, "minimum": 0.0116, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00585, "minimum": 0.0055, "nominal": null}, "K": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.002925, "minimum": null, "nominal": null}}, "family": "ep", "familySubtype": null, "magneticCircuit": "open", "name": "EP 17", "type": "standard"}, "type": "half set"}], "manufacturerInfo": {"cost": null, "datasheetUrl": "https://ferroxcube.com/upload/media/product/file/Pr_ds/EP17.pdf", "family": null, "name": "Ferroxcube", "orderCode": null, "reference": "EP17-3C91-A100", "status": "production"}, "name": "EP 17 - 3C91 - Gapped 0.414 mm", "processedDescription": {"columns": [{"area": 2.6e-05, "coordinates": [0, 0, 0], "depth": 0.005675, "height": 0.0113, "minimumDepth": null, "minimumWidth": null, "shape": "round", "type": "central", "width": 0.005675}, {"area": 0.000163, "coordinates": [0, 0, -0.007038], "depth": 0.002075, "height": 0.0113, "minimumDepth": null, "minimumWidth": 0.003001, "shape": "irregular", "type": "lateral", "width": 0.078555}], "depth": 0.011, "effectiveParameters": {"effectiveArea": 3.4461049794381896e-05, "effectiveLength": 0.02849715978271948, "effectiveVolume": 9.820420422707531e-07, "minimumArea": 2.5790801226066944e-05}, "height": 0.016800000000000002, "width": 0.018000000000000002, "windingWindows": [{"angle": null, "area": 3.573625e-05, "coordinates": [0.0028374999999999997, 0], "height": 0.0113, "radialHeight": null, "sectionsAlignment": null, "sectionsOrientation": null, "shape": null, "width": 0.0031625000000000004}]}}})"));
        OpenMagnetics::OperatingPoint operatingPoint(json::parse(R"({"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42}, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746.40291779551, "effectiveFrequency": 110746.40291779551, "peak": 5, "rms": 2.8874560332150576, "thd": 0.12151487440704967}, "harmonics": {"amplitudes": [1.1608769501236793e-14, 4.05366124583194, 1.787369544444173e-15, 0.4511310569983995, 9.749053004706756e-16, 0.16293015292554872, 4.036157626725542e-16, 0.08352979924600704, 3.4998295008010614e-16, 0.0508569581336163, 3.1489164048780735e-16, 0.034320410449418075, 3.142469873118059e-16, 0.024811988673843106, 2.3653352035940994e-16, 0.018849001010678823, 2.9306524147249266e-16, 0.014866633059596499, 1.796485796132569e-16, 0.012077180559557785, 1.6247782523152451e-16, 0.010049063750920609, 1.5324769149805092e-16, 0.008529750975091871, 1.0558579733068502e-16, 0.007363501410705499, 7.513269775674661e-17, 0.006450045785294609, 5.871414177162291e-17, 0.005722473794997712, 9.294731722001391e-17, 0.005134763398167541, 1.194820309200107e-16, 0.004654430423785411, 8.2422739080512e-17, 0.004258029771397032, 9.5067306351894e-17, 0.0039283108282380024, 1.7540347128474968e-16, 0.0036523670873925395, 9.623794010508822e-17, 0.0034204021424253787, 1.4083470894369491e-16, 0.0032248884817922927, 1.4749333016985644e-16, 0.0030599828465501895, 1.0448590642474364e-16, 0.002921112944200383, 7.575487373767413e-18, 0.002804680975178716, 7.419510610361002e-17, 0.0027078483284668376, 3.924741709073613e-17, 0.0026283777262804953, 2.2684279102637236e-17, 0.0025645167846443107, 8.997077625295079e-17, 0.0025149120164513483, 7.131074184849219e-17, 0.0024785457043284276, 9.354417496250849e-17, 0.0024546904085875065, 1.2488589642405877e-16, 0.0024428775264784264], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485.5360118389, "effectiveFrequency": 553357.3374711702, "peak": 70.5, "rms": 51.572309672924284, "thd": 0.4833151484524849}, "harmonics": {"amplitudes": [24.2890625, 57.92076613061847, 1.421875, 19.27588907896988, 1.421875, 11.528257939603122, 1.421875, 8.194467538528329, 1.421875, 6.331896912839248, 1.421875, 5.137996046859012, 1.421875, 4.304077056139349, 1.421875, 3.6860723299088454, 1.421875, 3.207698601961777, 1.421875, 2.8247804703632298, 1.421875, 2.509960393415185, 1.421875, 2.2453859950684323, 1.421875, 2.01890737840567, 1.421875, 1.8219644341144872, 1.421875, 1.6483482744897402, 1.421875, 1.4934420157473332, 1.421875, 1.3537375367153817, 1.421875, 1.2265178099275544, 1.421875, 1.1096421410704556, 1.421875, 1.0013973584174929, 1.421875, 0.9003924136274832, 1.421875, 0.8054822382572133, 1.421875, 0.7157117294021269, 1.421875, 0.6302738400635857, 1.421875, 0.5484777114167545, 1.421875, 0.46972405216147894, 1.421875, 0.3934858059809043, 1.421875, 0.31929270856030145, 1.421875, 0.24671871675852053, 1.421875, 0.17537155450693565, 1.421875, 0.10488380107099537, 1.421875, 0.034905072061178544], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"ancillaryLabel": null, "data": [-5, 5, -5], "numberPeriods": null, "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746.40291779551, "effectiveFrequency": 110746.40291779551, "peak": 5, "rms": 2.8874560332150576, "thd": 0.12151487440704967}, "harmonics": {"amplitudes": [1.1608769501236793e-14, 4.05366124583194, 1.787369544444173e-15, 0.4511310569983995, 9.749053004706756e-16, 0.16293015292554872, 4.036157626725542e-16, 0.08352979924600704, 3.4998295008010614e-16, 0.0508569581336163, 3.1489164048780735e-16, 0.034320410449418075, 3.142469873118059e-16, 0.024811988673843106, 2.3653352035940994e-16, 0.018849001010678823, 2.9306524147249266e-16, 0.014866633059596499, 1.796485796132569e-16, 0.012077180559557785, 1.6247782523152451e-16, 0.010049063750920609, 1.5324769149805092e-16, 0.008529750975091871, 1.0558579733068502e-16, 0.007363501410705499, 7.513269775674661e-17, 0.006450045785294609, 5.871414177162291e-17, 0.005722473794997712, 9.294731722001391e-17, 0.005134763398167541, 1.194820309200107e-16, 0.004654430423785411, 8.2422739080512e-17, 0.004258029771397032, 9.5067306351894e-17, 0.0039283108282380024, 1.7540347128474968e-16, 0.0036523670873925395, 9.623794010508822e-17, 0.0034204021424253787, 1.4083470894369491e-16, 0.0032248884817922927, 1.4749333016985644e-16, 0.0030599828465501895, 1.0448590642474364e-16, 0.002921112944200383, 7.575487373767413e-18, 0.002804680975178716, 7.419510610361002e-17, 0.0027078483284668376, 3.924741709073613e-17, 0.0026283777262804953, 2.2684279102637236e-17, 0.0025645167846443107, 8.997077625295079e-17, 0.0025149120164513483, 7.131074184849219e-17, 0.0024785457043284276, 9.354417496250849e-17, 0.0024546904085875065, 1.2488589642405877e-16, 0.0024428775264784264], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-50, 50, 50, -50, -50], "numberPeriods": null, "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485.536011839, "effectiveFrequency": 591449.4202715514, "peak": 50, "rms": 50, "thd": 0.48331514845248497}, "harmonics": {"amplitudes": [0.78125, 63.64919355013018, 1.5625, 21.18229569117569, 1.5625, 12.668415318245188, 1.5625, 9.004909382998164, 1.5625, 6.958128475647527, 1.5625, 5.646149502042871, 1.5625, 4.729755006746538, 1.5625, 4.050628933965765, 1.5625, 3.524943518639316, 1.5625, 3.104154363036517, 1.5625, 2.7581982345221827, 1.5625, 2.467457137437843, 1.5625, 2.2185795367095267, 1.5625, 2.0021587188071255, 1.5625, 1.8113717302085082, 1.5625, 1.6411450722498175, 1.5625, 1.487623666720196, 1.5625, 1.3478217691511587, 1.5625, 1.2193869682092893, 1.5625, 1.100436657601639, 1.5625, 0.9894422127774558, 1.5625, 0.8851453167661671, 1.5625, 0.7864964059364037, 1.5625, 0.6926086154544899, 1.5625, 0.60272275979863, 1.5625, 0.5161802771005264, 1.5625, 0.43240198459440116, 1.5625, 0.3508711083080249, 1.5625, 0.27111946896540395, 1.5625, 0.192715993963664, 1.5625, 0.11525692425384548, 1.5625, 0.03835722204524927], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}]})"));

        OpenMagnetics::Painter painter(outFile);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

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

    // TEST(Course0) {
    //     for (int64_t index = 1; index < 20; index++) {

    //         std::vector<int64_t> numberTurns = {index, index * 10};
    //         std::vector<int64_t> numberParallels = {1, 1};
    //         uint8_t interleavingLevel = 1;
    //         int64_t numberStacks = 1;
    //         std::string coreShape = "PQ 26/25";
    //         std::string coreMaterial = "3C97";
    //         auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

    //         auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel);
    //         auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

    //         auto outFile = outputFilePath;
    //         outFile.append("Course" + std::to_string(index) + ".svg");
    //         std::filesystem::remove(outFile);
    //         OpenMagnetics::Painter painter(outFile);
    //         OpenMagnetics::Magnetic magnetic;
    //         magnetic.set_core(core);
    //         magnetic.set_coil(coil);

    //         painter.paint_core(magnetic);
    //         painter.paint_bobbin(magnetic);
    //         painter.paint_coil_turns(magnetic);

    //         painter.export_svg();
    //         std::this_thread::sleep_for(std::chrono::milliseconds(200));
    //         CHECK(std::filesystem::exists(outFile));
    //         settings->reset();
    //     }
    // }
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
            outFile.append("Test_Wire_Painter_Litz_Simple_Insulated.svg");
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

    TEST(Test_Wire_Painter_Web_0) {
        OpenMagnetics::clear_databases(); 
        settings->set_painter_simple_litz(false);
        settings->set_painter_advanced_litz(true);

        OpenMagnetics::WireWrapper wire(json::parse(R"({"standard": "IEC 60317", "type": "litz", "strand": "Round 0.2 - Grade 1", "numberConductors": 17, "coating": {"breakdownVoltage": null, "grade": null, "material": null, "numberLayers": 1, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "served"}})"));

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
        double frequency = 250000;
        auto wire = OpenMagnetics::find_wire_by_name("Round 3.55 - Grade 1");
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Round_Enamelled_Grade_1_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Current_Density_Litz_Single_Served_Simple) {
        OpenMagnetics::clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000000;
        auto wire = OpenMagnetics::find_wire_by_name("Litz 10x0.02 - Grade 2 - Single Served");
        settings->set_painter_simple_litz(true);
        // settings->set_painter_advanced_litz(false);

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Litz_Single_Served_Simple_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Current_Density_Litz_Single_Served_Normal) {
        OpenMagnetics::clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 100e6;
        auto wire = OpenMagnetics::find_wire_by_name("Litz 8x0.1 - Grade 1 - Single Served");
        settings->set_painter_simple_litz(false);
        // settings->set_painter_advanced_litz(false);

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Litz_Single_Served_Normal_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Current_Density_Litz_Single_Served_Normal_Many_Strands) {
        OpenMagnetics::clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000000;
        auto wire = OpenMagnetics::find_wire_by_name("Litz 1000x0.02 - Grade 1 - Single Served");
        settings->set_painter_simple_litz(false);
        settings->set_painter_advanced_litz(true);

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Litz_Single_Served_Normal_Many_Strands_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Current_Density_Litz_Single_Served_Advanced) {
        OpenMagnetics::clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000000;
        auto wire = OpenMagnetics::find_wire_by_name("Litz 10x0.02 - Grade 2 - Single Served");
        settings->set_painter_simple_litz(false);
        settings->set_painter_advanced_litz(true);

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Litz_Single_Served_Advanced_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Current_Density_Rectangular_Enamelled_Grade_1) {
        OpenMagnetics::clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000;
        auto wire = OpenMagnetics::find_wire_by_name("Rectangular 2x0.80 - Grade 1");
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Rectangular_Enamelled_Grade_1_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings->reset();
    }
    TEST(Test_Wire_Painter_Current_Density_Foil) {
        OpenMagnetics::clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000;
        auto wire = OpenMagnetics::find_wire_by_name("Foil 0.5");
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithTolerance);
        wire.set_outer_width(wire.get_conducting_width().value());
        wire.set_outer_height(dimensionWithTolerance);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency, 0.001, 25, OpenMagnetics::WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Foil_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Current_Density_Web_0) {
        OpenMagnetics::clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 100e6;
        OpenMagnetics::WireWrapper wire(json::parse(R"({"coating": {"breakdownVoltage": 6000, "grade": null, "material": "TCA", "numberLayers": 3, "temperatureRating": 155, "thickness": null, "thicknessLayers": 7.62e-05, "type": "insulated"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00026900000000000003, "minimum": 0.000261, "nominal": 0.000265}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.265 - Grade 1", "numberConductors": 92, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000297, "minimum": 0.000283, "nominal": 0.0018383994787140629}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": "0.265 mm", "strand": {"coating": {"breakdownVoltage": 500, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00010300000000000001, "minimum": 9.7e-05, "nominal": 0.0001}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.1 - Grade 1", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00011700000000000001, "minimum": 0.000108, "nominal": null}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": "0.1 mm", "strand": null, "type": "round"}, "type": "litz"})"));
        OpenMagnetics::OperatingPoint operatingPoint(json::parse(R"({"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42}, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000000, "current": {"waveform": {"ancillaryLabel": null, "data": [-5, 5, -5], "numberPeriods": null, "time": [0, 5e-09, 1e-08]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746402.91779621, "effectiveFrequency": 110746402.91779621, "peak": 5, "rms": 2.887456033215047, "thd": 0.12151487440705}, "harmonics": {"amplitudes": [1.9671764217576992e-15, 4.053661245831925, 9.528877815088396e-15, 0.4511310569983986, 3.867859480230825e-15, 0.16293015292554885, 2.6154487633755876e-15, 0.08352979924600734, 1.9651652434930284e-15, 0.05085695813361692, 1.6136051601223988e-15, 0.034320410449418734, 1.3808009023999193e-15, 0.024811988673843956, 1.1386794830042506e-15, 0.018849001010679767, 1.0370440806127272e-15, 0.014866633059597558, 8.851494177631955e-16, 0.012077180559558507, 9.326505898457337e-16, 0.010049063750921247, 8.13280272912139e-16, 0.008529750975092641, 7.240094812312869e-16, 0.007363501410706191, 6.482105482773775e-16, 0.006450045785295375, 6.066430502657454e-16, 0.0057224737949983576, 5.440719337350025e-16, 0.005134763398167652, 6.25059167873935e-16, 0.004654430423786248, 4.892419283154727e-16, 0.004258029771397705, 5.345312821060954e-16, 0.0039283108282387675, 4.672784117529907e-16, 0.003652367087393178, 3.3833965320162966e-16, 0.003420402142426004, 4.2496275203954786e-16, 0.0032248884817931253, 3.8603919165942636e-16, 0.0030599828465507064, 3.809779240559984e-16, 0.0029211129442011045, 4.459588218866508e-16, 0.002804680975179318, 3.8344856196395527e-16, 0.0027078483284674916, 3.809203930919387e-16, 0.0026283777262812655, 4.0712748931445545e-16, 0.0025645167846449803, 3.9966373222574786e-16, 0.0025149120164521047, 4.572347145426278e-16, 0.0024785457043292464, 3.6473608081573336e-16, 0.0024546904085880616, 3.509577233239641e-16, 0.0024428775264793146], "frequencies": [0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000, 900000000, 1000000000, 1100000000, 1200000000, 1300000000, 1400000000, 1500000000, 1600000000, 1700000000, 1800000000, 1900000000, 2000000000, 2100000000, 2200000000, 2300000000, 2400000000, 2500000000, 2600000000, 2700000000, 2800000000, 2900000000, 3000000000, 3100000000, 3200000000, 3300000000, 3400000000, 3500000000, 3600000000, 3700000000, 3800000000, 3900000000, 4000000000, 4100000000, 4200000000, 4300000000, 4400000000, 4500000000, 4600000000, 4700000000, 4800000000, 4900000000, 5000000000, 5100000000, 5200000000, 5300000000, 5400000000, 5500000000, 5600000000, 5700000000, 5800000000, 5900000000, 6000000000, 6100000000, 6200000000, 6300000000]}}, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-50, 50, 50, -50, -50], "numberPeriods": null, "time": [0, 0, 5e-09, 5e-09, 1e-08]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485536.0118392, "effectiveFrequency": 591449420.2715513, "peak": 50, "rms": 50, "thd": 0.48331514845248497}, "harmonics": {"amplitudes": [0.78125, 63.64919355013018, 1.5625, 21.18229569117569, 1.5625, 12.668415318245188, 1.5625, 9.004909382998164, 1.5625, 6.958128475647527, 1.5625, 5.646149502042871, 1.5625, 4.729755006746538, 1.5625, 4.050628933965765, 1.5625, 3.524943518639316, 1.5625, 3.104154363036517, 1.5625, 2.7581982345221827, 1.5625, 2.467457137437843, 1.5625, 2.2185795367095267, 1.5625, 2.0021587188071255, 1.5625, 1.8113717302085082, 1.5625, 1.6411450722498175, 1.5625, 1.487623666720196, 1.5625, 1.3478217691511587, 1.5625, 1.2193869682092893, 1.5625, 1.100436657601639, 1.5625, 0.9894422127774558, 1.5625, 0.8851453167661671, 1.5625, 0.7864964059364037, 1.5625, 0.6926086154544899, 1.5625, 0.60272275979863, 1.5625, 0.5161802771005264, 1.5625, 0.43240198459440116, 1.5625, 0.3508711083080249, 1.5625, 0.27111946896540395, 1.5625, 0.192715993963664, 1.5625, 0.11525692425384548, 1.5625, 0.03835722204524927], "frequencies": [0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000, 900000000, 1000000000, 1100000000, 1200000000, 1300000000, 1400000000, 1500000000, 1600000000, 1700000000, 1800000000, 1900000000, 2000000000, 2100000000, 2200000000, 2300000000, 2400000000, 2500000000, 2600000000, 2700000000, 2800000000, 2900000000, 3000000000, 3100000000, 3200000000, 3300000000, 3400000000, 3500000000, 3600000000, 3700000000, 3800000000, 3900000000, 4000000000, 4100000000, 4200000000, 4300000000, 4400000000, 4500000000, 4600000000, 4700000000, 4800000000, 4900000000, 5000000000, 5100000000, 5200000000, 5300000000, 5400000000, 5500000000, 5600000000, 5700000000, 5800000000, 5900000000, 6000000000, 6100000000, 6200000000, 6300000000]}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}]})"));
        settings->set_painter_simple_litz(false);
        // settings->set_painter_advanced_litz(false);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Web_0.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire_with_current_density(wire, operatingPoint);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Current_Density_Web_1) {
        OpenMagnetics::clear_databases();

        OpenMagnetics::WireWrapper wire(json::parse(R"({"coating": {"breakdownVoltage": 13500, "grade": null, "material": "FEP", "numberLayers": 3, "temperatureRating": 155, "thickness": null, "thicknessLayers": 7.62e-05, "type": "insulated"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 1.0602875205865554e-06}, "conductingDiameter": null, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": null, "name": "Litz 15x0.3 - Grade 1 - Unserved", "numberConductors": 18, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0016300000000000002, "minimum": 0.0011557, "nominal": 0.012585990719707143}, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": null, "strand": {"coating": {"breakdownVoltage": 2500, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.004115}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Nearson", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 6.0 - Single Build", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.004186}, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": "6 AWG", "strand": null, "type": "round"}, "type": "litz"})"));
        settings->set_painter_simple_litz(false);
        // settings->set_painter_advanced_litz(false);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Web_1.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Wire_Painter_Current_Density_Web_2) {
        OpenMagnetics::clear_databases();

        // OpenMagnetics::WireWrapper wire(json::parse(R"({"coating": {"breakdownVoltage": null, "grade": null, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "bare"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 1.0602875205865554e-06}, "conductingDiameter": null, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": null, "name": "Litz 15x0.3 - Grade 1 - Unserved", "numberConductors": 15, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0016300000000000002, "minimum": 0.0011557, "nominal": null}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": null, "strand": {"coating": {"breakdownVoltage": 2200, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000303999999999, "minimum": 0.000296, "nominal": 0.00030000000000000003}, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.3 - Grade 1", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00033400000000000004, "minimum": 0.000319, "nominal": null}, "standard": "IEC 60317", "standardName": "0.3 mm", "type": "round"}, "type": "litz"})"));
        OpenMagnetics::WireWrapper wire(json::parse(R"({"coating": {"breakdownVoltage": null, "grade": null, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "bare"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 1.0602875205865554e-06}, "conductingDiameter": null, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": null, "name": "Litz 15x0.3 - Grade 1 - Unserved", "numberConductors": 15, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0016300000000000002, "minimum": 0.0011557, "nominal": 0.0015933066187962695}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": null, "strand": {"coating": {"breakdownVoltage": 2200, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000303999999999, "minimum": 0.000296, "nominal": 0.00030000000000000003}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.3 - Grade 1", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00033400000000000004, "minimum": 0.000319, "nominal": null}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": "0.3 mm", "strand": null, "type": "round"}, "type": "litz"})"));
        settings->set_painter_simple_litz(false);
        // settings->set_painter_advanced_litz(false);

        {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Web_2.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings->reset();
    }

}