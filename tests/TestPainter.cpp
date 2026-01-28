#include "support/Painter.h"
#include "json.hpp"
#include "TestingUtils.h"
#include <source_location>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fstream>
#include <string>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <chrono>
#include <thread>

using namespace MAS;
using namespace OpenMagnetics;

namespace { 
    const auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");

    TEST_CASE("Test_Painter_Contour_Many_Turns", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        OpenMagneticsTesting::PainterTestConfig config;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Contour_Many_Turns_Logarithmic_Scale", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        OpenMagneticsTesting::PainterTestConfig config;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_Logarithmic_Scale.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(true);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Contour_Many_Turns_No_Fringing", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        OpenMagneticsTesting::PainterTestConfig config;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);
 
        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_No_Fringing.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Contour_Many_Turns_Limit_Scale", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        OpenMagneticsTesting::PainterTestConfig config;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_Limit_Scale.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(25500);
        settings.set_painter_minimum_value_colorbar(0);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Contour_One_Turn", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
        OpenMagneticsTesting::PainterTestConfig config;
        config.numberTurns = {1};
        config.numberParallels = {1};
        config.interleavingLevel = 1;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_One_Turn.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(true);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Quiver_Many_Turns", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        OpenMagneticsTesting::PainterTestConfig config;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Quiver_One_Turn", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
        OpenMagneticsTesting::PainterTestConfig config;
        config.numberTurns = {1};
        config.numberParallels = {1};
        config.interleavingLevel = 1;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_One_Turn.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(true);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Quiver_Many_Turns_No_Fringing", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        OpenMagneticsTesting::PainterTestConfig config;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_No_Fringing.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_include_fringing(false);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Quiver_Many_Turns_Logarithmic_Scale", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        OpenMagneticsTesting::PainterTestConfig config;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_Logarithmic_Scale.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(true);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Quiver_Many_Turns_Limit_Scale", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        OpenMagneticsTesting::PainterTestConfig config;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_Limit_Scale.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(2500);
        settings.set_painter_minimum_value_colorbar(0);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Quiver_One_Turn_Rectangular", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
        OpenMagneticsTesting::PainterTestConfig config;
        config.numberTurns = {1};
        config.numberParallels = {1};
        config.interleavingLevel = 1;
        config.wireNames = {"Rectangular 2.36x1.12 - Grade 1"};
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_One_Turn_Rectangular.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Contour_One_Turn_Rectangular", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
        OpenMagneticsTesting::PainterTestConfig config;
        config.numberTurns = {1};
        config.numberParallels = {1};
        config.interleavingLevel = 1;
        config.wireNames = {"Rectangular 2.36x1.12 - Grade 1"};
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_One_Turn_Rectangular.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Quiver_Many_Turns_Rectangular", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
        OpenMagneticsTesting::PainterTestConfig config;
        config.numberTurns = {10};
        config.numberParallels = {1};
        config.interleavingLevel = 1;
        config.turnsAlignment = CoilAlignment::SPREAD;
        config.wireNames = {"Rectangular 4.50x0.90 - Grade 1"};
        config.compactCoil = false;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_Rectangular.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Contour_Many_Turns_Rectangular", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
        OpenMagneticsTesting::PainterTestConfig config;
        config.numberTurns = {10};
        config.numberParallels = {1};
        config.interleavingLevel = 1;
        config.turnsAlignment = CoilAlignment::SPREAD;
        config.wireNames = {"Rectangular 4.50x0.90 - Grade 1"};
        config.compactCoil = false;
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_Rectangular.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Quiver_One_Turn_Foil", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
        // Custom wire setup - modify Foil wire dimensions
        OpenMagnetics::Wire wire = find_wire_by_name("Foil 0.15");
        DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithTolerance);
        wire.set_outer_width(wire.get_conducting_width().value());
        wire.set_outer_height(dimensionWithTolerance);

        OpenMagneticsTesting::PainterTestConfig config;
        config.numberTurns = {1};
        config.numberParallels = {1};
        config.interleavingLevel = 1;
        config.customWires = {wire};
        auto [magnetic, inputs] = OpenMagneticsTesting::prepare_painter_test(config);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_One_Turn_Foil.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Contour_One_Turn_Foil", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
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
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        OpenMagnetics::Wire wire = find_wire_by_name("Foil 0.15");
        DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithTolerance);
        wire.set_outer_width(wire.get_conducting_width().value());
        wire.set_outer_height(dimensionWithTolerance);
        auto wires = std::vector<OpenMagnetics::Wire>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_One_Turn_Foil.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Quiver_Many_Turns_Foil", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
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
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        OpenMagnetics::Wire wire = find_wire_by_name("Foil 0.15");
        DimensionWithTolerance dimensionWithToleranceHeight;
        dimensionWithToleranceHeight.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithToleranceHeight);
        DimensionWithTolerance dimensionWithToleranceWidth;
        dimensionWithToleranceWidth.set_nominal(0.2e-3);
        wire.set_outer_width(dimensionWithToleranceWidth);
        wire.set_outer_height(dimensionWithToleranceHeight);
        auto wires = std::vector<OpenMagnetics::Wire>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Quiver_Many_Turns_Foil.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Contour_Many_Turns_Foil", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
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
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        OpenMagnetics::Wire wire = find_wire_by_name("Foil 0.15");
        DimensionWithTolerance dimensionWithToleranceHeight;
        dimensionWithToleranceHeight.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithToleranceHeight);
        DimensionWithTolerance dimensionWithToleranceWidth;
        dimensionWithToleranceWidth.set_nominal(0.2e-3);
        wire.set_outer_width(dimensionWithToleranceWidth);
        wire.set_outer_height(dimensionWithToleranceHeight);
        auto wires = std::vector<OpenMagnetics::Wire>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Contour_Many_Turns_Foil.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Text_Color", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
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
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        OpenMagnetics::Wire wire = find_wire_by_name("Foil 0.15");
        DimensionWithTolerance dimensionWithToleranceHeight;
        dimensionWithToleranceHeight.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithToleranceHeight);
        DimensionWithTolerance dimensionWithToleranceWidth;
        dimensionWithToleranceWidth.set_nominal(0.2e-3);
        wire.set_outer_width(dimensionWithToleranceWidth);
        wire.set_outer_height(dimensionWithToleranceHeight);
        auto wires = std::vector<OpenMagnetics::Wire>({wire});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Text_Color.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        settings.set_painter_color_text("0xFF0000");
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Coil_Magnetic_Field_Painter_Web_0", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
        clear_databases();

        auto json_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_magnetic_field_painter_web_0.json");
        std::ifstream json_file(json_path);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file));
        auto json_path_609 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_magnetic_field_painter_web_0_609.json");
        std::ifstream json_file_609(json_path_609);
        MAS::OperatingPoint operatingPoint(json::parse(json_file_609));
        // settings.set_painter_simple_litz(false);
        // settings.set_painter_advanced_litz(false);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Magnetic_Field_Painter_Web_0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);

            painter.paint_magnetic_field(operatingPoint, magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
        }
        settings.reset();
    }

    TEST_CASE("Test_Coil_Magnetic_Field_Basic_Painter", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {14, 16};
        std::vector<int64_t> numberParallels = {1, 1};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Round 0.80 - Grade 3")});
        wires.push_back({find_wire_by_name("Round 0.80 - Grade 3")});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Magnetic_Field_Basic_Painter.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);

            settings.set_painter_number_points_x(50);
            settings.set_painter_number_points_y(100);
            settings.set_painter_include_fringing(true);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
        }
        settings.reset();
    }

    TEST_CASE("Test_Coil_Magnetic_Field_Basic_Painter_Planar", "[support][painter][magnetic-field-painter][rectangular-winding-window][smoke-test]") {

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "leakage_inductance_planar.json");
        OpenMagnetics::Mas mas;
        OpenMagnetics::from_file(path, mas);
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Magnetic_Field_Basic_Painter_Planar.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);

            settings.set_painter_number_points_x(50);
            settings.set_painter_number_points_y(100);
            settings.set_painter_include_fringing(false);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
        }
        settings.reset();
    }

    TEST_CASE("Test_Coil_Electric_Field_Basic_Painter", "[support][painter][electric-field-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {14, 16};
        std::vector<int64_t> numberParallels = {1, 1};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Round 0.80 - Grade 3")});
        wires.push_back({find_wire_by_name("Round 0.80 - Grade 3")});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Electric_Field_Basic_Painter.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);

            settings.set_painter_number_points_x(50);
            settings.set_painter_number_points_y(100);
            settings.set_painter_include_fringing(false);
            painter.paint_electric_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
        }
        settings.reset();
    }

    TEST_CASE("Test_Coil_Electric_Field_Basic_Painter_Planar", "[support][painter][electric-field-painter][rectangular-winding-window][smoke-test]") {

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "leakage_inductance_planar.json");
        OpenMagnetics::Mas mas;
        OpenMagnetics::from_file(path, mas);
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Electric_Field_Basic_Painter_Planar.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);

            settings.set_painter_number_points_x(50);
            settings.set_painter_number_points_y(100);
            settings.set_painter_include_fringing(false);
            painter.paint_electric_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
        }
        settings.reset();
    }

    TEST_CASE("Test_Coil_Wire_Losses_Basic_Painter", "[support][painter][wire-losses-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {14, 16};
        std::vector<int64_t> numberParallels = {1, 1};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Round 0.80 - Grade 3")});
        wires.push_back({find_wire_by_name("Round 0.80 - Grade 3")});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Wire_Losses_Basic_Painter.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);

            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.paint_wire_losses(magnetic, std::nullopt, inputs.get_operating_point(0));
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
        }
        settings.reset();
    }


    TEST_CASE("Test_Coil_Wire_Losses_Basic_Painter_Planar", "[support][painter][wire-losses-painter][rectangular-winding-window][smoke-test]") {

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "leakage_inductance_planar.json");
        OpenMagnetics::Mas mas;
        OpenMagnetics::from_file(path, mas);
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Wire_Losses_Basic_Painter_Planar.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);

            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.paint_wire_losses(magnetic, std::nullopt, inputs.get_operating_point(0));
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
        }
        settings.reset();
    }

    TEST_CASE("Test_Painter_Toroid_Round_Wires", "[support][painter][magnetic-field-painter][round-winding-window]") {
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({100, 5});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 1e-3;
        std::string shapeName = "T 20/10/7";

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
                                                                                         offset,
                                                                                         turnsRatios);

        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::INNER_OR_TOP;
        auto sectionsAlignment = CoilAlignment::INNER_OR_TOP;

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


        settings.reset();

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_Toroid_Round_Wires.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            // settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_mode(PainterModes::QUIVER);
            settings.set_painter_logarithmic_scale(false);
            settings.set_painter_include_fringing(true);
            settings.set_painter_number_points_x(50);
            settings.set_painter_number_points_y(50);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }

    TEST_CASE("Test_Painter_Toroid_Quiver_One_Turn_Rectangular", "[support][painter][magnetic-field-painter][round-winding-window]") {
        clear_databases();
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_One_Turn_Rectangular.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        settings.set_painter_number_points_x(100);
        settings.set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Toroid_Quiver_One_Turn_Rectangular_Inner", "[support][painter][magnetic-field-painter][round-winding-window]") {
        clear_databases();
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_One_Turn_Rectangular_Inner.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        settings.set_painter_number_points_x(100);
        settings.set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Toroid_Quiver_Four_Turns_Rectangular_Inner", "[support][painter][magnetic-field-painter][round-winding-window]") {
        clear_databases();
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_Four_Turns_Rectangular_Inner.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        settings.set_painter_number_points_x(100);
        settings.set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Toroid_Quiver_Four_Turns_Rectangular_Spread", "[support][painter][magnetic-field-painter][round-winding-window]") {
        clear_databases();
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_Four_Turns_Rectangular_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        settings.set_painter_number_points_x(100);
        settings.set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Toroid_Quiver_Two_Turn_Rectangular", "[support][painter][magnetic-field-painter][round-winding-window]") {
        clear_databases();
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "T 20/10/7";

        std::string coreMaterial = "3C97";
        auto emptyGapping = json::array();
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;

        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Toroid_Quiver_Two_Turn_Rectangular.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        settings.set_painter_mode(PainterModes::QUIVER);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(false);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        settings.set_painter_number_points_x(100);
        settings.set_painter_number_points_y(100);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        // painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Toroid_Rectangular_Wires", "[support][painter][magnetic-field-painter][round-winding-window]") {
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns = {11, 90};
        std::vector<int64_t> numberParallels = {1, 1};
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        uint8_t interleavingLevel = 1;
        std::string coreShape = "T 20/10/7";
        auto emptyGapping = json::array();
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
        std::vector<OpenMagnetics::Wire> wires;

        wires.push_back({find_wire_by_name("Rectangular 2.50x1.18 - Grade 1")});
        wires.push_back({find_wire_by_name("Round 0.335 - Grade 1")});


        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 1e-3;
        std::string shapeName = "T 20/10/7";

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


        settings.reset();

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_Toroid_Rectangular_Wires.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            // settings.set_painter_mode(PainterModes::CONTOUR);
            settings.set_painter_mode(PainterModes::QUIVER);
            settings.set_painter_logarithmic_scale(false);
            settings.set_painter_include_fringing(true);
            settings.set_painter_number_points_x(50);
            settings.set_painter_number_points_y(50);
            settings.set_painter_maximum_value_colorbar(std::nullopt);
            settings.set_painter_minimum_value_colorbar(std::nullopt);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }

    TEST_CASE("Test_Painter_T_Core", "[support][painter][magnetic-painter][round-winding-window][smoke-test]") { 
        clear_databases();
        settings.set_coil_try_rewind(false);
        std::vector<int64_t> numberTurns = {12, 12};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "T 80/20/50";
        std::string coreMaterial = "3C97"; 
        settings.set_coil_delimit_and_compact(false);
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel,
                                                        WindingOrientation::OVERLAPPING,
                                                        WindingOrientation::OVERLAPPING,
                                                        CoilAlignment::SPREAD,
                                                        CoilAlignment::SPREAD
                                                        );
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_T_Core.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        // painter.paint_coil_sections(magnetic);
        // painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_T_Core_Overlapping", "[support][painter][magnetic-painter][round-winding-window][smoke-test]") {
        clear_databases();
        std::vector<int64_t> numberTurns = {2, 2};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "T 80/20/50";
        std::string coreMaterial = "3C97"; 
        settings.set_coil_delimit_and_compact(false);
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, WindingOrientation::OVERLAPPING);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_T_Core_Overlapping.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_T_Core_Contiguous", "[support][painter][magnetic-painter][round-winding-window][smoke-test]") {
        clear_databases();
        std::vector<int64_t> numberTurns = {72, 72};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "T 80/20/50";
        std::string coreMaterial = "3C97";
        // settings.set_coil_delimit_and_compact(false);
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, WindingOrientation::CONTIGUOUS);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_T_Core_Contiguous_turns.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, false);
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
            Painter painter(outFile, false, false, false);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_layers(magnetic);
            // painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
            settings.reset();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_T_Core_Contiguous_sections.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, false);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_coil_sections(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
            settings.reset();
        }
    }

    TEST_CASE("Test_Painter_T_Core_Contiguous_Sections_With_Margin", "[support][painter][magnetic-painter][round-winding-window][smoke-test]") {
        clear_databases();
        std::vector<int64_t> numberTurns = {2, 2};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "T 80/20/50";
        std::string coreMaterial = "3C97"; 
        settings.set_coil_delimit_and_compact(false);
        settings.set_painter_draw_spacer(false);
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, WindingOrientation::CONTIGUOUS);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        double margin = 0.005;
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 1.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        coil.add_margin_to_section_by_index(3, std::vector<double>{margin * 3.5, margin * 3.5});

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_T_Core_Contiguous_Sections_With_Margin.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_sections(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_T_Core_Contiguous_Sections_With_Spacer", "[support][painter][magnetic-painter][round-winding-window][smoke-test]") {
        clear_databases();
        std::vector<int64_t> numberTurns = {2, 2};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "T 80/20/50";
        std::string coreMaterial = "3C97"; 
        settings.set_coil_delimit_and_compact(false);
        settings.set_painter_draw_spacer(true);
        auto emptyGapping = json::array();

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, WindingOrientation::CONTIGUOUS);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, emptyGapping, numberStacks, coreMaterial);

        double margin = 0.001;
        coil.add_margin_to_section_by_index(0, std::vector<double>{margin, margin});
        coil.add_margin_to_section_by_index(1, std::vector<double>{margin * 1.5, margin * 0.5});
        coil.add_margin_to_section_by_index(2, std::vector<double>{margin * 0.5, margin * 3.5});
        coil.add_margin_to_section_by_index(3, std::vector<double>{margin * 3.5, margin * 3.5});

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_T_Core_Contiguous_Sections_With_Spacer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_sections(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Coil_Painter_Web_0", "[support][painter][magnetic-painter][round-winding-window][smoke-test]") {
        clear_databases();

        OpenMagnetics::Magnetic magnetic(json::parse(R"({"coil": {"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.005974999999999999, "columnShape": "rectangular", "columnThickness": 0.0013999999999999993, "columnWidth": 0.005999999999999999, "coordinates": [0, 0, 0], "pins": null, "wallThickness": 0.0011204000000000006, "windingWindows": [{"angle": null, "area": 0.00011625151999999999, "coordinates": [0.008799999999999999, 0, 0], "height": 0.0207592, "radialHeight": null, "sectionsOrientation": "overlapping", "shape": "rectangular", "width": 0.0056}]}}, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "primary", "numberParallels": 1, "numberTurns": 20, "wire": {"coating": {"breakdownVoltage": 18000, "grade": null, "material": "PFA", "numberLayers": 3, "temperatureRating": 180, "thickness": null, "thicknessLayers": 7.62e-05, "type": "insulated"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 6.532502100168472e-07}, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.000912}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Rubadue", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round T19A01PXXX-3", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.001369}, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": "19 AWG", "strand": null, "type": "round"}}], "layersDescription": [{"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 8.673617379884035e-19], "dimensions": [0.0045850000000000005, 0.020700000000000003], "fillingFactor": 0.9971482523411308, "insulationMaterial": null, "name": "primary section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [1], "winding": "primary"}], "section": "primary section 0", "turnsAlignment": "centered", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}], "sectionsDescription": [{"coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0], "dimensions": [0.0045850000000000005, 0.020700000000000003], "fillingFactor": 0.9624513055163947, "layersAlignment": null, "layersOrientation": "overlapping", "margin": [0, 0], "name": "primary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [1], "winding": "primary"}], "type": "conduction", "windingStyle": "windByConsecutiveTurns"}], "turnsDescription": [{"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.009832500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 0", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.008797500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 1", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.007762500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 2", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.006727500000000003], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 3", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.005692500000000002], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 4", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.004657500000000002], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 5", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.0036225000000000016], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 6", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.0025875000000000013], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 7", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.0015525000000000011], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 8", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, 0.000517500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 9", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0005174999999999991], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 10", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0015524999999999992], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 11", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0025874999999999995], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 12", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0036225], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 13", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0046575], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 14", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.0056925000000000005], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 15", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.006727500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 16", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.007762500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 17", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.008797500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 18", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.008292499999999998, -0.009832500000000001], "dimensions": [0.0045850000000000005, 0.0010350000000000001], "layer": "primary section 0 layer 0", "length": 0.06230420231670919, "name": "primary parallel 0 turn 19", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "primary section 0", "winding": "primary"}]}, "core": {"distributorsInfo": [{"cost": 1.87, "country": "USA", "distributedArea": "International", "email": null, "link": "https://www.digikey.com/en/products/detail/epcos-tdk-electronics/B66229G0500X127/3914560", "name": "Digi-Key", "phone": null, "quantity": 660, "reference": "495-B66229G0500X127-ND", "updatedAt": "16/04/2024"}, {"cost": 1.86, "country": "USA", "distributedArea": "International", "email": null, "link": "https://www.mouser.es/ProductDetail/EPCOS-TDK/B66229G0500X127?qs=%2FsLciWRBLmC6djd1SV6djQ%3D%3D", "name": "Mouser", "phone": null, "quantity": 355, "reference": "871-B66229G0500X127", "updatedAt": "16/04/2024"}], "functionalDescription": {"coating": null, "gapping": [{"area": 8.5e-05, "coordinates": [0, 0.00025, 0], "distanceClosestNormalSurface": 0.011, "distanceClosestParallelSurface": 0.006999999999999999, "length": 0.0005, "sectionDimensions": [0.0092, 0.00915], "shape": "rectangular", "type": "subtractive"}, {"area": 4.1e-05, "coordinates": [0.013825, 0, 0], "distanceClosestNormalSurface": 0.011498, "distanceClosestParallelSurface": 0.006999999999999999, "length": 5e-06, "sectionDimensions": [0.004451, 0.00915], "shape": "rectangular", "type": "residual"}, {"area": 4.1e-05, "coordinates": [-0.013825, 0, 0], "distanceClosestNormalSurface": 0.011498, "distanceClosestParallelSurface": 0.006999999999999999, "length": 5e-06, "sectionDimensions": [0.004451, 0.00915], "shape": "rectangular", "type": "residual"}], "material": "N27", "numberStacks": 1, "shape": {"aliases": ["E 32/9", "EF 32"], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0329, "minimum": 0.0313, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0164, "minimum": 0.0158, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0095, "minimum": 0.0088, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0118, "minimum": 0.0112, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0237, "minimum": 0.0227, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0095, "minimum": 0.0089, "nominal": null}}, "family": "e", "familySubtype": null, "magneticCircuit": "open", "name": "E 32/16/9", "type": "standard"}, "type": "two-piece set"}, "geometricalDescription": null, "manufacturerInfo": {"cost": null, "datasheetUrl": "https://product.tdk.com/system/files/dam/doc/product/ferrite/ferrite/ferrite-acc/data_sheet/80/db/fer/e_32_16_9.pdf", "family": null, "name": "TDK", "orderCode": null, "reference": "B66229G0500X127 (N27", "status": "production"}, "name": "E 32/16/9 - N27 - Gapped 0.500 mm", "processedDescription": {"columns": [{"area": 8.5e-05, "coordinates": [0, 0, 0], "depth": 0.00915, "height": 0.023, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "central", "width": 0.0092}, {"area": 4.1e-05, "coordinates": [0.013825, 0, 0], "depth": 0.00915, "height": 0.023, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.004451}, {"area": 4.1e-05, "coordinates": [-0.013825, 0, 0], "depth": 0.00915, "height": 0.023, "minimumDepth": null, "minimumWidth": null, "shape": "rectangular", "type": "lateral", "width": 0.004451}], "depth": 0.00915, "effectiveParameters": {"effectiveArea": 8.316165619067432e-05, "effectiveLength": 0.07431657444154591, "effectiveVolume": 6.180289412976495e-06, "minimumArea": 8.143500000000005e-05}, "height": 0.032200000000000006, "width": 0.032100000000000004, "windingWindows": [{"angle": null, "area": 0.00016099999999999998, "coordinates": [0.0046, 0], "height": 0.023, "radialHeight": null, "sectionsOrientation": null, "shape": null, "width": 0.006999999999999999}]}}, "distributorsInfo": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "", "recommendations": null, "reference": "E 32/16/9 - N27 - 0.5 mm, Turns: 20, Order: 0, Non-Interleaved, Margin Taped 00", "status": null}, "rotation": null})"));
        // settings.set_painter_simple_litz(false);
        // settings.set_painter_advanced_litz(false);

        {
            auto outFile = outputFilePath;
            outFile.append("Test_Coil_Painter_Web_0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);

            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(std::filesystem::exists(outFile));
        }
        settings.reset();
    }

    TEST_CASE("Test_Coil_Painter_Web_1", "[support][painter][magnetic-painter][round-winding-window][smoke-test]") {
        clear_databases();

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_1.svg");
        auto json_path_1429 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_painter_web_1_1429.json");
        std::ifstream json_file_1429(json_path_1429);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_1429));

        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings.reset();
    }

    TEST_CASE("Test_Coil_Painter_Web_2", "[support][painter][magnetic-painter][round-winding-window][smoke-test]") {
        clear_databases();
        // settings.set_coil_delimit_and_compact(false);

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_2.svg");
        auto json_path_1446 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_painter_web_2_1446.json");
        std::ifstream json_file_1446(json_path_1446);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_1446));

        magnetic.get_mutable_coil().set_layers_orientation(WindingOrientation::OVERLAPPING, "Primary section 0");
        magnetic.get_mutable_coil().set_layers_orientation(WindingOrientation::CONTIGUOUS, "Secondary section 0");
        magnetic.get_mutable_coil().set_turns_alignment(CoilAlignment::INNER_OR_TOP, "Primary section 0");
        magnetic.get_mutable_coil().set_turns_alignment(CoilAlignment::OUTER_OR_BOTTOM, "Secondary section 0");

        magnetic.get_mutable_coil().clear();
        magnetic.get_mutable_coil().wind();
        // magnetic.get_mutable_coil().delimit_and_compact();

        Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings.reset();
    }

    TEST_CASE("Test_Coil_Painter_Web_3", "[support][painter][magnetic-painter][round-winding-window][smoke-test]") {
        clear_databases();
        // settings.set_coil_delimit_and_compact(false);

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_3.svg");
        auto json_path_1473 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_painter_web_3_1473.json");
        std::ifstream json_file_1473(json_path_1473);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_1473));
        OperatingPoint operatingPoint(json::parse(R"({"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42}, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746.40291779551, "effectiveFrequency": 110746.40291779551, "peak": 5, "rms": 2.8874560332150576, "thd": 0.12151487440704967}, "harmonics": {"amplitudes": [1.1608769501236793e-14, 4.05366124583194, 1.787369544444173e-15, 0.4511310569983995, 9.749053004706756e-16, 0.16293015292554872, 4.036157626725542e-16, 0.08352979924600704, 3.4998295008010614e-16, 0.0508569581336163, 3.1489164048780735e-16, 0.034320410449418075, 3.142469873118059e-16, 0.024811988673843106, 2.3653352035940994e-16, 0.018849001010678823, 2.9306524147249266e-16, 0.014866633059596499, 1.796485796132569e-16, 0.012077180559557785, 1.6247782523152451e-16, 0.010049063750920609, 1.5324769149805092e-16, 0.008529750975091871, 1.0558579733068502e-16, 0.007363501410705499, 7.513269775674661e-17, 0.006450045785294609, 5.871414177162291e-17, 0.005722473794997712, 9.294731722001391e-17, 0.005134763398167541, 1.194820309200107e-16, 0.004654430423785411, 8.2422739080512e-17, 0.004258029771397032, 9.5067306351894e-17, 0.0039283108282380024, 1.7540347128474968e-16, 0.0036523670873925395, 9.623794010508822e-17, 0.0034204021424253787, 1.4083470894369491e-16, 0.0032248884817922927, 1.4749333016985644e-16, 0.0030599828465501895, 1.0448590642474364e-16, 0.002921112944200383, 7.575487373767413e-18, 0.002804680975178716, 7.419510610361002e-17, 0.0027078483284668376, 3.924741709073613e-17, 0.0026283777262804953, 2.2684279102637236e-17, 0.0025645167846443107, 8.997077625295079e-17, 0.0025149120164513483, 7.131074184849219e-17, 0.0024785457043284276, 9.354417496250849e-17, 0.0024546904085875065, 1.2488589642405877e-16, 0.0024428775264784264], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485.5360118389, "effectiveFrequency": 553357.3374711702, "peak": 70.5, "rms": 51.572309672924284, "thd": 0.4833151484524849}, "harmonics": {"amplitudes": [24.2890625, 57.92076613061847, 1.421875, 19.27588907896988, 1.421875, 11.528257939603122, 1.421875, 8.194467538528329, 1.421875, 6.331896912839248, 1.421875, 5.137996046859012, 1.421875, 4.304077056139349, 1.421875, 3.6860723299088454, 1.421875, 3.207698601961777, 1.421875, 2.8247804703632298, 1.421875, 2.509960393415185, 1.421875, 2.2453859950684323, 1.421875, 2.01890737840567, 1.421875, 1.8219644341144872, 1.421875, 1.6483482744897402, 1.421875, 1.4934420157473332, 1.421875, 1.3537375367153817, 1.421875, 1.2265178099275544, 1.421875, 1.1096421410704556, 1.421875, 1.0013973584174929, 1.421875, 0.9003924136274832, 1.421875, 0.8054822382572133, 1.421875, 0.7157117294021269, 1.421875, 0.6302738400635857, 1.421875, 0.5484777114167545, 1.421875, 0.46972405216147894, 1.421875, 0.3934858059809043, 1.421875, 0.31929270856030145, 1.421875, 0.24671871675852053, 1.421875, 0.17537155450693565, 1.421875, 0.10488380107099537, 1.421875, 0.034905072061178544], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"ancillaryLabel": null, "data": [-5, 5, -5], "numberPeriods": null, "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746.40291779551, "effectiveFrequency": 110746.40291779551, "peak": 5, "rms": 2.8874560332150576, "thd": 0.12151487440704967}, "harmonics": {"amplitudes": [1.1608769501236793e-14, 4.05366124583194, 1.787369544444173e-15, 0.4511310569983995, 9.749053004706756e-16, 0.16293015292554872, 4.036157626725542e-16, 0.08352979924600704, 3.4998295008010614e-16, 0.0508569581336163, 3.1489164048780735e-16, 0.034320410449418075, 3.142469873118059e-16, 0.024811988673843106, 2.3653352035940994e-16, 0.018849001010678823, 2.9306524147249266e-16, 0.014866633059596499, 1.796485796132569e-16, 0.012077180559557785, 1.6247782523152451e-16, 0.010049063750920609, 1.5324769149805092e-16, 0.008529750975091871, 1.0558579733068502e-16, 0.007363501410705499, 7.513269775674661e-17, 0.006450045785294609, 5.871414177162291e-17, 0.005722473794997712, 9.294731722001391e-17, 0.005134763398167541, 1.194820309200107e-16, 0.004654430423785411, 8.2422739080512e-17, 0.004258029771397032, 9.5067306351894e-17, 0.0039283108282380024, 1.7540347128474968e-16, 0.0036523670873925395, 9.623794010508822e-17, 0.0034204021424253787, 1.4083470894369491e-16, 0.0032248884817922927, 1.4749333016985644e-16, 0.0030599828465501895, 1.0448590642474364e-16, 0.002921112944200383, 7.575487373767413e-18, 0.002804680975178716, 7.419510610361002e-17, 0.0027078483284668376, 3.924741709073613e-17, 0.0026283777262804953, 2.2684279102637236e-17, 0.0025645167846443107, 8.997077625295079e-17, 0.0025149120164513483, 7.131074184849219e-17, 0.0024785457043284276, 9.354417496250849e-17, 0.0024546904085875065, 1.2488589642405877e-16, 0.0024428775264784264], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-50, 50, 50, -50, -50], "numberPeriods": null, "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485.536011839, "effectiveFrequency": 591449.4202715514, "peak": 50, "rms": 50, "thd": 0.48331514845248497}, "harmonics": {"amplitudes": [0.78125, 63.64919355013018, 1.5625, 21.18229569117569, 1.5625, 12.668415318245188, 1.5625, 9.004909382998164, 1.5625, 6.958128475647527, 1.5625, 5.646149502042871, 1.5625, 4.729755006746538, 1.5625, 4.050628933965765, 1.5625, 3.524943518639316, 1.5625, 3.104154363036517, 1.5625, 2.7581982345221827, 1.5625, 2.467457137437843, 1.5625, 2.2185795367095267, 1.5625, 2.0021587188071255, 1.5625, 1.8113717302085082, 1.5625, 1.6411450722498175, 1.5625, 1.487623666720196, 1.5625, 1.3478217691511587, 1.5625, 1.2193869682092893, 1.5625, 1.100436657601639, 1.5625, 0.9894422127774558, 1.5625, 0.8851453167661671, 1.5625, 0.7864964059364037, 1.5625, 0.6926086154544899, 1.5625, 0.60272275979863, 1.5625, 0.5161802771005264, 1.5625, 0.43240198459440116, 1.5625, 0.3508711083080249, 1.5625, 0.27111946896540395, 1.5625, 0.192715993963664, 1.5625, 0.11525692425384548, 1.5625, 0.03835722204524927], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}]})"));

        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings.reset();
    }

    TEST_CASE("Test_Coil_Painter_Web_4", "[support][painter][magnetic-painter][round-winding-window]") {
        clear_databases();
        // settings.set_coil_delimit_and_compact(false);

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_4.svg");
        auto json_path_1493 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_painter_web_4_1493.json");
        std::ifstream json_file_1493(json_path_1493);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_1493));
        auto json_path_1494 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_painter_web_4_1494.json");
        std::ifstream json_file_1494(json_path_1494);
        OperatingPoint operatingPoint(json::parse(json_file_1494));

        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings.reset();
    }

    TEST_CASE("Test_Coil_Painter_Web_5", "[support][painter][magnetic-painter][round-winding-window]") {
        clear_databases();
        // settings.set_coil_delimit_and_compact(false);

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_5.svg");
        std::string file_path_1513 = std::source_location::current().file_name();
        auto json_path_1513 = file_path_1513.substr(0, file_path_1513.rfind("/")).append("/testData/test_coil_painter_web_5_1513.json");
        std::ifstream json_file_1513(json_path_1513);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_1513));

        Painter painter(outFile, false);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings.reset();
    }

    TEST_CASE("Test_Coil_Painter_Web_6", "[support][painter][magnetic-painter][round-winding-window]") {
        clear_databases();
        // settings.set_coil_delimit_and_compact(false);

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_6.svg");
        auto json_path_1530 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coil_painter_web_6_1530.json");
        std::ifstream json_file_1530(json_path_1530);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_1530));
        OperatingPoint operatingPoint(json::parse(R"({"conditions": {"ambientTemperature": 25}, "excitationsPerWinding": [{"frequency": 132000.0, "current": {"processed": {"label": "Triangular", "rms": 1.8989304812834227, "peakToPeak": 6.578088147248217, "offset": 0, "dutyCycle": 0.5}, "waveform": {"ancillaryLabel": null, "data": [-3.2890440736241087, 3.2890440736241087, -3.2890440736241087], "numberPeriods": null, "time": [0.0, 3.787878787878788e-06, 7.575757575757576e-06]}}, "voltage": {"harmonics": {"amplitudes": [62.952303569165444, 5128.785093504464, 125.90460713833089, 1706.847115062678, 125.90460713833089, 1020.8075863736791, 125.90460713833089, 725.6061300368981, 125.90460713833089, 560.6786765724385, 125.90460713833089, 454.9607903353537, 125.90460713833089, 381.1186854303848, 125.90460713833089, 326.39542054023394, 125.90460713833089, 284.0362424954169, 125.90460713833089, 250.12949476790286, 125.90460713833089, 222.25271368137805, 125.90460713833089, 198.82510177266022, 125.90460713833089, 178.77080638371407, 125.90460713833089, 161.331844428797, 125.90460713833089, 145.95842948696347, 125.90460713833089, 132.2417443703173, 125.90460713833089, 119.87115092997756, 125.90460713833089, 108.60606101597877, 125.90460713833089, 98.25691979647422, 125.90460713833089, 88.67202883580953, 125.90460713833089, 79.72821317492881, 125.90460713833089, 71.32407895537791, 125.90460713833089, 63.37505344328454, 125.90460713833089, 55.809674002829496, 125.90460713833089, 48.56676626289743, 125.90460713833089, 41.5932640005737, 125.90460713833089, 34.84249727756338, 125.90460713833089, 28.27282499053598, 125.90460713833089, 21.84652174569095, 125.90460713833089, 15.528852165931653, 125.90460713833089, 9.2872817716177, 125.90460713833089, 3.0907846224163222], "frequencies": [0.0, 132000.0, 264000.0, 396000.0, 528000.0, 660000.0, 792000.0, 924000.0, 1056000.0, 1188000.0, 1320000.0, 1452000.0, 1584000.0, 1716000.0, 1848000.0, 1980000.0, 2112000.0, 2244000.0, 2376000.0, 2508000.0, 2640000.0, 2772000.0, 2904000.0, 3036000.0, 3168000.0, 3300000.0, 3432000.0, 3564000.0, 3696000.0, 3828000.0, 3960000.0, 4092000.0, 4224000.0, 4356000.0, 4488000.0, 4620000.0, 4752000.0, 4884000.0, 5016000.0, 5148000.0, 5280000.0, 5412000.0, 5544000.0, 5676000.0, 5808000.0, 5940000.0, 6072000.0, 6204000.0, 6336000.0, 6468000.0, 6600000.0, 6732000.0, 6864000.0, 6996000.0, 7128000.0, 7260000.0, 7392000.0, 7524000.0, 7656000.0, 7788000.0, 7920000.0, 8052000.0, 8184000.0, 8316000.0]}, "processed": {"acEffectiveFrequency": 780760.9075356274, "average": 62.952303569165224, "dutyCycle": 0.51, "effectiveFrequency": 780713.2347584474, "label": "Custom", "offset": 0.0, "peak": 4028.9474284265884, "peakToPeak": 8057.894856853177, "phase": null, "rms": 4028.947428426593, "thd": 0.4833151484524849}, "waveform": {"ancillaryLabel": null, "data": [4028.9474284265884, 4028.9474284265884, -4028.9474284265884, -4028.9474284265884, 4028.947428426589], "numberPeriods": null, "time": [0.0, 3.787878787878788e-06, 3.787878787878788e-06, 7.575757575757576e-06, 7.575757575757576e-06]}}}, {"frequency": 132000.0, "current": {"processed": {"label": "Triangular", "rms": 1.9251336898395721, "peakToPeak": 6.668858724329366, "offset": 0, "dutyCycle": 0.5}, "waveform": {"ancillaryLabel": null, "data": [-3.334429362164683, 3.334429362164683, -3.334429362164683], "numberPeriods": null, "time": [0.0, 3.787878787878788e-06, 7.575757575757576e-06]}}, "voltage": {"harmonics": {"amplitudes": [63.82097799183204, 5199.556839373717, 127.64195598366408, 1730.399778717443, 127.64195598366408, 1034.893638677906, 127.64195598366408, 735.6187181449823, 127.64195598366408, 568.4154423150599, 127.64195598366408, 461.2387623788434, 127.64195598366408, 386.37771544618, 127.64195598366408, 330.8993280610653, 127.64195598366408, 287.9556386886795, 127.64195598366408, 253.58101412685173, 127.64195598366408, 225.31956329286425, 127.64195598366408, 201.56867541018784, 127.64195598366408, 181.23765220539858, 127.64195598366408, 163.55805123730477, 127.64195598366408, 147.97249962068966, 127.64195598366408, 134.0665389279477, 127.64195598366408, 121.52524453616411, 127.64195598366408, 110.10470843636273, 127.64195598366408, 99.61276014286315, 127.64195598366408, 89.89560794393546, 127.64195598366408, 80.82837719790015, 127.64195598366408, 72.30827491956077, 127.64195598366408, 64.24956136182033, 127.64195598366408, 56.579787780958924, 127.64195598366408, 49.23693566500447, 127.64195598366408, 42.167206533952594, 127.64195598366408, 35.32328645430246, 127.64195598366408, 28.66295972005929, 127.64195598366408, 22.147980367358002, 127.64195598366408, 15.743133708069479, 127.64195598366408, 9.415436321549691, 127.64195598366408, 3.1334341426909305], "frequencies": [0.0, 132000.0, 264000.0, 396000.0, 528000.0, 660000.0, 792000.0, 924000.0, 1056000.0, 1188000.0, 1320000.0, 1452000.0, 1584000.0, 1716000.0, 1848000.0, 1980000.0, 2112000.0, 2244000.0, 2376000.0, 2508000.0, 2640000.0, 2772000.0, 2904000.0, 3036000.0, 3168000.0, 3300000.0, 3432000.0, 3564000.0, 3696000.0, 3828000.0, 3960000.0, 4092000.0, 4224000.0, 4356000.0, 4488000.0, 4620000.0, 4752000.0, 4884000.0, 5016000.0, 5148000.0, 5280000.0, 5412000.0, 5544000.0, 5676000.0, 5808000.0, 5940000.0, 6072000.0, 6204000.0, 6336000.0, 6468000.0, 6600000.0, 6732000.0, 6864000.0, 6996000.0, 7128000.0, 7260000.0, 7392000.0, 7524000.0, 7656000.0, 7788000.0, 7920000.0, 8052000.0, 8184000.0, 8316000.0]}, "processed": {"acEffectiveFrequency": 780760.9075356275, "average": 63.820977991832244, "dutyCycle": 0.51, "effectiveFrequency": 780713.2347584475, "label": "Custom", "offset": 0.0, "peak": 4084.5425914772504, "peakToPeak": 8169.085182954501, "phase": null, "rms": 4084.5425914772554, "thd": 0.48331514845248497}, "waveform": {"ancillaryLabel": null, "data": [4084.5425914772504, 4084.5425914772504, -4084.5425914772504, -4084.5425914772504, 4084.542591477251], "numberPeriods": null, "time": [0.0, 3.787878787878788e-06, 3.787878787878788e-06, 7.575757575757576e-06, 7.575757575757576e-06]}}}, {"frequency": 132000.0, "current": {"processed": {"label": "Triangular", "rms": 0.358288770053476, "peakToPeak": 1.2411487070279656, "offset": 0, "dutyCycle": 0.5}, "waveform": {"ancillaryLabel": null, "data": [-0.6205743535139828, 0.6205743535139828, -0.6205743535139828], "numberPeriods": null, "time": [0.0, 3.787878787878788e-06, 7.575757575757576e-06]}}, "voltage": {"harmonics": {"amplitudes": [11.87779312625763, 967.6953006612196, 23.75558625251526, 322.04662548352417, 23.75558625251526, 192.60520497616582, 23.75558625251526, 136.906816988094, 23.75558625251526, 105.7884295419695, 23.75558625251526, 85.84165855384032, 23.75558625251526, 71.90918593026129, 23.75558625251526, 61.58404161136487, 23.75558625251526, 53.59174386705979, 23.75558625251526, 47.194244295830735, 23.75558625251526, 41.93447427950528, 23.75558625251526, 37.51417014578499, 23.75558625251526, 33.73034082711585, 23.75558625251526, 30.439970646942818, 23.75558625251526, 27.539326318295032, 23.75558625251526, 24.951272522701416, 23.75558625251526, 22.61719828867497, 23.75558625251526, 20.491709625656366, 23.75558625251526, 18.5390414710329, 23.75558625251526, 16.730571478454657, 23.75558625251526, 15.043059089609208, 23.75558625251526, 13.45737338780712, 23.75558625251526, 11.957557253449902, 23.75558625251526, 10.530127170345182, 23.75558625251526, 9.163540804320238, 23.75558625251526, 7.847785660485627, 23.75558625251526, 6.574056090106289, 23.75558625251526, 5.334495281233224, 23.75558625251526, 4.12198523503605, 23.75558625251526, 2.9299721067795446, 23.75558625251526, 1.7523173153996083, 23.75558625251526, 0.583166909889826], "frequencies": [0.0, 132000.0, 264000.0, 396000.0, 528000.0, 660000.0, 792000.0, 924000.0, 1056000.0, 1188000.0, 1320000.0, 1452000.0, 1584000.0, 1716000.0, 1848000.0, 1980000.0, 2112000.0, 2244000.0, 2376000.0, 2508000.0, 2640000.0, 2772000.0, 2904000.0, 3036000.0, 3168000.0, 3300000.0, 3432000.0, 3564000.0, 3696000.0, 3828000.0, 3960000.0, 4092000.0, 4224000.0, 4356000.0, 4488000.0, 4620000.0, 4752000.0, 4884000.0, 5016000.0, 5148000.0, 5280000.0, 5412000.0, 5544000.0, 5676000.0, 5808000.0, 5940000.0, 6072000.0, 6204000.0, 6336000.0, 6468000.0, 6600000.0, 6732000.0, 6864000.0, 6996000.0, 7128000.0, 7260000.0, 7392000.0, 7524000.0, 7656000.0, 7788000.0, 7920000.0, 8052000.0, 8184000.0, 8316000.0]}, "processed": {"acEffectiveFrequency": 780760.9075356274, "average": 11.877793126257593, "dutyCycle": 0.51, "effectiveFrequency": 780713.2347584474, "label": "Custom", "offset": 0.0, "peak": 760.1787600804884, "peakToPeak": 1520.3575201609767, "phase": null, "rms": 760.1787600804879, "thd": 0.48331514845248497}, "waveform": {"ancillaryLabel": null, "data": [760.1787600804884, 760.1787600804884, -760.1787600804884, -760.1787600804884, 760.1787600804885], "numberPeriods": null, "time": [0.0, 3.787878787878788e-06, 3.787878787878788e-06, 7.575757575757576e-06, 7.575757575757576e-06]}}}]})"));

        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings.reset();
    }

    TEST_CASE("Test_Coil_Painter_Web_7", "[support][painter][magnetic-painter][round-winding-window]") {
        clear_databases();

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "error_C_shape.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);
        auto magnetic = mas.get_magnetic();

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_7.svg");
        Painter painter(outFile, false);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings.reset();
    }

    TEST_CASE("Test_Coil_Painter_Web_8", "[support][painter][magnetic-painter][round-winding-window]") {
        clear_databases();

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "error_paint_field.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        auto outFile = outputFilePath;
        outFile.append("Test_Coil_Painter_Web_8.svg");
        Painter painter(outFile, false);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();

        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Distributed_Gap", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Distributed_Gap_Many", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Grinded_Gap", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_U_Core_Distributed_Gap", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_U_Core_Grinded_Gap", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Bobbin", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Section", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Bobbin_And_Section", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Bobbin_And_Two_Sections", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Epx_Core_Grinded_Gap", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Epx_Core_Spacer_Gap", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        // painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_P_Core_Grinded_Gap", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_U80_Core_Grinded_Gap", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Ep_Core_Grinded_Gap", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_All_Cores", "[support][painter][magnetic-painter][rectangular-winding-window]") {

        for (std::string& shapeName : get_core_shape_names()){
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
            Painter painter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);

            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_sections(magnetic);

            painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        }
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Grinded_Gap_Layers_No_Interleaving", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Grinded_Gap_Turns_No_Interleaving", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Top_Alignment", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Top_Alignment.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Bottom_Alignment", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Bottom_Alignment.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Spread_Alignment", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Spread_Alignment.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Vertical_Sections", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Vertical_Sections_Vectical_Layers", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {3, 3};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Vectical_Layers.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Vertical_Sections_Horizontal_Layers", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {1, 1};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Vertical_Sections_Horizontal_Layers_Spread_Turns", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Spread_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Vertical_Sections_Horizontal_Layers_Inner_Turns", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
        CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Inner_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Vertical_Sections_Horizontal_Layers_Outer_Turns", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
        CoilAlignment turnsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Outer_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Vertical_Sections_Horizontal_Layers_Centered_Turns", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {35, 35};
        std::vector<int64_t> numberParallels = {4, 4};
        uint8_t interleavingLevel = 3;
        int64_t numberStacks = 1;
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Centered_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Foil_Centered", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        OpenMagnetics::Wire wire;
        wire.set_nominal_value_outer_height(0.014);
        wire.set_nominal_value_outer_width(0.0002);
        wire.set_nominal_value_conducting_height(0.014);
        wire.set_nominal_value_conducting_width(0.0002);
        wire.set_type(WireType::FOIL);
        auto wires = std::vector<OpenMagnetics::Wire>({wire});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Foil_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Foil_Top", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

        OpenMagnetics::Wire wire;
        wire.set_nominal_value_outer_height(0.014);
        wire.set_nominal_value_outer_width(0.0002);
        wire.set_nominal_value_conducting_height(0.014);
        wire.set_nominal_value_conducting_width(0.0002);
        wire.set_type(WireType::FOIL);
        auto wires = std::vector<OpenMagnetics::Wire>({wire});


        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Foil_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Litz", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

        auto wire = find_wire_by_name("Litz TXXL07/28FXXX-2(MWXX)");
        auto wires = std::vector<OpenMagnetics::Wire>({wire});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Litz.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Foil_With_Insulation_Centered", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {4};
        std::vector<int64_t> numberParallels = {1};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        OpenMagnetics::Wire wire;
        wire.set_nominal_value_outer_height(0.014);
        wire.set_nominal_value_outer_width(0.0002);
        wire.set_nominal_value_conducting_height(0.0139);
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_type(WireType::FOIL);
        auto wires = std::vector<OpenMagnetics::Wire>({wire});

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Foil_With_Insulation_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Horizontal_Centered", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::INNER_OR_TOP;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Vertical_Centered", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Horizontal_Horizontal_Centered", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Horizontal_Centered.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Vertical_Top", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Top.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Horizontal_Inner", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Inner.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Horizontal_Outer", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Outer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Vertical_Bottom", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::OUTER_OR_BOTTOM;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Bottom.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Vertical_Spread", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Vertical_Spread_Two_Sections", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Spread_Two_Sections.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Vertical_Spread_One_Section", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23};
        std::vector<int64_t> numberParallels = {2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Spread_One_Section.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_Two_Sections", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23, 23};
        std::vector<int64_t> numberParallels = {2, 2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_Two_Sections.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_One_Section", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {23};
        std::vector<int64_t> numberParallels = {2};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_One_Section.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Bobbin_Vertical_Sections_And_Insulation", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);
        WindingOrientation sectionOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::CONTIGUOUS;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Vertical_Sections_And_Insulation.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Bobbin_Horizontal_Sections_And_Insulation", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Horizontal_Sections_And_Insulation.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_sections(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Bobbin_Layers_And_Insulation", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Layers_And_Insulation.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_layers(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Pq_Core_Bobbin_Turns_And_Insulation", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {
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
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.set_inputs(inputs);
        coil.wind();
        coil.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Turns_And_Insulation.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Turns_Not_Fitting", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {

        settings.set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {42, 42};
        std::vector<int64_t> numberParallels = {6, 6};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_outer_diameter(0.0015);
        wire.set_type(WireType::ROUND);
        auto wires = std::vector<OpenMagnetics::Wire>({wire, wire});

        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Turns_Not_Fitting.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);
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
            REQUIRE(false);
        }
        catch (const std::exception &e)
        {
            painter.export_svg();
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            REQUIRE(true);
            std::cerr << e.what() << std::endl;
        }
        settings.reset();
    }

    TEST_CASE("Test_Painter_Planar", "[support][painter][magnetic-painter][rectangular-winding-window][smoke-test]") {

        settings.set_coil_wind_even_if_not_fit(false);
        settings.set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {20, 5};
        std::vector<int64_t> numberParallels = {4, 4};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        double voltagePeakToPeak = 2000;
        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
        std::vector<size_t> stackUp = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        auto core = OpenMagneticsTesting::get_quick_core("ELP 38/8/25", json::parse("[]"), 1, "3C95");
        auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0008);
        wire.set_nominal_value_conducting_height(0.000076);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);
        wire.set_nominal_value_conducting_width(0.0032);
        wire.set_nominal_value_conducting_height(0.000076);
        wires.push_back(wire);

        OpenMagnetics::Coil coil;
        for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
            OpenMagnetics::Winding coilFunctionalDescription; 
            coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
            coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
            coilFunctionalDescription.set_name(OpenMagnetics::to_string(isolationSides[windingIndex]));
            coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
            coilFunctionalDescription.set_wire(wires[windingIndex]);
            coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
        }
        coil.set_bobbin(bobbin);
        coil.set_strict(false);

        coil.wind_by_planar_sections(stackUp, {{{0, 1}, 0.0001}}, 0.0001);
        coil.wind_by_planar_layers();
        coil.wind_by_planar_turns(0.0002, {{0, 0.0002}, {1, 0.0002}});
        coil.delimit_and_compact();

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(125000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        {
            OpenMagnetics::Mas mas;
            mas.set_inputs(inputs);
            mas.set_magnetic(magnetic);
            mas.set_outputs({});
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_Planar.mas.json");
            OpenMagnetics::to_file(outFile, mas);
        }
        {

            auto outFile = outputFilePath;
            outFile.append("Test_Painter_Planar.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            // settings.set_painter_number_points_x(300);
            // settings.set_painter_number_points_y(300);
            settings.set_painter_mode(PainterModes::SCATTER);
            // settings.set_painter_logarithmic_scale(false);
            settings.set_painter_include_fringing(false);
            // settings.set_painter_maximum_value_colorbar(std::nullopt);
            // settings.set_painter_minimum_value_colorbar(std::nullopt);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            // painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
            REQUIRE(std::filesystem::exists(outFile));
            settings.reset();
        }
    }

    TEST_CASE("Test_Field_Painter_Web_0", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        clear_databases();
        settings.set_painter_color_bobbin("0x7F539796");

        auto outFile = outputFilePath;
        outFile.append("Test_Field_Painter_Web_0.svg");
        auto json_path_3410 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_field_painter_web_0_3410.json");
        std::ifstream json_file_3410(json_path_3410);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_3410));
        OperatingPoint operatingPoint(json::parse(R"({"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42}, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"ancillaryLabel": null, "data": [-0.5, 0.5, -0.5], "numberPeriods": null, "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 1, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746.40291779555, "effectiveFrequency": 110746.40291779555, "peak": 0.5, "rms": 0.28874560332150573, "thd": 0.12151487440704967}, "harmonics": {"amplitudes": [1.1639994523804376e-15, 0.405366124583194, 1.7625133133414615e-16, 0.04511310569983995, 9.794112399879062e-17, 0.016293015292554884, 4.7988858344639954e-17, 0.008352979924600703, 4.095487581202406e-17, 0.005085695813361643, 3.44883176864717e-17, 0.003432041044941819, 3.006242608568973e-17, 0.002481198867384312, 2.7043620343491673e-17, 0.001884900101067885, 2.3307495124652547e-17, 0.001486663305959666, 2.2626256707089782e-17, 0.0012077180559557953, 1.7697356001223085e-17, 0.0010049063750920616, 2.036097837850594e-17, 0.0008529750975091948, 1.2438706838922496e-17, 0.0007363501410705534, 1.422708446941111e-17, 0.0006450045785294559, 1.022638379745176e-17, 0.0005722473794997664, 8.850209619122964e-18, 0.0005134763398167083, 1.2327526502681126e-17, 0.00046544304237854605, 6.116630655891947e-18, 0.0004258029771397049, 1.0990156230869393e-17, 0.0003928310828238053, 1.5090528294941328e-17, 0.00036523670873923176, 1.1048170322215263e-17, 0.0003420402142425299, 1.2778884268149239e-17, 0.00032248884817922554, 1.1578435604880218e-17, 0.000305998284655021, 8.597983809652515e-18, 0.0002921112944200549, 2.4395966302795302e-18, 0.00028046809751785943, 1.0171613695358528e-17, 0.00027078483284668445, 5.9415411850148865e-18, 0.0002628377726280502, 4.740823506837856e-18, 0.0002564516784644257, 1.0500990956496667e-17, 0.00025149120164513414, 7.790761372538887e-18, 0.00024785457043285074, 1.1991920621603288e-17, 0.00024546904085874163, 1.0014371763736106e-17, 0.0002442877526479259], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-500, 500, 500, -500, -500], "numberPeriods": null, "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 1000, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485.5360118389, "effectiveFrequency": 591449.4202715511, "peak": 500, "rms": 500, "thd": 0.4833151484524849}, "harmonics": {"amplitudes": [7.8125, 636.4919355013018, 15.625, 211.8229569117569, 15.625, 126.68415318245187, 15.625, 90.04909382998163, 15.625, 69.58128475647527, 15.625, 56.461495020428714, 15.625, 47.29755006746536, 15.625, 40.50628933965767, 15.625, 35.249435186393164, 15.625, 31.04154363036516, 15.625, 27.581982345221817, 15.625, 24.674571374378406, 15.625, 22.18579536709526, 15.625, 20.021587188071276, 15.625, 18.11371730208506, 15.625, 16.411450722498188, 15.625, 14.876236667201972, 15.625, 13.4782176915116, 15.625, 12.19386968209292, 15.625, 11.004366576016427, 15.625, 9.894422127774538, 15.625, 8.85145316766167, 15.625, 7.864964059364016, 15.625, 6.926086154544901, 15.625, 6.0272275979863, 15.625, 5.1618027710052665, 15.625, 4.324019845944026, 15.625, 3.50871108308025, 15.625, 2.7111946896540857, 15.625, 1.9271599396366668, 15.625, 1.1525692425384477, 15.625, 0.3835722204524359], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"ancillaryLabel": null, "data": [-5, 5, -5], "numberPeriods": null, "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746.40291779551, "effectiveFrequency": 110746.40291779551, "peak": 5, "rms": 2.8874560332150576, "thd": 0.12151487440704967}, "harmonics": {"amplitudes": [1.1608769501236793e-14, 4.05366124583194, 1.787369544444173e-15, 0.4511310569983995, 9.749053004706756e-16, 0.16293015292554872, 4.036157626725542e-16, 0.08352979924600704, 3.4998295008010614e-16, 0.0508569581336163, 3.1489164048780735e-16, 0.034320410449418075, 3.142469873118059e-16, 0.024811988673843106, 2.3653352035940994e-16, 0.018849001010678823, 2.9306524147249266e-16, 0.014866633059596499, 1.796485796132569e-16, 0.012077180559557785, 1.6247782523152451e-16, 0.010049063750920609, 1.5324769149805092e-16, 0.008529750975091871, 1.0558579733068502e-16, 0.007363501410705499, 7.513269775674661e-17, 0.006450045785294609, 5.871414177162291e-17, 0.005722473794997712, 9.294731722001391e-17, 0.005134763398167541, 1.194820309200107e-16, 0.004654430423785411, 8.2422739080512e-17, 0.004258029771397032, 9.5067306351894e-17, 0.0039283108282380024, 1.7540347128474968e-16, 0.0036523670873925395, 9.623794010508822e-17, 0.0034204021424253787, 1.4083470894369491e-16, 0.0032248884817922927, 1.4749333016985644e-16, 0.0030599828465501895, 1.0448590642474364e-16, 0.002921112944200383, 7.575487373767413e-18, 0.002804680975178716, 7.419510610361002e-17, 0.0027078483284668376, 3.924741709073613e-17, 0.0026283777262804953, 2.2684279102637236e-17, 0.0025645167846443107, 8.997077625295079e-17, 0.0025149120164513483, 7.131074184849219e-17, 0.0024785457043284276, 9.354417496250849e-17, 0.0024546904085875065, 1.2488589642405877e-16, 0.0024428775264784264], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-50, 50, 50, -50, -50], "numberPeriods": null, "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485.536011839, "effectiveFrequency": 591449.4202715514, "peak": 50, "rms": 50, "thd": 0.48331514845248497}, "harmonics": {"amplitudes": [0.78125, 63.64919355013018, 1.5625, 21.18229569117569, 1.5625, 12.668415318245188, 1.5625, 9.004909382998164, 1.5625, 6.958128475647527, 1.5625, 5.646149502042871, 1.5625, 4.729755006746538, 1.5625, 4.050628933965765, 1.5625, 3.524943518639316, 1.5625, 3.104154363036517, 1.5625, 2.7581982345221827, 1.5625, 2.467457137437843, 1.5625, 2.2185795367095267, 1.5625, 2.0021587188071255, 1.5625, 1.8113717302085082, 1.5625, 1.6411450722498175, 1.5625, 1.487623666720196, 1.5625, 1.3478217691511587, 1.5625, 1.2193869682092893, 1.5625, 1.100436657601639, 1.5625, 0.9894422127774558, 1.5625, 0.8851453167661671, 1.5625, 0.7864964059364037, 1.5625, 0.6926086154544899, 1.5625, 0.60272275979863, 1.5625, 0.5161802771005264, 1.5625, 0.43240198459440116, 1.5625, 0.3508711083080249, 1.5625, 0.27111946896540395, 1.5625, 0.192715993963664, 1.5625, 0.11525692425384548, 1.5625, 0.03835722204524927], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}]})"));

        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Field_Painter_Web_1", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        clear_databases();

        auto outFile = outputFilePath;
        outFile.append("Test_Field_Painter_Web_1.svg");
        auto json_path_3429 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_field_painter_web_1_3429.json");
        std::ifstream json_file_3429(json_path_3429);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_3429));
        OperatingPoint operatingPoint(json::parse(R"({"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42}, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"ancillaryLabel": null, "data": [-0.5, 0.5, -0.5], "numberPeriods": null, "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 1, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746.40291779555, "effectiveFrequency": 110746.40291779555, "peak": 0.5, "rms": 0.28874560332150573, "thd": 0.12151487440704967}, "harmonics": {"amplitudes": [1.1639994523804376e-15, 0.405366124583194, 1.7625133133414615e-16, 0.04511310569983995, 9.794112399879062e-17, 0.016293015292554884, 4.7988858344639954e-17, 0.008352979924600703, 4.095487581202406e-17, 0.005085695813361643, 3.44883176864717e-17, 0.003432041044941819, 3.006242608568973e-17, 0.002481198867384312, 2.7043620343491673e-17, 0.001884900101067885, 2.3307495124652547e-17, 0.001486663305959666, 2.2626256707089782e-17, 0.0012077180559557953, 1.7697356001223085e-17, 0.0010049063750920616, 2.036097837850594e-17, 0.0008529750975091948, 1.2438706838922496e-17, 0.0007363501410705534, 1.422708446941111e-17, 0.0006450045785294559, 1.022638379745176e-17, 0.0005722473794997664, 8.850209619122964e-18, 0.0005134763398167083, 1.2327526502681126e-17, 0.00046544304237854605, 6.116630655891947e-18, 0.0004258029771397049, 1.0990156230869393e-17, 0.0003928310828238053, 1.5090528294941328e-17, 0.00036523670873923176, 1.1048170322215263e-17, 0.0003420402142425299, 1.2778884268149239e-17, 0.00032248884817922554, 1.1578435604880218e-17, 0.000305998284655021, 8.597983809652515e-18, 0.0002921112944200549, 2.4395966302795302e-18, 0.00028046809751785943, 1.0171613695358528e-17, 0.00027078483284668445, 5.9415411850148865e-18, 0.0002628377726280502, 4.740823506837856e-18, 0.0002564516784644257, 1.0500990956496667e-17, 0.00025149120164513414, 7.790761372538887e-18, 0.00024785457043285074, 1.1991920621603288e-17, 0.00024546904085874163, 1.0014371763736106e-17, 0.0002442877526479259], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-500, 500, 500, -500, -500], "numberPeriods": null, "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 1000, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485.5360118389, "effectiveFrequency": 591449.4202715511, "peak": 500, "rms": 500, "thd": 0.4833151484524849}, "harmonics": {"amplitudes": [7.8125, 636.4919355013018, 15.625, 211.8229569117569, 15.625, 126.68415318245187, 15.625, 90.04909382998163, 15.625, 69.58128475647527, 15.625, 56.461495020428714, 15.625, 47.29755006746536, 15.625, 40.50628933965767, 15.625, 35.249435186393164, 15.625, 31.04154363036516, 15.625, 27.581982345221817, 15.625, 24.674571374378406, 15.625, 22.18579536709526, 15.625, 20.021587188071276, 15.625, 18.11371730208506, 15.625, 16.411450722498188, 15.625, 14.876236667201972, 15.625, 13.4782176915116, 15.625, 12.19386968209292, 15.625, 11.004366576016427, 15.625, 9.894422127774538, 15.625, 8.85145316766167, 15.625, 7.864964059364016, 15.625, 6.926086154544901, 15.625, 6.0272275979863, 15.625, 5.1618027710052665, 15.625, 4.324019845944026, 15.625, 3.50871108308025, 15.625, 2.7111946896540857, 15.625, 1.9271599396366668, 15.625, 1.1525692425384477, 15.625, 0.3835722204524359], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"ancillaryLabel": null, "data": [-5, 5, -5], "numberPeriods": null, "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746.40291779551, "effectiveFrequency": 110746.40291779551, "peak": 5, "rms": 2.8874560332150576, "thd": 0.12151487440704967}, "harmonics": {"amplitudes": [1.1608769501236793e-14, 4.05366124583194, 1.787369544444173e-15, 0.4511310569983995, 9.749053004706756e-16, 0.16293015292554872, 4.036157626725542e-16, 0.08352979924600704, 3.4998295008010614e-16, 0.0508569581336163, 3.1489164048780735e-16, 0.034320410449418075, 3.142469873118059e-16, 0.024811988673843106, 2.3653352035940994e-16, 0.018849001010678823, 2.9306524147249266e-16, 0.014866633059596499, 1.796485796132569e-16, 0.012077180559557785, 1.6247782523152451e-16, 0.010049063750920609, 1.5324769149805092e-16, 0.008529750975091871, 1.0558579733068502e-16, 0.007363501410705499, 7.513269775674661e-17, 0.006450045785294609, 5.871414177162291e-17, 0.005722473794997712, 9.294731722001391e-17, 0.005134763398167541, 1.194820309200107e-16, 0.004654430423785411, 8.2422739080512e-17, 0.004258029771397032, 9.5067306351894e-17, 0.0039283108282380024, 1.7540347128474968e-16, 0.0036523670873925395, 9.623794010508822e-17, 0.0034204021424253787, 1.4083470894369491e-16, 0.0032248884817922927, 1.4749333016985644e-16, 0.0030599828465501895, 1.0448590642474364e-16, 0.002921112944200383, 7.575487373767413e-18, 0.002804680975178716, 7.419510610361002e-17, 0.0027078483284668376, 3.924741709073613e-17, 0.0026283777262804953, 2.2684279102637236e-17, 0.0025645167846443107, 8.997077625295079e-17, 0.0025149120164513483, 7.131074184849219e-17, 0.0024785457043284276, 9.354417496250849e-17, 0.0024546904085875065, 1.2488589642405877e-16, 0.0024428775264784264], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-50, 50, 50, -50, -50], "numberPeriods": null, "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485.536011839, "effectiveFrequency": 591449.4202715514, "peak": 50, "rms": 50, "thd": 0.48331514845248497}, "harmonics": {"amplitudes": [0.78125, 63.64919355013018, 1.5625, 21.18229569117569, 1.5625, 12.668415318245188, 1.5625, 9.004909382998164, 1.5625, 6.958128475647527, 1.5625, 5.646149502042871, 1.5625, 4.729755006746538, 1.5625, 4.050628933965765, 1.5625, 3.524943518639316, 1.5625, 3.104154363036517, 1.5625, 2.7581982345221827, 1.5625, 2.467457137437843, 1.5625, 2.2185795367095267, 1.5625, 2.0021587188071255, 1.5625, 1.8113717302085082, 1.5625, 1.6411450722498175, 1.5625, 1.487623666720196, 1.5625, 1.3478217691511587, 1.5625, 1.2193869682092893, 1.5625, 1.100436657601639, 1.5625, 0.9894422127774558, 1.5625, 0.8851453167661671, 1.5625, 0.7864964059364037, 1.5625, 0.6926086154544899, 1.5625, 0.60272275979863, 1.5625, 0.5161802771005264, 1.5625, 0.43240198459440116, 1.5625, 0.3508711083080249, 1.5625, 0.27111946896540395, 1.5625, 0.192715993963664, 1.5625, 0.11525692425384548, 1.5625, 0.03835722204524927], "frequencies": [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000, 3600000, 3700000, 3800000, 3900000, 4000000, 4100000, 4200000, 4300000, 4400000, 4500000, 4600000, 4700000, 4800000, 4900000, 5000000, 5100000, 5200000, 5300000, 5400000, 5500000, 5600000, 5700000, 5800000, 5900000, 6000000, 6100000, 6200000, 6300000]}}}]})"));

        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Field_Painter_Web_2", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        clear_databases();

        auto outFile = outputFilePath;
        outFile.append("Test_Field_Painter_Web_2.svg");

        OpenMagnetics::Magnetic magnetic(json::parse(R"({"coil": {"bobbin": {"distributorsInfo": null, "functionalDescription": null, "manufacturerInfo": null, "name": null, "processedDescription": {"columnDepth": 0.0015375000000000002, "columnShape": "round", "columnThickness": 0.0006750000000000003, "columnWidth": 0.0015375000000000002, "coordinates": [0, 0, 0], "pins": null, "wallThickness": 0.0006328382137062499, "windingWindows": [{"angle": null, "area": 2.115668411496562e-06, "coordinates": [0.0018750000000000001, 0, 0], "height": 0.0031343235725874996, "radialHeight": null, "sectionsAlignment": "inner or top", "sectionsOrientation": "overlapping", "shape": "rectangular", "width": 0.0006749999999999998}]}}, "functionalDescription": [{"connections": null, "isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 10, "wire": {"coating": {"breakdownVoltage": 4500, "grade": null, "material": "ETFE", "numberLayers": 3, "temperatureRating": 155, "thickness": null, "thicknessLayers": 2.54e-05, "type": "insulated"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.00224}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": null, "material": "copper", "name": "Round 2.24 - Grade 1", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.002316}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": "2.24 mm", "strand": null, "type": "round"}}], "layersDescription": [{"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.0026955000000000017, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 0", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.005011500000000002, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 1", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.0073275000000000015, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 2", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.009643500000000001, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 3", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.011959500000000001, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 4", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.014275500000000002, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 5", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.016591500000000002, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 6", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.0189075, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 7", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.0212235, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 8", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}, {"additionalCoordinates": null, "coordinateSystem": "cartesian", "coordinates": [0.023539500000000005, 0], "dimensions": [0.002316, 0.002316], "fillingFactor": 0.738915413920732, "insulationMaterial": null, "name": "Primary section 0 layer 9", "orientation": "overlapping", "partialWindings": [{"connections": null, "parallelsProportion": [0.1], "winding": "Primary"}], "section": "Primary section 0", "turnsAlignment": "spread", "type": "conduction", "windingStyle": "windByConsecutiveTurns"}], "sectionsDescription": [{"coordinateSystem": "cartesian", "coordinates": [0.0131175, 0], "dimensions": [0.02316, 0.002316], "fillingFactor": 25.35300886874689, "layersAlignment": null, "layersOrientation": "overlapping", "margin": [0, 0], "name": "Primary section 0", "partialWindings": [{"connections": null, "parallelsProportion": [1], "winding": "Primary"}], "type": "conduction", "windingStyle": "windByConsecutiveTurns"}], "turnsDescription": [{"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0026955000000000017, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 0", "length": 0.016936325995502585, "name": "Primary parallel 0 turn 0", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.005011500000000002, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 1", "length": 0.03148818316693051, "name": "Primary parallel 0 turn 1", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0073275000000000015, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 2", "length": 0.04604004033835843, "name": "Primary parallel 0 turn 2", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.009643500000000001, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 3", "length": 0.060591897509786344, "name": "Primary parallel 0 turn 3", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.011959500000000001, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 4", "length": 0.07514375468121427, "name": "Primary parallel 0 turn 4", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.014275500000000002, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 5", "length": 0.08969561185264219, "name": "Primary parallel 0 turn 5", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.016591500000000002, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 6", "length": 0.10424746902407012, "name": "Primary parallel 0 turn 6", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0189075, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 7", "length": 0.11879932619549803, "name": "Primary parallel 0 turn 7", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.0212235, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 8", "length": 0.13335118336692595, "name": "Primary parallel 0 turn 8", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}, {"additionalCoordinates": null, "angle": null, "coordinateSystem": "cartesian", "coordinates": [0.023539500000000005, 0], "dimensions": [0.002316, 0.002316], "layer": "Primary section 0 layer 9", "length": 0.1479030405383539, "name": "Primary parallel 0 turn 9", "orientation": "clockwise", "parallel": 0, "rotation": 0, "section": "Primary section 0", "winding": "Primary"}]}, "core": {"distributorsInfo": null, "functionalDescription": {"coating": null, "gapping": [{"area": 3e-06, "coordinates": [0, 0, 0], "distanceClosestNormalSurface": 0.0021975000000000002, "distanceClosestParallelSurface": 0.00135, "length": 5e-06, "sectionDimensions": [0.001725, 0.001725], "shape": "round", "type": "residual"}, {"area": 1.9e-05, "coordinates": [0, 0, -0.002545], "distanceClosestNormalSurface": 0.0021975000000000002, "distanceClosestParallelSurface": 0.00135, "length": 5e-06, "sectionDimensions": [0.028658, 0.000663], "shape": "irregular", "type": "residual"}], "material": "44", "numberStacks": 1, "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00615, "minimum": 0.00585, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00305, "minimum": 0.00295, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0039, "minimum": 0.00365, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0023, "minimum": 0.0021, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00455, "minimum": 0.0043, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0018, "minimum": 0.00165, "nominal": null}, "K": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0009, "minimum": null, "nominal": null}}, "family": "ep", "familySubtype": null, "magneticCircuit": "open", "name": "EP 6", "type": "standard"}, "type": "two-piece set"}, "geometricalDescription": [{"coordinates": [0, 0, 0], "dimensions": null, "insulationMaterial": null, "machining": null, "material": "44", "rotation": [3.141592653589793, 3.141592653589793, 0], "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00615, "minimum": 0.00585, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00305, "minimum": 0.00295, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0039, "minimum": 0.00365, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0023, "minimum": 0.0021, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00455, "minimum": 0.0043, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0018, "minimum": 0.00165, "nominal": null}, "K": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0009, "minimum": null, "nominal": null}}, "family": "ep", "familySubtype": null, "magneticCircuit": "open", "name": "EP 6", "type": "standard"}, "type": "half set"}, {"coordinates": [0, 0, 0], "dimensions": null, "insulationMaterial": null, "machining": null, "material": "44", "rotation": [0, 0, 0], "shape": {"aliases": [], "dimensions": {"A": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00615, "minimum": 0.00585, "nominal": null}, "B": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00305, "minimum": 0.00295, "nominal": null}, "C": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0039, "minimum": 0.00365, "nominal": null}, "D": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0023, "minimum": 0.0021, "nominal": null}, "E": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00455, "minimum": 0.0043, "nominal": null}, "F": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0018, "minimum": 0.00165, "nominal": null}, "K": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0009, "minimum": null, "nominal": null}}, "family": "ep", "familySubtype": null, "magneticCircuit": "open", "name": "EP 6", "type": "standard"}, "type": "half set"}], "manufacturerInfo": null, "name": "Custom", "processedDescription": {"columns": [{"area": 3e-06, "coordinates": [0, 0, 0], "depth": 0.001725, "height": 0.0044, "minimumDepth": null, "minimumWidth": null, "shape": "round", "type": "central", "width": 0.001725}, {"area": 1.9e-05, "coordinates": [0, 0, -0.002545], "depth": 0.000663, "height": 0.0044, "minimumDepth": null, "minimumWidth": 0.000788, "shape": "irregular", "type": "lateral", "width": 0.028658}], "depth": 0.0037749999999999997, "effectiveParameters": {"effectiveArea": 3.148296314833282e-06, "effectiveLength": 0.010354750378339139, "effectiveVolume": 3.2599822457143645e-08, "minimumArea": 2.401737909959532e-06}, "height": 0.006, "width": 0.006, "windingWindows": [{"angle": null, "area": 5.94e-06, "coordinates": [0.0008625, 0], "height": 0.004399999999999999, "radialHeight": null, "sectionsAlignment": null, "sectionsOrientation": null, "shape": null, "width": 0.00135}]}}, "manufacturerInfo": {"name": "OpenMagnetics", "reference": "My custom magnetic"}})"));
        OperatingPoint operatingPoint(json::parse(R"({"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42}, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}]})"));

        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Field_Painter_Web_3", "[support][painter][magnetic-field-painter][rectangular-winding-window]") {
        clear_databases();

        auto outFile = outputFilePath;
        outFile.append("Test_Field_Painter_Web_3.svg");

        auto json_path_3469 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_field_painter_web_3_3469.json");
        std::ifstream json_file_3469(json_path_3469);
        OpenMagnetics::Magnetic magnetic(json::parse(json_file_3469));

        std::filesystem::remove(outFile);
        Painter painter(outFile, true);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Round_Enamelled_Grade_1", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Round 0.335 - Grade 1");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_Enamelled_Grade_1.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Round_Enamelled_Grade_2", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Round 0.335 - Grade 2");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_Enamelled_Grade_2.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Round_Enamelled_Grade_9", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Round 0.071 - FIW 9");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_Enamelled_Grade_9.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Round_TIW_PFA", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Round T17A01PXXX-1.5");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_TIW_PFA.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Round_DIW_ETFE", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Round D32A01TXX-3");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_DIW_ETFE.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Round_SIW_PEF", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Round S24A01FX-2");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_SIW_PEF.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Round_SIW_TCA", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Round TCA1 31 AWG");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Round_SIW_TCA.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Rectangular_Enamelled_Grade_1", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Rectangular 2x0.80 - Grade 1");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Rectangular_Enamelled_Grade_1.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Rectangular_Enamelled_Grade_2", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Rectangular 2x0.80 - Grade 2");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Rectangular_Enamelled_Grade_2.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Foil_Bare", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Foil 0.7");
        DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithTolerance);
        wire.set_outer_width(wire.get_conducting_width().value());
        wire.set_outer_height(dimensionWithTolerance);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Foil_Bare.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Litz_Simple_TIW", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Litz TXXL07/28FXXX-2(MWXX)");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_Insulation.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Litz_Simple_Bare", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Litz 75x0.3 - Grade 1 - Unserved");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_Bare.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Litz_Simple_Single_Served", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();

        auto wire = find_wire_by_name("Litz 75x0.3 - Grade 1 - Single Served");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_Single_Served.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Litz_Simple_Double_Served", "[support][painter][wire-painter][smoke-test]") {
        clear_databases(); 
        settings.set_painter_simple_litz(false);
        settings.set_painter_advanced_litz(true);

        auto wire = find_wire_by_name("Litz 270x0.02 - Grade 1 - Single Served");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_Double_Served.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Litz_Simple_Insulated", "[support][painter][wire-painter][smoke-test]") {
        clear_databases(); 
        settings.set_painter_simple_litz(false);
        settings.set_painter_advanced_litz(true);

        auto wire = find_wire_by_name("Litz DXXL550/44TXX-3(MWXX)");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Simple_Insulated.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Litz_Few_Strands", "[support][painter][wire-painter][smoke-test]") {
        clear_databases(); 
        settings.set_painter_simple_litz(false);
        settings.set_painter_advanced_litz(true);

        auto wire = find_wire_by_name("Litz 4x0.071 - Grade 2 - Double Served");

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Few_Strands.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Web_0", "[support][painter][wire-painter][smoke-test]") {
        clear_databases(); 
        settings.set_painter_simple_litz(false);
        settings.set_painter_advanced_litz(true);

        OpenMagnetics::Wire wire(json::parse(R"({"standard": "IEC 60317", "type": "litz", "strand": "Round 0.2 - Grade 1", "numberConductors": 17, "coating": {"breakdownVoltage": null, "grade": null, "material": null, "numberLayers": 1, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "served"}})"));

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Litz_Few_Strands.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Rectangular_Wire_Tiny", "[support][painter][wire-painter][smoke-test]") {
        clear_databases();


        {
            OpenMagnetics::Wire wire;


            DimensionWithTolerance dimensionWithTolerance;
            wire.set_type(WireType::RECTANGULAR);
            wire.set_nominal_value_conducting_height(2e-6);
            wire.set_nominal_value_conducting_width(4e-6);
            wire.set_nominal_value_outer_height(OpenMagnetics::Wire::get_outer_height_rectangular(2e-6, 2, WireStandard::IEC_60317));
            wire.set_nominal_value_outer_width(OpenMagnetics::Wire::get_outer_height_rectangular(4e-6, 2, WireStandard::IEC_60317));
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Rectangular_Wire_Tiny.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }

        {
            OpenMagnetics::Wire wire;


            DimensionWithTolerance dimensionWithTolerance;
            wire.set_type(WireType::RECTANGULAR);
            wire.set_nominal_value_conducting_height(2e-6);
            wire.set_nominal_value_conducting_width(5e-6);
            wire.set_nominal_value_outer_height(OpenMagnetics::Wire::get_outer_height_rectangular(2e-6, 2, WireStandard::IEC_60317));
            wire.set_nominal_value_outer_width(OpenMagnetics::Wire::get_outer_height_rectangular(5e-6, 2, WireStandard::IEC_60317));
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Rectangular_Wire_Tiny2.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Round_Enamelled_Grade_1", "[support][painter][wire-density-painter][smoke-test]") {
        clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000;
        auto wire = find_wire_by_name("Round 3.55 - Grade 1");
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency, 0.001, 25, WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Round_Enamelled_Grade_1_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Litz_Single_Served_Simple", "[support][painter][wire-density-painter][smoke-test]") {
        clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000000;
        auto wire = find_wire_by_name("Litz 10x0.02 - Grade 2 - Single Served");
        settings.set_painter_simple_litz(true);
        // settings.set_painter_advanced_litz(false);

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency, 0.001, 25, WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Litz_Single_Served_Simple_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Litz_Single_Served_Normal", "[support][painter][wire-density-painter][smoke-test]") {
        clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 100e6;
        auto wire = find_wire_by_name("Litz 8x0.1 - Grade 1 - Single Served");
        settings.set_painter_simple_litz(false);
        // settings.set_painter_advanced_litz(false);

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency, 0.001, 25, WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Litz_Single_Served_Normal_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Litz_Single_Served_Normal_Many_Strands", "[support][painter][wire-density-painter]") {
        clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000000;
        auto wire = find_wire_by_name("Litz 1000x0.02 - Grade 1 - Single Served");
        settings.set_painter_simple_litz(false);
        settings.set_painter_advanced_litz(true);

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency, 0.001, 25, WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Litz_Single_Served_Normal_Many_Strands_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Litz_Single_Served_Advanced", "[support][painter][wire-density-painter][smoke-test]") {
        clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000000;
        auto wire = find_wire_by_name("Litz 10x0.02 - Grade 2 - Single Served");
        settings.set_painter_simple_litz(false);
        settings.set_painter_advanced_litz(true);

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency, 0.001, 25, WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Litz_Single_Served_Advanced_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Rectangular_Enamelled_Grade_1", "[support][painter][wire-density-painter][smoke-test]") {
        clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000;
        auto wire = find_wire_by_name("Rectangular 2x0.80 - Grade 1");
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency, 0.001, 25, WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Rectangular_Enamelled_Grade_1_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Foil", "[support][painter][wire-density-painter][smoke-test]") {
        clear_databases();

        double currentPeakToPeak = 1000;
        double frequency = 250000;
        auto wire = find_wire_by_name("Foil 0.5");
        DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(0.010);
        wire.set_conducting_height(dimensionWithTolerance);
        wire.set_outer_width(wire.get_conducting_width().value());
        wire.set_outer_height(dimensionWithTolerance);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency, 0.001, 25, WaveformLabel::TRIANGULAR, currentPeakToPeak, 0.5, 0);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Foil_" + std::to_string(frequency) + ".svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);
            painter.paint_wire_with_current_density(wire, inputs.get_operating_point(0));
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Web_0", "[support][painter][wire-density-painter][smoke-test]") {
        clear_databases();

        OpenMagnetics::Wire wire(json::parse(R"({"coating": {"breakdownVoltage": 6000, "grade": null, "material": "TCA", "numberLayers": 3, "temperatureRating": 155, "thickness": null, "thicknessLayers": 7.62e-05, "type": "insulated"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00026900000000000003, "minimum": 0.000261, "nominal": 0.000265}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.265 - Grade 1", "numberConductors": 92, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000297, "minimum": 0.000283, "nominal": 0.0018383994787140629}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": "0.265 mm", "strand": {"coating": {"breakdownVoltage": 500, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00010300000000000001, "minimum": 9.7e-05, "nominal": 0.0001}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.1 - Grade 1", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00011700000000000001, "minimum": 0.000108, "nominal": null}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": "0.1 mm", "strand": null, "type": "round"}, "type": "litz"})"));
        OperatingPoint operatingPoint(json::parse(R"({"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42}, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 100000000, "current": {"waveform": {"ancillaryLabel": null, "data": [-5, 5, -5], "numberPeriods": null, "time": [0, 5e-09, 1e-08]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular", "acEffectiveFrequency": 110746402.91779621, "effectiveFrequency": 110746402.91779621, "peak": 5, "rms": 2.887456033215047, "thd": 0.12151487440705}, "harmonics": {"amplitudes": [1.9671764217576992e-15, 4.053661245831925, 9.528877815088396e-15, 0.4511310569983986, 3.867859480230825e-15, 0.16293015292554885, 2.6154487633755876e-15, 0.08352979924600734, 1.9651652434930284e-15, 0.05085695813361692, 1.6136051601223988e-15, 0.034320410449418734, 1.3808009023999193e-15, 0.024811988673843956, 1.1386794830042506e-15, 0.018849001010679767, 1.0370440806127272e-15, 0.014866633059597558, 8.851494177631955e-16, 0.012077180559558507, 9.326505898457337e-16, 0.010049063750921247, 8.13280272912139e-16, 0.008529750975092641, 7.240094812312869e-16, 0.007363501410706191, 6.482105482773775e-16, 0.006450045785295375, 6.066430502657454e-16, 0.0057224737949983576, 5.440719337350025e-16, 0.005134763398167652, 6.25059167873935e-16, 0.004654430423786248, 4.892419283154727e-16, 0.004258029771397705, 5.345312821060954e-16, 0.0039283108282387675, 4.672784117529907e-16, 0.003652367087393178, 3.3833965320162966e-16, 0.003420402142426004, 4.2496275203954786e-16, 0.0032248884817931253, 3.8603919165942636e-16, 0.0030599828465507064, 3.809779240559984e-16, 0.0029211129442011045, 4.459588218866508e-16, 0.002804680975179318, 3.8344856196395527e-16, 0.0027078483284674916, 3.809203930919387e-16, 0.0026283777262812655, 4.0712748931445545e-16, 0.0025645167846449803, 3.9966373222574786e-16, 0.0025149120164521047, 4.572347145426278e-16, 0.0024785457043292464, 3.6473608081573336e-16, 0.0024546904085880616, 3.509577233239641e-16, 0.0024428775264793146], "frequencies": [0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000, 900000000, 1000000000, 1100000000, 1200000000, 1300000000, 1400000000, 1500000000, 1600000000, 1700000000, 1800000000, 1900000000, 2000000000, 2100000000, 2200000000, 2300000000, 2400000000, 2500000000, 2600000000, 2700000000, 2800000000, 2900000000, 3000000000, 3100000000, 3200000000, 3300000000, 3400000000, 3500000000, 3600000000, 3700000000, 3800000000, 3900000000, 4000000000, 4100000000, 4200000000, 4300000000, 4400000000, 4500000000, 4600000000, 4700000000, 4800000000, 4900000000, 5000000000, 5100000000, 5200000000, 5300000000, 5400000000, 5500000000, 5600000000, 5700000000, 5800000000, 5900000000, 6000000000, 6100000000, 6200000000, 6300000000]}}, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-50, 50, 50, -50, -50], "numberPeriods": null, "time": [0, 0, 5e-09, 5e-09, 1e-08]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 591485536.0118392, "effectiveFrequency": 591449420.2715513, "peak": 50, "rms": 50, "thd": 0.48331514845248497}, "harmonics": {"amplitudes": [0.78125, 63.64919355013018, 1.5625, 21.18229569117569, 1.5625, 12.668415318245188, 1.5625, 9.004909382998164, 1.5625, 6.958128475647527, 1.5625, 5.646149502042871, 1.5625, 4.729755006746538, 1.5625, 4.050628933965765, 1.5625, 3.524943518639316, 1.5625, 3.104154363036517, 1.5625, 2.7581982345221827, 1.5625, 2.467457137437843, 1.5625, 2.2185795367095267, 1.5625, 2.0021587188071255, 1.5625, 1.8113717302085082, 1.5625, 1.6411450722498175, 1.5625, 1.487623666720196, 1.5625, 1.3478217691511587, 1.5625, 1.2193869682092893, 1.5625, 1.100436657601639, 1.5625, 0.9894422127774558, 1.5625, 0.8851453167661671, 1.5625, 0.7864964059364037, 1.5625, 0.6926086154544899, 1.5625, 0.60272275979863, 1.5625, 0.5161802771005264, 1.5625, 0.43240198459440116, 1.5625, 0.3508711083080249, 1.5625, 0.27111946896540395, 1.5625, 0.192715993963664, 1.5625, 0.11525692425384548, 1.5625, 0.03835722204524927], "frequencies": [0, 100000000, 200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000, 900000000, 1000000000, 1100000000, 1200000000, 1300000000, 1400000000, 1500000000, 1600000000, 1700000000, 1800000000, 1900000000, 2000000000, 2100000000, 2200000000, 2300000000, 2400000000, 2500000000, 2600000000, 2700000000, 2800000000, 2900000000, 3000000000, 3100000000, 3200000000, 3300000000, 3400000000, 3500000000, 3600000000, 3700000000, 3800000000, 3900000000, 4000000000, 4100000000, 4200000000, 4300000000, 4400000000, 4500000000, 4600000000, 4700000000, 4800000000, 4900000000, 5000000000, 5100000000, 5200000000, 5300000000, 5400000000, 5500000000, 5600000000, 5700000000, 5800000000, 5900000000, 6000000000, 6100000000, 6200000000, 6300000000]}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}, {"name": "Primary winding excitation", "frequency": 100000, "current": {"waveform": {"data": [-5, 5, -5], "time": [0, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 10, "offset": 0, "label": "Triangular"}}, "voltage": {"waveform": {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "time": [0, 0, 5e-06, 5e-06, 1e-05]}, "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "Rectangular"}}}]})"));
        settings.set_painter_simple_litz(false);
        // settings.set_painter_advanced_litz(false);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Web_0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);
            painter.paint_wire_with_current_density(wire, operatingPoint);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Web_1", "[support][painter][wire-density-painter][smoke-test]") {
        clear_databases();

        OpenMagnetics::Wire wire(json::parse(R"({"coating": {"breakdownVoltage": 13500, "grade": null, "material": "FEP", "numberLayers": 3, "temperatureRating": 155, "thickness": null, "thicknessLayers": 7.62e-05, "type": "insulated"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 1.0602875205865554e-06}, "conductingDiameter": null, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": null, "name": "Litz 15x0.3 - Grade 1 - Unserved", "numberConductors": 18, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0016300000000000002, "minimum": 0.0011557, "nominal": 0.012585990719707143}, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": null, "strand": {"coating": {"breakdownVoltage": 2500, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.004115}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Nearson", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 6.0 - Single Build", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 0.004186}, "outerHeight": null, "outerWidth": null, "standard": "NEMA MW 1000 C", "standardName": "6 AWG", "strand": null, "type": "round"}, "type": "litz"})"));
        settings.set_painter_simple_litz(false);
        // settings.set_painter_advanced_litz(false);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Web_1.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);
            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Wire_Painter_Current_Density_Web_2", "[support][painter][wire-density-painter][smoke-test]") {
        clear_databases();

        // OpenMagnetics::Wire wire(json::parse(R"({"coating": {"breakdownVoltage": null, "grade": null, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "bare"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 1.0602875205865554e-06}, "conductingDiameter": null, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": null, "name": "Litz 15x0.3 - Grade 1 - Unserved", "numberConductors": 15, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0016300000000000002, "minimum": 0.0011557, "nominal": null}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": null, "strand": {"coating": {"breakdownVoltage": 2200, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000303999999999, "minimum": 0.000296, "nominal": 0.00030000000000000003}, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.3 - Grade 1", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00033400000000000004, "minimum": 0.000319, "nominal": null}, "standard": "IEC 60317", "standardName": "0.3 mm", "type": "round"}, "type": "litz"})"));
        OpenMagnetics::Wire wire(json::parse(R"({"coating": {"breakdownVoltage": null, "grade": null, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "bare"}, "conductingArea": {"excludeMaximum": null, "excludeMinimum": null, "maximum": null, "minimum": null, "nominal": 1.0602875205865554e-06}, "conductingDiameter": null, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": null, "name": "Litz 15x0.3 - Grade 1 - Unserved", "numberConductors": 15, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.0016300000000000002, "minimum": 0.0011557, "nominal": 0.0015933066187962695}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": null, "strand": {"coating": {"breakdownVoltage": 2200, "grade": 1, "material": null, "numberLayers": null, "temperatureRating": null, "thickness": null, "thicknessLayers": null, "type": "enamelled"}, "conductingArea": null, "conductingDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.000303999999999, "minimum": 0.000296, "nominal": 0.00030000000000000003}, "conductingHeight": null, "conductingWidth": null, "edgeRadius": null, "manufacturerInfo": {"cost": null, "datasheetUrl": null, "family": null, "name": "Elektrisola", "orderCode": null, "reference": null, "status": null}, "material": "copper", "name": "Round 0.3 - Grade 1", "numberConductors": 1, "outerDiameter": {"excludeMaximum": null, "excludeMinimum": null, "maximum": 0.00033400000000000004, "minimum": 0.000319, "nominal": null}, "outerHeight": null, "outerWidth": null, "standard": "IEC 60317", "standardName": "0.3 mm", "strand": null, "type": "round"}, "type": "litz"})"));
        settings.set_painter_simple_litz(false);
        // settings.set_painter_advanced_litz(false);

        {
            auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Wire_Painter_Current_Density_Web_2.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, false, true);

            painter.paint_wire(wire);
            painter.export_svg();
        }
        settings.reset();
    }

    TEST_CASE("Test_Painter_Albach_2D_Field_With_Gap", "[support][painter][magnetic-field-painter][albach-2d]") {
        // Test case to visualize the magnetic field using the Albach 2D boundary value solver
        // Uses a PQ core with a center leg gap
        clear_databases();
        
        std::vector<int64_t> numberTurns = {8};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 2000;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        
        // Create a gap of 1mm in the center leg
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::SPREAD;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(100000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();
        
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        
        // Plot using BasicPainter (useAdvancedPainter = false)
        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Albach_2D_Field_With_Gap.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);  // Use BasicPainter
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);  // Enable fringing to show gap effects
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Albach_2D_Transformer_With_Gap", "[support][painter][magnetic-field-painter][albach-2d]") {
        // Test case with a transformer (two windings) and center leg gap
        clear_databases();
        
        std::vector<int64_t> numberTurns = {12, 6};
        std::vector<int64_t> numberParallels = {1, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 2;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 1500;
        std::string coreShape = "E 42/21/15";
        std::string coreMaterial = "3C95";
        
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0005);
        
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        
        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Round 0.80 - Grade 2")});
        wires.push_back({find_wire_by_name("Round 1.00 - Grade 2")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(150000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();
        
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        
        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Albach_2D_Transformer_With_Gap.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Albach_2D_Planar_Transformer", "[support][painter][magnetic-field-painter][albach-2d]") {
        // Test case with planar transformer and gap
        clear_databases();
        settings.set_coil_wind_even_if_not_fit(false);
        settings.set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns = {16, 4};
        std::vector<int64_t> numberParallels = {2, 2};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        double voltagePeakToPeak = 1000;
        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY, IsolationSide::SECONDARY};
        std::vector<size_t> stackUp = {0, 1, 0, 1, 0, 1, 0, 1};
        
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0003);
        auto core = OpenMagneticsTesting::get_quick_core("ELP 32/6/20", gapping, 1, "3C95");
        auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0006);
        wire.set_nominal_value_conducting_height(0.000070);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);
        wire.set_nominal_value_conducting_width(0.0024);
        wire.set_nominal_value_conducting_height(0.000070);
        wires.push_back(wire);

        OpenMagnetics::Coil coil;
        for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
            OpenMagnetics::Winding coilFunctionalDescription; 
            coilFunctionalDescription.set_number_turns(numberTurns[windingIndex]);
            coilFunctionalDescription.set_number_parallels(numberParallels[windingIndex]);
            coilFunctionalDescription.set_name(OpenMagnetics::to_string(isolationSides[windingIndex]));
            coilFunctionalDescription.set_isolation_side(isolationSides[windingIndex]);
            coilFunctionalDescription.set_wire(wires[windingIndex]);
            coil.get_mutable_functional_description().push_back(coilFunctionalDescription);
        }
        coil.set_bobbin(bobbin);
        coil.set_strict(false);

        coil.wind_by_planar_sections(stackUp, {{{0, 1}, 0.0001}}, 0.0001);
        coil.wind_by_planar_layers();
        coil.wind_by_planar_turns(0.0002, {{0, 0.0002}, {1, 0.0002}});
        coil.delimit_and_compact();

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(200000, 0.001, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Albach_2D_Planar_Transformer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Albach_2D_RM_Core_Inductor", "[support][painter][magnetic-field-painter][albach-2d]") {
        // Test case with RM core (pot-like shape) and larger gap
        clear_databases();
        
        std::vector<int64_t> numberTurns = {24};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 800;
        std::string coreShape = "RM 10";
        std::string coreMaterial = "3C90";
        
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.002);
        
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(75000, 0.002, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();
        
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        
        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Albach_2D_RM_Core_Inductor.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_RM_Core_Distributed_Gap_3", "[support][painter][magnetic-painter][pot-core][distributed-gap][smoke-test]") {
        // Test case with RM core (pot-like shape) with 3 distributed gaps
        clear_databases();
        
        std::vector<int64_t> numberTurns = {24};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 400;
        std::string coreShape = "RM 10";
        std::string coreMaterial = "3C90";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.002, 3);

        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(75000, 0.002, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_RM_Core_Distributed_Gap_3.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_magnetic_field_strength_model(MagneticFieldStrengthModels::ALBACH_2D);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_RM_Core_Single_Gap_Single_Turn", "[support][painter][magnetic-painter][pot-core][single-gap][albach-2d]") {
        // Simplified test: RM core with 1 gap and 1 turn to debug ALBACH_2D
        // Expected H_gap = N*I / l_gap = 1 * 50 / 0.001 = 50,000 A/m
        clear_databases();
        
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        std::string coreShape = "RM 10";
        std::string coreMaterial = "3C90";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);  // 1mm single gap

        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        // High current: 50A peak, 0A DC
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(75000, 50.0, 25.0, WaveformLabel::SINUSOIDAL, 0.5, 0.0, 0.0, turnsRatios);
        coil.delimit_and_compact();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_RM_Core_Single_Gap_Single_Turn.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_magnetic_field_strength_model(MagneticFieldStrengthModels::ALBACH_2D);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);

        painter.export_svg();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Albach_2D_Litz_Wire_Transformer", "[support][painter][magnetic-field-painter][albach-2d]") {
        // Test case with Litz wire transformer
        clear_databases();;
        
        std::vector<int64_t> numberTurns = {10, 5};
        std::vector<int64_t> numberParallels = {1, 1};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 500;
        std::string coreShape = "ETD 29/16/10";
        std::string coreMaterial = "3C95";
        
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0008);
        
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::INNER_OR_TOP;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        
        std::vector<OpenMagnetics::Wire> wires;
        wires.push_back({find_wire_by_name("Litz 75x0.3 - Grade 1 - Single Served")});
        wires.push_back({find_wire_by_name("Litz 270x0.02 - Grade 1 - Single Served")});
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(500000, 0.0005, 25, WaveformLabel::TRIANGULAR, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();
        
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        
        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Albach_2D_Litz_Wire_Transformer.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

    TEST_CASE("Test_Painter_Albach_2D_Rectangular_Wire", "[support][painter][magnetic-field-painter][albach-2d][rectangular]") {
        // Test case to visualize the magnetic field using ALBACH_2D with rectangular wire
        // Uses filamentary subdivision for rectangular conductor representation
        clear_databases();
        
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};
        std::vector<double> turnsRatios = {};
        uint8_t interleavingLevel = 1;
        int64_t numberStacks = 1;
        double voltagePeakToPeak = 500;
        std::string coreShape = "P 14/8";
        std::string coreMaterial = "3C97";
        
        // Create a small gap
        auto gapping = OpenMagneticsTesting::get_ground_gap(0);
        
        WindingOrientation sectionOrientation = WindingOrientation::OVERLAPPING;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::SPREAD;
        
        // Create rectangular wire
        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.001);
        wire.set_nominal_value_conducting_height(0.0003);
        wire.set_nominal_value_outer_height(0.00031);
        wire.set_nominal_value_outer_width(0.00104);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, coreShape, interleavingLevel, 
                                                          sectionOrientation, layersOrientation, turnsAlignment, 
                                                          sectionsAlignment, wires, false);
        auto core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(100000, 0.001, 25, WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        coil.delimit_and_compact();
        
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        
        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Albach_2D_Rectangular_Wire.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, false, false);
        settings.set_painter_mode(PainterModes::CONTOUR);
        settings.set_painter_logarithmic_scale(false);
        settings.set_painter_include_fringing(true);
        settings.set_painter_maximum_value_colorbar(std::nullopt);
        settings.set_painter_minimum_value_colorbar(std::nullopt);
        painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        REQUIRE(std::filesystem::exists(outFile));
        settings.reset();
    }

}  // namespace
