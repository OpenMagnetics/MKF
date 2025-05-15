#include "physical_models/LeakageInductance.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "constructive_models/CoreWrapper.h"
#include "processors/InputsWrapper.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>
#include <cmath>


SUITE(LeakageInductance) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double maximumError = 0.2;
    auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
    bool plot = false;

    TEST(Test_Leakage_Inductance_E_0) {
        settings->reset();
        std::vector<int64_t> numberTurns({69, 69});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "E 42/33/20";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000352);
        wire.set_number_conductors(25);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 6.7e-6;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_E_1) {

        // settings->set_coil_wind_even_if_not_fit(true);
        // settings->set_coil_try_rewind(false);

        std::vector<int64_t> numberTurns({64, 20});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "E 42/33/20";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000352);
        wire.set_number_conductors(25);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);
        wire.set_number_conductors(225);
        wire.set_nominal_value_outer_diameter(0.001056);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 13e-6;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_E_2) {

        std::vector<int64_t> numberTurns({16, 6});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "E 42/33/20";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(370);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(370) * 0.000055);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        strandConductingDiameter.set_nominal(0.0001);
        strandOuterDiameter.set_nominal(0.00011);
        wire.set_strand(strand);
        wire.set_number_conductors(666);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(666) * 0.00011);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 4e-6;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Leakage_Inductance_E_2.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_E_3) {

        std::vector<int64_t> numberTurns({36, 26});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "E 65/32/27";
        uint8_t interleavingLevel = 2;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(650);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(650) * 0.000055);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        strandConductingDiameter.set_nominal(0.000071);
        strandOuterDiameter.set_nominal(0.000073);
        wire.set_strand(strand);
        wire.set_number_conductors(700);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(700) * 0.000073);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 8e-6;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_Parallels_No_Interleaving) {

        std::vector<int64_t> numberTurns({24, 6});
        std::vector<int64_t> numberParallels({2, 4});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "PQ 32/25";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(75);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(75) * 0.000055);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        wire.set_strand(strand);
        wire.set_number_conductors(225);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(225) * 0.000055);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 9e-6;

        if (plot) {
            auto outFile = outputFilePath;
            outFile.append("Test_Leakage_Inductance_Parallels_No_Interleaving.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_Parallels_Interleaving) {

        std::vector<int64_t> numberTurns({24, 6});
        std::vector<int64_t> numberParallels({2, 4});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "PQ 32/25";
        uint8_t interleavingLevel = 2;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(75);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(75) * 0.000055);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        wire.set_strand(strand);
        wire.set_number_conductors(225);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(225) * 0.000055);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 5e-6;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
        auto leakageMagneticField = OpenMagnetics::LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency);
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Leakage_Inductance_Parallels_Interleaving.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile, true);
            painter.paint_magnetic_field(OpenMagnetics::OperatingPoint(), magnetic, 1, leakageMagneticField);
            painter.paint_core(magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }

        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_ETD_0) {

        std::vector<int64_t> numberTurns({60, 59});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "ETD 39";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.000071);
        strandOuterDiameter.set_nominal(0.000073);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(135);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(135) * 0.000073);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        wire.set_number_conductors(69);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(69) * 0.000073);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 40e-6;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_PQ_26_0) {

        settings->set_coil_try_rewind(false);
        std::vector<int64_t> numberTurns({27, 3});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "PQ 26/25";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.000071);
        strandOuterDiameter.set_nominal(0.000073);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(60);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(60) * 0.000073);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        wire.set_strand(strand);
        wire.set_number_conductors(370);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(370) * 0.000055);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;

        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 86e-6;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();

        if (plot) {
            settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);

            auto leakageMagneticField = OpenMagnetics::LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency, 0, 1);

            auto outFile = outputFilePath;
            outFile.append("Test_Leakage_Inductance_PQ_26_0.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile, true);

            painter.paint_magnetic_field(OpenMagnetics::OperatingPoint(), magnetic, 1, leakageMagneticField);
            painter.paint_core(magnetic);
            painter.paint_bobbin(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_PQ_40_Horizontal) {

        std::vector<int64_t> numberTurns({20, 2});
        std::vector<int64_t> numberParallels({1, 3});
        std::vector<double> turnsRatios({10});
        std::string shapeName = "PQ 40/40";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.0021);
        wire.set_number_conductors(800);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);
        wire.set_nominal_value_outer_diameter(0.002469);
        wire.set_number_conductors(1000);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 9.9e-6;
        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_PQ_40_Vertical) {

        settings->set_coil_try_rewind(false);
        std::vector<int64_t> numberTurns({20, 2});
        std::vector<int64_t> numberParallels({1, 3});
        std::vector<double> turnsRatios({10});
        std::string shapeName = "PQ 40/40";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.0021);
        wire.set_number_conductors(800);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);
        wire.set_nominal_value_outer_diameter(0.002469);
        wire.set_number_conductors(1000);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 40e-6;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_Three_Windings) {
        settings->reset();
        std::vector<int64_t> numberTurns({50, 100, 25});
        std::vector<int64_t> numberParallels({1, 1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1], double(numberTurns[0]) / numberTurns[2]});
        std::string shapeName = "E 42/21/15";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000752);
        wire.set_number_conductors(350);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);
        wires.push_back(wire);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        double frequency = 100000;

        auto leakageInductance_01 = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
        auto leakageInductance_10 = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 1, 0).get_leakage_inductance_per_winding()[0].get_nominal().value();

        CHECK_CLOSE(leakageInductance_01, leakageInductance_10 * pow(double(numberTurns[0]) / numberTurns[1], 2), leakageInductance_01 * 0.01);

        auto leakageInductance_02 = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 2).get_leakage_inductance_per_winding()[0].get_nominal().value();
        auto leakageInductance_20 = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 2, 0).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(leakageInductance_02, leakageInductance_20 * pow(double(numberTurns[0]) / numberTurns[2], 2), leakageInductance_02 * 0.01);

        auto leakageInductance_12 = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 1, 2).get_leakage_inductance_per_winding()[0].get_nominal().value();
        auto leakageInductance_21 = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 2, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(leakageInductance_12, leakageInductance_21 * pow(double(numberTurns[1]) / numberTurns[2], 2), leakageInductance_12 * 0.01);
        settings->reset();
    }

    TEST(Test_Leakage_Inductance_T_0) {
        settings->reset();
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns({10, 200});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "T 48/28/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        settings->set_coil_try_rewind(false);

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;

        wire.set_nominal_value_outer_diameter(0.001);
        wire.set_nominal_value_conducting_diameter(0.00095);
        wire.set_material("copper");
        wire.set_number_conductors(1);
        wire.set_type(OpenMagnetics::WireType::ROUND);
        wires.push_back(wire);

        wire.set_nominal_value_outer_diameter(0.0008);
        wire.set_nominal_value_conducting_diameter(0.00075);
        wire.set_material("copper");
        wire.set_number_conductors(1);
        wire.set_type(OpenMagnetics::WireType::ROUND);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        std::vector<double> proportionPerWinding = {16.185 / (222.42 + 16.185), 222.42 / (222.42 + 16.185)};
        std::vector<size_t> pattern = {0, 1};
        coil.wind(proportionPerWinding, pattern);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 0.00143;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 1, 0).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        if (plot) {

            settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
            auto leakageMagneticField = OpenMagnetics::LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency, 0, 1);

            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Leakage_Inductance_T_0.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile, true);
            painter.paint_magnetic_field(OpenMagnetics::OperatingPoint(), magnetic, 1, leakageMagneticField);
            painter.paint_core(magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }


        settings->reset();
    }

    TEST(Test_Leakage_Inductance_T_1) {
        settings->reset();
        OpenMagnetics::clear_databases();
        std::vector<int64_t> numberTurns({10, 5});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "T 48/28/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        settings->set_coil_try_rewind(false);

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0038);
        wire.set_nominal_value_conducting_height(0.002);
        wire.set_nominal_value_outer_width(0.0039);
        wire.set_nominal_value_outer_height(0.002076);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::RECTANGULAR);
        wires.push_back(wire);

        wire.set_nominal_value_conducting_width(0.0038);
        wire.set_nominal_value_conducting_height(0.001);
        wire.set_nominal_value_outer_width(0.0039);
        wire.set_nominal_value_outer_height(0.001076);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::RECTANGULAR);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        double frequency = 100000;
        double expectedLeakageInductance = 5e-6;

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
        auto leakageMagneticField = OpenMagnetics::LeakageInductance().calculate_leakage_magnetic_field(magnetic, frequency, 0, 1);
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
        if (true) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Leakage_Inductance_T_1.svg");
            std::filesystem::remove(outFile);
            OpenMagnetics::Painter painter(outFile, true);
            painter.paint_magnetic_field(OpenMagnetics::OperatingPoint(), magnetic, 1, leakageMagneticField);
            painter.paint_core(magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }


        settings->reset();
    }
}
