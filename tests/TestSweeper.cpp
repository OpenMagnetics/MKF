#include "processors/Sweeper.h"
#include "support/Painter.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include <UnitTest++.h>

using namespace MAS;
using namespace OpenMagnetics;


SUITE(Sweeper) {
    auto settings = Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

    TEST(Test_Sweeper_Impedance_Over_Frequency_Many_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {110, 110};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 0.15 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         1,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "A07";
        std::vector<CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 400000, 1000);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Q_Factor_Over_Frequency_Many_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {110, 110};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 0.15 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         1,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "A07";
        std::vector<CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_q_factor_over_frequency(magnetic, 1000, 400000, 1000);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Q_Factor_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Impedance_Over_Frequency_Few_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {18, 18};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 0.425 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         1,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "A07";
        std::vector<CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 10000);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Few_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Few_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {9, 9};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 36/23/15";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 2.50 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         1,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "A07";
        std::vector<CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Few_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Many_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {17, 17};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 36/23/15";
        std::vector<OpenMagnetics::Wire> wires;
        auto wire = find_wire_by_name("Round 1.40 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        WindingOrientation windingOrientation = WindingOrientation::CONTIGUOUS;
        WindingOrientation layersOrientation = WindingOrientation::OVERLAPPING;
        CoilAlignment sectionsAlignment = CoilAlignment::CENTERED;
        CoilAlignment turnsAlignment = CoilAlignment::CENTERED;
        
        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         1,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "A07";
        std::vector<CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Winding_Resistance_Over_Frequency_Many_Turns) {
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_winding_resistance_over_frequency(magnetic, 0.1, 1000000, 100, 0);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Winding_Resistance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Magnetizing_Inductance_Over_Frequency_Many_Turns) {
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "Kool Mµ 75";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_magnetizing_inductance_over_frequency(magnetic, 10000, 10000000, 100, 25);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Magnetizing_Inductance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Magnetizing_Inductance_Over_Temperature_Many_Turns) {
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "Kool Mµ 75";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_magnetizing_inductance_over_temperature(magnetic, -40, 200, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Magnetizing_Inductance_Over_Temperature_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, false);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Magnetizing_Inductance_Over_DC_Bias_Powder) {
        std::vector<int64_t> numberTurns = {12, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "Kool Mµ 75";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_magnetizing_inductance_over_dc_bias(magnetic, 0, 200, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Magnetizing_Inductance_Over_DC_Bias_Powder.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, false);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Magnetizing_Inductance_Over_DC_Bias_Ferrite) {
        std::vector<int64_t> numberTurns = {12, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_magnetizing_inductance_over_dc_bias(magnetic, 0, 260, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Magnetizing_Inductance_Over_DC_Bias_Ferrite.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, false);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Resistance_Over_Frequency_Many_Turns) {
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = Sweeper().sweep_resistance_over_frequency(magnetic, 10000, 10000000, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Resistance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Core_Resistance_Over_Frequency_Many_Turns) {
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "N87";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto coreSweep = Sweeper().sweep_core_resistance_over_frequency(magnetic, 10000, 1200000, 20);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Core_Resistance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(coreSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Core_Losses_Over_Frequency_Many_Turns) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "N87";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);


        double frequency = 57026;
        double magnetizingInductance = 100e-6;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 10;
        double dutyCycle = 0.5;
        double dcCurrent = 0;


        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        auto coreSweep = Sweeper().sweep_core_losses_over_frequency(magnetic, inputs.get_operating_points()[0], 10000, 1200000, 20);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Core_Losses_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(coreSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Core_Losses_Over_Frequency_Non_Steinmetz) {
        double temperature = 100;
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "A121";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);


        double frequency = 57026;
        double magnetizingInductance = 100e-6;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 10;
        double dutyCycle = 0.5;
        double dcCurrent = 0;


        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        auto coreSweep = Sweeper().sweep_core_losses_over_frequency(magnetic, inputs.get_operating_points()[0], 10000, 1200000, 20);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Core_Losses_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(coreSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Winding_Losses_Over_Frequency_Many_Turns) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "N87";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);


        double frequency = 57026;
        double magnetizingInductance = 100e-6;
        WaveformLabel waveShape = WaveformLabel::TRIANGULAR;
        double peakToPeak = 10;
        double dutyCycle = 0.5;
        double dcCurrent = 0;


        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        auto coreSweep = Sweeper().sweep_winding_losses_over_frequency(magnetic, inputs.get_operating_points()[0], 10000, 1200000, 20);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Winding_Losses_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(coreSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Resistance_Over_Frequency_Web_0) {
        OpenMagnetics::Magnetic magnetic = json::parse(R"({"coil":{"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.00635,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0079375,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.0007917304360898403,"coordinates":[0.015875,0,0],"height":null,"radialHeight":0.015875,"sectionsAlignment":"spread","sectionsOrientation":"contiguous","shape":"round","width":null}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":4,"wire":{"coating":{"breakdownVoltage":1000,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000505,"minimum":0.000495,"nominal":0.0005},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 0.5 - Grade 1","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000544,"minimum":0.0005239999999990001,"nominal":null},"outerHeight":null,"outerWidth":null,"standard":"IEC 60317","standardName":"0.5 mm","strand":null,"type":"round"}},{"connections":null,"isolationSide":"secondary","name":"Secondary","numberParallels":1,"numberTurns":8,"wire":{"coating":{"breakdownVoltage":5000,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001151},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 17.0 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.001191},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"17 AWG","strand":null,"type":"round"}}],"layersDescription":[{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.00026700000000000074,89.97744260415998],"dimensions":[0.0005339999999995018,179.95488520831995],"fillingFactor":0.04433090220254894,"insulationMaterial":null,"name":"Primary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Primary"}],"section":"Primary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0079375,179.909770417],"dimensions":[0.015875,0.09022958336009401],"fillingFactor":1,"insulationMaterial":null,"name":"Insulation between Primary and Primary section 1 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Primary and Primary section 1","turnsAlignment":"spread","type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0005955000000000005,270.02255739584],"dimensions":[0.0011910000000000035,179.95488520831995],"fillingFactor":0.20659336246431567,"insulationMaterial":null,"name":"Secondary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Secondary"}],"section":"Secondary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"polar","coordinates":[0.0079375,359.954885208],"dimensions":[0.015875,0.09022958336009401],"fillingFactor":1,"insulationMaterial":null,"name":"Insulation between Secondary and Secondary section 3 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Secondary and Secondary section 3","turnsAlignment":"spread","type":"insulation","windingStyle":null}],"sectionsDescription":[{"coordinateSystem":"polar","coordinates":[0.00026700000000000074,89.97744260415998],"dimensions":[0.0005339999999995018,179.95488520831995],"fillingFactor":0.002882066680966784,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0,0],"name":"Primary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"polar","coordinates":[0.0079375,179.9097704166399],"dimensions":[0.015875,0.09022958336009401],"fillingFactor":1,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Primary and Primary section 1","partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"polar","coordinates":[0.0005955000000000005,270.02255739584],"dimensions":[0.0011910000000000035,179.95488520831995],"fillingFactor":0.028673125080251508,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0,0],"name":"Secondary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Secondary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"polar","coordinates":[0.0079375,359.95488520832],"dimensions":[0.015875,0.09022958336009401],"fillingFactor":1,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Secondary and Secondary section 3","partialWindings":[],"type":"insulation","windingStyle":null}],"turnsDescription":[{"additionalCoordinates":[[0.02958105679173332,0.012249463991711979]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.014420499559776792,0.005971503700616564],"dimensions":[0.0005339999999995,0.0005339999999995],"layer":"Primary section 0 layer 0","length":0.05882761047701695,"name":"Primary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":22.494360651,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[0.012261109128712661,0.029576231875169615]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.00597718060033567,0.014418147456277829],"dimensions":[0.0005339999999995,0.0005339999999995],"layer":"Primary section 0 layer 0","length":0.058827610477016956,"name":"Primary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":67.48308195300001,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[-0.012237816956040178,0.029585877123223095]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.005965825875312337,0.014422849428093389],"dimensions":[0.0005339999999995,0.0005339999999995],"layer":"Primary section 0 layer 0","length":0.05882761047701695,"name":"Primary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":112.47180325500001,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[-0.029571402374279854,0.012272752365237236]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.01441579311796108,0.005982856573589742],"dimensions":[0.0005339999999995,0.0005339999999995],"layer":"Primary section 0 layer 0","length":0.05882761047701695,"name":"Primary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":157.460524557,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":[[-0.031719323454203245,-0.006333710591559409]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.014983704154163592,-0.0029919441957531033],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":191.292295117,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[-0.026882798159381348,-0.01798739929428273]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.01269900649166862,-0.00849696147893812],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.06089163685042544,"name":"Secondary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":213.786655768,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[-0.01795563541598219,-0.02690402436547971]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.008481956727782842,-0.012709033413994134],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":236.281016419,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[-0.006296242411054827,-0.031726782089446674]],"angle":null,"coordinateSystem":"cartesian","coordinates":[-0.0029742448229185584,-0.01498722749488184],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":258.77537707,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[0.006321222177019618,-0.031721814583008254]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.002986044867254216,-0.014984880923809327],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.06089163685042546,"name":"Secondary parallel 0 turn 4","orientation":"clockwise","parallel":0,"rotation":281.269737721,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[0.01797681412055799,-0.02688987772982441]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.008491961211144235,-0.012702350768201205],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 5","orientation":"clockwise","parallel":0,"rotation":303.76409837200003,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[0.026896953132013648,-0.017966226160890365]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.012705693075716947,-0.008486959627315216],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 6","orientation":"clockwise","parallel":0,"rotation":326.25845902300006,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":[[0.031724300794811054,-0.00630873278324186]],"angle":null,"coordinateSystem":"cartesian","coordinates":[0.014986055370741387,-0.0029801450761788814],"dimensions":[0.001191,0.001191],"layer":"Secondary section 0 layer 0","length":0.060891636850425444,"name":"Secondary parallel 0 turn 7","orientation":"clockwise","parallel":0,"rotation":348.7528196740001,"section":"Secondary section 0","winding":"Secondary"}]},"core":{"distributorsInfo":null,"functionalDescription":{"coating":null,"gapping":[],"material":"A07","numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0635},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03175},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 64/32/12.7","type":"standard"},"type":"toroidal","magneticCircuit":"closed"},"geometricalDescription":[{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"A07","rotation":[1.5707963267948966,1.5707963267948966,0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0635},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.03175},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 64/32/12.7","type":"standard"},"type":"toroidal"}],"manufacturerInfo":null,"name":"Custom","processedDescription":{"columns":[{"area":0.000202,"coordinates":[0,0,0],"depth":0.0127,"height":0.1496183501272139,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"central","width":0.015875}],"depth":0.0127,"effectiveParameters":{"effectiveArea":0.00020161250000000003,"effectiveLength":0.14961835012721392,"effectiveVolume":0.000030164929615022918,"minimumArea":0.0002016125},"height":0.0635,"width":0.0635,"windingWindows":[{"angle":360,"area":0.0007917304360898403,"coordinates":[0.015875,0],"height":null,"radialHeight":0.015875,"sectionsAlignment":null,"sectionsOrientation":null,"shape":null,"width":null}]}},"manufacturerInfo":{"name":"OpenMagnetics","reference":"My custom magnetic"}})");

        auto sweep = Sweeper().sweep_winding_resistance_over_frequency(magnetic, 1000, 4e+06, 1000, 0, 25);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Resistance_Over_Frequency_Web_0.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(sweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Winding_Losses_Over_Frequency_Web_0) {

        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/bug_winding_losses_sweep.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();
        auto operatingPoint = inputs.get_operating_point(0);

        auto sweep = Sweeper::sweep_winding_losses_over_frequency(magnetic, operatingPoint, 1000, 4e+06, 200, 25, "log");

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Winding_Losses_Over_Frequency_Web_0.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(sweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Winding_Losses_Over_Frequency_Web_1) {

        OpenMagnetics::Magnetic magnetic = json::parse(R"({"coil":{"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.005200500000000001,"columnShape":"round","columnThickness":0.0010000000000000009,"columnWidth":0.005200500000000001,"coordinates":[0,0,0],"pins":null,"wallThickness":0.0012375000000000008,"windingWindows":[{"angle":null,"area":0.000029798124999999987,"coordinates":[0.006938000000000001,0,0],"height":0.008575,"radialHeight":null,"sectionsAlignment":"inner or top","sectionsOrientation":"overlapping","shape":"rectangular","width":0.003474999999999999}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":40,"wire":{"coating":{"breakdownVoltage":1340,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.8322475214082728e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000484999999999,"minimum":0.00047799999999900006,"nominal":0.000483},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 24.5 - Single Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000522,"minimum":0.000503,"nominal":0.000513},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"24.5 AWG","strand":null,"type":"round"}},{"connections":null,"isolationSide":"secondary","name":"Secondary","numberParallels":1,"numberTurns":12,"wire":{"coating":{"breakdownVoltage":5000,"grade":2,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":3.247220852585116e-7},"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000643},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Nearson","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 22.0 - Heavy Build","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.000701},"outerHeight":null,"outerWidth":null,"standard":"NEMA MW 1000 C","standardName":"22 AWG","strand":null,"type":"round"}}],"layersDescription":[{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,0],"dimensions":[0.000513,0.0084755],"fillingFactor":0.8375510204081633,"insulationMaterial":null,"name":"Primary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[0.35],"winding":"Primary"}],"section":"Primary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,2.999999999721975e-9],"dimensions":[0.000513,0.00842838],"fillingFactor":0.777725947521866,"insulationMaterial":null,"name":"Primary section 0 layer 1","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[0.325],"winding":"Primary"}],"section":"Primary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,2.999999999721975e-9],"dimensions":[0.000513,0.00842838],"fillingFactor":0.777725947521866,"insulationMaterial":null,"name":"Primary section 0 layer 2","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[0.325],"winding":"Primary"}],"section":"Primary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.006752,0],"dimensions":[0.000025,0.008575],"fillingFactor":1,"insulationMaterial":null,"name":"Insulation between Primary and Secondary section 1 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Primary and Secondary section 1","turnsAlignment":"spread","type":"insulation","windingStyle":null},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,-4.336808689942018e-19],"dimensions":[0.000701,0.008561413],"fillingFactor":0.9809912536443149,"insulationMaterial":null,"name":"Secondary section 0 layer 0","orientation":"overlapping","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Secondary"}],"section":"Secondary section 0","turnsAlignment":"spread","type":"conduction","windingStyle":"windByConsecutiveTurns"},{"additionalCoordinates":null,"coordinateSystem":"cartesian","coordinates":[0.007478000000000001,0],"dimensions":[0.000025,0.008575],"fillingFactor":1,"insulationMaterial":null,"name":"Insulation between Secondary and Primary section 3 layer 0","orientation":"overlapping","partialWindings":[],"section":"Insulation between Secondary and Primary section 3","turnsAlignment":"spread","type":"insulation","windingStyle":null}],"sectionsDescription":[{"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,0],"dimensions":[0.001539,0.0084755],"fillingFactor":0.5146674334458955,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0,0],"name":"Primary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Primary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"cartesian","coordinates":[0.006752000000000002,0],"dimensions":[0.000025,0.008575],"fillingFactor":1,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Primary and Secondary section 1","partialWindings":[],"type":"insulation","windingStyle":null},{"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,0],"dimensions":[0.000701,0.008561413],"fillingFactor":0.6458557114859496,"layersAlignment":null,"layersOrientation":"overlapping","margin":[0,0],"name":"Secondary section 0","partialWindings":[{"connections":null,"parallelsProportion":[1],"winding":"Secondary"}],"type":"conduction","windingStyle":"windByConsecutiveTurns"},{"coordinateSystem":"cartesian","coordinates":[0.007478000000000003,0],"dimensions":[0.000025,0.008575],"fillingFactor":1,"layersAlignment":null,"layersOrientation":"overlapping","margin":null,"name":"Insulation between Secondary and Primary section 3","partialWindings":[],"type":"insulation","windingStyle":null}],"turnsDescription":[{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,0.00398125],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,0.00336875],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,0.00275625],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,0.00214375],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,0.00153125],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 4","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,0.0009187500000000001],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 5","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,0.0003062500000000001],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 6","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,-0.0003062499999999999],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 7","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,-0.0009187499999999999],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 8","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,-0.0015312499999999998],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 9","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,-0.00214375],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 10","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,-0.00275625],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 11","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,-0.00336875],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 12","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005457000000000002,-0.00398125],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 0","length":0.034287342221279014,"name":"Primary parallel 0 turn 13","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,0.003957693],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 14","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,0.003298078],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 15","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,0.0026384629999999997],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 16","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,0.0019788479999999996],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 17","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,0.0013192329999999995],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 18","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,0.0006596179999999995],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 19","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,2.9999999995051346e-9],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 20","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,-0.0006596120000000005],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 21","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,-0.0013192270000000005],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 22","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,-0.0019788420000000006],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 23","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,-0.0026384570000000007],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 24","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,-0.003298072000000001],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 25","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.005970000000000001,-0.0039576870000000005],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 1","length":0.03751061628386214,"name":"Primary parallel 0 turn 26","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,0.003957693],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 27","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,0.003298078],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 28","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,0.0026384629999999997],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 29","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,0.0019788479999999996],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 30","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,0.0013192329999999995],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 31","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,0.0006596179999999995],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 32","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,2.9999999995051346e-9],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 33","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,-0.0006596120000000005],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 34","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,-0.0013192270000000005],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 35","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,-0.0019788420000000006],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 36","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,-0.0026384570000000007],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 37","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,-0.003298072000000001],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 38","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.006483000000000002,-0.0039576870000000005],"dimensions":[0.000513,0.000513],"layer":"Primary section 0 layer 2","length":0.04073389034644527,"name":"Primary parallel 0 turn 39","orientation":"clockwise","parallel":0,"rotation":0,"section":"Primary section 0","winding":"Primary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,0.003930206499999999],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 0","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,0.0032156234999999992],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 1","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,0.0025010404999999993],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 2","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,0.0017864574999999993],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 3","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,0.0010718744999999993],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 4","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,0.00035729149999999934],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 5","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,-0.00035729150000000064],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 6","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,-0.0010718745000000006],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 7","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,-0.0017864575000000006],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 8","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,-0.0025010405000000006],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 9","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,-0.0032156235000000005],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 10","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"},{"additionalCoordinates":null,"angle":null,"coordinateSystem":"cartesian","coordinates":[0.007115000000000002,-0.0039302065],"dimensions":[0.000701,0.000701],"layer":"Secondary section 0 layer 0","length":0.04470486346058277,"name":"Secondary parallel 0 turn 11","orientation":"clockwise","parallel":0,"rotation":0,"section":"Secondary section 0","winding":"Secondary"}]},"core":{"distributorsInfo":[{"cost":2.25,"country":"USA","distributedArea":"International","email":null,"link":"https://www.mouser.es/ProductDetail/EPCOS-TDK/B65811J65A87?qs=fWY6JreybtknTb3trPTKFA%3D%3D","name":"Mouser","phone":null,"quantity":539,"reference":"871-B65811J65A87","updatedAt":"03/02/2025"}],"functionalDescription":{"coating":null,"gapping":[{"area":0.000056,"coordinates":[0,0.000775,0],"distanceClosestNormalSurface":0.003976,"distanceClosestParallelSurface":0.004475,"length":0.00155,"sectionDimensions":[0.008401,0.008401],"shape":"round","type":"subtractive"},{"area":0.000033,"coordinates":[0.010026,0,0],"distanceClosestNormalSurface":0.005521,"distanceClosestParallelSurface":0.004475,"length":0.00001,"sectionDimensions":[0.0027,0.012223],"shape":"irregular","type":"residual"},{"area":0.000033,"coordinates":[-0.010026,0,0],"distanceClosestNormalSurface":0.005521,"distanceClosestParallelSurface":0.004475,"length":0.00001,"sectionDimensions":[0.0027,0.012223],"shape":"irregular","type":"residual"}],"material":"3C97","numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0232,"minimum":0.0223,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00825,"minimum":0.00815,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.011,"minimum":0.0106,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00565,"minimum":0.0054,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0177,"minimum":0.017,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00855,"minimum":0.00825,"nominal":null},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":0.0095,"nominal":null},"H":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0046,"minimum":0.0044,"nominal":null},"J":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0197,"minimum":0.0189,"nominal":null},"R":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0003,"minimum":null,"nominal":null}},"family":"rm","familySubtype":"3","magneticCircuit":"open","name":"RM 8","type":"standard"},"type":"two-piece set","magneticCircuit":"open"},"geometricalDescription":[{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":[{"coordinates":[0,0.000775,0],"length":0.00155}],"material":"3C97","rotation":[3.141592653589793,3.141592653589793,0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0232,"minimum":0.0223,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00825,"minimum":0.00815,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.011,"minimum":0.0106,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00565,"minimum":0.0054,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0177,"minimum":0.017,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00855,"minimum":0.00825,"nominal":null},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":0.0095,"nominal":null},"H":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0046,"minimum":0.0044,"nominal":null},"J":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0197,"minimum":0.0189,"nominal":null},"R":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0003,"minimum":null,"nominal":null}},"family":"rm","familySubtype":"3","magneticCircuit":"open","name":"RM 8","type":"standard"},"type":"half set"},{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"3C97","rotation":[0,0,0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0232,"minimum":0.0223,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00825,"minimum":0.00815,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.011,"minimum":0.0106,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00565,"minimum":0.0054,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0177,"minimum":0.017,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00855,"minimum":0.00825,"nominal":null},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":0.0095,"nominal":null},"H":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0046,"minimum":0.0044,"nominal":null},"J":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0197,"minimum":0.0189,"nominal":null},"R":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0003,"minimum":null,"nominal":null}},"family":"rm","familySubtype":"3","magneticCircuit":"open","name":"RM 8","type":"standard"},"type":"half set"}],"manufacturerInfo":null,"name":"Custom","processedDescription":{"columns":[{"area":0.000056,"coordinates":[0,0,0],"depth":0.008401,"height":0.01105,"minimumDepth":null,"minimumWidth":null,"shape":"round","type":"central","width":0.008401},{"area":0.000033,"coordinates":[0.010026,0,0],"depth":0.012223,"height":0.01105,"minimumDepth":null,"minimumWidth":null,"shape":"irregular","type":"lateral","width":0.0027},{"area":0.000033,"coordinates":[-0.010026,0,0],"depth":0.012223,"height":0.01105,"minimumDepth":null,"minimumWidth":null,"shape":"irregular","type":"lateral","width":0.0027}],"depth":0.01735,"effectiveParameters":{"effectiveArea":0.00005202273495154223,"effectiveLength":0.03542796663543134,"effectiveVolume":0.000001843059718147126,"minimumArea":0.00003951338160052563},"height":0.016399999999999998,"thermalResistance":null,"width":0.02275,"windingWindows":[{"angle":null,"area":0.00004944875,"coordinates":[0.004200000000000001,0],"height":0.01105,"radialHeight":null,"sectionsAlignment":null,"sectionsOrientation":null,"shape":null,"width":0.004475}]}},"manufacturerInfo":{"name":"OpenMagnetics","reference":"My custom magnetic"}})");
        OperatingPoint operatingPoint = json::parse(R"({"name":"Op. Point No. 1","conditions":{"ambientTemperature":100},"excitationsPerWinding":[{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,0.4511310569983995,0.16293015292554872,0.08352979924600704],"frequencies":[0,100000,300000,500000,700000]}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular","acEffectiveFrequency":591485.5360118389,"effectiveFrequency":553357.3374711702,"peak":70.5,"rms":51.572309672924284,"thd":0.4833151484524849},"harmonics":{"amplitudes":[24.2890625,57.92076613061847,19.27588907896988,11.528257939603122,8.194467538528329,6.331896912839248,5.137996046859012,4.304077056139349,3.6860723299088454,3.207698601961777,2.8247804703632298],"frequencies":[0,100000,300000,500000,700000,900000,1100000,1300000,1500000,1700000,1900000]}}},{"current":{"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":110746.40291779551,"average":-1.2073675392798577e-14,"dutyCycle":0.4999999999999997,"effectiveFrequency":110746.40291779551,"label":"Triangular","offset":-7.993605777301127e-15,"peak":5,"peakToPeak":9.999999999999984,"phase":null,"rms":2.8874560332150576,"thd":0.12151487440704967},"waveform":{"ancillaryLabel":null,"data":[-5,5,-5],"numberPeriods":null,"time":[0,0.000005,0.00001]}},"frequency":100000,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[50,63.66836927259825,0,21.239845973980135,0,12.764409619153927,0,9.139463837994452,0,7.131406462516141,0,5.858362437526288,0,4.981163385580113,0,4.341543620727607,0,3.855727306176788,0,3.475223814311625,0,3.1700258281158034,0,2.9205737407387846,0,2.713577235072968,0,2.5396940337952887,0,2.3921692655409155,0,2.2660016324287593,0,2.1574129006256673,0,2.0634993994178332,0,1.981996677151262,0,1.9111167121328483,0,1.8494329245544314,0,1.7957974501020666,0,1.7492806654596298,0,1.7091263687047202,0,1.6747181778374791,0,1.6455541098570983,0,1.6212272284541696,0,1.601410873150707,0,1.5858474127267403,0,1.5743397677532718,0,1.5667451638950334,0,1.5629707375645523],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]},"processed":{"acEffectiveFrequency":599606.589028566,"average":-50,"dutyCycle":0.5,"effectiveFrequency":489576.7298435604,"label":"Custom","offset":-50,"peak":100,"peakToPeak":100,"phase":null,"rms":70.71067811865476,"thd":0.4831695829659298},"waveform":{"ancillaryLabel":null,"data":[0,0,-100,-100,0],"numberPeriods":null,"time":[0,0.000005,0.000005,0.00001,0.00001]}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}}]})");

        auto sweep = Sweeper::sweep_winding_losses_over_frequency(magnetic, operatingPoint, 1000, 4e+06, 200, 25, "log");

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Winding_Losses_Over_Frequency_Web_1.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_curve(sweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }
}

