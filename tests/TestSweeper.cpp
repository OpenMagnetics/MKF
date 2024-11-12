#include "Sweeper.h"
#include "Painter.h"
#include "Settings.h"
#include "TestingUtils.h"
#include <UnitTest++.h>


SUITE(Sweeper) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

    TEST(Test_Sweeper_Impedance_Over_Frequency_Many_Turns) {
        settings->set_coil_wind_even_if_not_fit(true);
        std::vector<int64_t> numberTurns = {110, 110};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::WireWrapper> wires;
        auto wire = OpenMagnetics::find_wire_by_name("Round 0.15 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        std::vector<OpenMagnetics::CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 400000, 1000);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
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
        std::vector<OpenMagnetics::WireWrapper> wires;
        auto wire = OpenMagnetics::find_wire_by_name("Round 0.425 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        std::vector<OpenMagnetics::CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 10000);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Few_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
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
        std::vector<OpenMagnetics::WireWrapper> wires;
        auto wire = OpenMagnetics::find_wire_by_name("Round 2.50 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        std::vector<OpenMagnetics::CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Few_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
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
        std::vector<OpenMagnetics::WireWrapper> wires;
        auto wire = OpenMagnetics::find_wire_by_name("Round 1.40 - Grade 1");
        wires.push_back(wire);
        wires.push_back(wire);

        OpenMagnetics::WindingOrientation windingOrientation = OpenMagnetics::WindingOrientation::CONTIGUOUS;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        
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
        std::vector<OpenMagnetics::CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 100);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Impedance_Over_Frequency_Larger_Core_Many_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }

    TEST(Test_Sweeper_Resistance_Over_Frequency_Many_Turns) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round 0.25 - FIW 6");
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

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_resistance_over_frequency(magnetic, 0.1, 1000000, 100, 0);

        auto outFile = outputFilePath;

        outFile.append("Test_Sweeper_Resistance_Over_Frequency_Many_Turns.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep, true);
        painter.export_svg();
        CHECK(std::filesystem::exists(outFile));

        settings->reset();
    }
}

