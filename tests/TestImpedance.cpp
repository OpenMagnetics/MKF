#include "support/Painter.h"
#include "processors/Sweeper.h"
#include "physical_models/Impedance.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "support/Utils.h"
#include <UnitTest++.h>

using namespace MAS;
using namespace OpenMagnetics;


SUITE(Impedance) {
    double maximumError = 0.25;
    auto settings = Settings::GetInstance();
    TEST(Test_Impedance_0) {

        std::vector<int64_t> numberTurns = {54, 54};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 17/10.7/6.8";
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
        std::string coreMaterial = "80";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        double expectedSelfResonantFrequency = 1400000;
        settings->_debug = true;
        auto selfResonantFrequency = Impedance().calculate_self_resonant_frequency(magnetic);
        CHECK_CLOSE(expectedSelfResonantFrequency, selfResonantFrequency, expectedSelfResonantFrequency * maximumError);
        settings->_debug = false;

        {
            auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 4000000, 1000);

            auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
            auto outFile = outputFilePath;

            outFile.append("Test_Impedance_0.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_curve(impedanceSweep, true);
            painter.export_svg();
            CHECK(std::filesystem::exists(outFile));

        }
        {
            auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_Impedance_0_magnetic.svg";
            outFile.append(filename);
            settings->set_painter_include_fringing(false);
            Painter painter(outFile, false);

            painter.paint_core(magnetic);
            // painter.paint_bobbin(magnetic);
            // painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }

    TEST(Test_Impedance_Many_Turns) {

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
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        std::map<double, double> expectedImpedances = {
            {2000, 558},
            {5000, 1350},
            {10000, 2690},
            {25000, 6900},
            {50000, 15900}
        };
        settings->_debug = true;

        for (auto [frequency, expectedImpedance] : expectedImpedances) {
            auto impedance = Impedance().calculate_impedance(magnetic, frequency);
            CHECK_CLOSE(expectedImpedance, abs(impedance), expectedImpedance * maximumError);
        }

        // {
        //     auto impedanceSweep = Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 400000, 1000);

        //     auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        //     auto outFile = outputFilePath;

        //     outFile.append("Test_Impedance_Many_Turns.svg");
        //     std::filesystem::remove(outFile);
        //     Painter painter(outFile, false, true);
        //     painter.paint_curve(impedanceSweep, true);
        //     painter.export_svg();
        //     CHECK(std::filesystem::exists(outFile));

        // }
    }

    TEST(Test_Self_Resonant_Frequency_Many_Turns) {

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
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        double expectedSelfResonantFrequency = 180000;
        auto selfResonantFrequency = Impedance().calculate_self_resonant_frequency(magnetic);
        CHECK_CLOSE(expectedSelfResonantFrequency, selfResonantFrequency, expectedSelfResonantFrequency * maximumError);

    }
    TEST(Test_Impedance_Few_Turns) {

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
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        std::map<double, double> expectedImpedances = {
            {2000, 12.7},
            {5000, 31.8},
            {10000, 62.6},
            {25000, 153},
            {50000, 305},
        };

        for (auto [frequency, expectedImpedance] : expectedImpedances) {
            auto impedance = Impedance().calculate_impedance(magnetic, frequency);
            CHECK_CLOSE(expectedImpedance, abs(impedance), expectedImpedance * maximumError);
        }

    }

    TEST(Test_Impedance_Many_Turns_Larger_Core) {

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
        std::string coreMaterial = "A05";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        std::map<double, double> expectedImpedances = {
            {2000, 7.49},
            {5000, 19},
            {10000, 37.9},
            {25000, 93.9},
            {50000, 188},
        };

        for (auto [frequency, expectedImpedance] : expectedImpedances) {
            auto impedance = Impedance().calculate_impedance(magnetic, frequency);
            CHECK_CLOSE(expectedImpedance, abs(impedance), expectedImpedance * maximumError);
        }
    }

    TEST(Test_Impedance_Few_Turns_Larger_Core) {

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
        std::string coreMaterial = "A05";
        std::vector<CoreGap> gapping = {};
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        std::map<double, double> expectedImpedances = {
            {2000, 21.6},
            {5000, 54.1},
            {10000, 108},
            {25000, 300},
            {50000, 600},
        };

        for (auto [frequency, expectedImpedance] : expectedImpedances) {
            auto impedance = Impedance().calculate_impedance(magnetic, frequency);
            CHECK_CLOSE(expectedImpedance, abs(impedance), expectedImpedance * maximumError);
        }

    }

}

