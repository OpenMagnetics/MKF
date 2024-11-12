#include "Impedance.h"
#include "TestingUtils.h"
#include "Utils.h"
#include <UnitTest++.h>


SUITE(Impedance) {
    double maximumError = 0.25;
    TEST(Test_Impedance_Many_Turns) {

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
        std::vector<OpenMagnetics::CoreGap> gapping = {};
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

        for (auto [frequency, expectedImpedance] : expectedImpedances) {
            auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
            CHECK_CLOSE(expectedImpedance, abs(impedance), expectedImpedance * maximumError);
        }

    }
    TEST(Test_Impedance_Few_Turns) {

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
        std::vector<OpenMagnetics::CoreGap> gapping = {};
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
            auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
            CHECK_CLOSE(expectedImpedance, abs(impedance), expectedImpedance * maximumError);
        }

    }

    TEST(Test_Impedance_Many_Turns_Larger_Core) {

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
        std::string coreMaterial = "A05";
        std::vector<OpenMagnetics::CoreGap> gapping = {};
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
            auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
            CHECK_CLOSE(expectedImpedance, abs(impedance), expectedImpedance * maximumError);
        }
    }

    TEST(Test_Impedance_Few_Turns_Larger_Core) {

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
        std::string coreMaterial = "A05";
        std::vector<OpenMagnetics::CoreGap> gapping = {};
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        std::map<double, double> expectedImpedances = {
            {2000, 21.6},
            {5000, 54.1},
            {10000, 108},
            {25000, 268},
            {50000, 537},
        };

        for (auto [frequency, expectedImpedance] : expectedImpedances) {
            auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
            CHECK_CLOSE(expectedImpedance, abs(impedance), expectedImpedance * maximumError);
        }

    }

}

