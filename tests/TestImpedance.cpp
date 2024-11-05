#include "Impedance.h"
#include "TestingUtils.h"
#include "Utils.h"
#include <UnitTest++.h>


SUITE(Impedance) {
    TEST(Test_Impedance) {

        std::vector<int64_t> numberTurns = {110, 1};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::WireWrapper> wires;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName);

        int64_t numberStacks = 1;
        std::string coreMaterial = "A07";
        std::vector<OpenMagnetics::CoreGap> gapping = {};
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        double frequency = 2000;
        double expectedImpedance = 73800;

        auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
        std::cout << "expectedImpedance: " << expectedImpedance << std::endl;

    }

}

