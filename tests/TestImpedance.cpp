#include "Impedance.h"
#include "TestingUtils.h"
#include <UnitTest++.h>


SUITE(Impedance) {
    TEST(Test_Impedance) {

        std::vector<int64_t> numberTurns = {16, 1};
        std::vector<int64_t> numberParallels = {2, 1};
        std::string shapeName = "PQ 35/35";
        std::vector<OpenMagnetics::WireWrapper> wires;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C90";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        double frequency = 1e6;

        auto impedance = OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency);
        std::cout << impedance.real() << std::endl;
        std::cout << impedance.imag() << std::endl;

    }

}

