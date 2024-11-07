#include "Sweeper.h"
#include "Painter.h"
#include "TestingUtils.h"
#include <UnitTest++.h>


SUITE(Sweeper) {
    auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");

    TEST(Test_Sweeper) {
        std::vector<int64_t> numberTurns = {110, 110};
        std::vector<int64_t> numberParallels = {1, 1};
        std::string shapeName = "T 12.5/7.5/5";
        std::vector<OpenMagnetics::WireWrapper> wires;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName);

        int64_t numberStacks = 1;
        std::string coreMaterial = "A07";
        std::vector<OpenMagnetics::CoreGap> gapping;
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto impedanceSweep = OpenMagnetics::Sweeper().sweep_impedance_over_frequency(magnetic, 1000, 200000, 1000);
        std::cout << impedanceSweep.get_x_points().size() << std::endl;
        std::cout << impedanceSweep.get_y_points().size() << std::endl;


        auto outFile = outputFilePath;

        outFile.append("Test_Sweep_Impedance_Over_Frequency.svg");
        std::filesystem::remove(outFile);
        OpenMagnetics::Painter painter(outFile, false, true);
        painter.paint_curve(impedanceSweep);
        painter.export_svg();

    }

}

