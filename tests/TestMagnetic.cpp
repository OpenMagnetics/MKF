#include "physical_models/MagneticEnergy.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "support/Utils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <typeinfo>
#include <vector>

using json = nlohmann::json;

SUITE(Magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double max_error = 0.05;

    TEST(Test_Magnetic_Saturation_Current) {

        std::vector<int64_t> numberTurns = {18};
        std::vector<int64_t> numberParallels = {1};
        std::string shapeName = "PQ 65/44";

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName);

        // coil.wind({0, 1, 0}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::MagneticWrapper magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto saturationCurrentAt20 = magnetic.calculate_saturation_current(20);
        std::cout << "saturationCurrentAt20: " << saturationCurrentAt20 << std::endl;
        auto saturationCurrentAt100 = magnetic.calculate_saturation_current(100);
        std::cout << "saturationCurrentAt100: " << saturationCurrentAt100 << std::endl;
        CHECK(saturationCurrentAt100 < saturationCurrentAt20);
    }

}