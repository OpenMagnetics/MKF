#include "MasWrapper.h"
#include "Settings.h"
#include "TestingUtils.h"
#include "Utils.h"
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
    TEST(Test_Expand_Magnetic) {

        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/example_basic.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        CHECK(!mas.get_magnetic().get_core().get_processed_description());

        auto magnetic = OpenMagnetics::MasWrapper::expand_magnetic(mas.get_magnetic());
        auto inputs = OpenMagnetics::MasWrapper::expand_inputs(magnetic, mas.get_inputs());
        CHECK(magnetic.get_core().get_processed_description());

    }

}