#include "constructive_models/Mas.h"
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

using namespace MAS;
using namespace OpenMagnetics;

using json = nlohmann::json;

SUITE(Mas) {
    TEST(Test_Expand_Magnetic) {

        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/example_basic.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        CHECK(!mas.get_magnetic().get_core().get_processed_description());

        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
        CHECK(magnetic.get_core().get_processed_description());

    }

}