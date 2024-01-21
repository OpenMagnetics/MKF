#include "Constants.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
using json = nlohmann::json;

SUITE(Constants) {
    TEST(AccessConstants) {
        auto constants = OpenMagnetics::Constants();
        CHECK_CLOSE(constants.residualGap, 5.0e-6, 1e-8);
    }
}
