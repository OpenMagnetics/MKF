#include "support/Utils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
using json = nlohmann::json;

SUITE(Constants) {
    TEST(AccessConstants) {
        CHECK_CLOSE(OpenMagnetics::constants.residualGap, 5.0e-6, 1e-8);
    }
}
