#include <UnitTest++.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>
#include "Constants.h"
#include "json.hpp"
using json = nlohmann::json;


SUITE(Constants)
{
    TEST(AccessConstants)
    {
        auto constants = OpenMagnetics::Constants();
        CHECK_CLOSE(constants.residualGap, 5.0e-6, 1e-8);
    }
}

