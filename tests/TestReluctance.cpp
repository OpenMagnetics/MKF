#include <UnitTest++.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>
#include "Reluctance.h"


SUITE(Reluctance)
{
    TEST(ZhangModel)
    {
        auto constants = OpenMagnetics::Constants();
        CHECK_CLOSE(constants.residualGap, 5.0e-6, 1e-8);
    }
}

