#include <UnitTest++.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>
#include "json.hpp"

using json = nlohmann::json;


int main(int, const char *[])
{
    return UnitTest::RunAllTests();
}
