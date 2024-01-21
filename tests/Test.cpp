#include "TestReporterStdout.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

using json = nlohmann::json;

namespace UnitTest {
class ListFilter {
    char** list;
    int n;

  public:
    ListFilter(char** list_, int count) {
        list = list_;
        n = count;
    }
    bool operator()(const Test* const t) const {
        for (int i = 0; i < n; ++i) {
            std::string dot_cpp_appended = std::string(list[i]) + ".cpp";
            if (!strcmp(t->m_details.testName, list[i]) || !strcmp(t->m_details.suiteName, list[i]) ||
                !strcmp(t->m_details.filename, list[i]) ||
                !strcmp(t->m_details.filename, dot_cpp_appended.c_str())) {
                // erring on the side of matching more tests
                return true;
            }
        }
        return false;
    }
};

int RunTheseTests(char** list, int n) {
    TestReporterStdout reporter;
    TestRunner runner(reporter);
    return runner.RunTestsIf(Test::GetTestList(), NULL, ListFilter(list, n), 0);
}
} // namespace UnitTest

int main(int argc, char* argv[]) {
    if (argc == 1) {
        UnitTest::RunAllTests();
    }
    else {
        UnitTest::RunTheseTests(argv + 1, argc - 1); // skip the program's name itself.
    }
}