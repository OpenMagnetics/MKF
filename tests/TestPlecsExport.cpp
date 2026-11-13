// PLECS exporter output syntax.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "processors/CircuitSimulatorInterface.h"
#include "constructive_models/Magnetic.h"
#include "TestingUtils.h"

using namespace OpenMagnetics;

// ABT #2369 (web bug report #195): a Connection's Points were written as one bracket per point,
// "Points [405, 210]; [440, 210]", which PLECS refuses ("syntax error before ;"). PLECS reads them
// as one matrix with a row per point: "Points [405, 210; 440, 210]".
TEST_CASE("PLECS export writes connection Points as one matrix", "[circuit][export][plecs][abt-2369][smoke-test]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 39", OpenMagneticsTesting::get_residual_gap(),
                                                             std::vector<int64_t>{10, 5}, 1, "3C97");
    auto coil = magnetic.get_coil();
    coil.wind();
    magnetic.set_coil(coil);

    CircuitSimulatorExporter exporter(CircuitSimulatorExporterModels::PLECS);
    std::string model = exporter.export_magnetic_as_subcircuit(magnetic, 100000.0, 25.0);
    REQUIRE(!model.empty());

    CHECK(model.find("]; [") == std::string::npos);

    const std::regex matrix(R"(^\s*Points\s+\[-?\d+, -?\d+(; -?\d+, -?\d+)+\]\s*$)");
    size_t pointsLines = 0;
    std::istringstream lines(model);
    std::string line;
    while (std::getline(lines, line)) {
        if (line.find("Points") == std::string::npos) {
            continue;
        }
        ++pointsLines;
        INFO(line);
        CHECK(std::regex_match(line, matrix));
    }
    REQUIRE(pointsLines > 0);
}
