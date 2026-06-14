#include "converter_models/DifferentialModeChoke.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <string>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

// Captured from the WebFrontend DMC wizard "Help me with the design" flow
// (2026-06-14): singlePhaseBalanced, operatingCurrent = 3 A, no peakCurrent.
// propose_design() -> process_operating_points() must size the choke without
// requiring an explicit peakCurrent (it derives it from the operating current).
const std::string DMC_HELP_MODE_INPUT = R"JSON(
{"configuration":"singlePhaseBalanced","inputVoltage":{"nominal":230,"tolerance":0.1},"operatingCurrent":3,"lineFrequency":50,"switchingFrequency":100000,"ambientTemperature":25,"filterCapacitance":1e-06,"minimumImpedance":[{"frequency":150000,"impedance":{"magnitude":500}}],"numberOfPeriods":2}
)JSON";

} // namespace

TEST_CASE("Test_Dmc_HelpMode_DerivesPeakCurrent", "[converter-model][dmc-topology][frontend-repro]") {
    // "Help me with the design": operatingCurrent given, peakCurrent omitted.
    // The model must derive the peak current and produce a design, not throw.
    OpenMagnetics::DifferentialModeChoke dmc(json::parse(DMC_HELP_MODE_INPUT));
    REQUIRE_NOTHROW(dmc.propose_design());
}
