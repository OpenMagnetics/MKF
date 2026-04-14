#include "converter_models/Dab.h"
#include "converter_models/Llc.h"
#include "advisers/CoreAdviser.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <iostream>

using namespace OpenMagnetics;

TEST_CASE("Dab_Debug_CompareWithLlc", "[debug]") {
    // LLC inputs that work with CoreAdviser
    json llcJson;
    llcJson["inputVoltage"]["nominal"] = 400.0;
    llcJson["bridgeType"] = "Half Bridge";
    llcJson["minSwitchingFrequency"] = 80000;
    llcJson["maxSwitchingFrequency"] = 120000;
    llcJson["resonantFrequency"] = 100000;
    llcJson["qualityFactor"] = 0.4;
    llcJson["inductanceRatio"] = 5.0;
    llcJson["integratedResonantInductor"] = true;
    llcJson["operatingPoints"] = json::array();
    json llcOp;
    llcOp["ambientTemperature"] = 25.0;
    llcOp["outputVoltages"] = {48.0};
    llcOp["outputCurrents"] = {10.4};
    llcOp["switchingFrequency"] = 100000;
    llcJson["operatingPoints"].push_back(llcOp);
    Llc llc(llcJson);
    auto llcInputs = llc.process();

    json llcReqJson; to_json(llcReqJson, llcInputs.get_design_requirements());
    std::cerr << "=== LLC magnetizingInductance: " << llcReqJson["magnetizingInductance"].dump() << "\n";

    // DAB inputs
    json dabJson;
    dabJson["inputVoltage"]["nominal"] = 400.0;
    dabJson["seriesInductance"] = 50e-6;
    dabJson["useLeakageInductance"] = false;
    dabJson["operatingPoints"] = json::array();
    json dabOp;
    dabOp["ambientTemperature"] = 25.0;
    dabOp["outputVoltages"] = {200.0};
    dabOp["outputCurrents"] = {10.0};
    dabOp["phaseShift"] = 25.0;
    dabOp["switchingFrequency"] = 100000;
    dabJson["operatingPoints"].push_back(dabOp);
    Dab dab(dabJson);
    auto req = dab.process_design_requirements();
    json dabReqJson; to_json(dabReqJson, req);
    std::cerr << "=== DAB magnetizingInductance: " << dabReqJson["magnetizingInductance"].dump() << "\n";
    std::cerr << "=== DAB turnsRatios: " << dabReqJson["turnsRatios"].dump() << "\n";
    std::cerr << "=== DAB insulationType: " << dabReqJson.value("insulationType", "<unset>") << "\n";
    
    CHECK(true);
}
