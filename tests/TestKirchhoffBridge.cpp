// KH<->MKF integration test — exercises the KirchhoffBridge across EVERY Kirchhoff topology.
//
// When MKF's converter_models were deleted, their per-topology tests went with them (the converter physics
// is now SPICE-gated inside Kirchhoff). What MKF must still guarantee is that the INTEGRATION works: for
// every topology, KH designs the converter and MKF gets a well-formed magnetic Inputs to feed the adviser.
// That is this test — one case per topology, driving the two bridge entry points:
//     design_converter_tas(topology, spec)  -> the assembled TAS
//     design_requirements_from_tas(tas)      -> the adviser's OpenMagnetics::Inputs
// plus a simulate_tas() smoke on a couple of them (proves the ngspice-via-KH path MKF now relies on).
//
// It needs no core/material DB (it stops at the magnetic Inputs), so it is fast.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "converter_models/KirchhoffBridge.h"   // design_converter_tas / design_requirements_from_tas / simulate_tas
#include "processors/Inputs.h"

using json = nlohmann::json;
using namespace OpenMagnetics;

namespace {

// A DC converter design spec (the same shape KH's own requirement fixtures use).
json dc_spec(double vin, double vout, double power, double fsw) {
    json s;
    s["designRequirements"]["efficiency"] = 1.0;
    s["designRequirements"]["inputVoltage"] = {{"minimum", vin * 0.9}, {"nominal", vin}, {"maximum", vin * 1.1}};
    s["designRequirements"]["switchingFrequency"]["nominal"] = fsw;
    s["designRequirements"]["outputs"] = json::array({ {{"name", "out"}, {"voltage", {{"nominal", vout}}}} });
    s["operatingPoints"] = json::array({ {{"inputVoltage", vin}, {"outputs", json::array({{{"power", power}}})}} });
    return s;
}

// A dual-output DC spec (isolated_buck / isolated_buck_boost: primary buck rail + isolated secondary).
json dual_output_spec(double vin, double vpri, double ppri, double vsec, double psec, double fsw) {
    json s;
    s["designRequirements"]["efficiency"] = 1.0;
    s["designRequirements"]["inputVoltage"] = {{"minimum", vin * 0.9}, {"nominal", vin}, {"maximum", vin * 1.1}};
    s["designRequirements"]["switchingFrequency"]["nominal"] = fsw;
    s["designRequirements"]["outputs"] = json::array({
        {{"name", "pri"}, {"voltage", {{"nominal", vpri}}}},
        {{"name", "sec"}, {"voltage", {{"nominal", vsec}}}}});
    s["operatingPoints"] = json::array({ {{"inputVoltage", vin},
        {"outputs", json::array({ {{"power", ppri}}, {{"power", psec}} })}} });
    return s;
}

// An AC-input PFC spec (single- or three-phase).
json ac_spec(const std::string& inputType, double vrms, double vout, double power, double fline, double fsw) {
    json s;
    s["designRequirements"]["efficiency"] = 1.0;
    s["designRequirements"]["inputType"] = inputType;
    s["designRequirements"]["inputVoltage"]["nominal"] = vrms;
    s["designRequirements"]["lineFrequency"]["nominal"] = fline;
    s["designRequirements"]["switchingFrequency"]["nominal"] = fsw;
    s["designRequirements"]["outputs"] = json::array({ {{"name", "out"}, {"voltage", {{"nominal", vout}}}} });
    s["operatingPoints"] = json::array({ {{"inputVoltage", vrms}, {"outputs", json::array({{{"power", power}}})}} });
    return s;
}

struct Case { std::string topology; json spec; };

// One representative operating point per topology (KH's own requirement fixtures).
std::vector<Case> all_topologies() {
    return {
        {"buck",                 dc_spec(12,   5,   2,    500e3)},
        {"boost",                dc_spec(5,    9,   2,    400e3)},
        {"flyback",              dc_spec(24,   6,   0.2,  250e3)},
        {"forward",              dc_spec(48,   5,   10,   200e3)},
        {"two_switch_forward",   dc_spec(48,   12,  8,    250e3)},
        {"acf",                  dc_spec(48,   3.3, 30,   250e3)},   // active-clamp forward
        {"push_pull",            dc_spec(5,    3.3, 0.35, 410e3)},
        {"sepic",                dc_spec(5,    12,  0.5,  600e3)},
        {"cuk",                  dc_spec(25,   25,  1,    100e3)},
        {"zeta",                 dc_spec(12,   5,   1,    600e3)},
        {"fsbb",                 dc_spec(24,   12,  8,    350e3)},
        {"isolated_buck",        dual_output_spec(24, 12, 1.0, 5,  0.5, 350e3)},
        {"isolated_buck_boost",  dual_output_spec(12, 5,  0.5, 12, 0.5, 250e3)},
        {"weinberg",             dc_spec(50,   150, 10,   50e3)},
        {"psfb",                 dc_spec(400,  12,  50,   100e3)},
        {"pshb",                 dc_spec(400,  12,  50,   100e3)},
        {"ahb",                  dc_spec(100,  5,   20,   200e3)},
        {"dab",                  dc_spec(800,  500, 20,   100e3)},
        {"llc",                  dc_spec(400,  12,  10,   100e3)},
        {"cllc",                 dc_spec(400,  48,  10,   200e3)},
        {"clllc",                dc_spec(400,  48,  20,   100e3)},
        {"src",                  dc_spec(400,  48,  10,   110e3)},
        {"pfc",                  ac_spec("acSinglePhase", 230, 400, 100, 50, 20e3)},
        {"vienna",               ac_spec("acThreePhase",  230, 800, 300, 50, 40e3)},
    };
}

// Assert the bridge produced a well-formed magnetic Inputs the adviser can consume.
void check_inputs(const OpenMagnetics::Inputs& inputs, const std::string& topology) {
    INFO("topology: " << topology);
    const auto& ops = inputs.get_operating_points();
    REQUIRE_FALSE(ops.empty());
    const auto& excs = ops.at(0).get_excitations_per_winding();
    REQUIRE_FALSE(excs.empty());
    // the main winding must carry a processed current (what the loss/coil filters key off)
    REQUIRE(excs.at(0).get_current().has_value());
    REQUIRE(excs.at(0).get_current()->get_processed().has_value());
}

}  // namespace

TEST_CASE("Kirchhoff bridge designs every topology into an adviser-ready Inputs", "[kirchhoff][bridge][integration]") {
    for (const auto& c : all_topologies()) {
        DYNAMIC_SECTION("topology = " << c.topology) {
            json tas = design_converter_tas(c.topology, c.spec);          // KH designs the converter
            REQUIRE(tas.contains("topology"));
            REQUIRE(tas.at("topology").contains("stages"));

            OpenMagnetics::Inputs inputs = design_requirements_from_tas(tas);            // MKF gets the magnetic Inputs
            check_inputs(inputs, c.topology);
        }
    }
}

TEST_CASE("Kirchhoff bridge simulates through KH's ngspice", "[kirchhoff][bridge][ngspice]") {
    // The ngspice-via-KH path MKF now relies on (its own ngspice was removed). Smoke a couple of topologies.
    for (const std::string topo : {"buck", "llc"}) {
        DYNAMIC_SECTION("simulate " << topo) {
            json tas = design_converter_tas(topo, topo == "buck" ? dc_spec(12, 5, 2, 500e3)
                                                                 : dc_spec(400, 12, 10, 100e3));
            json result = simulate_tas(tas, json{{"origin", "REQUIREMENTS"}});
            // success=false is acceptable ONLY if KH was built without libngspice; here it is built with it.
            INFO("simulate result: " << result.dump());
            REQUIRE(result.contains("success"));
            CHECK(result.at("success").get<bool>());
            if (result.at("success").get<bool>())
                CHECK(result.at("points").get<size_t>() > 1000);
        }
    }
}
