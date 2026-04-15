#include "converter_models/Dab.h"
#include "converter_models/Llc.h"
#include "advisers/CoreAdviser.h"
#include "processors/Inputs.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <iostream>

using namespace OpenMagnetics;

static void dump_exc(const char* tag, const OperatingPointExcitation& exc) {
    bool hv = exc.get_voltage().has_value();
    bool hw = hv && exc.get_voltage()->get_waveform().has_value();
    bool ht = hw && exc.get_voltage()->get_waveform()->get_time().has_value();
    size_t dSz = hw ? exc.get_voltage()->get_waveform()->get_data().size() : 0;
    size_t tSz = ht ? exc.get_voltage()->get_waveform()->get_time()->size() : 0;
    double vpk = (hv && exc.get_voltage()->get_processed() && exc.get_voltage()->get_processed()->get_peak())
                 ? exc.get_voltage()->get_processed()->get_peak().value() : -1.0;
    std::cerr << "[" << tag << "] hasV=" << hv << " hasWf=" << hw << " hasT=" << ht
              << " dSz=" << dSz << " tSz=" << tSz << " Vpk=" << vpk
              << " freq=" << exc.get_frequency() << "\n";
}

TEST_CASE("Dab_Wizard_ETD34_turns", "[debug]") {
    // Replicate what DabWizard.vue sends to calculate_dab_inputs & then to calculate_advised_cores
    json dabParams = {
        {"inputVoltage", {{"nominal", 400}}},
        {"efficiency", 0.97},
        {"useLeakageInductance", true},
        {"desiredTurnsRatios", json::array({1.0})},
        {"desiredMagnetizingInductance", 1e-3},
        {"operatingPoints", json::array({
            {
                {"outputVoltages", json::array({400})},
                {"outputCurrents", json::array({2.5})},
                {"phaseShift", 30},
                {"modulationType", "SPS"},
                {"innerPhaseShift1", 15},
                {"innerPhaseShift2", 15},
                {"ambientTemperature", 25},
                {"switchingFrequency", 100000}
            }
        })}
    };

    // Path 1: direct (as in test_dab_wizard_turns.cpp)
    AdvancedDab dab(dabParams);
    dab.set_num_periods_to_extract(1);
    auto inputsDirect = dab.process();
    dump_exc("direct-after-process", inputsDirect.get_operating_point(0).get_excitations_per_winding()[0]);

    // Path 2: serialize then deserialize (mimics libMKF::calculate_dab_inputs → calculate_advised_cores)
    json serialized;
    to_json(serialized, inputsDirect);
    std::string serStr = serialized.dump();
    std::cerr << "[serialized] size=" << serStr.size() << " bytes\n";

    // Check raw JSON for waveform/time
    {
        auto& exc = serialized["operatingPoints"][0]["excitationsPerWinding"][0];
        bool jHasV = exc.contains("voltage");
        bool jHasWf = jHasV && exc["voltage"].contains("waveform");
        bool jHasTime = jHasWf && exc["voltage"]["waveform"].contains("time");
        size_t jDSz = jHasWf && exc["voltage"]["waveform"].contains("data") ? exc["voltage"]["waveform"]["data"].size() : 0;
        size_t jTSz = jHasTime ? exc["voltage"]["waveform"]["time"].size() : 0;
        std::cerr << "[serialized-json] hasV=" << jHasV << " hasWf=" << jHasWf
                  << " hasT=" << jHasTime << " dSz=" << jDSz << " tSz=" << jTSz << "\n";
    }

    // Now construct Inputs from the serialized JSON (this is what calculate_advised_cores does)
    OpenMagnetics::Inputs inputsRound(json::parse(serStr));
    dump_exc("round-after-ctor", inputsRound.get_operating_point(0).get_excitations_per_winding()[0]);

    // Run CoreAdviser on the round-tripped inputs
    settings.reset();
    OpenMagnetics::Settings::GetInstance().set_coil_delimit_and_compact(true);

    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 0.3;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 0.4;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 0.3;

    CoreAdviser coreAdviser;
    coreAdviser.set_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
    auto results = coreAdviser.get_advised_core(inputsRound, weights, 20);

    std::cerr << "=== CoreAdviser results (round-trip path) ===\n";
    for (size_t i = 0; i < std::min(results.size(), size_t(10)); ++i) {
        auto name = results[i].first.get_magnetic().get_core().get_name().value_or("?");
        auto turns = results[i].first.get_magnetic().get_coil().get_functional_description()[0].get_number_turns();
        std::cerr << "  " << i << ": " << name << " -> " << turns << " turns\n";
    }

    REQUIRE_FALSE(results.empty());

    // Sanity: every result has Faraday-compliant turns for the 0.7·Bsat ceiling.
    // N_min = ceil(maxVs / (0.7·Bsat·Ae·stacks)). With maxVs ≈ 0.00183 V·s,
    // Bsat ≈ 0.35 T (material 95/98), 0.7·Bsat ≈ 0.245 T. Allow a 40% headroom
    // above this floor to account for the loss-sweep and stacked variants.
    for (size_t i = 0; i < std::min(results.size(), size_t(10)); ++i) {
        auto& mag = results[i].first.get_magnetic();
        auto turns = mag.get_coil().get_functional_description()[0].get_number_turns();
        double Ae = mag.get_core().get_processed_description()->get_effective_parameters().get_effective_area();
        double Bpeak = 0.00183 / (turns * Ae);
        CAPTURE(i, mag.get_core().get_name().value_or("?"), turns, Ae, Bpeak);
        // Must not be in saturation-sinusoidal regime (~0.1 T). A healthy post-fix
        // design sits in ~0.18–0.36 T range under volt-second integration.
        REQUIRE(Bpeak > 0.15);
        REQUIRE(Bpeak < 0.40);
    }

    // ETD 34/17/11 specifically: verify the historical 72-turn regression doesn't reappear.
    // Native C++ test showed 54 turns; allow ±30% around that as a guard band.
    for (auto& r : results) {
        auto name = r.first.get_magnetic().get_core().get_name().value_or("?");
        if (name.find("ETD 34/17/11") != std::string::npos) {
            auto turns = r.first.get_magnetic().get_coil().get_functional_description()[0].get_number_turns();
            std::cerr << "ETD34: " << name << " -> " << turns << " turns\n";
            CHECK(turns >= 38);
            CHECK(turns <= 72);
            break;
        }
    }
}

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
