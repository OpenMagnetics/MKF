#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "json.hpp"
#include "converter_models/Llc.h"
#include <cmath>
#include <vector>
#include <iostream>

using json = nlohmann::json;
using namespace OpenMagnetics;

// Helper function to validate waveform data
bool validateWaveformData(const std::vector<double>& data, const std::string& name) {
    if (data.empty()) {
        std::cerr << "ERROR: " << name << " waveform data is empty" << std::endl;
        return false;
    }
    
    int nanCount = 0;
    int infCount = 0;
    for (size_t i = 0; i < data.size(); ++i) {
        if (std::isnan(data[i])) {
            nanCount++;
            if (nanCount <= 3) {
                std::cerr << "ERROR: NaN found in " << name << " at index " << i << std::endl;
            }
        }
        if (std::isinf(data[i])) {
            infCount++;
            if (infCount <= 3) {
                std::cerr << "ERROR: Inf found in " << name << " at index " << i << std::endl;
            }
        }
    }
    
    if (nanCount > 0 || infCount > 0) {
        std::cerr << "ERROR: " << name << " has " << nanCount << " NaN and " << infCount << " Inf values out of " << data.size() << " total" << std::endl;
        return false;
    }
    
    return true;
}

// ═══════════════════════════════════════════════════════════════════════════
// LLC Tests
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("LLC Analytical Waveforms", "[topology][llc][analytical]") {
    json llcJson;
    llcJson["inputVoltage"]["nominal"] = 400.0;
    llcJson["inputVoltage"]["tolerance"] = 0.1;
    llcJson["minSwitchingFrequency"] = 80000.0;
    llcJson["maxSwitchingFrequency"] = 120000.0;
    llcJson["resonantFrequency"] = 100000.0;
    llcJson["qualityFactor"] = 0.4;
    llcJson["integratedResonantInductor"] = true;
    llcJson["magnetizingInductance"] = 200e-6;
    llcJson["turnsRatio"] = 4.0;
    llcJson["ambientTemperature"] = 25.0;
    llcJson["efficiency"] = 0.97;
    llcJson["insulationType"] = "Basic";
    
    json opJson;
    opJson["outputVoltages"] = json::array({48.0});
    opJson["outputCurrents"] = json::array({500.0 / 48.0});
    opJson["switchingFrequency"] = 100000.0;
    opJson["ambientTemperature"] = 25.0;
    llcJson["operatingPoints"] = json::array({opJson});
    
    Llc llc(llcJson);
    
    auto designRequirements = llc.process_design_requirements();
    REQUIRE_FALSE(designRequirements.get_turns_ratios().empty());
    
    std::vector<double> turnsRatios;
    for (const auto& tr : designRequirements.get_turns_ratios()) {
        if (tr.get_nominal().has_value()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
    }
    
    double magnetizingInductance = 200e-6;
    auto operatingPoints = llc.process_operating_points(turnsRatios, magnetizingInductance);
    
    REQUIRE_FALSE(operatingPoints.empty());
    
    // Validate waveforms
    for (const auto& op : operatingPoints) {
        REQUIRE_FALSE(op.get_excitations_per_winding().empty());
        
        for (const auto& excitation : op.get_excitations_per_winding()) {
            // Check current waveform
            if (excitation.get_current() && excitation.get_current().value().get_waveform()) {
                auto waveform = excitation.get_current().value().get_waveform().value();
                REQUIRE(validateWaveformData(waveform.get_data(), "LLC Current"));
            }
            
            // Check voltage waveform
            if (excitation.get_voltage() && excitation.get_voltage().value().get_waveform()) {
                auto waveform = excitation.get_voltage().value().get_waveform().value();
                REQUIRE(validateWaveformData(waveform.get_data(), "LLC Voltage"));
            }
        }
    }
}

TEST_CASE("LLC Simulated Waveforms", "[topology][llc][simulation]") {
    // Skip this test if ngspice is not enabled
    #ifndef ENABLE_NGSPICE
    SKIP("ngspice is not enabled, skipping simulation test");
    #endif
    
    json llcJson;
    llcJson["inputVoltage"]["nominal"] = 400.0;
    llcJson["inputVoltage"]["tolerance"] = 0.1;
    llcJson["minSwitchingFrequency"] = 80000.0;
    llcJson["maxSwitchingFrequency"] = 120000.0;
    llcJson["resonantFrequency"] = 100000.0;
    llcJson["qualityFactor"] = 0.4;
    llcJson["integratedResonantInductor"] = true;
    llcJson["magnetizingInductance"] = 200e-6;
    llcJson["turnsRatio"] = 4.0;
    llcJson["ambientTemperature"] = 25.0;
    llcJson["efficiency"] = 0.97;
    llcJson["insulationType"] = "Basic";
    
    json opJson;
    opJson["outputVoltages"] = json::array({48.0});
    opJson["outputCurrents"] = json::array({500.0 / 48.0});
    opJson["switchingFrequency"] = 100000.0;
    opJson["ambientTemperature"] = 25.0;
    llcJson["operatingPoints"] = json::array({opJson});
    
    Llc llc(llcJson);
    
    auto designRequirements = llc.process_design_requirements();
    std::vector<double> turnsRatios;
    for (const auto& tr : designRequirements.get_turns_ratios()) {
        if (tr.get_nominal().has_value()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
    }
    
    double magnetizingInductance = 200e-6;
    
    auto operatingPoints = llc.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
    auto converterWaveforms = llc.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
    
    REQUIRE_FALSE(operatingPoints.empty());
    REQUIRE_FALSE(converterWaveforms.empty());
    
    // Validate simulated waveforms
    for (const auto& op : operatingPoints) {
        for (const auto& excitation : op.get_excitations_per_winding()) {
            if (excitation.get_current() && excitation.get_current().value().get_waveform()) {
                auto waveform = excitation.get_current().value().get_waveform().value();
                REQUIRE(validateWaveformData(waveform.get_data(), "LLC Simulated Current"));
            }
        }
    }
}
