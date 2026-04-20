#include <catch2/catch_test_macros.hpp>
#include "json.hpp"
#include "converter_models/Flyback.h"
#include "advisers/CoreAdviser.h"
#include <iostream>

using json = nlohmann::json;
using namespace OpenMagnetics;

TEST_CASE("Flyback Core Adviser - Default Values", "[core-adviser][flyback-topology][default-values]") {
    // Default flyback values from WebFrontend/WebSharedComponents/assets/js/defaults.js
    json flybackJson;
    flybackJson["inputVoltage"]["minimum"] = 120.0;
    flybackJson["inputVoltage"]["maximum"] = 375.0;
    flybackJson["diodeVoltageDrop"] = 0.7;
    flybackJson["maximumDrainSourceVoltage"] = 600.0;
    flybackJson["maximumDutyCycle"] = 0.5;
    flybackJson["currentRippleRatio"] = 1.0;
    flybackJson["inductance"] = 200e-6;  // 200 µH
    flybackJson["dutyCycleRange"][0] = 0.5;
    flybackJson["dutyCycleRange"][1] = 0.3;
    flybackJson["deadTime"] = 0.0;
    flybackJson["efficiency"] = 0.85;
    flybackJson["numberOfOutputs"] = 1;
    
    // Output parameters
    json outputParams;
    outputParams["outputVoltages"] = json::array({12.0});
    outputParams["outputCurrents"] = json::array({5.0});
    outputParams["switchingFrequency"] = 100000.0;  // 100 kHz
    outputParams["ambientTemperature"] = 25.0;
    flybackJson["operatingPoints"] = json::array({outputParams});
    
    // Insulation settings
    flybackJson["insulationType"] = "Basic";
    flybackJson["insulation"]["cti"] = "Group II";
    flybackJson["insulation"]["pollutionDegree"] = "P2";
    flybackJson["insulation"]["overvoltageCategory"] = "OVC-III";
    flybackJson["insulation"]["altitude"]["maximum"] = 2000.0;
    flybackJson["insulation"]["mainSupplyVoltage"]["maximum"] = 400.0;
    flybackJson["insulation"]["standards"] = json::array({"IEC 60664-1"});
    
    // Create Flyback object
    OpenMagnetics::Flyback flyback(flybackJson);
    
    // Process to get Inputs
    auto inputs = flyback.process();
    
    std::cout << "[TEST] Operating points generated: " << inputs.get_operating_points().size() << std::endl;
    
    REQUIRE_FALSE(inputs.get_operating_points().empty());
    
    // Check the first operating point
    auto& firstOp = inputs.get_operating_points()[0];
    std::cout << "[TEST] Excitations per winding: " << firstOp.get_excitations_per_winding().size() << std::endl;
    
    REQUIRE_FALSE(firstOp.get_excitations_per_winding().empty());
    
    // Check frequency in first excitation
    auto& firstExc = firstOp.get_excitations_per_winding()[0];
    double frequency = firstExc.get_frequency();
    std::cout << "[TEST] Frequency in excitation: " << frequency << std::endl;
    
    // Frequency should be valid
    REQUIRE(std::isfinite(frequency));
    REQUIRE(frequency > 0);
    REQUIRE(frequency == 100000.0);  // Should be 100 kHz
    
    // Now test the core adviser
    CoreAdviser coreAdviser;
    coreAdviser.set_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
    
    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 0.4;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 0.3;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 0.3;
    
    std::cout << "[TEST] Calling get_advised_core with default flyback inputs..." << std::endl;
    
    try {
        auto results = coreAdviser.get_advised_core(inputs, weights, 5);
        
        std::cout << "[TEST] Core adviser returned " << results.size() << " results" << std::endl;
        
        // Should find at least one core
        REQUIRE_FALSE(results.empty());
        
        // Log the first few results
        for (size_t i = 0; i < std::min(results.size(), size_t(3)); ++i) {
            auto& magnetic = results[i].first.get_magnetic();
            if (magnetic.get_core().get_name()) {
                std::cout << "[TEST] Core " << i << ": " << magnetic.get_core().get_name().value() << std::endl;
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "[TEST ERROR] Exception in core adviser: " << e.what() << std::endl;
        FAIL("Core adviser threw exception: " << std::string(e.what()));
    }
}

TEST_CASE("Flyback Inputs - Frequency Validation", "[flyback-topology][frequency][validation]") {
    json flybackJson;
    flybackJson["inputVoltage"]["minimum"] = 120.0;
    flybackJson["inputVoltage"]["maximum"] = 375.0;
    flybackJson["diodeVoltageDrop"] = 0.7;
    flybackJson["maximumDrainSourceVoltage"] = 600.0;
    flybackJson["maximumDutyCycle"] = 0.5;
    flybackJson["currentRippleRatio"] = 1.0;
    flybackJson["inductance"] = 200e-6;
    
    json outputParams;
    outputParams["outputVoltages"] = json::array({12.0});
    outputParams["outputCurrents"] = json::array({5.0});
    outputParams["switchingFrequency"] = 100000.0;
    outputParams["ambientTemperature"] = 25.0;
    flybackJson["operatingPoints"] = json::array({outputParams});
    
    flybackJson["insulationType"] = "Basic";
    flybackJson["efficiency"] = 0.85;
    
    OpenMagnetics::Flyback flyback(flybackJson);
    auto inputs = flyback.process();
    
    // Serialize to JSON and check frequency
    json inputsJson;
    to_json(inputsJson, inputs);
    
    std::cout << "[TEST] Serialized inputs JSON keys in first operating point:" << std::endl;
    if (inputsJson.contains("operatingPoints") && !inputsJson["operatingPoints"].empty()) {
        auto& firstOp = inputsJson["operatingPoints"][0];
        for (auto& [key, value] : firstOp.items()) {
            std::cout << "  " << key << std::endl;
        }
        
        if (firstOp.contains("excitationsPerWinding") && !firstOp["excitationsPerWinding"].empty()) {
            auto& firstExc = firstOp["excitationsPerWinding"][0];
            std::cout << "[TEST] First excitation keys:" << std::endl;
            for (auto& [key, value] : firstExc.items()) {
                std::cout << "  " << key << ": " << value << std::endl;
            }
            
            if (firstExc.contains("frequency")) {
                double freq = firstExc["frequency"].get<double>();
                std::cout << "[TEST] Frequency in JSON: " << freq << std::endl;
                REQUIRE(std::isfinite(freq));
                REQUIRE(freq > 0);
            }
            else {
                std::cerr << "[TEST ERROR] No frequency key in excitation!" << std::endl;
                FAIL("Frequency key missing from excitation JSON");
            }
        }
    }
    
    // Deserialize back and check
    OpenMagnetics::Inputs inputs2(inputsJson);
    auto& firstOp2 = inputs2.get_operating_points()[0];
    auto& firstExc2 = firstOp2.get_excitations_per_winding()[0];
    double freq2 = firstExc2.get_frequency();
    
    std::cout << "[TEST] Frequency after deserialization: " << freq2 << std::endl;
    REQUIRE(std::isfinite(freq2));
    REQUIRE(freq2 > 0);
    REQUIRE(freq2 == 100000.0);
}