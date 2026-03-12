#include "physical_models/MagneticEnergy.h"
#include "TestingUtils.h"
#include "support/Utils.h"
#include "json.hpp"
#include "advisers/MagneticAdviser.h"
#include "advisers/CoreAdviser.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <typeinfo>
#include <vector>
#include <chrono>

using namespace MAS;
using namespace OpenMagnetics;

using json = nlohmann::json;

namespace {
    TEST_CASE("Test_Magnetic_Saturation_Current", "[constructive-model][magnetic][smoke-test]") {

        std::vector<int64_t> numberTurns = {18};
        std::vector<int64_t> numberParallels = {1};
        std::string shapeName = "PQ 65/44";

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName);

        // coil.wind({0, 1, 0}, 1);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto saturationCurrentAt20 = magnetic.calculate_saturation_current(20);
        auto saturationCurrentAt100 = magnetic.calculate_saturation_current(100);
        REQUIRE(saturationCurrentAt100 < saturationCurrentAt20);
    }

    void run_adviser_with_mode(const std::string& modeName, CoreAdviser::CoreAdviserModes mode) {
        std::cout << "\n\n========================================" << std::endl;
        std::cout << "Testing with mode: " << modeName << std::endl;
        std::cout << "========================================" << std::endl;
        
        auto json_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "magnetic_adviser_inductor_bug.json");
        std::ifstream json_file(json_path);
        json masJson = json::parse(json_file);

        OpenMagnetics::Inputs inputs(masJson["inputs"]);
        
        MagneticAdviser magneticAdviser;
        magneticAdviser.set_core_mode(mode);
        
        // Measure runtime
        auto start = std::chrono::high_resolution_clock::now();
        auto masMagnetics = magneticAdviser.get_advised_magnetic(inputs, 5);
        auto end = std::chrono::high_resolution_clock::now();
        auto runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        REQUIRE(masMagnetics.size() > 0);
        
        std::cout << "\n=== Magnetic Adviser Recommendations (" << modeName << ") ===" << std::endl;
        std::cout << "Number of recommendations: " << masMagnetics.size() << std::endl;
        std::cout << "Runtime: " << runtime << " ms" << std::endl;
        
        for (size_t i = 0; i < masMagnetics.size(); ++i) {
            auto& mas = masMagnetics[i].first;
            double scoring = masMagnetics[i].second;
            
            std::cout << "\n--- Recommendation " << (i + 1) << " ---" << std::endl;
            std::cout << "Scoring: " << scoring << std::endl;
            
            // Get core info
            auto coreName = mas.get_magnetic().get_core().get_name();
            if (coreName.has_value()) {
                std::cout << "Core: " << coreName.value() << std::endl;
            }
            
            // Get coil info
            auto& coil = mas.get_magnetic().get_coil();
            auto windings = coil.get_functional_description();
            std::cout << "Number of windings: " << windings.size() << std::endl;
            for (size_t w = 0; w < windings.size(); ++w) {
                std::cout << "  Winding " << (w + 1) << ": " 
                          << windings[w].get_number_turns() << " turns" << std::endl;
            }
            
            // Get gap info
            auto& core = mas.get_magnetic().get_core();
            if (!core.get_functional_description().get_gapping().empty()) {
                double gapLength = core.get_functional_description().get_gapping()[0].get_length();
                std::cout << "Gap: " << gapLength * 1000 << " mm" << std::endl;
            }
            
            // Get manufacturer info
            auto mfgInfo = mas.get_magnetic().get_manufacturer_info();
            if (mfgInfo.has_value()) {
                auto ref = mfgInfo.value().get_reference();
                if (ref.has_value()) {
                    std::cout << "Reference: " << ref.value() << std::endl;
                }
            }
            
            // Check outputs
            if (!mas.get_outputs().empty()) {
                auto outputs = mas.get_outputs()[0];
                auto inductance = outputs.get_inductance();
                if (inductance.has_value()) {
                    auto magInd = inductance.value().get_magnetizing_inductance().get_magnetizing_inductance().get_nominal();
                    if (magInd.has_value()) {
                        std::cout << "Magnetizing Inductance: " 
                                  << magInd.value() 
                                  << " H" << std::endl;
                    }
                }
                
                // Get saturation info from excitation
                if (!inputs.get_operating_points().empty()) {
                    auto& op = inputs.get_operating_points()[0];
                    auto excitation = OpenMagnetics::Inputs::get_primary_excitation(op);
                    auto bFieldOpt = excitation.get_magnetic_flux_density();
                    if (bFieldOpt.has_value()) {
                        auto bField = bFieldOpt.value();
                        auto processedOpt = bField.get_processed();
                        if (processedOpt.has_value() && processedOpt->get_peak().has_value()) {
                            double bPeak = processedOpt->get_peak().value();
                            auto& core = mas.get_mutable_magnetic().get_mutable_core();
                            double temperature = op.get_conditions().get_ambient_temperature();
                            double bSat = core.get_magnetic_flux_density_saturation(temperature, true);
                            double saturationRatio = (bSat > 0) ? (bPeak / bSat * 100.0) : 0.0;
                            std::cout << "Bpeak: " << bPeak << " T, Bsat: " << bSat 
                                      << " T, Saturation: " << saturationRatio << "%" << std::endl;
                        }
                    }
                }
            }
        }
    }

    TEST_CASE("Test_Magnetic_Inductor_Bug_Standard", "[magnetic][adviser][inductor][bug]") {
        clear_databases();
        settings.reset();
        // Test with STANDARD_CORES (default)
        run_adviser_with_mode("STANDARD_CORES", CoreAdviser::CoreAdviserModes::STANDARD_CORES);
        settings.reset();
    }
    
    TEST_CASE("Test_Magnetic_Inductor_Bug_Available", "[magnetic][adviser][inductor][bug][available]") {
        clear_databases();
        settings.reset();
        // Test with AVAILABLE_CORES
        run_adviser_with_mode("AVAILABLE_CORES", CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
        settings.reset();
    }

    TEST_CASE("Test_Magnetic_Inductor_Bug_Standard_GoldenSection", "[magnetic][adviser][inductor][bug][golden]") {
        clear_databases();
        settings.reset();
        settings.set_gapping_strategy(GappingOptimizationStrategy::GOLDEN_SECTION);
        // Test with STANDARD_CORES using Golden Section optimization
        run_adviser_with_mode("STANDARD_CORES_GOLDEN", CoreAdviser::CoreAdviserModes::STANDARD_CORES);
        settings.reset();
    }

}  // namespace
