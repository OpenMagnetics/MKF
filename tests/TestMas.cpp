#include <source_location>
#include "constructive_models/Mas.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/StrayCapacitance.h"
#include "processors/MagneticSimulator.h"
#include "processors/Sweeper.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <typeinfo>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

using json = nlohmann::json;

namespace {
    TEST_CASE("Test_Expand_Magnetic", "[constructive-model][mas][smoke-test]") {

        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "example_basic.json");
        auto mas = OpenMagneticsTesting::mas_loader(path.string());

        REQUIRE(!mas.get_magnetic().get_core().get_processed_description());

        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        auto inputs = OpenMagnetics::inputs_autocomplete(mas.get_inputs(), magnetic);
        REQUIRE(magnetic.get_core().get_processed_description());

    }

    TEST_CASE("Test_Load_Mas_With_Microsign_In_Wire_Name", "[constructive-model][mas][encoding][smoke-test]") {
        // This test verifies that files containing the micro sign (µ) character in wire names
        // are loaded correctly across platforms (Windows/Linux encoding differences)
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_microsign_encoding.json");
        
        std::ifstream json_file(path);
        REQUIRE(json_file.is_open());
        
        auto jsonData = json::parse(json_file);
        json_file.close();
        
        // Verify the JSON was parsed successfully
        REQUIRE(jsonData.contains("magnetic"));
        REQUIRE(jsonData["magnetic"].contains("coil"));
        REQUIRE(jsonData["magnetic"]["coil"].contains("functionalDescription"));
        
        // Get the wire name from the first winding
        auto& windings = jsonData["magnetic"]["coil"]["functionalDescription"];
        REQUIRE(windings.size() > 0);
        REQUIRE(windings[0].contains("wire"));
        REQUIRE(windings[0]["wire"].contains("name"));
        
        std::string wireName = windings[0]["wire"]["name"];
        
        // The wire name should contain "Planar" and should be non-empty
        REQUIRE(!wireName.empty());
        REQUIRE(wireName.find("Planar") != std::string::npos);
        
        // Verify we can construct a Mas object from the JSON
        OpenMagnetics::Mas mas(jsonData);
        
        // Get the wire name from the loaded Mas object
        // get_wire() returns std::variant<Wire, std::string>
        const auto& wireVariant = mas.get_magnetic().get_coil().get_functional_description()[0].get_wire();
        REQUIRE(std::holds_alternative<OpenMagnetics::Wire>(wireVariant));
        const auto& wire = std::get<OpenMagnetics::Wire>(wireVariant);
        auto loadedWireName = wire.get_name();
        REQUIRE(loadedWireName.has_value());
        REQUIRE(loadedWireName.value().find("Planar") != std::string::npos);
        
        // Verify the magnetic can be processed
        auto magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic());
        REQUIRE(magnetic.get_core().get_processed_description());
    }

// =============================================================================
// ABT #278: real-geometry physics battery over every MAS example.
// Each example is wound with REAL winding geometry (per-edge lead rows, toroidal
// angular corridors — ABT #229/#187) and then pushed through the full physics
// chain: magnetizing inductance, leakage inductance (multi-winding), stray
// capacitance matrix, core + winding losses (MagneticSimulator, per operating
// point) and the impedance-over-frequency sweep. Every produced number must be
// finite and physical. Not a smoke test: this is the example-battery gate.
// =============================================================================
TEST_CASE("Test_All_Examples_Real_Geometry_Physics", "[constructive-model][mas][example-battery]") {
    namespace fs = std::filesystem;
    auto examplesDir = fs::path{std::source_location::current().file_name()}.parent_path().append("..").append("MAS").append("examples");
    std::vector<fs::path> files;
    for (auto& entry : fs::directory_iterator(examplesDir)) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());
    REQUIRE(files.size() >= 24);

    for (auto& file : files) {
        std::string name = file.filename().string();
        DYNAMIC_SECTION(name) {
            settings.reset();
            auto mas = OpenMagneticsTesting::mas_loader(file.string());

            // Expand through the canonical path (resolves shape/bobbin, processes the core, winds
            // the coil) with REAL winding geometry ON, so the wind inside builds the lead rows /
            // toroidal corridors (ABT #229/#187).
            settings.set_coil_use_real_winding_geometry(true);

            // ABT #492 owner ruling: planar wires are PCBs — the real-winding connection model
            // (leads, markers, dragbacks) is for WOUND magnetics only, and real winding for planar
            // has not been started, so this combination must THROW (no via model, no fallback, no
            // silent skip). The planar example (09_planar_xfmr_er2510_3c94) therefore asserts the
            // loud gate here; the physics battery below applies to the wound examples.
            bool anyPlanarWire = false;
            for (const auto& winding : mas.get_magnetic().get_coil().get_functional_description()) {
                if (std::holds_alternative<OpenMagnetics::Wire>(winding.get_wire())
                    && std::get<OpenMagnetics::Wire>(winding.get_wire()).get_type() == WireType::PLANAR) {
                    anyPlanarWire = true;
                }
            }
            if (anyPlanarWire) {
                // The planar example ships fully wound, so autocomplete does not re-wind it; the
                // gate fires at the first real-winding machinery a planar coil can reach: a
                // re-wind, and the connection-resistance path the loss chain uses.
                OpenMagnetics::Magnetic planarMagnetic;
                REQUIRE_NOTHROW(planarMagnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic()));
                REQUIRE_THROWS_WITH(planarMagnetic.get_mutable_coil().wind(),
                                    Catch::Matchers::ContainsSubstring("not implemented for planar"));
                REQUIRE_THROWS_WITH(WindingOhmicLosses::calculate_connection_resistance_per_winding_per_parallel(
                                        planarMagnetic.get_coil(), 25.0),
                                    Catch::Matchers::ContainsSubstring("not implemented for planar"));
                settings.reset();
                continue;
            }

            OpenMagnetics::Magnetic magnetic;
            REQUIRE_NOTHROW(magnetic = OpenMagnetics::magnetic_autocomplete(mas.get_magnetic()));
            auto& core = magnetic.get_mutable_core();
            auto& coil = magnetic.get_mutable_coil();
            REQUIRE(coil.get_turns_description());

            // Magnetizing inductance.
            double magnetizingInductance = 0;
            REQUIRE_NOTHROW(magnetizingInductance = MagnetizingInductance()
                .calculate_inductance_from_number_turns_and_gapping(core, coil)
                .get_magnetizing_inductance().get_nominal().value());
            UNSCOPED_INFO("L = " << magnetizingInductance);
            CHECK(std::isfinite(magnetizingInductance));
            CHECK(magnetizingInductance > 0);

            // Leakage inductance (transformers only).
            if (coil.get_functional_description().size() >= 2) {
                double frequency = Defaults().measurementFrequency;
                if (!mas.get_inputs().get_operating_points().empty()) {
                    frequency = mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_frequency();
                }
                double leakageInductance = 0;
                REQUIRE_NOTHROW(leakageInductance = LeakageInductance()
                    .calculate_leakage_inductance(magnetic, frequency)
                    .get_leakage_inductance_per_winding()[0].get_nominal().value());
                UNSCOPED_INFO("Llk = " << leakageInductance);
                CHECK(std::isfinite(leakageInductance));
                CHECK(leakageInductance > 0);
                CHECK(leakageInductance < magnetizingInductance);
            }

            // Stray capacitance matrix.
            {
                StrayCapacitanceOutput capacitanceOutput;
                REQUIRE_NOTHROW(capacitanceOutput = StrayCapacitance().calculate_capacitance(coil));
                auto maxwellCapacitanceMatrix = capacitanceOutput.get_maxwell_capacitance_matrix().value();
                CHECK(!maxwellCapacitanceMatrix.empty());
                for (auto& matrixAtFrequency : maxwellCapacitanceMatrix) {
                    for (auto& [firstKey, row] : matrixAtFrequency.get_magnitude()) {
                        for (auto& [secondKey, capacitanceWithTolerance] : row) {
                            auto capacitance = OpenMagnetics::resolve_dimensional_values(capacitanceWithTolerance);
                            UNSCOPED_INFO("C[" << firstKey << "][" << secondKey << "] = " << capacitance);
                            CHECK(std::isfinite(capacitance));
                        }
                    }
                }
            }

            // Core + winding losses per operating point (full simulate, fast mode).
            {
                OpenMagnetics::Mas simulated;
                REQUIRE_NOTHROW(simulated = MagneticSimulator().simulate(mas.get_inputs(), magnetic, true));
                REQUIRE(!simulated.get_outputs().empty());
                for (size_t operatingPointIndex = 0; operatingPointIndex < simulated.get_outputs().size(); ++operatingPointIndex) {
                    auto& output = simulated.get_outputs()[operatingPointIndex];
                    REQUIRE(output.get_core_losses());
                    double coreLosses = output.get_core_losses()->get_core_losses();
                    UNSCOPED_INFO("op " << operatingPointIndex << " coreLosses = " << coreLosses);
                    CHECK(std::isfinite(coreLosses));
                    CHECK(coreLosses >= 0);
                    REQUIRE(output.get_winding_losses());
                    double windingLosses = output.get_winding_losses()->get_winding_losses();
                    UNSCOPED_INFO("op " << operatingPointIndex << " windingLosses = " << windingLosses);
                    CHECK(std::isfinite(windingLosses));
                    CHECK(windingLosses > 0);
                }
            }

            // Impedance sweep (the S-parameter-style frequency response). A material without
            // complex permeability AND without frequency-dependent initial permeability cannot
            // honestly produce an impedance curve — ComplexPermeability throws the loud
            // MaterialDataMissingException per the no-fallback rule. That exact data gap (tracked
            // as a MAS data task under ABT #278: XFlux 60, TP4A, PC95) is reported and tolerated
            // here; ANY other failure still fails the battery.
            {
                Curve2D impedanceCurve;
                bool materialLacksComplexPermeability = false;
                try {
                    impedanceCurve = Sweeper::sweep_impedance_over_frequency(magnetic, 1e3, 1e7, 15);
                }
                catch (const MaterialDataMissingException& e) {
                    materialLacksComplexPermeability = true;
                    WARN(name << ": impedance sweep skipped — " << e.what()
                              << " (MAS data gap, tracked in ABT #278)");
                }
                if (!materialLacksComplexPermeability) {
                    REQUIRE(!impedanceCurve.get_y_points().empty());
                    for (size_t pointIndex = 0; pointIndex < impedanceCurve.get_y_points().size(); ++pointIndex) {
                        double impedance = impedanceCurve.get_y_points()[pointIndex];
                        UNSCOPED_INFO("Z[" << impedanceCurve.get_x_points()[pointIndex] << " Hz] = " << impedance);
                        CHECK(std::isfinite(impedance));
                        CHECK(impedance > 0);
                    }
                }
            }

            settings.reset();
        }
    }
}

}  // namespace
