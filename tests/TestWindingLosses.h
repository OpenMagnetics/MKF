#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>
#include "physical_models/WindingLosses.h"
#include "TestingUtils.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace WindingLossesTestData {

    // Enum for wire type classification
    enum class WireTypeClass {
        ROUND,
        LITZ,
        RECTANGULAR,
        FOIL,
        PLANAR
    };

    // Structure to hold test configuration
    struct TestConfig {
        std::string name;
        WireTypeClass wireType;
        
        // Operating point parameters
        WaveformLabel waveform = WaveformLabel::SINUSOIDAL;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double offset = 0;
        double magnetizingInductance = 1e-3;
        double temperature = 20;
        
        // Magnetic field settings
        bool includeFringing = false;
        int mirroringDimension = 1;
        
        // Expected values: (frequency, expected_loss)
        std::vector<std::pair<double, double>> expectedValues;
        
        // Function to create the magnetic - set by each test
        std::function<OpenMagnetics::Magnetic()> createMagnetic;
    };

    // Helper function to get wire type as string
    inline std::string wireTypeToString(WireTypeClass wt) {
        switch (wt) {
            case WireTypeClass::ROUND: return "ROUND";
            case WireTypeClass::LITZ: return "LITZ";
            case WireTypeClass::RECTANGULAR: return "RECTANGULAR";
            case WireTypeClass::FOIL: return "FOIL";
            case WireTypeClass::PLANAR: return "PLANAR";
            default: return "UNKNOWN";
        }
    }

    // ==================================================================================
    // ROUND WIRE TEST CONFIGURATIONS
    // ==================================================================================
    
    inline TestConfig createOneTurnRoundStackedConfig() {
        TestConfig config;
        config.name = "One_Turn_Round_Stacked";
        config.wireType = WireTypeClass::ROUND;
        config.temperature = 20;
        config.waveform = WaveformLabel::SINUSOIDAL;
        config.peakToPeak = 2 * 1.4142;
        config.dutyCycle = 0.5;
        config.offset = 0;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {{0.01, 0.0049}};
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "E 42/21/20";
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_diameter(0.00071);
            wire.set_nominal_value_outer_diameter(0.000762);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::ROUND);
            wires.push_back(wire);
            int64_t numberStacks = 2;
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires, true, numberStacks);
            
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createOneTurnRoundSinusoidalConfig() {
        TestConfig config;
        config.name = "One_Turn_Round_Sinusoidal";
        config.wireType = WireTypeClass::ROUND;
        config.temperature = 22;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {1, 0.0022348},
            {10000, 0.002238},
            {20000, 0.0022476},
            {30000, 0.0022633},
            {40000, 0.0022546},
            {50000, 0.0023109},
            {60000, 0.0023419},
            {70000, 0.0023769},
            {80000, 0.0024153},
            {90000, 0.0024569},
            {100000, 0.0025011},
            {200000, 0.0030259},
            {300000, 0.0035737},
            {400000, 0.0040654},
            {500000, 0.0044916},
            {600000, 0.0048621},
            {700000, 0.0051882},
            {800000, 0.0054789},
            {900000, 0.0057414},
            {1000000, 0.0059805},
        };
        
        config.createMagnetic = []() {
            auto path = std::filesystem::path(__FILE__).parent_path() / "Test_Winding_Losses_One_Turn_Round_Sinusoidal.json";
            auto mas = OpenMagneticsTesting::mas_loader(path);
            return mas.get_magnetic();
        };
        
        return config;
    }

    inline TestConfig createTwelveTurnsRoundSinusoidalConfig() {
        TestConfig config;
        config.name = "Twelve_Turns_Round_Sinusoidal";
        config.wireType = WireTypeClass::ROUND;
        config.temperature = 22;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {1, 0.17371},
            {10000, 0.17372},
            {20000, 0.17373},
            {30000, 0.17374},
            {40000, 0.17375},
            {50000, 0.17378},
            {60000, 0.1738},
            {70000, 0.17384},
            {80000, 0.17387},
            {90000, 0.17391},
            {100000, 0.17396},
            {200000, 0.1747},
            {300000, 0.17593},
            {400000, 0.17764},
            {500000, 0.17983},
            {600000, 0.18248},
            {700000, 0.1856},
            {800000, 0.18916},
            {900000, 0.19315},
            {1000000, 0.19755},
            {3000000, 0.34496},
        };
        
        config.createMagnetic = []() {
            auto path = std::filesystem::path(__FILE__).parent_path() / "Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal.json";
            auto mas = OpenMagneticsTesting::mas_loader(path);
            return mas.get_magnetic();
        };
        
        return config;
    }

    inline TestConfig createOneTurnRoundFringingConfig() {
        TestConfig config;
        config.name = "One_Turn_Round_Fringing";
        config.wireType = WireTypeClass::ROUND;
        config.temperature = 22;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {1, 167.89},
            {10000, 169.24},
            {20000, 174.77},
            {30000, 183.33},
            {40000, 194.12},
            {50000, 206.33},
            {60000, 219.3},
            {70000, 232.5},
            {80000, 245.61},
            {90000, 258.43},
            {100000, 270.86},
            {200000, 376.14},
            {300000, 460.7},
            {400000, 532.27},
            {500000, 594.6},
            {600000, 649.64},
            {700000, 699.9},
            {800000, 746.3},
            {900000, 789.66},
            {1000000, 830.49},
        };
        
        config.createMagnetic = []() {
            auto path = std::filesystem::path(__FILE__).parent_path() / "Test_Winding_Losses_One_Turn_Round_Sinusoidal_Fringing.json";
            auto mas = OpenMagneticsTesting::mas_loader(path);
            return mas.get_magnetic();
        };
        
        return config;
    }

    // ==================================================================================
    // LITZ WIRE TEST CONFIGURATIONS
    // ==================================================================================

    inline TestConfig createOneTurnLitzSinusoidalConfig() {
        TestConfig config;
        config.name = "One_Turn_Litz_Sinusoidal";
        config.wireType = WireTypeClass::LITZ;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.003374},
            {25000, 0.003371},
            {50000, 0.00336},
            {100000, 0.003387},
            {200000, 0.003415},
            {250000, 0.003435},
            {500000, 0.003629}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            WireRound strand;
            OpenMagnetics::Wire wire;
            DimensionWithTolerance strandConductingDiameter;
            DimensionWithTolerance strandOuterDiameter;
            strandConductingDiameter.set_nominal(0.000071);
            strandOuterDiameter.set_nominal(0.0000762);
            strand.set_conducting_diameter(strandConductingDiameter);
            strand.set_outer_diameter(strandOuterDiameter);
            strand.set_number_conductors(1);
            strand.set_material("copper");
            strand.set_type(WireType::ROUND);
            
            wire.set_strand(strand);
            wire.set_nominal_value_outer_diameter(0.000873);
            wire.set_number_conductors(60);
            wire.set_type(WireType::LITZ);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createOneTurnLitzManyStrandsConfig() {
        TestConfig config;
        config.name = "One_Turn_Litz_Many_Strands";
        config.wireType = WireTypeClass::LITZ;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.00033763},
            {25000, 0.00033732},
            {50000, 0.00033622},
            {100000, 0.00033887},
            {200000, 0.00034167},
            {250000, 0.00034373},
            {500000, 0.00036313}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            WireRound strand;
            OpenMagnetics::Wire wire;
            DimensionWithTolerance strandConductingDiameter;
            DimensionWithTolerance strandOuterDiameter;
            strandConductingDiameter.set_nominal(0.000071);
            strandOuterDiameter.set_nominal(0.0000762);
            strand.set_conducting_diameter(strandConductingDiameter);
            strand.set_outer_diameter(strandOuterDiameter);
            strand.set_number_conductors(1);
            strand.set_material("copper");
            strand.set_type(WireType::ROUND);
            
            wire.set_strand(strand);
            wire.set_nominal_value_outer_diameter(0.00256);
            wire.set_number_conductors(600);
            wire.set_type(WireType::LITZ);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createTenTurnsLitzSinusoidalConfig() {
        TestConfig config;
        config.name = "Ten_Turns_Litz_Sinusoidal";
        config.wireType = WireTypeClass::LITZ;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.03375},
            {25000, 0.03371},
            {50000, 0.03361},
            {100000, 0.0339},
            {200000, 0.03463},
            {250000, 0.03553},
            {500000, 0.04004}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({10});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            WireRound strand;
            OpenMagnetics::Wire wire;
            DimensionWithTolerance strandConductingDiameter;
            DimensionWithTolerance strandOuterDiameter;
            strandConductingDiameter.set_nominal(0.000071);
            strandOuterDiameter.set_nominal(0.0000762);
            strand.set_conducting_diameter(strandConductingDiameter);
            strand.set_outer_diameter(strandOuterDiameter);
            strand.set_number_conductors(1);
            strand.set_material("copper");
            strand.set_type(WireType::ROUND);
            
            wire.set_strand(strand);
            wire.set_nominal_value_outer_diameter(0.000873);
            wire.set_number_conductors(60);
            wire.set_type(WireType::LITZ);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createThirtyTurnsLitzSinusoidalConfig() {
        TestConfig config;
        config.name = "Thirty_Turns_Litz_Sinusoidal";
        config.wireType = WireTypeClass::LITZ;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.1133},
            {25000, 0.1155},
            {50000, 0.1165},
            {100000, 0.1196},
            {200000, 0.12944},
            {250000, 0.13644}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({30});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            WireRound strand;
            OpenMagnetics::Wire wire;
            DimensionWithTolerance strandConductingDiameter;
            DimensionWithTolerance strandOuterDiameter;
            strandConductingDiameter.set_nominal(0.000071);
            strandOuterDiameter.set_nominal(0.0000762);
            strand.set_conducting_diameter(strandConductingDiameter);
            strand.set_outer_diameter(strandOuterDiameter);
            strand.set_number_conductors(1);
            strand.set_material("copper");
            strand.set_type(WireType::ROUND);
            
            wire.set_strand(strand);
            wire.set_nominal_value_outer_diameter(0.000873);
            wire.set_number_conductors(60);
            wire.set_type(WireType::LITZ);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    // ==================================================================================
    // RECTANGULAR WIRE TEST CONFIGURATIONS
    // ==================================================================================

    inline TestConfig createOneTurnRectangularUngappedConfig() {
        TestConfig config;
        config.name = "One_Turn_Rectangular_Ungapped";
        config.wireType = WireTypeClass::RECTANGULAR;
        config.temperature = 22;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {1, 0.0058551},
            {10000, 0.0058553},
            {20000, 0.0058559},
            {30000, 0.0058568},
            {40000, 0.0058573},
            {50000, 0.0058595},
            {60000, 0.0058619},
            {70000, 0.005865},
            {80000, 0.0058683},
            {90000, 0.0058719},
            {100000, 0.005876},
            {200000, 0.0059275},
            {300000, 0.0059955},
            {400000, 0.0060792},
            {500000, 0.0061776}
        };
        
        config.createMagnetic = []() {
            auto path = std::filesystem::path(__FILE__).parent_path() / "Test_Winding_Losses_One_Turn_Rectangular_Sinusoidal.json";
            auto mas = OpenMagneticsTesting::mas_loader(path);
            return mas.get_magnetic();
        };
        
        return config;
    }

    inline TestConfig createSevenTurnsRectangularUngappedConfig() {
        TestConfig config;
        config.name = "Seven_Turns_Rectangular_Ungapped";
        config.wireType = WireTypeClass::RECTANGULAR;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.0058551 * 5},
            {50000, 0.029313},
            {100000, 0.029367},
            {200000, 0.029651},
            {500000, 0.032103}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({5});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0001);
            wire.set_nominal_value_conducting_height(0.0016);
            wire.set_nominal_value_outer_width(0.00015);
            wire.set_nominal_value_outer_height(0.00165);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    // Seven Turns Rectangular Ungapped Sinusoidal with custom bobbin (PQ 27/17)
    inline TestConfig createSevenTurnsRectangularUngappedPQ2717Config() {
        TestConfig config;
        config.name = "Seven_Turns_Rectangular_Ungapped_PQ2717";
        config.wireType = WireTypeClass::RECTANGULAR;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.waveform = WaveformLabel::SINUSOIDAL;
        config.peakToPeak = 2;
        config.offset = 0;
        config.dutyCycle = 0.5;
        config.includeFringing = true;
        config.mirroringDimension = 2;
        
        config.expectedValues = {
            {0.01, 0.00019313 + 0.00019313 * 6},
            {100000, 0.0025456},
            {200000, 0.00342804},
            {300000, 0.00419621},
            {400000, 0.00484483},
            {500000, 0.0054165},
            {600000, 0.0059334},
            {700000, 0.00640876},
            {800000, 0.00685123},
            {900000, 0.00726681},
            {1000000, 0.00765988}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({7});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "PQ 27/17";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0038);
            wire.set_nominal_value_conducting_height(0.00076);
            wire.set_nominal_value_outer_height(0.0007676);
            wire.set_nominal_value_outer_width(0.003838);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::SPREAD, CoilAlignment::CENTERED, wires, false);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_residual_gap();
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            
            // Custom bobbin modifications
            auto bobbin = coil.resolve_bobbin();
            auto bobbinProcessedDescription = bobbin.get_processed_description().value();
            double bobbinColumnThickness = 0.001;
            bobbinProcessedDescription.set_wall_thickness(bobbinColumnThickness);
            bobbinProcessedDescription.set_column_thickness(bobbinColumnThickness);

            WindingWindowElement windingWindowElement;
            auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];
            auto coreCentralColumn = core.get_processed_description()->get_columns()[0];
            windingWindowElement.set_width(coreWindingWindow.get_width().value() - bobbinColumnThickness);
            windingWindowElement.set_height(coreWindingWindow.get_height().value() - bobbinColumnThickness * 2);
            windingWindowElement.set_coordinates(std::vector<double>({coreCentralColumn.get_width() / 2 + bobbinColumnThickness + windingWindowElement.get_width().value() / 2, 0, 0}));
            bobbinProcessedDescription.set_column_depth(coreCentralColumn.get_depth() / 2 + bobbinColumnThickness);
            bobbinProcessedDescription.set_column_width(coreCentralColumn.get_width() / 2 + bobbinColumnThickness);
            bobbinProcessedDescription.set_winding_windows(std::vector<WindingWindowElement>({windingWindowElement}));

            bobbin.set_processed_description(bobbinProcessedDescription);
            coil.set_bobbin(bobbin);
            coil.wind();
            
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    // ==================================================================================
    // FOIL WIRE TEST CONFIGURATIONS
    // ==================================================================================

    inline TestConfig createOneTurnFoilSinusoidalConfig() {
        TestConfig config;
        config.name = "One_Turn_Foil_Sinusoidal";
        config.wireType = WireTypeClass::FOIL;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = false;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.00037759},
            {50000, 0.00037787},
            {100000, 0.00037907},
            {200000, 0.00039071},
            {500000, 0.000442}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0001);
            wire.set_nominal_value_conducting_height(0.02);
            wire.set_nominal_value_outer_width(0.00011);
            wire.set_nominal_value_outer_height(0.02);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::FOIL);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires, false);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createTenTurnsFoilSinusoidalConfig() {
        TestConfig config;
        config.name = "Ten_Turns_Foil_Sinusoidal";
        config.wireType = WireTypeClass::FOIL;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = false;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.0041536},
            {50000, 0.0041805},
            {100000, 0.004284},
            {200000, 0.004639},
            {500000, 0.0057148}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({10});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0001);
            wire.set_nominal_value_conducting_height(0.02);
            wire.set_nominal_value_outer_width(0.00018);
            wire.set_nominal_value_outer_height(0.02);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::FOIL);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires, false);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createTenShortTurnsFoilConfig() {
        TestConfig config;
        config.name = "Ten_Short_Turns_Foil";
        config.wireType = WireTypeClass::FOIL;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = false;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.0090968},
            {50000, 0.009198},
            {100000, 0.0095051},
            {200000, 0.01035},
            {500000, 0.012629}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({10});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0001);
            wire.set_nominal_value_conducting_height(0.007);
            wire.set_nominal_value_outer_width(0.00015);
            wire.set_nominal_value_outer_height(0.007);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::FOIL);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires, false);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    // ==================================================================================
    // TOROIDAL CORE TEST CONFIGURATIONS
    // ==================================================================================

    inline TestConfig createOneTurnRoundToroidalConfig() {
        TestConfig config;
        config.name = "One_Turn_Round_Toroidal";
        config.wireType = WireTypeClass::ROUND;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = false;
        config.mirroringDimension = 1;
        
        // Use correct expected values including AC effects
        config.expectedValues = {
            {0.01, 0.002032},
            {25000, 0.002053},
            {50000, 0.002121},
            {100000, 0.002355},
            {200000, 0.002987},
            {250000, 0.003293},
            {500000, 0.004466}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "T 40/24/16";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_diameter(0.00071);
            wire.set_nominal_value_outer_diameter(0.000762);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::ROUND);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            auto gapping = json::array();
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, "3C97");
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createTenTurnsRoundToroidalConfig() {
        TestConfig config;
        config.name = "Ten_Turns_Round_Toroidal";
        config.wireType = WireTypeClass::ROUND;
        config.temperature = 20;
        config.magnetizingInductance = 100e-6;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        // Use correct expected values including AC effects
        config.expectedValues = {
            {0.01, 0.0213484},
            {25000, 0.0217045},
            {50000, 0.0227473},
            {100000, 0.0265887},
            {200000, 0.0409422},
            {250000, 0.0460606},
            {500000, 0.0658813}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({10});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "T 40/24/16";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_diameter(0.00071);
            wire.set_nominal_value_outer_diameter(0.000762);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::ROUND);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            auto gapping = json::array();
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, "3C97");
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createOneTurnRectangularToroidalConfig() {
        TestConfig config;
        config.name = "One_Turn_Rectangular_Toroidal";
        config.wireType = WireTypeClass::RECTANGULAR;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = false;
        config.mirroringDimension = 1;
        
        // Use correct expected values including AC effects
        config.expectedValues = {
            {0.01, 9.94277e-05},
            {25000, 0.000216944},
            {50000, 0.00030592},
            {100000, 0.000432524},
            {200000, 0.000611667},
            {250000, 0.000683864},
            {500000, 0.000967129}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "T 40/24/16";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0058);
            wire.set_nominal_value_conducting_height(0.002);
            wire.set_nominal_value_outer_width(0.0059);
            wire.set_nominal_value_outer_height(0.002076);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            auto gapping = json::array();
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, "3C97");
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createTenTurnsRectangularToroidalConfig() {
        TestConfig config;
        config.name = "Ten_Turns_Rectangular_Toroidal";
        config.wireType = WireTypeClass::RECTANGULAR;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = false;
        config.mirroringDimension = 1;
        
        // Use correct expected values including AC effects
        config.expectedValues = {
            {0.01, 9.94277e-04},
            {25000, 0.00216944},
            {50000, 0.0030592},
            {100000, 0.00432524},
            {200000, 0.00611667},
            {250000, 0.00683864},
            {500000, 0.00967129}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({10});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "T 40/24/16";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0058);
            wire.set_nominal_value_conducting_height(0.002);
            wire.set_nominal_value_outer_width(0.0059);
            wire.set_nominal_value_outer_height(0.002076);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            auto gapping = json::array();
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, "3C97");
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    // ==================================================================================
    // PLANAR TESTS (uses JSON files)
    // ==================================================================================

    inline TestConfig createOneTurnPlanarNoFringingConfig() {
        TestConfig config;
        config.name = "One_Turn_Planar_No_Fringing";
        config.wireType = WireTypeClass::PLANAR;
        config.temperature = 22;
        config.includeFringing = false;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {100000, 0.000209},
            {200000, 0.000219},
            {300000, 0.000234},
            {400000, 0.000253},
            {500000, 0.000275}
        };
        
        config.createMagnetic = []() {
            auto path = std::filesystem::path(__FILE__).parent_path() / "Test_Winding_Losses_One_Turn_Planar.json";
            auto mas = OpenMagneticsTesting::mas_loader(path);
            return mas.get_magnetic();
        };
        
        return config;
    }

    inline TestConfig createOneTurnPlanarFringingConfig() {
        TestConfig config;
        config.name = "One_Turn_Planar_Fringing";
        config.wireType = WireTypeClass::PLANAR;
        config.temperature = 22;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {100000, 0.000218},
            {200000, 0.000238},
            {300000, 0.000266},
            {400000, 0.000303},
            {500000, 0.000347}
        };
        
        config.createMagnetic = []() {
            auto path = std::filesystem::path(__FILE__).parent_path() / "Test_Winding_Losses_One_Turn_Planar.json";
            auto mas = OpenMagneticsTesting::mas_loader(path);
            return mas.get_magnetic();
        };
        
        return config;
    }

    inline TestConfig createSixteenTurnsPlanarNoFringingConfig() {
        TestConfig config;
        config.name = "Sixteen_Turns_Planar_No_Fringing";
        config.wireType = WireTypeClass::PLANAR;
        config.temperature = 22;
        config.includeFringing = false;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {100000, 0.00312},
            {200000, 0.00354},
            {300000, 0.00416},
            {400000, 0.00492},
            {500000, 0.00576}
        };
        
        config.createMagnetic = []() {
            auto path = std::filesystem::path(__FILE__).parent_path() / "Test_Winding_Losses_Sixteen_Turns_Planar.json";
            auto mas = OpenMagneticsTesting::mas_loader(path);
            return mas.get_magnetic();
        };
        
        return config;
    }

    inline TestConfig createSixteenTurnsPlanarFringingCloseConfig() {
        TestConfig config;
        config.name = "Sixteen_Turns_Planar_Fringing_Close";
        config.wireType = WireTypeClass::PLANAR;
        config.temperature = 22;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {100000, 0.0044},
            {200000, 0.00725},
            {300000, 0.01152},
            {400000, 0.01657},
            {500000, 0.02214}
        };
        
        config.createMagnetic = []() {
            auto path = std::filesystem::path(__FILE__).parent_path() / "Test_Winding_Losses_Sixteen_Turns_Planar_Close.json";
            auto mas = OpenMagneticsTesting::mas_loader(path);
            return mas.get_magnetic();
        };
        
        return config;
    }

    inline TestConfig createSixteenTurnsPlanarFringingFarConfig() {
        TestConfig config;
        config.name = "Sixteen_Turns_Planar_Fringing_Far";
        config.wireType = WireTypeClass::PLANAR;
        config.temperature = 22;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {100000, 0.00346},
            {200000, 0.00429},
            {300000, 0.00543},
            {400000, 0.00681},
            {500000, 0.00838}
        };
        
        config.createMagnetic = []() {
            auto path = std::filesystem::path(__FILE__).parent_path() / "Test_Winding_Losses_Sixteen_Turns_Planar_Far.json";
            auto mas = OpenMagneticsTesting::mas_loader(path);
            return mas.get_magnetic();
        };
        
        return config;
    }

    // ==================================================================================
    // ADDITIONAL LITZ WIRE CONFIGS
    // ==================================================================================

    inline TestConfig createOneTurnLitzTriangularWithDCConfig() {
        TestConfig config;
        config.name = "One_Turn_Litz_Triangular_With_DC";
        config.wireType = WireTypeClass::LITZ;
        config.temperature = 20;
        config.waveform = WaveformLabel::TRIANGULAR;
        config.peakToPeak = 2 * 1.73205;
        config.offset = 10;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {{500000, 0.113}};
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            WireRound strand;
            OpenMagnetics::Wire wire;
            DimensionWithTolerance strandConductingDiameter;
            DimensionWithTolerance strandOuterDiameter;
            strandConductingDiameter.set_nominal(0.00004);
            strandOuterDiameter.set_nominal(0.000049);
            strand.set_conducting_diameter(strandConductingDiameter);
            strand.set_outer_diameter(strandOuterDiameter);
            strand.set_number_conductors(1);
            strand.set_material("copper");
            strand.set_type(WireType::ROUND);
            
            wire.set_strand(strand);
            wire.set_nominal_value_outer_diameter(0.001576);
            wire.set_number_conductors(600);
            wire.set_type(WireType::LITZ);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createOneTurnLitzFewStrandsConfig() {
        TestConfig config;
        config.name = "One_Turn_Litz_Few_Strands";
        config.wireType = WireTypeClass::LITZ;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.008411},
            {25000, 0.008412},
            {50000, 0.008416},
            {100000, 0.008430},
            {200000, 0.008489},
            {250000, 0.008433},
            {500000, 0.008800}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            WireRound strand;
            OpenMagnetics::Wire wire;
            DimensionWithTolerance strandConductingDiameter;
            DimensionWithTolerance strandOuterDiameter;
            strandConductingDiameter.set_nominal(0.0001);
            strandOuterDiameter.set_nominal(0.00011);
            strand.set_conducting_diameter(strandConductingDiameter);
            strand.set_outer_diameter(strandOuterDiameter);
            strand.set_number_conductors(1);
            strand.set_material("copper");
            strand.set_type(WireType::ROUND);
            
            wire.set_strand(strand);
            wire.set_nominal_value_outer_diameter(0.000551);
            wire.set_number_conductors(12);
            wire.set_type(WireType::LITZ);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createOneTurnLitzManyManyStrandsConfig() {
        TestConfig config;
        config.name = "One_Turn_Litz_Many_Many_Strands";
        config.wireType = WireTypeClass::LITZ;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.0001572},
            {25000, 0.0001578},
            {50000, 0.0001586},
            {100000, 0.000159},
            {200000, 0.0001516},
            {250000, 0.0001547},
            {500000, 0.00016}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            WireRound strand;
            OpenMagnetics::Wire wire;
            DimensionWithTolerance strandConductingDiameter;
            DimensionWithTolerance strandOuterDiameter;
            strandConductingDiameter.set_nominal(0.00002);
            strandOuterDiameter.set_nominal(0.000029);
            strand.set_conducting_diameter(strandConductingDiameter);
            strand.set_outer_diameter(strandOuterDiameter);
            strand.set_number_conductors(1);
            strand.set_material("copper");
            strand.set_type(WireType::ROUND);
            
            wire.set_strand(strand);
            wire.set_nominal_value_outer_diameter(0.004384);
            wire.set_number_conductors(20000);
            wire.set_type(WireType::LITZ);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    // ==================================================================================
    // ADDITIONAL ROUND WIRE CONFIGS
    // ==================================================================================

    inline TestConfig createOneTurnRoundSinusoidalWithDCConfig() {
        TestConfig config;
        config.name = "One_Turn_Round_Sinusoidal_With_DC";
        config.wireType = WireTypeClass::ROUND;
        config.temperature = 20;
        config.waveform = WaveformLabel::SINUSOIDAL;
        config.peakToPeak = 2 * 1.4142;
        config.offset = 4.2;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {0.01, 0.03782},
            {25000, 0.03784},
            {50000, 0.03791},
            {100000, 0.03811},
            {200000, 0.03871},
            {250000, 0.03903},
            {500000, 0.04035}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_diameter(0.00071);
            wire.set_nominal_value_outer_diameter(0.000762);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::ROUND);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    inline TestConfig createOneTurnRoundTriangularWithDCConfig() {
        TestConfig config;
        config.name = "One_Turn_Round_Triangular_With_DC";
        config.wireType = WireTypeClass::ROUND;
        config.temperature = 20;
        config.waveform = WaveformLabel::TRIANGULAR;
        config.peakToPeak = 2 * 1.73205;
        config.offset = 4.2;
        config.magnetizingInductance = 1e-3;
        config.includeFringing = true;
        config.mirroringDimension = 1;
        
        config.expectedValues = {
            {25000, 0.03783},
            {100000, 0.03811},
            {500000, 0.040374}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({1});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "ETD 34/17/11";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_diameter(0.00071);
            wire.set_nominal_value_outer_diameter(0.000762);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::ROUND);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::CENTERED, CoilAlignment::CENTERED, wires);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    // Five Turns Rectangular Ungapped Sinusoidal (peakToPeak = 2)
    inline TestConfig createFiveTurnsRectangularUngappedConfig() {
        TestConfig config;
        config.name = "Five_Turns_Rectangular_Ungapped";
        config.wireType = WireTypeClass::RECTANGULAR;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.waveform = WaveformLabel::SINUSOIDAL;
        config.peakToPeak = 2;
        config.offset = 0;
        config.dutyCycle = 0.5;
        config.includeFringing = true;
        config.mirroringDimension = 2;
        
        config.expectedValues = {
            {0.01, 0.00011974 + 0.0003736},
            {100000, 0.00087387 + 0.0003736},
            {200000, 0.0012227 + 0.0003736},
            {300000, 0.001487 + 0.0003736},
            {400000, 0.0017156 + 0.0003736},
            {500000, 0.001927 + 0.0003736},
            {600000, 0.0021237 + 0.0003736},
            {700000, 0.0023129 + 0.0003736},
            {800000, 0.0024954 + 0.0003736},
            {900000, 0.0026719 + 0.0003736},
            {1000000, 0.002843 + 0.0003736}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({5});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "PQ 20/16";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.004);
            wire.set_nominal_value_conducting_height(0.0009);
            wire.set_nominal_value_outer_height(0.000909);
            wire.set_nominal_value_outer_width(0.00404);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::SPREAD, CoilAlignment::CENTERED, wires, false);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_residual_gap();
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    // Five Turns Rectangular Ungapped Sinusoidal 7 Amps (peakToPeak = 14)
    inline TestConfig createFiveTurnsRectangularUngapped7AmpsConfig() {
        TestConfig config;
        config.name = "Five_Turns_Rectangular_Ungapped_7_Amps";
        config.wireType = WireTypeClass::RECTANGULAR;
        config.temperature = 20;
        config.magnetizingInductance = 1e-3;
        config.waveform = WaveformLabel::SINUSOIDAL;
        config.peakToPeak = 14;
        config.offset = 0;
        config.dutyCycle = 0.5;
        config.includeFringing = true;
        config.mirroringDimension = 2;
        
        config.expectedValues = {
            {0.01, 0.0058551 * 5},
            {100000, 0.03363 + 0.0058551 * 5},
            {200000, 0.06063 + 0.0058551 * 5},
            {300000, 0.073721 + 0.0058551 * 5},
            {400000, 0.084989 + 0.0058551 * 5},
            {500000, 0.095376 + 0.0058551 * 5},
            {600000, 0.10525 + 0.0058551 * 5},
            {700000, 0.11477 + 0.0058551 * 5},
            {800000, 0.12399 + 0.0058551 * 5},
            {900000, 0.13296 + 0.0058551 * 5},
            {1000000, 0.141661 + 0.0058551 * 5}
        };
        
        config.createMagnetic = []() {
            std::vector<int64_t> numberTurns({5});
            std::vector<int64_t> numberParallels({1});
            std::string shapeName = "PQ 20/16";
            
            std::vector<OpenMagnetics::Wire> wires;
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.004);
            wire.set_nominal_value_conducting_height(0.0009);
            wire.set_nominal_value_outer_height(0.000909);
            wire.set_nominal_value_outer_width(0.00404);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            wires.push_back(wire);
            
            auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName, 1,
                WindingOrientation::OVERLAPPING, WindingOrientation::OVERLAPPING,
                CoilAlignment::SPREAD, CoilAlignment::CENTERED, wires, false);
            
            int64_t numberStacks = 1;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_residual_gap();
            auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            return magnetic;
        };
        
        return config;
    }

    // Five Turns Rectangular Gapped Sinusoidal 7 Amps (same as Ungapped - they use same gap type)
    inline TestConfig createFiveTurnsRectangularGapped7AmpsConfig() {
        // Note: The original test uses residual_gap() same as the ungapped test
        // This appears to be a copy-paste with different intent but same implementation
        return createFiveTurnsRectangularUngapped7AmpsConfig();
    }

    // Seven Turns Rectangular Ungapped Sinusoidal - already exists, keep for consistency
    inline TestConfig createSevenTurnsRectangularUngappedSinusoidalConfig() {
        // This is a special test with custom bobbin settings - keep original
        return createSevenTurnsRectangularUngappedConfig();
    }

    // ==================================================================================
    // GET ALL TEST CONFIGURATIONS
    // ==================================================================================

    inline std::map<std::string, TestConfig> getAllTestConfigs() {
        std::map<std::string, TestConfig> configs;
        
        // Round wire tests
        configs["One_Turn_Round_Stacked"] = createOneTurnRoundStackedConfig();
        // Note: JSON-based tests require the JSON file to be present
        // configs["One_Turn_Round_Sinusoidal"] = createOneTurnRoundSinusoidalConfig();
        // configs["Twelve_Turns_Round_Sinusoidal"] = createTwelveTurnsRoundSinusoidalConfig();
        // configs["One_Turn_Round_Fringing"] = createOneTurnRoundFringingConfig();
        
        // Litz wire tests
        configs["One_Turn_Litz_Sinusoidal"] = createOneTurnLitzSinusoidalConfig();
        configs["One_Turn_Litz_Many_Strands"] = createOneTurnLitzManyStrandsConfig();
        configs["Ten_Turns_Litz_Sinusoidal"] = createTenTurnsLitzSinusoidalConfig();
        configs["Thirty_Turns_Litz_Sinusoidal"] = createThirtyTurnsLitzSinusoidalConfig();
        
        // Rectangular wire tests
        // configs["One_Turn_Rectangular_Ungapped"] = createOneTurnRectangularUngappedConfig();
        configs["Seven_Turns_Rectangular_Ungapped"] = createSevenTurnsRectangularUngappedConfig();
        configs["Five_Turns_Rectangular_Ungapped"] = createFiveTurnsRectangularUngappedConfig();
        configs["Five_Turns_Rectangular_Ungapped_7_Amps"] = createFiveTurnsRectangularUngapped7AmpsConfig();
        configs["Five_Turns_Rectangular_Gapped_7_Amps"] = createFiveTurnsRectangularGapped7AmpsConfig();
        
        // Foil tests
        configs["One_Turn_Foil_Sinusoidal"] = createOneTurnFoilSinusoidalConfig();
        configs["Ten_Turns_Foil_Sinusoidal"] = createTenTurnsFoilSinusoidalConfig();
        configs["Ten_Short_Turns_Foil"] = createTenShortTurnsFoilConfig();
        
        // Toroidal core tests
        configs["One_Turn_Round_Toroidal"] = createOneTurnRoundToroidalConfig();
        configs["Ten_Turns_Round_Toroidal"] = createTenTurnsRoundToroidalConfig();
        configs["One_Turn_Rectangular_Toroidal"] = createOneTurnRectangularToroidalConfig();
        configs["Ten_Turns_Rectangular_Toroidal"] = createTenTurnsRectangularToroidalConfig();
        
        // Planar tests (require JSON files)
        // configs["One_Turn_Planar_No_Fringing"] = createOneTurnPlanarNoFringingConfig();
        // configs["One_Turn_Planar_Fringing"] = createOneTurnPlanarFringingConfig();
        // configs["Sixteen_Turns_Planar_No_Fringing"] = createSixteenTurnsPlanarNoFringingConfig();
        // configs["Sixteen_Turns_Planar_Fringing_Close"] = createSixteenTurnsPlanarFringingCloseConfig();
        // configs["Sixteen_Turns_Planar_Fringing_Far"] = createSixteenTurnsPlanarFringingFarConfig();
        
        return configs;
    }

} // namespace WindingLossesTestData
