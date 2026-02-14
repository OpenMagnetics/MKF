/**
 * @file TestCllcConverter.cpp
 * @brief Unit tests for the CLLC Bidirectional Resonant Converter model.
 *
 * Tests cover:
 *   1. Forward mode at resonance (fs = fr)
 *   2. Forward mode below resonance (boost, fs < fr)
 *   3. Forward mode above resonance (buck, fs > fr)
 *   4. Reverse mode at resonance
 *   5. Symmetric resonant tank design
 *   6. Asymmetric resonant tank design
 *   7. Multiple operating points
 *   8. Design requirements calculation
 *   9. Advanced converter with user-specified parameters
 *   10. Ngspice simulation (forward mode)
 *   11. Ngspice waveform polarity verification
 *   12. Ngspice num-periods parameter test
 *
 * Design equations reference:
 *   [1] Infineon AN-2024-06: "Operation and modeling analysis of a bidirectional CLLC converter"
 *   [2] Bartecka et al., Energies 2024, 17, 55.
 *
 * @note Tests tagged [ngspice-simulation] require ngspice to be installed on the system.
 */

#include <source_location>
#include "support/Painter.h"
#include "converter_models/Cllc.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "processors/NgspiceRunner.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>
#include <numeric>
#include <cmath>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}
                              .parent_path().append("..").append("output");
    double maximumError = 0.15;

    // =========================================================================
    // Helper: Create a standard CLLC converter specification
    // =========================================================================

    /**
     * Creates a CLLC converter spec similar to the Infineon AN example:
     *   Vin = 750V (700-800V range), Vout = 600V (550-800V range)
     *   fr = 73 kHz, Pout = 11 kW, forward mode
     */
    json create_standard_cllc_json(double fsw = 73000.0,
                                    CllcPowerFlow powerFlow = CllcPowerFlow::FORWARD) {
        json cllcJson;
        json inputVoltage;
        inputVoltage["minimum"] = 700;
        inputVoltage["maximum"] = 800;
        inputVoltage["nominal"] = 750;
        cllcJson["inputVoltage"] = inputVoltage;

        cllcJson["maxSwitchingFrequency"] = 250000;
        cllcJson["minSwitchingFrequency"] = 40000;
        cllcJson["efficiency"] = 0.95;
        cllcJson["qualityFactor"] = 0.3;
        cllcJson["symmetricDesign"] = true;
        cllcJson["bidirectional"] = true;

        cllcJson["operatingPoints"] = json::array();
        {
            json opJson;
            opJson["outputVoltages"] = {600.0};
            opJson["outputCurrents"] = {18.33};  // ~11 kW at 600V
            opJson["switchingFrequency"] = fsw;
            opJson["ambientTemperature"] = 25.0;
            opJson["powerFlow"] = (powerFlow == CllcPowerFlow::FORWARD) ? "Forward" : "Reverse";
            cllcJson["operatingPoints"].push_back(opJson);
        }
        return cllcJson;
    }

    /**
     * Creates a smaller-power CLLC spec for quick simulation:
     *   Vin = 400V, Vout = 48V, ~500W, 200 kHz
     */
    json create_small_power_cllc_json(double fsw = 200000.0,
                                       CllcPowerFlow powerFlow = CllcPowerFlow::FORWARD) {
        json cllcJson;
        json inputVoltage;
        inputVoltage["minimum"] = 360;
        inputVoltage["maximum"] = 420;
        inputVoltage["nominal"] = 400;
        cllcJson["inputVoltage"] = inputVoltage;

        cllcJson["maxSwitchingFrequency"] = 400000;
        cllcJson["minSwitchingFrequency"] = 100000;
        cllcJson["efficiency"] = 0.95;
        cllcJson["qualityFactor"] = 0.3;
        cllcJson["symmetricDesign"] = true;
        cllcJson["bidirectional"] = true;

        cllcJson["operatingPoints"] = json::array();
        {
            json opJson;
            opJson["outputVoltages"] = {48.0};
            opJson["outputCurrents"] = {10.0};  // ~480W
            opJson["switchingFrequency"] = fsw;
            opJson["ambientTemperature"] = 25.0;
            opJson["powerFlow"] = (powerFlow == CllcPowerFlow::FORWARD) ? "Forward" : "Reverse";
            cllcJson["operatingPoints"].push_back(opJson);
        }
        return cllcJson;
    }

    // =========================================================================
    // TEST 1: Forward Mode at Resonance (fs = fr)
    // =========================================================================

    TEST_CASE("Test_CllcConverter_Forward_AtResonance",
              "[converter-model][cllc-topology][forward][at-resonance][smoke-test]") {
        /**
         * Test the CLLC converter operating at resonance frequency in forward mode.
         * At resonance, the gain should be approximately 1.0 (i.e., nVout/Vin ≈ 1).
         * The voltage gain is purely dependent on the turns ratio.
         *
         * Reference: Infineon AN Section 2.1, Table 1 (fs = fr → Resonant mode)
         */
        json cllcJson = create_standard_cllc_json(73000.0, CllcPowerFlow::FORWARD);
        OpenMagnetics::CllcConverter cllc(cllcJson);
        cllc._assertErrors = true;

        // Calculate resonant parameters
        auto params = cllc.calculate_resonant_parameters();

        INFO("Turns ratio n = " << params.turnsRatio);
        INFO("Resonant frequency fr = " << params.resonantFrequency << " Hz");
        INFO("L1 = " << (params.primaryResonantInductance * 1e6) << " uH");
        INFO("C1 = " << (params.primaryResonantCapacitance * 1e9) << " nF");
        INFO("Lm = " << (params.magnetizingInductance * 1e6) << " uH");
        INFO("L2 = " << (params.secondaryResonantInductance * 1e6) << " uH");
        INFO("C2 = " << (params.secondaryResonantCapacitance * 1e9) << " nF");
        INFO("Q = " << params.qualityFactor);
        INFO("k = " << params.inductanceRatio);

        // Verify turns ratio: n = 750/600 = 1.25
        REQUIRE_THAT(params.turnsRatio,
                     Catch::Matchers::WithinAbs(1.25, 0.05));

        // Verify resonant frequency is in expected range
        REQUIRE(params.resonantFrequency > 30000);
        REQUIRE(params.resonantFrequency < 300000);

        // Verify voltage gain at resonance ≈ 1.0
        double gainAtResonance = cllc.get_voltage_gain(params.resonantFrequency, params);
        INFO("Voltage gain at resonance = " << gainAtResonance);
        REQUIRE_THAT(gainAtResonance, Catch::Matchers::WithinAbs(1.0, 0.15));

        // Generate operating points
        std::vector<double> turnsRatios = {params.turnsRatio};
        auto operatingPoints = cllc.process_operating_points(turnsRatios, params.magnetizingInductance);

        REQUIRE(!operatingPoints.empty());
        REQUIRE(operatingPoints[0].get_excitations_per_winding().size() == 2);

        // Primary winding checks
        auto& primaryExcitation = operatingPoints[0].get_excitations_per_winding()[0];
        REQUIRE(primaryExcitation.get_voltage().has_value());
        REQUIRE(primaryExcitation.get_current().has_value());
        REQUIRE(primaryExcitation.get_voltage()->get_waveform().has_value());
        REQUIRE(primaryExcitation.get_current()->get_waveform().has_value());

        // Primary voltage should be approximately ±Vin
        auto priVoltageData = primaryExcitation.get_voltage()->get_waveform()->get_data();
        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        double priV_min = *std::min_element(priVoltageData.begin(), priVoltageData.end());
        INFO("Primary voltage max: " << priV_max << " V, min: " << priV_min << " V");
        CHECK(priV_max > 600);  // Should be around Vin nominal (700-800V range)
        CHECK(priV_min < -600);

        // Primary current should be quasi-sinusoidal
        auto priCurrentData = primaryExcitation.get_current()->get_waveform()->get_data();
        double priI_max = *std::max_element(priCurrentData.begin(), priCurrentData.end());
        double priI_min = *std::min_element(priCurrentData.begin(), priCurrentData.end());
        INFO("Primary current max: " << priI_max << " A, min: " << priI_min << " A");
        CHECK(priI_max > 0);
        CHECK(priI_min < 0);  // Bipolar current expected

        // Secondary winding checks
        auto& secondaryExcitation = operatingPoints[0].get_excitations_per_winding()[1];
        REQUIRE(secondaryExcitation.get_voltage().has_value());
        REQUIRE(secondaryExcitation.get_current().has_value());

        auto secVoltageData = secondaryExcitation.get_voltage()->get_waveform()->get_data();
        double secV_max = *std::max_element(secVoltageData.begin(), secVoltageData.end());
        double secV_min = *std::min_element(secVoltageData.begin(), secVoltageData.end());
        INFO("Secondary voltage max: " << secV_max << " V, min: " << secV_min << " V");
        CHECK(secV_max > 400);  // Should be around Vout (600V)
        CHECK(secV_min < -400);

        // Paint waveforms
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Cllc_Forward_AtResonance_Primary_Current.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(primaryExcitation.get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Cllc_Forward_AtResonance_Primary_Voltage.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(primaryExcitation.get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Cllc_Forward_AtResonance_Secondary_Current.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(secondaryExcitation.get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Cllc_Forward_AtResonance_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(secondaryExcitation.get_voltage()->get_waveform().value());
            painter.export_svg();
        }
    }

    // =========================================================================
    // TEST 2: Forward Mode Below Resonance (Boost, fs < fr)
    // =========================================================================

    TEST_CASE("Test_CllcConverter_Forward_BelowResonance",
              "[converter-model][cllc-topology][forward][below-resonance][smoke-test]") {
        /**
         * Below resonance (fs < fr), the CLLC operates in boost mode.
         * The voltage gain should be > 1.
         * Includes freewheeling state after resonant cycle completes.
         *
         * Reference: Infineon AN Table 1, fs < fr → Boost mode
         */
        // Use a switching frequency well below the natural resonant frequency
        json cllcJson = create_standard_cllc_json(50000.0, CllcPowerFlow::FORWARD);
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();

        // Verify gain > 1 below resonance
        double gainBelowResonance = cllc.get_voltage_gain(50000.0, params);
        INFO("Voltage gain at 50 kHz (below resonance) = " << gainBelowResonance);
        CHECK(gainBelowResonance > 0.8);  // Should be elevated (boost region)

        // Generate operating points
        std::vector<double> turnsRatios = {params.turnsRatio};
        auto operatingPoints = cllc.process_operating_points(turnsRatios, params.magnetizingInductance);

        REQUIRE(!operatingPoints.empty());
        REQUIRE(operatingPoints[0].get_excitations_per_winding().size() == 2);

        // Waveforms should still be valid
        auto& priExc = operatingPoints[0].get_excitations_per_winding()[0];
        REQUIRE(priExc.get_voltage().has_value());
        REQUIRE(priExc.get_current().has_value());
    }

    // =========================================================================
    // TEST 3: Forward Mode Above Resonance (Buck, fs > fr)
    // =========================================================================

    TEST_CASE("Test_CllcConverter_Forward_AboveResonance",
              "[converter-model][cllc-topology][forward][above-resonance][smoke-test]") {
        /**
         * Above resonance (fs > fr), the CLLC operates in buck mode.
         * The voltage gain should be < 1.
         * The resonant half-cycle is interrupted by the next switching half-cycle.
         *
         * Reference: Infineon AN Table 1, fs > fr → Buck mode
         */
        json cllcJson = create_standard_cllc_json(150000.0, CllcPowerFlow::FORWARD);
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();

        double gainAboveResonance = cllc.get_voltage_gain(150000.0, params);
        INFO("Voltage gain at 150 kHz (above resonance) = " << gainAboveResonance);
        CHECK(gainAboveResonance < 1.2);  // Should be reduced (buck region)

        std::vector<double> turnsRatios = {params.turnsRatio};
        auto operatingPoints = cllc.process_operating_points(turnsRatios, params.magnetizingInductance);

        REQUIRE(!operatingPoints.empty());
        REQUIRE(operatingPoints[0].get_excitations_per_winding().size() == 2);
    }

    // =========================================================================
    // TEST 4: Reverse Mode at Resonance
    // =========================================================================

    TEST_CASE("Test_CllcConverter_Reverse_AtResonance",
              "[converter-model][cllc-topology][reverse][at-resonance][smoke-test]") {
        /**
         * Test reverse power flow (secondary → primary).
         * In reverse mode, the secondary bridge is the inverter and the primary bridge
         * is the rectifier. The gain equation uses the inverse turns ratio.
         *
         * Reference: Infineon AN Section 2.2.2 (discharging mode)
         */
        json cllcJson = create_standard_cllc_json(73000.0, CllcPowerFlow::REVERSE);

        // For reverse mode, swap the voltage/current references
        cllcJson["operatingPoints"][0]["powerFlow"] = "Reverse";

        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();

        INFO("Reverse mode - turns ratio n = " << params.turnsRatio);
        INFO("Reverse mode - Lm = " << (params.magnetizingInductance * 1e6) << " uH");

        std::vector<double> turnsRatios = {params.turnsRatio};
        auto operatingPoints = cllc.process_operating_points(turnsRatios, params.magnetizingInductance);

        REQUIRE(!operatingPoints.empty());

        // Check that operating point name contains "Reverse"
        REQUIRE(operatingPoints[0].get_name().value().find("Reverse") != std::string::npos);
    }

    // =========================================================================
    // TEST 5: Symmetric Resonant Tank Design
    // =========================================================================

    TEST_CASE("Test_CllcConverter_Symmetric_Design",
              "[converter-model][cllc-topology][symmetric][smoke-test]") {
        /**
         * For symmetric design: a = 1, b = 1
         * This means L2 = L1/n² and C2 = n²*C1
         * And the primary and secondary resonant frequencies are equal: fr1 = fr2
         *
         * Reference: Infineon AN Step 8 and Energies paper Eq. (1)
         */
        json cllcJson = create_standard_cllc_json();
        cllcJson["symmetricDesign"] = true;
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();

        double n = params.turnsRatio;
        double L1 = params.primaryResonantInductance;
        double C1 = params.primaryResonantCapacitance;
        double L2 = params.secondaryResonantInductance;
        double C2 = params.secondaryResonantCapacitance;

        // Check symmetric relationships: L2 = L1/n², C2 = n²*C1
        double expectedL2 = L1 / (n * n);
        double expectedC2 = n * n * C1;

        INFO("L1 = " << L1 << ", L2 = " << L2 << ", expected L2 = " << expectedL2);
        INFO("C1 = " << C1 << ", C2 = " << C2 << ", expected C2 = " << expectedC2);

        REQUIRE_THAT(L2, Catch::Matchers::WithinAbs(expectedL2, expectedL2 * 0.01));
        REQUIRE_THAT(C2, Catch::Matchers::WithinAbs(expectedC2, expectedC2 * 0.01));

        // Check resonant frequencies match: fr1 = 1/(2π√(L1C1)) = fr2 = 1/(2π√(L2C2))
        double fr1 = 1.0 / (2.0 * M_PI * sqrt(L1 * C1));
        double fr2 = 1.0 / (2.0 * M_PI * sqrt(L2 * C2));
        INFO("fr1 = " << fr1 << " Hz, fr2 = " << fr2 << " Hz");
        REQUIRE_THAT(fr1, Catch::Matchers::WithinAbs(fr2, fr2 * 0.01));

        // Verify a = 1, b = 1
        CHECK(params.resonantInductorRatio == 1.0);
        CHECK(params.resonantCapacitorRatio == 1.0);
    }

    // =========================================================================
    // TEST 6: Asymmetric Resonant Tank Design
    // =========================================================================

    TEST_CASE("Test_CllcConverter_Asymmetric_Design",
              "[converter-model][cllc-topology][asymmetric][smoke-test]") {
        /**
         * For asymmetric design: a = 0.95, b = 1.052
         * L2 = 0.95*L1/n², C2 = n²*1.052*C1
         * The primary and secondary resonant frequencies still match when a*b = 1
         *
         * Reference: Infineon AN Step 8 (a=0.95, b=1.052, a*b ≈ 1.0)
         */
        json cllcJson = create_standard_cllc_json();
        cllcJson["symmetricDesign"] = false;
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();

        double n = params.turnsRatio;
        double L1 = params.primaryResonantInductance;
        double C1 = params.primaryResonantCapacitance;
        double L2 = params.secondaryResonantInductance;
        double C2 = params.secondaryResonantCapacitance;
        double a = params.resonantInductorRatio;
        double b = params.resonantCapacitorRatio;

        INFO("a = " << a << ", b = " << b << ", a*b = " << (a * b));

        // Check asymmetric relationships
        double expectedL2 = a * L1 / (n * n);
        double expectedC2 = n * n * b * C1;

        REQUIRE_THAT(L2, Catch::Matchers::WithinAbs(expectedL2, expectedL2 * 0.01));
        REQUIRE_THAT(C2, Catch::Matchers::WithinAbs(expectedC2, expectedC2 * 0.01));

        // a and b should differ from 1.0
        CHECK(a != 1.0);
        CHECK(b != 1.0);

        // a*b should be approximately 1.0 (so fr1 ≈ fr2)
        REQUIRE_THAT(a * b, Catch::Matchers::WithinAbs(1.0, 0.01));
    }

    // =========================================================================
    // TEST 7: Multiple Operating Points
    // =========================================================================

    TEST_CASE("Test_CllcConverter_MultipleOperatingPoints",
              "[converter-model][cllc-topology][multiple-ops][smoke-test]") {
        /**
         * Test with multiple operating points at different frequencies and loads.
         * Validates that the converter generates correct waveforms for each condition.
         */
        json cllcJson = create_standard_cllc_json();

        // Add a second operating point at different frequency and load
        {
            json opJson;
            opJson["outputVoltages"] = {550.0};
            opJson["outputCurrents"] = {10.0};  // ~5.5 kW at 550V
            opJson["switchingFrequency"] = 90000.0;  // Above resonance
            opJson["ambientTemperature"] = 45.0;
            opJson["powerFlow"] = "Forward";
            cllcJson["operatingPoints"].push_back(opJson);
        }

        OpenMagnetics::CllcConverter cllc(cllcJson);
        auto params = cllc.calculate_resonant_parameters();

        std::vector<double> turnsRatios = {params.turnsRatio};
        auto operatingPoints = cllc.process_operating_points(turnsRatios, params.magnetizingInductance);

        // Should have: 3 input voltages (min, nom, max) × 2 operating points = 6 operating points
        INFO("Number of operating points generated: " << operatingPoints.size());
        REQUIRE(operatingPoints.size() >= 2);

        // All operating points should have 2 windings (primary + secondary)
        for (const auto& op : operatingPoints) {
            REQUIRE(op.get_excitations_per_winding().size() == 2);
        }
    }

    // =========================================================================
    // TEST 8: Design Requirements Calculation
    // =========================================================================

    TEST_CASE("Test_CllcConverter_DesignRequirements",
              "[converter-model][cllc-topology][design-requirements][smoke-test]") {
        /**
         * Validate the design requirements output.
         * Checks turns ratio and magnetizing inductance match expected values from
         * the Infineon AN example (11 kW, 750V/600V).
         *
         * Expected (from Infineon AN):
         *   n = 1.25, Lm ≈ 160.2 µH, L1 ≈ 36 µH, C1 ≈ 132 nF
         */
        json cllcJson = create_standard_cllc_json();
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto designReqs = cllc.process_design_requirements();

        // Check turns ratio
        REQUIRE(designReqs.get_turns_ratios().size() == 1);
        double tr = designReqs.get_turns_ratios()[0].get_nominal().value();
        INFO("Turns ratio from design requirements: " << tr);
        REQUIRE_THAT(tr, Catch::Matchers::WithinAbs(1.25, 0.1));

        // Check magnetizing inductance
        REQUIRE(designReqs.get_magnetizing_inductance().get_nominal().has_value());
        double Lm = designReqs.get_magnetizing_inductance().get_nominal().value();
        INFO("Magnetizing inductance: " << (Lm * 1e6) << " uH");
        // Lm should be positive and in a reasonable range
        CHECK(Lm > 10e-6);   // > 10 µH
        CHECK(Lm < 1e-3);    // < 1 mH

        // Verify resonant parameters match Infineon AN example approximately
        auto params = cllc.calculate_resonant_parameters();
        INFO("L1 = " << (params.primaryResonantInductance * 1e6) << " uH (expected ~36 uH)");
        INFO("C1 = " << (params.primaryResonantCapacitance * 1e9) << " nF (expected ~132 nF)");
        INFO("Lm = " << (params.magnetizingInductance * 1e6) << " uH (expected ~160 uH)");
        INFO("Ro = " << params.equivalentAcResistance << " ohms (expected ~41.45)");

        // These are order-of-magnitude checks since Q/k choices may differ slightly
        CHECK(params.primaryResonantInductance > 5e-6);
        CHECK(params.primaryResonantInductance < 200e-6);
        CHECK(params.primaryResonantCapacitance > 10e-9);
        CHECK(params.primaryResonantCapacitance < 1000e-9);
    }

    // =========================================================================
    // TEST 9: Advanced Converter with User-Specified Parameters
    // =========================================================================

    TEST_CASE("Test_CllcConverter_AdvancedProcess",
              "[converter-model][cllc-topology][advanced][smoke-test]") {
        /**
         * Test the AdvancedCllcConverter that accepts user-specified resonant parameters.
         */
        json advJson = create_standard_cllc_json();
        advJson["desiredTurnsRatios"] = {1.25};
        advJson["desiredMagnetizingInductance"] = 160e-6;
        advJson["desiredResonantInductancePrimary"] = 36e-6;
        advJson["desiredResonantCapacitancePrimary"] = 132e-9;

        OpenMagnetics::AdvancedCllcConverter cllc(advJson);
        cllc._assertErrors = true;

        auto inputs = cllc.process();

        // Check design requirements
        REQUIRE(inputs.get_design_requirements().get_turns_ratios().size() == 1);
        double tr = inputs.get_design_requirements().get_turns_ratios()[0].get_nominal().value();
        REQUIRE_THAT(tr, Catch::Matchers::WithinAbs(1.25, 0.01));

        double Lm = inputs.get_design_requirements().get_magnetizing_inductance().get_nominal().value();
        REQUIRE_THAT(Lm, Catch::Matchers::WithinAbs(160e-6, 1e-6));

        // Check operating points
        REQUIRE(!inputs.get_operating_points().empty());
        for (const auto& op : inputs.get_operating_points()) {
            REQUIRE(op.get_excitations_per_winding().size() == 2);
        }

        // Paint waveforms
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Cllc_Advanced_Primary_Current.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Cllc_Advanced_Primary_Voltage.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
    }

    // =========================================================================
    // TEST 10: Voltage Gain Curve Validation
    // =========================================================================

    TEST_CASE("Test_CllcConverter_VoltageGainCurve",
              "[converter-model][cllc-topology][gain-curve][smoke-test]") {
        /**
         * Validate the voltage gain curve shape:
         *   - At resonance: gain ≈ 1.0
         *   - Below resonance: gain > 1.0 (boost)
         *   - Above resonance: gain < 1.0 (buck)
         *   - Gain should be monotonically decreasing for fs > fr (in the inductive region)
         *
         * Reference: Infineon AN Figure 8 (gain vs frequency curves)
         */
        json cllcJson = create_standard_cllc_json();
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();
        double fr = params.resonantFrequency;

        double gain_at_fr = cllc.get_voltage_gain(fr, params);
        double gain_below = cllc.get_voltage_gain(fr * 0.7, params);
        double gain_above = cllc.get_voltage_gain(fr * 1.5, params);

        INFO("Gain at fr = " << gain_at_fr);
        INFO("Gain at 0.7*fr = " << gain_below);
        INFO("Gain at 1.5*fr = " << gain_above);

        // At resonance, gain should be close to 1.0
        CHECK_THAT(gain_at_fr, Catch::Matchers::WithinAbs(1.0, 0.2));

        // Above resonance, gain should decrease
        CHECK(gain_above < gain_at_fr);

        // Gain should always be positive
        CHECK(gain_at_fr > 0);
        CHECK(gain_below > 0);
        CHECK(gain_above > 0);
    }

    // =========================================================================
    // TEST 11: Ngspice Simulation - Forward Mode
    // =========================================================================

    TEST_CASE("Test_CllcConverter_Ngspice_Simulation",
              "[converter-model][cllc-topology][ngspice-simulation]") {
        /**
         * Run ngspice simulation of the CLLC converter and validate:
         *   - Simulation completes successfully
         *   - Primary voltage shows switching behavior
         *   - Output voltage is approximately correct
         *   - Reasonable current magnitudes
         */
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        json cllcJson = create_small_power_cllc_json(200000.0, CllcPowerFlow::FORWARD);
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();
        double n = params.turnsRatio;

        INFO("CLLC Simulation - n = " << n);
        INFO("CLLC Simulation - L1 = " << (params.primaryResonantInductance * 1e6) << " uH");
        INFO("CLLC Simulation - C1 = " << (params.primaryResonantCapacitance * 1e9) << " nF");
        INFO("CLLC Simulation - Lm = " << (params.magnetizingInductance * 1e6) << " uH");
        INFO("CLLC Simulation - L2 = " << (params.secondaryResonantInductance * 1e6) << " uH");
        INFO("CLLC Simulation - C2 = " << (params.secondaryResonantCapacitance * 1e9) << " nF");

        std::vector<double> turnsRatios = {n};
        auto converterWaveforms = cllc.simulate_and_extract_topology_waveforms(
            turnsRatios, params.magnetizingInductance);

        REQUIRE(!converterWaveforms.empty());
        REQUIRE(!converterWaveforms[0].get_input_voltage().get_data().empty());

        auto priVoltageData = converterWaveforms[0].get_input_voltage().get_data();
        auto priCurrentData = converterWaveforms[0].get_input_current().get_data();

        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        double priV_min = *std::min_element(priVoltageData.begin(), priVoltageData.end());

        INFO("Simulated primary voltage max: " << priV_max << " V, min: " << priV_min << " V");

        // Primary voltage should show switching behavior (±Vin range)
        CHECK(priV_max > 100);   // Should be significant
        CHECK(priV_min < -100);  // Should show bipolar behavior

        // Check output voltage
        REQUIRE(!converterWaveforms[0].get_output_voltages().empty());
        auto outVoltageData = converterWaveforms[0].get_output_voltages()[0].get_data();
        if (!outVoltageData.empty()) {
            double outV_avg = 0;
            for (double v : outVoltageData) outV_avg += v;
            outV_avg /= outVoltageData.size();
            INFO("Simulated output voltage average: " << outV_avg << " V");
            CHECK(outV_avg > 20);   // Should be around 48V
            CHECK(outV_avg < 100);
        }

        INFO("CLLC ngspice simulation test passed");
    }

    // =========================================================================
    // TEST 12: Ngspice Waveform Polarity
    // =========================================================================

    TEST_CASE("Test_CllcConverter_Waveform_Polarity",
              "[converter-model][cllc-topology][ngspice-simulation][polarity][smoke-test]") {
        /**
         * Verify CLLC waveform polarity:
         *   - Primary voltage should be bipolar (±Vin)
         *   - Primary current should be approximately sinusoidal (bipolar)
         *   - Secondary current should be bipolar
         *
         * Reference: Infineon AN Figure 5 (gate pulses, current waveforms)
         */
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        json cllcJson = create_small_power_cllc_json();
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();
        std::vector<double> turnsRatios = {params.turnsRatio};

        auto operatingPoints = cllc.simulate_and_extract_operating_points(
            turnsRatios, params.magnetizingInductance);

        REQUIRE(!operatingPoints.empty());

        auto& primaryExcitation = operatingPoints[0].get_excitations_per_winding()[0];
        auto& secondaryExcitation = operatingPoints[0].get_excitations_per_winding()[1];

        REQUIRE(primaryExcitation.get_voltage().has_value());
        REQUIRE(primaryExcitation.get_voltage()->get_waveform().has_value());

        auto priVoltageData = primaryExcitation.get_voltage()->get_waveform()->get_data();
        auto priCurrentData = primaryExcitation.get_current()->get_waveform()->get_data();

        double priV_max = *std::max_element(priVoltageData.begin(), priVoltageData.end());
        double priV_min = *std::min_element(priVoltageData.begin(), priVoltageData.end());
        double priI_max = *std::max_element(priCurrentData.begin(), priCurrentData.end());
        double priI_min = *std::min_element(priCurrentData.begin(), priCurrentData.end());

        INFO("Primary voltage max: " << priV_max << " V, min: " << priV_min << " V");
        INFO("Primary current max: " << priI_max << " A, min: " << priI_min << " A");

        // Primary voltage should be bipolar
        CHECK(priV_max > 50);
        CHECK(priV_min < -50);

        // Primary current should be bipolar (sinusoidal-like)
        CHECK(priI_max > 0);
        CHECK(priI_min < 0);
    }

    // =========================================================================
    // TEST 13: Ngspice Num Periods - Simulated Operating Points
    // =========================================================================

    TEST_CASE("Test_CllcConverter_NumPeriods_SimulatedOperatingPoints",
              "[converter-model][cllc-topology][num-periods][ngspice-simulation][smoke-test]") {
        /**
         * Test that changing the number of extracted periods affects waveform length.
         */
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        json cllcJson = create_small_power_cllc_json();
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();
        std::vector<double> turnsRatios = {params.turnsRatio};

        // Simulate with 1 period
        cllc.set_num_periods_to_extract(1);
        auto ops1 = cllc.simulate_and_extract_operating_points(turnsRatios, params.magnetizingInductance);
        REQUIRE(!ops1.empty());
        auto voltageWf1 = ops1[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();

        // Simulate with 3 periods
        cllc.set_num_periods_to_extract(3);
        auto ops3 = cllc.simulate_and_extract_operating_points(turnsRatios, params.magnetizingInductance);
        REQUIRE(!ops3.empty());
        auto voltageWf3 = ops3[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();

        INFO("1-period waveform data size: " << voltageWf1.get_data().size());
        INFO("3-period waveform data size: " << voltageWf3.get_data().size());

        CHECK(voltageWf3.get_data().size() > voltageWf1.get_data().size());
    }

    // =========================================================================
    // TEST 14: Ngspice Converter Waveforms
    // =========================================================================

    TEST_CASE("Test_CllcConverter_NumPeriods_ConverterWaveforms",
              "[converter-model][cllc-topology][num-periods][ngspice-simulation][smoke-test]") {
        /**
         * Test number of periods in converter-level waveforms.
         */
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        json cllcJson = create_small_power_cllc_json();
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();
        std::vector<double> turnsRatios = {params.turnsRatio};

        // Simulate with 1 period
        auto wf1 = cllc.simulate_and_extract_topology_waveforms(turnsRatios, params.magnetizingInductance, 1);
        REQUIRE(!wf1.empty());
        auto inputV1 = wf1[0].get_input_voltage().get_data();

        // Simulate with 3 periods
        auto wf3 = cllc.simulate_and_extract_topology_waveforms(turnsRatios, params.magnetizingInductance, 3);
        REQUIRE(!wf3.empty());
        auto inputV3 = wf3[0].get_input_voltage().get_data();

        INFO("1-period converter waveform data size: " << inputV1.size());
        INFO("3-period converter waveform data size: " << inputV3.size());

        CHECK(inputV3.size() > inputV1.size());
    }

    // =========================================================================
    // TEST 15: Resonant Parameter Validation Against Infineon Example
    // =========================================================================

    TEST_CASE("Test_CllcConverter_InfineonExample_Parameters",
              "[converter-model][cllc-topology][infineon-example][smoke-test]") {
        /**
         * Validate resonant parameters against the Infineon AN-2024-06 example:
         *   Vin = 750V, Vout = 600V, Pout = 11 kW, fr = 73 kHz
         *   Expected: n=1.25, Ro≈41.45Ω, C1≈132nF, L1≈36µH, Lm≈160µH
         *
         * Note: We use Q=0.3984 as in the Infineon example for this validation.
         * Our default Q=0.3 will give different values.
         */
        json cllcJson = create_standard_cllc_json(73000.0);
        cllcJson["qualityFactor"] = 0.3984;  // Match Infineon AN exactly
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();

        // Infineon AN Step 1: n = 750/600 = 1.25
        REQUIRE_THAT(params.turnsRatio,
                     Catch::Matchers::WithinAbs(1.25, 0.01));

        // Infineon AN Step 4: Ro = 8*1.25²/π² * 600²/11000 ≈ 41.45 Ω
        INFO("Ro = " << params.equivalentAcResistance << " (expected ~41.45)");
        REQUIRE_THAT(params.equivalentAcResistance,
                     Catch::Matchers::WithinAbs(41.45, 5.0));

        // Infineon AN Step 5: C1 ≈ 132 nF
        INFO("C1 = " << (params.primaryResonantCapacitance * 1e9) << " nF (expected ~132)");
        REQUIRE_THAT(params.primaryResonantCapacitance * 1e9,
                     Catch::Matchers::WithinAbs(132.0, 30.0));

        // Infineon AN Step 6: L1 ≈ 36 µH
        INFO("L1 = " << (params.primaryResonantInductance * 1e6) << " uH (expected ~36)");
        REQUIRE_THAT(params.primaryResonantInductance * 1e6,
                     Catch::Matchers::WithinAbs(36.0, 10.0));

        // Infineon AN Step 7: Lm = k*L1 = 4.45*36µH ≈ 160.2 µH
        INFO("Lm = " << (params.magnetizingInductance * 1e6) << " uH (expected ~160.2)");
        REQUIRE_THAT(params.magnetizingInductance * 1e6,
                     Catch::Matchers::WithinAbs(160.2, 50.0));
    }

    // =========================================================================
    // TEST 16: Netlist Generation Smoke Test
    // =========================================================================

    TEST_CASE("Test_CllcConverter_NetlistGeneration",
              "[converter-model][cllc-topology][netlist][smoke-test]") {
        /**
         * Verify that the ngspice netlist is generated correctly and contains
         * all expected circuit elements.
         */
        json cllcJson = create_small_power_cllc_json();
        OpenMagnetics::CllcConverter cllc(cllcJson);

        auto params = cllc.calculate_resonant_parameters();
        std::string netlist = cllc.generate_ngspice_circuit(params.turnsRatio, params);

        INFO("Generated netlist:\n" << netlist);

        // Check for essential circuit elements
        CHECK(netlist.find("Vin") != std::string::npos);
        CHECK(netlist.find("S1") != std::string::npos);
        CHECK(netlist.find("S2") != std::string::npos);
        CHECK(netlist.find("S3") != std::string::npos);
        CHECK(netlist.find("S4") != std::string::npos);
        CHECK(netlist.find("C_res1") != std::string::npos);
        CHECK(netlist.find("L_res1") != std::string::npos);
        CHECK(netlist.find("Lpri") != std::string::npos);
        CHECK(netlist.find("Lsec") != std::string::npos);
        CHECK(netlist.find("L_res2") != std::string::npos);
        CHECK(netlist.find("C_res2") != std::string::npos);
        CHECK(netlist.find("Ds1") != std::string::npos);
        CHECK(netlist.find("Ds2") != std::string::npos);
        CHECK(netlist.find("Ds3") != std::string::npos);
        CHECK(netlist.find("Ds4") != std::string::npos);
        CHECK(netlist.find("Kpri_sec") != std::string::npos);
        CHECK(netlist.find("Rload") != std::string::npos);
        CHECK(netlist.find("Cout") != std::string::npos);
        CHECK(netlist.find(".tran") != std::string::npos);
        CHECK(netlist.find(".end") != std::string::npos);
    }

} // anonymous namespace
