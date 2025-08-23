#include "support/Painter.h"
#define private public
#define protected public
#include "converter_models/Topology.h"
#undef private
#undef protected
#include "json.hpp"
#include <UnitTest++.h>
#include <filesystem>
#include <vector>
#include <cmath>

using namespace MAS;
using namespace OpenMagnetics;
using json = nlohmann::json;

SUITE(TwoLevelInverter) {
    auto outputFilePath = std::filesystem::path{__FILE__}.parent_path().append("..").append("output");

    TEST(Test_Inverter_SPWM_3phase) {
        std::filesystem::create_directories(outputFilePath);

        json invJson;
        invJson["dcBusVoltage"]["nominal"] = 400.0;
        invJson["dcBusCapacitor"]["capacitance"] = 1e-3;
        invJson["numberOfLegs"] = 3; // single-phase
        invJson["inverterLegPowerRating"] = 1000.0;
        invJson["lineRmsCurrent"]["nominal"] = 5.0;
        invJson["modulation"]["switchingFrequency"] = 10000.0;
        invJson["modulation"]["pwmType"] = "triangular";
        invJson["modulation"]["modulationStrategy"] = "SPWM";
        invJson["modulation"]["modulationDepth"] = 0.5;
        json op;
        op["fundamentalFrequency"] = 50.0;
        op["outputPower"] = 100.0;
        op["powerFactor"] = 1.0;
        op["load"]["loadType"] = "grid";
        op["load"]["gridFrequency"] = 50;
        op["load"]["gridResistance"] = 0.0001;
        op["load"]["gridInductance"] = 1e-7;
        op["currentPhaseAngle"] = 0.0;
        invJson["operatingPoints"] = json::array({op});
        invJson["downstreamFilter"]["filterTopology"] = "L";
        invJson["downstreamFilter"]["inductor"]["desiredInductance"]["nominal"] = 1e-3; // 1 mH
        invJson["downstreamFilter"]["inductor"]["resistance"] = 0.01; // 10 mΩ

        // Construct inverter
        MyInverter inverter(invJson);

        // Call normal workflow → this internally computes harmonics + plots if DEBUG_PLOTS
        auto results = inverter.process_operating_points();

        // Just check results exist
        CHECK(!results.empty());

        // Check that debug plots were generated
    #ifdef DEBUG_PLOTS
        CHECK(std::filesystem::exists("debug_plots/carrier_vs_refs.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_short.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_fundamental.png"));
        CHECK(std::filesystem::exists("debug_plots/power.png"));
        CHECK(std::filesystem::exists("debug_plots/vdc_ripple.png"));
        CHECK(std::filesystem::exists("debug_plots/fft_vl1_il1.png"));
        CHECK(std::filesystem::exists("debug_plots/final_fft_vl1_il1.png"));
    #endif
    }
    TEST(Test_Inverter_SVPWM_3phase) {
        std::filesystem::create_directories(outputFilePath);

        json invJson;
        invJson["dcBusVoltage"]["nominal"] = 400.0;
        invJson["dcBusCapacitor"]["capacitance"] = 1e-3;
        invJson["numberOfLegs"] = 3; // single-phase
        invJson["inverterLegPowerRating"] = 1000.0;
        invJson["lineRmsCurrent"]["nominal"] = 5.0;
        invJson["modulation"]["switchingFrequency"] = 10000.0;
        invJson["modulation"]["pwmType"] = "triangular";
        invJson["modulation"]["modulationStrategy"] = "SVPWM";
        invJson["modulation"]["modulationDepth"] = 0.5;
        json op;
        op["fundamentalFrequency"] = 50.0;
        op["outputPower"] = 100.0;
        op["powerFactor"] = 1.0;
        op["load"]["loadType"] = "grid";
        op["load"]["gridFrequency"] = 50;
        op["load"]["gridResistance"] = 0.0001;
        op["load"]["gridInductance"] = 1e-7;
        op["currentPhaseAngle"] = 0.0;
        invJson["operatingPoints"] = json::array({op});
        invJson["downstreamFilter"]["filterTopology"] = "L";
        invJson["downstreamFilter"]["inductor"]["desiredInductance"]["nominal"] = 1e-3; // 1 mH
        invJson["downstreamFilter"]["inductor"]["resistance"] = 0.01; // 10 mΩ

        // Construct inverter
        MyInverter inverter(invJson);

        // Call normal workflow → this internally computes harmonics + plots if DEBUG_PLOTS
        auto results = inverter.process_operating_points();

        // Just check results exist
        CHECK(!results.empty());

        // Check that debug plots were generated
    #ifdef DEBUG_PLOTS
        CHECK(std::filesystem::exists("debug_plots/carrier_vs_refs.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_short.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_fundamental.png"));
        CHECK(std::filesystem::exists("debug_plots/power.png"));
        CHECK(std::filesystem::exists("debug_plots/vdc_ripple.png"));
        CHECK(std::filesystem::exists("debug_plots/fft_vl1_il1.png"));
        CHECK(std::filesystem::exists("debug_plots/final_fft_vl1_il1.png"));
    #endif
    }
    TEST(Test_Inverter_THIPWM_3phase) {
        std::filesystem::create_directories(outputFilePath);

        json invJson;
        invJson["dcBusVoltage"]["nominal"] = 400.0;
        invJson["dcBusCapacitor"]["capacitance"] = 1e-3;
        invJson["numberOfLegs"] = 3; // single-phase
        invJson["inverterLegPowerRating"] = 1000.0;
        invJson["lineRmsCurrent"]["nominal"] = 5.0;
        invJson["modulation"]["switchingFrequency"] = 10000.0;
        invJson["modulation"]["pwmType"] = "triangular";
        invJson["modulation"]["modulationStrategy"] = "THIPWM";
        invJson["modulation"]["thirdHarmonicInjectionCoefficient"] = 0.166;
        invJson["modulation"]["modulationDepth"] = 0.5;
        json op;
        op["fundamentalFrequency"] = 50.0;
        op["outputPower"] = 100.0;
        op["powerFactor"] = 1.0;
        op["load"]["loadType"] = "grid";
        op["load"]["gridFrequency"] = 50;
        op["load"]["gridResistance"] = 0.0001;
        op["load"]["gridInductance"] = 1e-7;
        op["currentPhaseAngle"] = 0.0;
        invJson["operatingPoints"] = json::array({op});
        invJson["downstreamFilter"]["filterTopology"] = "L";
        invJson["downstreamFilter"]["inductor"]["desiredInductance"]["nominal"] = 1e-3; // 1 mH
        invJson["downstreamFilter"]["inductor"]["resistance"] = 0.01; // 10 mΩ

        // Construct inverter
        MyInverter inverter(invJson);

        // Call normal workflow → this internally computes harmonics + plots if DEBUG_PLOTS
        auto results = inverter.process_operating_points();

        // Just check results exist
        CHECK(!results.empty());

        // Check that debug plots were generated
    #ifdef DEBUG_PLOTS
        CHECK(std::filesystem::exists("debug_plots/carrier_vs_refs.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_short.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_fundamental.png"));
        CHECK(std::filesystem::exists("debug_plots/power.png"));
        CHECK(std::filesystem::exists("debug_plots/vdc_ripple.png"));
        CHECK(std::filesystem::exists("debug_plots/fft_vl1_il1.png"));
        CHECK(std::filesystem::exists("debug_plots/final_fft_vl1_il1.png"));
    #endif
    }
    TEST(Test_Inverter_SPWM_2phase) {
        std::filesystem::create_directories(outputFilePath);

        json invJson;
        invJson["dcBusVoltage"]["nominal"] = 400.0;
        invJson["dcBusCapacitor"]["capacitance"] = 1e-3;
        invJson["numberOfLegs"] = 2; // single-phase
        invJson["inverterLegPowerRating"] = 1000.0;
        invJson["lineRmsCurrent"]["nominal"] = 5.0;
        invJson["modulation"]["switchingFrequency"] = 10000.0;
        invJson["modulation"]["pwmType"] = "triangular";
        invJson["modulation"]["modulationStrategy"] = "SPWM";
        invJson["modulation"]["modulationDepth"] = 0.5;
        json op;
        op["fundamentalFrequency"] = 50.0;
        op["outputPower"] = 100.0;
        op["powerFactor"] = 1.0;
        op["load"]["loadType"] = "grid";
        op["load"]["gridFrequency"] = 50;
        op["load"]["gridResistance"] = 0.0001;
        op["load"]["gridInductance"] = 1e-7;
        op["currentPhaseAngle"] = 0.0;
        invJson["operatingPoints"] = json::array({op});
        invJson["downstreamFilter"]["filterTopology"] = "L";
        invJson["downstreamFilter"]["inductor"]["desiredInductance"]["nominal"] = 1e-3; // 1 mH
        invJson["downstreamFilter"]["inductor"]["resistance"] = 0.01; // 10 mΩ

        // Construct inverter
        MyInverter inverter(invJson);

        // Call normal workflow → this internally computes harmonics + plots if DEBUG_PLOTS
        auto results = inverter.process_operating_points();

        // Just check results exist
        CHECK(!results.empty());

        // Check that debug plots were generated
    #ifdef DEBUG_PLOTS
        CHECK(std::filesystem::exists("debug_plots/carrier_vs_refs.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_short.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_fundamental.png"));
        CHECK(std::filesystem::exists("debug_plots/power.png"));
        CHECK(std::filesystem::exists("debug_plots/vdc_ripple.png"));
        CHECK(std::filesystem::exists("debug_plots/fft_vl1_il1.png"));
        CHECK(std::filesystem::exists("debug_plots/final_fft_vl1_il1.png"));
    #endif
    }
    TEST(Test_Inverter_SVPWM_2phase) {
        std::filesystem::create_directories(outputFilePath);

        json invJson;
        invJson["dcBusVoltage"]["nominal"] = 400.0;
        invJson["dcBusCapacitor"]["capacitance"] = 1e-3;
        invJson["numberOfLegs"] = 2; // single-phase
        invJson["inverterLegPowerRating"] = 1000.0;
        invJson["lineRmsCurrent"]["nominal"] = 5.0;
        invJson["modulation"]["switchingFrequency"] = 10000.0;
        invJson["modulation"]["pwmType"] = "triangular";
        invJson["modulation"]["modulationStrategy"] = "SVPWM";
        invJson["modulation"]["modulationDepth"] = 0.5;
        json op;
        op["fundamentalFrequency"] = 50.0;
        op["outputPower"] = 100.0;
        op["powerFactor"] = 1.0;
        op["load"]["loadType"] = "grid";
        op["load"]["gridFrequency"] = 50;
        op["load"]["gridResistance"] = 0.0001;
        op["load"]["gridInductance"] = 1e-7;
        op["currentPhaseAngle"] = 0.0;
        invJson["operatingPoints"] = json::array({op});
        invJson["downstreamFilter"]["filterTopology"] = "L";
        invJson["downstreamFilter"]["inductor"]["desiredInductance"]["nominal"] = 1e-3; // 1 mH
        invJson["downstreamFilter"]["inductor"]["resistance"] = 0.01; // 10 mΩ

        // Construct inverter
        MyInverter inverter(invJson);

        // Call normal workflow → this internally computes harmonics + plots if DEBUG_PLOTS
        auto results = inverter.process_operating_points();

        // Just check results exist
        CHECK(!results.empty());

        // Check that debug plots were generated
    #ifdef DEBUG_PLOTS
        CHECK(std::filesystem::exists("debug_plots/carrier_vs_refs.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_short.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_fundamental.png"));
        CHECK(std::filesystem::exists("debug_plots/power.png"));
        CHECK(std::filesystem::exists("debug_plots/vdc_ripple.png"));
        CHECK(std::filesystem::exists("debug_plots/fft_vl1_il1.png"));
        CHECK(std::filesystem::exists("debug_plots/final_fft_vl1_il1.png"));
    #endif
    }
    TEST(Test_Inverter_THIPWM_2phase) {
        std::filesystem::create_directories(outputFilePath);

        json invJson;
        invJson["dcBusVoltage"]["nominal"] = 400.0;
        invJson["dcBusCapacitor"]["capacitance"] = 1e-3;
        invJson["numberOfLegs"] = 2; // single-phase
        invJson["inverterLegPowerRating"] = 1000.0;
        invJson["lineRmsCurrent"]["nominal"] = 5.0;
        invJson["modulation"]["switchingFrequency"] = 10000.0;
        invJson["modulation"]["pwmType"] = "triangular";
        invJson["modulation"]["modulationStrategy"] = "THIPWM";
        invJson["modulation"]["modulationDepth"] = 0.5;
        json op;
        op["fundamentalFrequency"] = 50.0;
        op["outputPower"] = 100.0;
        op["powerFactor"] = 1.0;
        op["load"]["loadType"] = "grid";
        op["load"]["gridFrequency"] = 50;
        op["load"]["gridResistance"] = 0.0001;
        op["load"]["gridInductance"] = 1e-7;
        op["currentPhaseAngle"] = 0.0;
        invJson["operatingPoints"] = json::array({op});
        invJson["downstreamFilter"]["filterTopology"] = "L";
        invJson["downstreamFilter"]["inductor"]["desiredInductance"]["nominal"] = 1e-3; // 1 mH
        invJson["downstreamFilter"]["inductor"]["resistance"] = 0.01; // 10 mΩ

        // Construct inverter
        MyInverter inverter(invJson);

        // Call normal workflow → this internally computes harmonics + plots if DEBUG_PLOTS
        auto results = inverter.process_operating_points();

        // Just check results exist
        CHECK(!results.empty());

        // Check that debug plots were generated
    #ifdef DEBUG_PLOTS
        CHECK(std::filesystem::exists("debug_plots/carrier_vs_refs.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_short.png"));
        CHECK(std::filesystem::exists("debug_plots/va_vb_vc_fundamental.png"));
        CHECK(std::filesystem::exists("debug_plots/power.png"));
        CHECK(std::filesystem::exists("debug_plots/vdc_ripple.png"));
        CHECK(std::filesystem::exists("debug_plots/fft_vl1_il1.png"));
        CHECK(std::filesystem::exists("debug_plots/final_fft_vl1_il1.png"));
    #endif
    }
}