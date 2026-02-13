#include <source_location>
#include "support/Painter.h"
#include "converter_models/PhaseShiftedFullBridge.h"
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

using namespace OpenMagnetics;

namespace {
auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}
    .parent_path().append("..").append("output");
double maximumError = 0.1;


// =========================================================================
// Helper: create a typical PSFB JSON (400V -> 12V, 50A, 100 kHz, center-tapped)
// =========================================================================
json make_psfb_json(double Vin_nom = 400.0, double Vin_min = 370.0, double Vin_max = 410.0,
                    double Vo = 12.0, double Io = 50.0, double Fs = 100000.0,
                    double phaseShift = 126.0,
                    std::string rectType = "Center Tapped") {
    json j;
    json inputVoltage;
    inputVoltage["nominal"] = Vin_nom;
    inputVoltage["minimum"] = Vin_min;
    inputVoltage["maximum"] = Vin_max;
    j["inputVoltage"] = inputVoltage;
    j["rectifierType"] = rectType;

    j["operatingPoints"] = json::array();
    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {Vo};
        op["outputCurrents"] = {Io};
        op["switchingFrequency"] = Fs;
        op["phaseShift"] = phaseShift;
        j["operatingPoints"].push_back(op);
    }
    return j;
}


// =========================================================================
// TEST 1: PSFB Design - Center Tapped, 400V -> 12V, 600W
// =========================================================================
TEST_CASE("Test_Psfb_CenterTapped_Design", "[converter-model][psfb-topology][smoke-test]") {
    auto psfbJson = make_psfb_json();
    OpenMagnetics::Psfb psfb(psfbJson);

    SECTION("Input validation") {
        CHECK(psfb.run_checks(false) == true);
    }

    SECTION("Turns ratio") {
        auto req = psfb.process_design_requirements();
        double n = resolve_dimensional_values(req.get_turns_ratios()[0]);
        // n = Vin * Deff / (Vo + Vd)
        // Deff = 126/180 = 0.7
        // n = 400 * 0.7 / (12 + 0.6) = 280 / 12.6 ≈ 22.22
        double expectedN = 400.0 * (126.0 / 180.0) / (12.0 + 0.6);
        REQUIRE_THAT(n, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.05));
        CHECK(n > 0);
    }

    SECTION("Magnetizing inductance is positive and reasonable") {
        auto req = psfb.process_design_requirements();
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
        CHECK(Lm > 0);
        CHECK(Lm > 10e-6);
        CHECK(Lm < 500e-3);
    }

    SECTION("Computed series inductance is positive") {
        psfb.process_design_requirements();
        double Lr = psfb.get_computed_series_inductance();
        CHECK(Lr > 0);
        CHECK(Lr < 1e-3);
    }

    SECTION("Computed output inductance is positive") {
        psfb.process_design_requirements();
        double Lo = psfb.get_computed_output_inductance();
        CHECK(Lo > 0);
    }

    SECTION("Effective duty cycle") {
        psfb.process_design_requirements();
        double Deff = psfb.get_computed_effective_duty_cycle();
        CHECK(Deff > 0);
        CHECK(Deff < 1.0);
        REQUIRE_THAT(Deff, Catch::Matchers::WithinAbs(0.7, 0.01));
    }
}


// =========================================================================
// TEST 2: PSFB Operating Points Generation
// =========================================================================
TEST_CASE("Test_Psfb_OperatingPoints_Generation", "[converter-model][psfb-topology][smoke-test]") {
    auto psfbJson = make_psfb_json();
    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    SECTION("Multiple input voltages") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        // Vin has nominal, min, max → 3 operating points
        CHECK(ops.size() == 3);
    }

    SECTION("Waveform structure") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        REQUIRE(!ops.empty());

        auto& op = ops[0];
        // Primary + 1 secondary = 2 windings
        CHECK(op.get_excitations_per_winding().size() == 2);

        auto& priExc = op.get_excitations_per_winding()[0];
        REQUIRE(priExc.get_current().has_value());
        REQUIRE(priExc.get_voltage().has_value());
        CHECK(priExc.get_frequency() == 100e3);

        // Waveforms should have 2*256+1 = 513 samples
        auto currentWfm = priExc.get_current()->get_waveform().value();
        CHECK(currentWfm.get_data().size() == 513);

        auto voltageWfm = priExc.get_voltage()->get_waveform().value();
        CHECK(voltageWfm.get_data().size() == 513);
    }

    SECTION("Primary voltage is 3-level") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        auto& priExc = ops[0].get_excitations_per_winding()[0];
        auto vData = priExc.get_voltage()->get_waveform().value().get_data();

        double vMax = *std::max_element(vData.begin(), vData.end());
        double vMin = *std::min_element(vData.begin(), vData.end());

        // Primary voltage should swing to approximately +Vin and -Vin
        double Vin_min = 370.0;
        REQUIRE_THAT(vMax, Catch::Matchers::WithinAbs(Vin_min, Vin_min * 0.05));
        REQUIRE_THAT(vMin, Catch::Matchers::WithinAbs(-Vin_min, Vin_min * 0.05));

        // Should have zero-voltage intervals (freewheeling)
        int zeroCount = 0;
        for (double v : vData) {
            if (std::abs(v) < 1.0) zeroCount++;
        }
        CHECK(zeroCount > 0);
    }

    SECTION("Primary current antisymmetry") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        auto& priExc = ops[0].get_excitations_per_winding()[0];
        auto iData = priExc.get_current()->get_waveform().value().get_data();

        // PSFB primary current should be antisymmetric: i(t + T/2) = -i(t)
        size_t N = iData.size();
        size_t half = N / 2;
        double asymmetrySum = 0;
        for (size_t k = 1; k < half; k++) {
            asymmetrySum += std::abs(iData[k] + iData[half + k]);
        }
        double avgAsymmetry = asymmetrySum / half;
        double iPeak = *std::max_element(iData.begin(), iData.end());
        // Antisymmetry error should be small relative to peak
        CHECK(avgAsymmetry / iPeak < 0.05);
    }

    SECTION("Secondary winding excitation exists") {
        auto ops = psfb.process_operating_points(turnsRatios, Lm);
        auto& secExc = ops[0].get_excitations_per_winding()[1];
        REQUIRE(secExc.get_current().has_value());
        REQUIRE(secExc.get_voltage().has_value());
    }
}


// =========================================================================
// TEST 3: PSFB Waveform Plotting
// =========================================================================
TEST_CASE("Test_Psfb_Waveform_Plotting", "[converter-model][psfb-topology][visualization]") {
    auto psfbJson = make_psfb_json();
    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto ops = psfb.process_operating_points(turnsRatios, Lm);
    REQUIRE(!ops.empty());

    SECTION("Primary current waveform plot") {
        auto outFile = outputFilePath;
        outFile.append("Test_Psfb_Primary_Current_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
        painter.export_svg();
    }

    SECTION("Primary voltage waveform plot") {
        auto outFile = outputFilePath;
        outFile.append("Test_Psfb_Primary_Voltage_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
        painter.export_svg();
    }

    SECTION("Secondary current waveform plot") {
        auto outFile = outputFilePath;
        outFile.append("Test_Psfb_Secondary_Current_Waveform.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile, false, true);
        painter.paint_waveform(ops[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
        painter.export_svg();
    }
}


// =========================================================================
// TEST 4: PSFB SPICE Netlist Generation
// =========================================================================
TEST_CASE("Test_Psfb_Spice_Netlist", "[converter-model][psfb-topology][spice]") {
    auto psfbJson = make_psfb_json();
    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    SECTION("Netlist is non-empty and contains key elements") {
        std::string netlist = psfb.generate_ngspice_circuit(turnsRatios, Lm);
        CHECK(!netlist.empty());
        CHECK(netlist.find("Phase-Shifted Full Bridge") != std::string::npos);
        CHECK(netlist.find("L_pri") != std::string::npos);
        CHECK(netlist.find("L_sec") != std::string::npos);
        CHECK(netlist.find("K_trafo") != std::string::npos);
        CHECK(netlist.find("L_out") != std::string::npos);
        CHECK(netlist.find("R_load") != std::string::npos);
        CHECK(netlist.find(".tran") != std::string::npos);
        // Full bridge: 4 switches (SA, SB, SC, SD)
        CHECK(netlist.find("SA ") != std::string::npos);
        CHECK(netlist.find("SB ") != std::string::npos);
        CHECK(netlist.find("SC ") != std::string::npos);
        CHECK(netlist.find("SD ") != std::string::npos);
    }

    SECTION("Netlist saved to file") {
        std::string netlist = psfb.generate_ngspice_circuit(turnsRatios, Lm);
        auto outFile = outputFilePath;
        outFile.append("Test_Psfb_Netlist.cir");
        std::ofstream ofs(outFile);
        ofs << netlist;
        ofs.close();
        CHECK(std::filesystem::exists(outFile));
    }
}


// =========================================================================
// TEST 5: PSFB Multiple Outputs
// =========================================================================
TEST_CASE("Test_Psfb_Multiple_Outputs", "[converter-model][psfb-topology][smoke-test]") {
    json psfbJson;
    json inputVoltage;
    inputVoltage["nominal"] = 400.0;
    psfbJson["inputVoltage"] = inputVoltage;
    psfbJson["rectifierType"] = "Center Tapped";
    psfbJson["operatingPoints"] = json::array();

    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0, 5.0};
        op["outputCurrents"] = {50.0, 10.0};
        op["switchingFrequency"] = 100000;
        op["phaseShift"] = 126.0;
        psfbJson["operatingPoints"].push_back(op);
    }

    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    CHECK(req.get_turns_ratios().size() == 2);

    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());

    auto ops = psfb.process_operating_points(turnsRatios, Lm);
    REQUIRE(!ops.empty());

    // Primary + 2 secondaries = 3 windings
    CHECK(ops[0].get_excitations_per_winding().size() == 3);
}


// =========================================================================
// TEST 6: PSFB Static calculations
// =========================================================================
TEST_CASE("Test_Psfb_Static_Calculations", "[converter-model][psfb-topology][unit]") {

    SECTION("Effective duty cycle") {
        CHECK(Psfb::compute_effective_duty_cycle(0.0) == 0.0);
        CHECK(Psfb::compute_effective_duty_cycle(90.0) == 0.5);
        CHECK(Psfb::compute_effective_duty_cycle(180.0) == 1.0);
        REQUIRE_THAT(Psfb::compute_effective_duty_cycle(126.0),
                     Catch::Matchers::WithinAbs(0.7, 1e-6));
    }

    SECTION("Output voltage - center tapped") {
        // Vo = Vin * Deff / n - Vd
        double Vo = Psfb::compute_output_voltage(400.0, 0.7, 22.0, 0.6,
                                                  PsfbRectifierType::CENTER_TAPPED);
        double expected = 400.0 * 0.7 / 22.0 - 0.6;
        REQUIRE_THAT(Vo, Catch::Matchers::WithinAbs(expected, 0.01));
    }

    SECTION("Turns ratio round-trip") {
        double n = Psfb::compute_turns_ratio(400.0, 12.0, 0.7, 0.6,
                                              PsfbRectifierType::CENTER_TAPPED);
        double Vo_check = Psfb::compute_output_voltage(400.0, 0.7, n, 0.6,
                                                        PsfbRectifierType::CENTER_TAPPED);
        REQUIRE_THAT(Vo_check, Catch::Matchers::WithinAbs(12.0, 0.01));
    }
}


// =========================================================================
// TEST 7: PSFB Current Doubler rectifier
// =========================================================================
TEST_CASE("Test_Psfb_CurrentDoubler_Design", "[converter-model][psfb-topology][smoke-test]") {
    auto psfbJson = make_psfb_json(400.0, 370.0, 410.0, 12.0, 50.0, 100000.0,
                                   126.0, "Current Doubler");
    OpenMagnetics::Psfb psfb(psfbJson);
    auto req = psfb.process_design_requirements();

    double n = resolve_dimensional_values(req.get_turns_ratios()[0]);
    // Current doubler: n = Vin * Deff / (2*(Vo+Vd))
    double expectedN = 400.0 * 0.7 / (2.0 * (12.0 + 0.6));
    REQUIRE_THAT(n, Catch::Matchers::WithinAbs(expectedN, expectedN * 0.05));
}


// =========================================================================
// TEST 8: AdvancedPsfb JSON round-trip
// =========================================================================
TEST_CASE("Test_AdvancedPsfb_Process", "[converter-model][psfb-topology][advanced]") {
    json advJson;
    json inputVoltage;
    inputVoltage["nominal"] = 400.0;
    inputVoltage["minimum"] = 370.0;
    inputVoltage["maximum"] = 410.0;
    advJson["inputVoltage"] = inputVoltage;
    advJson["rectifierType"] = "Center Tapped";
    advJson["operatingPoints"] = json::array();

    {
        json op;
        op["ambientTemperature"] = 25.0;
        op["outputVoltages"] = {12.0};
        op["outputCurrents"] = {50.0};
        op["switchingFrequency"] = 100000;
        op["phaseShift"] = 126.0;
        advJson["operatingPoints"].push_back(op);
    }

    advJson["desiredTurnsRatios"] = {22.0};
    advJson["desiredMagnetizingInductance"] = 1e-3;

    OpenMagnetics::AdvancedPsfb advPsfb(advJson);
    auto inputs = advPsfb.process();

    CHECK(inputs.get_design_requirements().get_turns_ratios().size() == 1);
    double n = resolve_dimensional_values(inputs.get_design_requirements().get_turns_ratios()[0]);
    REQUIRE_THAT(n, Catch::Matchers::WithinAbs(22.0, 0.01));

    double Lm = resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance());
    REQUIRE_THAT(Lm, Catch::Matchers::WithinAbs(1e-3, 1e-6));

    CHECK(!inputs.get_operating_points().empty());
}

} // anonymous namespace
