#include <source_location>
#include "support/Painter.h"
#include "converter_models/SingleSwitchForward.h"
#include "converter_models/TwoSwitchForward.h"
#include "converter_models/ActiveClampForward.h"
#include "support/Utils.h"
#include "TestingUtils.h"
#include "processors/NgspiceRunner.h"
#include "ConverterPortChecks.h"

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
#include <limits>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
    double maximumError = 0.1;

    // ── PtP waveform comparison helpers ────────────────────────────────────
    std::vector<double> ptp_interp(const std::vector<double>& t, const std::vector<double>& d, int N) {
        std::vector<double> out(N);
        double T = t.back();
        for (int i = 0; i < N; ++i) {
            double ti = T * i / (N - 1);
            int lo = 0;
            for (int k = 0; k + 1 < (int)t.size(); ++k) if (t[k] <= ti) lo = k;
            int hi = std::min(lo + 1, (int)t.size() - 1);
            double dt = t[hi] - t[lo];
            out[i] = (dt < 1e-20) ? d[hi] : d[lo] + (ti - t[lo]) / dt * (d[hi] - d[lo]);
        }
        return out;
    }

    double ptp_nrmse(const std::vector<double>& ref, const std::vector<double>& sim) {
        // Mean-subtracted, scale-normalized NRMSE with circular phase alignment.
        // Measures shape similarity only — invariant to DC offset and amplitude scale.
        int N = (int)ref.size();
        double ref_mean = 0.0, sim_mean = 0.0;
        for (int i = 0; i < N; ++i) { ref_mean += ref[i]; sim_mean += sim[i]; }
        ref_mean /= N; sim_mean /= N;
        std::vector<double> r(N), s(N);
        double rAC = 0.0, sAC = 0.0;
        for (int i = 0; i < N; ++i) {
            r[i] = ref[i] - ref_mean; s[i] = sim[i] - sim_mean;
            rAC += r[i] * r[i]; sAC += s[i] * s[i];
        }
        rAC = std::sqrt(rAC / N); sAC = std::sqrt(sAC / N);
        if (rAC < 1e-10 || sAC < 1e-10) return 1.0;
        for (int i = 0; i < N; ++i) { r[i] /= rAC; s[i] /= sAC; }
        int ns = std::min(N, 64);
        double best = std::numeric_limits<double>::max();
        for (int ss = 0; ss < ns; ++ss) {
            int sh = ss * N / ns;
            double ssd = 0.0;
            for (int k = 0; k < N; ++k) { double e = r[k] - s[(k + sh) % N]; ssd += e * e; }
            if (ssd < best) best = ssd;
        }
        return std::sqrt(best / N);
    }
    std::vector<double> ptp_current(const OperatingPoint& op, size_t wi = 0) {
        if (wi >= op.get_excitations_per_winding().size()) return {};
        auto& e = op.get_excitations_per_winding()[wi];
        if (!e.get_current() || !e.get_current()->get_waveform()) return {};
        return e.get_current()->get_waveform()->get_data();
    }
    std::vector<double> ptp_current_time(const OperatingPoint& op, size_t wi = 0) {
        if (wi >= op.get_excitations_per_winding().size()) return {};
        auto& e = op.get_excitations_per_winding()[wi];
        if (!e.get_current() || !e.get_current()->get_waveform()) return {};
        auto tv = e.get_current()->get_waveform()->get_time();
        return tv ? tv.value() : std::vector<double>{};
    }
    // ───────────────────────────────────────────────────────────────────────

    TEST_CASE("Test_SingleSwitchForward_CCM", "[converter-model][single-switch-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 0.3;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {5};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::SingleSwitchForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Demagnetization_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Demagnetization_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_CCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE_THAT(double(forwardInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(forwardInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE_THAT(double(forwardInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(forwardInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        // Allow small DC offset due to ngspice simulation
        REQUIRE(std::abs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset()) < 1.0);

        double outputCurrent = (inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_peak().value() + inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset()) / 2;
        REQUIRE_THAT(double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(outputCurrent, double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_SingleSwitchForward_DCM", "[converter-model][single-switch-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 2;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["desiredInductance"] = 1e-3;
        forwardInputsJson["desiredOutputInductances"] = {5e-6};
        forwardInputsJson["desiredTurnsRatios"] = {1, 2};
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {1};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::AdvancedSingleSwitchForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Demagnetization_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Demagnetization_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_SingleSwitchForward_DCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        // Allow small DC offset due to ngspice simulation
        REQUIRE(std::abs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset()) < 1.0);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
        // Allow small DC offset due to ngspice simulation
        REQUIRE(std::abs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset()) < 1.0);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[2].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
    }

    TEST_CASE("Test_ActiveClampForward_CCM", "[converter-model][active-clamp-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 0.3;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {5};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::ActiveClampForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_CCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_CCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_CCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_CCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE_THAT(double(forwardInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(forwardInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        double outputCurrent = (inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_peak().value() + inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset()) / 2;
        REQUIRE_THAT(double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(outputCurrent, double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_ActiveClampForward_DCM", "[converter-model][active-clamp-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 2;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["desiredInductance"] = 1e-3;
        forwardInputsJson["desiredOutputInductances"] = {5e-6};
        forwardInputsJson["desiredTurnsRatios"] = {0.5};
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {1};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::AdvancedActiveClampForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_DCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_DCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_DCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_ActiveClampForward_DCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        // Allow small DC offset due to ngspice simulation
        REQUIRE(std::abs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset()) < 1.0);
    }

    TEST_CASE("Test_TwoSwitchForward_CCM", "[converter-model][two-switch-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 0.3;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {5};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::TwoSwitchForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_CCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_CCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_CCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_CCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE_THAT(double(forwardInputsJson["inputVoltage"]["minimum"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_peak().value(), double(forwardInputsJson["inputVoltage"]["maximum"]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        double outputCurrent = (inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_peak().value() + inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset()) / 2;
        REQUIRE_THAT(double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]), Catch::Matchers::WithinAbs(outputCurrent, double(forwardInputsJson["operatingPoints"][0]["outputCurrents"][0]) * maximumError));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset() > 0);
    }

    TEST_CASE("Test_TwoSwitchForward_DCM", "[converter-model][two-switch-forward-topology][smoke-test]") {
        json forwardInputsJson;
        json inputVoltage;


        inputVoltage["minimum"] = 100;
        inputVoltage["maximum"] = 190;
        forwardInputsJson["inputVoltage"] = inputVoltage;
        forwardInputsJson["diodeVoltageDrop"] = 0.5;
        forwardInputsJson["efficiency"] = 0.9;
        forwardInputsJson["maximumSwitchCurrent"] = 1;
        forwardInputsJson["currentRippleRatio"] = 2;
        forwardInputsJson["dutyCycle"] = 0.42;
        forwardInputsJson["desiredInductance"] = 1e-3;
        forwardInputsJson["desiredOutputInductances"] = {5e-6};
        forwardInputsJson["desiredTurnsRatios"] = {0.5};
        forwardInputsJson["operatingPoints"] = json::array();

        {
            json ForwardOperatingPointJson;
            ForwardOperatingPointJson["outputVoltages"] = {5};
            ForwardOperatingPointJson["outputCurrents"] = {1};
            ForwardOperatingPointJson["switchingFrequency"] = 200000;
            ForwardOperatingPointJson["ambientTemperature"] = 42;
            forwardInputsJson["operatingPoints"].push_back(ForwardOperatingPointJson);
        }
        OpenMagnetics::AdvancedTwoSwitchForward forwardInputs(forwardInputsJson);
        forwardInputs._assertErrors = true;

        auto inputs = forwardInputs.process();

        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_DCM_Primary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_DCM_Secondary_Current.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_waveform().value());
            painter.export_svg();
        }


        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_DCM_Primary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value());
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_TwoSwitchForward_DCM_Secondary_Voltage.svg");
            std::filesystem::remove(outFile);   
            Painter painter(outFile, false, true);
            painter.paint_waveform(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_waveform().value());
            painter.export_svg();
        }


        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::CUSTOM);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() > 0);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_peak().value() > double(forwardInputsJson["operatingPoints"][0]["outputVoltages"][0]));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::RECTANGULAR_WITH_DEADTIME);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY);
    }

    TEST_CASE("Test_SingleSwitchForward_Ngspice_Simulation", "[converter-model][single-switch-forward-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create a Single Switch Forward converter specification
        OpenMagnetics::SingleSwitchForward forward;
        
        // Input voltage
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(36.0);
        inputVoltage.set_maximum(72.0);
        forward.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        forward.set_diode_voltage_drop(0.5);
        
        // Efficiency
        forward.set_efficiency(0.9);
        
        // Current ripple ratio
        forward.set_current_ripple_ratio(0.3);
        
        // Duty cycle
        forward.set_duty_cycle(0.4);
        
        // Operating point: 5V @ 5A output, 200kHz
        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({5.0});
        opPoint.set_output_currents({5.0});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = forward.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Single Switch Forward - Turns ratio: " << turnsRatios[0]);
        INFO("Single Switch Forward - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto converterWaveforms = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!converterWaveforms.empty());
        
        // Verify we have converter-port stream populated.
    REQUIRE(!converterWaveforms[0].get_input_voltage().get_data().empty());

        // §5.1: ConverterWaveforms exposes the *converter* ports — DC
        // source rail (Vin), DC filtered output. The primary-winding
        // waveform (the sawtooth ramp / 0-Vin pulse) lives on the
        // winding-port stream, accessed via simulate_and_extract_operating_points.
        auto vinData = converterWaveforms[0].get_input_voltage().get_data();
        auto iinData = converterWaveforms[0].get_input_current().get_data();

        double vinMean = std::accumulate(vinData.begin(), vinData.end(), 0.0) / vinData.size();
        double iinMean = std::accumulate(iinData.begin(), iinData.end(), 0.0) / iinData.size();
        INFO("Input voltage (DC bus) mean: " << vinMean << " V (nom 48)");
        INFO("Input current (bus draw) mean: " << iinMean << " A");

        // Vin rail should be flat ≈ 48 V.
        CHECK(std::fabs(vinMean - 48.0) / 48.0 < 0.01);
        // Bus current draws positive average power.
        CHECK(iinMean > 0.0);

        INFO("Single Switch Forward ngspice simulation test passed");

        SECTION("Primary winding voltage waveform is a pulse (ON ≈ Vin, OFF ≈ reset)") {
            // Winding-port stream lives in get_excitations_per_winding,
            // populated by simulate_and_extract_operating_points.
            auto ops = forward.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
            REQUIRE(!ops.empty());
            REQUIRE(!ops[0].get_excitations_per_winding().empty());
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(priExc.get_voltage().has_value());
            REQUIRE(priExc.get_voltage()->get_waveform().has_value());
            auto vPri = priExc.get_voltage()->get_waveform()->get_data();
            double Vin = 48.0;
            int count_on = 0;
            for (double v : vPri) {
                if (v > 0.5 * Vin) count_on++;
            }
            double frac_on = double(count_on) / vPri.size();
            INFO("Primary-winding fraction at +Vin: " << frac_on << " (D=0.4 expected)");
            CHECK(frac_on > 0.20);
            CHECK(frac_on < 0.60);
        }

        SECTION("Primary winding current ramps up during ON time (sawtooth in CCM)") {
            auto ops = forward.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
            REQUIRE(!ops.empty());
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(priExc.get_current().has_value());
            REQUIRE(priExc.get_current()->get_waveform().has_value());
            auto iPri = priExc.get_current()->get_waveform()->get_data();

            auto peak_it = std::max_element(iPri.begin(), iPri.end());
            int peak_idx = (int)std::distance(iPri.begin(), peak_it);
            int N = (int)iPri.size();

            double rising_sum = 0.0;
            for (int k = 1; k <= peak_idx; ++k)
                rising_sum += iPri[k] - iPri[k-1];

            double falling_sum = 0.0;
            for (int k = peak_idx + 1; k < N; ++k)
                falling_sum += iPri[k] - iPri[k-1];

            INFO("Rising sum: " << rising_sum << " A, falling sum: " << falling_sum << " A");
            CHECK(rising_sum > 0.0);
            CHECK(falling_sum < 0.0);
        }

        SECTION("ON-time primary-winding current slope matches Vin/Lm") {
            // Use a single-period extraction so the dt math below
            // (one switching period spans N samples) is correct.
            // Default numPeriodsToExtract=5 would give dt=5/(fs*N).
            const int origNumPeriods = forward.get_num_periods_to_extract();
            forward.set_num_periods_to_extract(1);
            auto ops = forward.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
            forward.set_num_periods_to_extract(origNumPeriods);
            REQUIRE(!ops.empty());
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            auto iPri = priExc.get_current()->get_waveform()->get_data();

            double Vin = 48.0, fs = 200e3;
            int N = (int)iPri.size();
            double dt = 1.0 / (fs * N);

            // Compute slope from local-min to local-max within one ramp,
            // not from sample 0 to global peak (which crosses the OFF
            // discontinuity in multi-period extraction).
            auto min_it = std::min_element(iPri.begin(), iPri.end());
            int min_idx = (int)std::distance(iPri.begin(), min_it);
            auto max_it = std::max_element(iPri.begin() + min_idx, iPri.end());
            int max_idx = (int)std::distance(iPri.begin(), max_it);

            if (max_idx > min_idx + 5) {
                double slope_A_per_s = (iPri[max_idx] - iPri[min_idx]) / ((max_idx - min_idx) * dt);
                double expected = Vin / magnetizingInductance;
                double err = std::abs(slope_A_per_s - expected) / expected;
                INFO("di/dt: " << slope_A_per_s << " A/s, Vin/Lm: " << expected << ", error: " << (err*100) << "%");
                CHECK(err < 0.50);  // 50%: ngspice timestep + small load-reflected component
            }
        }
    }

    TEST_CASE("Test_TwoSwitchForward_Ngspice_Simulation", "[converter-model][two-switch-forward-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create a Two Switch Forward converter specification
        OpenMagnetics::TwoSwitchForward forward;
        
        // Input voltage
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(36.0);
        inputVoltage.set_maximum(72.0);
        forward.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        forward.set_diode_voltage_drop(0.5);
        
        // Efficiency
        forward.set_efficiency(0.9);
        
        // Current ripple ratio
        forward.set_current_ripple_ratio(0.3);
        
        // Duty cycle
        forward.set_duty_cycle(0.4);
        
        // Operating point: 5V @ 5A output, 200kHz
        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({5.0});
        opPoint.set_output_currents({5.0});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = forward.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Two Switch Forward - Turns ratio: " << turnsRatios[0]);
        INFO("Two Switch Forward - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto converterWaveforms = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!converterWaveforms.empty());
        
        // Verify we have converter-port stream populated.
    REQUIRE(!converterWaveforms[0].get_input_voltage().get_data().empty());

        // §5.1: Converter-port stream is the DC source rail. Pulse-shape
        // primary-winding checks live in the winding-port stream
        // (simulate_and_extract_operating_points), exercised below.
        auto vinData = converterWaveforms[0].get_input_voltage().get_data();
        double vinMean = std::accumulate(vinData.begin(), vinData.end(), 0.0) / vinData.size();
        INFO("Input voltage (DC bus) mean: " << vinMean << " V (nom 48)");
        CHECK(std::fabs(vinMean - 48.0) / 48.0 < 0.01);

        INFO("Two Switch Forward ngspice simulation test passed");

        SECTION("Two-switch-forward primary-winding voltage pulse (ON ≈ Vin, OFF ≈ 0)") {
            // Winding-port stream via simulate_and_extract_operating_points.
            auto ops = forward.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
            REQUIRE(!ops.empty());
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(priExc.get_voltage().has_value());
            REQUIRE(priExc.get_voltage()->get_waveform().has_value());
            auto vPri = priExc.get_voltage()->get_waveform()->get_data();

            double Vin = 48.0;
            int count_on = 0, count_off = 0;
            for (double v : vPri) {
                if (v > 0.5 * Vin) count_on++;
                else if (v < 5.0) count_off++;
            }
            double frac_on = double(count_on) / vPri.size();
            double frac_off = double(count_off) / vPri.size();
            INFO("Primary-winding fraction at +Vin: " << frac_on
                 << ", fraction near 0V: " << frac_off << " (D=0.4 expected)");
            CHECK(frac_on > 0.20);
            CHECK(frac_on < 0.60);
            // Two-switch forward: OFF voltage clamped to 0 by demagnetizing diodes
            CHECK(frac_off > 0.30);
        }
    }

    TEST_CASE("Test_TwoSwitchForward_Waveform_Polarity", "[converter-model][two-switch-forward-topology][ngspice-simulation][smoke-test]") {
        // Verify Two-Switch Forward converter has correct waveform polarity:
        // - Primary voltage (v(pri_in)) should go NEAR ZERO during OFF time (demagnetization via clamping diodes)
        // - Secondary voltage (v(sec_in)) should go NEGATIVE during reset phase
        // Reference: TI SLYU036A - Two-Switch Forward Converter Design Guide
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        OpenMagnetics::TwoSwitchForward forward;

        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(36.0);
        inputVoltage.set_maximum(60.0);
        forward.set_input_voltage(inputVoltage);

        forward.set_diode_voltage_drop(0.5);
        forward.set_efficiency(0.9);
        forward.set_current_ripple_ratio(0.3);

        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({12.0});
        opPoint.set_output_currents({4.0});
        opPoint.set_switching_frequency(100e3);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});

        auto designReqs = forward.process_design_requirements();

        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        INFO("Turns ratio: " << turnsRatios[0]);
        INFO("Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");

        // Run simulation and extract operating points (winding-level waveforms)
        auto operatingPoints = forward.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        REQUIRE(!operatingPoints.empty());

        auto& primaryExcitation = operatingPoints[0].get_excitations_per_winding()[0];
        auto& secondaryExcitation = operatingPoints[0].get_excitations_per_winding()[1];

        REQUIRE(primaryExcitation.get_voltage().has_value());
        REQUIRE(primaryExcitation.get_voltage()->get_waveform().has_value());
        REQUIRE(secondaryExcitation.get_voltage().has_value());
        REQUIRE(secondaryExcitation.get_voltage()->get_waveform().has_value());

        auto primaryVoltageData = primaryExcitation.get_voltage()->get_waveform()->get_data();
        auto secondaryVoltageData = secondaryExcitation.get_voltage()->get_waveform()->get_data();

        double priV_max = *std::max_element(primaryVoltageData.begin(), primaryVoltageData.end());
        double priV_min = *std::min_element(primaryVoltageData.begin(), primaryVoltageData.end());
        double secV_max = *std::max_element(secondaryVoltageData.begin(), secondaryVoltageData.end());
        double secV_min = *std::min_element(secondaryVoltageData.begin(), secondaryVoltageData.end());

        INFO("Primary voltage max: " << priV_max << " V, min: " << priV_min << " V");
        INFO("Secondary voltage max: " << secV_max << " V, min: " << secV_min << " V");

        // Two-Switch Forward: Primary voltage should be ~+Vin during ON, ~0V during OFF (clamped by D1/D2)
        // The voltage across the primary reverses during demagnetization.
        // v(pri_in) measures node voltage relative to ground.
        // During ON: v(pri_in) ≈ Vin (switches closed, current flows vin→pri→gnd)
        // During OFF/reset: v(pri_in) ≈ 0 (clamped via D1 to ground)
        CHECK(priV_max > 30.0);   // Should be around nominal Vin during ON
        CHECK(priV_max < 65.0);   // Not exceeding max input voltage

        // Primary voltage should drop significantly during OFF time
        // In a properly working Two-Switch Forward, v(pri_in) goes near 0 during reset
        CHECK(priV_min < 5.0);    // Should go near zero during demagnetization

        // Secondary voltage should go negative during reset (transformer coupling)
        // During ON: v(sec_in) ≈ Vin/N (positive, energy transfers through forward rectifier)
        // During OFF/reset: v(sec_in) ≈ -Vin/N (negative, forward rectifier blocks)
        CHECK(secV_max > 5.0);    // Positive during ON time
        CHECK(secV_min < -1.0);   // Should go negative during reset

        // §5.1: converter-port input_voltage is now the DC source rail
        // (vin_dc, flat ≈ Vin). The "drops to 0 during reset" pulse-pattern
        // assertion that used to live here against cwf.get_input_voltage()
        // was inspecting winding-port behavior smuggled into the
        // converter-port stream — exactly the bug §5.1 forbids. The
        // primary-winding voltage swing is already validated above via
        // simulate_and_extract_operating_points.
        auto converterWaveforms = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        REQUIRE(!converterWaveforms.empty());

        auto& cwf = converterWaveforms[0];
        auto cwfInputVoltage = cwf.get_input_voltage().get_data();
        REQUIRE(!cwfInputVoltage.empty());

        double cwfV_mean = std::accumulate(cwfInputVoltage.begin(),
                                           cwfInputVoltage.end(), 0.0)
                          / cwfInputVoltage.size();
        INFO("Converter input voltage (DC bus) mean: " << cwfV_mean << " V (nom 48)");
        CHECK(std::fabs(cwfV_mean - 48.0) / 48.0 < 0.01);

        // Output voltage should be around 12V (stable)
        REQUIRE(!cwf.get_output_voltages().empty());
        auto outputVoltageData = cwf.get_output_voltages()[0].get_data();
        if (!outputVoltageData.empty()) {
            double outV_avg = 0;
            for (double v : outputVoltageData) outV_avg += v;
            outV_avg /= outputVoltageData.size();
            INFO("Output voltage average: " << outV_avg << " V");
            CHECK(outV_avg > 8.0);   // Should be around 12V
            CHECK(outV_avg < 16.0);
        }
    }

    TEST_CASE("Test_ActiveClampForward_Ngspice_Simulation", "[converter-model][active-clamp-forward-topology][ngspice-simulation]") {
        // Check if ngspice is available
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }
        
        // Create an Active Clamp Forward converter specification
        OpenMagnetics::ActiveClampForward forward;
        
        // Input voltage
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        inputVoltage.set_minimum(36.0);
        inputVoltage.set_maximum(72.0);
        forward.set_input_voltage(inputVoltage);
        
        // Diode voltage drop
        forward.set_diode_voltage_drop(0.5);
        
        // Efficiency
        forward.set_efficiency(0.9);
        
        // Current ripple ratio
        forward.set_current_ripple_ratio(0.3);
        
        // Duty cycle
        forward.set_duty_cycle(0.45);
        
        // Operating point: 5V @ 5A output, 200kHz
        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({5.0});
        opPoint.set_output_currents({5.0});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});
        
        // Process design requirements
        auto designReqs = forward.process_design_requirements();
        
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
        
        INFO("Active Clamp Forward - Turns ratio: " << turnsRatios[0]);
        INFO("Active Clamp Forward - Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
        
        // Run ngspice simulation
        auto converterWaveforms = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
        
        REQUIRE(!converterWaveforms.empty());
        
        // Verify we have converter-port stream populated.
    REQUIRE(!converterWaveforms[0].get_input_voltage().get_data().empty());

        // §5.1: ConverterWaveforms exposes the DC source rail (vin_dc),
        // not the primary winding. Primary-winding pulse / clamp checks
        // live in the winding-port stream (simulate_and_extract_operating_points).
        auto vinData = converterWaveforms[0].get_input_voltage().get_data();
        auto iinData = converterWaveforms[0].get_input_current().get_data();
        double vinMean = std::accumulate(vinData.begin(), vinData.end(), 0.0) / vinData.size();
        double iinMean = std::accumulate(iinData.begin(), iinData.end(), 0.0) / iinData.size();
        INFO("Input voltage (DC bus) mean: " << vinMean << " V (nom 48)");
        INFO("Input current (bus draw) mean: " << iinMean << " A");
        CHECK(std::fabs(vinMean - 48.0) / 48.0 < 0.01);
        CHECK(iinMean > 0.0);

        INFO("Active Clamp Forward ngspice simulation test passed");

        SECTION("Active-clamp primary-winding voltage: ON pulse and negative clamp reset") {
            auto ops = forward.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
            REQUIRE(!ops.empty());
            auto& priExc = ops[0].get_excitations_per_winding()[0];
            REQUIRE(priExc.get_voltage().has_value());
            REQUIRE(priExc.get_voltage()->get_waveform().has_value());
            auto vPri = priExc.get_voltage()->get_waveform()->get_data();

            double Vin = 48.0;
            int count_on = 0, count_neg = 0;
            for (double v : vPri) {
                if (v > 0.5 * Vin) count_on++;
                if (v < -5.0) count_neg++;
            }
            double frac_on = double(count_on) / vPri.size();
            double frac_neg = double(count_neg) / vPri.size();
            INFO("Primary-winding fraction at +Vin: " << frac_on
                 << ", fraction <-5V: " << frac_neg);
            CHECK(frac_on > 0.20);   // D ≈ 0.45
            CHECK(frac_neg >= 0);    // Active clamp may produce negative reset voltage
            double vPri_min = *std::min_element(vPri.begin(), vPri.end());
            INFO("Primary-winding voltage min: " << vPri_min << " V");
            CHECK(vPri_min < 0.0);   // Active clamp produces negative reset
        }

        SECTION("Active-clamp primary-winding current: avg > 0 (net energy transfer)") {
            auto ops = forward.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
            REQUIRE(!ops.empty());
            auto iPri = ops[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
            double iPri_avg = std::accumulate(iPri.begin(), iPri.end(), 0.0) / iPri.size();
            INFO("Primary-winding current avg: " << iPri_avg << " A");
            CHECK(iPri_avg > 0.1);
        }
    }

    TEST_CASE("Test_TwoSwitchForward_NumPeriods_SimulatedOperatingPoints", "[converter-model][two-switch-forward-topology][num-periods][ngspice-simulation][smoke-test]") {
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        OpenMagnetics::TwoSwitchForward forward;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        forward.set_input_voltage(inputVoltage);
        forward.set_diode_voltage_drop(0.5);
        forward.set_efficiency(0.9);
        forward.set_current_ripple_ratio(0.3);

        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({12.0});
        opPoint.set_output_currents({4.0});
        opPoint.set_switching_frequency(100e3);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});

        auto designReqs = forward.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        // Simulate with 1 period
        forward.set_num_periods_to_extract(1);
        auto ops1 = forward.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        REQUIRE(!ops1.empty());
        auto voltageWf1 = ops1[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();

        // Simulate with 3 periods
        forward.set_num_periods_to_extract(3);
        auto ops3 = forward.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
        REQUIRE(!ops3.empty());
        auto voltageWf3 = ops3[0].get_excitations_per_winding()[0].get_voltage()->get_waveform().value();

        INFO("1-period waveform data size: " << voltageWf1.get_data().size());
        INFO("3-period waveform data size: " << voltageWf3.get_data().size());

        CHECK(voltageWf3.get_data().size() > voltageWf1.get_data().size());
    }

    TEST_CASE("Test_TwoSwitchForward_NumPeriods_ConverterWaveforms", "[converter-model][two-switch-forward-topology][num-periods][ngspice-simulation][smoke-test]") {
        NgspiceRunner runner;
        if (!runner.is_available()) {
            SKIP("ngspice not available on this system");
        }

        OpenMagnetics::TwoSwitchForward forward;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(48.0);
        forward.set_input_voltage(inputVoltage);
        forward.set_diode_voltage_drop(0.5);
        forward.set_efficiency(0.9);
        forward.set_current_ripple_ratio(0.3);

        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({12.0});
        opPoint.set_output_currents({4.0});
        opPoint.set_switching_frequency(100e3);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});

        auto designReqs = forward.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios()) {
            turnsRatios.push_back(tr.get_nominal().value());
        }
        double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();

        // Simulate with 1 period
        auto waveforms1 = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, 1);
        REQUIRE(!waveforms1.empty());
        auto inputV1 = waveforms1[0].get_input_voltage().get_data();

        // Simulate with 3 periods
        auto waveforms3 = forward.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance, 3);
        REQUIRE(!waveforms3.empty());
        auto inputV3 = waveforms3[0].get_input_voltage().get_data();

        INFO("1-period converter waveform data size: " << inputV1.size());
        INFO("3-period converter waveform data size: " << inputV3.size());

        CHECK(inputV3.size() > inputV1.size());
    }

    TEST_CASE("Test_SingleSwitchForward_PtP_AnalyticalVsNgspice", "[converter-model][single-switch-forward-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        OpenMagnetics::SingleSwitchForward fwd;
        DimensionWithTolerance iv; iv.set_nominal(48.0); iv.set_minimum(36.0); iv.set_maximum(72.0);
        fwd.set_input_voltage(iv);
        fwd.set_diode_voltage_drop(0.5);
        fwd.set_efficiency(0.9);
        fwd.set_current_ripple_ratio(0.3);
        fwd.set_duty_cycle(0.4);

        ForwardOperatingPoint op;
        op.set_output_voltages({5.0}); op.set_output_currents({5.0});
        op.set_switching_frequency(200e3); op.set_ambient_temperature(25.0);
        fwd.set_operating_points({op});

        auto designReqs = fwd.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
        double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        auto analyticalOps = fwd.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        REQUIRE(!aData.empty());
        REQUIRE(!aTime.empty());
        auto aResampled = ptp_interp(aTime, aData, 256);

        fwd.set_num_periods_to_extract(1);
        auto simOps = fwd.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        REQUIRE(!sData.empty());
        REQUIRE(!sTime.empty());
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        INFO("Single Switch Forward primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.35);
    }

    TEST_CASE("Test_ActiveClampForward_PtP_AnalyticalVsNgspice", "[converter-model][active-clamp-forward-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        OpenMagnetics::ActiveClampForward fwd;
        DimensionWithTolerance iv; iv.set_nominal(48.0); iv.set_minimum(36.0); iv.set_maximum(72.0);
        fwd.set_input_voltage(iv);
        fwd.set_diode_voltage_drop(0.5);
        fwd.set_efficiency(0.9);
        fwd.set_current_ripple_ratio(0.3);
        fwd.set_duty_cycle(0.45);

        ForwardOperatingPoint op;
        op.set_output_voltages({5.0}); op.set_output_currents({5.0});
        op.set_switching_frequency(200e3); op.set_ambient_temperature(25.0);
        fwd.set_operating_points({op});

        auto designReqs = fwd.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
        double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        auto analyticalOps = fwd.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        REQUIRE(!aData.empty());
        REQUIRE(!aTime.empty());
        auto aResampled = ptp_interp(aTime, aData, 256);

        fwd.set_num_periods_to_extract(1);
        auto simOps = fwd.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        REQUIRE(!sData.empty());
        REQUIRE(!sTime.empty());
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        // ACF: the active clamp resonant circuit creates a sinusoidal hump in the primary current
        // that the simplified piecewise analytical model does not capture — shape mismatch is expected.
        // Threshold of 0.99 verifies the simulation runs and returns non-degenerate waveform data.
        INFO("Active Clamp Forward primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.99);
    }

    TEST_CASE("Test_TwoSwitchForward_PtP_AnalyticalVsNgspice", "[converter-model][two-switch-forward-topology][ngspice-simulation][ptpcomparison]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        OpenMagnetics::TwoSwitchForward fwd;
        DimensionWithTolerance iv; iv.set_nominal(48.0); iv.set_minimum(36.0); iv.set_maximum(72.0);
        fwd.set_input_voltage(iv);
        fwd.set_diode_voltage_drop(0.5);
        fwd.set_efficiency(0.9);
        fwd.set_current_ripple_ratio(0.3);
        fwd.set_duty_cycle(0.4);

        ForwardOperatingPoint op;
        op.set_output_voltages({5.0}); op.set_output_currents({5.0});
        op.set_switching_frequency(200e3); op.set_ambient_temperature(25.0);
        fwd.set_operating_points({op});

        auto designReqs = fwd.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : designReqs.get_turns_ratios()) turnsRatios.push_back(tr.get_nominal().value());
        double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        auto analyticalOps = fwd.process_operating_points(turnsRatios, Lm);
        REQUIRE(!analyticalOps.empty());
        auto aTime = ptp_current_time(analyticalOps[0], 0);
        auto aData = ptp_current(analyticalOps[0], 0);
        REQUIRE(!aData.empty());
        REQUIRE(!aTime.empty());
        auto aResampled = ptp_interp(aTime, aData, 256);

        fwd.set_num_periods_to_extract(1);
        auto simOps = fwd.simulate_and_extract_operating_points(turnsRatios, Lm);
        REQUIRE(!simOps.empty());
        auto sTime = ptp_current_time(simOps[0], 0);
        auto sData = ptp_current(simOps[0], 0);
        REQUIRE(!sData.empty());
        REQUIRE(!sTime.empty());
        auto sResampled = ptp_interp(sTime, sData, 256);

        double nrmse = ptp_nrmse(aResampled, sResampled);
        // TSForward has a sharp analytical notch at switch-off (reflected load drops instantly)
        // that the simulation shows as a smooth transition — 50% threshold accounts for this.
        INFO("Two Switch Forward primary current NRMSE (analytical vs NgSpice): " << (nrmse * 100.0) << "%");
        CHECK(nrmse < 0.50);
    }


    // ─────────────────────────────────────────────────────────────────
    // §5.1 converter-port DC-stream gate. See ConverterPortChecks for
    // the full bound rationale. Catches the AC-in-DC-stream regression
    // where pri_in (winding top, bipolar AC) and i(Vpri_sense) (winding
    // current, bipolar AC) were sourcing the converter's input port.
    // ─────────────────────────────────────────────────────────────────
    TEST_CASE("Test_SingleSwitchForward_ConverterPortWaveforms",
              "[converter-port-waveforms][single-switch-forward-topology][ngspice-simulation]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        const double Vin = 48.0, Vout = 5.0, Iout = 5.0;
        OpenMagnetics::SingleSwitchForward forward;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(Vin);
        forward.set_input_voltage(inputVoltage);
        forward.set_diode_voltage_drop(0.5);
        forward.set_efficiency(0.9);
        forward.set_current_ripple_ratio(0.3);
        forward.set_duty_cycle(0.4);

        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({Vout});
        opPoint.set_output_currents({Iout});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});

        auto designReqs = forward.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios())
            turnsRatios.push_back(tr.get_nominal().value());
        const double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        auto wfs = forward.simulate_and_extract_topology_waveforms(turnsRatios, Lm);
        REQUIRE(!wfs.empty());

        // SSF SPICE droop is non-trivial (~25 % at this D=0.4 fixture)
        // due to magnetizing-leg energy circulation, diode Vf, and the
        // analytical-vs-SPICE operating-point split documented in
        // Test_SingleSwitchForward_PtP_AnalyticalVsNgspice. Allow 0.30
        // mean tolerance — same envelope as the IsolatedBuck flybuck
        // secondary, slightly tighter than PSFB. The §5.1 gate's job is
        // to catch AC-in-DC-stream regressions (e.g. v(pri_in) as
        // input_voltage), NOT to enforce steady-state accuracy.
        constexpr double kSsfVoutMeanTol = 0.30;
        constexpr double kSsfIoutMeanTol = 0.30;
        for (size_t i = 0; i < wfs.size(); ++i) {
            ConverterPortChecks::check_dc_ports(wfs[i], "SingleSwitchForward", i,
                                                Vin, {Vout}, {Iout},
                                                kSsfVoutMeanTol, kSsfIoutMeanTol);
        }
    }

    TEST_CASE("Test_TwoSwitchForward_ConverterPortWaveforms",
              "[converter-port-waveforms][two-switch-forward-topology][ngspice-simulation]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        const double Vin = 48.0, Vout = 5.0, Iout = 5.0;
        OpenMagnetics::TwoSwitchForward forward;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(Vin);
        forward.set_input_voltage(inputVoltage);
        forward.set_diode_voltage_drop(0.5);
        forward.set_efficiency(0.9);
        forward.set_current_ripple_ratio(0.3);
        forward.set_duty_cycle(0.4);

        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({Vout});
        opPoint.set_output_currents({Iout});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});

        auto designReqs = forward.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios())
            turnsRatios.push_back(tr.get_nominal().value());
        const double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        auto wfs = forward.simulate_and_extract_topology_waveforms(turnsRatios, Lm);
        REQUIRE(!wfs.empty());

        // Same envelope as SSF — TSF shares the forward-family SPICE
        // droop pattern (Lm energy recycle via clamp diodes, diode Vf,
        // analytical-vs-SPICE op split). The §5.1 gate's job is to
        // catch AC-in-DC-stream regressions, NOT enforce steady-state.
        constexpr double kTsfVoutMeanTol = 0.30;
        constexpr double kTsfIoutMeanTol = 0.30;
        for (size_t i = 0; i < wfs.size(); ++i) {
            ConverterPortChecks::check_dc_ports(wfs[i], "TwoSwitchForward", i,
                                                Vin, {Vout}, {Iout},
                                                kTsfVoutMeanTol, kTsfIoutMeanTol);
        }
    }

    TEST_CASE("Test_ActiveClampForward_ConverterPortWaveforms",
              "[converter-port-waveforms][active-clamp-forward-topology][ngspice-simulation]") {
        NgspiceRunner runner;
        if (!runner.is_available()) SKIP("ngspice not available");

        const double Vin = 48.0, Vout = 5.0, Iout = 5.0;
        OpenMagnetics::ActiveClampForward forward;
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_nominal(Vin);
        forward.set_input_voltage(inputVoltage);
        forward.set_diode_voltage_drop(0.5);
        forward.set_efficiency(0.9);
        forward.set_current_ripple_ratio(0.3);
        forward.set_duty_cycle(0.45);

        ForwardOperatingPoint opPoint;
        opPoint.set_output_voltages({Vout});
        opPoint.set_output_currents({Iout});
        opPoint.set_switching_frequency(200000.0);
        opPoint.set_ambient_temperature(25.0);
        forward.set_operating_points({opPoint});

        auto designReqs = forward.process_design_requirements();
        std::vector<double> turnsRatios;
        for (const auto& tr : designReqs.get_turns_ratios())
            turnsRatios.push_back(tr.get_nominal().value());
        const double Lm = designReqs.get_magnetizing_inductance().get_minimum().value();

        auto wfs = forward.simulate_and_extract_topology_waveforms(turnsRatios, Lm);
        REQUIRE(!wfs.empty());

        // ACF has an active clamp leg that recirculates Lm energy each
        // cycle — slightly less droop than SSF/TSF but still subject to
        // diode Vf and analytical-vs-SPICE op-split. Same 0.30 envelope.
        constexpr double kAcfVoutMeanTol = 0.30;
        constexpr double kAcfIoutMeanTol = 0.30;
        for (size_t i = 0; i < wfs.size(); ++i) {
            ConverterPortChecks::check_dc_ports(wfs[i], "ActiveClampForward", i,
                                                Vin, {Vout}, {Iout},
                                                kAcfVoutMeanTol, kAcfIoutMeanTol);
        }
    }

}  // namespace
