// Tests for OpenMagnetics::WaveformProcessor — the pure converter-waveform DSP
// extracted from Inputs/Topology. These exercise the ported functions
// (create_waveform, calculate_sampled_waveform, calculate_harmonics_data,
// calculate_processed_data, calculate_basic_processed_data, complete_excitation)
// and the ported helpers (try_guess_duty_cycle, try_guess_waveform_label,
// compress_waveform, is_waveform_sampled). They were moved here from
// TestInputs.cpp so the tests live with the logic, and they call
// WaveformProcessor:: directly rather than going through Inputs.

#include "processors/WaveformProcessor.h"
#include "processors/Inputs.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <magic_enum.hpp>
#include <chrono>
#include <cmath>
#include <numbers>
#include <random>
#include <vector>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

TEST_CASE("Test_Get_Duty_Cycle_Triangular", "[processor][waveform-processor][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["waveform"]["data"] = {-5, 5, -5};
    windingExcitation["current"]["waveform"]["time"] = {0, 0.0000025, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    auto dutyCycle = WaveformProcessor::try_guess_duty_cycle(excitation.get_current().value().get_waveform().value());
    REQUIRE(dutyCycle == 0.25);
}

TEST_CASE("Test_Get_Duty_Cycle_Rectangular", "[processor][waveform-processor][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["waveform"]["data"] = {7.5, 7.5, -2.5, -2.5, 7.5};
    windingExcitation["current"]["waveform"]["time"] = {0, 0.0000075, 0.0000075, 0.00001, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    auto dutyCycle = WaveformProcessor::try_guess_duty_cycle(excitation.get_current().value().get_waveform().value());
    double max_error = 0.02;

    REQUIRE_THAT(0.75, Catch::Matchers::WithinAbs(dutyCycle, max_error * 0.75));
}

TEST_CASE("Test_Get_Duty_Cycle_Custom", "[processor][waveform-processor][smoke-test]") {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"] = json();
    operatingPoint["conditions"]["ambientTemperature"] = 42;

    json windingExcitation = json();
    windingExcitation["frequency"] = 100000;
    windingExcitation["current"]["waveform"]["data"] = {0, 3, 8, 3, 0};
    windingExcitation["current"]["waveform"]["time"] = {0, 0.0000025, 0.0000042, 0.0000075, 0.00001};
    operatingPoint["excitationsPerWinding"] = json::array();
    operatingPoint["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPoint);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    OpenMagnetics::Inputs inputs(inputsJson);

    auto excitation = inputs.get_operating_points()[0].get_excitations_per_winding()[0];
    auto dutyCycle = WaveformProcessor::try_guess_duty_cycle(excitation.get_current().value().get_waveform().value());

    REQUIRE(dutyCycle == 0.42);
}

TEST_CASE("Test_Try_Guess_Duty_Cycle", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double frequency = 100000;

    Waveform voltageWaveform;
    voltageWaveform.set_time(std::vector<double>{0, dutyCycle / frequency, (dutyCycle - 0.01) / frequency, 1.0 / frequency, 1.0 / frequency});
    voltageWaveform.set_data(std::vector<double>{peakToPeak * (1 - dutyCycle), peakToPeak * (1 - dutyCycle), -peakToPeak * dutyCycle, -peakToPeak * dutyCycle, peakToPeak * (1 - dutyCycle)});

    double expectedValue = 0.25;

    auto dutyCycleGuessed = WaveformProcessor::try_guess_duty_cycle(voltageWaveform);
    REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(dutyCycleGuessed, max_error * expectedValue));
}

TEST_CASE("Test_Try_Guess_Sinusoidal", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.5;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::SINUSOIDAL;

    ProcessedWaveform processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = WaveformProcessor::create_waveform(processed, frequency);
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);

    auto calculatedProcessed = WaveformProcessor::calculate_basic_processed_data(waveform);

    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(processed.get_offset(), Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Triangular", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::TRIANGULAR;

    ProcessedWaveform processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = WaveformProcessor::create_waveform(processed, frequency);
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);

    auto calculatedProcessed = WaveformProcessor::calculate_basic_processed_data(waveform);

    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(processed.get_offset(), Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Unipolar_Triangular", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::UNIPOLAR_TRIANGULAR;

    ProcessedWaveform processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = WaveformProcessor::create_waveform(processed, frequency);
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);

    auto calculatedProcessed = WaveformProcessor::calculate_processed_data(waveform, frequency);

    double expectedAverage = (peakToPeak * dutyCycle / 2);
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(0, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE_THAT(expectedAverage, Catch::Matchers::WithinAbs(calculatedProcessed.get_average().value(), max_error * expectedAverage));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Unipolar_Rectangular", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::UNIPOLAR_RECTANGULAR;

    ProcessedWaveform processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = WaveformProcessor::create_waveform(processed, frequency);
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);

    auto calculatedProcessed = WaveformProcessor::calculate_processed_data(waveform, frequency);

    double expectedAverage = peakToPeak * dutyCycle;
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(0, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE_THAT(expectedAverage, Catch::Matchers::WithinAbs(calculatedProcessed.get_average().value(), max_error * expectedAverage));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Rectangular", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.25;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::RECTANGULAR;

    ProcessedWaveform processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = WaveformProcessor::create_waveform(processed, frequency);
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);

    auto calculatedProcessed = WaveformProcessor::calculate_basic_processed_data(waveform);

    double expectedOffset = 0;
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(expectedOffset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Bipolar_Rectangular", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.2;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.5;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::BIPOLAR_RECTANGULAR;

    ProcessedWaveform processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = WaveformProcessor::create_waveform(processed, frequency);
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);

    auto calculatedProcessed = WaveformProcessor::calculate_basic_processed_data(waveform);

    double expectedOffset = 0;
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(expectedOffset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Bipolar_Triangular", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.2;
    double peakToPeak = 53.3333;
    double dutyCycle = 0.5;
    double offset = 0;
    double frequency = 100000;
    auto label = WaveformLabel::BIPOLAR_TRIANGULAR;

    ProcessedWaveform processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = WaveformProcessor::create_waveform(processed, frequency);
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);

    auto calculatedProcessed = WaveformProcessor::calculate_basic_processed_data(waveform);

    double expectedOffset = 0;
    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE_THAT(expectedOffset, Catch::Matchers::WithinAbs(calculatedProcessed.get_offset(), max_error));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Flyback_Primary", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 50;
    double dutyCycle = 0.25;
    double offset = 10;
    double frequency = 100000;
    auto label = WaveformLabel::FLYBACK_PRIMARY;

    ProcessedWaveform processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = WaveformProcessor::create_waveform(processed, frequency);
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);

    auto calculatedProcessed = WaveformProcessor::calculate_basic_processed_data(waveform);

    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Try_Guess_Flyback_Secondary", "[processor][waveform-processor][smoke-test]") {
    double max_error = 0.1;
    double peakToPeak = 50;
    double dutyCycle = 0.25;
    double offset = 10;
    double frequency = 100000;
    auto label = WaveformLabel::FLYBACK_SECONDARY;

    ProcessedWaveform processed;
    processed.set_peak_to_peak(peakToPeak);
    processed.set_duty_cycle(dutyCycle);
    processed.set_offset(offset);
    processed.set_label(label);

    auto waveform = WaveformProcessor::create_waveform(processed, frequency);
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);

    auto calculatedProcessed = WaveformProcessor::calculate_basic_processed_data(waveform);

    REQUIRE(magic_enum::enum_name(label) == magic_enum::enum_name(guessedLabel));
    REQUIRE_THAT(processed.get_peak_to_peak().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_peak_to_peak().value(), max_error * processed.get_peak_to_peak().value()));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(calculatedProcessed.get_duty_cycle().value(), max_error * processed.get_duty_cycle().value()));
    REQUIRE(magic_enum::enum_name(calculatedProcessed.get_label()) == magic_enum::enum_name(processed.get_label()));
}

TEST_CASE("Test_Sample_And_Compress_Sinusoidal", "[processor][waveform-processor][smoke-test]") {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double peakToPeak = 50;
    double offset = 0;
    size_t numberPointsSampledWaveforms = 69;
    for (size_t i = 0; i < numberPointsSampledWaveforms; ++i) {
        double angle = i * 2 * M_PI / numberPointsSampledWaveforms;
        time.push_back(angle);
        data.push_back((sin(angle) * peakToPeak / 2) + offset);
    }
    waveform.set_data(data);
    waveform.set_time(time);

    REQUIRE(!WaveformProcessor::is_waveform_sampled(waveform));

    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform);
    REQUIRE(WaveformProcessor::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = WaveformProcessor::compress_waveform(sampledWaveform);
    REQUIRE(115UL == compressedWaveform.get_data().size());
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::SINUSOIDAL) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Rectangular", "[processor][waveform-processor][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({-2.5, 7.5, 7.5, -2.5, -2.5}));
    waveform.set_time(std::vector<double>({0, 0, 0.0000025, 0.0000025, 0.00001}));

    REQUIRE(!WaveformProcessor::is_waveform_sampled(waveform));

    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform);
    REQUIRE(WaveformProcessor::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = WaveformProcessor::compress_waveform(sampledWaveform);
    REQUIRE(5UL == compressedWaveform.get_data().size());
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::RECTANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Triangular", "[processor][waveform-processor][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({-5, 5, -5}));
    waveform.set_time(std::vector<double>({0, 0.0000025, 0.00001}));

    REQUIRE(!WaveformProcessor::is_waveform_sampled(waveform));

    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform);
    REQUIRE(WaveformProcessor::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = WaveformProcessor::compress_waveform(sampledWaveform);
    REQUIRE(3UL == compressedWaveform.get_data().size());
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::TRIANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Unipolar_Triangular", "[processor][waveform-processor][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({10, 60, 10, 10}));
    waveform.set_time(std::vector<double>({0, 0.0000025, 0.0000025, 0.00001}));

    REQUIRE(!WaveformProcessor::is_waveform_sampled(waveform));

    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform);
    REQUIRE(WaveformProcessor::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = WaveformProcessor::compress_waveform(sampledWaveform);
    REQUIRE(4UL == compressedWaveform.get_data().size());
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::UNIPOLAR_TRIANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Unipolar_Rectangular", "[processor][waveform-processor][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({10, 60, 60, 10, 10}));
    waveform.set_time(std::vector<double>({0, 0, 0.0000025, 0.0000025, 0.00001}));

    REQUIRE(!WaveformProcessor::is_waveform_sampled(waveform));

    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform);
    REQUIRE(WaveformProcessor::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = WaveformProcessor::compress_waveform(sampledWaveform);
    REQUIRE(5UL == compressedWaveform.get_data().size());
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::UNIPOLAR_RECTANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Bipolar_Rectangular", "[processor][waveform-processor][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({0, 0, 25, 25, 0, 0, -25, -25, 0, 0}));
    waveform.set_time(std::vector<double>({0, 1.25e-6, 1.25e-6, 3.75e-6, 3.75e-6, 6.25e-6, 6.25e-6, 8.75e-6, 8.75e-6, 10e-6}));

    REQUIRE(!WaveformProcessor::is_waveform_sampled(waveform));

    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform);
    REQUIRE(WaveformProcessor::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = WaveformProcessor::compress_waveform(sampledWaveform);
    REQUIRE(10UL == compressedWaveform.get_data().size());
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::BIPOLAR_RECTANGULAR) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Flyback_Primary", "[processor][waveform-processor][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({0, 30, 80, 0, 0}));
    waveform.set_time(std::vector<double>({0, 0, 1.4e-6, 1.4e-6, 0.00001}));

    REQUIRE(!WaveformProcessor::is_waveform_sampled(waveform));

    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform);
    REQUIRE(WaveformProcessor::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = WaveformProcessor::compress_waveform(sampledWaveform);
    REQUIRE(5UL == compressedWaveform.get_data().size());
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::FLYBACK_PRIMARY) == magic_enum::enum_name(guessedLabel));
}

TEST_CASE("Test_Sample_And_Compress_Flyback_Secondary", "[processor][waveform-processor][smoke-test]") {
    Waveform waveform;
    waveform.set_data(std::vector<double>({0, 0, 80, 30, 0}));
    waveform.set_time(std::vector<double>({0, 1.4e-6, 1.4e-6, 0.00001, 0.00001}));

    REQUIRE(!WaveformProcessor::is_waveform_sampled(waveform));

    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform);
    REQUIRE(WaveformProcessor::is_waveform_sampled(sampledWaveform));
    REQUIRE(128UL == sampledWaveform.get_data().size());

    auto compressedWaveform = WaveformProcessor::compress_waveform(sampledWaveform);
    REQUIRE(5UL == compressedWaveform.get_data().size());
    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(compressedWaveform);
    REQUIRE(magic_enum::enum_name(WaveformLabel::FLYBACK_SECONDARY) == magic_enum::enum_name(guessedLabel));
}

// Direct smoke test for complete_excitation (ported from Topology). It builds
// primary current/voltage waveforms with create_waveform and confirms
// complete_excitation populates a full OperatingPointExcitation (sampled
// waveform + harmonics + processed) for both signals.
TEST_CASE("Test_Complete_Excitation_Smoke", "[processor][waveform-processor][smoke-test]") {
    double frequency = 100000;

    auto currentWaveform = WaveformProcessor::create_waveform(WaveformLabel::TRIANGULAR, 10, frequency, 0.5);
    auto voltageWaveform = WaveformProcessor::create_waveform(WaveformLabel::RECTANGULAR, 100, frequency, 0.5);

    auto excitation = WaveformProcessor::complete_excitation(currentWaveform, voltageWaveform, frequency, "Primary");

    REQUIRE(excitation.get_frequency() == frequency);
    REQUIRE(excitation.get_name().value() == "Primary");

    REQUIRE(excitation.get_current());
    REQUIRE(excitation.get_current()->get_waveform());
    REQUIRE(excitation.get_current()->get_harmonics());
    REQUIRE(excitation.get_current()->get_processed());
    REQUIRE(WaveformProcessor::is_waveform_sampled(excitation.get_current()->get_waveform().value()));

    REQUIRE(excitation.get_voltage());
    REQUIRE(excitation.get_voltage()->get_waveform());
    REQUIRE(excitation.get_voltage()->get_harmonics());
    REQUIRE(excitation.get_voltage()->get_processed());
    REQUIRE(WaveformProcessor::is_waveform_sampled(excitation.get_voltage()->get_waveform().value()));
}

// --- ABT #602 -------------------------------------------------------------
// An imported (SPICE-simulated) waveform used to be classified SINUSOIDAL no
// matter its shape, because the sine-fit error was divided by N twice: `area`
// accumulated a SUM while `error` was averaged. At N = 512 the 5% acceptance
// threshold became an effective 2560% relative error, so CUSTOM was
// unreachable. try_guess_duty_cycle then returned the canonical 0.5 for the
// SINUSOIDAL label without ever looking at the samples — a user importing a
// FlyBuck switch-node voltage with a ~21% duty was shown 50%.
//
// The rectangle below mimics the reported waveform: +15.086 V for 20.7% of the
// period and -3.92 V for the rest, with finite (2-sample) edges so it looks
// like real simulator output rather than an ideal step.
namespace {
// `ripple` reproduces what makes real simulator output different from an ideal
// rectangle: the plateaus are not flat, so compress_waveform cannot collapse the
// signal to the 5 vertices the RECTANGULAR branch needs (the reported file
// compressed to 61 points). With ripple = 0 this builds the ideal rectangle,
// which SHOULD still be recognised as RECTANGULAR.
Waveform build_imported_rectangle(size_t numberPoints, double dutyCycle, double high, double low,
                                  double period, double ripple) {
    std::vector<double> data;
    std::vector<double> time;
    size_t pointsHigh = static_cast<size_t>(numberPoints * dutyCycle);
    for (size_t i = 0; i < numberPoints; ++i) {
        double value = (i < pointsHigh) ? high : low;
        // Two-sample linear edges, as a finite-slew simulator export would have.
        if (i == pointsHigh || i == pointsHigh + 1) {
            value = high + (low - high) * (i - pointsHigh + 1) / 3.0;
        }
        // Deterministic plateau ripple — no rand(), so the test cannot flake.
        value += ripple * sin(2 * 3.14159265358979323846 * i * 7 / numberPoints);
        data.push_back(value);
        time.push_back(period * i / numberPoints);
    }
    Waveform waveform;
    waveform.set_data(data);
    waveform.set_time(time);
    return waveform;
}
}  // namespace

TEST_CASE("Test_Ideal_Rectangle_Still_Recognised_As_Rectangular", "[processor][waveform-processor][smoke-test]") {
    // The clean case must keep using the analytical vertex path, duty and all.
    auto waveform = build_imported_rectangle(512, 0.207, 15.086, -3.92, 1.0 / 300000, 0.0);

    auto processed = WaveformProcessor::calculate_basic_processed_data(waveform);
    REQUIRE(magic_enum::enum_name(processed.get_label()) == magic_enum::enum_name(WaveformLabel::RECTANGULAR));
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(0.207, 0.01));
}

TEST_CASE("Test_Imported_Rectangle_Is_Custom_Not_Sinusoidal", "[processor][waveform-processor][smoke-test]") {
    // 0.25 V of plateau ripple on a 19 V swing (~1.3%), matching the reported file.
    auto waveform = build_imported_rectangle(512, 0.207, 15.086, -3.92, 1.0 / 300000, 0.25);

    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);
    REQUIRE(magic_enum::enum_name(guessedLabel) == magic_enum::enum_name(WaveformLabel::CUSTOM));

    auto processed = WaveformProcessor::calculate_basic_processed_data(waveform);
    REQUIRE(magic_enum::enum_name(processed.get_label()) == magic_enum::enum_name(WaveformLabel::CUSTOM));
    // Measured, not assumed: the duty must track the real high time, not be 0.5.
    REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(0.207, 0.01));
}

TEST_CASE("Test_Imported_Sine_Still_Classifies_Sinusoidal", "[processor][waveform-processor][smoke-test]") {
    // The other side of the fix: tightening the error metric must not cost us
    // the genuine sine, or every analytical path would start reporting CUSTOM.
    std::vector<double> data;
    std::vector<double> time;
    const size_t numberPoints = 512;
    const double period = 1.0 / 300000;
    for (size_t i = 0; i < numberPoints; ++i) {
        data.push_back(3.0 * sin(2 * 3.14159265358979323846 * i / numberPoints) + 1.0);
        time.push_back(period * i / numberPoints);
    }
    Waveform waveform;
    waveform.set_data(data);
    waveform.set_time(time);

    auto guessedLabel = WaveformProcessor::try_guess_waveform_label(waveform);
    REQUIRE(magic_enum::enum_name(guessedLabel) == magic_enum::enum_name(WaveformLabel::SINUSOIDAL));
}

TEST_CASE("Test_Phase_Shifted_Non_Sine_Stays_Custom", "[processor][waveform-processor][smoke-test]") {
    // Guarding the other direction: making the comparison phase-aware must not
    // let a triangle or a square through as a sine. Both sit around 0.27-0.36
    // mean relative error against their own best-phase sine, far above the 0.05
    // the classifier accepts, at every phase.
    const double pi = 3.14159265358979323846;
    const size_t numberPoints = 512;
    const double period = 1.0 / 100000;

    for (double shift : {0.0, 0.125, 0.25, 0.5, 0.8}) {
        std::vector<double> triangle;
        std::vector<double> square;
        std::vector<double> time;
        for (size_t i = 0; i < numberPoints; ++i) {
            double positionInPeriod = std::fmod(double(i) / numberPoints + shift, 1.0);
            triangle.push_back(positionInPeriod < 0.5 ? 4 * positionInPeriod - 1 : 3 - 4 * positionInPeriod);
            square.push_back(positionInPeriod < 0.5 ? 1.0 : -1.0);
            time.push_back(period * i / numberPoints);
        }
        // Sampled at 512 points these are imported data, not the 3- and 5-point
        // analytical shapes the vertex tests recognise, so CUSTOM is the honest
        // answer — the point here is that it is never SINUSOIDAL.
        for (const auto& data : {triangle, square}) {
            Waveform waveform;
            waveform.set_data(data);
            waveform.set_time(time);

            INFO("shift " << shift << " of a period");
            REQUIRE(magic_enum::enum_name(WaveformProcessor::try_guess_waveform_label(waveform)) !=
                    magic_enum::enum_name(WaveformLabel::SINUSOIDAL));
        }
    }
}

TEST_CASE("Test_Triangular_Recognised_At_Every_Duty", "[processor][waveform-processor][smoke-test]") {
    // A triangle's apex only lands ON a sample when the duty divides the sample
    // count. At any other duty compression cannot pick one vertex, emits both
    // neighbours, and the 3-point test missed a mathematically perfect triangle:
    // an ideal 1024-point triangle was TRIANGULAR at duty 0.5 and CUSTOM at 0.1,
    // 0.2, 0.3, 0.7 and 0.8.
    //
    // The duty came out wrong with it. Falling through to CUSTOM sent the
    // measurement to the mid-level crossing, which a triangle crosses halfway up
    // and halfway down whatever its ramps — so every duty reported 0.5.
    const size_t numberPoints = 1024;
    const double period = 1.0 / 100000;

    for (double duty : {0.1, 0.2, 0.3, 0.5, 0.7, 0.8}) {
        std::vector<double> data;
        std::vector<double> time;
        for (size_t i = 0; i < numberPoints; ++i) {
            double position = double(i) / numberPoints;
            data.push_back(position < duty ? (position / duty) * 2 - 1
                                           : 1 - ((position - duty) / (1 - duty)) * 2);
            time.push_back(period * i / numberPoints);
        }
        Waveform waveform;
        waveform.set_data(data);
        waveform.set_time(time);

        INFO("duty " << duty);
        auto processed = WaveformProcessor::calculate_basic_processed_data(waveform);
        REQUIRE(magic_enum::enum_name(processed.get_label()) ==
                magic_enum::enum_name(WaveformLabel::TRIANGULAR));
        REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(duty, 0.01));
    }
}

TEST_CASE("Test_Rectangular_Duty_Unaffected_By_Triangle_Path", "[processor][waveform-processor][smoke-test]") {
    // The control for the test above: rectangles already reported their duty
    // correctly, and the triangle work must not disturb either the label or the
    // measurement. A rectangle's duty IS its mid-level crossing.
    const size_t numberPoints = 1024;
    const double period = 1.0 / 100000;

    for (double duty : {0.2, 0.5, 0.8}) {
        std::vector<double> data;
        std::vector<double> time;
        for (size_t i = 0; i < numberPoints; ++i) {
            data.push_back(double(i) / numberPoints < duty ? 1.0 : -1.0);
            time.push_back(period * i / numberPoints);
        }
        Waveform waveform;
        waveform.set_data(data);
        waveform.set_time(time);

        INFO("duty " << duty);
        auto processed = WaveformProcessor::calculate_basic_processed_data(waveform);
        REQUIRE(magic_enum::enum_name(processed.get_label()) ==
                magic_enum::enum_name(WaveformLabel::RECTANGULAR));
        REQUIRE_THAT(processed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(duty, 0.01));
    }
}

TEST_CASE("Test_Noisy_Triangle_Still_Triangular_And_Noisy_Sine_Still_Sinusoidal",
          "[processor][waveform-processor][smoke-test]") {
    // Measured data never survives the vertex tests: noise above about 1e-6
    // relative leaves a real triangle compressing to 27-32 points, never 3, so
    // TRIANGULAR was unreachable for anything measured. A quarter of the
    // CoreDataX exchange is triangular and was reaching CUSTOM this way.
    //
    // The tolerance path recognises it, and must not do so by stealing sines.
    const double pi = 3.14159265358979323846;
    const size_t numberPoints = 1024;
    const double period = 1.0 / 100000;
    std::mt19937 generator(20260905);

    for (double noise : {1e-5, 1e-4, 1e-3}) {
        std::normal_distribution<double> jitter(0.0, noise);
        std::vector<double> triangle;
        std::vector<double> sine;
        std::vector<double> time;
        for (size_t i = 0; i < numberPoints; ++i) {
            double position = double(i) / numberPoints;
            triangle.push_back((position < 0.3 ? (position / 0.3) * 2 - 1
                                               : 1 - ((position - 0.3) / 0.7) * 2) + jitter(generator));
            sine.push_back(sin(2 * pi * position) + jitter(generator));
            time.push_back(period * i / numberPoints);
        }

        INFO("relative noise " << noise);
        Waveform noisyTriangle;
        noisyTriangle.set_data(triangle);
        noisyTriangle.set_time(time);
        auto triangleProcessed = WaveformProcessor::calculate_basic_processed_data(noisyTriangle);
        REQUIRE(magic_enum::enum_name(triangleProcessed.get_label()) ==
                magic_enum::enum_name(WaveformLabel::TRIANGULAR));
        REQUIRE_THAT(triangleProcessed.get_duty_cycle().value(), Catch::Matchers::WithinAbs(0.3, 0.02));

        Waveform noisySine;
        noisySine.set_data(sine);
        noisySine.set_time(time);
        REQUIRE(magic_enum::enum_name(WaveformProcessor::try_guess_waveform_label(noisySine)) ==
                magic_enum::enum_name(WaveformLabel::SINUSOIDAL));
    }
}

TEST_CASE("Test_Sampled_Sine_Classifies_Sinusoidal_At_Any_Phase", "[processor][waveform-processor][smoke-test]") {
    // The reference sine was built at phase zero, so a perfect sine that did not
    // start at a rising zero crossing was compared against a rotated reference and
    // came back CUSTOM — a cosine was never recognised as a sine at all. Every
    // analytical waveform MKF generates itself starts at a rising zero crossing,
    // so nothing here presented a shifted one.
    //
    // Goes through calculate_basic_processed_data on purpose. An earlier attempt
    // at this was verified only through try_guess_waveform_label on the raw
    // samples, passed, and still broke every real import — because
    // calculate_basic_processed_data used to classify from the COMPRESSED
    // waveform, which no phase or distortion measure can be taken on.
    const double pi = 3.14159265358979323846;
    const size_t numberPoints = 1024;
    const double period = 1.0 / 100000;

    for (double phaseDegrees : {0.0, 45.0, 90.0, 180.0, 270.0, 290.0, 327.0}) {
        std::vector<double> data;
        std::vector<double> time;
        for (size_t i = 0; i < numberPoints; ++i) {
            data.push_back(0.1 * sin(2 * pi * i / numberPoints + phaseDegrees * pi / 180) + 0.02);
            time.push_back(period * i / numberPoints);
        }
        Waveform waveform;
        waveform.set_data(data);
        waveform.set_time(time);

        INFO("phase " << phaseDegrees << " degrees");
        REQUIRE(magic_enum::enum_name(WaveformProcessor::calculate_basic_processed_data(waveform).get_label()) ==
                magic_enum::enum_name(WaveformLabel::SINUSOIDAL));
    }
}

TEST_CASE("Test_Phase_Shifted_Non_Sine_Never_Sinusoidal", "[processor][waveform-processor][smoke-test]") {
    // Guarding the other direction: a phase-aware comparison must not start
    // admitting triangles and squares. Which named label they land on is the
    // vertex tests' business; the invariant here is that neither is ever called a
    // sine, at any shift.
    const size_t numberPoints = 1024;
    const double period = 1.0 / 100000;

    for (double shift : {0.0, 0.125, 0.25, 0.375, 0.5}) {
        std::vector<double> triangle;
        std::vector<double> square;
        std::vector<double> time;
        for (size_t i = 0; i < numberPoints; ++i) {
            double positionInPeriod = std::fmod(double(i) / numberPoints + shift, 1.0);
            triangle.push_back(positionInPeriod < 0.5 ? 4 * positionInPeriod - 1 : 3 - 4 * positionInPeriod);
            square.push_back(positionInPeriod < 0.5 ? 1.0 : -1.0);
            time.push_back(period * i / numberPoints);
        }
        for (const auto& data : {triangle, square}) {
            Waveform waveform;
            waveform.set_data(data);
            waveform.set_time(time);

            INFO("shift " << shift << " of a period");
            REQUIRE(magic_enum::enum_name(WaveformProcessor::calculate_basic_processed_data(waveform).get_label()) !=
                    magic_enum::enum_name(WaveformLabel::SINUSOIDAL));
        }
    }
}

// An LTspice export with a fine time step carries hundreds of thousands of points in one switching
// period. Sampling restarted its segment search from the first point for every sample, so the cost
// was O(N*M): 400k points spent over two minutes in the web engine, whose watchdog then killed the
// import. With the search resuming where the previous sample stopped this is milliseconds; the bound
// is generous so machine load cannot flake it, and still two orders of magnitude under the old cost.
TEST_CASE("Test_Sampled_Waveform_Dense_Import_Is_Linear", "[processor][waveform-processor][smoke-test]") {
    const size_t numberPoints = 400000;
    const double period = 1.0 / 45000;
    std::vector<double> time(numberPoints);
    std::vector<double> data(numberPoints);
    for (size_t i = 0; i < numberPoints; ++i) {
        time[i] = period * static_cast<double>(i) / static_cast<double>(numberPoints - 1);
        double phase = time[i] / period;
        data[i] = phase < 0.4? 2 + phase / 0.4 : 3 - (phase - 0.4) / 0.6;
    }
    Waveform waveform;
    waveform.set_time(time);
    waveform.set_data(data);

    auto start = std::chrono::steady_clock::now();
    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform, 45000);
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    REQUIRE(sampledWaveform.get_data().size() == 524288);
    REQUIRE(seconds < 5);
    auto sampledTime = sampledWaveform.get_time().value();
    for (size_t i = 0; i < sampledTime.size(); i += 4099) {
        double phase = sampledTime[i] / period;
        double expected = phase < 0.4? 2 + phase / 0.4 : 3 - (phase - 0.4) / 0.6;
        CHECK_THAT(sampledWaveform.get_data()[i], Catch::Matchers::WithinAbs(expected, 1e-6));
    }
}

// Resuming the search must not change which segment an instant lands in: still the FIRST one that
// holds it, including on zero-length segments (repeated timestamps, as in FLYBACK_PRIMARY) where the
// value on the left of the step is the one taken. Compared against the full search it replaced.
TEST_CASE("Test_Sampled_Waveform_Matches_Full_Search", "[processor][waveform-processor][smoke-test]") {
    std::vector<double> time = {0, 0, 1e-6, 2.5e-6, 2.5e-6, 2.5e-6, 4e-6, 7e-6, 7e-6, 1e-5};
    std::vector<double> data = {-3, 5, 6, 7.5, -2.5, 1, 0.5, -1, 4, -3};
    Waveform waveform;
    waveform.set_time(time);
    waveform.set_data(data);

    auto sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform, 1e5, 1024);
    auto sampledTime = sampledWaveform.get_time().value();
    REQUIRE(sampledTime.size() == 1024);

    for (size_t i = 0; i < sampledTime.size(); ++i) {
        std::optional<double> expected;
        for (size_t k = 0; k + 1 < time.size() && !expected; ++k) {
            if (time[k + 1] == time[k]) {
                if (sampledTime[i] == time[k]) {
                    expected = data[k];
                }
            }
            else if (time[k] <= sampledTime[i] && sampledTime[i] <= time[k + 1]) {
                double proportion = (sampledTime[i] - time[k]) / (time[k + 1] - time[k]);
                expected = data[k] + (data[k + 1] - data[k]) * proportion;
            }
        }
        REQUIRE(expected);
        REQUIRE_THAT(sampledWaveform.get_data()[i], Catch::Matchers::WithinAbs(expected.value(), 1e-12));
    }
}

// ABT #1325: a caller can cap how many samples a dense waveform keeps. Capped, the samples are
// still exact on a piecewise-linear waveform; the cap never goes below what was requested, and a
// cap that the FFT could not take is refused.
TEST_CASE("Test_Sampled_Waveform_Capped_At_Maximum_Number_Points", "[processor][waveform-processor][smoke-test]") {
    const size_t numberPoints = 400000;
    const double period = 1.0 / 45000;
    std::vector<double> time(numberPoints);
    std::vector<double> data(numberPoints);
    for (size_t i = 0; i < numberPoints; ++i) {
        time[i] = period * static_cast<double>(i) / static_cast<double>(numberPoints - 1);
        double phase = time[i] / period;
        data[i] = phase < 0.4? 2 + phase / 0.4 : 3 - (phase - 0.4) / 0.6;
    }
    Waveform waveform;
    waveform.set_time(time);
    waveform.set_data(data);

    auto capped = WaveformProcessor::calculate_sampled_waveform(waveform, 45000, std::nullopt, 128, 8192);
    REQUIRE(capped.get_data().size() == 8192);
    auto cappedTime = capped.get_time().value();
    for (size_t i = 0; i < cappedTime.size(); i += 97) {
        double phase = cappedTime[i] / period;
        double expected = phase < 0.4? 2 + phase / 0.4 : 3 - (phase - 0.4) / 0.6;
        CHECK_THAT(capped.get_data()[i], Catch::Matchers::WithinAbs(expected, 1e-6));
    }

    auto requestedAboveCap = WaveformProcessor::calculate_sampled_waveform(waveform, 45000, 16384, 128, 8192);
    REQUIRE(requestedAboveCap.get_data().size() == 16384);

    REQUIRE_THROWS_AS(WaveformProcessor::calculate_sampled_waveform(waveform, 45000, std::nullopt, 128, 6000), std::invalid_argument);
}

// ABT #1325: calculate_instantaneous_power averaged |v*i| over the first 128 samples whatever the
// waveforms carried, so an imported waveform (8192 samples) was averaged over 1.6 % of its period.
// Square voltage (+24 V for 40 %, -12 V after) across a triangular current (2 -> 3 -> 2 A):
// mean |v*i| = 0.4 * 24 * 2.5 + 0.6 * 12 * 2.5 = 42 W.
TEST_CASE("Test_Instantaneous_Power_Averages_The_Whole_Period", "[processor][inputs][smoke-test]") {
    const double frequency = 45000;
    const double period = 1 / frequency;
    const size_t numberPoints = 8192;
    std::vector<double> time, voltage, current;
    for (size_t i = 0; i < numberPoints; ++i) {
        double t = period * static_cast<double>(i) / numberPoints;
        double phase = t / period;
        time.push_back(t);
        voltage.push_back(phase < 0.4? 24 : -12);
        current.push_back(phase < 0.4? 2 + phase / 0.4 : 3 - (phase - 0.4) / 0.6);
    }
    OperatingPointExcitation excitation;
    excitation.set_frequency(frequency);
    SignalDescriptor voltageSignal, currentSignal;
    Waveform voltageWaveform, currentWaveform;
    voltageWaveform.set_time(time);
    voltageWaveform.set_data(voltage);
    currentWaveform.set_time(time);
    currentWaveform.set_data(current);
    voltageSignal.set_waveform(voltageWaveform);
    currentSignal.set_waveform(currentWaveform);
    excitation.set_voltage(voltageSignal);
    excitation.set_current(currentSignal);

    REQUIRE_THAT(OpenMagnetics::Inputs::calculate_instantaneous_power(excitation), Catch::Matchers::WithinRel(42.0, 1e-3));
}

// --- ABT #1460 ------------------------------------------------------------
// Harmonics come from the exact Fourier series of the piecewise-linear knots, not from a
// 128-point DFT that folds everything above harmonic 64 back onto 1..63. A switched current
// decays as 1/k, so the aliased high harmonics read up to ~1.5x high and inflated litz
// proximity loss by 11-17% on flyback waveforms.
namespace {
std::vector<double> exact_amplitudes(const std::vector<double>& time, const std::vector<double>& data, size_t numberHarmonics) {
    auto coefficients = WaveformProcessor::calculate_exact_fourier_coefficients(time, data, numberHarmonics);
    std::vector<double> amplitudes{std::abs(coefficients[0])};
    for (size_t k = 1; k < coefficients.size(); ++k) {
        amplitudes.push_back(2 * std::abs(coefficients[k]));
    }
    return amplitudes;
}
double sinc_pi(double x) { return x == 0 ? 1.0 : std::sin(std::numbers::pi * x) / (std::numbers::pi * x); }
}  // namespace

TEST_CASE("Test_Exact_Harmonics_Analytic_Series", "[processor][waveform-processor][harmonics][smoke-test]") {
    const double period = 1e-5;
    const double amplitude = 3.0;
    const size_t numberHarmonics = 64;

    SECTION("Square wave: 4A/(pi k) at odd k, 0 at even k") {
        auto a = exact_amplitudes({0, period / 2, period / 2, period}, {amplitude, amplitude, -amplitude, -amplitude}, numberHarmonics);
        REQUIRE(a.size() == numberHarmonics + 1);
        CHECK(std::fabs(a[0]) < 1e-12);
        for (size_t k = 1; k <= numberHarmonics; ++k) {
            double expected = (k % 2 == 1) ? 4 * amplitude / (std::numbers::pi * static_cast<double>(k)) : 0.0;
            CHECK(std::fabs(a[k] - expected) < 1e-9);
        }
    }
    SECTION("Triangle: 8A/(pi^2 k^2) at odd k, 0 at even k") {
        auto a = exact_amplitudes({0, period / 2, period}, {-amplitude, amplitude, -amplitude}, numberHarmonics);
        CHECK(std::fabs(a[0]) < 1e-12);
        for (size_t k = 1; k <= numberHarmonics; ++k) {
            double kDouble = static_cast<double>(k);
            double expected = (k % 2 == 1) ? 8 * amplitude / (std::numbers::pi * std::numbers::pi * kDouble * kDouble) : 0.0;
            CHECK(std::fabs(a[k] - expected) < 1e-9);
        }
    }
    SECTION("Trapezoid: box (rise + top) convolved with box (rise)") {
        const double rise = 0.07 * period;
        const double top = 0.31 * period;
        auto a = exact_amplitudes({0, rise, rise + top, 2 * rise + top, period}, {0, amplitude, amplitude, 0, 0}, numberHarmonics);
        const double width = (rise + top) / period;
        CHECK(std::fabs(a[0] - amplitude * width) < 1e-9);
        for (size_t k = 1; k <= numberHarmonics; ++k) {
            double kDouble = static_cast<double>(k);
            double expected = 2 * amplitude * width * std::fabs(sinc_pi(kDouble * width) * sinc_pi(kDouble * rise / period));
            CHECK(std::fabs(a[k] - expected) < 1e-9);
        }
    }
}

TEST_CASE("Test_Exact_Harmonics_Flyback_DCM_Against_Fine_DFT", "[processor][waveform-processor][harmonics][smoke-test]") {
    // A DCM flyback-like current: a ramp that steps down, a step up into a falling ramp, an
    // idle stretch, and a step back at the period boundary.
    const double period = 1e-5;
    // The step sits on the reference's sample grid (0.375 = 24576 / 65536) so its mid-step
    // sample is exact; the kink at 0.85 T is off-grid.
    const std::vector<double> time = {0, 0.375 * period, 0.375 * period, 0.85 * period, period};
    const std::vector<double> data = {0.2, 2.0, 3.1, 0, 0};
    auto exact = exact_amplitudes(time, data, 64);

    // Reference: the plain DFT of 65536 point samples (data only, so the legacy path), taking
    // the mid-step value where a sample lands on a step. Its aliasing error is ~A/(pi N).
    const size_t numberSamples = 65536;
    std::vector<double> samples(numberSamples);
    for (size_t n = 0; n < numberSamples; ++n) {
        double t = period * static_cast<double>(n) / static_cast<double>(numberSamples);
        if (n == 0) {
            samples[n] = 0.5 * (data.back() + data.front());
        }
        else if (n == 24576) {
            samples[n] = 0.5 * (data[1] + data[2]);
        }
        else if (t < time[1]) {
            samples[n] = data[0] + (data[1] - data[0]) * t / time[1];
        }
        else if (t < time[3]) {
            samples[n] = data[2] + (data[3] - data[2]) * (t - time[2]) / (time[3] - time[2]);
        }
        else {
            samples[n] = 0;
        }
    }
    Waveform sampled;
    sampled.set_data(samples);
    auto reference = WaveformProcessor::calculate_harmonics_data(sampled, 1 / period, false);
    for (size_t k = 1; k <= 64; ++k) {
        INFO("harmonic " << k << ": exact " << exact[k] << ", 65536-point DFT " << reference.get_amplitudes()[k]);
        CHECK(std::fabs(exact[k] - reference.get_amplitudes()[k]) <= 1e-3 * reference.get_amplitudes()[k]);
    }
    CHECK(std::fabs(exact[0] - reference.get_amplitudes()[0]) <= 1e-3 * reference.get_amplitudes()[0]);

    // And the 128-point path this replaces was aliased well beyond that at the top harmonics.
    auto sampled128 = WaveformProcessor::calculate_sampled_waveform([&] { Waveform w; w.set_time(time); w.set_data(data); return w; }(), 1 / period);
    std::vector<double> samples128 = sampled128.get_data();
    Waveform dataOnly128;
    dataOnly128.set_data(samples128);
    auto aliased = WaveformProcessor::calculate_harmonics_data(dataOnly128, 1 / period, false);
    CHECK(aliased.get_amplitudes()[63] > 1.2 * reference.get_amplitudes()[63]);
}

TEST_CASE("Test_Exact_Harmonics_Knot_Invariance_And_Errors", "[processor][waveform-processor][harmonics][smoke-test]") {
    const double period = 1e-5;
    const std::vector<double> time = {0, 0, 0.3 * period, 0.3 * period, period};
    const std::vector<double> data = {0, 1.5, 4.0, 0, 0};
    auto reference = WaveformProcessor::calculate_exact_fourier_coefficients(time, data, 64);

    SECTION("Repeated knots (same instant, same value) change nothing") {
        const std::vector<double> repeatedTime = {0, 0, 0, 0.1 * period, 0.1 * period, 0.3 * period, 0.3 * period, 0.3 * period, 0.7 * period, period, period};
        // The ramp is also split at 0.1 T by a repeated knot on the line.
        const std::vector<double> repeatedData = {0, 0, 1.5, 1.5 + 2.5 * (0.1 / 0.3), 1.5 + 2.5 * (0.1 / 0.3), 4.0, 4.0, 0, 0, 0, 0};
        auto repeated = WaveformProcessor::calculate_exact_fourier_coefficients(repeatedTime, repeatedData, 64);
        REQUIRE(repeated.size() == reference.size());
        for (size_t k = 0; k < reference.size(); ++k) {
            CHECK(std::abs(repeated[k] - reference[k]) < 1e-12);
        }
    }
    SECTION("Uniform samples: the sinc^2-corrected DFT is the closed form of their interpolant") {
        Waveform knots;
        knots.set_time(time);
        knots.set_data(data);
        auto sampled = WaveformProcessor::calculate_sampled_waveform(knots, 1 / period);
        REQUIRE(WaveformProcessor::is_waveform_sampled(sampled));
        auto viaSamples = WaveformProcessor::calculate_harmonics_data(sampled, 1 / period, false);
        auto closedTime = sampled.get_time().value();
        auto closedData = sampled.get_data();
        closedTime.push_back(closedTime.back() + (closedTime[1] - closedTime[0]));
        closedData.push_back(closedData.front());
        auto closedForm = exact_amplitudes(closedTime, closedData, 64);
        REQUIRE(viaSamples.get_amplitudes().size() == 65);
        for (size_t k = 0; k <= 64; ++k) {
            CHECK(std::fabs(viaSamples.get_amplitudes()[k] - closedForm[k]) < 1e-9);
        }
    }
    SECTION("A knot waveform gets harmonics 0..64 from its knots") {
        Waveform knots;
        knots.set_time(time);
        knots.set_data(data);
        auto harmonics = WaveformProcessor::calculate_harmonics_data(knots, 1 / period, false);
        REQUIRE(harmonics.get_amplitudes().size() == 65);
        REQUIRE(harmonics.get_frequencies().size() == 65);
        CHECK(harmonics.get_frequencies()[64] == 64 / period);
        for (size_t k = 1; k <= 64; ++k) {
            CHECK(harmonics.get_amplitudes()[k] == 2 * std::abs(reference[k]));
        }
    }
    SECTION("Time running backwards throws") {
        CHECK_THROWS(WaveformProcessor::calculate_exact_fourier_coefficients({0, 0.5 * period, 0.4 * period, period}, {0, 1, 2, 0}, 64));
    }
    SECTION("128 points on an uneven time axis are knots, not samples") {
        // mas_autocomplete's compressed magnetizing current is exactly this: 128 uneven knots.
        std::vector<double> unevenTime, unevenData;
        for (size_t n = 0; n < 128; ++n) {
            double x = static_cast<double>(n) / 128.0;
            unevenTime.push_back(period * x * x);
            unevenData.push_back(std::sin(2 * std::numbers::pi * x));
        }
        Waveform uneven;
        uneven.set_time(unevenTime);
        uneven.set_data(unevenData);
        auto harmonics = WaveformProcessor::calculate_harmonics_data(uneven, 1 / period, false);
        auto exact = WaveformProcessor::calculate_exact_fourier_coefficients(unevenTime, unevenData, 64);
        REQUIRE(harmonics.get_amplitudes().size() == 65);
        CHECK(harmonics.get_amplitudes()[0] == std::abs(exact[0]));
        for (size_t k = 1; k <= 64; ++k) {
            CHECK(harmonics.get_amplitudes()[k] == 2 * std::abs(exact[k]));
        }
    }
}

TEST_CASE("Test_Rectangle_Label_Volt_Second_Balance_Is_Relative_To_Amplitude", "[processor][waveform-processor][smoke-test]") {
    // ABT #1641: the RECTANGULAR vs UNIPOLAR_RECTANGULAR decision compared a volt-second sum
    // (V*s) against an absolute tolerance equal to the period (s). Sampled at 128 points, a
    // balanced rectangle loses one sample of its high plateau (the sampler takes the pre-edge
    // value at t = 0), an imbalance of peak-to-peak * period / 128: at -250/+750 V and
    // 100 kHz that is 7.8e-5 V*s against a "tolerance" of 1e-5, and it was labelled
    // UNIPOLAR_RECTANGULAR with its minimum (-250 V) as offset. A genuinely unipolar 0/1 V
    // one, 2.5e-6 V*s, passed as balanced.
    const double frequency = 100000;
    const double period = 1 / frequency;
    const double dutyCycle = 0.25;
    for (double amplitude : {1e-3, 1.0, 1000.0, 1e6}) {
        INFO("peak-to-peak " << amplitude);
        Waveform balanced;
        balanced.set_data({-0.25 * amplitude, 0.75 * amplitude, 0.75 * amplitude, -0.25 * amplitude, -0.25 * amplitude});
        balanced.set_time(std::vector<double>{0, 0, dutyCycle * period, dutyCycle * period, period});
        auto sampledBalanced = WaveformProcessor::calculate_sampled_waveform(balanced, frequency);
        REQUIRE(WaveformProcessor::is_waveform_sampled(sampledBalanced));
        CHECK(WaveformProcessor::try_guess_waveform_label(balanced) == WaveformLabel::RECTANGULAR);
        CHECK(WaveformProcessor::try_guess_waveform_label(sampledBalanced) == WaveformLabel::RECTANGULAR);
        auto processedBalanced = WaveformProcessor::calculate_processed_data(sampledBalanced, frequency);
        CHECK(processedBalanced.get_label() == WaveformLabel::RECTANGULAR);
        CHECK(processedBalanced.get_offset() == 0);

        Waveform unipolar;
        unipolar.set_data({0, amplitude, amplitude, 0, 0});
        unipolar.set_time(std::vector<double>{0, 0, dutyCycle * period, dutyCycle * period, period});
        auto sampledUnipolar = WaveformProcessor::calculate_sampled_waveform(unipolar, frequency);
        CHECK(WaveformProcessor::try_guess_waveform_label(unipolar) == WaveformLabel::UNIPOLAR_RECTANGULAR);
        CHECK(WaveformProcessor::try_guess_waveform_label(sampledUnipolar) == WaveformLabel::UNIPOLAR_RECTANGULAR);
    }
}

