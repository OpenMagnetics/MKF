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
#include <cmath>
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
