#include "InputsWrapper.h"

#include "Constants.h"
#include "Utils.h"
#include "json.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <string>
#include <valarray>
#include <vector>

using json = nlohmann::json;

namespace OpenMagnetics {

// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
// !!! Warning : in some cases this code make result different from not optimased version above (need to fix bug)
// The bug is now fixed @2017/05/30
// from https://rosettacode.org/wiki/Fast_Fourier_transform#C++
void fft(std::vector<std::complex<double>>& x) {
    // DFT
    size_t N = x.size(), k = N, n;
    double thetaT = std::numbers::pi / N;
    std::complex<double> phiT = std::complex<double>(cos(thetaT), -sin(thetaT)), T;
    while (k > 1) {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0L;
        for (size_t l = 0; l < k; l++) {
            for (size_t a = l; a < N; a += n) {
                size_t b = a + k;
                std::complex<double> t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    // Decimate
    size_t m = (size_t) log2(N);
    for (size_t a = 0; a < N; a++) {
        size_t b = a;
        // Reverse bits
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a) {
            std::complex<double> t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }
}

double InputsWrapper::try_get_duty_cycle(Waveform waveform, double frequency) {
    if (waveform.get_time()) {
        std::vector<double> data = waveform.get_data();
        std::vector<double> time = waveform.get_time().value();
        if (data.size() == 3) {
            if (roundFloat(1 / (time[2] - time[0]), 9) == frequency)
                if (data[0] == data[2])
                    return roundFloat((time[1] - time[0]) / (time[2] - time[0]), 3);
        }
        else if (data.size() == 5) {
            if (roundFloat(1 / (time[4] - time[0]), 9) == frequency)
                if (data[0] == data[1])
                    if (time[1] == time[2])
                        if (data[2] == data[3])
                            if (data[0] == data[4]) {
                                return roundFloat((time[1] - time[0]) / (time[4] - time[0]), 3);
                            }
        }
        else if (data.size() == 10) {
            if (roundFloat(1 / (time[9] - time[0]), 9) == frequency)
                if (data[0] == data[1])
                    if (time[1] == time[2])
                        if (data[2] == data[3])
                            if (time[3] == time[4])
                                if (data[4] == data[5])
                                    if (time[5] == time[6])
                                        if (data[6] == data[7])
                                            if (time[7] == time[8])
                                                if (data[8] == data[9])
                                                    if (data[0] == data[9])
                                                        return roundFloat((time[1] - time[0]) / (time[7] - time[0]), 3);
        }
    }

    auto sampledWaveform = InputsWrapper::get_sampled_waveform(waveform, frequency);
    std::vector<double> data = sampledWaveform.get_data();
    std::vector<double> diff_data;
    std::vector<double> diff_diff_data;

    data.push_back(data[0]);
    data.push_back(data[1]);
    for (size_t i = 0; i < data.size() - 1; ++i) {
        diff_data.push_back(roundFloat(data[i + 1] - data[i], 9));
    }
    for (size_t i = 0; i < diff_data.size() - 1; ++i) {
        diff_diff_data.push_back(fabs(roundFloat(diff_data[i + 1] - diff_data[i], 9)));
    }

    auto constants = Constants();
    auto maximum = std::max_element(diff_diff_data.rbegin(), diff_diff_data.rend()).base();
    auto maximum_index = std::distance(diff_diff_data.begin(), std::prev(maximum));
    auto dutyCycle = roundFloat((maximum_index + 1) / constants.number_points_samples_waveforms, 2);
    return dutyCycle;
}
std::vector<double> linear_spaced_array(double startingValue, double endingValue, std::size_t numberPoints) {
    double h = (endingValue - startingValue) / static_cast<double>(numberPoints);
    std::vector<double> xs(numberPoints);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = startingValue; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}

// In case the waveform comes defined with processed data only, we create a MAS format waveform from it, as the rest of
// the code depends on it.
SignalDescriptor InputsWrapper::standarize_waveform(SignalDescriptor parameter, double frequency) {
    SignalDescriptor standardized_parameter(parameter);
    if (!parameter.get_waveform()) {
        Waveform waveform;
        double period = 1 / frequency;
        auto processed = parameter.get_processed().value();

        auto peakToPeak = processed.get_peak_to_peak().value();
        auto offset = processed.get_offset();
        double dutyCycle = 0.5;
        if (processed.get_duty_cycle()) {
            dutyCycle = processed.get_duty_cycle().value();
        }
        std::vector<double> data;
        std::vector<double> time;
        switch (processed.get_label()) {
            case OpenMagnetics::WaveformLabel::TRIANGULAR:
                data = {-peakToPeak / 2 + offset, peakToPeak / 2 + offset, -peakToPeak / 2 + offset};
                time = {0., dutyCycle * period, period};
                break;
            case OpenMagnetics::WaveformLabel::RECTANGULAR: {
                double max = peakToPeak * (1 - dutyCycle);
                double min = -peakToPeak * dutyCycle;
                double dc = dutyCycle * period;
                data = {max, max, min, min, max};
                time = {0, dc, dc, period, period};
                break;
            }
            case OpenMagnetics::WaveformLabel::BIPOLAR_RECTANGULAR: {
                double max = +peakToPeak / 2;
                double min = -peakToPeak / 2;
                double dc = dutyCycle * period;
                data = {0, 0, max, max, 0, 0, min, min, 0, 0};
                time = {0,
                        0.25 * period - dc / 2,
                        0.25 * period - dc / 2,
                        0.25 * period + dc / 2,
                        0.25 * period + dc / 2,
                        0.75 * period - dc / 2,
                        0.75 * period - dc / 2,
                        0.75 * period + dc / 2,
                        0.75 * period + dc / 2,
                        1};
                break;
            }
            case OpenMagnetics::WaveformLabel::SINUSOIDAL: {
                auto constants = Constants();
                for (size_t i = 0; i < constants.number_points_samples_waveforms; ++i) {
                    double angle = i * 2 * std::numbers::pi / constants.number_points_samples_waveforms;
                    time.push_back(angle * period);
                    data.push_back((sin(angle) * peakToPeak / 2) + offset);
                }
                break;
            }
            default:
                break;
        }

        waveform.set_data(data);
        if (processed.get_label() == OpenMagnetics::WaveformLabel::SINUSOIDAL) {
            waveform.set_time(std::nullopt);
        }
        else {
            waveform.set_time(time);
        }
        standardized_parameter.set_waveform(waveform);
    }

    return standardized_parameter;
}

bool is_size_power_of_2(std::vector<double> data) {
    return data.size() > 0 && ((data.size() & (data.size() - 1)) == 0);
}

bool InputsWrapper::is_waveform_sampled(Waveform waveform) {
    auto constants = Constants();
    if (!waveform.get_time()) {
        return false;
    }
    else {
        return waveform.get_data().size() == constants.number_points_samples_waveforms;
    }

}

Waveform InputsWrapper::get_sampled_waveform(Waveform waveform, double frequency) {
    auto constants = Constants();

    std::vector<double> time;
    auto data = waveform.get_data();
    if (!waveform.get_time()) { // This means the waveform is equidistant
        time = linear_spaced_array(0, 1. / frequency, data.size());
    }
    else {
        time = waveform.get_time().value();
        // Check if we need to guess the frequency from the time vector
        if (frequency == 0) {
            frequency = 1 / (time[waveform.get_data().size() - 1] - time[0]);
        }
    }

    auto sampledTime = linear_spaced_array(0, 1. / frequency, constants.number_points_samples_waveforms);

    std::vector<double> sampledData;

    for (int i = 0; i < constants.number_points_samples_waveforms; i++) {
        for (size_t interpIndex = 0; interpIndex < data.size() - 1; interpIndex++) {
            if (time[interpIndex] <= sampledTime[i] && sampledTime[i] <= time[interpIndex + 1]) {
                double proportion = (sampledTime[i] - time[interpIndex]) / (time[interpIndex + 1] - time[interpIndex]);
                double interpPoint = std::lerp(data[interpIndex], data[interpIndex + 1], proportion);
                sampledData.push_back(interpPoint);
                break;
            }
        }
    }

    Waveform sampledWaveform;
    sampledWaveform.set_data(sampledData);
    sampledWaveform.set_time(sampledTime);
    return sampledWaveform;
}

SignalDescriptor InputsWrapper::get_induced_voltage(OperatingPointExcitation& excitation,
                                                            double magnetizingInductance) {
    auto constants = Constants();
    Waveform sourceWaveform = excitation.get_current().value().get_waveform().value();
    std::vector<double> source = sourceWaveform.get_data();
    bool isWaveformSampled =
        !bool(sourceWaveform.get_time()) || source.size() == constants.number_points_samples_waveforms;
    std::vector<double> time = sourceWaveform.get_time().value();
    std::vector<double> derivative;
    std::vector<double> derivativeTime;
    std::vector<double> voltageData;
    Waveform voltageWaveform;
    SignalDescriptor voltageSignalDescriptor;
    Waveform resultWaveform(sourceWaveform);

    if (isWaveformSampled) {
        source.push_back(source[0]);
        double difference = time[time.size() - 1] - time[time.size() - 2];
        time.push_back(time[time.size() - 1] + difference);
    }
    else {
        source.push_back(source[1]);
        time.push_back(time[time.size() - 1] + time[1]);
    }


    std::adjacent_difference(source.begin(), source.end(), source.begin());
    derivative = std::vector<double>(source.begin() + 1, source.end());
    std::adjacent_difference(time.begin(), time.end(), time.begin());
    derivativeTime = std::vector<double>(time.begin() + 1, time.end());

    for (size_t i = 0; i < derivative.size(); ++i) {
        voltageData.push_back(magnetizingInductance * derivative[i] / derivativeTime[i]);
    }

    voltageWaveform.set_data(voltageData);
    voltageWaveform.set_time(time);
    voltageSignalDescriptor.set_waveform(voltageWaveform);
    auto sampledWaveform = InputsWrapper::get_sampled_waveform(voltageWaveform, excitation.get_frequency());
    voltageSignalDescriptor.set_harmonics(get_harmonics_data(sampledWaveform, excitation.get_frequency()));
    voltageSignalDescriptor.set_processed(
        get_processed_data(voltageSignalDescriptor, sampledWaveform, true, true));

    resultWaveform.set_data(derivative);
    return voltageSignalDescriptor;
}

SignalDescriptor InputsWrapper::add_offset_to_excitation(SignalDescriptor signalDescriptor,
                                                                 double offset,
                                                                 double frequency) {
    auto waveform = signalDescriptor.get_waveform().value();
    std::vector<double> source = waveform.get_data();
    std::vector<double> modified_data;

    for (auto& elem : source) {
        modified_data.push_back(elem + offset);
    }

    waveform.set_data(modified_data);
    signalDescriptor.set_waveform(waveform);
    auto sampledWaveform = InputsWrapper::get_sampled_waveform(waveform, frequency);
    signalDescriptor.set_harmonics(get_harmonics_data(sampledWaveform, frequency));
    signalDescriptor.set_processed(get_processed_data(signalDescriptor, sampledWaveform, true, true));
    return signalDescriptor;
}

Waveform get_magnetizing_current_waveform(Waveform sourceWaveform,
                                          double frequency,
                                          double magnetizingInductance,
                                          double dcCurrent = 0) {
    std::vector<double> source = sourceWaveform.get_data();
    std::vector<double> integration;
    Waveform resultWaveform(sourceWaveform);

    double integral = 0;
    integration.push_back(integral);
    double timePerPoint = 1 / frequency / source.size();
    for (auto& point : source) {
        integral += point / magnetizingInductance * timePerPoint;
        integration.push_back(integral);
    }

    integration = std::vector<double>(integration.begin(), integration.end() - 1);

    double integrationAverage =
        std::accumulate(std::begin(integration), std::end(integration), 0.0) / integration.size();
    for (size_t i = 0; i < integration.size(); ++i) {
        integration[i] = integration[i] - integrationAverage + dcCurrent;
    }

    resultWaveform.set_data(integration);
    return resultWaveform;
}

SignalDescriptor InputsWrapper::reflect_waveform(SignalDescriptor primarySignalDescriptor,
                                                         double ratio) {
    SignalDescriptor reflected_waveform;
    auto primaryWaveform = primarySignalDescriptor.get_waveform().value();
    Waveform waveform = primaryWaveform;

    waveform.set_data({});
    for (auto& datum : primaryWaveform.get_data()) {
        waveform.get_mutable_data().push_back(datum * ratio);
    }
    // waveform.set_time(primarySignalDescriptor.get_waveform().value().get_time());
    reflected_waveform.set_waveform(waveform);

    return reflected_waveform;
}

std::pair<bool, std::string> InputsWrapper::check_integrity() {
    auto operatingPoints = get_mutable_operating_points();
    auto turnsRatios = get_design_requirements().get_turns_ratios();
    auto magnetizingInductance = resolve_dimensional_values(get_design_requirements().get_magnetizing_inductance());
    std::pair<bool, std::string> result;
    result.first = true;
    result.second = "";

    for (auto& operating_point : operatingPoints) {
        if (operating_point.get_excitations_per_winding().size() == 0) {
            result.first = false;
            throw std::invalid_argument("Missing excitation for primary");
        }
    }

    for (size_t i = 0; i < operatingPoints.size(); ++i) {
        std::vector<OperatingPointExcitation> processed_excitations_per_winding;


        for (auto& excitation : operatingPoints[i].get_mutable_excitations_per_winding()) {
            // Here we processed this excitation voltage
            if (excitation.get_voltage()) {
                auto voltage_excitation = excitation.get_voltage().value();
                voltage_excitation = standarize_waveform(voltage_excitation, excitation.get_frequency());
                excitation.set_voltage(voltage_excitation);
            }
            // Here we processed this excitation current
            if (excitation.get_current()) {
                auto current_excitation = excitation.get_current().value();
                current_excitation = standarize_waveform(current_excitation, excitation.get_frequency());
                excitation.set_current(current_excitation);
            }
            else {
                auto voltageWaveform = excitation.get_voltage().value().get_waveform().value();
                auto sampledWaveform = get_sampled_waveform(voltageWaveform, excitation.get_frequency());
                excitation.set_current(
                    get_magnetizing_current(excitation, sampledWaveform, magnetizingInductance));
            }
            processed_excitations_per_winding.push_back(excitation);
        }
        operatingPoints[i].set_excitations_per_winding(processed_excitations_per_winding);
    }

    for (size_t i = 0; i < operatingPoints.size(); ++i) {

        if (turnsRatios.size() > operatingPoints[i].get_excitations_per_winding().size() - 1) {
            if (turnsRatios.size() == 1 && operatingPoints[i].get_excitations_per_winding().size() == 1) {
                // We are missing excitation only for secondary
                for (size_t turnRatioIndex = 0; turnRatioIndex < turnsRatios.size(); ++turnRatioIndex) {
                    if (turnRatioIndex >= operatingPoints[i].get_excitations_per_winding().size() - 1) {
                        double turnRatio = resolve_dimensional_values(turnsRatios[turnRatioIndex]);
                        auto excitationOfPrimaryWinding = operatingPoints[i].get_excitations_per_winding()[0];
                        OperatingPointExcitation excitationOfThisWinding(excitationOfPrimaryWinding);

                        excitationOfThisWinding.set_voltage(
                            reflect_waveform(excitationOfPrimaryWinding.get_voltage().value(), 1 / turnRatio));


                        excitationOfThisWinding.set_current(
                            reflect_waveform(excitationOfPrimaryWinding.get_current().value(), turnRatio));
                        operatingPoints[i].get_mutable_excitations_per_winding().push_back(excitationOfThisWinding);
                    }
                }
                result.second = "Had to create the excitations of some windings based on primary";
            }
            else {
                throw std::invalid_argument("Missing excitation for more than one secondary. Only one can be guessed");
            }
        }
    }

    set_operating_points(operatingPoints);

    return result;
}

Processed InputsWrapper::get_processed_data(SignalDescriptor excitation,
                                            Waveform sampledWaveform,
                                            bool force = false,
                                            bool includeAdvancedData = true) {
    auto constants = Constants();
    Processed processed;
    std::vector<double> dataToProcess;
    auto harmonics = excitation.get_harmonics().value();
    auto sampledDataToProcess = sampledWaveform.get_data();

    if (sampledWaveform.get_time() && sampledWaveform.get_data().size() < constants.number_points_samples_waveforms) {
        sampledDataToProcess = get_sampled_waveform(sampledWaveform, 0).get_data();
    }

    if (!excitation.get_processed() || force) {
        auto waveform = excitation.get_waveform().value();
        dataToProcess = waveform.get_data();

        for (size_t i = 0; i < dataToProcess.size(); ++i) {
            if (std::isnan(dataToProcess[i])) {
                throw std::invalid_argument("Waveform data contains NaN");
            }
        }

        std::string labelString;
        if (waveform.get_ancillary_label()) {
            labelString = waveform.get_ancillary_label().value();
        }
        else {
            labelString = "custom";
        }
        std::transform(labelString.begin(), labelString.end(), labelString.begin(),
                       [](unsigned char c) { return std::toupper(c); });
        std::transform(labelString.begin(), labelString.end(), labelString.begin(),
                       [](unsigned char c) { return (c=='-' || c==' ')? '_' : c; });
        std::optional<WaveformLabel> label;
        try {
            label = magic_enum::enum_cast<WaveformLabel>(labelString);
            processed.set_label(label.value());

            processed.set_offset(std::accumulate(std::begin(sampledDataToProcess), std::end(sampledDataToProcess), 0.0) /
                                 sampledDataToProcess.size());
            processed.set_peak_to_peak(*max_element(dataToProcess.begin(), dataToProcess.end()) -
                                       *min_element(dataToProcess.begin(), dataToProcess.end()));
            processed.set_peak(*max_element(dataToProcess.begin(), dataToProcess.end()));
        }
        catch (...) {
            throw std::invalid_argument("Unknown excitation label: " + labelString);
        }
    }
    else {
        processed = excitation.get_processed().value();
    }

    if (includeAdvancedData) {
        if (!processed.get_effective_frequency() || force) {
            double effectiveFrequency = 0;
            std::vector<double> dividend;
            std::vector<double> divisor;

            for (size_t i = 0; i < harmonics.get_amplitudes().size(); ++i) {
                double dataSquared = pow(harmonics.get_amplitudes()[i], 2);
                double timeSquared = pow(harmonics.get_frequencies()[i], 2);
                dividend.push_back(dataSquared * timeSquared);
                divisor.push_back(dataSquared);
            }
            double sumDivisor = std::reduce(divisor.begin(), divisor.end());
            if (sumDivisor > 0)
                effectiveFrequency = sqrt(std::reduce(dividend.begin(), dividend.end()) / sumDivisor);

            processed.set_effective_frequency(effectiveFrequency);
        }
        if (!processed.get_ac_effective_frequency() || force) {
            double effectiveFrequency = 0;
            std::vector<double> dividend;
            std::vector<double> divisor;

            for (size_t i = 1; i < harmonics.get_amplitudes().size(); ++i) {
                double dataSquared = pow(harmonics.get_amplitudes()[i], 2);
                double timeSquared = pow(harmonics.get_frequencies()[i], 2);
                dividend.push_back(dataSquared * timeSquared);
                divisor.push_back(dataSquared);
            }
            double sumDivisor = std::reduce(divisor.begin(), divisor.end());
            if (sumDivisor > 0)
                effectiveFrequency = sqrt(std::reduce(dividend.begin(), dividend.end()) / sumDivisor);

            processed.set_ac_effective_frequency(effectiveFrequency);
        }

        if (!processed.get_rms() || force) {
            double rms = 0.0;
            for (int i = 0; i < int(sampledDataToProcess.size()); ++i) {
                rms += sampledDataToProcess[i] * sampledDataToProcess[i];
            }
            rms /= sampledDataToProcess.size();
            rms = sqrt(rms);
            processed.set_rms(rms);
        }
        if (!processed.get_thd() || force) {
            std::vector<double> dividend;
            double divisor = harmonics.get_amplitudes()[1];
            double thd = 0;
            for (size_t i = 2; i < harmonics.get_amplitudes().size(); ++i) {
                dividend.push_back(pow(harmonics.get_amplitudes()[i], 2));
            }

            if (divisor > 0)
                thd = sqrt(std::reduce(dividend.begin(), dividend.end())) / divisor;
            processed.set_thd(thd);
        }
    }

    return processed;
}

Harmonics InputsWrapper::get_harmonics_data(Waveform waveform, double frequency) {
    auto dataToProcess = waveform.get_data();
    Harmonics harmonics;

    std::vector<std::complex<double>> data;
    for (std::size_t i = 0; i < dataToProcess.size(); ++i) {
        data.emplace_back(dataToProcess[i]);
    }

    if (data.size() > 0 && ((data.size() & (data.size() - 1)) != 0)) {

        throw std::invalid_argument("Data vector size is not a power of 2");
    }
    fft(data);

    harmonics.get_mutable_amplitudes().push_back(abs(data[0] / static_cast<double>(data.size())));
    for (size_t i = 1; i < data.size() / 2; ++i) {
        harmonics.get_mutable_amplitudes().push_back(abs(2. * data[i] / static_cast<double>(data.size())));
    }
    for (size_t i = 0; i < data.size() / 2; ++i) {
        harmonics.get_mutable_frequencies().push_back(frequency * i);
    }

    return harmonics;
}

SignalDescriptor InputsWrapper::get_magnetizing_current(OperatingPointExcitation& excitation,
                                                                Waveform sampledWaveform,
                                                                double magnetizingInductance) {
    double dcCurrent = 0;
    if (excitation.get_current()) {
        if (!excitation.get_current().value().get_processed()) {
            auto currentExcitation = excitation.get_current().value();
            auto currentExcitationWaveform = currentExcitation.get_waveform().value();
            auto sampledCurrentWaveform = get_sampled_waveform(currentExcitationWaveform, excitation.get_frequency());
            currentExcitation.set_harmonics(get_harmonics_data(sampledCurrentWaveform, excitation.get_frequency()));
            currentExcitation.set_processed(get_processed_data(currentExcitation, sampledCurrentWaveform));
            excitation.set_current(currentExcitation);
        }
        dcCurrent = excitation.get_current().value().get_processed().value().get_offset();
    }

    auto sampledMagnetizingCurrentWaveform =
        get_magnetizing_current_waveform(sampledWaveform, excitation.get_frequency(), magnetizingInductance, dcCurrent);
    SignalDescriptor magnetizing_current_excitation;

    magnetizing_current_excitation.set_waveform(sampledMagnetizingCurrentWaveform);
    magnetizing_current_excitation.set_harmonics(
        get_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
    magnetizing_current_excitation.set_processed(
        get_processed_data(magnetizing_current_excitation, sampledMagnetizingCurrentWaveform));
    return magnetizing_current_excitation;
}

SignalDescriptor InputsWrapper::get_magnetizing_current(OperatingPointExcitation& excitation,
                                                                double magnetizingInductance) {
    double dcCurrent = 0;
    if (excitation.get_current()) {
        dcCurrent = excitation.get_current().value().get_processed().value().get_offset();
    }

    auto voltage_excitation = excitation.get_voltage().value();
    voltage_excitation = standarize_waveform(voltage_excitation, excitation.get_frequency());
    auto waveform = voltage_excitation.get_waveform().value();
    auto sampledWaveform = get_sampled_waveform(waveform, excitation.get_frequency());

    auto sampledMagnetizingCurrentWaveform =
        get_magnetizing_current_waveform(sampledWaveform, excitation.get_frequency(), magnetizingInductance, dcCurrent);
    SignalDescriptor magnetizing_current_excitation;

    magnetizing_current_excitation.set_waveform(sampledMagnetizingCurrentWaveform);
    magnetizing_current_excitation.set_harmonics(
        get_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
    magnetizing_current_excitation.set_processed(
        get_processed_data(magnetizing_current_excitation, sampledMagnetizingCurrentWaveform));
    return magnetizing_current_excitation;
}

OperatingPoint InputsWrapper::process_operating_point(OperatingPoint operatingPoint, double magnetizingInductance) {
    std::vector<OperatingPointExcitation> processed_excitations_per_winding;
    for (auto& excitation : operatingPoint.get_mutable_excitations_per_winding()) {
        // Here we processed this excitation current
        if (excitation.get_current()) {
            auto current_excitation = excitation.get_current().value();
            current_excitation = standarize_waveform(current_excitation, excitation.get_frequency());
            auto waveform = current_excitation.get_waveform().value();
            auto sampledWaveform = get_sampled_waveform(waveform, excitation.get_frequency());
            current_excitation.set_harmonics(get_harmonics_data(sampledWaveform, excitation.get_frequency()));
            current_excitation.set_processed(get_processed_data(current_excitation, sampledWaveform));
            excitation.set_current(current_excitation);
        }
        // Here we processed this excitation voltage
        if (excitation.get_voltage()) {
            auto voltage_excitation = excitation.get_voltage().value();
            voltage_excitation = standarize_waveform(voltage_excitation, excitation.get_frequency());
            auto waveform = voltage_excitation.get_waveform().value();
            auto sampledWaveform = get_sampled_waveform(waveform, excitation.get_frequency());
            voltage_excitation.set_harmonics(get_harmonics_data(sampledWaveform, excitation.get_frequency()));
            voltage_excitation.set_processed(get_processed_data(voltage_excitation, sampledWaveform));
            excitation.set_voltage(voltage_excitation);

            if (!excitation.get_magnetizing_current()) {
                excitation.set_magnetizing_current(
                    get_magnetizing_current(excitation, sampledWaveform, magnetizingInductance));
            }
        }
        processed_excitations_per_winding.push_back(excitation);
    }
    operatingPoint.set_excitations_per_winding(processed_excitations_per_winding);
    return operatingPoint;
}

void InputsWrapper::process_waveforms() {
    auto operating_points = get_mutable_operating_points();
    std::vector<OperatingPoint> processed_operating_points;
    for (auto& operating_point : operating_points) {
        processed_operating_points.push_back(process_operating_point(
            operating_point, resolve_dimensional_values(get_design_requirements().get_magnetizing_inductance())));
    }
    set_operating_points(processed_operating_points);
}

InputsWrapper InputsWrapper::create_quick_operating_point(double frequency,
                                                          double magnetizingInductance,
                                                          double temperature,
                                                          WaveformLabel waveShape,
                                                          double peakToPeak,
                                                          double dutyCycle,
                                                          double dcCurrent,
                                                          std::vector<double> turnsRatios) {
    json inputsJson;

    inputsJson["operatingPoints"] = json::array();
    json operatingPointJson = json();
    operatingPointJson["name"] = "Nominal";
    operatingPointJson["conditions"] = json();
    operatingPointJson["conditions"]["ambientTemperature"] = temperature;

    json windingExcitation = json();
    windingExcitation["winding"] = "primary";
    windingExcitation["frequency"] = frequency;
    windingExcitation["voltage"]["processed"]["dutyCycle"] = dutyCycle;
    windingExcitation["voltage"]["processed"]["label"] = waveShape;
    windingExcitation["voltage"]["processed"]["offset"] = 0;
    windingExcitation["voltage"]["processed"]["peakToPeak"] = peakToPeak;
    operatingPointJson["excitationsPerWinding"] = json::array();
    operatingPointJson["excitationsPerWinding"].push_back(windingExcitation);
    inputsJson["operatingPoints"].push_back(operatingPointJson);

    inputsJson["designRequirements"] = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["minimum"] = magnetizingInductance * 0.8;
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = magnetizingInductance;
    inputsJson["designRequirements"]["magnetizingInductance"]["maximum"] = magnetizingInductance * 1.2;
    json turnsRatiosJson = json::array();
    for (auto& turnsRatio : turnsRatios) {
        json turnsRatioJson = json();
        turnsRatioJson["nominal"] = turnsRatio;
        turnsRatiosJson.push_back(turnsRatioJson);
    }
    inputsJson["designRequirements"]["turnsRatios"] = turnsRatiosJson;

    InputsWrapper inputs(inputsJson);

    auto voltage =
        inputs.get_mutable_operating_points()[0].get_mutable_excitations_per_winding()[0].get_voltage().value();

    auto operatingPoint = inputs.get_mutable_operating_points()[0];

    std::vector<OperatingPointExcitation> excitationsWithDcCurrent;
    for (auto& excitation : operatingPoint.get_mutable_excitations_per_winding()) {
        auto magnetizingCurrentSignalDescriptor = get_magnetizing_current(excitation, magnetizingInductance);
        excitation.set_current(
            InputsWrapper::add_offset_to_excitation(magnetizingCurrentSignalDescriptor, dcCurrent, frequency));
        excitation.set_magnetizing_current(
            InputsWrapper::add_offset_to_excitation(magnetizingCurrentSignalDescriptor, dcCurrent, frequency));
        excitationsWithDcCurrent.push_back(excitation);
    }
    operatingPoint.set_excitations_per_winding(excitationsWithDcCurrent);
    inputs.set_operating_point_by_index(operatingPoint, 0);

    return inputs;
}

OperatingPoint InputsWrapper::get_operating_point(size_t index) {
    return get_mutable_operating_points()[index];
}

OperatingPointExcitation InputsWrapper::get_winding_excitation(size_t operatingPointIndex, size_t windingIndex) {
    return get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex];
}

OperatingPointExcitation InputsWrapper::get_primary_excitation(size_t operatingPointIndex) {
    return get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0];
}

OperatingPointExcitation InputsWrapper::get_primary_excitation(OperatingPoint operatingPoint) {
    return operatingPoint.get_mutable_excitations_per_winding()[0];
}

void InputsWrapper::make_waveform_size_power_of_two(OperatingPoint* operatingPoint) {
    OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(*operatingPoint);
    double frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();

    // TODO: iterate over windings here in the future

    if (excitation.get_current()) {
        auto current = operatingPoint->get_excitations_per_winding()[0].get_current().value();
        auto currentWaveform = current.get_waveform().value();
        if (!is_size_power_of_2(currentWaveform.get_data())) {
            auto currentSampledWaveform = InputsWrapper::get_sampled_waveform(currentWaveform, frequency);
            current.set_waveform(currentSampledWaveform);
            current.set_harmonics(get_harmonics_data(currentSampledWaveform, frequency));
            current.set_processed(get_processed_data(current, currentSampledWaveform, true));
            operatingPoint->get_mutable_excitations_per_winding()[0].set_current(current);


            currentWaveform = current.get_waveform().value();

            currentWaveform = operatingPoint->get_excitations_per_winding()[0].get_current().value().get_waveform().value();
        }
    }
    if (excitation.get_voltage()) {
        auto voltage = operatingPoint->get_excitations_per_winding()[0].get_voltage().value();
        auto voltageWaveform = voltage.get_waveform().value();
        if (!is_size_power_of_2(voltageWaveform.get_data())) {
            auto voltageSampledWaveform = InputsWrapper::get_sampled_waveform(voltageWaveform, frequency);
            voltage.set_waveform(voltageSampledWaveform);
            operatingPoint->get_mutable_excitations_per_winding()[0].set_voltage(voltage);
        }
    }
}

double InputsWrapper::get_waveform_coefficient(OperatingPoint* operatingPoint) {
    auto constants = Constants();
    OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(*operatingPoint);
    double frequency = excitation.get_frequency();
    Waveform sampledWaveform = excitation.get_voltage().value().get_waveform().value();

    if (sampledWaveform.get_time() && sampledWaveform.get_data().size() < constants.number_points_samples_waveforms) {
        sampledWaveform = get_sampled_waveform(sampledWaveform, frequency);
    }

    std::vector<double> source = sampledWaveform.get_data();
    std::vector<double> integration;

    double integral = 0;
    integration.push_back(integral);
    double time_per_point = 1 / frequency / source.size();
    source = std::vector<double>(source.begin(), source.end() - source.size() / 2);

    for (size_t i = 0; i < source.size() - 1; i++) {
        integral += fabs(source[i + 1] - source[i]) / 2. + source[i];
    }
    integral *= time_per_point;

    double voltageRms = excitation.get_voltage().value().get_processed().value().get_rms().value();

    double waveformCoefficient = 2 * voltageRms / (frequency * integral);
    return waveformCoefficient;
}

} // namespace OpenMagnetics
