#include "Settings.h"
#include "InputsWrapper.h"

#include "Constants.h"
#include "Defaults.h"
#include "Utils.h"
#include "json.hpp"

#include <algorithm>
#include <cctype>
#include <cfloat>
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
#include <list>


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

double calculate_waveform_average(Waveform waveform) {
    double integration = 0;
    double period = waveform.get_time()->back() - waveform.get_time()->front();
    for (size_t i = 0; i < waveform.get_data().size() - 1; ++i)
    {
        double area = (waveform.get_data()[i + 1] + waveform.get_data()[i]) / 2 * (waveform.get_time().value()[i + 1] - waveform.get_time().value()[i]);
        integration += area;
    }
    return integration / period;
}

Waveform multiply_waveform(Waveform waveform, double scalarValue){
    Waveform scaledWaveform = waveform;
    scaledWaveform.get_mutable_data().clear();
    for (auto& datum : waveform.get_data()) {
        scaledWaveform.get_mutable_data().push_back(datum * scalarValue);
    }
    return scaledWaveform;
}

Waveform sum_waveform(Waveform waveform, double scalarValue){
    Waveform scaledWaveform = waveform;
    scaledWaveform.get_mutable_data().clear();
    for (auto& datum : waveform.get_data()) {
        scaledWaveform.get_mutable_data().push_back(datum + scalarValue);
    }
    return scaledWaveform;
}

double InputsWrapper::try_guess_duty_cycle(Waveform waveform, WaveformLabel label) {
    if (label != WaveformLabel::CUSTOM) {
        switch(label) {
            case WaveformLabel::TRIANGULAR: {
                return (waveform.get_time().value()[1] - waveform.get_time().value()[0]) / (waveform.get_time().value()[2] - waveform.get_time().value()[0]);
            }
            case WaveformLabel::UNIPOLAR_TRIANGULAR: {
                return (waveform.get_time().value()[1] - waveform.get_time().value()[0]) / (waveform.get_time().value()[3] - waveform.get_time().value()[0]);
            }
            case WaveformLabel::RECTANGULAR: {
                return (waveform.get_time().value()[2] - waveform.get_time().value()[0]) / (waveform.get_time().value()[4] - waveform.get_time().value()[0]);
            }
            case WaveformLabel::UNIPOLAR_RECTANGULAR: {
                return (waveform.get_time().value()[2] - waveform.get_time().value()[0]) / (waveform.get_time().value()[4] - waveform.get_time().value()[0]);
            }
            case WaveformLabel::BIPOLAR_RECTANGULAR: {
                return (waveform.get_time().value()[3] - waveform.get_time().value()[2]) / (waveform.get_time().value()[9] - waveform.get_time().value()[0]);
            }
            case WaveformLabel::BIPOLAR_TRIANGULAR: {
                return (waveform.get_time().value()[2] - waveform.get_time().value()[1]) / (waveform.get_time().value()[5] - waveform.get_time().value()[0]);
            }
            case WaveformLabel::FLYBACK_PRIMARY:{
                return (waveform.get_time().value()[2] - waveform.get_time().value()[0]) / (waveform.get_time().value()[4] - waveform.get_time().value()[0]);
            }
            case WaveformLabel::FLYBACK_SECONDARY:{
                return (waveform.get_time().value()[2] - waveform.get_time().value()[0]) / (waveform.get_time().value()[4] - waveform.get_time().value()[0]);
            }
            case WaveformLabel::SINUSOIDAL: {
                return 0.5;
            }
            default:
                break;
            }
    }

    Waveform sampledWaveform;
    if (!is_waveform_sampled(waveform)) {
        sampledWaveform = InputsWrapper::calculate_sampled_waveform(waveform, 0);
    }
    else {
        sampledWaveform = waveform;
    }
 
    std::vector<double> data = sampledWaveform.get_data();
    std::vector<double> diff_data;
    std::vector<double> diff_diff_data;

    for (size_t i = 0; i < data.size() - 1; ++i) {
        diff_data.push_back(roundFloat(data[i + 1] - data[i], 9));
    }
    for (size_t i = 0; i < diff_data.size() - 1; ++i) {
        diff_diff_data.push_back(fabs(roundFloat(diff_data[i + 1] - diff_data[i], 9)));
    }
    for (size_t i = 0; i < diff_data.size() - 1; ++i) {
    }

    auto settings = OpenMagnetics::Settings::GetInstance();
    double maximum = *max_element(diff_diff_data.begin(), diff_diff_data.end());
    size_t maximum_index = 0;
    size_t distanceToMiddle = settings->get_inputs_number_points_sampled_waveforms();
    for (size_t i = 0; i < diff_diff_data.size(); ++i)
    {
        if (diff_diff_data[i] == maximum) {
            if (fabs(double(settings->get_inputs_number_points_sampled_waveforms()) / 2 - i) < distanceToMiddle) {
                distanceToMiddle = fabs(double(settings->get_inputs_number_points_sampled_waveforms()) / 2 - i);
                maximum_index = i;
            }
        }
    }
    auto dutyCycle = roundFloat((maximum_index + 1.0) / settings->get_inputs_number_points_sampled_waveforms(), 2);
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
SignalDescriptor InputsWrapper::standarize_waveform(SignalDescriptor signal, double frequency) {
    SignalDescriptor standardized_signal(signal);
    if (!signal.get_waveform()) {
        if (!signal.get_processed() && !signal.get_harmonics()) {
            throw std::runtime_error("Signal is not processed");
        }
        if (signal.get_processed()) {
            auto waveform = create_waveform(signal.get_processed().value(), frequency);
            standardized_signal.set_waveform(waveform);
        }
        else {
            auto waveform = reconstruct_signal(signal.get_harmonics().value(), frequency);
            standardized_signal.set_waveform(waveform);
        }
    }

    if (standardized_signal.get_waveform() && !standardized_signal.get_waveform()->get_time()) {
        auto time = linear_spaced_array(0, 1. / roundFloat(frequency, 9), standardized_signal.get_waveform()->get_data().size());
        auto waveform = standardized_signal.get_waveform().value();
        waveform.set_time(time);
        standardized_signal.set_waveform(waveform);
    }

    return standardized_signal;
}

Waveform InputsWrapper::reconstruct_signal(Harmonics harmonics, double frequency) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    std::vector<double> data = std::vector<double>(settings->get_inputs_number_points_sampled_waveforms(), 0);
    for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_frequencies().size(); ++harmonicIndex) {
        auto amplitude = harmonics.get_amplitudes()[harmonicIndex];
        auto frequencyMultiplier = harmonics.get_frequencies()[harmonicIndex] / frequency;
        auto totalAngle = 2 * std::numbers::pi / (settings->get_inputs_number_points_sampled_waveforms() - 1) * frequencyMultiplier;
        for (size_t i = 0; i < settings->get_inputs_number_points_sampled_waveforms(); ++i) {
            double angle = i * totalAngle;
            data[i] += sin(angle) * amplitude;
        }
    }
    std::vector<double> time = linear_spaced_array(0, 1. / roundFloat(frequency, 9), settings->get_inputs_number_points_sampled_waveforms());
    Waveform waveform;
    waveform.set_data(data);
    waveform.set_time(time);
    return waveform;
}

Waveform InputsWrapper::create_waveform(Processed processed, double frequency) {
    Waveform waveform;
    double period = 1 / frequency;
    if (!processed.get_peak_to_peak()) {
        throw std::runtime_error("Signal is missing peak to peak");
    }

    auto peakToPeak = processed.get_peak_to_peak().value();
    auto offset = processed.get_offset();
    double dutyCycle = 0.5;
    if (processed.get_duty_cycle()) {
        dutyCycle = processed.get_duty_cycle().value();
    }
    std::vector<double> data;
    std::vector<double> time;

    switch (processed.get_label()) {
        case WaveformLabel::TRIANGULAR: {
            double max = peakToPeak / 2 + offset;
            double min = -peakToPeak / 2 + offset;
            double dc = dutyCycle * period;
            data = {min, max, min};
            time = {0, dc, period};
            break;
        }
        case WaveformLabel::UNIPOLAR_TRIANGULAR: {
            double max = peakToPeak + offset;
            double min = offset;
            double dc = dutyCycle * period;
            data = {min, max, min, min};
            time = {0, dc, dc, period};
            break;
        }
        case WaveformLabel::RECTANGULAR: {
            double max = peakToPeak * (1 - dutyCycle);
            double min = -peakToPeak * dutyCycle;
            double dc = dutyCycle * period;
            data = {min, max, max, min, min};
            time = {0, 0, dc, dc, period};
            break;
        }
        case WaveformLabel::UNIPOLAR_RECTANGULAR: {
            double max = peakToPeak + offset;
            double min = offset;
            double dc = std::min(0.5, dutyCycle) * period;
            data = {min, max, max, min, min};
            time = {0, 0, dc, dc, period};
            break;
        }
        case WaveformLabel::BIPOLAR_RECTANGULAR: {
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
                    period};
            break;
        }
        case WaveformLabel::BIPOLAR_TRIANGULAR: {
            double max = +peakToPeak / 2;
            double min = -peakToPeak / 2;
            double dc = std::min(0.5, dutyCycle) * period;
            data = {min, min, max, max, min, min};
            time = {0,
                    0.25 * period - dc / 2,
                    0.25 * period + dc / 2,
                    0.75 * period - dc / 2,
                    0.75 * period + dc / 2,
                    period};
            break;
        }
        case WaveformLabel::FLYBACK_PRIMARY:{
            double max = peakToPeak + offset;
            double min = offset;
            double dc = dutyCycle * period;
            data = {0, min, max, 0, 0};
            time = {0, 0, dc, dc, period};
            break;
        }
        case WaveformLabel::FLYBACK_SECONDARY:{
            double max = peakToPeak + offset;
            double min = offset;
            double dc = dutyCycle * period;
            data = {0, 0, max, min, 0};
            time = {0, dc, dc, period, period};
            break;
        }
        case WaveformLabel::SINUSOIDAL: {
            auto settings = OpenMagnetics::Settings::GetInstance();
            for (size_t i = 0; i < settings->get_inputs_number_points_sampled_waveforms(); ++i) {
                double angle = i * 2 * std::numbers::pi / (settings->get_inputs_number_points_sampled_waveforms() - 1);
                time.push_back(i * period / (settings->get_inputs_number_points_sampled_waveforms() - 1));
                data.push_back((sin(angle) * peakToPeak / 2) + offset);
            }
            break;
        }
        default:
            break;
    }

    waveform.set_data(data);
    waveform.set_time(time);

    return waveform;
}

bool InputsWrapper::is_waveform_sampled(Waveform waveform) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    if (!waveform.get_time()) {
        return false;
    }
    else {
        return waveform.get_data().size() == settings->get_inputs_number_points_sampled_waveforms();
    }

}

bool InputsWrapper::is_waveform_imported(Waveform waveform) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    if (!waveform.get_time()) {
        return false;
    }
    else {
        return waveform.get_data().size() > settings->get_inputs_number_points_sampled_waveforms();
    }

}

Waveform InputsWrapper::calculate_sampled_waveform(Waveform waveform, double frequency) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    std::vector<double> time;
    auto data = waveform.get_data();

    if (!waveform.get_time()) { // This means the waveform is equidistant
        time = linear_spaced_array(0, 1. / roundFloat(frequency, 9), data.size());
    }
    else {
        time = waveform.get_time().value();
        // Check if we need to guess the frequency from the time vector
        double period = (time.back() - time.front());
        if (frequency == 0) {
            frequency = 1.0 / period;
        }
        else {
            if (fabs((1.0 / period) - frequency) / frequency > 0.01)
                throw std::invalid_argument("Frequency: " + std::to_string(frequency) + " is not matching waveform time info with calculated frequency of: " + std::to_string(1.0 / period));
        }
    }

    size_t numberPoints = settings->get_inputs_number_points_sampled_waveforms();

    if (data.size() > numberPoints) {
        if (is_size_power_of_2(data)) {
            numberPoints = data.size();
        }
        else {
            numberPoints = round_up_size_to_power_of_2(data);
        }
    }

    auto sampledTime = linear_spaced_array(0, 1. / roundFloat(frequency, 9), numberPoints);

    std::vector<double> sampledData;

    for (size_t i = 0; i < numberPoints; i++) {
        bool found = false;
        for (size_t interpIndex = 0; interpIndex < data.size() - 1; interpIndex++) {
            if (time[interpIndex] <= sampledTime[i] && sampledTime[i] <= time[interpIndex + 1]) {
                if (time[interpIndex + 1] == time[interpIndex]) {
                    sampledData.push_back(data[interpIndex]);
                }
                else {
                    double proportion = (sampledTime[i] - time[interpIndex]) / (time[interpIndex + 1] - time[interpIndex]);
                    double interpPoint = std::lerp(data[interpIndex], data[interpIndex + 1], proportion);
                    sampledData.push_back(interpPoint);
                }
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::invalid_argument("Error while sampling waveform in point: " + std::to_string(i));
        }
    }

    if (sampledData.size() != numberPoints) {
        throw std::invalid_argument("Wrong number of sampled points");
    }

    Waveform sampledWaveform;
    sampledWaveform.set_data(sampledData);
    sampledWaveform.set_time(sampledTime);
    return sampledWaveform;
}

SignalDescriptor InputsWrapper::calculate_induced_voltage(OperatingPointExcitation& excitation,
                                                            double magnetizingInductance) {
    if (!excitation.get_current()->get_waveform()) {
        throw std::runtime_error("Current waveform is missing");
    }
    Waveform sourceWaveform = excitation.get_current()->get_waveform().value();
    std::vector<double> source = sourceWaveform.get_data();
    bool isWaveformSampled = is_waveform_sampled(sourceWaveform);
    bool isWaveformImported = is_waveform_imported(sourceWaveform);
    std::vector<double> time = sourceWaveform.get_time().value();
    std::vector<double> derivative;
    std::vector<double> derivativeTime;
    std::vector<double> voltageData;
    Waveform voltageWaveform;
    SignalDescriptor voltageSignalDescriptor;
    Waveform resultWaveform(sourceWaveform);
    auto originalTime = time;

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

    if (isWaveformSampled || isWaveformImported) {
        for (size_t i = 0; i < derivative.size(); ++i) {
            if (derivativeTime[i] == 0) {
                if (i > 0) {
                    voltageData.push_back(voltageData[i - 1]);
                }
                else {
                    voltageData.push_back(0);
                }
            }
            else { 
                voltageData.push_back(magnetizingInductance * derivative[i] / derivativeTime[i]);
            }
        }
        voltageWaveform.set_time(originalTime);
    }
    else {
        std::vector<double> finalTime;
        for (size_t i = 0; i < derivative.size() - 1; ++i) {
            if (derivativeTime[i] == 0) {
                continue;
            }
            voltageData.push_back(magnetizingInductance * derivative[i] / derivativeTime[i]);
            voltageData.push_back(magnetizingInductance * derivative[i] / derivativeTime[i]);
            finalTime.push_back(originalTime[i]);
            finalTime.push_back(originalTime[i + 1]);
        }
        if (derivativeTime[derivative.size() - 1] != 0) {
            finalTime.push_back(originalTime[derivative.size() - 1]);
            voltageData.push_back(magnetizingInductance * derivative[derivative.size() - 1] / derivativeTime[derivative.size() - 1]);
        }
        voltageWaveform.set_time(finalTime);
    }

    voltageWaveform.set_data(voltageData);
    voltageSignalDescriptor.set_waveform(voltageWaveform);
    auto sampledWaveform = InputsWrapper::calculate_sampled_waveform(voltageWaveform, excitation.get_frequency());
    voltageSignalDescriptor.set_harmonics(calculate_harmonics_data(sampledWaveform, excitation.get_frequency()));
    voltageSignalDescriptor.set_processed(
        calculate_processed_data(voltageSignalDescriptor, sampledWaveform, true));

    resultWaveform.set_data(derivative);
    return voltageSignalDescriptor;
}

Waveform InputsWrapper::calculate_derivative_waveform(Waveform waveform) {
    std::vector<double> sourceData = waveform.get_data();
    std::vector<double> sourceTime = waveform.get_time().value();

    std::vector<double> tempData = sourceData;
    std::vector<double> tempTime = sourceTime;


    std::vector<double> derivative;
    std::vector<double> derivativeTime;
    std::vector<double> data;
    Waveform derivativeWaveform;
    SignalDescriptor voltageSignalDescriptor;
    auto originalTime = tempTime;

    if (is_waveform_sampled(waveform)) {
        tempData.push_back(tempData[0]);
        double difference = tempTime[tempTime.size() - 1] - tempTime[tempTime.size() - 2];
        tempTime.push_back(tempTime[tempTime.size() - 1] + difference);
    }
    else {
        tempData.push_back(tempData[1]);
        tempTime.push_back(tempTime[tempTime.size() - 1] + tempTime[1]);
    }

    std::adjacent_difference(tempData.begin(), tempData.end(), tempData.begin());
    derivative = std::vector<double>(tempData.begin() + 1, tempData.end());

    std::adjacent_difference(tempTime.begin(), tempTime.end(), tempTime.begin());
    derivativeTime = std::vector<double>(tempTime.begin() + 1, tempTime.end());

    if (is_waveform_sampled(waveform)) {
        for (size_t i = 0; i < derivative.size(); ++i) {
            data.push_back(derivative[i] / derivativeTime[i]);
        }
        derivativeWaveform.set_time(originalTime);
    }
    else {
        std::vector<double> finalTime;
        for (size_t i = 0; i < derivative.size() - 1; ++i) {
            if (derivativeTime[i] == 0) {
                continue;
            }
            data.push_back(derivative[i] / derivativeTime[i]);
            data.push_back(derivative[i] / derivativeTime[i]);
            finalTime.push_back(originalTime[i]);
            finalTime.push_back(originalTime[i + 1]);
        }
        if (derivativeTime[derivative.size() - 1] != 0) {
            finalTime.insert(finalTime.begin(), finalTime[0]);
            data.insert(data.begin(), data[data.size() - 1]);
        }
        derivativeWaveform.set_time(finalTime);
    }

    derivativeWaveform.set_data(data);
    return derivativeWaveform;
}


Waveform InputsWrapper::calculate_integral_waveform(Waveform waveform) {
    std::vector<double> data = waveform.get_data();
    std::vector<double> time = waveform.get_time().value();
    std::vector<double> integration;
    Waveform resultWaveform(waveform);

    double integral = 0;
    integration.push_back(integral);
    for (size_t i = 0; i < time.size() - 1; ++i) {
        double timePerPoint = time[i + 1] - time[i];
        integral += data[i] * timePerPoint;
        integration.push_back(integral);
    }
    resultWaveform.set_data(integration);

    double integrationAverage = calculate_waveform_average(resultWaveform);
    resultWaveform = sum_waveform(resultWaveform, -integrationAverage);

    std::vector<double> distinctData;
    std::vector<double> distinctTime;
    for (size_t i = 0; i < resultWaveform.get_data().size(); ++i){
        if (distinctData.size() != 0) {
            if (resultWaveform.get_data()[i] == distinctData.back() && resultWaveform.get_time().value()[i] == distinctTime.back())
                continue;
        }

        distinctData.push_back(resultWaveform.get_data()[i]);
        distinctTime.push_back(resultWaveform.get_time().value()[i]);
    }

    resultWaveform.set_data(distinctData);
    resultWaveform.set_time(distinctTime);

    return resultWaveform;
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
    auto sampledWaveform = InputsWrapper::calculate_sampled_waveform(waveform, frequency);
    signalDescriptor.set_harmonics(calculate_harmonics_data(sampledWaveform, frequency));
    signalDescriptor.set_processed(calculate_processed_data(signalDescriptor, sampledWaveform, true, signalDescriptor.get_processed()));
    return signalDescriptor;
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
    // waveform.set_time(primarySignalDescriptor.get_waveform()->get_time());
    reflected_waveform.set_waveform(waveform);

    return reflected_waveform;
}

SignalDescriptor InputsWrapper::reflect_waveform(SignalDescriptor signal,
                                                 double ratio,
                                                 WaveformLabel label) {
    
    if (label == WaveformLabel::CUSTOM) {
        return reflect_waveform(signal, ratio);
    }
    Processed processed;
    if (!signal.get_processed()) {
        auto waveform = signal.get_waveform().value();
        if (is_waveform_sampled(waveform)) {
            waveform = compress_waveform(waveform);
        }
        processed = calculate_basic_processed_data(waveform);
    }
    else {
        processed = signal.get_processed().value();
    }
    Waveform newWaveform;

    double period = signal.get_waveform()->get_time()->back() - signal.get_waveform()->get_time()->front();
    double frequency = 1.0 / period;
    double peakToPeak = processed.get_peak_to_peak().value() * ratio;
    double offset = processed.get_offset() * ratio;
    double dutyCycle = processed.get_duty_cycle().value();


    switch(label) {
        case WaveformLabel::FLYBACK_PRIMARY:
            processed.set_label(WaveformLabel::FLYBACK_SECONDARY);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            newWaveform = create_waveform(processed, frequency);
            break;
        case WaveformLabel::FLYBACK_SECONDARY:
            processed.set_label(WaveformLabel::FLYBACK_PRIMARY);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            newWaveform = create_waveform(processed, frequency);
            break;
        case WaveformLabel::UNIPOLAR_TRIANGULAR: {
            double max = peakToPeak * dutyCycle / (1 - dutyCycle) + offset;
            double min = offset;
            double dc = dutyCycle * period;
            std::vector<double> data = {min, min, max, min};
            std::vector<double> time = {0, dc, dc, period};
            newWaveform.set_data(data);
            newWaveform.set_time(time);
            break;
        }
        case WaveformLabel::UNIPOLAR_RECTANGULAR: {
            double max = peakToPeak * dutyCycle / (1 - dutyCycle) + offset;
            double min = offset;
            double dc = dutyCycle * period;
            std::vector<double> data = {-min, -min, -max, -max, -min};
            std::vector<double> time = {0, dc, dc, period, period};
            newWaveform.set_data(data);
            newWaveform.set_time(time);
            break;
        }
        default:
            return reflect_waveform(signal, ratio);
            break;
    }

    SignalDescriptor newSignal;
    newSignal.set_waveform(newWaveform);
    return newSignal;
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
                auto voltageWaveform = excitation.get_voltage()->get_waveform().value();
                auto sampledWaveform = calculate_sampled_waveform(voltageWaveform, excitation.get_frequency());
                excitation.set_current(
                    calculate_magnetizing_current(excitation, sampledWaveform, magnetizingInductance));
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

Processed InputsWrapper::calculate_processed_data(SignalDescriptor excitation,
                                            Waveform sampledWaveform,
                                            bool includeAdvancedData,
                                            std::optional<Processed> processed) {

    auto harmonics = excitation.get_harmonics().value();
    return InputsWrapper::calculate_processed_data(harmonics, sampledWaveform, includeAdvancedData, processed);
}

double calculate_offset(Waveform waveform, WaveformLabel label) {
    switch (label) {
        case WaveformLabel::TRIANGULAR:
            return (waveform.get_data()[0] + waveform.get_data()[1]) / 2;
        case WaveformLabel::UNIPOLAR_TRIANGULAR:
            return *min_element(waveform.get_data().begin(), waveform.get_data().end());
        case WaveformLabel::RECTANGULAR:
            return 0;
        case WaveformLabel::UNIPOLAR_RECTANGULAR:
            return *min_element(waveform.get_data().begin(), waveform.get_data().end());
        case WaveformLabel::BIPOLAR_RECTANGULAR:
            return 0;
        case WaveformLabel::BIPOLAR_TRIANGULAR:
            return 0;
        case WaveformLabel::SINUSOIDAL:
            return (*max_element(waveform.get_data().begin(), waveform.get_data().end()) + *min_element(waveform.get_data().begin(), waveform.get_data().end())) / 2;
        case WaveformLabel::CUSTOM:
            return (*max_element(waveform.get_data().begin(), waveform.get_data().end()) + *min_element(waveform.get_data().begin(), waveform.get_data().end())) / 2;
        case WaveformLabel::FLYBACK_PRIMARY:
            return waveform.get_data()[1];
        case WaveformLabel::FLYBACK_SECONDARY:
            return waveform.get_data()[3];
    }
    return 0;
}

Processed InputsWrapper::calculate_basic_processed_data(Waveform waveform) {
    Processed processed;
    std::vector<double> dataToProcess;
    auto sampledWaveform = waveform;
    auto compressedWaveform = waveform;

    for (size_t i = 0; i < waveform.get_data().size(); ++i) {
        if (std::isnan(waveform.get_data()[i])) {
            throw std::invalid_argument("Waveform data contains NaN");
        }
    }

    if (!is_waveform_sampled(waveform)) {
        sampledWaveform = calculate_sampled_waveform(waveform);
    }
    else {
        compressedWaveform = compress_waveform(waveform);
    }

    WaveformLabel label;

    label = try_guess_waveform_label(compressedWaveform);
    processed.set_label(label);

    if (is_waveform_sampled(waveform)) {
        processed.set_average(std::accumulate(std::begin(sampledWaveform.get_data()), std::end(sampledWaveform.get_data()), 0.0) /
                         sampledWaveform.get_data().size());
    }
    else {
        auto average = calculate_waveform_average(compressedWaveform);
        processed.set_average(average);
    }

    double offset = calculate_offset(compressedWaveform, label);
    processed.set_offset(offset);

    processed.set_peak_to_peak(*max_element(compressedWaveform.get_data().begin(), compressedWaveform.get_data().end()) -
                               *min_element(compressedWaveform.get_data().begin(), compressedWaveform.get_data().end()));

    if (label == WaveformLabel::FLYBACK_PRIMARY ||
        label == WaveformLabel::FLYBACK_SECONDARY ||
        label == WaveformLabel::UNIPOLAR_TRIANGULAR ||
        label == WaveformLabel::UNIPOLAR_RECTANGULAR) {
        processed.set_peak_to_peak(processed.get_peak_to_peak().value() - offset);
    }

    processed.set_peak(*max_element(compressedWaveform.get_data().begin(), compressedWaveform.get_data().end()));

    processed.set_duty_cycle(try_guess_duty_cycle(compressedWaveform, label));

    return processed;
}

Processed InputsWrapper::calculate_processed_data(Harmonics harmonics,
                                            Waveform waveform,
                                            bool includeAdvancedData,
                                            std::optional<Processed> processed) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto sampledDataToProcess = waveform;

    if (waveform.get_time() && waveform.get_data().size() < settings->get_inputs_number_points_sampled_waveforms()) {
        auto frequency = harmonics.get_frequencies()[1];
        sampledDataToProcess = calculate_sampled_waveform(waveform, frequency);
    }

    Processed processedResult;
    if (processed) {
        processedResult = processed.value();
    }
    else {
        processedResult = calculate_basic_processed_data(sampledDataToProcess);
    }

    if (includeAdvancedData) {
         {
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

            processedResult.set_effective_frequency(effectiveFrequency);
        }
        {
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

            processedResult.set_ac_effective_frequency(effectiveFrequency);
        }
        {
            double rms = 0.0;
            for (int i = 0; i < int(sampledDataToProcess.get_data().size()); ++i) {
                rms += sampledDataToProcess.get_data()[i] * sampledDataToProcess.get_data()[i];
            }
            rms /= sampledDataToProcess.get_data().size();
            rms = sqrt(rms);
            processedResult.set_rms(rms);
        }
        {
            std::vector<double> dividend;
            double divisor = harmonics.get_amplitudes()[1];
            double thd = 0;
            for (size_t i = 2; i < harmonics.get_amplitudes().size(); ++i) {
                dividend.push_back(pow(harmonics.get_amplitudes()[i], 2));
            }

            if (divisor > 0)
                thd = sqrt(std::reduce(dividend.begin(), dividend.end())) / divisor;
            processedResult.set_thd(thd);
        }
    }

    return processedResult;
}

Harmonics InputsWrapper::calculate_harmonics_data(Waveform waveform, double frequency) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    bool trimHarmonics = settings->get_inputs_trim_harmonics();
    bool isWaveformImported = is_waveform_imported(waveform);
    Harmonics harmonics;

    std::vector<std::complex<double>> data;
    for (std::size_t i = 0; i < waveform.get_data().size(); ++i) {
        data.emplace_back(waveform.get_data()[i]);
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


    if (isWaveformImported && trimHarmonics) {
        auto mainHarmonicIndexes = get_main_harmonic_indexes(harmonics, Defaults().importedWaveformHarmonicAmplitudeThreshold);
        Harmonics reducedHarmonics;
        reducedHarmonics.get_mutable_amplitudes().push_back(harmonics.get_amplitudes()[0]);
        reducedHarmonics.get_mutable_amplitudes().push_back(harmonics.get_amplitudes()[1]);
        reducedHarmonics.get_mutable_frequencies().push_back(harmonics.get_frequencies()[0]);
        reducedHarmonics.get_mutable_frequencies().push_back(harmonics.get_frequencies()[1]);
        for (auto harmonicIndex : mainHarmonicIndexes) {
            if (harmonicIndex == 0 || harmonicIndex == 1) {
                continue;
            }
            reducedHarmonics.get_mutable_amplitudes().push_back(harmonics.get_amplitudes()[harmonicIndex]);
            reducedHarmonics.get_mutable_frequencies().push_back(harmonics.get_frequencies()[harmonicIndex]);
        }
        harmonics = reducedHarmonics;
    }

    return harmonics;
}

Waveform InputsWrapper::compress_waveform(Waveform waveform) {
    Waveform compressedWaveform;
    auto data = waveform.get_data();
    data.push_back(data[0]);
    auto time = waveform.get_time().value();
    time.push_back(time[time.size() - 1] + (time[time.size() - 1] - time[time.size() - 2]));
    std::vector<double> compressedData;
    std::vector<double> compressedTime;
    double previousSlope = DBL_MAX;
    for (size_t i = 0; i < data.size() - 1; ++i){
        auto slope = (data[i + 1] - data[i]) / (time[i + 1] - time[i]);
        if (fabs((slope - previousSlope) / previousSlope) > 0.01) {
            compressedData.push_back(data[i]);
            compressedTime.push_back(time[i]);
        }
        previousSlope = slope;
    }
    compressedData.push_back(data.back());
    compressedTime.push_back(time.back());
    waveform.set_data(compressedData);
    waveform.set_time(compressedTime);
    return waveform;
}

SignalDescriptor InputsWrapper::calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                                Waveform sampledWaveform,
                                                                double magnetizingInductance,
                                                                bool compress,
                                                                double offset) {
    double dcCurrent = 0;
    if (magnetizingInductance <= 0) {
        throw std::runtime_error("magnetizingInductance cannot be zero or negative");
    }
    if (excitation.get_current()) {

        if (!excitation.get_current()->get_processed()) {
            auto currentExcitation = excitation.get_current().value();
            auto currentExcitationWaveform = currentExcitation.get_waveform().value();
            auto sampledCurrentWaveform = calculate_sampled_waveform(currentExcitationWaveform, excitation.get_frequency());
            currentExcitation.set_harmonics(calculate_harmonics_data(sampledCurrentWaveform, excitation.get_frequency()));
            currentExcitation.set_processed(calculate_processed_data(currentExcitation, sampledCurrentWaveform, true, currentExcitation.get_processed()));
            excitation.set_current(currentExcitation);
        }

        dcCurrent = excitation.get_current()->get_processed()->get_offset();
    }
    else {
        dcCurrent = offset;
    }

    SignalDescriptor magnetizingCurrentExcitation;
    Waveform sampledMagnetizingCurrentWaveform;

    if (excitation.get_current() && 
        (excitation.get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_PRIMARY || 
        excitation.get_current()->get_processed()->get_label() == WaveformLabel::FLYBACK_SECONDARY || 
        excitation.get_current()->get_processed()->get_label() == WaveformLabel::UNIPOLAR_TRIANGULAR)) {

        double offset = excitation.get_current()->get_processed()->get_offset();
        double peakToPeak = excitation.get_current()->get_processed()->get_peak_to_peak().value();
        Processed triangularProcessed;
        triangularProcessed.set_label(WaveformLabel::TRIANGULAR);
        triangularProcessed.set_offset(offset + peakToPeak / 2);
        triangularProcessed.set_peak_to_peak(peakToPeak);
        auto newWaveform = create_waveform(triangularProcessed, excitation.get_frequency());
        sampledMagnetizingCurrentWaveform = calculate_sampled_waveform(newWaveform, excitation.get_frequency());            
    }
    else {
        sampledMagnetizingCurrentWaveform = calculate_integral_waveform(sampledWaveform);
        SignalDescriptor magnetizingCurrentExcitation;

        sampledMagnetizingCurrentWaveform = multiply_waveform(sampledMagnetizingCurrentWaveform, 1.0 / magnetizingInductance);
        sampledMagnetizingCurrentWaveform = sum_waveform(sampledMagnetizingCurrentWaveform, dcCurrent);
    }
 
    if (compress){
        magnetizingCurrentExcitation.set_waveform(compress_waveform(sampledMagnetizingCurrentWaveform));
    }
    else {
        magnetizingCurrentExcitation.set_waveform(sampledMagnetizingCurrentWaveform);
    }

    magnetizingCurrentExcitation.set_harmonics(
        calculate_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
    magnetizingCurrentExcitation.set_processed(
        calculate_processed_data(magnetizingCurrentExcitation, sampledMagnetizingCurrentWaveform));
    return magnetizingCurrentExcitation;
}

SignalDescriptor InputsWrapper::calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                                double magnetizingInductance,
                                                                bool compress,
                                                                double offset) {
    if (!excitation.get_voltage()) {
        throw std::invalid_argument("Missing voltage signal");
    }
    auto voltage_excitation = excitation.get_voltage().value();
    voltage_excitation = standarize_waveform(voltage_excitation, excitation.get_frequency());
    auto waveform = voltage_excitation.get_waveform().value();

    auto sampledWaveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
    return calculate_magnetizing_current(excitation, sampledWaveform, magnetizingInductance, compress, offset);
}

OperatingPoint InputsWrapper::process_operating_point(OperatingPoint operatingPoint, double magnetizingInductance) {
    std::vector<OperatingPointExcitation> processed_excitations_per_winding;
    for (auto& excitation : operatingPoint.get_mutable_excitations_per_winding()) {
        // Here we processed this excitation current
        if (excitation.get_current()) {
            auto current_excitation = excitation.get_current().value();
            current_excitation = standarize_waveform(current_excitation, excitation.get_frequency());
            auto waveform = current_excitation.get_waveform().value();
            auto sampledWaveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
            current_excitation.set_harmonics(calculate_harmonics_data(sampledWaveform, excitation.get_frequency()));
            current_excitation.set_processed(calculate_processed_data(current_excitation, sampledWaveform, true, current_excitation.get_processed()));
            excitation.set_current(current_excitation);
        }
        // Here we processed this excitation voltage
        if (excitation.get_voltage()) {
            auto voltage_excitation = excitation.get_voltage().value();
            voltage_excitation = standarize_waveform(voltage_excitation, excitation.get_frequency());
            auto waveform = voltage_excitation.get_waveform().value();
            auto sampledWaveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
            voltage_excitation.set_harmonics(calculate_harmonics_data(sampledWaveform, excitation.get_frequency()));
            voltage_excitation.set_processed(calculate_processed_data(voltage_excitation, sampledWaveform));
            excitation.set_voltage(voltage_excitation);

            if (!excitation.get_magnetizing_current() && magnetizingInductance > 0) {
                excitation.set_magnetizing_current(
                    calculate_magnetizing_current(excitation, sampledWaveform, magnetizingInductance));
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

    InputsWrapper inputs;

    DesignRequirements designRequirements;
    DimensionWithTolerance magnetizingInductanceWithTolerance;
    magnetizingInductanceWithTolerance.set_minimum(magnetizingInductance * 0.8);
    magnetizingInductanceWithTolerance.set_nominal(magnetizingInductance);
    magnetizingInductanceWithTolerance.set_maximum(magnetizingInductance * 1.2);

    OpenMagnetics::InsulationRequirements insulationRequirements;
    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    auto cti = OpenMagnetics::Cti::GROUP_I;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
    auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = OpenMagnetics::InsulationType::BASIC;

    insulationRequirements.set_altitude(altitude);
    insulationRequirements.set_cti(cti);
    insulationRequirements.set_insulation_type(insulationType);
    insulationRequirements.set_main_supply_voltage(mainSupplyVoltage);
    insulationRequirements.set_overvoltage_category(overvoltageCategory);
    insulationRequirements.set_pollution_degree(pollutionDegree);
    insulationRequirements.set_standards(standards);
    designRequirements.set_insulation(insulationRequirements);

    designRequirements.set_magnetizing_inductance(magnetizingInductanceWithTolerance);
    for (auto& turnsRatio : turnsRatios) {
        DimensionWithTolerance turnsRatiosWithTolerance;
        turnsRatiosWithTolerance.set_nominal(turnsRatio);
        designRequirements.get_mutable_turns_ratios().push_back(turnsRatiosWithTolerance);
    }
    inputs.set_design_requirements(designRequirements);

    OperatingPoint operatingPoint;
    OperatingConditions conditions;
    conditions.set_ambient_temperature(temperature);
    conditions.set_ambient_relative_humidity(std::nullopt);
    conditions.set_cooling(std::nullopt);
    conditions.set_name(std::nullopt);
    operatingPoint.set_conditions(conditions);
    {
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        SignalDescriptor voltage;
        Processed processed;
        processed.set_label(waveShape);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(0);
        voltage.set_processed(processed);
        voltage = standarize_waveform(voltage, frequency);

        excitation.set_voltage(voltage);
        if (magnetizingInductance > 0) {
            auto current = calculate_magnetizing_current(excitation, magnetizingInductance, true, dcCurrent);
            excitation.set_current(current);
            excitation.set_magnetizing_current(current);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }
    for (auto& turnsRatio : turnsRatios) {
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        SignalDescriptor voltage;
        Processed processed;
        processed.set_label(waveShape);
        processed.set_peak_to_peak(peakToPeak * turnsRatio);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(0);
        voltage.set_processed(processed);
        voltage = standarize_waveform(voltage, frequency);
        excitation.set_voltage(voltage);
        if (magnetizingInductance > 0) {
            auto current = calculate_magnetizing_current(excitation, magnetizingInductance, true, dcCurrent);
            excitation.set_current(current);
            excitation.set_magnetizing_current(current);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    inputs.get_mutable_operating_points().push_back(operatingPoint);
    inputs.process_waveforms();
    return inputs;

}
InputsWrapper InputsWrapper::create_quick_operating_point_only_current(double frequency,
                                                          double magnetizingInductance,
                                                          double temperature,
                                                          WaveformLabel waveShape,
                                                          double peakToPeak,
                                                          double dutyCycle,
                                                          double dcCurrent,
                                                          std::vector<double> turnsRatios) {

    InputsWrapper inputs;

    DesignRequirements designRequirements;
    DimensionWithTolerance magnetizingInductanceWithTolerance;
    magnetizingInductanceWithTolerance.set_minimum(magnetizingInductance * 0.8);
    magnetizingInductanceWithTolerance.set_nominal(magnetizingInductance);
    magnetizingInductanceWithTolerance.set_maximum(magnetizingInductance * 1.2);
    OpenMagnetics::InsulationRequirements insulationRequirements;

    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    auto cti = OpenMagnetics::Cti::GROUP_I;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
    auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = OpenMagnetics::InsulationType::BASIC;

    insulationRequirements.set_altitude(altitude);
    insulationRequirements.set_cti(cti);
    insulationRequirements.set_insulation_type(insulationType);
    insulationRequirements.set_main_supply_voltage(mainSupplyVoltage);
    insulationRequirements.set_overvoltage_category(overvoltageCategory);
    insulationRequirements.set_pollution_degree(pollutionDegree);
    insulationRequirements.set_standards(standards);
    designRequirements.set_insulation(insulationRequirements);
    designRequirements.set_magnetizing_inductance(magnetizingInductanceWithTolerance);
    for (auto& turnsRatio : turnsRatios) {
        DimensionWithTolerance turnsRatiosWithTolerance;
        turnsRatiosWithTolerance.set_nominal(turnsRatio);
        designRequirements.get_mutable_turns_ratios().push_back(turnsRatiosWithTolerance);
    }
    designRequirements.set_wiring_technology(WiringTechnology::WOUND);
    inputs.set_design_requirements(designRequirements);

    OperatingPoint operatingPoint;
    OperatingConditions conditions;
    conditions.set_ambient_temperature(temperature);
    conditions.set_ambient_relative_humidity(std::nullopt);
    conditions.set_cooling(std::nullopt);
    conditions.set_name(std::nullopt);
    operatingPoint.set_conditions(conditions);
    {
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        SignalDescriptor current;
        Processed processed;
        processed.set_label(waveShape);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(dcCurrent);
        current.set_processed(processed);
        current = standarize_waveform(current, frequency);
        excitation.set_current(current);
        if (magnetizingInductance > 0) {
            auto voltage = calculate_induced_voltage(excitation, magnetizingInductance);
            excitation.set_voltage(voltage);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }
    for (auto& turnsRatio : turnsRatios) {
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        SignalDescriptor current;
        Processed processed;
        processed.set_label(waveShape);
        processed.set_peak_to_peak(peakToPeak * turnsRatio);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(dcCurrent);
        current.set_processed(processed);
        current = standarize_waveform(current, frequency);
        excitation.set_current(current);
        if (magnetizingInductance > 0) {
            auto voltage = calculate_induced_voltage(excitation, magnetizingInductance / pow(turnsRatio, 2));
            excitation.set_voltage(voltage);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    inputs.get_mutable_operating_points().push_back(operatingPoint);
    inputs.process_waveforms();

    return inputs;
}
InputsWrapper InputsWrapper::create_quick_operating_point_only_current(double frequency,
                                                          double magnetizingInductance,
                                                          double temperature,
                                                          WaveformLabel waveShape,
                                                          std::vector<double> peakToPeaks,
                                                          double dutyCycle,
                                                          double dcCurrent,
                                                          std::vector<double> turnsRatios) {

    InputsWrapper inputs;

    DesignRequirements designRequirements;
    DimensionWithTolerance magnetizingInductanceWithTolerance;
    magnetizingInductanceWithTolerance.set_minimum(magnetizingInductance * 0.8);
    magnetizingInductanceWithTolerance.set_nominal(magnetizingInductance);
    magnetizingInductanceWithTolerance.set_maximum(magnetizingInductance * 1.2);
    OpenMagnetics::InsulationRequirements insulationRequirements;

    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    auto cti = OpenMagnetics::Cti::GROUP_I;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
    auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = OpenMagnetics::InsulationType::BASIC;

    insulationRequirements.set_altitude(altitude);
    insulationRequirements.set_cti(cti);
    insulationRequirements.set_insulation_type(insulationType);
    insulationRequirements.set_main_supply_voltage(mainSupplyVoltage);
    insulationRequirements.set_overvoltage_category(overvoltageCategory);
    insulationRequirements.set_pollution_degree(pollutionDegree);
    insulationRequirements.set_standards(standards);
    designRequirements.set_insulation(insulationRequirements);
    designRequirements.set_magnetizing_inductance(magnetizingInductanceWithTolerance);
    for (auto& turnsRatio : turnsRatios) {
        DimensionWithTolerance turnsRatiosWithTolerance;
        turnsRatiosWithTolerance.set_nominal(turnsRatio);
        designRequirements.get_mutable_turns_ratios().push_back(turnsRatiosWithTolerance);
    }
    designRequirements.set_wiring_technology(WiringTechnology::WOUND);
    inputs.set_design_requirements(designRequirements);

    OperatingPoint operatingPoint;
    OperatingConditions conditions;
    conditions.set_ambient_temperature(temperature);
    conditions.set_ambient_relative_humidity(std::nullopt);
    conditions.set_cooling(std::nullopt);
    conditions.set_name(std::nullopt);
    operatingPoint.set_conditions(conditions);
    for (size_t i = 0; i < peakToPeaks.size(); i++) {
        double peakToPeak = peakToPeaks[i];
        double turnsRatio = (i==0) ? 1 : turnsRatios[i - 1];
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        SignalDescriptor current;
        Processed processed;
        processed.set_label(waveShape);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(dcCurrent);
        current.set_processed(processed);
        current = standarize_waveform(current, frequency);
        excitation.set_current(current);
        if (magnetizingInductance > 0) {
            auto voltage = calculate_induced_voltage(excitation, magnetizingInductance / pow(turnsRatio, 2));
            excitation.set_voltage(voltage);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    inputs.get_mutable_operating_points().push_back(operatingPoint);
    inputs.process_waveforms();

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

            auto currentSampledWaveform = InputsWrapper::calculate_sampled_waveform(currentWaveform, frequency);
            current.set_waveform(currentSampledWaveform);
            current.set_harmonics(calculate_harmonics_data(currentSampledWaveform, frequency));
            current.set_processed(calculate_processed_data(current, currentSampledWaveform, true, current.get_processed()));
            operatingPoint->get_mutable_excitations_per_winding()[0].set_current(current);


            currentWaveform = current.get_waveform().value();

            currentWaveform = operatingPoint->get_excitations_per_winding()[0].get_current()->get_waveform().value();
        }
    }
    if (excitation.get_voltage()) {
        auto voltage = operatingPoint->get_excitations_per_winding()[0].get_voltage().value();
        auto voltageWaveform = voltage.get_waveform().value();
        if (!is_size_power_of_2(voltageWaveform.get_data())) {
            auto voltageSampledWaveform = InputsWrapper::calculate_sampled_waveform(voltageWaveform, frequency);
            voltage.set_waveform(voltageSampledWaveform);
            operatingPoint->get_mutable_excitations_per_winding()[0].set_voltage(voltage);
        }
    }
}

double InputsWrapper::calculate_waveform_coefficient(OperatingPoint* operatingPoint) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(*operatingPoint);
    double frequency = excitation.get_frequency();
    Waveform sampledWaveform = excitation.get_voltage()->get_waveform().value();

    if (sampledWaveform.get_time() && sampledWaveform.get_data().size() < settings->get_inputs_number_points_sampled_waveforms()) {
        sampledWaveform = calculate_sampled_waveform(sampledWaveform, frequency);
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

    double voltageRms = excitation.get_voltage()->get_processed()->get_rms().value();

    double waveformCoefficient = 2 * voltageRms / (frequency * integral);
    return waveformCoefficient;
}

double InputsWrapper::calculate_instantaneous_power(OperatingPointExcitation excitation) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double frequency = excitation.get_frequency();
    if (!excitation.get_voltage()) {
        throw std::runtime_error("Voltage signal is missing");
    }
    if (!excitation.get_voltage()->get_waveform()) {
        throw std::runtime_error("Voltage waveform is missing");
    }
    if (!excitation.get_current()) {
        throw std::runtime_error("Current signal is missing");
    }
    if (!excitation.get_current()->get_waveform()) {
        throw std::runtime_error("Current waveform is missing");
    }
    Waveform voltageSampledWaveform = excitation.get_voltage()->get_waveform().value();
    Waveform currentSampledWaveform = excitation.get_current()->get_waveform().value();

    if (voltageSampledWaveform.get_time() && voltageSampledWaveform.get_data().size() != settings->get_inputs_number_points_sampled_waveforms()) {
        voltageSampledWaveform = calculate_sampled_waveform(voltageSampledWaveform, frequency);
    }

    if (currentSampledWaveform.get_time() && currentSampledWaveform.get_data().size() != settings->get_inputs_number_points_sampled_waveforms()) {
        currentSampledWaveform = calculate_sampled_waveform(currentSampledWaveform, frequency);
    }

    std::vector<double> powerPoints;

    for (size_t i = 0; i < settings->get_inputs_number_points_sampled_waveforms(); i++) {
        powerPoints.push_back(fabs(voltageSampledWaveform.get_data()[i] * currentSampledWaveform.get_data()[i]));
    }

    double instantaneousPower = std::accumulate(std::begin(powerPoints), std::end(powerPoints), 0.0) /
                             powerPoints.size();
    return instantaneousPower;
}

bool is_closed_enough(double x, double y, double error){
    return fabs(x - y) <= error;
}

WaveformLabel InputsWrapper::try_guess_waveform_label(Waveform waveform) {;
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto compressedWaveform = waveform;
    if (is_waveform_sampled(waveform))
        compressedWaveform = compress_waveform(waveform);
    double period = compressedWaveform.get_time()->back() - compressedWaveform.get_time()->front();

    if (compressedWaveform.get_data().size() == 3 && 
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[2]) {
            return WaveformLabel::TRIANGULAR;
    }
    else if (compressedWaveform.get_data().size() == 4 &&
        is_closed_enough(compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[2], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[2] == compressedWaveform.get_data()[3] &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[3]) {
            return WaveformLabel::UNIPOLAR_TRIANGULAR;
    }
    else if (compressedWaveform.get_data().size() == 5 &&
        !is_closed_enough((compressedWaveform.get_time().value()[2] - compressedWaveform.get_time().value()[0]) * compressedWaveform.get_data()[2] + (compressedWaveform.get_time().value()[4] - compressedWaveform.get_time().value()[2]) * compressedWaveform.get_data()[4], 0 , period) &&
        is_closed_enough(compressedWaveform.get_time().value()[0], compressedWaveform.get_time().value()[1], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[1] == compressedWaveform.get_data()[2] &&
        is_closed_enough(compressedWaveform.get_time().value()[2], compressedWaveform.get_time().value()[3], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[3] == compressedWaveform.get_data()[4] &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
            return WaveformLabel::UNIPOLAR_RECTANGULAR;
    }
    else if (compressedWaveform.get_data().size() == 5 &&
        is_closed_enough((compressedWaveform.get_time().value()[2] - compressedWaveform.get_time().value()[0]) * compressedWaveform.get_data()[2] + (compressedWaveform.get_time().value()[4] - compressedWaveform.get_time().value()[2]) * compressedWaveform.get_data()[4], 0 , period) &&
        is_closed_enough(compressedWaveform.get_time().value()[0], compressedWaveform.get_time().value()[1], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[1] == compressedWaveform.get_data()[2] &&
        is_closed_enough(compressedWaveform.get_time().value()[2], compressedWaveform.get_time().value()[3], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[3] == compressedWaveform.get_data()[4] &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
            return WaveformLabel::RECTANGULAR;
    }
    else if (compressedWaveform.get_data().size() == 5 &&
        is_closed_enough((compressedWaveform.get_time().value()[1] - compressedWaveform.get_time().value()[0]) * compressedWaveform.get_data()[1] + (compressedWaveform.get_time().value()[3] - compressedWaveform.get_time().value()[2]) * compressedWaveform.get_data()[3], 0 , period) &&
        is_closed_enough(compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[2], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[1] &&
        is_closed_enough(compressedWaveform.get_time().value()[3], compressedWaveform.get_time().value()[4], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[2] == compressedWaveform.get_data()[3] &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
            return WaveformLabel::RECTANGULAR;
    }
    else if (compressedWaveform.get_data().size() == 10 &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[1] &&
        is_closed_enough(compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[2], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[2] == compressedWaveform.get_data()[3] &&
        is_closed_enough(compressedWaveform.get_time().value()[3], compressedWaveform.get_time().value()[4], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[4] == compressedWaveform.get_data()[5] &&
        is_closed_enough(compressedWaveform.get_time().value()[5], compressedWaveform.get_time().value()[6], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[6] == compressedWaveform.get_data()[7] &&
        is_closed_enough(compressedWaveform.get_time().value()[7], compressedWaveform.get_time().value()[8], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[8] == compressedWaveform.get_data()[9] &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[9]) {
            return WaveformLabel::BIPOLAR_RECTANGULAR;
    }
    else if (compressedWaveform.get_data().size() == 6 &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[1] &&
        is_closed_enough(compressedWaveform.get_time().value()[2] - compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[4] - compressedWaveform.get_time().value()[3], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[2] == compressedWaveform.get_data()[3] &&
        compressedWaveform.get_data()[4] == compressedWaveform.get_data()[5] &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[5]) {
            return WaveformLabel::BIPOLAR_TRIANGULAR;
    }
    else if (compressedWaveform.get_data().size() == 5 &&
        is_closed_enough(compressedWaveform.get_time().value()[0], compressedWaveform.get_time().value()[1], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[1] < compressedWaveform.get_data()[2] &&
        is_closed_enough(compressedWaveform.get_time().value()[2], compressedWaveform.get_time().value()[3], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[3] == compressedWaveform.get_data()[4] &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
            return WaveformLabel::FLYBACK_PRIMARY;
    }
    else if (compressedWaveform.get_data().size() == 5 &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[1] &&
        is_closed_enough(compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[2], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[2] > compressedWaveform.get_data()[3] &&
        is_closed_enough(compressedWaveform.get_time().value()[3], compressedWaveform.get_time().value()[4], 1.5 * period / settings->get_inputs_number_points_sampled_waveforms()) &&
        compressedWaveform.get_data()[0] == waveform.get_data()[4]) {
            return WaveformLabel::FLYBACK_SECONDARY;
    }
    else {
        auto settings = OpenMagnetics::Settings::GetInstance();
        double error = 0;
        double area = 0;
        double maximum = *max_element(waveform.get_data().begin(), waveform.get_data().end());
        double minimum = *min_element(waveform.get_data().begin(), waveform.get_data().end());

        double peakToPeak = maximum - minimum;
        double offset = (maximum - minimum) / 2;

        for (size_t i = 0; i < waveform.get_data().size(); ++i) {
            double angle = i * 2 * std::numbers::pi / settings->get_inputs_number_points_sampled_waveforms();
            double calculated_data = (sin(angle) * peakToPeak / 2) + offset;
            area += fabs(waveform.get_data()[i]);
            error += fabs(calculated_data - waveform.get_data()[i]);
        }
        error /= waveform.get_data().size();
        error /= area;
        if (error < 0.05) {
            return WaveformLabel::SINUSOIDAL;
        }
        else {
            return WaveformLabel::CUSTOM;
        }
    }
}

void InputsWrapper::scale_time_to_frequency(InputsWrapper& inputs, double newFrequency){
    for (auto& operatingPoint : inputs.get_mutable_operating_points()) {
        OpenMagnetics::InputsWrapper::scale_time_to_frequency(operatingPoint, newFrequency);
    }}

void InputsWrapper::scale_time_to_frequency(OperatingPoint& operatingPoint, double newFrequency){
    for (auto& excitation : operatingPoint.get_mutable_excitations_per_winding()) {
        scale_time_to_frequency(excitation, newFrequency);
    }}

void InputsWrapper::scale_time_to_frequency(OperatingPointExcitation& excitation, double newFrequency){
    excitation.set_frequency(newFrequency);
    if (excitation.get_current() && excitation.get_current()->get_waveform()) {
        auto current = excitation.get_current().value();
        current.set_waveform(scale_time_to_frequency(current.get_waveform().value(), newFrequency));
        excitation.set_current(current);
    }
    if (excitation.get_voltage() && excitation.get_voltage()->get_waveform()) {
        auto voltage = excitation.get_voltage().value();
        voltage.set_waveform(scale_time_to_frequency(voltage.get_waveform().value(), newFrequency));
        excitation.set_voltage(voltage);
    }
    if (excitation.get_magnetizing_current() && excitation.get_magnetizing_current()->get_waveform()) {
        auto magnetizingCurrent = excitation.get_magnetizing_current().value();
        magnetizingCurrent.set_waveform(scale_time_to_frequency(magnetizingCurrent.get_waveform().value(), newFrequency));
        excitation.set_magnetizing_current(magnetizingCurrent);
    }
    if (excitation.get_magnetic_flux_density() && excitation.get_magnetic_flux_density()->get_waveform()) {
        auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();
        magneticFluxDensity.set_waveform(scale_time_to_frequency(magneticFluxDensity.get_waveform().value(), newFrequency));
        excitation.set_magnetic_flux_density(magneticFluxDensity);
    }
    if (excitation.get_magnetic_field_strength() && excitation.get_magnetic_field_strength()->get_waveform()) {
        auto magneticFieldStrength = excitation.get_magnetic_field_strength().value();
        magneticFieldStrength.set_waveform(scale_time_to_frequency(magneticFieldStrength.get_waveform().value(), newFrequency));
        excitation.set_magnetic_field_strength(magneticFieldStrength);
    }
}
Waveform InputsWrapper::scale_time_to_frequency(Waveform waveform, double newFrequency){
    std::vector<double> scaledTime;
    std::vector<double> time = waveform.get_time().value();
    double oldFrequency = 1.0 / (time.back() - time.front());
    for (auto& timePoint : time) {
        scaledTime.push_back(timePoint * oldFrequency / newFrequency);
    }
    waveform.set_time(scaledTime);
    return waveform;
}


double InputsWrapper::get_maximum_voltage_peak() {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumVoltage = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        for (auto& excitation : operatingPoint.get_excitations_per_winding()) {
            if (!excitation.get_voltage()) 
                throw std::invalid_argument("Missing voltage in excitation");
            if (!excitation.get_voltage()->get_processed()) 
                throw std::invalid_argument("Voltage has not been processed");
            maximumVoltage = std::max(maximumVoltage, excitation.get_voltage()->get_processed()->get_peak().value());
        }
    }

    return maximumVoltage;
}

double InputsWrapper::get_maximum_voltage_rms() {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumVoltage = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        for (auto& excitation : operatingPoint.get_excitations_per_winding()) {
            if (!excitation.get_voltage()) 
                throw std::invalid_argument("Missing voltage in excitation");
            if (!excitation.get_voltage()->get_processed()) 
                throw std::invalid_argument("Voltage has not been processed");
            maximumVoltage = std::max(maximumVoltage, excitation.get_voltage()->get_processed()->get_rms().value());
        }
    }

    return maximumVoltage;
}

double InputsWrapper::get_maximum_current_rms() {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumCurrentRms = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        for (auto& excitation : operatingPoint.get_excitations_per_winding()) {
            if (!excitation.get_current()) 
                throw std::invalid_argument("Missing current in excitation");
            if (!excitation.get_current()->get_processed()) 
                throw std::invalid_argument("Current has not been processed");
            maximumCurrentRms = std::max(maximumCurrentRms, excitation.get_current()->get_processed()->get_rms().value());
        }
    }
    return maximumCurrentRms;
}

double InputsWrapper::get_maximum_current_peak() {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumCurrentPeak = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        for (auto& excitation : operatingPoint.get_excitations_per_winding()) {
            if (!excitation.get_current()) 
                throw std::invalid_argument("Missing current in excitation");
            if (!excitation.get_current()->get_processed()) 
                throw std::invalid_argument("Current has not been processed");
            maximumCurrentPeak = std::max(maximumCurrentPeak, excitation.get_current()->get_processed()->get_peak().value());
        }
    }
    return maximumCurrentPeak;
}

double InputsWrapper::get_maximum_voltage_peak(size_t windingIndex) {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumVoltage = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        if (!excitation.get_voltage()) 
            throw std::invalid_argument("Missing voltage in excitation");
        if (!excitation.get_voltage()->get_processed()) 
            throw std::invalid_argument("Voltage has not been processed");
        maximumVoltage = std::max(maximumVoltage, excitation.get_voltage()->get_processed()->get_peak().value());
    }

    return maximumVoltage;
}

double InputsWrapper::get_maximum_voltage_rms(size_t windingIndex) {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumVoltage = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        if (!excitation.get_voltage()) 
            throw std::invalid_argument("Missing voltage in excitation");
        if (!excitation.get_voltage()->get_processed()) 
            throw std::invalid_argument("Voltage has not been processed");
        maximumVoltage = std::max(maximumVoltage, excitation.get_voltage()->get_processed()->get_rms().value());
    }

    return maximumVoltage;
}

double InputsWrapper::get_maximum_current_rms(size_t windingIndex) {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumCurrentRms = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        if (!excitation.get_current()) 
            throw std::invalid_argument("Missing current in excitation");
        if (!excitation.get_current()->get_processed()) 
            throw std::invalid_argument("Current has not been processed");
        maximumCurrentRms = std::max(maximumCurrentRms, excitation.get_current()->get_processed()->get_rms().value());
    }
    return maximumCurrentRms;
}

double InputsWrapper::get_maximum_current_peak(size_t windingIndex) {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumCurrentPeak = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        if (!excitation.get_current()) 
            throw std::invalid_argument("Missing current in excitation");
        if (!excitation.get_current()->get_processed()) 
            throw std::invalid_argument("Current has not been processed");
        maximumCurrentPeak = std::max(maximumCurrentPeak, excitation.get_current()->get_processed()->get_peak().value());
    }
    return maximumCurrentPeak;
}

double InputsWrapper::get_maximum_current_effective_frequency() {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumCurrentEffectiveFrequency = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        for (auto& excitation : operatingPoint.get_excitations_per_winding()) {
            if (!excitation.get_current()) 
                throw std::invalid_argument("Missing current in excitation");
            if (!excitation.get_current()->get_processed()) 
                throw std::invalid_argument("Current has not been processed");
            maximumCurrentEffectiveFrequency = std::max(maximumCurrentEffectiveFrequency, excitation.get_current()->get_processed()->get_effective_frequency().value());
        }
    }
    return maximumCurrentEffectiveFrequency;
}

double InputsWrapper::get_maximum_frequency() {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumFrequency = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        for (auto& excitation : operatingPoint.get_excitations_per_winding()) {
            maximumFrequency = std::max(maximumFrequency, excitation.get_frequency());
        }
    }

    return maximumFrequency;
}

DimensionWithTolerance InputsWrapper::get_altitude() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_altitude())
        throw std::invalid_argument("Missing altitude in insulation requirements");

    return get_design_requirements().get_insulation()->get_altitude().value();
}

Cti InputsWrapper::get_cti() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_cti())
        throw std::invalid_argument("Missing cti in insulation requirements");

    return get_design_requirements().get_insulation()->get_cti().value();
}

InsulationType InputsWrapper::get_insulation_type() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_insulation_type())
        throw std::invalid_argument("Missing insulation_type in insulation requirements");

    return get_design_requirements().get_insulation()->get_insulation_type().value();
}

DimensionWithTolerance InputsWrapper::get_main_supply_voltage() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_main_supply_voltage())
        throw std::invalid_argument("Missing main_supply_voltage in insulation requirements");

    return get_design_requirements().get_insulation()->get_main_supply_voltage().value();
}

OvervoltageCategory InputsWrapper::get_overvoltage_category() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_overvoltage_category())
        throw std::invalid_argument("Missing overvoltage_category in insulation requirements");

    return get_design_requirements().get_insulation()->get_overvoltage_category().value();
}

PollutionDegree InputsWrapper::get_pollution_degree() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_pollution_degree())
        throw std::invalid_argument("Missing pollution_degree in insulation requirements");

    return get_design_requirements().get_insulation()->get_pollution_degree().value();
}

std::vector<InsulationStandards> InputsWrapper::get_standards() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_standards())
        throw std::invalid_argument("Missing standards in insulation requirements");

    return get_design_requirements().get_insulation()->get_standards().value();
}

void InputsWrapper::set_current_as_magnetizing_current(OperatingPoint* operatingPoint) {
    OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(*operatingPoint);
    auto current = excitation.get_current().value();

    auto currentExcitation = excitation.get_current().value();
    auto currentExcitationWaveform = currentExcitation.get_waveform().value();
    auto sampledCurrentWaveform = InputsWrapper::calculate_sampled_waveform(currentExcitationWaveform, excitation.get_frequency());

    if (sampledCurrentWaveform.get_data().size() > 0 && ((sampledCurrentWaveform.get_data().size() & (sampledCurrentWaveform.get_data().size() - 1)) != 0)) {
        throw std::invalid_argument("sampledCurrentWaveform vector size is not a power of 2");
    }

    currentExcitation.set_harmonics(InputsWrapper::calculate_harmonics_data(sampledCurrentWaveform, excitation.get_frequency()));
    currentExcitation.set_processed(InputsWrapper::calculate_processed_data(currentExcitation, sampledCurrentWaveform, true));
    excitation.set_current(currentExcitation);
    excitation.set_magnetizing_current(excitation.get_current().value());
    operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
}

double InputsWrapper::get_switching_frequency(OperatingPointExcitation excitation) {
    if (excitation.get_current()) {
        if (excitation.get_current()->get_waveform()) {
            auto waveform = excitation.get_current()->get_waveform().value();
            if (waveform.get_data().size() > Constants().numberPointsSampledWaveforms) {
                if (excitation.get_current()->get_harmonics()) {
                    auto harmonics = excitation.get_current()->get_harmonics().value();
                    double mainHarmonicAmplitude = harmonics.get_amplitudes()[1];
                    double strongestHarmonicAmplitudeAfterMain = 0;
                    double strongestHarmonicAmplitudeAfterMainFrequency = harmonics.get_frequencies()[1];
                    for (size_t harmonicIndex = 2; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
                        if (harmonics.get_amplitudes()[harmonicIndex] > 0.01 * mainHarmonicAmplitude) {
                            if (harmonics.get_amplitudes()[harmonicIndex] > strongestHarmonicAmplitudeAfterMain) {
                                strongestHarmonicAmplitudeAfterMain = harmonics.get_amplitudes()[harmonicIndex];
                                strongestHarmonicAmplitudeAfterMainFrequency = harmonics.get_frequencies()[harmonicIndex];
                            }
                        }
                    }
                    return strongestHarmonicAmplitudeAfterMainFrequency;
                }
            }
        }
    }
    return excitation.get_frequency();
}

double InputsWrapper::get_magnetic_flux_density_peak(OperatingPointExcitation excitation, double switchingFrequency) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();

    if (excitation.get_frequency() != switchingFrequency) {
        if (!excitation.get_magnetic_flux_density()->get_harmonics()) {
            auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value();
            auto sampledWaveform = InputsWrapper::calculate_sampled_waveform(magneticFluxDensityWaveform, excitation.get_frequency());
            magneticFluxDensity.set_harmonics(calculate_harmonics_data(sampledWaveform, excitation.get_frequency()));
            excitation.set_magnetic_flux_density(magneticFluxDensity);
        }
        auto harmonics = excitation.get_magnetic_flux_density()->get_harmonics().value();
        for (size_t harmonicIndex = 2; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            if (harmonics.get_frequencies()[harmonicIndex] == switchingFrequency) {
                return harmonics.get_amplitudes()[harmonicIndex];
            }
        }
    }

    return magneticFluxDensity.get_processed().value().get_peak().value();
}

double InputsWrapper::get_magnetic_flux_density_peak_to_peak(OperatingPointExcitation excitation, double switchingFrequency) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();

    if (excitation.get_frequency() != switchingFrequency) {
        if (!excitation.get_magnetic_flux_density()->get_harmonics()) {
            auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value();
            auto sampledWaveform = InputsWrapper::calculate_sampled_waveform(magneticFluxDensityWaveform, excitation.get_frequency());
            magneticFluxDensity.set_harmonics(calculate_harmonics_data(sampledWaveform, excitation.get_frequency()));
            excitation.set_magnetic_flux_density(magneticFluxDensity);
        }
        auto harmonics = excitation.get_magnetic_flux_density()->get_harmonics().value();
        for (size_t harmonicIndex = 2; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            if (harmonics.get_frequencies()[harmonicIndex] == switchingFrequency) {
                return harmonics.get_amplitudes()[harmonicIndex] * 2;
            }
        }
    }

    return magneticFluxDensity.get_processed().value().get_peak_to_peak().value();
}

} // namespace OpenMagnetics
