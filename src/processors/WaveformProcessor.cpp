#include "processors/WaveformProcessor.h"
#include "support/Exceptions.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdint>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

// Kept C++17-compatible (no <numbers>/std::lerp) so the Kirchhoff converter-model library, which builds
// at C++17, can reuse this file; MKF (C++23) compiles it unchanged.
namespace { constexpr double kWaveformPi = 3.14159265358979323846264338327950288; }
static inline double kWaveformLerp(double a, double b, double t) { return a + t * (b - a); }

// ---------------------------------------------------------------------------
// Pure / shared helpers.
//
// These are duplicated here (in an anonymous namespace, internal linkage) so
// that WaveformProcessor stays self-contained and free of support/Utils.h and
// the magnetics stack. The originals still live in support/Utils.* and
// processors/Inputs.cpp for the rest of MKF; these copies are byte-for-byte
// equivalent so behavior is identical.
// ---------------------------------------------------------------------------
namespace {

// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// from https://rosettacode.org/wiki/Fast_Fourier_transform#C++
void fft(std::vector<std::complex<double>>& x) {
    // DFT
    unsigned int N = x.size(), k = N, n;
    double thetaT = kWaveformPi / N;
    std::complex<double> phiT = std::complex<double>(cos(thetaT), -sin(thetaT)), T;
    while (k > 1) {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0L;
        for (unsigned int l = 0; l < k; l++) {
            for (long a = l; a < N; a += n) {
                long b = a + k;
                std::complex<double> t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    // Decimate
    unsigned int m = (unsigned int) log2(N);
    for (unsigned int a = 0; a < N; a++) {
        unsigned int b = a;
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

double roundFloat(double value, int64_t decimals = 9) {
    return round(value * pow(10, decimals)) / pow(10, decimals);
}

bool is_size_power_of_2(std::vector<double> data) {
    return data.size() > 0 && ((data.size() & (data.size() - 1)) == 0);
}

size_t round_up_size_to_power_of_2(std::vector<double> data) {
    return pow(2, ceil(log(data.size()) / log(2)));
}

std::vector<double> linear_spaced_array(double a, double b, size_t N) {
    double h = (b - a) / static_cast<double>(N-1);
    std::vector<double> xs(N);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}

std::vector<size_t> get_main_harmonic_indexes(Harmonics harmonics, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex = std::nullopt) {
    double mainOrMaximumHarmonicAmplitudeTimesRootFrequency = 0;
    if (!mainHarmonicIndex) {
        for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            mainOrMaximumHarmonicAmplitudeTimesRootFrequency = std::max(harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex]), mainOrMaximumHarmonicAmplitudeTimesRootFrequency);
        }
    }
    else {
        mainOrMaximumHarmonicAmplitudeTimesRootFrequency = harmonics.get_amplitudes()[mainHarmonicIndex.value()] * sqrt(harmonics.get_frequencies()[mainHarmonicIndex.value()]);
    }

    std::vector<size_t> mainHarmonicIndexes;
    for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
        if ((harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex])) < mainOrMaximumHarmonicAmplitudeTimesRootFrequency * windingLossesHarmonicAmplitudeThreshold) {
            continue;
        }

        if (std::find(mainHarmonicIndexes.begin(), mainHarmonicIndexes.end(), harmonicIndex) == mainHarmonicIndexes.end()) {
            mainHarmonicIndexes.push_back(harmonicIndex);
        }
    }
    return mainHarmonicIndexes;
}

bool is_close_enough(double x, double y, double error){
    return fabs(x - y) <= error;
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
            if (waveform.get_data().size() == 4) {
                return waveform.get_data()[0];
            }
            else {
                return waveform.get_data()[1];
            }
        case WaveformLabel::FLYBACK_SECONDARY:
            return waveform.get_data()[3];
    // FIXED: BUG-06 — Handle missing waveform labels
    case WaveformLabel::RECTANGULAR_WITH_DEADTIME:
    case WaveformLabel::SECONDARY_RECTANGULAR:
    case WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME:
        return 0; // Rectangular variants: zero offset like RECTANGULAR
    case WaveformLabel::TRIANGULAR_WITH_DEADTIME:
    case WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME:
        return (*max_element(waveform.get_data().begin(), waveform.get_data().end()) +
                *min_element(waveform.get_data().begin(), waveform.get_data().end())) / 2;
    }
    return 0;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Waveform-shape helpers (public; ported from Inputs).
// ---------------------------------------------------------------------------

bool WaveformProcessor::is_waveform_imported(Waveform waveform, size_t numberPointsSampledWaveforms) {
    if (!waveform.get_time()) {
        return false;
    }
    else {
        return waveform.get_data().size() > numberPointsSampledWaveforms;
    }
}

bool WaveformProcessor::is_waveform_sampled(Waveform waveform, size_t numberPointsSampledWaveforms) {
    if (!waveform.get_time()) {
        return false;
    }
    else if (is_waveform_imported(waveform, numberPointsSampledWaveforms)) {
        return is_size_power_of_2(waveform.get_data());
    }
    else {
        return waveform.get_data().size() == numberPointsSampledWaveforms;
    }
}

double WaveformProcessor::calculate_waveform_average(Waveform waveform) {
    double integration = 0;
    double period = waveform.get_time()->back() - waveform.get_time()->front();
    for (size_t i = 0; i < waveform.get_data().size() - 1; ++i)
    {
        double area = (waveform.get_data()[i + 1] + waveform.get_data()[i]) / 2 * (waveform.get_time().value()[i + 1] - waveform.get_time().value()[i]);
        integration += area;
    }
    return integration / period;
}

Waveform WaveformProcessor::compress_waveform(const Waveform& waveform) {
    Waveform compressedWaveform;
    auto data = waveform.get_data();
    if (data.size() < 2) {
        return waveform;
    }
    data.push_back(data[0]);
    auto time = waveform.get_time().value();
    if (time.size() < 2) {
        return waveform;
    }
    time.push_back(time[time.size() - 1] + (time[time.size() - 1] - time[time.size() - 2]));
    std::vector<double> compressedData;
    std::vector<double> compressedTime;
    double previousSlope = DBL_MAX;
    for (size_t i = 0; i < data.size() - 1; ++i){
        // FIXED: BUG-01 — Guard against division by zero on identical time points
        double slope;
        if (time[i + 1] == time[i]) {
            slope = (data[i + 1] >= data[i]) ? DBL_MAX : -DBL_MAX;
        } else {
            slope = (data[i + 1] - data[i]) / (time[i + 1] - time[i]);
        }

        // Emit point if: (1) discontinuity, OR (2) slope changes significantly
        bool isDiscontinuity = (slope == DBL_MAX || slope == -DBL_MAX);
        bool previousWasDiscontinuity = (previousSlope == DBL_MAX || previousSlope == -DBL_MAX);

        // Detect slope change: use relative change if possible, absolute if previous slope is near zero
        bool slopeChanged = false;
        if (!isDiscontinuity && !previousWasDiscontinuity) {
            if (fabs(previousSlope) > 1e-9) {
                // Non-zero previous slope: use relative change
                slopeChanged = fabs((slope - previousSlope) / previousSlope) > 0.01;
            } else {
                // Previous slope near zero: use absolute change
                slopeChanged = fabs(slope) > 1e-9;
            }
        }
        if (isDiscontinuity || previousWasDiscontinuity || slopeChanged) {
            compressedData.push_back(data[i]);
            compressedTime.push_back(time[i]);
        }
        previousSlope = slope;
    }
    compressedData.push_back(data.back());
    compressedTime.push_back(time.back());
    compressedWaveform.set_data(compressedData);
    compressedWaveform.set_time(compressedTime);
    if (waveform.get_ancillary_label()) {
        compressedWaveform.set_ancillary_label(waveform.get_ancillary_label().value());
    }
    return compressedWaveform;
}

double WaveformProcessor::try_guess_duty_cycle(Waveform waveform, WaveformLabel label, double frequency, size_t numberPointsSampledWaveforms) {
    if (label != WaveformLabel::CUSTOM) {
        switch(label) {
            case WaveformLabel::TRIANGULAR: {
                if (waveform.get_time()->size() == 3) {
                    return (waveform.get_time().value()[1] - waveform.get_time()->front()) / (waveform.get_time()->back() - waveform.get_time()->front());
                }
                else if (waveform.get_time()->size() == 4) {
                    return ((waveform.get_time().value()[1] + waveform.get_time().value()[2]) / 2 - waveform.get_time()->front()) / (waveform.get_time()->back() - waveform.get_time()->front());
                }
                break;
            }
            case WaveformLabel::UNIPOLAR_TRIANGULAR: {
                return (waveform.get_time().value()[1] - waveform.get_time()->front()) / (waveform.get_time()->back() - waveform.get_time()->front());
            }
            case WaveformLabel::RECTANGULAR: {
                return (waveform.get_time().value()[2] - waveform.get_time()->front()) / (waveform.get_time()->back() - waveform.get_time()->front());
            }
            case WaveformLabel::UNIPOLAR_RECTANGULAR: {
                return (waveform.get_time().value()[2] - waveform.get_time()->front()) / (waveform.get_time()->back() - waveform.get_time()->front());
            }
            case WaveformLabel::BIPOLAR_RECTANGULAR: {
                return (waveform.get_time().value()[3] - waveform.get_time().value()[2]) / (waveform.get_time()->back() - waveform.get_time()->front());
            }
            case WaveformLabel::BIPOLAR_TRIANGULAR: {
                return (waveform.get_time().value()[2] - waveform.get_time().value()[1]) / (waveform.get_time()->back() - waveform.get_time()->front());
            }
            case WaveformLabel::FLYBACK_PRIMARY:{
                if (waveform.get_time()->size() == 4) {
                    return (waveform.get_time().value()[1] - waveform.get_time()->front()) / (waveform.get_time()->back() - waveform.get_time()->front());
                }
                else if (waveform.get_time()->size() == 5) {
                    return (waveform.get_time().value()[2] - waveform.get_time()->front()) / (waveform.get_time()->back() - waveform.get_time()->front());
                }
                break;
            }
            case WaveformLabel::FLYBACK_SECONDARY:{
                if (waveform.get_time()->size() == 4) {
                    return (waveform.get_time().value()[1] - waveform.get_time()->front()) / (waveform.get_time()->back() - waveform.get_time()->front());
                }
                else if (waveform.get_time()->size() == 5) {
                    return (waveform.get_time().value()[2] - waveform.get_time()->front()) / (waveform.get_time()->back() - waveform.get_time()->front());
                }
                break;
            }
            case WaveformLabel::SINUSOIDAL: {
                return 0.5;
            }
            default:
                break;
            }
    }

    Waveform sampledWaveform;
    if (!is_waveform_sampled(waveform, numberPointsSampledWaveforms)) {
        if (frequency > 0) {
            sampledWaveform = WaveformProcessor::calculate_sampled_waveform(waveform, frequency, std::nullopt, numberPointsSampledWaveforms);
        }
        else if (waveform.get_time() && waveform.get_time()->size() >= 3) {
            // Can compute duty cycle directly from the waveform time points without sampling
            auto timeVec = waveform.get_time().value();
            auto dataVec = waveform.get_data();
            double totalPeriod = timeVec.back() - timeVec.front();
            if (totalPeriod <= 0) {
                return 0.5;
            }
            if (dataVec.size() == 3) {
                // 3-point triangular: peak at middle time point
                return roundFloat((timeVec[1] - timeVec[0]) / totalPeriod, 2);
            }
            else if (dataVec.size() == 4) {
                // 4-point waveform: duty cycle from average of middle two points
                return roundFloat(((timeVec[1] + timeVec[2]) / 2 - timeVec[0]) / totalPeriod, 2);
            }
            else if (dataVec.size() == 5) {
                // 5-point rectangular: time from point 0 to 2 over total
                return roundFloat((timeVec[2] - timeVec[0]) / totalPeriod, 2);
            }
            else {
                // Cannot compute duty cycle without frequency for arbitrary waveforms
                return 0.5;
            }
        }
        else {
            // Cannot sample without frequency, return default duty cycle
            return 0.5;
        }
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

    double maximum = *max_element(diff_diff_data.begin(), diff_diff_data.end());
    size_t maximum_index = 0;
    size_t distanceToMiddle = numberPointsSampledWaveforms;
    for (size_t i = 0; i < diff_diff_data.size(); ++i)
    {
        if (diff_diff_data[i] == maximum) {
            if (fabs(double(numberPointsSampledWaveforms) / 2 - i) < distanceToMiddle) {
                distanceToMiddle = fabs(double(numberPointsSampledWaveforms) / 2 - i);
                maximum_index = i;
            }
        }
    }
    auto dutyCycle = roundFloat((maximum_index + 1.0) / numberPointsSampledWaveforms, 2);

    if (dutyCycle <= 0.03 || dutyCycle >= 0.97) {

        double maximum = *max_element(data.begin(), data.end());
        double threshold = maximum * 0.05;

        double numberPointsOn = 0;
        double numberPointsOff = 0;
        for (size_t i = 0; i < data.size() - 1; ++i) {
            if (data[i] < threshold) {
                numberPointsOff++;
            }
            else {
                numberPointsOn++;
            }
        }

        dutyCycle = numberPointsOn / data.size();
    }

    return dutyCycle;
}

WaveformLabel WaveformProcessor::try_guess_waveform_label(Waveform waveform, size_t numberPointsSampledWaveforms) {
    if (waveform.get_ancillary_label()) {
        return waveform.get_ancillary_label().value();
    }

    auto compressedWaveform = waveform;
    if (is_waveform_sampled(waveform, numberPointsSampledWaveforms))
        compressedWaveform = compress_waveform(waveform);
    double period = 0;
    if (compressedWaveform.get_time()) {
        period = compressedWaveform.get_time()->back() - compressedWaveform.get_time()->front();
    }

    if (compressedWaveform.get_data().size() == 3 &&
        compressedWaveform.get_data()[0] == compressedWaveform.get_data()[2]) {
            return WaveformLabel::TRIANGULAR;
    }
    else if (compressedWaveform.get_time()) {
        if (compressedWaveform.get_data().size() == 4 &&
            is_close_enough(compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[2], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[2] == compressedWaveform.get_data()[3] &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[3]) {
                return WaveformLabel::UNIPOLAR_TRIANGULAR;
        }
        else if (compressedWaveform.get_data().size() == 5 &&
            !is_close_enough((compressedWaveform.get_time().value()[2] - compressedWaveform.get_time().value()[0]) * compressedWaveform.get_data()[2] + (compressedWaveform.get_time().value()[4] - compressedWaveform.get_time().value()[2]) * compressedWaveform.get_data()[4], 0 , period) &&
            is_close_enough(compressedWaveform.get_time().value()[0], compressedWaveform.get_time().value()[1], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[1] == compressedWaveform.get_data()[2] &&
            is_close_enough(compressedWaveform.get_time().value()[2], compressedWaveform.get_time().value()[3], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[3] == compressedWaveform.get_data()[4] &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
                return WaveformLabel::UNIPOLAR_RECTANGULAR;
        }
        else if (compressedWaveform.get_data().size() == 5 &&
            !is_close_enough((compressedWaveform.get_time().value()[2] - compressedWaveform.get_time().value()[0]) * compressedWaveform.get_data()[2] + (compressedWaveform.get_time().value()[4] - compressedWaveform.get_time().value()[2]) * compressedWaveform.get_data()[4], 0 , period) &&
            is_close_enough(compressedWaveform.get_time().value()[0], compressedWaveform.get_time().value()[1], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[1] == compressedWaveform.get_data()[2] &&
            is_close_enough(compressedWaveform.get_time().value()[2], compressedWaveform.get_time().value()[3], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[3] == compressedWaveform.get_data()[4] &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
                return WaveformLabel::UNIPOLAR_RECTANGULAR;
        }
        else if (compressedWaveform.get_data().size() == 5 &&
            is_close_enough((compressedWaveform.get_time().value()[2] - compressedWaveform.get_time().value()[0]) * compressedWaveform.get_data()[2] + (compressedWaveform.get_time().value()[4] - compressedWaveform.get_time().value()[2]) * compressedWaveform.get_data()[4], 0 , period) &&
            is_close_enough(compressedWaveform.get_time().value()[0], compressedWaveform.get_time().value()[1], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[1] == compressedWaveform.get_data()[2] &&
            is_close_enough(compressedWaveform.get_time().value()[2], compressedWaveform.get_time().value()[3], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[3] == compressedWaveform.get_data()[4] &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
                return WaveformLabel::RECTANGULAR;
        }
        else if (compressedWaveform.get_data().size() == 5 &&
            is_close_enough((compressedWaveform.get_time().value()[1] - compressedWaveform.get_time().value()[0]) * compressedWaveform.get_data()[1] + (compressedWaveform.get_time().value()[3] - compressedWaveform.get_time().value()[2]) * compressedWaveform.get_data()[3], 0 , period) &&
            is_close_enough(compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[2], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[1] &&
            is_close_enough(compressedWaveform.get_time().value()[3], compressedWaveform.get_time().value()[4], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[2] == compressedWaveform.get_data()[3] &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
                return WaveformLabel::RECTANGULAR;
        }
        else if (compressedWaveform.get_data().size() == 10 &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[1] &&
            is_close_enough(compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[2], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[2] == compressedWaveform.get_data()[3] &&
            is_close_enough(compressedWaveform.get_time().value()[3], compressedWaveform.get_time().value()[4], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[4] == compressedWaveform.get_data()[5] &&
            is_close_enough(compressedWaveform.get_time().value()[5], compressedWaveform.get_time().value()[6], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[6] == compressedWaveform.get_data()[7] &&
            is_close_enough(compressedWaveform.get_time().value()[7], compressedWaveform.get_time().value()[8], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[8] == compressedWaveform.get_data()[9] &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[9]) {
                return WaveformLabel::BIPOLAR_RECTANGULAR;
        }
        else if (compressedWaveform.get_data().size() == 6 &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[1] &&
            is_close_enough(compressedWaveform.get_time().value()[2] - compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[4] - compressedWaveform.get_time().value()[3], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[2] == compressedWaveform.get_data()[3] &&
            compressedWaveform.get_data()[4] == compressedWaveform.get_data()[5] &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[5]) {
                return WaveformLabel::BIPOLAR_TRIANGULAR;
        }
        else if (compressedWaveform.get_data().size() == 5 &&
            is_close_enough(compressedWaveform.get_time().value()[0], compressedWaveform.get_time().value()[1], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[1] < compressedWaveform.get_data()[2] &&
            is_close_enough(compressedWaveform.get_time().value()[2], compressedWaveform.get_time().value()[3], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[3] == compressedWaveform.get_data()[4] &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
                return WaveformLabel::FLYBACK_PRIMARY;
        }
        else if (compressedWaveform.get_data().size() == 5 &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[1] &&
            is_close_enough(compressedWaveform.get_time().value()[1], compressedWaveform.get_time().value()[2], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[2] > compressedWaveform.get_data()[3] &&
            is_close_enough(compressedWaveform.get_time().value()[3], compressedWaveform.get_time().value()[4], 1.5 * period / numberPointsSampledWaveforms) &&
            compressedWaveform.get_data()[0] == compressedWaveform.get_data()[4]) {
                return WaveformLabel::FLYBACK_SECONDARY;
        }
        else {
            double error = 0;
            double area = 0;
            double maximum = *max_element(waveform.get_data().begin(), waveform.get_data().end());
            double minimum = *min_element(waveform.get_data().begin(), waveform.get_data().end());

            double peakToPeak = maximum - minimum;
            double offset = (maximum + minimum) / 2; // FIXED: BUG-08

            for (size_t i = 0; i < waveform.get_data().size(); ++i) {
                double angle = i * 2 * kWaveformPi / numberPointsSampledWaveforms;
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
    else {
        double error = 0;
        double area = 0;
        double maximum = *max_element(waveform.get_data().begin(), waveform.get_data().end());
        double minimum = *min_element(waveform.get_data().begin(), waveform.get_data().end());

        double peakToPeak = maximum - minimum;
        double offset = (maximum + minimum) / 2; // FIXED: BUG-08

        for (size_t i = 0; i < waveform.get_data().size(); ++i) {
            double angle = i * 2 * kWaveformPi / numberPointsSampledWaveforms;
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

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

Waveform WaveformProcessor::create_waveform(ProcessedWaveform processed, double frequency, size_t numberPointsSampledWaveforms) {
    if (!processed.get_peak_to_peak()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Signal is missing peak to peak");
    }

    auto label = processed.get_label();
    auto peakToPeak = processed.get_peak_to_peak().value();
    auto offset = processed.get_offset();
    double dutyCycle = 0.5;
    if (processed.get_duty_cycle()) {
        dutyCycle = processed.get_duty_cycle().value();
    }
    double deadTime = 0;
    if (processed.get_dead_time()) {
        deadTime = processed.get_dead_time().value();
    }

    return create_waveform(label, peakToPeak, frequency, dutyCycle, offset, deadTime, 0, 0, numberPointsSampledWaveforms);
}

Waveform WaveformProcessor::create_waveform(WaveformLabel label, double peakToPeak, double frequency, double dutyCycle, double offset, double deadTime, double skew, double phase, size_t numberPointsSampledWaveforms) {
    Waveform waveform;
    std::vector<double> data;
    std::vector<double> time;
    double period = 1 / frequency;

    switch (label) {
        case WaveformLabel::TRIANGULAR: {
            double max = peakToPeak / 2 + offset;
            double min = -peakToPeak / 2 + offset;
            double dc = dutyCycle * period;
            data = {min, max, min};
            time = {0, dc, period};
            break;
        }
        case WaveformLabel::TRIANGULAR_WITH_DEADTIME: {
            double max = peakToPeak / 2 + offset;
            double min = -peakToPeak / 2 + offset;
            double dc = dutyCycle * period;
            data = {min, max, min, 0};
            time = {0, dc, period - deadTime, period};
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
            double max = peakToPeak * (1 - dutyCycle) + offset;
            double min = -peakToPeak * dutyCycle + offset;
            double dc = dutyCycle * period;
            data = {min, max, max, min, min};
            time = {0, 0, dc, dc, period};
            break;
        }
        case WaveformLabel::RECTANGULAR_WITH_DEADTIME: {
            double max = peakToPeak * (1 - dutyCycle) + offset;
            double min = -peakToPeak * dutyCycle + offset;
            double dc = dutyCycle * period;
            data = {0, max, max, min, min, 0, 0};
            time = {0, 0, dc, dc, period - deadTime, period - deadTime, period};
            break;
        }
        case WaveformLabel::SECONDARY_RECTANGULAR: {
            double max = -peakToPeak * (1 - dutyCycle) + offset;
            double min = peakToPeak * dutyCycle + offset;
            double dc = dutyCycle * period;
            data = {min, max, max, min, min};
            time = {0, 0, dc, dc, period};
            break;
        }
        case WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME: {
            double max = -peakToPeak * (1 - dutyCycle) + offset;
            double min = peakToPeak * dutyCycle + offset;
            double dc = dutyCycle * period;
            data = {0, max, max, min, min, 0, 0};
            time = {0, 0, dc, dc, period - deadTime, period - deadTime, period};
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
        case WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME:{
            double max = peakToPeak + offset;
            double min = offset;
            double dc = dutyCycle * period;
            data = {0, 0, max, min, 0, 0};
            time = {0, dc, dc, period - deadTime, period - deadTime, period};
            break;
        }
        case WaveformLabel::SINUSOIDAL: {
            for (size_t i = 0; i < numberPointsSampledWaveforms; ++i) {
                double angle = i * 2 * kWaveformPi / (numberPointsSampledWaveforms - 1);
                time.push_back(i * period / (numberPointsSampledWaveforms - 1));
                data.push_back((sin(angle + phase) * peakToPeak / 2) + offset);
            }
            break;
        }
        default:
            break;
    }

    waveform.set_ancillary_label(label);
    waveform.set_data(data);
    waveform.set_time(time);

    if (skew > 0) {
        auto sampledWaveform = calculate_sampled_waveform(waveform, frequency, std::nullopt, numberPointsSampledWaveforms);
        auto timeStep = period / sampledWaveform.get_data().size();
        size_t stepNeededForSkew = round(skew / timeStep);

        auto auxData = sampledWaveform.get_data();
        // FIXED: BUG-05 — Removed period subtraction
        std::rotate(auxData.rbegin(), auxData.rbegin() + stepNeededForSkew, auxData.rend());
        waveform = sampledWaveform;
        waveform.set_data(auxData);
        waveform = compress_waveform(waveform);
    }

    return waveform;
}

Waveform WaveformProcessor::calculate_sampled_waveform(Waveform waveform, double frequency, std::optional<size_t> numberPoints, size_t numberPointsSampledWaveforms) {
    std::vector<double> time;
    auto data = waveform.get_data();

    // Validate input data
    if (data.size() < 2) {
        throw std::invalid_argument("Waveform must have at least 2 data points");
    }

    if (!waveform.get_time()) { // This means the waveform is equidistant
        // Without time data, frequency must be provided
        if (!std::isfinite(frequency) || frequency <= 0) {
            throw std::invalid_argument("Invalid frequency: " + std::to_string(frequency));
        }
        // Reject unreasonably small frequencies (likely uninitialized garbage values)
        // Minimum reasonable frequency is 1 Hz (period = 1 second)
        if (frequency < 1.0) {
            throw std::invalid_argument("Frequency too small (likely uninitialized): " + std::to_string(frequency));
        }
        time = linear_spaced_array(0, 1. / roundFloat(frequency, 9), data.size());
    }
    else {
        time = waveform.get_time().value();
        // The waveform time data defines the actual period to sample over
        // Always use the frequency derived from the waveform time data for sampling
        // This ensures the sampled time points fall within the waveform's time range
        double period = (time.back() - time.front());
        if (period <= 0) {
            throw std::invalid_argument("Invalid waveform period: " + std::to_string(period));
        }
        // Use frequency from waveform time data for sampling purposes
        // The passed frequency parameter is only used as fallback when time data is missing
        frequency = 1.0 / period;
        // Validate frequency after inference
        if (!std::isfinite(frequency) || frequency <= 0) {
            throw std::invalid_argument("Invalid frequency derived from waveform: " + std::to_string(frequency));
        }
    }

    size_t numberPointsForSampling = numberPointsSampledWaveforms;
    if (numberPoints) {
        numberPointsForSampling = numberPoints.value();
    }

    if (data.size() > numberPointsForSampling) {
        if (is_size_power_of_2(data)) {
            numberPointsForSampling = data.size();
        }
        else {
            numberPointsForSampling = round_up_size_to_power_of_2(data);
        }
    }

    auto sampledTime = linear_spaced_array(0, 1. / roundFloat(frequency, 9), numberPointsForSampling + 1);

    std::vector<double> sampledData;

    for (size_t i = 0; i < numberPointsForSampling; i++) {
        bool found = false;
        for (size_t interpIndex = 0; interpIndex < data.size() - 1; interpIndex++) {
            // Skip zero-length segments (where time[i] == time[i+1])
            // These occur in waveforms like FLYBACK_PRIMARY: time = {0, 0, dc, dc, period}
            if (time[interpIndex + 1] == time[interpIndex]) {
                // If sampled time is exactly at this zero-length segment, use the data value
                if (sampledTime[i] == time[interpIndex]) {
                    sampledData.push_back(data[interpIndex]);
                    found = true;
                    break;
                }
                // Otherwise skip this zero-length segment
                continue;
            }

            // For non-zero-length segments, check if sampled time falls within [start, end]
            // Use inclusive checks to handle edge cases properly
            if (time[interpIndex] <= sampledTime[i] && sampledTime[i] <= time[interpIndex + 1]) {
                double proportion = (sampledTime[i] - time[interpIndex]) / (time[interpIndex + 1] - time[interpIndex]);
                double interpPoint = kWaveformLerp(data[interpIndex], data[interpIndex + 1], proportion);
                sampledData.push_back(interpPoint);
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::invalid_argument("Error while sampling waveform in point: " + std::to_string(i));
        }
    }

    if (sampledData.size() != numberPointsForSampling) {
        throw std::invalid_argument("Wrong number of sampled points");
    }

    sampledTime.pop_back();

    Waveform sampledWaveform;
    sampledWaveform.set_data(sampledData);
    sampledWaveform.set_time(sampledTime);
    if (waveform.get_ancillary_label()) {
        sampledWaveform.set_ancillary_label(waveform.get_ancillary_label().value());
    }
    return sampledWaveform;
}

Harmonics WaveformProcessor::calculate_harmonics_data(Waveform waveform, double frequency, bool trimHarmonics, double harmonicAmplitudeThreshold, size_t numberPointsSampledWaveforms) {
    bool isWaveformImported = is_waveform_imported(waveform, numberPointsSampledWaveforms);
    Harmonics harmonics;

    std::vector<std::complex<double>> data;
    for (std::size_t i = 0; i < waveform.get_data().size(); ++i) {
        data.emplace_back(waveform.get_data()[i]);
    }

    if (data.size() > 0 && ((data.size() & (data.size() - 1)) != 0)) {
        throw std::invalid_argument("Data vector size is not a power of 2: " + std::to_string(data.size()));
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
        auto mainHarmonicIndexes = get_main_harmonic_indexes(harmonics, harmonicAmplitudeThreshold);
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

ProcessedWaveform WaveformProcessor::calculate_processed_data(SignalDescriptor excitation,
                                            Waveform sampledWaveform,
                                            bool includeAdvancedData,
                                            std::optional<ProcessedWaveform> processed,
                                            size_t numberPointsSampledWaveforms) {

    auto harmonics = excitation.get_harmonics().value();
    return calculate_processed_data(harmonics, sampledWaveform, includeAdvancedData, processed, numberPointsSampledWaveforms);
}

ProcessedWaveform WaveformProcessor::calculate_basic_processed_data(Waveform waveform, size_t numberPointsSampledWaveforms) {
    ProcessedWaveform processed;
    std::vector<double> dataToProcess;
    auto sampledWaveform = waveform;
    auto compressedWaveform = waveform;

    for (size_t i = 0; i < waveform.get_data().size(); ++i) {
        if (std::isnan(waveform.get_data()[i])) {
            throw std::invalid_argument("Waveform data contains NaN");
        }
    }

    if (is_waveform_sampled(waveform, numberPointsSampledWaveforms)) {
        compressedWaveform = compress_waveform(waveform);
    }

    WaveformLabel label;

    label = try_guess_waveform_label(compressedWaveform, numberPointsSampledWaveforms);
    processed.set_label(label);

    if (is_waveform_sampled(waveform, numberPointsSampledWaveforms)) {
        processed.set_average(std::accumulate(std::begin(sampledWaveform.get_data()), std::end(sampledWaveform.get_data()), 0.0) /
                         sampledWaveform.get_data().size());
    }
    else {
        // auto average = calculate_waveform_average(compressedWaveform);
        // processed.set_average(average);
    }

    double offset = calculate_offset(compressedWaveform, label);
    processed.set_offset(offset);

    processed.set_peak_to_peak(*max_element(compressedWaveform.get_data().begin(), compressedWaveform.get_data().end()) -
                               *min_element(compressedWaveform.get_data().begin(), compressedWaveform.get_data().end()));

    if (label == WaveformLabel::FLYBACK_PRIMARY ||
        label == WaveformLabel::FLYBACK_SECONDARY) {
        // Flyback waveforms include a 0-level floor, so max - min overstates the
        // ramp peak-to-peak by the offset. UNIPOLAR_* waveforms never touch 0
        // (their minimum IS the offset), so max - min is already correct for them.
        processed.set_peak_to_peak(processed.get_peak_to_peak().value() - offset);
    }

    double positivePeak = *max_element(compressedWaveform.get_data().begin(), compressedWaveform.get_data().end());
    double negativePeak = *min_element(compressedWaveform.get_data().begin(), compressedWaveform.get_data().end());
    processed.set_peak(std::max(positivePeak, -negativePeak));
    processed.set_positive_peak(positivePeak);
    processed.set_negative_peak(negativePeak);

    processed.set_duty_cycle(try_guess_duty_cycle(compressedWaveform, label, 0, numberPointsSampledWaveforms));

    return processed;
}

ProcessedWaveform WaveformProcessor::calculate_processed_data(Waveform waveform,
                                            std::optional<double> frequency,
                                            bool includeAdvancedData,
                                            std::optional<ProcessedWaveform> processed,
                                            bool trimHarmonics,
                                            double harmonicAmplitudeThreshold,
                                            size_t numberPointsSampledWaveforms) {
    double frequencyValue;
    if (frequency) {
        frequencyValue = frequency.value();
    }
    else {
        if (!waveform.get_time()) {
            throw std::invalid_argument("Either frequency or time must be provided");
        }
        frequencyValue = 1.0 / (waveform.get_time()->back() - waveform.get_time()->front());
    }
    auto sampledWaveform = waveform;
    if (!is_size_power_of_2(waveform.get_data())) {
        sampledWaveform = calculate_sampled_waveform(waveform, frequencyValue, std::nullopt, numberPointsSampledWaveforms);
    }
    auto harmonics = calculate_harmonics_data(sampledWaveform, frequencyValue, trimHarmonics, harmonicAmplitudeThreshold, numberPointsSampledWaveforms);
    return calculate_processed_data(harmonics, waveform, includeAdvancedData, processed, numberPointsSampledWaveforms);
}

ProcessedWaveform WaveformProcessor::calculate_processed_data(Harmonics harmonics,
                                            Waveform waveform,
                                            bool includeAdvancedData,
                                            std::optional<ProcessedWaveform> processed,
                                            size_t numberPointsSampledWaveforms) {
    auto sampledDataToProcess = waveform;

    if (waveform.get_time() && waveform.get_data().size() < numberPointsSampledWaveforms) {
        auto frequency = harmonics.get_frequencies()[1];
        sampledDataToProcess = calculate_sampled_waveform(waveform, frequency, std::nullopt, numberPointsSampledWaveforms);
    }

    ProcessedWaveform processedResult;
    if (processed) {
        processedResult = processed.value();
    }
    else {
        processedResult = calculate_basic_processed_data(sampledDataToProcess, numberPointsSampledWaveforms);
    }

    {
        // calculate_waveform_average dereferences waveform.get_time() (UB for a data-only
        // waveform, which the schema allows) — for those, average the resampled uniform
        // waveform instead, same as the sampled branch.
        if (is_waveform_sampled(waveform, numberPointsSampledWaveforms) || !waveform.get_time()) {
            processedResult.set_average(std::accumulate(std::begin(sampledDataToProcess.get_data()), std::end(sampledDataToProcess.get_data()), 0.0) /
                             sampledDataToProcess.get_data().size());
        }
        else {
            auto average = calculate_waveform_average(waveform);
            processedResult.set_average(average);
        }
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

OperatingPointExcitation WaveformProcessor::complete_excitation(Waveform currentWaveform, Waveform voltageWaveform, double switchingFrequency, std::string name, bool trimHarmonics, double harmonicAmplitudeThreshold, size_t numberPointsSampledWaveforms) {
    if (switchingFrequency <= 0 || !std::isfinite(switchingFrequency)) {
        throw std::invalid_argument("complete_excitation: Invalid switchingFrequency: " + std::to_string(switchingFrequency));
    }

    OperatingPointExcitation excitation;
    excitation.set_frequency(switchingFrequency);

    SignalDescriptor current;
    auto currentProcessed = calculate_processed_data(currentWaveform, switchingFrequency, true, std::nullopt, trimHarmonics, harmonicAmplitudeThreshold, numberPointsSampledWaveforms);
    auto sampledCurrentWaveform = calculate_sampled_waveform(currentWaveform, switchingFrequency, std::nullopt, numberPointsSampledWaveforms);
    auto currentHarmonics = calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency, trimHarmonics, harmonicAmplitudeThreshold, numberPointsSampledWaveforms);
    // Store the resampled (power-of-2) waveform so downstream consumers
    // (MagnetizingInductance, harmonics derivation, FFT-based pipelines)
    // get the standardized waveform contract MKF expects. The raw
    // analytical waveform is not necessarily a power-of-2 length (e.g.
    // DAB/PSFB/PSHB use 2*N+1 = 513 samples, AHB has variable size due to
    // discontinuity duplicates), and the size-check gates in
    // MagnetizingInductance.cpp would otherwise throw.
    current.set_waveform(sampledCurrentWaveform);
    current.set_processed(currentProcessed);
    current.set_harmonics(currentHarmonics);
    excitation.set_current(current);
    SignalDescriptor voltage;
    auto voltageProcessed = calculate_processed_data(voltageWaveform, switchingFrequency, true, std::nullopt, trimHarmonics, harmonicAmplitudeThreshold, numberPointsSampledWaveforms);
    auto sampledVoltageWaveform = calculate_sampled_waveform(voltageWaveform, switchingFrequency, std::nullopt, numberPointsSampledWaveforms);
    auto voltageHarmonics = calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency, trimHarmonics, harmonicAmplitudeThreshold, numberPointsSampledWaveforms);
    voltage.set_waveform(sampledVoltageWaveform);
    voltage.set_processed(voltageProcessed);
    voltage.set_harmonics(voltageHarmonics);
    excitation.set_voltage(voltage);
    excitation.set_name(name);
    return excitation;
}

} // namespace OpenMagnetics
